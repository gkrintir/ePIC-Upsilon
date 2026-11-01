"""Independent selection, derived-coordinate and sparse-bin checks for final maps."""
import copy
import unittest

import numpy as np

from diagnostics.core import Accumulator
from diagnostics.final_maps import FinalMapAccumulator, final_selection_masks, PROTON_MASS_GEV
from original_selection_reference import selection_masks
from test_core import arrays, sample, settings


def accumulate(data, spec=None, cfg=None, chunks=None):
    acc = FinalMapAccumulator(spec or sample(), cfg or settings())
    if chunks is None:
        acc.consume(data)
    else:
        for part in chunks:
            acc.consume({k: v[part] for k, v in data.items()})
    return acc.finish()


def map_rows(result, name):
    return [r for r in result["rows"] if r["map"] == name]


class FinalMapTests(unittest.TestCase):
    def test_exact_producer_predicate_for_all_states_and_boundaries(self):
        for state, cap in [("1S", 1.), ("2S", .01), ("3S", .01)]:
            data = arrays(30, state)
            data["true_mass"][:4] += [-.2, .2, -.21, .21]
            data["true_Q2"][4:9] = [0, cap, np.nextafter(cap, np.inf), -999, np.nan]
            data["true_W"][9:13] = [20, 100, 19.999, np.inf]
            data["reco_mass"][13:19] = [8, 12, 7.99, 12.01, -999, np.nan]
            data["reco_W"][19:23] = [20, 100, -999, 100.01]
            data["sigma_epz"][23:] = [16.2, 19.8, np.nextafter(16.2, 18),
                                     np.nextafter(19.8, 18), -999, np.nan, 18]
            expected = selection_masks(data, state, 9., cap, 20., 100., .2)
            actual = final_selection_masks(data, sample(state, cap), settings())
            for a, b in zip(actual, expected):
                np.testing.assert_array_equal(a, b)
            result = accumulate(data, sample(state, cap))
            self.assertEqual(result["audit"]["eligible_truth"], int(expected[0].sum()))
            self.assertEqual(result["audit"]["selected"], int(expected[1].sum()))

    def test_pt_y_counts_identical_to_old_diagnostics(self):
        rng = np.random.default_rng(42872)
        data = arrays(1200)
        data["true_pT"] = rng.uniform(-.1, 3.1, 1200)
        data["true_y"] = rng.uniform(-1.2, 3.2, 1200)
        data["true_Q2"] = 10 ** rng.uniform(-7, 0, 1200)
        data["reco_mass"][::4] = -999
        data["reco_W"][::5] = 101
        data["sigma_epz"][::7] = 22
        data["true_pT"][:4] = [np.nan, -999, 3, 0]
        data["true_y"][4:8] = [np.inf, -999, -1, 3]
        old = Accumulator(sample(), settings())
        old.consume(data)
        old_rows = [r for r in old.finish()["tables"]["maps"]
                    if r["q2_band"] == "all" and r["y_variable"] == "true_y"]
        new_rows = map_rows(accumulate(data), "pt_y")
        self.assertEqual(len(old_rows), len(new_rows))
        for a, b in zip(old_rows, new_rows):
            for field in ("x_low", "x_high", "y_low", "y_high", "denominator", "selected",
                          "efficiency", "wilson_low", "wilson_high"):
                self.assertEqual(a[field], b[field])

    def test_x_definitions_and_nominal_mass_independent_of_stored_mass(self):
        data = arrays(2)
        data["true_Q2"][:] = [.01, .1]
        data["true_W"][:] = [20, 100]
        data["true_mass"] += [.1, -.1]
        data["reco_mass"][1] = -999
        cfg = settings()
        cfg["final_2d"] = dict(xb_edges=[1e-6, 1.5e-5, 3e-5],
                               xv_edges=[.008, .02, .25], q2_edges=[.001, .03, 1])
        result = accumulate(data, cfg=cfg)
        xb = [r for r in map_rows(result, "xb_q2") if r["denominator"]]
        xv = [r for r in map_rows(result, "xv_q2") if r["denominator"]]
        self.assertEqual(len(xb), 2)
        self.assertEqual(len(xv), 2)
        # Explicit independent expectations: W20 Q2=.01 is high-x, low-Q2;
        # W100 Q2=.1 is low-x, high-Q2. Mass-window fluctuations do not shift xV.
        for rows, split in [(xb, 1.5e-5), (xv, .02)]:
            low = next(r for r in rows if r["x_high"] == split)
            high = next(r for r in rows if r["x_low"] == split)
            self.assertEqual((low["denominator"], low["selected"], low["y_low"]), (1, 0, .03))
            self.assertEqual((high["denominator"], high["selected"], high["y_high"]), (1, 1, .03))
        self.assertAlmostEqual(result["audit"]["x_definitions"]["proton_mass_GeV"], PROTON_MASS_GEV)

    def test_zero_q2_and_positive_log_underflow_are_not_clipped(self):
        data = arrays(4)
        data["true_Q2"][:] = [0, 1e-12, 1e-10, 1]
        result = accumulate(data)
        self.assertEqual(result["audit"]["eligible_truth"], 4)
        self.assertEqual(result["audit"]["zero_q2_eligible_truth"], 1)
        self.assertEqual(sum(r["denominator"] for r in map_rows(result, "pt_y")), 4)
        flow = next(f for f in result["audit"]["histogram_flows"] if f["map"] == "pt_q2")
        self.assertEqual(flow["truth_in_histogram"], 2)
        self.assertEqual(flow["truth_nonpositive_log_coordinates"], 1)
        self.assertEqual(flow["truth_positive_log_underflow"], 1)
        self.assertEqual(flow["truth_valid_outside_any_axis"], 2)
        self.assertEqual(flow["truth_invalid_coordinates"], 0)

    def test_xv_uses_fixed_nominal_mass_and_retains_proton_mass_term(self):
        data = arrays(2)
        data["true_W"][:] = 20.
        data["true_Q2"][:] = .01
        data["true_mass"] += [-.1, .1]
        # Tight interval separates the requested proton-mass-corrected proxy
        # from both a W^2+Q2 denominator and fluctuating stored masses.
        nominal = (9.46030**2 + .01) / (400. + .01 - .9382720813**2)
        cfg = settings()
        cfg["final_2d"] = {"xv_edges": [nominal * .9999, nominal * 1.0001]}
        result = accumulate(data, cfg=cfg)
        self.assertEqual(sum(r["denominator"] for r in map_rows(result, "xv_q2")), 2)

    def test_x_maps_have_no_pt_or_rapidity_requirement(self):
        data = arrays(5)
        data["true_pT"][:] = [np.nan, -999, 4, -.1, 900]
        data["true_y"][:] = [-999, np.nan, 4, 900, -2]
        result = accumulate(data)
        self.assertEqual(sum(r["denominator"] for r in map_rows(result, "pt_y")), 0)
        for name in ["xb_q2", "xv_q2"]:
            self.assertEqual(sum(r["denominator"] for r in map_rows(result, name)), 5)
            self.assertEqual(sum(r["selected"] for r in map_rows(result, name)), 5)

    def test_high_state_q2_cap_is_inserted_and_included(self):
        cfg = settings()
        cfg["final_2d"] = {"q2_edges": [1e-8, .001, .1, 1]}
        data = arrays(3, "2S")
        data["true_Q2"][:] = [.001, .01, np.nextafter(.01, 1.)]
        result = accumulate(data, sample("2S", .01), cfg)
        self.assertEqual(result["audit"]["edges"]["true_Q2"], [1e-8, .001, .01])
        self.assertEqual(sum(r["denominator"] for r in map_rows(result, "pt_q2")), 2)
        self.assertEqual(result["audit"]["eligible_truth"], 2)

    def test_zero_success_and_empty_are_distinct(self):
        data = arrays(8)
        data["reco_mass"][:] = -999
        result = accumulate(data)
        populated = next(r for r in map_rows(result, "pt_y") if r["denominator"])
        empty = next(r for r in map_rows(result, "pt_y") if not r["denominator"])
        self.assertEqual(populated["efficiency"], 0)
        self.assertAlmostEqual(populated["wilson_high"], 1/9)
        self.assertIsNone(empty["efficiency"])
        self.assertIsNone(empty["wilson_high"])
        self.assertFalse(empty["sparse"])

    def test_chunking_and_extra_weights_do_not_change_counts(self):
        data = arrays(35)
        data["true_pT"] = np.linspace(-.1, 3.1, 35)
        data["true_Q2"] = np.geomspace(1e-12, 1, 35)
        data["true_W"] = np.linspace(20, 100, 35)
        data["reco_mass"][::3] = -999
        whole = accumulate(data)
        split = accumulate(data, chunks=[slice(0, 5), slice(5, 17), slice(17, 35)])
        self.assertEqual(whole["rows"], split["rows"])
        self.assertEqual(whole["audit"]["histogram_flows"], split["audit"]["histogram_flows"])
        weighted = copy.deepcopy(data)
        weighted["weight"] = np.linspace(-100, 100, 35)
        self.assertEqual(whole["rows"], accumulate(weighted)["rows"])

    def test_refuses_invalid_log_edges_and_post_finish_mutation(self):
        for override in [{"q2_edges": [0, .01, 1]}, {"xb_edges": [-1, .1]},
                         {"xv_edges": [.1, .1]}, {"q2_edges": [1, 2]}]:
            cfg = settings()
            cfg["final_2d"] = override
            with self.assertRaises(ValueError):
                FinalMapAccumulator(sample(), cfg)
        acc = FinalMapAccumulator(sample(), settings())
        result = acc.finish()
        self.assertEqual(result["audit"]["raw_rows"], 0)
        self.assertEqual(result["audit"]["status"], "EMPTY_INPUT")
        with self.assertRaises(RuntimeError):
            acc.consume(arrays())


if __name__ == "__main__":
    unittest.main()
