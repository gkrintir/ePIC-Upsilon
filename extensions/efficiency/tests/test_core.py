"""Independent count and legacy-selection checks; no production inputs needed."""
import copy
import unittest

import numpy as np

from diagnostics.core import Accumulator
from original_selection_reference import MASSES, selection_masks


def settings():
    return dict(proton_selection="inclusive", electron_energy=9., w_range=[20., 100.], true_mass_window=.2,
                q2_split=.01, coverage=.6826894921370859, min_denominator=50,
                pt_edges=np.linspace(0, 3, 16).tolist(),
                w_edges=np.linspace(20, 100, 11).tolist(),
                y_edges=np.linspace(-1, 3, 17).tolist(),
                q2_edges=[0, 1e-6, 1e-5, 1e-4, .001, .01, .03, .1, .3, 1],
                tail_edges=[0, 1, 1.6, 2.2, 2.6, 3],
                strata_w_edges=[20, 40, 60, 80, 100],
                strata_y_edges=[-1, 0, 1, 2, 3])


def sample(state="1S", cap=1.):
    return dict(key="test", label="Synthetic test", state=state,
                id_mode="synthetic", q2_max=cap, path="synthetic.root")


def arrays(n=8, state="1S"):
    return dict(true_mass=np.full(n, MASSES[state]), true_Q2=np.full(n, .001),
                true_W=np.full(n, 50.), true_y=np.full(n, 1.),
                true_pT=np.full(n, .1), reco_mass=np.full(n, MASSES[state]),
                reco_W=np.full(n, 50.), sigma_epz=np.full(n, 18.))


def run(data, spec=None, chunks=None):
    acc = Accumulator(spec or sample(), settings())
    if chunks is None:
        acc.consume(data)
    else:
        for idx in chunks:
            acc.consume({k: v[idx] for k, v in data.items()})
    return acc.finish()


def rows(result, table, band="all", **fields):
    return [r for r in result["tables"][table]
            if r["q2_band"] == band and all(r[k] == v for k, v in fields.items())]


class EfficiencyCoreTests(unittest.TestCase):
    def test_original_producer_selection_and_histograms(self):
        rng = np.random.default_rng(983142)
        for state, cap in [("1S", 1.), ("2S", .01), ("3S", .01)]:
            data = arrays(2000, state)
            data.update(true_mass=MASSES[state] + rng.normal(0, .19, 2000),
                        true_Q2=rng.uniform(-.001, 1.1 * cap, 2000),
                        true_W=rng.uniform(10, 110, 2000),
                        true_y=rng.uniform(-2, 4, 2000),
                        true_pT=rng.uniform(-.2, 3.5, 2000),
                        reco_mass=rng.uniform(7, 13, 2000),
                        reco_W=rng.uniform(10, 110, 2000),
                        sigma_epz=rng.uniform(15, 21, 2000))
            # Exercise exact endpoints, missing-value sentinels and NaNs.
            for key in data:
                data[key][0:3] = [-999., np.nan, np.inf]
            data["true_Q2"][3:7] = [0, .01, cap, np.nextafter(cap, np.inf)]
            data["true_W"][7:11] = [20, 100, np.nextafter(20., 0.), np.nextafter(100., np.inf)]
            data["reco_mass"][11:15] = [8, 12, np.nextafter(8., 0.), np.nextafter(12., np.inf)]
            truth, selected = selection_masks(data, state, 9., cap, 20., 100., .2)
            result = run(data, sample(state, cap))
            for variable, edge_key in [("true_pT", "pt_edges"), ("true_W", "w_edges"), ("true_y", "y_edges")]:
                valid = np.isfinite(data[variable]) & (np.abs(data[variable]) < 900)
                expected_n = np.histogram(data[variable][truth & valid], settings()[edge_key])[0]
                expected_k = np.histogram(data[variable][selected & valid], settings()[edge_key])[0]
                saved = sorted(rows(result, "efficiency", variable=variable), key=lambda x: x["low"])
                np.testing.assert_array_equal([r["denominator"] for r in saved], expected_n)
                np.testing.assert_array_equal([r["selected"] for r in saved], expected_k)

    def test_q2_partition_includes_split_only_in_low(self):
        data = arrays(5)
        data["true_Q2"][:] = [.0099, .01, np.nextafter(.01, 1.), .999, 1.]
        result = run(data)
        counts = {}
        for band in ["all", "low", "high"]:
            counts[band] = sum(r["denominator"] for r in rows(result, "efficiency", band, variable="true_pT"))
        self.assertEqual(counts, {"all": 5, "low": 2, "high": 3})
        for variable in ["true_pT", "true_W", "true_y"]:
            band_rows = {b: sorted(rows(result, "efficiency", b, variable=variable), key=lambda r: r["low"])
                         for b in counts}
            for a, lo, hi in zip(band_rows["all"], band_rows["low"], band_rows["high"]):
                for name in ["denominator", "selected"]:
                    self.assertEqual(a[name], lo[name] + hi[name])

    def test_chunking_preserves_all_numerical_tables(self):
        data = arrays(35)
        data["true_pT"] = np.linspace(-.1, 3.1, 35)
        data["true_Q2"] = np.geomspace(1e-5, 1, 35)
        data["reco_mass"][::3] = -999
        whole = run(data)
        split = run(data, chunks=[slice(0, 5), slice(5, 17), slice(17, 35)])
        self.assertEqual(whole["tables"], split["tables"])

    def test_cutflow_and_nminus1_with_overlapping_failures(self):
        data = arrays(8)
        data["reco_mass"][[1, 4, 5, 7]] = -999
        data["reco_W"][[2, 4, 6, 7]] = -999
        data["sigma_epz"][[3, 5, 6, 7]] = -999
        result = run(data)
        cutflow = {r["stage"]: r for r in rows(result, "cutflow") if r["pt_low"] == 0}
        expected = {"truth": (8, 8, 0), "reco_mass": (8, 4, 4),
                    "reco_W": (4, 2, 2), "sigma_epz": (2, 1, 1)}
        for name, (entered, passed, failed) in expected.items():
            self.assertEqual((cutflow[name]["entered"], cutflow[name]["passed"], cutflow[name]["failed"]),
                             (entered, passed, failed))
            self.assertEqual(cutflow[name]["baseline"], 8)
        check = [r for r in rows(result, "nminus1") if r["pt_low"] == 0 and r["cut"] != "proton_tag"]
        self.assertEqual(len(check), 3)
        for row in check:
            self.assertEqual(row["denominator"], 2)
            self.assertEqual(row["selected"], 1)
            self.assertEqual(row["efficiency"], .5)

    def test_empty_and_zero_success_intervals_are_distinct(self):
        data = arrays(8)
        data["reco_mass"][:] = -999
        result = run(data)
        saved = sorted(rows(result, "efficiency", variable="true_pT"), key=lambda r: r["low"])
        zero_success, empty = saved[0], saved[1]
        self.assertEqual(zero_success["denominator"], 8)
        self.assertEqual(zero_success["efficiency"], 0.)
        self.assertAlmostEqual(zero_success["wilson_high"], 1/9, places=10)
        self.assertEqual(empty["denominator"], 0)
        self.assertIsNone(empty["efficiency"])
        self.assertIsNone(empty["wilson_low"])
        self.assertIsNone(empty["wilson_high"])

    def test_histogram_range_and_last_edge(self):
        data = arrays(7)
        data["true_pT"][:] = [-999, np.nan, -.2, 0, 3, 3.01, 900]
        result = run(data)
        saved = sorted(rows(result, "efficiency", variable="true_pT"), key=lambda r: r["low"])
        self.assertEqual(sum(r["denominator"] for r in saved), 2)
        self.assertEqual(saved[0]["denominator"], 1)
        self.assertEqual(saved[-1]["denominator"], 1)
        flow = next(r for r in result["audit"]["histogram_flows"]
                    if r["q2_band"] == "all" and r["histogram"] == "efficiency/true_pT")
        self.assertEqual(flow["truth_total"], 7)
        self.assertEqual(flow["truth_invalid_coordinates"], 3)
        self.assertEqual(flow["truth_valid_below_any_axis"], 1)
        self.assertEqual(flow["truth_valid_above_any_axis"], 1)
        self.assertEqual(flow["truth_valid_outside_any_axis"], 2)
        self.assertEqual(flow["truth_in_histogram"], 2)

    def test_low_only_cap_never_creates_high_q2_population(self):
        data = arrays(4, "2S")
        data["true_Q2"][:] = [0, .009, .01, .01001]
        result = run(data, sample("2S", .01))
        self.assertEqual(result["audit"]["band_counts"]["all"]["eligible_truth"], 3)
        self.assertEqual(result["audit"]["band_counts"]["low"]["eligible_truth"], 3)
        self.assertNotIn("high", result["audit"]["band_counts"])
        self.assertFalse(rows(result, "efficiency", "high"))

    def test_q2_mixture_decomposes_selected_counts(self):
        data = arrays(5)
        data["true_Q2"][:] = [.001, .01, .02, .1, 1.]
        data["reco_mass"][[1, 4]] = -999
        result = run(data)
        row = next(r for r in rows(result, "q2_mixture", variable="true_pT") if r["low"] == 0)
        self.assertEqual((row["n_all"], row["n_low"], row["n_high"]), (5, 2, 3))
        self.assertEqual((row["k_all"], row["k_low"], row["k_high"]), (3, 1, 2))
        self.assertAlmostEqual(row["high_truth_fraction"], .6)
        self.assertAlmostEqual(row["high_selected_fraction"], 2/3)
        self.assertAlmostEqual(row["low_contribution"], .2)
        self.assertAlmostEqual(row["high_contribution"], .4)
        self.assertAlmostEqual(row["eff_all"], .6)
        self.assertAlmostEqual(row["mixture_residual"], 0.)
        empty = next(r for r in rows(result, "q2_mixture", variable="true_pT") if r["low"] == .2)
        self.assertIsNone(empty["eff_all"])

    def test_weights_do_not_enter_efficiency(self):
        data = arrays(8)
        result = run(data)
        weighted = copy.deepcopy(data)
        weighted["weight"] = np.array([0, -3, 1, 1e8, 4, 9, np.nan, np.inf])
        self.assertEqual(result["tables"], run(weighted)["tables"])

    def test_failure_reasons_partition_all_truth_with_overlaps(self):
        data = arrays(12)
        data["reco_mass"][:] = [-999, np.nan, np.inf, -np.inf, 7.9, 12.1,
                                8, 12, 10, -1, 900, 9.5]
        data["reco_W"][:] = [-999, np.nan, np.inf, -np.inf, 19, 101,
                             20, 100, 50, -1, 900, 55]
        data["sigma_epz"][:] = [-999, np.nan, np.inf, -np.inf, 16, 20,
                                np.nextafter(16.2, 18.), np.nextafter(19.8, 18.),
                                18, -1, 900, 17]
        result = run(data)
        expected = {"nonfinite": 3, "minus_999_sentinel": 1,
                    "finite_lower_rejection": 2, "finite_upper_rejection": 2,
                    "pass": 4}
        for cut in ["reco_mass", "reco_W", "sigma_epz"]:
            failure_rows = [r for r in rows(result, "failure_reasons", cut=cut) if r["pt_low"] == 0]
            self.assertEqual({r["reason"]: r["count"] for r in failure_rows}, expected)
            self.assertEqual(sum(r["count"] for r in failure_rows), 12)
            for row in failure_rows:
                self.assertEqual(row["baseline"], 12)
                self.assertEqual(row["denominator"], 12)
                self.assertAlmostEqual(row["fraction"], row["count"] / 12)
        # These are independent field diagnoses, not exclusive first failures.
        final = rows(result, "efficiency", variable="true_pT")[0]
        self.assertEqual(final["selected"], 4)

    def test_failure_reasons_use_original_epz_endpoint_predicate(self):
        data = arrays(4)
        data["sigma_epz"][:] = [16.2, 19.8, np.nextafter(16.2, 18.), np.nextafter(19.8, 18.)]
        _, selected = selection_masks(data, "1S", 9., 1., 20., 100., .2)
        result = run(data)
        failure_rows = [r for r in rows(result, "failure_reasons", cut="sigma_epz") if r["pt_low"] == 0]
        counts = {r["reason"]: r["count"] for r in failure_rows}
        self.assertEqual(counts["pass"], int(selected.sum()))
        self.assertEqual(sum(counts.values()), 4)

    def test_missing_or_mismatched_branch_rejected(self):
        data = arrays(8)
        del data["sigma_epz"]
        with self.assertRaises((ValueError, KeyError)):
            run(data)
        data = arrays(8)
        data["true_Q2"] = np.ones(7)
        with self.assertRaises((ValueError, TypeError)):
            run(data)


if __name__ == "__main__":
    unittest.main()
