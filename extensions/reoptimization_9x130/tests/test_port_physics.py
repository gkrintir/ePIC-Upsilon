"""Independent numerical and response-contract checks; no production ROOT is needed.

Run from the framework's parent with YPSILON_TEST_FRAMEWORK set if necessary.
These are software tests, not validation of the projected data or generator.
"""
from copy import deepcopy
import json
import math
import os
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad

sys.dont_write_bytecode = True
PACKAGE = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ.get(
    "YPSILON_TEST_FRAMEWORK",
    "tmp/note_update_review/production/Ypsilon-Fit-Framework-ProtonTagged",
)).resolve()
sys.path.insert(0, str(SOURCE / "scripts"))
sys.path.insert(0, str(SOURCE / "workflow"))
from reoptimization.data import integrated_flux, observable_cross_section_scaling, state_settings
from reoptimization.response import evaluate_edges, candidate_summary, optimize


def configuration():
    return json.loads((PACKAGE / "reoptimization_9x130.json").read_text())


def independent_flux(q2max):
    """Numerically integrate EPA in y and log(Q²), independent of dN/dW code."""
    alpha = 1.0 / 137.035999084
    electron_mass = 0.00051099895
    proton_mass = 0.9382720813
    denominator = 4.0 * 9.0 * 130.0
    lower_y = (20.0**2 - proton_mass**2) / denominator
    upper_y = 2.0 / (1.0 + math.sqrt(1.0 + 4.0 * electron_mass**2 / q2max))

    def integrated_virtuality(y):
        q2min = electron_mass**2 * y**2 / (1.0 - y)
        if q2min >= q2max:
            return 0.0
        return quad(
            lambda log_q2: alpha / (2.0 * math.pi) * (
                (1.0 + (1.0-y)**2) / y
                - 2.0 * electron_mass**2 * y / math.exp(log_q2)),
            math.log(q2min), math.log(q2max), epsabs=1e-12, epsrel=1e-11,
        )[0]

    return quad(integrated_virtuality, lower_y, upper_y,
                epsabs=1e-12, epsrel=1e-10, limit=200)[0]


def sample(truth, reco, selected=None, weight=1.0):
    truth, reco = np.asarray(truth, float), np.asarray(reco, float)
    size = len(truth)
    return {
        "truth_pt2": truth.copy(), "truth_abs_t": truth.copy(),
        "reco_pt2": reco.copy(), "weight": np.full(size, weight),
        "truth_selected": np.ones(size, dtype=bool),
        "reco_selected": np.ones(size, dtype=bool) if selected is None else np.asarray(selected, bool),
        "is_train": np.arange(size) % 2 == 0, "train_fraction": 0.5,
        "source_lumi": 10.0, "target_lumi": 1.0,
        "metadata": {"synthetic_test_input": True},
    }


class PortPhysicsTests(unittest.TestCase):
    def test_all_energy_normalization_settings_are_explicit(self):
        cfg = configuration()
        self.assertEqual((cfg["electron_energy"], cfg["proton_energy"]), (9.0, 130.0))
        self.assertEqual(cfg["w_range"], [20.0, 70.0])
        self.assertEqual(cfg["target_lumi_invfb"], 1.0)
        expected = {
            "1S": (119.685, .0238, 1.0),
            "2S": (29.297, .0191, .01),
            "3S": (23.447, .0218, .01),
        }
        for state, (rate, br, q2) in expected.items():
            with self.subTest(state=state):
                declared = cfg["states"][state]
                self.assertTrue({"state", "mass", "sigma", "visible_xsec_fb", "br_ee",
                                 "q2_max", "expected_source_lumi_invfb"}.issubset(declared))
                resolved = state_settings(declared)
                self.assertEqual((resolved["visible_xsec_fb"], resolved["br_ee"], resolved["q2_max"]),
                                 (rate, br, q2))
                self.assertEqual(declared["expected_source_lumi_invfb"], 10.0)
                # All-generated source normalization is reduced by luminosity only.
                self.assertAlmostEqual(rate * 10.0 * (cfg["target_lumi_invfb"] / 10.0), rate)
        self.assertEqual(cfg["optimization"]["fit_range"], [0.0, 1.0])
        self.assertEqual(cfg["optimization"]["reco_control_max"], 5.0)
        self.assertEqual(cfg["stress_test_modes"], ["shape_only", "tail_only", "combined"])

    def test_flux_matches_benchmark_and_independent_double_integration(self):
        cfg = configuration()
        for state, expected in (("1S", 0.06880106724446178),
                                ("2S", 0.04959731276197783),
                                ("3S", 0.04959731276197783)):
            with self.subTest(state=state):
                value = integrated_flux(cfg, cfg["states"][state])
                self.assertAlmostEqual(value, expected, delta=2e-11)
                self.assertAlmostEqual(value, independent_flux(cfg["states"][state]["q2_max"]),
                                       delta=2e-11)

    def test_w_flux_additivity_and_rapidity_primary_convention(self):
        cfg = configuration()
        for state in ("1S", "2S", "3S"):
            scfg = cfg["states"][state]
            flux = integrated_flux(cfg, scfg)
            w = observable_cross_section_scaling(cfg, scfg, "w", [20, 30, 40, 50, 60, 70])
            self.assertAlmostEqual(sum(w["bin_flux"]), flux, delta=2e-11)
            self.assertTrue(w["photon_flux_applied"])
            self.assertFalse(w["bin_width_applied"])
            y = observable_cross_section_scaling(cfg, scfg, "y", [-1, 0, 2, 3])
            self.assertFalse(y["photon_flux_applied"])
            self.assertTrue(y["bin_width_applied"])
            np.testing.assert_allclose(y["factors"],
                                       1 / (1000 * scfg["br_ee"] * np.array([1, 2, 1])))
            np.testing.assert_allclose(y["alternate_factors"], np.array(y["factors"]) / flux)

    def test_unsupported_fitted_truth_bin_is_rejected(self):
        cfg = configuration()
        truth = np.repeat([.1, .4, .8], 400)
        data = sample(truth, truth, truth < .6, weight=.1)
        result = evaluate_edges(data, "abs_t_via_pt2", [0, .3, .6, 1], cfg, split="all")
        self.assertEqual(result["metadata"]["all_missed_truth_columns"], [2])
        summary = candidate_summary(result, cfg, "abs_t_via_pt2")
        self.assertFalse(summary["passes"])
        self.assertIn("full_column_rank", summary["failure_reasons"])
        self.assertIn("conditioning", summary["failure_reasons"])
        cfg["optimization"]["candidate_edges"] = [0, .3, .6, 1]
        self.assertEqual(optimize(data, "abs_t_via_pt2", cfg)["status"],
                         "NO_FEASIBLE_SLOPE_BINNING")

    def test_outside_truth_retained_as_feed_in_and_reco_control(self):
        cfg = configuration()
        # Last three rows are outside the truth fit interval. Two migrate in;
        # the middle one supplies the reconstructed 1-5 control interval.
        data = sample([.1, .4, .8, 1.5, 2, 6], [.1, .4, .8, .1, 2, .2])
        result = evaluate_edges(data, "abs_t_via_pt2", [0, .3, .6, 1], cfg, split="all")
        np.testing.assert_array_equal(result["truth_expected"], [1, 1, 1])
        np.testing.assert_array_equal(result["feed_in"], [2, 0, 0, 1])
        np.testing.assert_array_equal(result["reco_expected"], [3, 1, 1, 1])
        np.testing.assert_array_equal(result["reco_edges"], [0, .3, .6, 1, 5])
        np.testing.assert_array_equal(result["response"].sum(axis=0), [1, 1, 1])
        self.assertEqual(result["metadata"]["control_reco_indices"], [3])
        self.assertEqual(result["metadata"]["generated_truth_expected_in_fit_range"], 3)
        # Unobserved tails stay outside the reduced measurement; no yield repair.
        dropped_tail = deepcopy(data)
        dropped_tail["reco_selected"][3:] = False
        reduced = evaluate_edges(dropped_tail, "abs_t_via_pt2", [0, .3, .6, 1], cfg, split="all")
        np.testing.assert_array_equal(reduced["truth_expected"], result["truth_expected"])
        np.testing.assert_array_equal(reduced["feed_in"], [0, 0, 0, 0])
        self.assertEqual(len(dropped_tail["weight"]), 6)

    def test_distinct_truth_pt2_and_true_t_are_not_interchanged(self):
        cfg = configuration()
        data = sample([.1, .4, .8], [.1, .4, .8])
        data["truth_abs_t"] = np.array([.4, .8, 2.0])
        pt2 = evaluate_edges(data, "pt2", [0, .3, .6, 1], cfg, split="all")
        abs_t = evaluate_edges(data, "abs_t_via_pt2", [0, .3, .6, 1], cfg, split="all")
        np.testing.assert_array_equal(pt2["truth_expected"], [1, 1, 1])
        np.testing.assert_array_equal(abs_t["truth_expected"], [0, 1, 1])
        np.testing.assert_array_equal(abs_t["feed_in"], [0, 0, 1, 0])


if __name__ == "__main__":
    unittest.main()
