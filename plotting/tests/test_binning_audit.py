import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy.stats import poisson

from audit_projected_binning import audit_bin_counts, audit_projected_binning
from upsilon_esr.common import Config, DEFAULT_COVERAGE


def build_inputs(directory):
    root = Path(directory)
    for state, total in (("1S", 40.), ("2S", 20.), ("3S", 10.)):
        base = root / state / "artifacts"
        perf, cross, diff = [base / p for p in ("performance/9x275", "cross_sections", "differential")]
        for p in (perf, cross, diff):
            p.mkdir(parents=True)
        q2 = 1. if state == "1S" else .01
        audit = dict(state=state, energy_label="9x275 GeV",
                     truth_selection=dict(Q2_GeV2=[0, q2], W_GeV=[20, 30]),
                     reco_selection={"mass_GeV": [8, 12], "fixture_only": True},
                     statistical_configuration={"active_luminosity_invfb": 10.})
        (perf / "performance_audit.json").write_text(json.dumps(audit))
        (perf / "mass_fit.json").write_text(json.dumps(dict(
            state=state, selected_projected_yield=total,
            fitted_yield=9999.,  # Must never be used as the selected-yield budget.
            statistical_configuration={"active_luminosity_invfb": 10.})))
        pd.DataFrame(dict(low=[8., 10.], high=[10., 12.], projected_yield=[total*.4, total*.6])).to_csv(perf / "mass_bins.csv", index=False)
        counts, flux = total*np.array([.4, .6]), np.array([.01, .02])
        truth, br = 2*counts, .0234
        sigma = truth/(1000*10*br*flux)
        pd.DataFrame(dict(state=state, W_low=[20., 25.], W_high=[25., 30.], W_center=[22.5, 27.5],
                          n_reco_weighted_visible=counts, n_true_weighted_visible=truth,
                          sigma_gamma_p_pb=sigma, sigma_gamma_p_truth_closure_pb=sigma,
                          flux_int=flux, reco_count_to_sigma_pb=sigma/counts,
                          yield_lumi_invfb=10., output_convention="inclusive_upsilon",
                          plotted=True)).to_csv(cross / f"9x275_sigma_gamma_p_{state.lower()}.csv", index=False)
        pd.DataFrame([dict(state=state, br_ee=br)]).to_csv(cross / f"9x275_sigma_gamma_p_{state.lower()}_normalization.csv", index=False)
        pd.DataFrame([dict(state=state, q2_max_flux=q2, q2_cut_events=q2, w_min=20., w_max=30., fixture_only=True)]).to_csv(diff / f"9x275_{state.lower()}_metadata.csv", index=False)
        for var, edges, fractions in (("pT2", [0., .5, 1., 5.], [.2, .3, .5]),
                                      ("y", [-1., 0., 1.], [.3, .7])):
            edges, counts = np.array(edges), total*np.array(fractions)
            factor = 1./np.diff(edges)
            pd.DataFrame(dict(state=state, x_low=edges[:-1], x_high=edges[1:],
                              x_center=.5*(edges[1:]+edges[:-1]), bin_width=np.diff(edges),
                              n_reco_visible=counts, dsigma_dx_pb=counts*factor,
                              reco_count_to_dsigma_dx_pb=factor, passes_min_eff=True,
                              yield_lumi_invfb=10., output_convention="inclusive_upsilon")).to_csv(
                                  diff / f"9x275_{state.lower()}_{state}_{var}.csv", index=False)
    return Config(runs={s:root/s for s in ("1S", "2S", "3S")}, lumi_invfb=2.5)


class ProjectedBinningAuditTests(unittest.TestCase):
    def test_common_selection_fixture_closes_but_is_not_certified(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            result = audit_projected_binning(cfg, write_outputs=False)
            self.assertEqual(result["status"], "DIAGNOSTICS_ONLY")
            self.assertFalse(result["errors"])
            self.assertFalse(any(c["status"] == "FAIL" for c in result["checks"]))
            self.assertEqual(result["cross_observable_selection_match"], "NOT_VERIFIED")
            for state, expected in (("1S", 10.), ("2S", 5.), ("3S", 2.5)):
                summary = result["states"][state]
                self.assertAlmostEqual(summary["mass_expected_selected_yield"], expected)
                for obs in ("W", "pT2", "y"):
                    self.assertAlmostEqual(summary["observables"][obs]["full_expected_reco_sum"], expected)
                    self.assertAlmostEqual(summary["observables"][obs]["full_to_mass_expected_yield_ratio"], 1.)
                self.assertAlmostEqual(summary["observables"]["pT2"]["displayed_expected_reco_sum"], expected*.5)
                fit = summary["observables"]["pT2"]["slope_fit_count_diagnostic"]
                self.assertFalse(fit["has_minimum_three_supported_bins"])
                self.assertAlmostEqual(fit["probability_zero_events"], np.exp(-expected*.5))
                self.assertFalse(fit["interval_coverage_verified"])

    def test_mass_source_luminosity_mismatch_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            p = cfg.runs["2S"] / "artifacts/performance/9x275/mass_fit.json"
            raw = json.loads(p.read_text())
            raw["statistical_configuration"]["active_luminosity_invfb"] = 20.
            p.write_text(json.dumps(raw))
            result = audit_projected_binning(cfg, write_outputs=False)
            mismatch = [c for c in result["checks"] if c["state"] == "2S" and c["check"] == "mass_fit_vs_performance_source_lumi"]
            self.assertEqual(mismatch[0]["status"], "FAIL")
            self.assertEqual(result["status"], "INPUT_ISSUES")
            summary = result["states"]["2S"]
            self.assertEqual(summary["mass_expected_selected_yield"], 2.5)
            self.assertEqual(summary["observables"]["W"]["full_expected_reco_sum"], 5.)
            self.assertEqual(summary["observables"]["W"]["cross_observable_selection_match"], "NOT_VERIFIED")

    def test_invalid_mass_edges_and_selected_yield_are_reported(self):
        for defect in ("nonfinite", "negative_width", "gap", "negative_selected", "nan_selected"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                cfg = build_inputs(directory)
                directory = cfg.runs["1S"] / "artifacts/performance/9x275"
                if defect.endswith("selected"):
                    p = directory / "mass_fit.json"
                    data = json.loads(p.read_text())
                    data["selected_projected_yield"] = -1. if defect == "negative_selected" else float("nan")
                    p.write_text(json.dumps(data))
                else:
                    p = directory / "mass_bins.csv"
                    data = pd.read_csv(p)
                    if defect == "nonfinite":
                        data.loc[0, "low"] = np.nan
                    elif defect == "negative_width":
                        data.loc[0, "high"] = 7.
                    else:
                        data.loc[1, "low"] = 10.1
                    data.to_csv(p, index=False)
                result = audit_projected_binning(cfg, write_outputs=False)
                self.assertEqual(result["status"], "INPUT_ISSUES")
                self.assertIsNone(result["states"]["1S"]["mass_expected_selected_yield"])
                self.assertIn("W", result["states"]["1S"]["observables"])
                self.assertTrue(any(e["state"] == "1S" and e["observable"] == "mass" for e in result["errors"]))

    def test_missing_w_selection_metadata_is_not_verified(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            result = audit_projected_binning(cfg, write_outputs=False)
            check = [c for c in result["checks"] if c["state"] == "1S" and c["observable"] == "W" and c["check"] == "recorded_Q2_W_vs_performance"][0]
            self.assertEqual(check["status"], "NOT_VERIFIED")
            self.assertIn("q2_cut_events", check["evidence"]["missing_normalization_columns"])

    def test_sparse_and_zero_count_probabilities_and_intervals(self):
        means = np.array([0., .25, 1., 5.])
        frame = pd.DataFrame(dict(W_low=[0., 1., 2., 3.], W_high=[1., 2., 3., 4.],
                                  W_center=[.5, 1.5, 2.5, 3.5],
                                  n_reco_weighted_visible=4*means,
                                  sigma_gamma_p_pb=8*means, reco_count_to_sigma_pb=2.,
                                  yield_lumi_invfb=10., plotted=True))
        result = audit_bin_counts(frame, state="1S", observable="W", target_lumi=2.5, coverage=DEFAULT_COVERAGE)
        np.testing.assert_allclose(result.expected_reco_count, means)
        np.testing.assert_allclose(result.probability_zero, [1., np.exp(-.25), np.exp(-1), np.exp(-5)])
        np.testing.assert_allclose(result.probability_below_5, poisson.cdf(4, means))
        self.assertAlmostEqual(result.count_interval_high.iloc[0], 1.841021645, places=8)
        self.assertAlmostEqual(result.upper_cross_section.iloc[0], 8*1.841021645, places=7)
        self.assertTrue(result.supported.iloc[0])
        self.assertFalse(result.displayed.iloc[0])

    def test_bad_count_factor_preserves_diagnostic_and_flags_renderer_error(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            p = cfg.runs["1S"] / "artifacts/cross_sections/9x275_sigma_gamma_p_1s.csv"
            raw = pd.read_csv(p)
            raw.loc[0, "reco_count_to_sigma_pb"] *= 2
            raw.to_csv(p, index=False)
            result = audit_projected_binning(cfg, write_outputs=False)
            self.assertEqual(result["status"], "INPUT_ISSUES")
            bad = [r for r in result["bin_records"] if r["state"] == "1S" and r["observable"] == "W" and r["bin_index"] == 0][0]
            self.assertEqual(bad["count_factor_closure"], "FAIL")
            self.assertFalse(bad["valid"])
            self.assertTrue(any("does not close" in e["message"] for e in result["errors"]))

    def test_missing_observable_does_not_abort_other_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            (cfg.runs["3S"] / "artifacts/differential/9x275_3s_3S_y.csv").unlink()
            output = Path(directory) / "audit"
            result = audit_projected_binning(cfg, output)
            self.assertEqual(result["status"], "INPUT_ISSUES")
            self.assertIn("W", result["states"]["3S"]["observables"])
            self.assertTrue(any(e["state"] == "3S" and e["observable"] == "y" for e in result["errors"]))
            self.assertTrue((output / "projected_bins.csv").is_file())
            saved = json.loads((output / "projected_binning_audit.json").read_text())
            self.assertEqual(saved["status"], "INPUT_ISSUES")

    def test_unknown_boolean_is_an_input_failure_not_a_silent_mask(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            p = cfg.runs["1S"] / "artifacts/cross_sections/9x275_sigma_gamma_p_1s.csv"
            raw = pd.read_csv(p)
            raw["plotted"] = raw["plotted"].astype(object)
            raw.loc[0, "plotted"] = "typo"
            raw.to_csv(p, index=False)
            result = audit_projected_binning(cfg, write_outputs=False)
            self.assertEqual(result["status"], "INPUT_ISSUES")
            check = [c for c in result["checks"] if c["state"] == "1S" and c["check"] == "strict_boolean_plotted"][0]
            self.assertEqual(check["status"], "FAIL")
            self.assertEqual(check["evidence"]["invalid_values"], ["typo"])

    def test_missing_explicit_factor_warns_that_closure_is_derived(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg = build_inputs(directory)
            p = cfg.runs["1S"] / "artifacts/differential/9x275_1s_1S_pT2.csv"
            pd.read_csv(p).drop(columns="reco_count_to_dsigma_dx_pb").to_csv(p, index=False)
            result = audit_projected_binning(cfg, write_outputs=False)
            check = [c for c in result["checks"] if c["state"] == "1S" and c["observable"] == "pT2" and c["check"] == "reco_count_factor_closure"][0]
            self.assertEqual(check["status"], "WARN")
            self.assertFalse(check["evidence"]["independent_factor_present"])
            rows = [r for r in result["bin_records"] if r["state"] == "1S" and r["observable"] == "pT2"]
            self.assertTrue(all(r["count_factor_closure"] == "DERIVED_NOT_INDEPENDENT" for r in rows))

    def test_corrupt_counts_fail_even_when_renderer_can_drop_them(self):
        for bad_value in (np.nan, -1.):
            with self.subTest(value=bad_value), tempfile.TemporaryDirectory() as directory:
                cfg = build_inputs(directory)
                p = cfg.runs["1S"] / "artifacts/cross_sections/9x275_sigma_gamma_p_1s.csv"
                raw = pd.read_csv(p)
                raw.loc[0, "n_reco_weighted_visible"] = bad_value
                raw.to_csv(p, index=False)
                result = audit_projected_binning(cfg, write_outputs=False)
                self.assertEqual(result["status"], "INPUT_ISSUES")
                check = [c for c in result["checks"] if c["state"] == "1S" and c["observable"] == "W" and c["check"] == "all_expected_reco_counts_finite_nonnegative"][0]
                self.assertEqual(check["status"], "FAIL")
                self.assertEqual(check["evidence"]["invalid_bin_indices"], [0])

    def test_invalid_correction_requires_explicit_mask_and_empty_closure_warns(self):
        for masked in (False, True):
            with self.subTest(masked=masked), tempfile.TemporaryDirectory() as directory:
                cfg = build_inputs(directory)
                p = cfg.runs["1S"] / "artifacts/cross_sections/9x275_sigma_gamma_p_1s.csv"
                raw = pd.read_csv(p)
                raw["reco_count_to_sigma_pb"] = np.nan
                raw["plotted"] = not masked
                raw.to_csv(p, index=False)
                result = audit_projected_binning(cfg, write_outputs=False)
                self.assertEqual(result["status"], "DIAGNOSTICS_ONLY" if masked else "INPUT_ISSUES")
                checks = {c["check"]: c for c in result["checks"] if c["state"] == "1S" and c["observable"] == "W"}
                self.assertEqual(checks["active_bins_have_valid_central_and_correction"]["status"], "PASS" if masked else "FAIL")
                self.assertEqual(checks["reco_count_factor_closure"]["status"], "WARN")
                self.assertEqual(checks["reco_count_factor_closure"]["evidence"]["n_tested_bins"], 0)
                self.assertEqual(checks["renderer_projection_agreement"]["status"], "WARN")


if __name__ == "__main__":
    unittest.main()
