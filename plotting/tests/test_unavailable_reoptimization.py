"""Synthetic tests for explicitly unavailable saved reoptimization cases.

An absent numerical extraction is accepted only when producer metadata says
NO_FEASIBLE_SLOPE_BINNING. Missing/corrupt results are never silently skipped.
"""
import json
import tempfile
import unittest
from unittest import mock
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import test_proton_inputs
from upsilon_esr.common import write_json
from upsilon_esr.reoptimization_inputs import load_case, preflight_reoptimization


UNAVAILABLE = "NO_FEASIBLE_SLOPE_BINNING"


def mark_unavailable(cfg, state="2S", observable="pt2"):
    """Model a complete no-feasible-binning decision with no extraction."""
    base = cfg.reoptimization_run / "artifacts/reoptimization"
    result_path = base / state / observable / "result.json"
    entry = json.loads(result_path.read_text())
    entry.update(status=UNAVAILABLE, overall_status=UNAVAILABLE, complete=True,
                 publication_ready=False, failed_checks=[])
    entry["optimization"] = {"status": UNAVAILABLE, "selected_edges": None,
                             "candidates": [], "reasons": []}
    for key in ("inference", "cross_section_scaling", "cross_section_factor",
                "model", "validation_model"):
        entry.pop(key, None)
    write_json(result_path, entry)
    path = base / "validation_status.csv"
    table = pd.read_csv(path)
    mask = (table.state == state) & (table.observable == observable)
    table.loc[mask, ["status", "overall_status"]] = UNAVAILABLE
    table.loc[mask, "stage_status"] = "NOT_RUN"
    table.loc[mask, "failed_checks"] = "[]"
    table.loc[mask, "optimization_status"] = UNAVAILABLE
    table.loc[mask, "inference_status"] = "NOT_RUN"
    table.loc[mask, "proton_selection"] = json.dumps(entry["proton_selection"])
    table.loc[mask, "proton_selection_mode"] = "required"
    table.to_csv(path, index=False)
    path = base / "differential_results.csv"
    table = pd.read_csv(path)
    table = table.loc[~((table.state == state) & (table.observable == observable))]
    table.to_csv(path, index=False)
    for path in base.glob(f"{state}_{observable}_*covariance*"):
        path.unlink()
    all_path = base / "all_results.json"
    if all_path.is_file():
        values = json.loads(all_path.read_text())
        values[state][observable] = entry
        write_json(all_path, values)
    return entry


class UnavailableReoptimizationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cfg = test_proton_inputs.ProtonInputsTests().reopt_fixture(self.temp.name)
        self.base = self.cfg.reoptimization_run / "artifacts/reoptimization"
        self.warning_context = warnings.catch_warnings()
        self.warning_context.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)
        self.addCleanup(self.warning_context.__exit__, None, None, None)
        self.addCleanup(plt.close, "all")

    def test_specific_exception_preserves_unavailable_contract(self):
        from upsilon_esr.reoptimization_inputs import SavedCaseUnavailable
        mark_unavailable(self.cfg)
        with self.assertRaises(SavedCaseUnavailable) as context:
            load_case(self.cfg, "2S", "pt2")
        metadata = context.exception.metadata
        self.assertEqual(metadata["overall_status"], UNAVAILABLE)
        self.assertEqual(metadata["stage_status"], "NOT_RUN")
        self.assertFalse(metadata["publication_ready"])
        self.assertEqual(metadata["proton_selection_mode"], "required")
        self.assertEqual(context.exception.entry["optimization"]["status"], UNAVAILABLE)
        self.assertTrue(context.exception.sources)
        self.assertTrue(all(len(source["sha256"]) == 64 for source in context.exception.sources))

    def test_preflight_retains_all_cases_and_marks_unavailable_bins(self):
        for state in ("2S", "3S"):
            mark_unavailable(self.cfg, state)
        table = preflight_reoptimization(self.cfg, "differential")
        self.assertEqual(len(table), 6)
        unavailable = table.loc[(table.observable == "pt2") & table.state.isin(["2S", "3S"])]
        self.assertEqual(set(unavailable.overall_status), {UNAVAILABLE})
        self.assertEqual(set(unavailable.stage_status), {"NOT_RUN"})
        self.assertEqual(unavailable.bins.tolist(), [0, 0])
        self.assertFalse(unavailable.available.any())
        rapidity = table.loc[table.observable == "y"]
        self.assertEqual(set(rapidity.overall_status), {"NOT_VALIDATED"})
        self.assertTrue(rapidity.available.all())

    def test_partial_plot_retains_available_points_without_refit_or_fallback(self):
        from upsilon_esr.differential import plot_differential
        for state in ("2S", "3S"):
            mark_unavailable(self.cfg, state)
        with mock.patch("upsilon_esr.differential.project_cross_section", side_effect=AssertionError("legacy reprojection")), \
             mock.patch("upsilon_esr.differential.fit_projected_exponential", side_effect=AssertionError("legacy refit")):
            fig, audit = plot_differential(self.cfg)
        plt.close(fig)
        points = pd.read_csv(self.cfg.output_dir / "upsilon_es_differential_points.csv")
        self.assertEqual(set(points.loc[points.variable == "pT2", "state"]), {"1S"})
        self.assertEqual(set(points.loc[points.variable == "y", "state"]), {"1S", "2S", "3S"})
        self.assertEqual(set(points.loc[points.variable == "y", "overall_status"]), {"NOT_VALIDATED"})
        for state in ("2S", "3S"):
            meta = audit["validation"][state]["pT2"]
            self.assertEqual(meta["overall_status"], UNAVAILABLE)
            self.assertEqual(meta["stage_status"], "NOT_RUN")
            self.assertFalse(meta["publication_ready"])
            self.assertEqual(audit["slope_fits"][state]["overall_status"], UNAVAILABLE)
            self.assertNotIn("b_slope_GeVminus2", audit["slope_fits"][state])
            self.assertFalse((self.cfg.output_dir / f"upsilon_es_{state.lower()}_pt2_fit_bins.csv").exists())
        self.assertEqual(audit["slope_fits"]["1S"]["b_slope_GeVminus2"], 2.7)
        slopes = pd.read_csv(self.cfg.output_dir / "upsilon_es_slopes.csv")
        missing = slopes.loc[slopes.state.isin(["2S", "3S"])]
        self.assertEqual(set(missing.overall_status), {UNAVAILABLE})
        self.assertEqual(set(missing.stage_status), {"NOT_RUN"})
        self.assertTrue(missing.b_slope_GeVminus2.isna().all())
        self.assertFalse(audit["publication_ready"])
        for extension in ("pdf", "png"):
            self.assertTrue((self.cfg.output_dir / ("upsilon_es_differential." + extension)).is_file())

    def test_missing_csv_for_available_case_is_still_fatal(self):
        path = self.base / "differential_results.csv"
        table = pd.read_csv(path)
        table = table.loc[~((table.state == "2S") & (table.observable == "pt2"))]
        table.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "missing 2S/pt2"):
            load_case(self.cfg, "2S", "pt2")

    def test_no_feasible_case_with_remaining_numerical_rows_is_rejected(self):
        from upsilon_esr.reoptimization_inputs import SavedCaseUnavailable
        path = self.base / "differential_results.csv"
        original = pd.read_csv(path)
        mark_unavailable(self.cfg)
        original.to_csv(path, index=False)
        try:
            load_case(self.cfg, "2S", "pt2")
        except SavedCaseUnavailable:
            self.fail("Conflicting numerical rows were silently treated as unavailable")
        except ValueError:
            pass
        else:
            self.fail("Conflicting numerical rows were accepted")

    def test_inconsistent_unavailable_metadata_is_never_skipped(self):
        from upsilon_esr.reoptimization_inputs import SavedCaseUnavailable
        for defect in ("status", "edges", "coarse", "scaling", "proton", "missing_proton", "status_proton", "incomplete"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as root:
                cfg = test_proton_inputs.ProtonInputsTests().reopt_fixture(root)
                base = cfg.reoptimization_run / "artifacts/reoptimization"
                entry = mark_unavailable(cfg)
                if defect == "status":
                    path = base / "validation_status.csv"
                    table = pd.read_csv(path)
                    table.loc[(table.state == "2S") & (table.observable == "pt2"), "overall_status"] = "NOT_VALIDATED"
                    table.to_csv(path, index=False)
                elif defect == "status_proton":
                    path = base / "validation_status.csv"
                    table = pd.read_csv(path)
                    table.loc[(table.state == "2S") & (table.observable == "pt2"), "proton_selection_mode"] = "inclusive"
                    table.to_csv(path, index=False)
                elif defect == "edges":
                    entry["optimization"]["selected_edges"] = [0., .3, 1.]
                elif defect == "coarse":
                    entry["inference"] = {"coarse": {"status": "PASS_CONDITIONAL", "fitted_truth": [1., 2.]}}
                elif defect == "scaling":
                    entry["cross_section_scaling"] = {"factors": [1., 2.]}
                elif defect == "missing_proton":
                    entry.pop("proton_selection")
                elif defect == "proton":
                    entry["proton_selection"]["mode"] = "inclusive"
                else:
                    entry["complete"] = False
                write_json(base / "2S/pt2/result.json", entry)
                try:
                    load_case(cfg, "2S", "pt2")
                except SavedCaseUnavailable:
                    self.fail(f"{defect} mismatch was silently skipped")
                except ValueError:
                    pass
                else:
                    self.fail(f"{defect} mismatch was accepted")

    def test_no_available_pt2_still_renders_available_rapidity(self):
        from upsilon_esr.differential import plot_differential
        for state in ("1S", "2S", "3S"):
            mark_unavailable(self.cfg, state)
        fig, audit = plot_differential(self.cfg)
        plt.close(fig)
        points = pd.read_csv(self.cfg.output_dir / "upsilon_es_differential_points.csv")
        self.assertEqual(set(points.variable), {"y"})
        self.assertEqual(set(points.state), {"1S", "2S", "3S"})
        self.assertTrue(all(row["overall_status"] == UNAVAILABLE for row in audit["slope_fits"].values()))

    def test_partial_preflight_does_not_hide_later_corrupt_rapidity(self):
        for state in ("2S", "3S"):
            mark_unavailable(self.cfg, state)
        (self.base / "3S/y/result.json").unlink()
        with self.assertRaises(FileNotFoundError):
            preflight_reoptimization(self.cfg, "differential")


if __name__ == "__main__":
    unittest.main()
