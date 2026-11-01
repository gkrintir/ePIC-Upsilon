"""Synthetic ROOT tests for the external support diagnostic; no production input needed."""
import ast
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import uproot

PACKAGE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("audit_t_support", PACKAGE / "audit_t_support.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
SOURCE = Path(os.environ.get("YPSILON_TEST_FRAMEWORK", "tmp/note_update_review/production/Ypsilon-Fit-Framework-ProtonTagged")).resolve()


@contextmanager
def policy_environment(mode="required"):
    old = os.environ.get("YPSILON_PROTON_SELECTION")
    os.environ["YPSILON_PROTON_SELECTION"] = mode
    try:
        yield
    finally:
        if old is None:
            os.environ.pop("YPSILON_PROTON_SELECTION", None)
        else:
            os.environ["YPSILON_PROTON_SELECTION"] = old


def fixture(n=2000):
    truth = np.resize(np.array([.1, .3, .7, 2.]), n)
    return {"weight": np.full(n, 119.685 * 10 / n),
            "true_mass": np.full(n, 9.46030), "true_t": truth,
            "true_W": np.full(n, 40.), "true_Q2": np.full(n, .001),
            "reco_mass": np.full(n, 9.46030), "reco_W": np.full(n, 40.),
            "sigma_epz": np.full(n, 18.), "reco_pT": np.sqrt(truth),
            "reco_topology": np.ones(n), "reco_tagged_proton_p": np.full(n, 100.)}


class SupportAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "tagged"
        (self.root / "workflow").mkdir(parents=True)
        (self.root / "scripts/unfolding").mkdir(parents=True)
        for rel in ["workflow/proton_selection.py", "workflow/proton_category_isolation.py",
                    "scripts/unfolding/unfold_t_spectrum.py"]:
            shutil.copy2(SOURCE / rel, self.root / rel)
        (self.root / ".proton_category.json").write_text(json.dumps({
            "schema_version": 1, "installation_status": "COMPLETE", "mode": "required",
            "source_root": str(Path(self.temp.name) / "original"), "category_root": str(self.root)}))
        self.input = Path(self.temp.name) / "normalized.root"
        self.out = self.root / "results/generated/runs/synthetic-support/artifacts/support_audit"

    def run_audit(self, arrays, lumi=1., name="synthetic-support"):
        with uproot.recreate(self.input) as target:
            target["UpsilonTree"] = arrays
        output = self.root / "results/generated/runs" / name / "artifacts/support_audit"
        with policy_environment():
            result = audit.main(["--framework-root", str(self.root), "--input", str(self.input),
                                 "--out-dir", str(output), "--projection-lumi", str(lumi)])
        self.assertEqual(result, 0)
        return json.loads((output / "audit.json").read_text())

    def test_zero_support_attributed_to_tag_and_luminosity_independent(self):
        arrays = fixture()
        tail = arrays["true_t"] > 1.
        arrays["reco_topology"][tail] = 0.
        arrays["reco_tagged_proton_p"][tail] = -999.
        one = self.run_audit(arrays)
        two = self.run_audit(arrays, 2.5, "synthetic-support-L2p5")
        self.assertEqual(one["responses"]["legacy_tagged"]["unsupported_truth_bins"], [3])
        self.assertEqual(one["responses"]["without_proton_tag"]["unsupported_truth_bins"], [])
        self.assertIn("proton requirement removes", one["diagnoses"][0]["observed_finite_MC_reason"])
        np.testing.assert_array_equal(one["responses"]["legacy_tagged"]["A_truth_rows_reco_columns"],
                                      two["responses"]["legacy_tagged"]["A_truth_rows_reco_columns"])
        np.testing.assert_allclose(np.array(one["responses"]["legacy_tagged"]["truth_test_target_yield"]) * 2.5,
                                   two["responses"]["legacy_tagged"]["truth_test_target_yield"])
        self.assertTrue(one["normalization"]["weight_sum_matches_assigned"])
        self.assertEqual(one["input_rows"], len(arrays["weight"]))

    def test_training_hole_is_not_assigned_to_proton_cut(self):
        arrays = fixture(12)
        arrays["true_t"][:] = .1
        arrays["reco_pT"][:] = np.sqrt(.1)
        test_indices = np.flatnonzero(np.random.default_rng(12345).random(12) >= .5)
        arrays["true_t"][test_indices[0]] = 2.
        arrays["reco_pT"][test_indices[0]] = np.sqrt(2.)
        result = self.run_audit(arrays)
        self.assertEqual(result["responses"]["legacy_tagged"]["unsupported_truth_bins"], [3])
        self.assertIn("No positive-weight training truth", result["diagnoses"][0]["observed_finite_MC_reason"])

    def test_nonproton_cut_hole_and_global_finite_filter(self):
        arrays = fixture()
        arrays["reco_mass"][arrays["true_t"] > 1.] = -999.
        arrays["reco_W"][0] = np.nan
        result = self.run_audit(arrays)
        self.assertEqual(result["rows_dropped_by_legacy_global_finite_filter"], 1)
        self.assertEqual(result["truth_fiducial_rows_dropped_by_legacy_global_finite_filter"], 1)
        self.assertIn("non-proton reconstruction cuts", result["diagnoses"][0]["observed_finite_MC_reason"])

    def test_negative_pt_leakage_and_beyond_five_feed_in(self):
        arrays = fixture()
        arrays["reco_pT"][0] = -.3  # Legacy square passes; physical selection rejects.
        arrays["reco_pT"][1] = -999.  # Legacy range already rejects this sentinel.
        arrays["true_t"][-1] = 6.
        arrays["reco_pT"][-1] = .5
        result = self.run_audit(arrays)
        self.assertEqual(result["negative_reco_pT_rows_passing_legacy_final"], 1)
        self.assertEqual(result["sentinel_reco_pT_minus999_rows"], 1)
        rows = [row for row in result["tail_feed_in"] if row["region"] == "truth_above_5"
                and row["selection"] == "reco_core_with_proton_tag"]
        self.assertEqual(rows[0]["rows"], 1)

    def test_response_matches_actual_framework_including_right_edge(self):
        tree = ast.parse((SOURCE / "scripts/unfolding/unfold_t_spectrum.py").read_text())
        nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name in {"hist_weighted", "response_matrix"}]
        namespace = {"np": np}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), "framework_response", "exec"), namespace)
        truth = np.array([.1, .3, .7, 5., .1, .3, .7, 5.])
        reco = truth.copy()
        weights = np.linspace(.1, .8, 8)
        train = np.arange(8) < 4
        selected = np.ones(8, dtype=bool)
        got = audit.response(truth, reco, weights, selected, train, ~train, audit.EDGES, .2)
        expected = namespace["response_matrix"](truth[train], reco[train], weights[train], audit.EDGES, selected[train])
        np.testing.assert_array_equal(got["H_truth_rows_reco_columns"], expected[0])
        np.testing.assert_array_equal(got["A_truth_rows_reco_columns"], expected[4])
        self.assertAlmostEqual(got["efficiency"][3], 1., places=14)

    def test_refuses_overwrite_and_inclusive_policy(self):
        arrays = fixture()
        self.run_audit(arrays)
        with self.assertRaises(SystemExit):
            with policy_environment():
                audit.main(["--framework-root", str(self.root), "--input", str(self.input), "--out-dir", str(self.out)])
        with self.assertRaises(ValueError):
            with policy_environment("inclusive"):
                audit.helpers(self.root, self.out.parent / "new-audit")

    def test_refuses_missing_category_marker(self):
        (self.root / ".proton_category.json").unlink()
        with self.assertRaisesRegex(ValueError, "completed isolated proton-tagged category marker"):
            with policy_environment():
                audit.helpers(self.root, self.out)


if __name__ == "__main__":
    unittest.main()
