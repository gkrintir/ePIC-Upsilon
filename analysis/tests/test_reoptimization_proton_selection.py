"""Synthetic regression checks for explicit proton selection, not physics validation."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "workflow"))
from reoptimization.data import STATE_DEFAULTS, load_sample
from reoptimization.response import evaluate_edges
from proton_selection import proton_selection_metadata

HAS_UPROOT = importlib.util.find_spec("uproot") is not None


def synthetic_input():
    rows = 2400
    state = dict(STATE_DEFAULTS["1S"], state="1S")
    weight = state["visible_xsec_fb"] * 10 / rows
    pt2 = np.repeat([.1, .3, .7], rows // 3)
    w = np.repeat([35., 55., 85.], rows // 3)
    y = np.repeat([-1., 0., 1.], rows // 3)
    tagged = np.arange(rows) % 2 == 0
    arrays = {
        "weight": np.full(rows, weight), "true_mass": np.full(rows, state["mass"]),
        "true_W": w, "true_Q2": np.full(rows, .1), "true_pT": np.sqrt(pt2),
        "true_t": pt2, "true_y": y,
        "reco_mass": np.full(rows, state["mass"]), "reco_W": w.copy(),
        "reco_pT": np.sqrt(pt2), "reco_y": y.copy(), "sigma_epz": np.full(rows, 18.),
        "reco_topology": np.where(tagged, 1, 3).astype(np.int32),
        "reco_tagged_proton_p": np.where(tagged, 270., -999.),
    }
    metadata = dict(entries=rows, generated_events_h_EventCount=rows,
        cross_section_visible_ee_fb=state["visible_xsec_fb"], luminosity_invfb=10.,
        generator_branching_fraction_ee=state["br_ee"], target_sum_weight=weight*rows,
        corrected_event_weight=weight, written_sum_weight=weight*rows,
        normalization_profile="jul26-purestates-estarlight-939b11a",
        normalization_scope="All generated events, before truth-state, Q2, W, and reconstruction cuts",
        normalization_source="EXPLICITLY SYNTHETIC PROTON-SELECTION TEST")
    return state, arrays, metadata


def write_root(path, arrays, metadata):
    import uproot
    with uproot.recreate(path) as out:
        tree = out.mktree("UpsilonTree", {key: value.dtype for key, value in arrays.items()})
        tree.extend(arrays)
        out["h_EventCount"] = (np.array([len(arrays["weight"])], float), np.array([0., 1.]))
        out["framework_normalization"] = json.dumps(metadata)


@unittest.skipUnless(HAS_UPROOT, "ROOT regression tests require uproot")
class ProtonSelectionInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "synthetic.root"
        self.state, self.arrays, self.meta = synthetic_input()
        self.config = dict(requested_observables=["pt2", "abs_t_via_pt2", "w", "y"])

    def tearDown(self):
        self.temp.cleanup()

    def load(self, mode="required", config=None):
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": mode}):
            return load_sample(self.path, self.state, self.config if config is None else config)

    def test_tag_is_numerator_only_for_all_four_observables(self):
        write_root(self.path, self.arrays, self.meta)
        file_before = self.path.read_bytes()
        inclusive, tagged = self.load("inclusive"), self.load()
        self.assertEqual(self.path.read_bytes(), file_before)
        for key in ("weight", "truth_selected", "is_train"):
            np.testing.assert_array_equal(tagged[key], inclusive[key])
        self.assertEqual(tagged["truth_selected"].sum(), 2400)
        self.assertEqual(tagged["reco_selected"].sum(), 1200)
        self.assertEqual(tagged["reco_selected_y"].sum(), 1200)
        self.assertAlmostEqual(tagged["weight"].sum(), self.state["visible_xsec_fb"]*2.5)
        audit = tagged["metadata"]["proton_selection_cutflow"]
        self.assertEqual(audit["truth_fiducial_reco_before_rows"], 2400)
        self.assertEqual(audit["truth_fiducial_reco_after_rows"], 1200)
        self.assertAlmostEqual(audit["truth_fiducial_reco_after_target_yield"],
                               audit["truth_fiducial_reco_before_target_yield"] / 2)
        self.assertEqual(tagged["metadata"]["proton_selection"]["mode"], "required")
        for observable, edges in (("pt2", [0,.2,.5,1]), ("abs_t_via_pt2", [0,.2,.5,1]),
                                  ("w", [20,45,70,100]), ("y", [-2,-.5,.5,2])):
            cfg = {"fit_range": [edges[0], edges[-1]], "reco_control_max": None}
            for split in ("train", "test"):
                with self.subTest(observable=observable, split=split):
                    before = evaluate_edges(inclusive, observable, edges, cfg, split=split)
                    after = evaluate_edges(tagged, observable, edges, cfg, split=split)
                    np.testing.assert_array_equal(before["truth_expected"], after["truth_expected"])
                    self.assertTrue(np.all(after["reco_expected"] < before["reco_expected"]))
                    self.assertTrue(np.all(after["metrics"]["efficiency_selected"] < 1))
                    self.assertTrue(np.all(after["raw_categories"][-1] > 0))
                    np.testing.assert_allclose(after["response"] @ after["truth_expected"] + after["feed_in"], after["reco_expected"])
                    self.assertEqual(after["metadata"]["proton_selection"]["mode"], "required")

    def test_missing_tag_branches_fail_required_but_allow_explicit_inclusive(self):
        for branch in ("reco_topology", "reco_tagged_proton_p"):
            with self.subTest(branch=branch):
                arrays = copy.deepcopy(self.arrays)
                arrays.pop(branch)
                write_root(self.path, arrays, self.meta)
                with self.assertRaisesRegex(ValueError, "missing required branches"):
                    self.load()
                inclusive = self.load("inclusive")
                self.assertEqual(inclusive["reco_selected"].sum(), 2400)
                self.assertEqual(inclusive["metadata"]["proton_selection"]["mode"], "inclusive")

    def test_inconsistent_tag_topology_fails(self):
        for topology, momentum in ((1, -999.), (3, 270.), (7, 270.), (2, np.nan)):
            with self.subTest(topology=topology, momentum=momentum):
                arrays = copy.deepcopy(self.arrays)
                arrays["reco_topology"][0] = topology
                arrays["reco_tagged_proton_p"][0] = momentum
                write_root(self.path, arrays, self.meta)
                with self.assertRaises(ValueError):
                    self.load()

    def test_policy_mismatch_fails_instead_of_relabeling(self):
        write_root(self.path, self.arrays, self.meta)
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
            policy = proton_selection_metadata()
        with self.assertRaisesRegex(ValueError, "Configured proton selection differs"):
            self.load(config=dict(self.config, proton_selection=policy))


class ProtonSelectionProvenanceTests(unittest.TestCase):
    def runner(self):
        spec = importlib.util.spec_from_file_location("proton_test_runner", ROOT / "workflow/run_reoptimization.py")
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        return runner

    def test_policy_and_helper_are_in_resume_fingerprint(self):
        runner = self.runner()
        fingerprints = {}
        for mode in ("required", "inclusive"):
            with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": mode}):
                args = runner.parse_args(["plan", "--states", "1S"])
                config, inputs = runner.prepare_config(args)
                self.assertEqual(config["proton_selection"]["mode"], mode)
                fingerprint = runner.source_fingerprint(config, inputs)
                self.assertIn("proton_selection.py", fingerprint["code_sha256"])
                fingerprints[mode] = runner.digest(fingerprint)
        self.assertNotEqual(fingerprints["required"], fingerprints["inclusive"])

    def test_report_csvs_keep_policy_and_cutflow(self):
        import csv
        from reoptimization.reporting import write_report
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            policy = proton_selection_metadata()
        cutflow = {"truth_fiducial_reco_before_rows": 10, "truth_fiducial_reco_after_rows": 5}
        entry = dict(status="NOT_VALIDATED", overall_status="NOT_VALIDATED",
            proton_selection=policy,
            selection_metadata={"proton_selection": policy, "proton_selection_cutflow": cutflow},
            state_config=dict(state="1S", q2_max=1., visible_xsec_fb=210.136, br_ee=.0238, sigma=.053),
            optimization=dict(status="SELECTED_TRAINING_ONLY", selected_edges=[0., .2, .5, 1.],
                              candidates=[], thresholds={}))
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            write_report(out, {"1S": {"pt2": entry}}, {"target_lumi_invfb": 2.5, "synthetic": True})
            for name in ("validation_status", "selected_bins", "selection_settings"):
                with (out / (name + ".csv")).open() as stream:
                    rows = list(csv.DictReader(stream))
                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0]["proton_selection_mode"], "required")
                self.assertEqual(json.loads(rows[0]["proton_selection"]), policy)
                self.assertEqual(rows[0]["overall_status"], "NOT_VALIDATED")
                if name == "selection_settings":
                    self.assertEqual(json.loads(rows[0]["proton_selection_cutflow"]), cutflow)


if __name__ == "__main__":
    unittest.main()
