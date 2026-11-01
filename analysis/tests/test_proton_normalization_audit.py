"""Independent W audit: proton selection changes reco only, never normalization."""
import importlib.util
import os
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "proton_normalization_physics", ROOT / "tools/normalization_audit/normalization_physics.py")
physics = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(physics)


def settings():
    return dict(state="1S", electron_energy=9., proton_energy=275., q2_cut=1., q2_max=1.,
                w_edges=[20., 40., 60.], lumi_invfb=10., br_ee=.0238, mass=9.46030,
                mass_sigma=.053, mass_window_sigma=3., true_mass_window=.2,
                sigma_epz_center=18., sigma_epz_halfwidth=1.8, min_eff=0.)


def events():
    # Every row is eligible truth; the first row has no valid reconstruction.
    # Last row is tagged but fails EPZ, to distinguish tagging from final selection.
    return dict(weight=np.array([1., 2., 3., 4., 5., 6.]),
                true_mass=np.full(6, 9.46030), true_Q2=np.full(6, .001),
                true_W=np.array([30., 30., 30., 50., 50., 50.]),
                reco_mass=np.array([-999., 9.46030, 9.46030, 9.46030, 9.46030, 9.46030]),
                reco_W=np.array([-999., 30., 30., 50., 50., 50.]),
                sigma_epz=np.array([-999., 18., 18., 18., 18., 0.]),
                reco_topology=np.array([0, 1, 2, 3, 4, 1]),
                reco_tagged_proton_p=np.array([-999., 260., 250., -999., -999., 240.]))


class ProtonNormalizationAuditTests(unittest.TestCase):
    def test_tag_preserves_generated_truth_and_flux(self):
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            tagged = physics.physics_summary(events(), settings())
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
            inclusive = physics.physics_summary(events(), settings())
        self.assertEqual(tagged["cutflow"]["all"]["count"], 6)
        self.assertEqual(tagged["cutflow"]["all"]["sumw"], 21.)
        self.assertEqual(tagged["cutflow"]["reco_W_before_proton"]["sumw"], 14.)
        self.assertEqual(tagged["cutflow"]["reco_W"]["sumw"], 5.)
        self.assertEqual(inclusive["cutflow"]["reco_W"]["sumw"], 14.)
        for a, b in zip(tagged["bins"], inclusive["bins"]):
            for field in ("truth_count", "truth_sumw", "truth_sumw2", "flux_int", "independent_flux"):
                self.assertEqual(a[field], b[field], field)
        self.assertEqual(tagged["bins"][0]["reco_count"], 2)
        self.assertAlmostEqual(tagged["bins"][0]["eff_reco_over_true"], 5./6.)
        self.assertEqual(tagged["bins"][1]["reco_count"], 0)
        self.assertIsNone(tagged["bins"][1]["corrected_visible"])
        self.assertEqual(tagged["settings"]["proton_selection"], "required")

    def test_schema_required_only_in_tagged_mode(self):
        arrays = events()
        del arrays["reco_topology"]
        del arrays["reco_tagged_proton_p"]
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            with self.assertRaisesRegex(ValueError, "Missing required physics branches"):
                physics.physics_summary(arrays, settings())
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
            result = physics.physics_summary(arrays, settings())
            self.assertEqual(result["cutflow"]["all"]["count"], 6)
            self.assertEqual(result["cutflow"]["reco_W"]["count"], 4)

    def test_contract_violations_fail_before_accounting(self):
        cases = [("reco_topology", 0, np.nan), ("reco_topology", 1, 1.5),
                 ("reco_topology", 2, 5), ("reco_tagged_proton_p", 0, 260.),
                 ("reco_tagged_proton_p", 1, -999.), ("reco_tagged_proton_p", 2, np.nan),
                 ("reco_tagged_proton_p", 3, 0.)]
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            for key, index, value in cases:
                with self.subTest(key=key, index=index, value=value):
                    arrays = {k: np.asarray(v, float).copy() for k, v in events().items()}
                    arrays[key][index] = value
                    accumulator = physics.PhysicsAccumulator(settings())
                    with self.assertRaisesRegex(ValueError, "topology|Proton"):
                        accumulator.update(arrays)
                    self.assertEqual(accumulator.cutflow["all"]["count"], 0)

    def test_policy_and_required_branches_frozen_across_chunks(self):
        arrays = events()
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            accumulator = physics.PhysicsAccumulator(settings())
            accumulator.update({k: v[:3] for k, v in arrays.items()})
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
            accumulator.update({k: v[3:] for k, v in arrays.items()})
        self.assertEqual(accumulator.cutflow["all"]["count"], 6)
        self.assertEqual(accumulator.cutflow["reco_W"]["sumw"], 5.)
        self.assertIn("reco_topology", accumulator.required_branches)

    def test_metadata_dictionary_mode_is_supported(self):
        configured = dict(settings(), proton_selection={"mode": "required"})
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
            result = physics.physics_summary(events(), configured)
        self.assertEqual(result["settings"]["proton_selection"], "required")
        self.assertEqual(result["proton_selection"]["mode"], "required")
        self.assertEqual(result["cutflow"]["reco_W"]["sumw"], 5.)
        configured["proton_selection"] = {"mode": "inclusive"}
        arrays = {k: v for k, v in events().items() if not k.startswith("reco_topology") and k != "reco_tagged_proton_p"}
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            result = physics.physics_summary(arrays, configured)
        self.assertEqual(result["cutflow"]["reco_W"]["sumw"], 14.)

    def test_unknown_policy_rejected(self):
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "automatic"}):
            with self.assertRaisesRegex(ValueError, "must be required or inclusive"):
                physics.PhysicsAccumulator(settings())


if __name__ == "__main__":
    unittest.main()
