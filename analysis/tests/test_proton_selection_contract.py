import os
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches


class ProtonSelectionContract(unittest.TestCase):
    def test_mask_does_not_change_full_generated_denominator(self):
        a = {"reco_topology": np.array([0, 1, 2, 3, 4]),
             "reco_tagged_proton_p": np.array([-999, 250, 260, -999, -999.])}
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            m = proton_tag_mask(a)
        np.testing.assert_array_equal(m, [False, True, True, False, False])
        self.assertEqual(len(a["reco_topology"]), 5)
        self.assertEqual(m.mean(), 0.4)

    def test_missing_or_inconsistent_schema_fails(self):
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            with self.assertRaises(ValueError):
                proton_tag_mask({"weight": np.ones(2)})
            for t, p in [(1, -999), (2, float("nan")), (3, 200), (0, 250), (5, -999), (1.5, 200)]:
                with self.subTest(t=t,p=p), self.assertRaises(ValueError):
                    proton_tag_mask({"reco_topology": [t], "reco_tagged_proton_p": [p]})

    def test_explicit_inclusive_category(self):
        with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
            self.assertEqual(proton_branches(), [])
            np.testing.assert_array_equal(proton_tag_mask({"weight": np.ones(3)}), [True]*3)


if __name__ == "__main__":
    unittest.main()
