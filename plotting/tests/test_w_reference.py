"""Small regression checks for published W points and normalization auditing."""
import unittest
import numpy as np
import pandas as pd

from upsilon_esr.w_cross_section import zeus_reference, _check_normalization


class WReferenceTests(unittest.TestCase):
    def test_published_reference_energies_and_asymmetric_totals(self):
        data = zeus_reference()
        np.testing.assert_array_equal(data.W_ref, [100., 180.])
        np.testing.assert_array_equal(data.sigma_pb, [160., 321.])
        np.testing.assert_allclose(data.err_low_pb, np.hypot([51., 88.], [21., 114.]))
        np.testing.assert_allclose(data.err_high_pb, np.hypot([51., 88.], [48., 46.]))
        self.assertFalse(np.allclose(data.err_low_pb, data.err_high_pb))

    def test_truth_cross_section_audit_rejects_extra_branching_factor(self):
        br = .0238
        frame = pd.DataFrame({"flux_int": [.01, .02],
                              "n_true_weighted_visible": [25., 75.],
                              "yield_lumi_invfb": [10., 10.]})
        frame["sigma_gamma_p_truth_closure_pb"] = (
            frame.n_true_weighted_visible / (1000 * 10 * br * frame.flux_int))
        normalization = pd.DataFrame([{"br_ee": br}])
        result = _check_normalization(frame, normalization, "1S")
        self.assertEqual(result["status"], "passed")
        frame["sigma_gamma_p_truth_closure_pb"] /= br
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            _check_normalization(frame, normalization, "1S")


if __name__ == "__main__":
    unittest.main()
