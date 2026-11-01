import unittest
import numpy as np
import pandas as pd
from upsilon_esr.common import expected_intervals, project_cross_section, source_luminosity


class StatisticalProjectionTests(unittest.TestCase):
    def fixture(self):
        return pd.DataFrame({"yield_lumi_invfb": [10., 10., 10.],
            "sigma_gamma_p_pb": [20., 2., 0.],
            "n_reco_weighted_visible": [100., 1., 0.],
            "reco_count_to_sigma_pb": [.2, 2., 2.], "plotted": [True, True, True]})

    def test_zero_poisson_upper_limit(self):
        ci = expected_intervals(np.array([0.]))
        self.assertEqual(ci["lower"][0], 0.)
        self.assertAlmostEqual(ci["upper"][0], 1.841021645, places=8)

    def test_luminosity_changes_counts_not_sigma(self):
        frame = self.fixture()
        low = project_cross_section(frame, 2.5, "sigma_gamma_p_pb")
        high = project_cross_section(frame, 10., "sigma_gamma_p_pb")
        np.testing.assert_allclose(low["mean"], high["mean"])
        np.testing.assert_allclose(low.mean_count, high.mean_count / 4)
        np.testing.assert_allclose(low.projection_variance, high.projection_variance * 4)
        self.assertGreater(low.upper.iloc[1]-low.lower.iloc[1], high.upper.iloc[1]-high.lower.iloc[1])

    def test_unsupported_zero_bin_is_not_invented(self):
        frame = self.fixture().drop(columns="reco_count_to_sigma_pb")
        projected = project_cross_section(frame, 2.5, "sigma_gamma_p_pb")
        self.assertFalse(projected.valid.iloc[2])
        self.assertTrue(np.isnan(projected.upper.iloc[2]))

    def test_correction_closure_catches_factor_error(self):
        frame = self.fixture()
        frame.loc[0, "sigma_gamma_p_pb"] *= 2
        with self.assertRaisesRegex(ValueError, "does not close"):
            project_cross_section(frame, 2.5, "sigma_gamma_p_pb")

    def test_unknown_lumi_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "absent"):
            source_luminosity(pd.DataFrame({"events": [1]}))

    def test_efficiency_mask_is_retained(self):
        frame = self.fixture()
        frame.loc[0, "plotted"] = False
        self.assertFalse(project_cross_section(frame, 2.5, "sigma_gamma_p_pb").valid.iloc[0])


if __name__ == "__main__":
    unittest.main()
