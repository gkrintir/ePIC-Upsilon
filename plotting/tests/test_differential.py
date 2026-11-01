"""Numerical checks on the finite-bin, expected-count slope forecast."""
import unittest
import numpy as np
import pandas as pd
from upsilon_esr.differential import exponential_bin_average, fit_projected_exponential


class DifferentialTests(unittest.TestCase):
    @staticmethod
    def sample(lumi=2.5):
        edges = np.array([0., .08, .2, .4, .7, 1., 5.])
        means = exponential_bin_average(edges[:-1], edges[1:], 3.2, 150.)
        factor = np.array([3., 2., 1.8, 2.5, 3., 4.]) * 2.5 / lumi
        return pd.DataFrame(dict(x_low=edges[:-1], x_high=edges[1:],
                                 mean_count=means/factor, mean=means,
                                 count_to_cross_section=factor, valid=True))

    def test_fit_exact_model_and_luminosity(self):
        low, _, _ = fit_projected_exponential(self.sample(2.5))
        high, _, _ = fit_projected_exponential(self.sample(10.))
        self.assertAlmostEqual(low["b_slope_GeVminus2"], 3.2, places=5)
        self.assertAlmostEqual(high["b_slope_GeVminus2"], 3.2, places=5)
        self.assertGreater(low["b_slope_err_low"], high["b_slope_err_low"])
        self.assertGreater(low["b_slope_err_high"], high["b_slope_err_high"])
        self.assertAlmostEqual(low["amplitude_pb_per_GeV2"], 150., places=3)
        self.assertEqual(low["included_bin_indices"], [0, 1, 2, 3, 4])

    def test_wide_bin_average_is_not_center_value(self):
        average = exponential_bin_average([1.], [5.], 3.2, 150.)[0]
        analytic = 150. * (np.exp(-3.2) - np.exp(-16.)) / (3.2*4.)
        self.assertAlmostEqual(average, analytic)
        self.assertGreater(average, 150. * np.exp(-3.2*3.))

    def test_supported_zero_count_is_kept(self):
        data = self.sample()
        data.loc[3, "mean_count"] = 0.
        result, _, _ = fit_projected_exponential(data)
        self.assertIn(3, result["included_bin_indices"])
        data.loc[3, "valid"] = False
        result, _, _ = fit_projected_exponential(data)
        self.assertNotIn(3, result["included_bin_indices"])

    def test_zero_expected_yield_rejects_fit(self):
        data = self.sample()
        data["mean_count"] = 0.
        with self.assertRaisesRegex(ValueError, "not identifiable"):
            fit_projected_exponential(data)

    def test_partial_bin_never_silently_split(self):
        frame = self.sample()
        original = frame.copy(deep=True)
        result, _, _ = fit_projected_exponential(frame, (.05, 1.))
        self.assertEqual(result["partly_overlapping_bins_excluded"], [0])
        self.assertNotIn(0, result["included_bin_indices"])
        pd.testing.assert_frame_equal(frame, original)


if __name__ == "__main__":
    unittest.main()
