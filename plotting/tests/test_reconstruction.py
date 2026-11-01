"""Numerical tests; the inputs below are synthetic, not production results."""
import unittest

import numpy as np

from upsilon_esr.reconstruction import MASSES, SHAPE_NAMES, fit_combined_mass, mass_probabilities


class CombinedMassFitTests(unittest.TestCase):
    def setUp(self):
        self.edges = np.linspace(8., 12., 81)
        self.calibrations, self.means = {}, {}
        self.yields = {"1S": 105., "2S": 32., "3S": 27.}
        for i, state in enumerate(MASSES):
            shape = [MASSES[state] - .006, .06 * (1 + .05*i), 1.4, 3.5, 1.8, 4.]
            self.calibrations[state] = {
                "converged": True,
                "parameters": {name: {"value": value} for name, value in zip(SHAPE_NAMES, shape)}}
            self.means[state] = self.yields[state] * mass_probabilities(self.edges, shape)

    def test_recovers_injected_three_state_mixture(self):
        _, result = fit_combined_mass(self.edges, self.means, self.calibrations)
        self.assertTrue(result["converged"])
        self.assertLess(result["poisson_deviance"], 1e-8)
        for state in MASSES:
            self.assertAlmostEqual(result["states"][state]["fitted_yield"], self.yields[state], places=4)
        self.assertAlmostEqual(result["common_mass_offset_GeV"], -.006, places=6)
        self.assertAlmostEqual(result["common_width_scale"], 1., places=6)

    def test_luminosity_changes_yields_not_core_widths(self):
        _, first = fit_combined_mass(self.edges, self.means, self.calibrations)
        _, second = fit_combined_mass(self.edges, {s: 4*x for s, x in self.means.items()}, self.calibrations)
        for state in MASSES:
            self.assertAlmostEqual(second["states"][state]["fitted_yield"] / first["states"][state]["fitted_yield"], 4., places=5)
            self.assertAlmostEqual(second["states"][state]["fitted_sigma_core_GeV"], first["states"][state]["fitted_sigma_core_GeV"], places=6)

    def test_exact_bin_integrals_survive_rebinning(self):
        fine = np.linspace(8., 12., 161)
        for state in MASSES:
            shape = [self.calibrations[state]["parameters"][name]["value"] for name in SHAPE_NAMES]
            fine_prob = mass_probabilities(fine, shape).reshape(-1, 2).sum(axis=1)
            np.testing.assert_allclose(fine_prob, mass_probabilities(self.edges, shape), rtol=1e-9, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
