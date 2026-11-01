"""Projected mass-parameter interval checks; all inputs are synthetic."""
import copy
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np
from scipy.optimize import brentq

from upsilon_esr.reconstruction import MASSES, SHAPE_NAMES, fit_combined_mass, mass_probabilities
from upsilon_esr.mass_fit_intervals import compute_mass_fit_intervals


def fixture(scale=1., width_scale=1., floated=True):
    edges = np.linspace(8., 12., 81)
    calibrations, means = {}, {}
    for i, state in enumerate(MASSES):
        shape = np.array([MASSES[state] - .006, .06 * (1 + .05*i), 1.4, 3.5, 1.8, 4.])
        calibrations[state] = {"converged": True, "parameters": {
            name: {"value": value} for name, value in zip(SHAPE_NAMES, shape)}}
        shape[1] *= width_scale
        means[state] = scale * [105., 32., 27.][i] * mass_probabilities(edges, shape)
    _, summary = fit_combined_mass(edges, means, calibrations, float_common_shape=floated)
    return edges, means, calibrations, summary


class MassIntervalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = fixture()
        cls.result = compute_mass_fit_intervals(*cls.inputs, luminosity_invfb=2.5)

    def test_nominal_fit_and_inputs_unchanged(self):
        edges, means, calibrations, summary = fixture()
        saved = copy.deepcopy((edges, means, calibrations, summary))
        result = compute_mass_fit_intervals(edges, means, calibrations, summary, luminosity_invfb=2.5)
        np.testing.assert_array_equal(edges, saved[0])
        for state in MASSES:
            np.testing.assert_array_equal(means[state], saved[1][state])
            self.assertEqual(result["states"][state]["fitted_yield"]["value"], summary["states"][state]["fitted_yield"])
            self.assertEqual(result["states"][state]["fitted_sigma_core_MeV"]["value"],
                             summary["states"][state]["fitted_sigma_core_GeV"] * 1000.)
        self.assertEqual(calibrations, saved[2])
        self.assertEqual(summary, saved[3])
        self.assertFalse(result["coverage_validated"])
        self.assertFalse(result["publication_ready"])
        self.assertAlmostEqual(result["delta_poisson_deviance_threshold"], 1., places=12)

    def test_profile_agrees_with_independent_fisher_locally(self):
        # Independent expected Fisher benchmark from finite differences of
        # exact model integrals (not least_squares' residual Jacobian).
        edges, means, calibrations, summary = self.inputs
        nominal = np.r_[[np.log(summary["states"][s]["fitted_yield"]) for s in MASSES],
                         summary["common_mass_offset_GeV"], np.log(summary["common_width_scale"])]
        def prediction(p):
            values = []
            for i, state in enumerate(MASSES):
                shape = np.array([calibrations[state]["parameters"][n]["value"] for n in SHAPE_NAMES])
                shape[0] = MASSES[state] + p[3]
                shape[1] *= np.exp(p[4])
                values.append(np.exp(p[i]) * mass_probabilities(edges, shape))
            return np.sum(values, axis=0)
        mean = prediction(nominal)
        jac = np.stack([(prediction(nominal + np.eye(5)[i]*1e-5) -
                         prediction(nominal - np.eye(5)[i]*1e-5))/(2e-5) for i in range(5)], axis=1)
        covariance = np.linalg.inv(jac.T @ (jac/mean[:, None]))
        for i, state in enumerate(MASSES):
            r = self.result["states"][state]["fitted_yield"]
            self.assertEqual(r["status"], "TWO_SIDED")
            symmetric = .5 * (r["error_low"] + r["error_high"])
            fisher = r["value"] * np.sqrt(covariance[i, i])
            self.assertLess(abs(symmetric/fisher - 1), .02)
            self.assertGreater(r["error_high"], r["error_low"])

    def test_shared_resolution_scale_and_luminosity_scaling(self):
        brighter = compute_mass_fit_intervals(*fixture(scale=4.), luminosity_invfb=10.)
        for state in MASSES:
            r = self.result["states"][state]
            b = brighter["states"][state]
            self.assertAlmostEqual(b["fitted_yield"]["value"]/r["fitted_yield"]["value"], 4., places=5)
            self.assertAlmostEqual(b["fitted_sigma_core_MeV"]["value"], r["fitted_sigma_core_MeV"]["value"], places=5)
            ey = lambda obj: .5*(obj["error_low"] + obj["error_high"])
            self.assertLess(abs(ey(b["fitted_yield"])/ey(r["fitted_yield"]) - 2.), .03)
            self.assertLess(abs(ey(b["fitted_sigma_core_MeV"])/ey(r["fitted_sigma_core_MeV"]) - .5), .01)
            for endpoint in ("lower", "upper"):
                self.assertAlmostEqual(r["fitted_sigma_core_MeV"][endpoint]/r["fitted_sigma_core_MeV"]["value"],
                                       self.result["common_width_scale"][endpoint]/self.result["common_width_scale"]["value"], places=12)

    def test_disjoint_templates_match_exact_poisson_likelihood_roots(self):
        edges = np.arange(4.)
        yields = [8., 3., 1.2]
        calibrations, means = {}, {}
        for i, state in enumerate(MASSES):
            shape = [MASSES[state], .06, 1.4, 3.5, 1.8, 4.]
            calibrations[state] = {"converged": True, "parameters": {
                name: {"value": value} for name, value in zip(SHAPE_NAMES, shape)}}
            means[state] = np.eye(3)[i] * yields[i]
        def separated(edges, shape):
            closest = int(np.argmin([abs(shape[0] - m) for m in MASSES.values()]))
            return np.eye(3)[closest]
        with mock.patch("upsilon_esr.reconstruction.mass_probabilities", side_effect=separated), \
             mock.patch("upsilon_esr.mass_fit_intervals.mass_probabilities", side_effect=separated):
            _, summary = fit_combined_mass(edges, means, calibrations, float_common_shape=False)
            result = compute_mass_fit_intervals(edges, means, calibrations, summary, luminosity_invfb=2.5)
        for state, observed in zip(MASSES, yields):
            delta = lambda mean: 2*(mean-observed+observed*np.log(observed/mean))-1
            low = brentq(delta, 1e-8, observed)
            high = brentq(delta, observed, 20*observed)
            interval = result["states"][state]["fitted_yield"]
            self.assertAlmostEqual(interval["lower"], low, places=6)
            self.assertAlmostEqual(interval["upper"], high, places=6)
            self.assertEqual(result["states"][state]["fitted_sigma_core_MeV"]["status"], "FIXED_PARAMETER")
            self.assertIsNone(result["states"][state]["fitted_sigma_core_MeV"]["error_low"])

    def test_bound_limited_resolution_has_no_displayable_errors(self):
        inputs = fixture(width_scale=2.)
        result = compute_mass_fit_intervals(*inputs, luminosity_invfb=2.5)
        for state in MASSES:
            interval = result["states"][state]["fitted_sigma_core_MeV"]
            self.assertIn("UPPER_BOUND_REACHED", interval["status"])
            self.assertIsNone(interval["error_low"])
            self.assertIsNone(interval["error_high"])
            self.assertIsNone(interval["upper"])

    def test_failed_optimization_cannot_supply_errors(self):
        fake = lambda *args, **kw: SimpleNamespace(success=False, x=np.asarray(args[1]))
        with mock.patch("upsilon_esr.mass_fit_intervals.least_squares", side_effect=fake):
            with self.assertRaisesRegex(RuntimeError, "All nuisance-fit starts failed"):
                compute_mass_fit_intervals(*self.inputs, luminosity_invfb=2.5)

    def test_mismatched_nominal_deviance_or_luminosity_is_rejected(self):
        edges, means, cals, summary = copy.deepcopy(self.inputs)
        summary["poisson_deviance"] = 100.
        with self.assertRaisesRegex(ValueError, "Nominal deviance differs"):
            compute_mass_fit_intervals(edges, means, cals, summary, luminosity_invfb=2.5)
        summary["target_luminosity_invfb"] = 10.
        with self.assertRaisesRegex(ValueError, "Target luminosity differs"):
            compute_mass_fit_intervals(edges, means, cals, summary, luminosity_invfb=2.5)

    def test_nonstationary_nominal_cannot_be_kept_and_given_errors(self):
        edges, means, cals, summary = copy.deepcopy(self.inputs)
        from upsilon_esr.reconstruction import poisson_residual
        summary["states"]["1S"]["fitted_yield"] *= 1.1
        model = np.zeros(len(edges)-1)
        for state in MASSES:
            shape = np.array([cals[state]["parameters"][n]["value"] for n in SHAPE_NAMES])
            shape[0] = summary["states"][state]["fitted_mean_GeV"]
            shape[1] = summary["states"][state]["fitted_sigma_core_GeV"]
            model += summary["states"][state]["fitted_yield"] * mass_probabilities(edges, shape)
        summary["poisson_deviance"] = float(np.sum(poisson_residual(np.sum(list(means.values()), axis=0), model)**2))
        with self.assertRaisesRegex(RuntimeError, "better mass-fit minimum"):
            compute_mass_fit_intervals(edges, means, cals, summary, luminosity_invfb=2.5)


if __name__ == "__main__":
    unittest.main()
