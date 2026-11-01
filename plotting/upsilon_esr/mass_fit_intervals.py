"""Conditional projected profile-likelihood intervals for the saved joint mass fit.

This additive module never changes the nominal fit or its displayed curves.
The interval threshold is asymptotic: a nominal 68.27% interval is not a claim
of toy-validated coverage for the sparse, possibly misspecified mass model.
"""
import copy

import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.stats import chi2

from .reconstruction import MASSES, SHAPE_NAMES, mass_probabilities, poisson_residual


def _deviance(observed, model):
    return float(np.sum(poisson_residual(observed, model) ** 2))


def _interval(value, low, high, low_status, high_status):
    good = low_status == high_status == "CROSSING"
    status = "TWO_SIDED" if good else ";".join(
        direction + "_" + flag for direction, flag in
        (("LOWER", low_status), ("UPPER", high_status)) if flag != "CROSSING")
    return {"value": float(value), "lower": None if low is None else float(low),
            "upper": None if high is None else float(high),
            "error_low": float(value - low) if good else None,
            "error_high": float(high - value) if good else None,
            "status": status}


class _Profile:
    """Bounded fit with exact original bin integrals and nuisance parameters."""
    def __init__(self, edges, expected, shapes, initial, lower, upper, active):
        self.edges, self.expected, self.shapes = edges, expected, shapes
        self.nominal = initial.copy()
        self.lower, self.upper, self.active = lower, upper, active
        self.cache = {}
        self.n_profile_fits = 0
        self.n_failed_starts = 0
        self.baseline = _deviance(expected, self.model(initial))
        self.allowed_improvement = max(2e-6, 2e-7 * max(1., self.baseline))

    def model(self, full):
        values = []
        for j, state in enumerate(MASSES):
            shape = self.shapes[state].copy()
            shape[0] = MASSES[state] + full[3]
            shape[1] *= np.exp(full[4])
            values.append(np.exp(full[j]) * mass_probabilities(self.edges, shape))
        return np.sum(values, axis=0)

    def optimize(self, fixed_index=None, fixed_value=None):
        free = [i for i in self.active if i != fixed_index]
        seeds = [self.nominal.copy()]
        # A second start changes both common-shape coordinates. Yields are
        # also initialized at the observed total so nuisance profiling does
        # not inherit a poor normalization from a distant scanned value.
        second = self.nominal.copy()
        if 3 in free:
            second[3] = np.clip(second[3] + .006, self.lower[3], self.upper[3])
        if 4 in free:
            second[4] = np.clip(second[4] - .08, self.lower[4], self.upper[4])
        for i in free:
            if i < 3:
                second[i] = np.clip(second[i] + .05, self.lower[i], self.upper[i])
        seeds.append(second)
        if fixed_index is not None:
            nearby = [(abs(x - fixed_value), p) for (j, x), (_, p) in self.cache.items()
                      if j == fixed_index]
            if nearby:
                seeds.append(min(nearby, key=lambda item: item[0])[1].copy())
        candidates = []
        for seed in seeds:
            if fixed_index is not None:
                seed[fixed_index] = fixed_value

            def expand(p):
                full = seed.copy()
                full[free] = p
                return full

            try:
                result = least_squares(
                    lambda p: poisson_residual(self.expected, self.model(expand(p))),
                    seed[free], bounds=(self.lower[free], self.upper[free]),
                    x_scale="jac", ftol=1e-10, xtol=1e-10, gtol=1e-9,
                    max_nfev=5000)
                full = expand(result.x)
                value = _deviance(self.expected, self.model(full))
                if result.success and np.isfinite(value) and np.all(np.isfinite(full)):
                    candidates.append((value, full))
                else:
                    self.n_failed_starts += 1
            except (ValueError, FloatingPointError, RuntimeError):
                self.n_failed_starts += 1
        self.n_profile_fits += 1
        if not candidates:
            raise RuntimeError("All nuisance-fit starts failed; no mass-fit interval may be quoted.")
        best = min(candidates, key=lambda item: item[0])
        if best[0] < self.baseline - self.allowed_improvement:
            raise RuntimeError("Profiling found a better mass-fit minimum; nominal fit must be reviewed before quoting errors.")
        return best

    def profile(self, index, value):
        key = (index, float(value))
        if key not in self.cache:
            self.cache[key] = self.optimize(index, float(value))
        return self.cache[key][0]

    def endpoints(self, index, threshold):
        nominal = float(self.nominal[index])
        result = []
        for sign, boundary in ((-1., self.lower[index]), (1., self.upper[index])):
            if abs(nominal - boundary) < 1e-7:
                result.append((None, "BOUND_REACHED"))
                continue
            # Geometric stepping brackets both broad low-count and narrow
            # high-count intervals without a Hessian approximation.
            step = .002 if index == 3 else .025
            previous = nominal
            found = False
            while True:
                trial = float(np.clip(nominal + sign * step, self.lower[index], self.upper[index]))
                delta = self.profile(index, trial) - self.baseline - threshold
                if delta >= 0:
                    left, right = sorted((previous, trial))
                    root, report = brentq(
                        lambda x: self.profile(index, x) - self.baseline - threshold,
                        left, right, xtol=2e-8, rtol=1e-9,
                        maxiter=100, full_output=True)
                    residual = self.profile(index, root) - self.baseline - threshold
                    if not report.converged or abs(residual) > max(2e-5, threshold * 2e-5):
                        raise RuntimeError("Mass-profile crossing did not converge; no interval may be quoted.")
                    result.append((float(root), "CROSSING"))
                    found = True
                    break
                if trial == boundary:
                    break
                previous = trial
                step *= 1.8
            if not found:
                result.append((None, "BOUND_REACHED"))
        return result


def compute_mass_fit_intervals(edges, means, calibrations, fit_summary, *,
                               luminosity_invfb, coverage=.6826894921370859):
    """Return projected asymmetric intervals without mutating the nominal fit.

    Inputs must be the same rebinned, luminosity-scaled per-state histograms
    used for ``fit_combined_mass``. ``means`` are expected event counts, not
    MC event counts or sumw2. The table can quote errors only when ``status``
    is ``TWO_SIDED``; bound-limited cases deliberately have null errors.
    """
    if not np.isfinite(luminosity_invfb) or luminosity_invfb <= 0:
        raise ValueError("A positive target luminosity is required.")
    if not 0 < coverage < 1:
        raise ValueError("coverage must lie strictly between zero and one.")
    if not fit_summary.get("converged", False):
        raise ValueError("The nominal simultaneous fit did not converge.")
    if "target_luminosity_invfb" in fit_summary and not np.isclose(
            fit_summary["target_luminosity_invfb"], luminosity_invfb, rtol=1e-10, atol=1e-12):
        raise ValueError("Target luminosity differs from the nominal fit.")
    edges = np.asarray(edges, float)
    if edges.ndim != 1 or len(edges) < 3 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError("Mass edges must be finite and strictly increasing.")
    counts, shapes = [], {}
    for state in MASSES:
        y = np.asarray(means[state], float)
        if y.shape != (len(edges)-1,) or np.any(~np.isfinite(y)) or np.any(y < 0) or y.sum() <= 0:
            raise ValueError(f"{state}: invalid projected mass histogram.")
        if not calibrations[state].get("converged", False):
            raise ValueError(f"{state}: upstream shape fit did not converge.")
        shapes[state] = np.array([calibrations[state]["parameters"][name]["value"]
                                  for name in SHAPE_NAMES], float)
        mass_probabilities(edges, shapes[state])
        if not np.isclose(y.sum(), fit_summary["states"][state]["expected_selected_yield"], rtol=1e-9, atol=1e-10):
            raise ValueError(f"{state}: expected yield differs from the nominal fit.")
        counts.append(y.copy())
    totals = np.array([y.sum() for y in counts])
    expected = np.sum(counts, axis=0)
    full = np.r_[[np.log(fit_summary["states"][s]["fitted_yield"]) for s in MASSES],
                 fit_summary["common_mass_offset_GeV"], np.log(fit_summary["common_width_scale"])]
    lower = np.r_[np.full(3, -20.), -.1, np.log(.5)]
    upper = np.r_[np.log(np.maximum(10 * totals, 100.)), .1, np.log(2.)]
    active = np.arange(5 if fit_summary["common_shape_floating"] else 3)
    if not np.all(np.isfinite(full)) or np.any(full < lower - 1e-10) or np.any(full > upper + 1e-10):
        raise ValueError("Nominal fit parameters lie outside the original fit bounds.")
    full = np.clip(full, lower, upper)
    for state in MASSES:
        nominal_state = fit_summary["states"][state]
        if not np.isclose(nominal_state["fitted_mean_GeV"], MASSES[state] + full[3], atol=1e-10, rtol=0):
            raise ValueError("Nominal fit does not use the fixed mass-difference convention.")
        if not np.isclose(nominal_state["fitted_sigma_core_GeV"], shapes[state][1] * np.exp(full[4]), rtol=1e-10, atol=1e-12):
            raise ValueError("Nominal fit does not use the shared core-width scale.")
    profile = _Profile(edges, expected, shapes, full, lower, upper, active)
    if not np.isclose(profile.baseline, fit_summary["poisson_deviance"], rtol=2e-7, atol=2e-6):
        raise ValueError("Nominal deviance differs from the supplied mass histogram/model.")
    profile.optimize()  # Validate the unchanged nominal optimum before profiling.
    threshold = float(chi2.ppf(coverage, 1))
    output = {"method": "conditional_asymptotic_profile_likelihood",
              "target_luminosity_invfb": float(luminosity_invfb),
              "nominal_coverage": float(coverage), "delta_poisson_deviance_threshold": threshold,
              "coverage_validated": False, "publication_ready": False,
              "interpretation": "Expected conditional statistical profile widths from projected histogram means; nominal asymptotic coverage, not an interval for observed data or validated low-count coverage.",
              "assumptions": ["Independent signal-only Poisson mass-bin counts at target luminosity",
                              "Three yields and (if enabled) common mass offset and common width scale profiled jointly",
                              "Fixed state mass differences, DSCB tails and relative widths from upstream MC calibration",
                              "No background, calibration-MC, luminosity, branching-fraction or other systematic uncertainty",
                              "Known template/DSCB shape mismatch is not covered by these statistical intervals"],
              "resolution_correlation": "All core widths share one scale: fractional width intervals are identical and their width-scale uncertainty is fully correlated.",
              "upstream_shape_mismatch_flags": {s: bool(calibrations[s].get("shape_mismatch_flag", False)) for s in MASSES},
              "states": {}, "common_width_scale": None}
    for index, state in enumerate(MASSES):
        endpoints = profile.endpoints(index, threshold)
        vals = [None if x is None else np.exp(x) for x, _ in endpoints]
        output["states"][state] = {"fitted_yield": _interval(
            fit_summary["states"][state]["fitted_yield"], vals[0], vals[1], endpoints[0][1], endpoints[1][1])}
    if 4 in active:
        endpoints = profile.endpoints(4, threshold)
        vals = [None if x is None else np.exp(x) for x, _ in endpoints]
        width = _interval(fit_summary["common_width_scale"], vals[0], vals[1], endpoints[0][1], endpoints[1][1])
    else:
        width = {"value": float(fit_summary["common_width_scale"]), "lower": None, "upper": None,
                 "error_low": None, "error_high": None, "status": "FIXED_PARAMETER"}
    output["common_width_scale"] = width
    for state in MASSES:
        resolution = copy.deepcopy(width)
        factor = float(shapes[state][1] * 1000.)
        for key in ("value", "lower", "upper", "error_low", "error_high"):
            if resolution[key] is not None:
                resolution[key] *= factor
        # Preserve the precise nominal central value, including its floating
        # point multiplication order, rather than reconstructing that value.
        resolution["value"] = float(fit_summary["states"][state]["fitted_sigma_core_GeV"] * 1000.)
        if resolution["status"] == "TWO_SIDED":
            resolution["error_low"] = resolution["value"] - resolution["lower"]
            resolution["error_high"] = resolution["upper"] - resolution["value"]
        output["states"][state]["fitted_sigma_core_MeV"] = resolution
    output["diagnostics"] = {"nominal_deviance": profile.baseline,
                             "profile_fit_evaluations": profile.n_profile_fits,
                             "failed_optimizer_starts": profile.n_failed_starts,
                             "nominal_parameter_on_bound": bool(np.any(np.isclose(full[active], lower[active], atol=1e-5, rtol=0) |
                                                                         np.isclose(full[active], upper[active], atol=1e-5, rtol=0))),
                             "profile_solutions_with_nuisance_on_bound": sum(
                                 bool(np.any(np.isclose(p[[i for i in active if i != j]], lower[[i for i in active if i != j]], atol=1e-5, rtol=0) |
                                             np.isclose(p[[i for i in active if i != j]], upper[[i for i in active if i != j]], atol=1e-5, rtol=0)))
                                 for (j, _), (_, p) in profile.cache.items()),
                             "nominal_values_preserved": True}
    return output
