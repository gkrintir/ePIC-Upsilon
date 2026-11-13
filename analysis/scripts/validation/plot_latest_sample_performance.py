#!/usr/bin/env python3
"""Signal mass, tuple-conditional selection efficiency, and stored detector diagnostics.

The signal-only double-sided Crystal Ball fit uses unweighted MC counts and a
binned Poisson likelihood. MC calibration is retained in both statistical modes.
Projected mode additionally forecasts signal-only data precision using expected
Poisson Fisher information, independently of finite-MC template uncertainty.
Efficiency denominators comprise eligible rows in UpsilonTree, not an independently
verified number of generated or simulated events. Stored histograms retain their
producer selections and are shown separately from the tuple selections.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

# Shared explicit reconstructed-proton selection; truth denominators are unchanged.
from pathlib import Path as _ProtonPath
import sys as _proton_sys
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
import uproot
from scipy.optimize import least_squares
from scipy.special import erf
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projected_intervals import (add_interval_arguments, expected_count_intervals,
                                 draw_absolute_intervals)

MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
# Retain the historical constant for callers of the original 1S helpers.
MASS = MASSES["1S"]


def crystal_ball_integral(z, alpha_l, n_l, alpha_r, n_r):
    """Unnormalized DSCB cumulative integral from minus infinity to z."""
    z = np.asarray(z, dtype=float)
    left_total = np.exp(-0.5 * alpha_l**2) * n_l / alpha_l / (n_l - 1)
    core_total = np.sqrt(np.pi / 2) * (erf(alpha_r / np.sqrt(2)) + erf(alpha_l / np.sqrt(2)))
    right_total = np.exp(-0.5 * alpha_r**2) * n_r / alpha_r / (n_r - 1)
    result = np.empty_like(z)
    left = z < -alpha_l
    right = z > alpha_r
    core = ~(left | right)
    b_l = n_l / alpha_l - alpha_l
    b_r = n_r / alpha_r - alpha_r
    result[left] = np.exp(n_l * np.log(n_l / alpha_l) - alpha_l**2 / 2
                          - np.log(n_l - 1) + (1 - n_l) * np.log(b_l - z[left]))
    result[core] = left_total + np.sqrt(np.pi / 2) * (erf(z[core] / np.sqrt(2)) + erf(alpha_l / np.sqrt(2)))
    result[right] = left_total + core_total + right_total - np.exp(
        n_r * np.log(n_r / alpha_r) - alpha_r**2 / 2 - np.log(n_r - 1)
        + (1 - n_r) * np.log(b_r + z[right]))
    return result


def mass_probabilities(edges, pars):
    mean, sigma, alpha_l, n_l, alpha_r, n_r = pars
    cdf = crystal_ball_integral((edges - mean) / sigma, alpha_l, n_l, alpha_r, n_r)
    p = np.maximum(np.diff(cdf), 1e-16)
    return p / p.sum()


def poisson_deviance_residual(observed, expected):
    expected = np.maximum(expected, 1e-14)
    term = expected - observed
    positive = observed > 0
    term[positive] += observed[positive] * np.log(observed[positive] / expected[positive])
    return np.sign(observed - expected) * np.sqrt(np.maximum(2 * term, 0))


def fit_mass(values, edges, state="1S"):
    """Fit one state's signal; the two-argument API retains its 1S meaning."""
    if state not in MASSES:
        raise ValueError(f"Unknown Upsilon state: {state}")
    counts, _ = np.histogram(values, edges)
    n = int(counts.sum())
    if n < 100:
        raise RuntimeError("At least 100 selected MC events are required for the mass fit")
    lower = np.array([9.2, 0.005, 0.1, 1.01, 0.1, 1.01])
    upper = np.array([9.7, 0.3, 8.0, 100.0, 8.0, 100.0])
    # Translate the original mean range instead of constraining every fit to 1S.
    lower[0] += MASSES[state] - MASS
    upper[0] += MASSES[state] - MASS
    initial = [float(np.median(values)), 0.05, 1.5, 3.0, 1.5, 3.0]
    result = least_squares(lambda p: poisson_deviance_residual(counts, n * mass_probabilities(edges, p)),
                           np.clip(initial, lower + 1e-5, upper - 1e-5), bounds=(lower, upper),
                           x_scale="jac", max_nfev=2500, ftol=1e-10, xtol=1e-10, gtol=1e-8)
    expected = n * mass_probabilities(edges, result.x)
    deviance = float(np.sum(poisson_deviance_residual(counts, expected)**2))
    ndof = len(counts) - len(result.x) - 1
    rank = int(np.linalg.matrix_rank(result.jac))
    bound = np.isclose(result.x, lower, atol=1e-4, rtol=1e-3) | np.isclose(result.x, upper, atol=1e-4, rtol=1e-3)
    covariance = np.linalg.pinv(result.jac.T @ result.jac)
    errors = np.sqrt(np.maximum(np.diag(covariance), 0))
    names = ["mean_GeV", "sigma_core_GeV", "alpha_left", "n_left", "alpha_right", "n_right"]
    summary = {
        "state": state, "nominal_mass_GeV": MASSES[state],
        "model": "signal-only double-sided Crystal Ball, normalized in fit interval",
        "likelihood": "binned Poisson likelihood of unweighted MC event counts",
        "fit_interval_GeV": [float(edges[0]), float(edges[-1])],
        "bin_width_GeV": float(edges[1] - edges[0]), "mc_events_fitted": n,
        "converged": bool(result.success), "optimizer_message": result.message,
        "poisson_deviance": deviance, "ndof": ndof, "deviance_per_ndof": deviance / ndof,
        "asymptotic_deviance_pvalue": float(chi2.sf(deviance, ndof)),
        "pvalue_caveat": "asymptotic approximation; sparse tail bins can invalidate chi-square calibration",
        "covariance_rank": rank,
        "mc_covariance_parameter_order": names,
        "mc_covariance": covariance.tolist(),
        "mc_covariance_method": "unscaled Gauss-Newton inverse of the Poisson-deviance residual Jacobian; no chi-square/ndf rescaling",
        "uncertainty_interpretation": "approximate local likelihood curvature from finite MC sample, not projected-data precision",
        "parameters": {name: {"value": float(value), "mc_stat_error": float(error),
                              "at_bound": bool(at_bound), "bounds": [float(lo), float(hi)]}
                       for name, value, error, at_bound, lo, hi in zip(names, result.x, errors, bound, lower, upper)},
        "model_adequacy": "DSCB is a descriptive baseline only; significant shape mismatch" if deviance / ndof > 2 else "inspect residuals; DSCB sigma describes its fitted core, not full tail resolution",
        "shape_mismatch_flag": bool(deviance / ndof > 2),
    }
    if not result.success or rank < 6 or np.any(bound):
        summary["fit_warning"] = "Convergence, parameter-bound, or covariance issue; inspect before interpreting fit errors."
    return counts, expected, summary


def selection_masks(arrays, state, electron_energy, q2_cut, w_min, w_max,
                    true_mass_window=0.2):
    """Select one pure-state tuple while retaining its full reconstructed mass tails."""
    mass = MASSES[state]
    truth = (np.isfinite(arrays["true_mass"])
             & (np.abs(arrays["true_mass"] - mass) <= true_mass_window)
             & np.isfinite(arrays["true_Q2"]) & (arrays["true_Q2"] >= 0)
             & (arrays["true_Q2"] <= q2_cut)
             & np.isfinite(arrays["true_W"]) & (arrays["true_W"] >= w_min)
             & (arrays["true_W"] <= w_max))
    reco = (truth & np.isfinite(arrays["reco_mass"])
            & (arrays["reco_mass"] >= 8) & (arrays["reco_mass"] <= 12)
            & np.isfinite(arrays["reco_W"]) & (arrays["reco_W"] >= w_min)
            & (arrays["reco_W"] <= w_max)
            & np.isfinite(arrays["sigma_epz"])
            & (np.abs(arrays["sigma_epz"] - 2 * electron_energy) <= 0.2 * electron_energy))
    reco &= proton_tag_mask(arrays)
    return truth, reco



def statistical_configuration(stat_mode="mc", source_lumi_invfb=10.0,
                              projection_lumi_invfb=None, projection_shape_mode="floating"):
    """Keep source normalization explicit and scale projection weights exactly once."""
    if stat_mode not in ("mc", "projected"):
        raise ValueError("stat_mode must be mc or projected")
    if projection_shape_mode not in ("floating", "fixed_tails"):
        raise ValueError("projection_shape_mode must be floating or fixed_tails")
    target = source_lumi_invfb if projection_lumi_invfb is None else projection_lumi_invfb
    if not np.isfinite(source_lumi_invfb) or source_lumi_invfb <= 0 or not np.isfinite(target) or target <= 0:
        raise ValueError("Source and projection luminosities must be finite and positive")
    scale = target / source_lumi_invfb if stat_mode == "projected" else 1.0
    return {"stat_mode": stat_mode, "source_weight_luminosity_invfb": float(source_lumi_invfb),
            "requested_projection_luminosity_invfb": float(target),
            "projection_luminosity_invfb": float(target), "active_luminosity_invfb": float(target if stat_mode == "projected" else source_lumi_invfb),
            "projection_luminosity_ignored_in_mc_mode": stat_mode == "mc",
            "yield_weight_scale": float(scale), "projection_shape_mode": projection_shape_mode,
            "displayed_yield_variance": "expected_count" if stat_mode == "projected" else "sum_squared_weights",
            "efficiency_uncertainty": "finite MC Wilson intervals, unaffected by statistical mode or luminosity",
            "resolution_calibration": "finite MC quantities, unaffected by statistical mode or luminosity",
            "components_combined": False,
            "normalization": "Input weights already represent source luminosity; only target/source scaling is applied"}


def mass_projection_forecast(edges, calibration, expected_yield, shape_mode="floating"):
    """Expected extended-Poisson Fisher covariance at the MC-calibrated DSCB.

    This is a separate signal-only forecast. It does not reinterpret the MC
    calibration covariance as a data uncertainty or constrain data shape
    parameters using that covariance. Fixed tails mean exactly fixed tails,
    without a hidden MC constraint. No mismatch-based covariance inflation is used.
    """
    shape_names = ["mean_GeV", "sigma_core_GeV", "alpha_left", "n_left", "alpha_right", "n_right"]
    names = ["yield_events"] + shape_names
    parameters = np.array([calibration["parameters"][name]["value"] for name in shape_names], float)
    total = float(expected_yield)
    active = list(range(7)) if shape_mode == "floating" else [0, 1, 2]
    summary = {
        "status": "not_evaluated", "valid_local_covariance": False,
        "model": "extended signal-only DSCB Poisson model, lambda_i = N * p_i(theta)",
        "asimov_definition": "Expected bin counts from the MC-calibrated DSCB at the projected selected yield",
        "shape_mode": shape_mode, "parameter_order": names,
        "floating_parameters": [names[i] for i in active],
        "fixed_parameters": [names[i] for i in range(7) if i not in active],
        "expected_selected_yield": total, "background_model": "none",
        "calibration_shape_mismatch_flag": bool(calibration.get("shape_mismatch_flag", False)),
        "calibration_covariance_added": False, "MC_shape_constraint_applied": False,
        "chi_square_rescaling_applied": False,
        "interpretation": "Conditional expected data-statistical precision assuming the calibrated signal model and no background; not total experimental precision",
        "limitations": ["Model and calibration uncertainty, finite-MC template uncertainty, background, efficiencies, luminosity, and detector systematics are excluded",
                        "Local Fisher approximation can be inaccurate at small projected yields or near parameter boundaries"],
    }
    if not np.isfinite(total) or total <= 0:
        summary.update(status="invalid_expected_yield", covariance=None)
        return summary
    if shape_mode not in ("floating", "fixed_tails"):
        raise ValueError("shape_mode must be floating or fixed_tails")
    probabilities = mass_probabilities(np.asarray(edges, float), parameters)
    expected = total * probabilities
    summary["asimov_expected_bin_counts"] = expected.tolist()
    summary["low_expected_statistics_flag"] = total < 100
    summary["calibration_boundary_parameters"] = [name for name in shape_names if calibration["parameters"][name].get("at_bound", False)]
    if not calibration.get("converged", False) or summary["calibration_boundary_parameters"]:
        summary.update(status="invalid_calibration_or_boundary", covariance=None)
        return summary
    # Dimensionless derivatives improve conditioning without changing the model.
    scales = np.array([total, parameters[1], parameters[1], max(parameters[2], 1),
                       max(parameters[3], 1), max(parameters[4], 1), max(parameters[5], 1)])
    derivative = np.empty((len(probabilities), 7))
    derivative[:, 0] = probabilities * scales[0]
    for j in range(6):
        step = 1e-4 * scales[j + 1]
        plus, minus = parameters.copy(), parameters.copy()
        plus[j] += step
        minus[j] -= step
        derivative[:, j + 1] = total * (mass_probabilities(edges, plus) - mass_probabilities(edges, minus)) / (2 * step) * scales[j + 1]
    jacobian = derivative[:, active]
    fisher = (jacobian.T / expected) @ jacobian
    fisher = (fisher + fisher.T) / 2
    singular = np.linalg.svd(fisher, compute_uv=False)
    tolerance = np.finfo(float).eps * len(active) * singular[0]
    rank = int(np.sum(singular > tolerance))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else None
    summary.update(fisher_parameter_order=[names[i] for i in active],
                   fisher_parameter_scales=scales[active].tolist(),
                   fisher_scaled=fisher.tolist(), fisher_rank=rank, fisher_condition_scaled=condition)
    if rank < len(active) or condition is None or condition > 1e12:
        summary.update(status="singular_or_ill_conditioned_fisher", covariance=None)
        return summary
    try:
        inverse = np.linalg.inv(fisher)
    except np.linalg.LinAlgError:
        summary.update(status="fisher_inversion_failed", covariance=None)
        return summary
    covariance_active = inverse * scales[active, None] * scales[None, active]
    if np.any(~np.isfinite(covariance_active)) or np.any(np.diag(covariance_active) <= 0):
        summary.update(status="invalid_fisher_covariance", covariance=None)
        return summary
    covariance = np.zeros((7, 7))
    covariance[np.ix_(active, active)] = covariance_active
    errors = np.sqrt(np.diag(covariance))
    values = np.r_[total, parameters]
    summary.update(status="ok_conditional_model", valid_local_covariance=True,
                   covariance=covariance.tolist(),
                   parameters={name: {"value": float(value), "projected_stat_error": float(error), "fixed": i not in active}
                               for i, (name, value, error) in enumerate(zip(names, values, errors))})
    if summary["calibration_shape_mismatch_flag"]:
        summary["limitations"].append("The MC calibration shows significant DSCB shape mismatch; this forecast does not validate that model")
    return summary


def mass_interval_label(method, coverage):
    if method == "median_garwood":
        return f"Median expected Garwood intervals ({100 * coverage:.2f}%)"
    return f"Poisson Gaussian bands ({100 * coverage:.2f}%)"


def projected_yield_interval(expected_yield, coverage, method):
    """Keep count confidence-limit forecasts separate from Fisher covariance."""
    interval = expected_count_intervals(np.array([expected_yield]), coverage=coverage, method=method)
    reference = interval["reference_count"]
    return {"central_expected_yield": float(expected_yield),
            "lower": float(interval["lower"][0]), "upper": float(interval["upper"][0]),
            "reference_count": float(reference[0]) if reference is not None and np.isfinite(reference[0]) else None,
            "method": interval["method"], "coverage": float(interval["coverage"]),
            "poisson_standard_deviation": float(np.sqrt(expected_yield)),
            "interpretation": ("Median expected endpoints of exact central Garwood confidence intervals over hypothetical integer-count experiments; the expected-mean marker is not an observed count"
                               if method == "median_garwood" else "Normal-approximation band using the Poisson standard-deviation scale, not an exact confidence interval")}


def write_csv(path, rows):
    if rows:
        with path.open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def finish(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def wilson(k, n):
    n = np.asarray(n, float)
    p = np.divide(k, n, out=np.zeros_like(n), where=n > 0)
    safe = np.where(n > 0, n, 1)
    center = (p + 0.5 / safe) / (1 + 1 / safe)
    half = np.sqrt(p * (1 - p) / safe + 0.25 / safe**2) / (1 + 1 / safe)
    return p, center - half, center + half


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--state", choices=list(MASSES), default="1S")
    parser.add_argument("--electron-energy", type=float, default=9.0)
    parser.add_argument("--energy-label", default="9x130 GeV")
    parser.add_argument("--q2-cut", type=float, default=1.0)
    parser.add_argument("--true-mass-window", type=float, default=0.2)
    parser.add_argument("--w-min", type=float, default=20.0)
    parser.add_argument("--w-max", type=float, default=70.0)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--stat-mode", choices=["mc", "projected"], default="mc")
    parser.add_argument("--lumi-invfb", type=float, default=10.0, help="Luminosity already represented by input weights")
    parser.add_argument("--projection-lumi-invfb", type=float)
    parser.add_argument("--projection-shape-mode", choices=["floating", "fixed_tails"], default="floating")
    add_interval_arguments(parser)
    args = parser.parse_args()
    if args.w_min >= args.w_max or args.q2_cut <= 0 or args.electron_energy <= 0 or args.true_mass_window <= 0:
        parser.error("Require w-min < w-max and positive q2-cut, electron-energy and true-mass-window")
    try:
        stat = statistical_configuration(args.stat_mode, args.lumi_invfb, args.projection_lumi_invfb, args.projection_shape_mode)
        expected_count_intervals(np.array([0.]), coverage=args.interval_coverage, method=args.count_interval)
    except ValueError as error:
        parser.error(str(error))
    output = Path(args.out_dir)
    stat["count_interval"] = {"method": args.count_interval, "coverage": args.interval_coverage,
                              "applied": args.stat_mode == "projected", "covariance_changed": False}
    output.mkdir(parents=True, exist_ok=True)
    branches = ["weight", "true_mass", "true_Q2", "true_W", "true_y", "true_pT",
                "reco_mass", "reco_W", "sigma_epz"]
    branches += proton_branches()
    histograms = {}
    with uproot.open(args.input) as source:
        tree = source[args.tree]
        arrays = {name: tree[name].array(library="np") for name in branches}
        for name in ["h_EventCount", "h_CutFlow", "h_EleID_CutFlow", "h_EleID_EpRatio", "h_EleID_IsoFraction"]:
            if name in source:
                hist = source[name]
                histograms[name] = {"counts": hist.values(), "edges": hist.axis().edges(),
                                    "labels": hist.axis().labels(), "flow": hist.values(flow=True)[[0, -1]]}
    truth, reco = selection_masks(arrays, args.state, args.electron_energy, args.q2_cut,
                                  args.w_min, args.w_max, args.true_mass_window)
    weights = np.asarray(arrays["weight"], float)
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise RuntimeError("This diagnostic requires finite, nonnegative normalized tuple weights")
    selected_weight = weights[reco] * stat["yield_weight_scale"]
    if not selected_weight.size or selected_weight.sum() <= 0:
        raise RuntimeError("No selected events with positive total weight")
    # A direct rate rescaling of an unweighted shape fit requires constant event weights.
    constant_weights = bool(np.allclose(selected_weight, selected_weight[0], rtol=1e-8, atol=1e-14))
    if args.stat_mode == "projected" and not constant_weights:
        raise RuntimeError("Projected mass-parameter forecasting requires constant normalized event weights")
    mass = arrays["reco_mass"][reco]
    residual = mass - arrays["true_mass"][reco]
    quantiles = np.quantile(residual, [0.025, 0.16, 0.5, 0.84, 0.975])
    edges = np.linspace(8, 12, 161)
    centers = (edges[:-1] + edges[1:]) / 2
    counts, expected, fit = fit_mass(mass, edges, state=args.state)
    projected, _ = np.histogram(mass, edges, weights=selected_weight)
    projected_sumw2, _ = np.histogram(mass, edges, weights=selected_weight**2)
    poisson_variance = projected.copy()
    shown_variance = poisson_variance if args.stat_mode == "projected" else projected_sumw2
    intervals = (expected_count_intervals(projected, coverage=args.interval_coverage, method=args.count_interval)
                 if args.stat_mode == "projected" else None)
    yield_interval = (projected_yield_interval(selected_weight.sum(), args.interval_coverage, args.count_interval)
                      if args.stat_mode == "projected" else None)
    forecast = mass_projection_forecast(edges, fit, selected_weight.sum(), args.projection_shape_mode) if args.stat_mode == "projected" else {"status": "not_requested", "valid_local_covariance": False}
    fit["projected_forecast"] = forecast
    fit["statistical_configuration"] = stat
    fit["projected_yield_interval"] = yield_interval
    pars = fit["parameters"]
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    axes[0].errorbar(centers, counts, yerr=np.sqrt(counts), fmt=".", ms=4, color="black", label="Selected signal MC")
    axes[0].stairs(expected, edges, color="tab:red", label="Double-sided Crystal Ball")
    axes[0].set_yscale("log")
    axes[0].set_ylim(0.4, max(counts.max(), expected.max()) * 6)
    axes[0].set_ylabel("MC events / 0.025 GeV")
    axes[0].set_title(f"Upsilon({args.state}) mass fit, {args.energy_label}")
    axes[0].legend(loc="upper left", frameon=False, fontsize=10)
    note = (f"Mean = {pars['mean_GeV']['value']:.5f} GeV\n"
            f"Core sigma = {1000 * pars['sigma_core_GeV']['value']:.2f} MeV\n"
            f"Poisson deviance / ndf = {fit['poisson_deviance']:.1f} / {fit['ndof']}\n"
            f"Residual central 68% half-width = {500 * (quantiles[3] - quantiles[1]):.1f} MeV\n"
            "Fit uncertainties describe finite MC statistics")
    if fit["shape_mismatch_flag"]:
        note += "\nDSCB baseline has significant shape mismatch"
    axes[0].text(0.97, 0.96, note, transform=axes[0].transAxes, ha="right", va="top", fontsize=9)
    axes[1].axhline(0, color="gray", lw=1)
    axes[1].plot(centers, poisson_deviance_residual(counts, expected), ".", color="black", ms=4)
    axes[1].set_ylabel("Deviance\nresidual")
    axes[1].set_xlabel(r"Reconstructed $m_{ee}$ [GeV]")
    axes[1].set_xlim(8, 12)
    finish(fig, output / "mass_fit_mc.png")
    fig, ax = plt.subplots(figsize=(9, 6))
    if intervals is not None:
        draw_absolute_intervals(ax, centers, projected, intervals["lower"], intervals["upper"],
                                fmt=".", ms=4, color="black", label="Expected mean; count intervals")
    else:
        ax.errorbar(centers, projected, yerr=np.sqrt(shown_variance), fmt=".", ms=4, color="black",
                    label="Normalized tuple weights (MC error)")
    if constant_weights:
        ax.stairs(expected * selected_weight[0], edges, color="tab:red", label="MC fit scaled by event weight")
    ax.set_ylabel("Projected events / 0.025 GeV")
    ax.set_xlabel(r"Reconstructed $m_{ee}$ [GeV]")
    ax.set_title(f"Upsilon({args.state}), {args.energy_label}, {stat['active_luminosity_invfb']:g} fb$^{{-1}}$")
    ax.set_xlim(8, 12)
    if args.stat_mode == "projected":
        ax.set_ylim(0, 1.45 * max(float(np.max(intervals["upper"])), float(projected.max()), 1.))
    else:
        ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=10, loc="upper left" if args.stat_mode == "projected" else "best")
    projection_note = f"Selected expected yield = {selected_weight.sum():.3f}\n"
    if yield_interval is not None:
        projection_note += f"Expected yield bounds = [{yield_interval['lower']:.2f}, {yield_interval['upper']:.2f}]\n"
    projection_note += mass_interval_label(args.count_interval, args.interval_coverage) if args.stat_mode == "projected" else "Error bars: finite MC uncertainty"
    if forecast.get("valid_local_covariance"):
        projection_note += (f"\nExpected mean error = {1000 * forecast['parameters']['mean_GeV']['projected_stat_error']:.2f} MeV"
                            f"\nExpected core-width error = {1000 * forecast['parameters']['sigma_core_GeV']['projected_stat_error']:.2f} MeV"
                            f"\nLocal signal-only Fisher precision ({args.projection_shape_mode})")
        if forecast.get("low_expected_statistics_flag"):
            projection_note += "\nLow yield: local shape approximation may be inaccurate"
    elif args.stat_mode == "projected":
        projection_note += "\nProjected fit covariance unavailable: " + forecast["status"]
    if args.stat_mode == "projected" and fit["shape_mismatch_flag"]:
        projection_note += "\nConditional on DSCB; MC shape mismatch remains"
    ax.text(0.97, 0.94, projection_note, transform=ax.transAxes, ha="right", va="top", fontsize=9)
    finish(fig, output / "mass_projected.png")
    write_csv(output / "mass_bins.csv", [dict(low=float(edges[i]), high=float(edges[i+1]), mc_count=int(counts[i]),
              fit_mc_expectation=float(expected[i]), projected_yield=float(projected[i]), projected_mc_variance=float(projected_sumw2[i]), projected_data_variance=float(poisson_variance[i]),
              displayed_variance=float(shown_variance[i]), stat_mode=args.stat_mode,
              projected_interval_low=float(intervals["lower"][i]) if intervals is not None else None,
              projected_interval_high=float(intervals["upper"][i]) if intervals is not None else None,
              interval_reference_count=float(intervals["reference_count"][i]) if intervals is not None and intervals["reference_count"] is not None and np.isfinite(intervals["reference_count"][i]) else None,
              count_interval_method=args.count_interval if intervals is not None else "not_applied_mc",
              interval_coverage=args.interval_coverage if intervals is not None else None,
              mc_effective_entries=float(projected[i]**2/projected_sumw2[i]) if projected_sumw2[i] > 0 else None,
              model_asimov_expectation=float(forecast["asimov_expected_bin_counts"][i]) if "asimov_expected_bin_counts" in forecast else None) for i in range(len(counts))])
    fit["mass_residual_unweighted_quantiles_GeV"] = dict(zip(["q025", "q16", "q50", "q84", "q975"], map(float, quantiles)))
    fit["mass_residual_central68_halfwidth_GeV"] = float((quantiles[3] - quantiles[1]) / 2)
    fit["mass_residual_rms_GeV"] = float(np.std(residual))
    fit["proton_selection"] = proton_selection_metadata()
    fit["constant_selected_weights"] = constant_weights
    fit["selected_projected_yield"] = float(selected_weight.sum())
    (output / "mass_fit.json").write_text(json.dumps(fit, indent=2) + "\n")
    parameter_rows = []
    for name, parameter in pars.items():
        projected_error = forecast.get("parameters", {}).get(name, {}).get("projected_stat_error")
        parameter_rows.append(dict(parameter=name, **parameter, projected_stat_error=projected_error,
                                   projected_parameter_fixed=forecast.get("parameters", {}).get(name, {}).get("fixed"),
                                   selected_stat_error=projected_error if args.stat_mode == "projected" else parameter["mc_stat_error"],
                                   stat_mode=args.stat_mode))
    write_csv(output / "mass_fit_parameters.csv", parameter_rows)
    # Efficiencies refer strictly to eligible stored tuple rows, independent of yield weights.
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    efficiency_rows = []
    variables = [("true_W", np.linspace(args.w_min, args.w_max, 11), r"True $W_{\gamma p}$ [GeV]"),
                 ("true_y", np.linspace(-1, 3, 17), r"True $y_\Upsilon$"),
                 ("true_pT", np.linspace(0, 3, 16), r"True $p_{T,\Upsilon}$ [GeV]")]
    variable_overflows = {}
    for ax, (branch, bins, xlabel) in zip(axes, variables):
        valid = np.isfinite(arrays[branch]) & (np.abs(arrays[branch]) < 900)
        denominator, _ = np.histogram(arrays[branch][truth & valid], bins)
        numerator, _ = np.histogram(arrays[branch][reco & valid], bins)
        efficiency, low, high = wilson(numerator, denominator)
        populated = denominator > 0
        x = (bins[:-1] + bins[1:]) / 2
        ax.errorbar(x[populated], efficiency[populated], yerr=[efficiency[populated] - low[populated], high[populated] - efficiency[populated]], fmt="o", ms=4, capsize=2)
        ax.set_xlabel(xlabel)
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.2)
        variable_overflows[branch] = {"truth_below_range": int(np.sum(truth & valid & (arrays[branch] < bins[0]))),
                                     "truth_above_range": int(np.sum(truth & valid & (arrays[branch] > bins[-1])))}
        for i in range(len(denominator)):
            efficiency_rows.append(dict(variable=branch, low=float(bins[i]), high=float(bins[i+1]), denominator=int(denominator[i]),
                                        selected=int(numerator[i]), efficiency=float(efficiency[i]) if denominator[i] else None,
                                        wilson_low=float(low[i]) if denominator[i] else None, wilson_high=float(high[i]) if denominator[i] else None))
    axes[0].set_ylabel("Selected / eligible tuple rows")
    fig.suptitle(f"Upsilon({args.state}), {args.energy_label}: tuple-conditional selection efficiency\nUnweighted MC counts, 68% Wilson intervals", fontsize=12)
    finish(fig, output / "selection_efficiency.png")
    write_csv(output / "selection_efficiency.csv", efficiency_rows)
    # Stored producer diagnostics are not forced to close with the selection above.
    stored_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    for ax, name in zip(axes, ["h_CutFlow", "h_EleID_CutFlow"]):
        if name not in histograms:
            ax.text(0.5, 0.5, f"{name} absent", transform=ax.transAxes, ha="center")
            ax.set_axis_off()
            continue
        hist = histograms[name]
        raw_labels = hist["labels"] or [str(i+1) for i in range(len(hist["counts"]))]
        labels = [label.replace("#Sigma", "Sigma") for label in raw_labels]
        if name == "h_CutFlow" and len(hist["counts"]) == 10:
            labels = [
                "1. Total processed",
                "2. Generated Upsilon candidate found",
                "3. Generated decay leptons in acceptance",
                "4. Identified central e+ available",
                "5. Identified central e+ and e- available",
                "6. Additional e- candidate available",
                "7. Valid assignment and mass [8, 12] GeV",
                "8. Producer longitudinal-balance selection",
                "9. Reconstructed W >= 0",
                "10. RP/OMD proton candidate",
            ]
            ax.axhline(2.5, color="gray", linestyle="--", linewidth=1)
        ax.barh(np.arange(len(labels)), hist["counts"], color="tab:blue")
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Stored raw counts")
        ax.set_title(name)
        for i, (label, count) in enumerate(zip(raw_labels, hist["counts"])):
            stored_rows.append(dict(histogram=name, bin=i+1, label=label, count=float(count)))
    fig.suptitle(f"{args.energy_label} | stored producer counters\nEvent and electron-ID histograms retain their producer counting units", fontsize=12)
    finish(fig, output / "stored_cutflows.png")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    for ax, name, xlabel in zip(axes, ["h_EleID_EpRatio", "h_EleID_IsoFraction"], [r"Electron candidate $E/p$", "Electron candidate isolation fraction"]):
        if name not in histograms:
            ax.text(0.5, 0.5, f"{name} absent", transform=ax.transAxes, ha="center")
            ax.set_axis_off()
            continue
        hist = histograms[name]
        ax.stairs(hist["counts"], hist["edges"], color="tab:blue")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Stored raw candidate counts / bin")
        ax.set_yscale("log")
        for i, count in enumerate(hist["counts"]):
            stored_rows.append(dict(histogram=name, bin=i+1, label=f"[{hist['edges'][i]:g}, {hist['edges'][i+1]:g})", count=float(count)))
    fig.suptitle(f"{args.energy_label}: stored electron-ID inputs, producer selections", fontsize=12)
    finish(fig, output / "stored_electron_id.png")
    write_csv(output / "stored_histogram_counts.csv", stored_rows)
    audit = {"proton_selection": proton_selection_metadata(),
             "truth_rows_with_proton_tag": int(np.sum(truth & proton_tag_mask(arrays))),
             "input": str(Path(args.input).resolve()), "tree": args.tree, "state": args.state,
             "energy_label": args.energy_label, "electron_energy_GeV": args.electron_energy,
             "truth_selection": {"mass_window_GeV": [MASSES[args.state]-args.true_mass_window, MASSES[args.state]+args.true_mass_window], "Q2_GeV2": [0, args.q2_cut], "W_GeV": [args.w_min, args.w_max]},
             "reco_selection": {"mass_GeV": [8, 12], "W_GeV": [args.w_min, args.w_max],
                                "sigma_epz_GeV": [1.8*args.electron_energy, 2.2*args.electron_energy]},
             "tuple_rows": len(weights), "eligible_truth_rows": int(truth.sum()), "selected_rows": int(reco.sum()),
             "conditional_efficiency": float(reco.sum()/truth.sum()) if truth.sum() else None,
             "denominator_caveat": "Eligible stored UpsilonTree rows only; not a claim of full generation or simulation efficiency",
             "statistical_configuration": stat,
             "projected_data_variance_definition": "Expected independent-event Poisson variance equals expected yield; not an observed-count confidence interval",
             "projected_yield_interval": yield_interval,
             "interval_display_note": "Projected count-interval endpoints do not alter the expected mean, variance, or local Fisher covariance; zero-mean bins remain at zero on a linear axis",
             "mass_projected_axis_scale": "linear" if args.stat_mode == "projected" else "log",
             "empty_bin_caveat": "A zero MC-template bin is not proof of zero physical expectation; model Asimov expectations are recorded separately",
             "selected_sum_weight": float(selected_weight.sum()), "selected_sum_weight2": float(np.square(selected_weight).sum()),
             "stored_histogram_caveat": "Raw producer counts preserved as stored, with producer selections and counting units; no luminosity rescaling",
             "histograms_available": list(histograms),
             "stored_histogram_sums": {name: float(h["counts"].sum()) for name, h in histograms.items()},
             "stored_histogram_flow_counts": {name: list(map(float, h['flow'])) for name, h in histograms.items()},
             "efficiency_variable_overflows": variable_overflows}
    (output / "performance_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"eligible_truth_rows": int(truth.sum()), "selected_rows": int(reco.sum()),
                      "mass_fit_converged": fit["converged"], "mass_sigma_core_GeV": pars["sigma_core_GeV"]["value"],
                      "out_dir": str(output)}, indent=2))


if __name__ == "__main__":
    main()
