#!/usr/bin/env python3
"""
Extract differential exclusive Upsilon photoproduction cross sections.

The script computes bin-by-bin efficiency-corrected differential cross sections
for one or more Upsilon states and variables such as pT^2 and rapidity:

  B_ee dσ(γp -> Υp)/dpT^2
  B_ee dσ(γp -> Υp)/dy

or, optionally, the inclusive Υ cross section after dividing by B_ee.

Default convention
------------------
The eSTARlight tuples used here are Υ(nS) -> e+e- samples and the tuple weights
are visible e+e- weights. Therefore the default output convention is
"visible_decay", i.e. B_ee dσ/dx. Use

  --output-convention inclusive_upsilon

to divide by the dielectron branching fraction.

Important command-line note
---------------------------
For negative bin edges, either write

  --y-bins=-1,0,1,2

or quote it. This v4 script also pre-processes the argv so quoted negative
lists work in most shells.

Examples
--------
All states at 18x275:

python3 extract_differential_gamma_p_swan_v4.py /eos/user/g/gkrintir/merged_elecs_test_tuples_fix2.root \
  --states 1S,2S,3S \
  --electron-energy 18 \
  --proton-energy 275 \
  --visible-xsec-1s-fb 327.584 \
  --visible-xsec-2s-fb 84.780 \
  --visible-xsec-3s-fb 68.929 \
  --w-min 37 \
  --w-max 114 \
  --q2-max 0.001 \
  --q2-cut 0.001 \
  --variables pT2,y \
  --y-bins=-1,0,1,2 \
  --output-convention visible_decay \
  --out-prefix diff_gamma_p_family_18x275_W37to114
"""

import argparse
from pathlib import Path as _BranchingPath
import sys as _branching_sys
_branching_sys.path.insert(0, str(_BranchingPath(__file__).resolve().parents[1]))
_branching_sys.path.insert(0, str(_BranchingPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata
from branching_fractions import (add_branching_fraction_argument, branching_fractions_from_args,
                                branching_fraction_metadata)
import csv
import math
import os
import sys
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join(tempfile.gettempdir(), f"{os.environ.get('USER', 'user')}_mplconfig")
)

import uproot

from statistical_components import (configure_statistics, metadata as stat_metadata, statistical_label,
                                    binned_components, covariance_rows, write_covariance,
                                    add_interval_arguments, binned_intervals, draw_interval_bars)
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813

STATES = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0234, "color": "royalblue"},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191, "color": "crimson"},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218, "color": "forestgreen"},
}

DEFAULT_VISIBLE_XSECS = {
    "18x275": {"1S": 327.584, "2S": 84.780, "3S": 68.929},
    "10x100": {"1S": 104.793, "2S": 29.297, "3S": 23.447},
}

VARIABLES = {
    "pT2": {
        "true_expr": "true_pT2",
        "reco_expr": "reco_pT2",
        "xlabel": r"$p_{T,\Upsilon}^{2}$ [GeV$^2$]",
        "ylabel_unit": r"pb/GeV$^2$",
        "default_bins": "0,0.05,0.10,0.20,0.35,0.50,0.75,1.0,1.5,2.0,3.0,5.0",
    },
    "y": {
        "true_expr": "true_y",
        "reco_expr": "reco_y",
        "xlabel": r"$y_{\Upsilon}$",
        "ylabel_unit": r"pb",
        "default_bins": "-1,0,1,2",
    },
    "W": {
        "true_expr": "true_W",
        "reco_expr": "reco_W",
        "xlabel": r"$W_{\gamma p}$ [GeV]",
        "ylabel_unit": r"pb/GeV",
        "default_bins": "20,30,40,50,60,70,80,90,100,115,130,150",
    },
    "t": {
        "true_expr": "true_t",
        "reco_expr": "reco_t",
        "xlabel": r"$|t|$ [GeV$^2$]",
        "ylabel_unit": r"pb/GeV$^2$",
        "default_bins": "0,0.5,1.0,5.0",
    },
}


def preprocess_negative_bin_args(argv):
    """Allow patterns like --y-bins "-4,-3,...".

    argparse may interpret a token beginning with '-' as an option. Convert
    [--y-bins, -4,-3,...] to [--y-bins=-4,-3,...].
    """
    bin_opts = {"--y-bins", "--pt2-bins", "--w-bins", "--t-bins"}
    out = []
    i = 0
    while i < len(argv):
        if argv[i] in bin_opts and i + 1 < len(argv):
            nxt = argv[i + 1]
            if nxt.startswith("-") and "," in nxt:
                out.append(f"{argv[i]}={nxt}")
                i += 2
                continue
        out.append(argv[i])
        i += 1
    return out


def parse_edges(s):
    return np.asarray([float(x) for x in str(s).split(",") if x.strip()], dtype=float)


def parse_states(s):
    states = [x.strip() for x in str(s).split(",") if x.strip()]
    for st in states:
        if st not in STATES:
            raise ValueError(f"Unknown state {st}; valid states are {list(STATES)}")
    return states


def photon_flux_dW(W, e_e, e_p, q2_max):
    s_minus_mp2 = 4.0 * e_e * e_p
    y = (W * W - M_P * M_P) / s_minus_mp2
    if y <= 0.0 or y >= 1.0:
        return 0.0

    q2_min = (M_E * M_E * y * y) / (1.0 - y)
    if q2_min <= 0.0 or q2_min >= q2_max:
        return 0.0

    log_term = math.log(q2_max / q2_min)
    dndy = (ALPHA_EM / (2.0 * math.pi)) * (
        ((1.0 + (1.0 - y) ** 2) / y) * log_term
        - (2.0 * (1.0 - y) / y) * (1.0 - q2_min / q2_max)
    )
    dydw = 2.0 * W / s_minus_mp2
    return dndy * dydw


def integrated_flux(w_low, w_high, e_e, e_p, q2_max, n_steps=600):
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.asarray([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:
        trapezoid = np.trapz
    return float(trapezoid(vals, grid))


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def branch_exists(path, tree, branch):
    with uproot.open(os.path.expanduser(path)) as source:
        return tree in source and branch in source[tree].keys()


def get_default_visible_xsec(state, electron_energy, proton_energy):
    if abs(electron_energy - 10.0) < 1e-6 and abs(proton_energy - 100.0) < 1e-6:
        return DEFAULT_VISIBLE_XSECS["10x100"][state]
    if abs(electron_energy - 18.0) < 1e-6 and abs(proton_energy - 275.0) < 1e-6:
        return DEFAULT_VISIBLE_XSECS["18x275"][state]
    return None


def get_visible_xsec_for_state(args, state):
    explicit = {
        "1S": args.visible_xsec_1s_fb,
        "2S": args.visible_xsec_2s_fb,
        "3S": args.visible_xsec_3s_fb,
    }[state]
    if explicit is not None:
        return explicit
    if args.visible_xsec_fb is not None and state == "1S":
        # backward-compatible single-state option
        return args.visible_xsec_fb
    return get_default_visible_xsec(state, args.electron_energy, args.proton_energy)


def weighted_exp_fit(x, y, yerr, covariance=None):
    """Log-space GLS, treating the supplied covariance as absolute, not chi2-scaled."""
    x, y, yerr = np.asarray(x), np.asarray(y), np.asarray(yerr)
    m = np.isfinite(x) & np.isfinite(y) & (y > 0) & np.isfinite(yerr) & (yerr > 0)
    if np.sum(m) < 2:
        return np.nan, np.nan, np.nan, np.nan
    covariance = np.diag(yerr**2) if covariance is None else np.asarray(covariance)
    log_cov = covariance[np.ix_(m, m)] / np.outer(y[m], y[m])
    try:
        precision = np.linalg.inv(log_cov)
        design = np.column_stack((x[m], np.ones(np.sum(m))))
        parameter_covariance = np.linalg.inv(design.T @ precision @ design)
        slope, intercept = parameter_covariance @ design.T @ precision @ np.log(y[m])
        slope_err, intercept_err = np.sqrt(np.diag(parameter_covariance))
        amplitude = math.exp(intercept)
    except (np.linalg.LinAlgError, ValueError, OverflowError):
        return np.nan, np.nan, np.nan, np.nan
    return amplitude, amplitude * intercept_err, -slope, slope_err


def _b_shape_fractions(b, edges):
    """Bin fractions for a shape proportional to b exp(-b t)."""
    edges = np.asarray(edges, dtype=float)
    if b <= 0:
        return np.zeros(len(edges) - 1)
    norm = math.exp(-b * edges[0]) - math.exp(-b * edges[-1])
    if norm <= 0:
        return np.zeros(len(edges) - 1)
    return (np.exp(-b * edges[:-1]) - np.exp(-b * edges[1:])) / norm


def fit_binned_exponential_slope(bin_edges, counts, errors, min_count=0.0, b_min=0.05, b_max=20.0, n_scan=4000, covariance=None):
    """Finite-bin GLS fit with free normalization and absolute input covariance.

    Profile Delta-chi2=1 intervals do not rescale errors by the reduced chi2.
    The supplied covariance includes all selected statistical components.
    """
    from scipy.optimize import minimize_scalar, brentq
    edges, counts, errors = np.asarray(bin_edges, float), np.asarray(counts, float), np.asarray(errors, float)
    valid = np.isfinite(counts) & np.isfinite(errors) & (errors > 0) & (counts > min_count)
    empty = {"b": np.nan, "b_err_low": np.nan, "b_err_high": np.nan,
             "norm_yield": np.nan, "chi2": np.nan, "ndof": 0,
             "valid_bins": int(np.sum(valid)), "covariance_rank": 0,
             "fit_uncertainty_method": "GLS profile Delta-chi2=1; absolute covariance; no reduced-chi2 scaling"}
    if np.sum(valid) < 2:
        return empty
    covariance = np.diag(errors**2) if covariance is None else np.asarray(covariance, float)
    cov = covariance[np.ix_(valid, valid)]
    if not np.isfinite(cov).all():
        raise ValueError("Nonfinite covariance in estimable b-fit bins")
    eigenvalues = np.linalg.eigvalsh((cov + cov.T) / 2)
    if eigenvalues[0] <= 0:
        raise ValueError("b-fit covariance must be positive definite in the selected bins")
    precision = np.linalg.inv(cov)
    y = counts[valid]

    def profile(b):
        f = _b_shape_fractions(b, edges)[valid]
        den = float(f @ precision @ f)
        norm = float(f @ precision @ y) / den if den > 0 else np.nan
        residual = y - norm * f
        chi2 = float(residual @ precision @ residual) if norm > 0 else np.inf
        return chi2, norm

    b_grid = np.linspace(b_min, b_max, max(int(n_scan), 3))
    chi2 = np.array([profile(b)[0] for b in b_grid])
    if not np.isfinite(chi2).any():
        return empty
    idx = int(np.nanargmin(chi2))
    lower, upper = b_grid[max(idx - 1, 0)], b_grid[min(idx + 1, len(b_grid) - 1)]
    optimum = minimize_scalar(lambda b: profile(b)[0], bounds=(lower, upper), method="bounded",
                              options={"xatol": 1e-11})
    candidates = [(float(b_grid[idx]), profile(b_grid[idx])[0]), (float(optimum.x), float(optimum.fun))]
    b_hat, chi2_min = min(candidates, key=lambda pair: pair[1])
    norm = profile(b_hat)[1]
    target = chi2_min + 1.0
    b_low = b_high = np.nan
    left = b_grid[(b_grid < b_hat) & (chi2 >= target)]
    right = b_grid[(b_grid > b_hat) & (chi2 >= target)]
    if left.size:
        b_low = brentq(lambda b: profile(b)[0] - target, float(left[-1]), b_hat)
    if right.size:
        b_high = brentq(lambda b: profile(b)[0] - target, b_hat, float(right[0]))
    return dict(empty, b=b_hat,
                b_err_low=b_hat-b_low if np.isfinite(b_low) else np.nan,
                b_err_high=b_high-b_hat if np.isfinite(b_high) else np.nan,
                norm_yield=norm, chi2=chi2_min, ndof=int(np.sum(valid)-2),
                covariance_rank=int(np.sum(valid)))


def fit_poisson_exponential_slope(bin_edges, reco_counts, corrections, valid=None,
                                  correction_covariance=None, b_min=0.05, b_max=20., n_scan=4000):
    """Reco-level Poisson likelihood, profiling a free corrected normalization.

    Fractional counts are Asimov expectations, not observed fractional events.
    Delta(-2 log L)=1 intervals have approximate asymptotic coverage. Optional
    MC calibration uses Gaussian log-correction nuisances with delta-method
    covariance; corrections are positive but may exceed one due to migration.
    """
    from scipy.optimize import minimize, minimize_scalar, brentq
    from scipy.special import logsumexp, xlogy
    edges = np.asarray(bin_edges, float)
    counts, correction = np.asarray(reco_counts, float), np.asarray(corrections, float)
    if np.any(np.isfinite(counts) & (counts < 0)):
        raise ValueError("Poisson fit counts must be nonnegative")
    valid = np.ones(counts.shape, bool) if valid is None else np.asarray(valid, bool).copy()
    valid &= np.isfinite(counts) & np.isfinite(correction) & (correction > 0)
    count = counts[valid]
    correction = correction[valid]
    total = float(count.sum())
    empty = dict(b=np.nan, b_err_low=np.nan, b_err_high=np.nan, norm_yield=np.nan,
                 chi2=np.nan, deviance=np.nan, constraint_penalty=0., fit_statistic_value=np.nan,
                 ndof=max(int(np.sum(valid))-2, 0), valid_bins=int(np.sum(valid)), covariance_rank=0,
                 calibration_nuisance_rank=0,
                 fit_statistic="poisson_deviance",
                 fit_uncertainty_method="Asimov reco-Poisson profile Delta(-2 log L)=1; approximate asymptotic coverage",
                 fit_interval_coverage=0.6826894921370859,
                 fit_calibration_model="Fixed bin-by-bin R/T corrections; MC-shape-dependent migration approximation",
                 fit_boundary_low=False, fit_boundary_high=False)
    if len(count) < 2 or total <= 0:
        return empty
    basis = np.zeros((len(count), 0))
    if correction_covariance is not None:
        cov = np.asarray(correction_covariance, float)[np.ix_(valid, valid)] / np.outer(correction, correction)
        if not np.isfinite(cov).all():
            raise ValueError("Nonfinite log-correction covariance in Poisson fit")
        eigenvalues, eigenvectors = np.linalg.eigh((cov+cov.T)/2)
        scale = max(float(np.max(np.abs(eigenvalues))), 1e-30)
        if np.min(eigenvalues) < -1e-10*scale:
            raise ValueError("Log-correction covariance must be positive semidefinite")
        keep = eigenvalues > 1e-12*scale
        basis = eigenvectors[:, keep]*np.sqrt(eigenvalues[keep])
        empty.update(fit_calibration_model="Correlated Gaussian log(R/T) nuisances; delta-method covariance; approximate MC calibration model",
                     fit_statistic="poisson_deviance_plus_gaussian_constraint", calibration_nuisance_rank=basis.shape[1])

    def profile(b):
        logshape = np.log(correction)+np.log(_b_shape_fractions(b, edges)[valid])
        def evaluate(z, gradient=False):
            logq = logshape+basis@z
            logsum = logsumexp(logq)
            mu = total*np.exp(logq-logsum)
            deviance = max(0., float(2*np.sum(mu-count+xlogy(count, count)-xlogy(count, mu))))
            penalty = float(z@z)
            value = deviance+penalty
            if gradient:
                return value, 2*basis.T@(mu-count)+2*z
            return value, total*np.exp(-logsum), deviance, penalty
        if basis.shape[1]:
            optimum = minimize(lambda z: evaluate(z, True), np.zeros(basis.shape[1]), jac=True,
                               method="BFGS", options={"gtol":1e-8})
            if not optimum.success and np.linalg.norm(optimum.jac, ord=np.inf) > 1e-5:
                raise RuntimeError("Poisson calibration-nuisance profiling did not converge")
            return evaluate(optimum.x)
        return evaluate(np.zeros(0))

    # A modest grid brackets the minimum; continuous minimization and root
    # crossings determine precision, keeping repeated low-count toy fits fast.
    grid_size = min(max(int(n_scan), 3), 101)
    grid = np.linspace(b_min, b_max, grid_size)
    values = np.array([profile(b)[0] for b in grid])
    index = int(np.argmin(values))
    lo, hi = grid[max(index-1, 0)], grid[min(index+1, len(grid)-1)]
    optimum = minimize_scalar(lambda b: profile(b)[0], bounds=(lo, hi), method="bounded",
                              options={"xatol":1e-10})
    b_hat, minimum = min([(float(grid[index]), float(values[index])), (float(optimum.x), float(optimum.fun))], key=lambda pair:pair[1])
    value, norm, deviance, penalty = profile(b_hat)
    threshold = minimum+1.
    left, right = grid[(grid < b_hat) & (values >= threshold)], grid[(grid > b_hat) & (values >= threshold)]
    low = brentq(lambda b:profile(b)[0]-threshold, float(left[-1]), b_hat) if left.size else np.nan
    high = brentq(lambda b:profile(b)[0]-threshold, b_hat, float(right[0])) if right.size else np.nan
    return dict(empty, b=b_hat, b_err_low=b_hat-low if np.isfinite(low) else np.nan,
                b_err_high=high-b_hat if np.isfinite(high) else np.nan, norm_yield=norm,
                chi2=value, deviance=deviance, constraint_penalty=penalty, fit_statistic_value=value,
                fit_boundary_low=not np.isfinite(low), fit_boundary_high=not np.isfinite(high))


def extract_b_slope_for_state(data, args, state, flux_int, lumi_pb):
    """ZEUS-style finite-bin b-slope extraction using pT2 or |t|."""
    props = dict(STATES[state], br_ee=branching_fractions_from_args(args)[state])
    br = props["br_ee"]
    visible_xsec = get_visible_xsec_for_state(args, state)
    if visible_xsec is None and args.state_normalize:
        raise RuntimeError(f"Need visible xsec for {state}; provide the per-state visible xsec or use --no-state-normalize.")

    truth_state_all = (
        (data["true_mass"] >= props["mass"] - args.truth_mass_window)
        & (data["true_mass"] <= props["mass"] + args.truth_mass_window)
    )

    raw_sumw_state = float(np.sum(data["weight"][truth_state_all]))
    if args.tuple_lumi_invfb <= 0:
        raise ValueError("--tuple-lumi-invfb must be positive")
    if raw_sumw_state <= 0:
        raise ValueError(
            f"No positive-weight truth events are available for Upsilon({state}); "
            "the requested b-slope extraction is undefined"
        )
    raw_visible_xsec_fb = raw_sumw_state / args.tuple_lumi_invfb if args.tuple_lumi_invfb > 0 else np.nan
    state_scale = visible_xsec / raw_visible_xsec_fb if (args.state_normalize and raw_visible_xsec_fb > 0) else 1.0
    weight = data["weight"] * state_scale

    sigma_epz_center = args.sigma_epz_center if args.sigma_epz_center is not None else 2.0 * args.electron_energy
    sigma_epz_halfwidth = args.sigma_epz_halfwidth if args.sigma_epz_halfwidth is not None else 0.2 * args.electron_energy

    if args.b_slope_zeus_mass_window and state == "1S":
        reco_mass_min = 9.33
        reco_mass_max = 9.66
    else:
        reco_mass_min = args.mass_min if args.mass_min is not None else props["mass"] - args.reco_mass_window_sigma * props["sigma"]
        reco_mass_max = args.mass_max if args.mass_max is not None else props["mass"] + args.reco_mass_window_sigma * props["sigma"]

    if args.b_slope_variable == "pT2":
        true_values = data["true_pT2"]
        reco_values = data["reco_pT2"]
    elif args.b_slope_variable == "t":
        true_values = data["true_t"]
        reco_values = data["reco_t"]
    else:
        raise ValueError(args.b_slope_variable)

    bins = parse_edges(args.b_slope_bins)
    selected_proton = proton_tag_mask(data)

    truth_base = (
        truth_state_all
        & (data["true_W"] >= args.w_min)
        & (data["true_W"] <= args.w_max)
        & (data["true_Q2"] >= 0.0)
        & (data["true_Q2"] <= args.q2_cut)
        & (true_values >= bins[0])
        & (true_values <= bins[-1])
    )

    reco_base = (
        truth_base
        & selected_proton
        & (data["reco_mass"] >= reco_mass_min)
        & (data["reco_mass"] <= reco_mass_max)
        & (np.abs(data["sigma_epz"] - sigma_epz_center) <= sigma_epz_halfwidth)
        & (data["reco_W"] >= args.w_min)
        & (data["reco_W"] <= args.w_max)
        & (reco_values >= bins[0])
        & (reco_values <= bins[-1])
    )

    stats = binned_components(true_values, reco_values, weight, truth_base, reco_base, bins, args)
    n_true, n_reco = stats["truth"], stats["reco"]
    eff, n_corr, err_corr = stats["efficiency"], stats["corrected"], stats["error"]

    valid = np.isfinite(eff) & (eff > 0) & (eff >= args.min_eff)
    if args.stat_mode == "projected":
        fit = fit_poisson_exponential_slope(
            bins, n_reco, eff, valid=valid,
            correction_covariance=stats["efficiency_mc_cov"] if args.include_efficiency_mc_uncertainty else None,
            b_min=args.b_slope_min, b_max=args.b_slope_max, n_scan=args.b_slope_n_scan)
    else:
        valid &= (n_corr > 0) & (err_corr > 0)
        fit = fit_binned_exponential_slope(
            bins, np.where(valid, n_corr, 0.), np.where(valid, err_corr, 0.),
            b_min=args.b_slope_min, b_max=args.b_slope_max, n_scan=args.b_slope_n_scan,
            covariance=stats["cov_selected"])
    intervals = binned_intervals(stats, args)

    br_corr = 1.0 / br if args.output_convention == "inclusive_upsilon" else 1.0
    sigma_total_pb = fit["norm_yield"] / (lumi_pb * flux_int) * br_corr if flux_int > 0 and np.isfinite(fit["norm_yield"]) else np.nan

    factor = br_corr / (lumi_pb * flux_int) if flux_int > 0 else np.nan
    write_covariance(f"{args.out_prefix}_{state}_b_slope_statistical_covariance.csv",
                     covariance_rows(stats, factor, bins, dict(state=state, variable=args.b_slope_variable, **stat_metadata(args))))
    bin_rows = []
    for i in range(len(bins) - 1):
        bin_rows.append({
            "proton_selection": proton_selection_metadata()["mode"],
            "state": state,
            "b_variable": args.b_slope_variable,
            "t_low": bins[i],
            "t_high": bins[i + 1],
            "t_center": 0.5 * (bins[i] + bins[i + 1]),
            "n_truth_visible": n_true[i],
            "n_reco_visible": n_reco[i],
            "efficiency": eff[i],
            "n_corrected_visible": n_corr[i],
            "n_corrected_visible_err": err_corr[i],
            "n_corrected_visible_interval_low": intervals["lower"][i],
            "n_corrected_visible_interval_high": intervals["upper"][i],
            "n_reco_count_interval_low": intervals["reco_lower"][i],
            "n_reco_count_interval_high": intervals["reco_upper"][i],
            "n_reco_interval_reference_count": intervals["reference_count"][i],
            "n_corrected_visible_err_mc": np.sqrt(stats["cov_mc"][i, i]),
            "n_corrected_visible_err_projected": np.sqrt(stats["cov_projected"][i, i]),
            "n_corrected_visible_err_correction_mc": np.sqrt(stats["cov_correction_mc"][i, i]),
            "used_in_b_fit": bool(valid[i]),
            **stat_metadata(args),
            **branching_fraction_metadata(args),
        })

    summary = {
        "proton_selection": proton_selection_metadata()["mode"],
        "state": state,
        "br_ee": br,
        "b_variable": args.b_slope_variable,
        "b_slope_GeVminus2": fit["b"],
        "b_slope_err_low": fit["b_err_low"],
        "b_slope_err_high": fit["b_err_high"],
        "fit_norm_corrected_visible_yield": fit["norm_yield"],
        "fit_integrated_sigma_pb": sigma_total_pb,
        "chi2": fit["chi2"],
        "ndof": fit["ndof"],
        "valid_bins": fit["valid_bins"],
        "mass_window_min": reco_mass_min,
        "mass_window_max": reco_mass_max,
        "b_fit_t_min": bins[0],
        "b_fit_t_max": bins[-1],
        "zeus_style_mass_window": bool(args.b_slope_zeus_mass_window and state == "1S"),
        "covariance_rank": fit["covariance_rank"],
        "fit_uncertainty_method": fit["fit_uncertainty_method"],
        "fit_statistic": fit.get("fit_statistic", "gaussian_gls_chi2"),
        "deviance": fit.get("deviance", np.nan),
        "constraint_penalty": fit.get("constraint_penalty", 0.),
        "fit_statistic_value": fit.get("fit_statistic_value", fit["chi2"]),
        "calibration_nuisance_rank": fit.get("calibration_nuisance_rank", 0),
        "fit_calibration_model": fit.get("fit_calibration_model", "Same-event full MC covariance"),
        "fit_interval_coverage": fit.get("fit_interval_coverage", 0.6826894921370859),
        "fit_boundary_low": fit.get("fit_boundary_low", not np.isfinite(fit["b_err_low"])),
        "fit_boundary_high": fit.get("fit_boundary_high", not np.isfinite(fit["b_err_high"])),
        **stat_metadata(args),
        **branching_fraction_metadata(args),
    }
    return summary, bin_rows


def extract_for_state(data, args, state, bins_by_var, flux_int, lumi_pb):
    props = dict(STATES[state], br_ee=branching_fractions_from_args(args)[state])
    br = props["br_ee"]

    visible_xsec = get_visible_xsec_for_state(args, state)
    if visible_xsec is None and args.state_normalize:
        raise RuntimeError(f"Need visible xsec for {state}; provide --visible-xsec-{state.lower()}-fb or use --no-state-normalize.")

    truth_state_all = (
        (data["true_mass"] >= props["mass"] - args.truth_mass_window)
        & (data["true_mass"] <= props["mass"] + args.truth_mass_window)
    )

    raw_sumw_state = float(np.sum(data["weight"][truth_state_all]))
    if args.tuple_lumi_invfb <= 0:
        raise ValueError("--tuple-lumi-invfb must be positive")
    if raw_sumw_state <= 0:
        raise ValueError(
            f"No positive-weight truth events are available for Upsilon({state}); "
            "refusing to emit zero or NaN differential cross sections"
        )
    raw_visible_xsec_fb = raw_sumw_state / args.tuple_lumi_invfb if args.tuple_lumi_invfb > 0 else np.nan

    state_scale = 1.0
    if args.state_normalize:
        state_scale = visible_xsec / raw_visible_xsec_fb

    weight = data["weight"] * state_scale

    sigma_epz_center = args.sigma_epz_center if args.sigma_epz_center is not None else 2.0 * args.electron_energy
    sigma_epz_halfwidth = args.sigma_epz_halfwidth if args.sigma_epz_halfwidth is not None else 0.2 * args.electron_energy

    reco_mass_min = args.mass_min if args.mass_min is not None else props["mass"] - args.reco_mass_window_sigma * props["sigma"]
    reco_mass_max = args.mass_max if args.mass_max is not None else props["mass"] + args.reco_mass_window_sigma * props["sigma"]

    truth_base = (
        truth_state_all
        & (data["true_W"] >= args.w_min)
        & (data["true_W"] <= args.w_max)
        & (data["true_Q2"] >= 0.0)
        & (data["true_Q2"] <= args.q2_cut)
    )

    reco_base = (
        truth_base
        & proton_tag_mask(data)
        & (data["reco_mass"] >= reco_mass_min)
        & (data["reco_mass"] <= reco_mass_max)
        & (np.abs(data["sigma_epz"] - sigma_epz_center) <= sigma_epz_halfwidth)
        & (data["reco_W"] >= args.w_min)
        & (data["reco_W"] <= args.w_max)
    )

    br_corr = 1.0 / br if args.output_convention == "inclusive_upsilon" else 1.0

    all_rows = []
    per_var_results = {}

    for var, bins in bins_by_var.items():
        widths = np.diff(bins)
        centers = 0.5 * (bins[:-1] + bins[1:])
        true_values = data[VARIABLES[var]["true_expr"]]
        reco_values = data[VARIABLES[var]["reco_expr"]]

        truth_mask = truth_base & (true_values >= bins[0]) & (true_values <= bins[-1])
        reco_mask = reco_base & (reco_values >= bins[0]) & (reco_values <= bins[-1])

        stats = binned_components(true_values, reco_values, weight, truth_mask, reco_mask, bins, args)
        intervals = binned_intervals(stats, args)
        n_true, n_reco = stats["truth"], stats["reco"]
        eff, n_corr, err_corr = stats["efficiency"], stats["corrected"], stats["error"]
        good = stats["valid"] & (eff >= args.min_eff) & (flux_int > 0)
        denom = lumi_pb * flux_int * widths
        factor = np.divide(br_corr, denom, out=np.full_like(denom, np.nan), where=denom > 0)
        dsigma, dsigma_err = n_corr * factor, err_corr * factor
        closure = np.divide(n_corr, n_true, out=np.full_like(n_corr, np.nan), where=n_true > 0)
        covariance = stats["cov_selected"] * np.outer(factor, factor)
        write_covariance(f"{args.out_prefix}_{state}_{var}_statistical_covariance.csv",
                         covariance_rows(stats, factor, bins, dict(state=state, variable=var, **stat_metadata(args))))

        rows = []
        for i in range(len(centers)):
            row = {
                "proton_selection": proton_selection_metadata()["mode"],
                "state": state,
                "variable": var,
                "x_low": bins[i],
                "x_high": bins[i + 1],
                "x_center": centers[i],
                "bin_width": widths[i],
                "n_truth_visible": n_true[i],
                "n_reco_visible": n_reco[i],
                "efficiency": eff[i],
                "n_corrected_visible": n_corr[i],
                "n_corrected_visible_err": err_corr[i],
                "flux_integral": flux_int,
                "lumi_pb_inv": lumi_pb,
                "output_convention": args.output_convention,
                "br_corr": br_corr,
                "dsigma_dx_pb": dsigma[i],
                "dsigma_dx_pb_err": dsigma_err[i],
                "reco_count_to_dsigma_dx_pb": factor[i]/eff[i] if stats["valid"][i] else np.nan,
                "dsigma_dx_pb_interval_low": intervals["lower"][i]*factor[i],
                "dsigma_dx_pb_interval_high": intervals["upper"][i]*factor[i],
                "dsigma_dx_pb_err_low": intervals["error_low"][i]*factor[i],
                "dsigma_dx_pb_err_high": intervals["error_high"][i]*factor[i],
                "dsigma_dx_pb_count_interval_low": intervals["count_lower"][i]*factor[i],
                "dsigma_dx_pb_count_interval_high": intervals["count_upper"][i]*factor[i],
                "n_reco_count_interval_low": intervals["reco_lower"][i],
                "n_reco_count_interval_high": intervals["reco_upper"][i],
                "n_reco_interval_reference_count": intervals["reference_count"][i],
                "dsigma_dx_pb_err_mc": np.sqrt(stats["cov_mc"][i, i]) * factor[i],
                "dsigma_dx_pb_err_projected": np.sqrt(stats["cov_projected"][i, i]) * factor[i],
                "dsigma_dx_pb_err_correction_mc": np.sqrt(stats["cov_correction_mc"][i, i]) * factor[i],
                "efficiency_err_mc": np.sqrt(stats["efficiency_mc_cov"][i, i]),
                **stat_metadata(args),
                **branching_fraction_metadata(args),
                "closure_corr_over_truth": closure[i],
                "passes_min_eff": bool(good[i]),
                "raw_visible_xsec_fb": raw_visible_xsec_fb,
                "target_visible_xsec_fb": visible_xsec,
                "state_scale": state_scale,
            }
            rows.append(row)
            all_rows.append(row)

        per_var_results[var] = {
            "bins": bins,
            "centers": centers,
            "widths": widths,
            "y": dsigma,
            "yerr": dsigma_err,
            "covariance": covariance,
            "interval_low": intervals["lower"]*factor,
            "interval_high": intervals["upper"]*factor,
            "eff": eff,
            "good": good,
            "rows": rows,
        }

    metadata = {
        "proton_selection": proton_selection_metadata()["mode"],
        "state": state,
        "br_ee": br,
        "raw_visible_xsec_fb": raw_visible_xsec_fb,
        "target_visible_xsec_fb": visible_xsec,
        "state_scale": state_scale,
        **stat_metadata(args),
        **branching_fraction_metadata(args),
    }
    return all_rows, per_var_results, metadata


def write_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path}")


def main():
    argv = preprocess_negative_bin_args(sys.argv[1:])

    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--state", choices=list(STATES.keys()), default=None,
                   help="Backward-compatible single state option.")
    p.add_argument("--states", default=None,
                   help="Comma-separated states to measure and plot, e.g. 1S,2S,3S.")
    p.add_argument("--electron-energy", type=float, default=18.0)
    p.add_argument("--proton-energy", type=float, default=275.0)
    p.add_argument("--lumi-invfb", type=float, default=10.0)
    p.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    p.add_argument("--stat-mode", choices=["mc", "projected"], default="mc")
    p.add_argument("--projection-lumi-invfb", type=float, default=None)
    p.add_argument("--include-efficiency-mc-uncertainty", type=int, choices=[0, 1], default=0,
                   help="Add finite-MC correction covariance to projected data statistics; projected mode only.")

    # Backward-compatible single-state xsec and per-state options.
    p.add_argument("--visible-xsec-fb", type=float, default=None)
    p.add_argument("--visible-xsec-1s-fb", type=float, default=None)
    p.add_argument("--visible-xsec-2s-fb", type=float, default=None)
    p.add_argument("--visible-xsec-3s-fb", type=float, default=None)

    p.add_argument("--state-normalize", action="store_true", default=True)
    p.add_argument("--no-state-normalize", dest="state_normalize", action="store_false")
    p.add_argument("--output-convention", choices=["visible_decay", "inclusive_upsilon"], default="visible_decay")
    p.add_argument("--q2-max", type=float, default=0.001)
    p.add_argument("--q2-cut", type=float, default=None)
    p.add_argument("--w-min", type=float, default=20.0)
    p.add_argument("--w-max", type=float, default=140.0)

    p.add_argument("--variables", default="pT2,y")
    p.add_argument("--pt2-bins", default=VARIABLES["pT2"]["default_bins"])
    p.add_argument("--y-bins", default=VARIABLES["y"]["default_bins"])
    p.add_argument("--y-plot-min", type=float, default=None)
    p.add_argument("--y-plot-max", type=float, default=None)
    p.add_argument("--w-bins", default=VARIABLES["W"]["default_bins"])
    p.add_argument("--t-bins", default=VARIABLES["t"]["default_bins"])

    p.add_argument("--truth-mass-window", type=float, default=0.20)
    p.add_argument("--reco-mass-window-sigma", type=float, default=3.0)
    p.add_argument("--mass-min", type=float, default=None)
    p.add_argument("--mass-max", type=float, default=None)
    p.add_argument("--sigma-epz-center", type=float, default=None)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    p.add_argument("--min-eff", type=float, default=0.05)

    # ZEUS-style b-slope extraction from the finite-bin pT2 or |t| spectrum.
    p.add_argument("--fit-b-slope", action="store_true", default=True)
    p.add_argument("--no-fit-b-slope", dest="fit_b_slope", action="store_false")
    p.add_argument("--b-slope-variable", choices=["pT2", "t"], default="pT2")
    p.add_argument("--b-slope-bins", default="0,0.2,0.5,1.0,5.0",
                   help="Finite bins used for the b-slope fit, ZEUS-like by default.")
    p.add_argument("--b-slope-min", type=float, default=0.05)
    p.add_argument("--b-slope-max", type=float, default=20.0)
    p.add_argument("--b-slope-n-scan", type=int, default=4000)
    p.add_argument("--b-slope-zeus-mass-window", action="store_true", default=False,
                   help="For 1S use the ZEUS Section 6 signal mass window 9.33--9.66 GeV.")

    p.add_argument("--out-prefix", default="diff_gamma_p")
    p.add_argument("--campaign", default="")
    p.add_argument("--energy-label", default="")
    add_interval_arguments(p)
    add_branching_fraction_argument(p)
    args = configure_statistics(p.parse_args(argv))
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)

    if args.states is None:
        args.states = args.state if args.state is not None else "1S"
    states = parse_states(args.states)

    if args.q2_cut is None:
        args.q2_cut = args.q2_max

    variables = [v.strip() for v in args.variables.split(",") if v.strip()]
    for v in variables:
        if v not in VARIABLES:
            raise ValueError(f"Unknown variable {v}. Available: {list(VARIABLES)}")

    bins_by_var = {
        "pT2": parse_edges(args.pt2_bins),
        "y": parse_edges(args.y_bins),
        "W": parse_edges(args.w_bins),
        "t": parse_edges(args.t_bins),
    }
    bins_by_var = {v: bins_by_var[v] for v in variables}

    filename = os.path.expanduser(args.input)

    cols = [
        "weight",
        "true_mass", "true_W", "true_Q2",
        "reco_mass", "reco_W", "sigma_epz",
        "true_pT", "reco_pT",
        "true_y", "reco_y",
        "true_t",
    ]
    cols += proton_branches()
    with uproot.open(filename) as source:
        tree = source[args.tree]
        if "reco_t_eX" in tree.keys():
            cols.append("reco_t_eX")
        arr = {name: tree[name].array(library="np") for name in cols}
    arr["true_pT2"] = arr.pop("true_pT") ** 2
    arr["reco_pT2"] = arr.pop("reco_pT") ** 2
    arr["reco_t"] = arr.pop("reco_t_eX") if "reco_t_eX" in arr else arr["reco_pT2"]
    data = {k: np.asarray(v, dtype=float) for k, v in arr.items()}

    # Validate the full producer contract before any existing finite-value filter.
    proton_tag_mask(data)
    finite = np.ones(len(data["weight"]), dtype=bool)
    for name, v in data.items():
        if name not in proton_branches():
            finite &= np.isfinite(v)
    data = {k: v[finite] for k, v in data.items()}

    flux_int = integrated_flux(args.w_min, args.w_max, args.electron_energy, args.proton_energy, args.q2_max)
    lumi_pb = args.lumi_invfb * 1000.0

    all_rows = []
    metadata_rows = []
    all_results = {}

    for state in states:
        rows, results, meta = extract_for_state(data, args, state, bins_by_var, flux_int, lumi_pb)
        all_rows.extend(rows)
        meta_full = dict(meta)
        meta_full.update({
            "proton_selection": proton_selection_metadata(),
            "input": filename,
            "w_min": args.w_min,
            "w_max": args.w_max,
            "q2_max_flux": args.q2_max,
            "q2_cut_events": args.q2_cut,
            "flux_integral": flux_int,
            "output_convention": args.output_convention,
            **stat_metadata(args),
            **branching_fraction_metadata(args),
        })
        metadata_rows.append(meta_full)
        all_results[state] = results

        if args.fit_b_slope:
            b_summary, b_bins = extract_b_slope_for_state(data, args, state, flux_int, lumi_pb)
            all_results[state]["b_slope"] = b_summary
            all_results[state]["b_slope_bins"] = b_bins
            write_rows(f"{args.out_prefix}_{state}_b_slope_bins.csv", b_bins)

        # state-specific output files
        state_rows = rows
        write_rows(f"{args.out_prefix}_{state}_all.csv", state_rows)
        for var in variables:
            var_rows = [r for r in rows if r["variable"] == var]
            write_rows(f"{args.out_prefix}_{state}_{var}.csv", var_rows)

    write_rows(f"{args.out_prefix}_all_states_all.csv", all_rows)
    write_rows(f"{args.out_prefix}_metadata.csv", metadata_rows)

    if args.fit_b_slope:
        b_summary_rows = [all_results[st]["b_slope"] for st in states if "b_slope" in all_results[st]]
        write_rows(f"{args.out_prefix}_b_slope_summary.csv", b_summary_rows)

    # Plot pT2 and y in side-by-side style when available.
    plot_vars = [v for v in ["pT2", "y"] if v in variables]
    if plot_vars:
        fig, axes = plt.subplots(1, len(plot_vars), figsize=(7.6 * len(plot_vars), 6.0))
        if len(plot_vars) == 1:
            axes = [axes]

        for ax, var in zip(axes, plot_vars):
            for state in states:
                res = all_results[state][var]
                x = res["centers"]
                xerr = 0.5 * res["widths"]
                y = res["y"]
                yerr = res["yerr"]
                good = res["good"] & (y > 0)
                color = STATES[state]["color"]

                draw_interval_bars(ax, x[good], y[good], res["interval_low"][good], res["interval_high"][good], xerr=xerr[good],
                            fmt="o", color=color, label=rf"$\Upsilon({state})$", capsize=3)

                if var == "pT2" and args.fit_b_slope and "b_slope" in all_results[state]:
                    bs = all_results[state]["b_slope"]
                    B = bs["b_slope_GeVminus2"]
                    Blo = bs["b_slope_err_low"]
                    Bhi = bs["b_slope_err_high"]
                    sig_tot = bs["fit_integrated_sigma_pb"]
                    if np.isfinite(B) and np.isfinite(sig_tot):
                        bedges = parse_edges(args.b_slope_bins)
                        norm_int = math.exp(-B * bedges[0]) - math.exp(-B * bedges[-1])
                        xx = np.linspace(bedges[0], bedges[-1], 300)
                        yy = sig_tot * B * np.exp(-B * xx) / norm_int if norm_int > 0 else np.zeros_like(xx)
                        ax.plot(xx, yy, color=color, linewidth=1.4, alpha=0.75,
                                label=rf"$\Upsilon({state})$ binned fit, $b={B:.2f}^{{+{Bhi:.2f}}}_{{-{Blo:.2f}}}$")
                elif var == "pT2" and args.stat_mode == "mc" and np.sum(good) >= 2:
                    A, Aerr, B, Berr = weighted_exp_fit(x[good], y[good], yerr[good], res["covariance"][np.ix_(good, good)])
                    if np.isfinite(A) and np.isfinite(B):
                        xx = np.linspace(0.0, max(x[good]) * 1.15, 300)
                        yy = A * np.exp(-B * xx)
                        ax.plot(xx, yy, color=color, linewidth=1.4, alpha=0.75,
                                label=rf"$\Upsilon({state})$ exp. fit, $B={B:.2f}\pm{Berr:.2f}$")

            if var == "pT2":
                ax.set_xlabel(r"$p_{T,\Upsilon}^{2}$ [GeV$^2$]")
                ylabel_prefix = r"$d\sigma_{\gamma p}/dp_T^2$"
                if args.output_convention == "visible_decay":
                    ylabel_prefix = r"$B_{ee}\,d\sigma_{\gamma p}/dp_T^2$"
                ax.set_ylabel(ylabel_prefix + r" [pb/GeV$^2$]")
                ax.set_yscale("log")
            elif var == "y":
                ax.set_xlabel(r"$y_{\Upsilon}$")
                if args.y_plot_min is not None or args.y_plot_max is not None:
                    ax.set_xlim(
                        args.y_plot_min if args.y_plot_min is not None else bins_by_var["y"][0],
                        args.y_plot_max if args.y_plot_max is not None else bins_by_var["y"][-1],
                    )
                else:
                    ax.set_xlim(bins_by_var["y"][0], bins_by_var["y"][-1])
                ylabel_prefix = r"$d\sigma_{\gamma p}/dy$"
                if args.output_convention == "visible_decay":
                    ylabel_prefix = r"$B_{ee}\,d\sigma_{\gamma p}/dy$"
                ax.set_ylabel(ylabel_prefix + r" [pb]")
                ax.set_ylim(bottom=0.0)

            ax.legend(frameon=False, fontsize=9)
            ax.text(0.04, 0.94, "ePIC Performance", transform=ax.transAxes,
                    fontsize=12, fontweight="bold", ha="left", va="top")
            label = args.energy_label if args.energy_label else rf"e+p, {args.electron_energy:g}$\times${args.proton_energy:g} GeV"
            ax.text(0.04, 0.88, label, transform=ax.transAxes, fontsize=10.5, ha="left", va="top")
            ax.text(0.04, 0.82, rf"$W_{{\gamma p}}\in[{args.w_min:g},{args.w_max:g}]$ GeV, $Q^2<{args.q2_cut:g}$ GeV$^2$",
                    transform=ax.transAxes, fontsize=10.5, ha="left", va="top")
            ax.text(0.04, 0.76, rf"$\mathcal{{L}}_{{ep}}={args.lumi_invfb:g}$ fb$^{{-1}}$",
                    transform=ax.transAxes, fontsize=10.5, ha="left", va="top")

        fig.suptitle(r"$\gamma p\to\Upsilon(nS)p$ differential cross sections", fontsize=16)
        fig.text(0.5, 0.012, statistical_label(args),
                 ha="center", va="bottom", fontsize=9)
        fig.tight_layout(rect=[0, 0.04, 1, 0.94])
        out_png = f"{args.out_prefix}_pt2_y.png"
        out_pdf = f"{args.out_prefix}_pt2_y.pdf"
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_pdf)
        print(f"Wrote {out_png}")
        print(f"Wrote {out_pdf}")

    print("\n=== Metadata ===")
    for row in metadata_rows:
        print(row)


if __name__ == "__main__":
    main()
