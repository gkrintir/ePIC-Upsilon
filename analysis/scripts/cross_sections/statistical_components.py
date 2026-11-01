"""Explicit MC and signal-only projected statistics for same-sample correction.

MC is Poissonized weighted-event sumw2 with the supplied normalization held fixed.
Truth and reconstruction contributions from the same event are correlated. For
R/(R/T), those correlations cancel R exactly, leaving the covariance of T.
Projected data are independent Poisson reconstructed counts; finite-MC correction
uncertainty is separate and is combined only on explicit request.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projected_intervals import DEFAULT_COVERAGE, add_interval_arguments, expected_count_intervals, draw_absolute_intervals
import numpy as np


def configure_statistics(args):
    args.count_interval = getattr(args, "count_interval", "median_garwood")
    args.interval_coverage = getattr(args, "interval_coverage", DEFAULT_COVERAGE)
    if not 0 < args.interval_coverage < 1:
        raise ValueError("--interval-coverage must be between zero and one")
    requested_projection = args.projection_lumi_invfb
    source = args.tuple_lumi_invfb
    if args.stat_mode == "mc":
        if not np.isclose(args.lumi_invfb, source, rtol=1e-12, atol=0):
            raise ValueError("MC mode requires --lumi-invfb to match --tuple-lumi-invfb; use projected mode to change luminosity")
        target = args.lumi_invfb
    else:
        target = args.lumi_invfb if requested_projection is None else requested_projection
    if not np.isfinite(target) or target <= 0 or not np.isfinite(source) or source <= 0:
        raise ValueError("Target and tuple luminosities must be finite and positive")
    include = bool(args.include_efficiency_mc_uncertainty) or bool(getattr(args, "include_eff_mc_unc", False))
    if include and args.stat_mode != "projected":
        raise ValueError("Efficiency MC uncertainty can be added only in projected mode; adding it to same-sample MC double-counts correlated fluctuations")
    args.requested_projection_lumi_invfb = requested_projection
    args.projection_lumi_invfb = float(target)
    args.lumi_invfb = float(target)
    args.include_efficiency_mc_uncertainty = int(include)
    args.stat_lumi_scale = float(target / source)
    return args


def metadata(args):
    return {
        "stat_mode": args.stat_mode,
        "count_interval": args.count_interval if args.stat_mode == "projected" else "mc_standard_deviation",
        "interval_coverage": args.interval_coverage,
        "interval_interpretation": (
            "Ensemble-median expected Garwood confidence limits; mean central value unchanged; not an interval for a fractional observation"
            if args.stat_mode == "projected" and args.count_interval == "median_garwood" else
            "Standard-deviation bars; not a confidence interval"),
        "interval_combination": ("quadrature_with_correction_mc_approximation" if args.include_efficiency_mc_uncertainty else "count_only")
            if args.stat_mode == "projected" else "mc_covariance_standard_deviation",
        "projection_lumi_invfb": args.projection_lumi_invfb,
        "tuple_lumi_invfb": args.tuple_lumi_invfb,
        "yield_lumi_invfb": args.lumi_invfb,
        "requested_projection_lumi_invfb": args.requested_projection_lumi_invfb,
        "yield_luminosity_scale": args.stat_lumi_scale,
        "include_efficiency_mc_uncertainty": bool(args.include_efficiency_mc_uncertainty),
        "uncertainty_model": ("same_sample_corrected_truth_sumw2" if args.stat_mode == "mc" else
                              "projected_reco_poisson_plus_correction_mc" if args.include_efficiency_mc_uncertainty else "projected_reco_poisson"),
        "mc_covariance_convention": "Poissonized weighted events; supplied normalization fixed; true/reco same-event covariance retained",
        "projection_assumption": "Independent signal-only Poisson reconstructed counts; no background, systematic, luminosity or branching-ratio uncertainty",
        "legacy_mc_change": "Corrected same-MC yield is R/(R/T)=T; its error is truth sumw2, replacing reco sumw2/efficiency with efficiency held fixed",
    }


def statistical_label(args):
    if args.stat_mode == "mc":
        return "Signal-only MC; same-event MC statistical covariance"
    suffix = " + correction-MC quadrature approximation" if args.include_efficiency_mc_uncertainty else "; MC correction held fixed"
    interval = (f"median expected Garwood limits ({100*args.interval_coverage:.2f}%)"
                if args.count_interval == "median_garwood" else "reconstructed-count standard deviations")
    return f"Signal-only projection at {args.projection_lumi_invfb:g} fb^-1; {interval}" + suffix


def _indices(values, edges):
    i = np.searchsorted(edges, values, side="right") - 1
    i = np.where(values == edges[-1], len(edges) - 2, i)
    return i


def _masked_covariance(cov, valid):
    result = np.asarray(cov, dtype=float).copy()
    result[~valid, :] = np.nan
    result[:, ~valid] = np.nan
    return result


def binned_components(true_values, reco_values, weights, truth_mask, reco_mask, edges, args):
    """Build yield and full correction covariance from aligned event arrays."""
    weights = np.asarray(weights, float) * args.stat_lumi_scale
    if not np.isfinite(weights).all():
        raise ValueError("Nonfinite input weights are not supported")
    if args.stat_mode == "projected" and np.any(weights < 0):
        raise ValueError("Projected Poisson statistics require nonnegative event weights")
    edges = np.asarray(edges, float)
    if len(edges) < 2 or not np.isfinite(edges).all() or not np.all(np.diff(edges) > 0):
        raise ValueError("Histogram edges must be finite and strictly increasing")
    ti, ri = _indices(np.asarray(true_values), edges), _indices(np.asarray(reco_values), edges)
    n = len(edges) - 1
    tm = np.asarray(truth_mask, bool) & np.isfinite(true_values) & (ti >= 0) & (ti < n)
    rm = np.asarray(reco_mask, bool) & np.isfinite(reco_values) & (ri >= 0) & (ri < n)
    t = np.bincount(ti[tm], weights=weights[tm], minlength=n)
    r = np.bincount(ri[rm], weights=weights[rm], minlength=n)
    vt = np.bincount(ti[tm], weights=weights[tm]**2, minlength=n)
    vr = np.bincount(ri[rm], weights=weights[rm]**2, minlength=n)
    both = tm & rm
    ctr = np.bincount(ti[both]*n + ri[both], weights=weights[both]**2, minlength=n*n).reshape(n, n)
    efficiency = np.divide(r, t, out=np.full(n, np.nan), where=t > 0)
    valid = (t > 0) & (r > 0) & np.isfinite(efficiency) & (efficiency > 0)
    d_r = np.divide(1.0, t, out=np.zeros(n), where=t > 0)
    d_t = np.divide(-r, t*t, out=np.zeros(n), where=t > 0)
    ce = (np.diag(d_r*d_r*vr + d_t*d_t*vt)
          + (d_t[:, None] * ctr) * d_r[None, :]
          + (d_r[:, None] * ctr.T) * d_t[None, :])
    ce = (ce + ce.T) / 2
    # Remove cancellation noise on the diagonal without treating migrating bins
    # as binomial probabilities (R/T may exceed one).
    np.fill_diagonal(ce, np.maximum(np.diag(ce), 0))
    corrected = np.divide(r, efficiency, out=np.full(n, np.nan), where=valid)
    projected_var = np.divide(r, efficiency**2, out=np.full(n, np.nan), where=valid)
    derivative = np.divide(-r, efficiency**2, out=np.zeros(n), where=valid)
    calibration = derivative[:, None] * ce * derivative[None, :]
    mc_cov = _masked_covariance(np.diag(vt), valid)
    projected_cov = _masked_covariance(np.diag(np.nan_to_num(projected_var)), valid)
    calibration = _masked_covariance(calibration, valid)
    if args.stat_mode == "mc":
        selected = mc_cov.copy()
    else:
        selected = projected_cov.copy()
        if args.include_efficiency_mc_uncertainty:
            selected += calibration
    return {
        "truth": t, "reco": r, "truth_mc_var": vt, "reco_mc_var": vr,
        "truth_reco_mc_cov": ctr, "efficiency": efficiency,
        "efficiency_mc_cov": _masked_covariance(ce, t > 0),
        "valid": valid, "corrected": corrected,
        "cov_mc": mc_cov, "cov_projected": projected_cov,
        "cov_correction_mc": calibration, "cov_selected": selected,
        "error": np.sqrt(np.maximum(np.diag(selected), 0)),
    }


def binned_intervals(stats, args):
    """Propagate expected reconstructed-count endpoints through fixed correction.

    Covariances remain standard deviations. Optional calibration errors are
    combined in quadrature with the two count-error lengths, an approximation
    with no exact combined-coverage claim. Pure count endpoints are also retained.
    """
    center = stats["corrected"]
    if args.stat_mode == "mc":
        lower, upper = center-stats["error"], center+stats["error"]
        return dict(lower=lower, upper=upper, count_lower=lower.copy(), count_upper=upper.copy(),
                    error_low=stats["error"], error_high=stats["error"],
                    reco_lower=stats["reco"]-np.sqrt(stats["reco_mc_var"]),
                    reco_upper=stats["reco"]+np.sqrt(stats["reco_mc_var"]),
                    reference_count=np.full_like(center, np.nan))
    count = expected_count_intervals(stats["reco"], coverage=args.interval_coverage, method=args.count_interval)
    low = np.divide(count["lower"], stats["efficiency"], out=np.full_like(center, np.nan), where=stats["valid"])
    high = np.divide(count["upper"], stats["efficiency"], out=np.full_like(center, np.nan), where=stats["valid"])
    lower, upper = low.copy(), high.copy()
    if args.include_efficiency_mc_uncertainty:
        calibration_variance = np.diag(stats["cov_correction_mc"])
        lower = center-np.sqrt((center-low)**2+calibration_variance)
        upper = center+np.sqrt((high-center)**2+calibration_variance)
    return dict(lower=lower, upper=upper, count_lower=low, count_upper=high,
                error_low=center-lower, error_high=upper-center,
                reco_lower=count["lower"], reco_upper=count["upper"], reference_count=count["reference_count"])


draw_interval_bars = draw_absolute_intervals


def normalized_shape(values, covariance):
    """Normalize only estimable bins and retain normalization-induced covariance."""
    values, covariance = np.asarray(values, float), np.asarray(covariance, float)
    valid = np.isfinite(values) & (values >= 0) & np.isfinite(np.diag(covariance))
    fraction = np.full(len(values), np.nan)
    result = np.full(covariance.shape, np.nan)
    total = float(values[valid].sum())
    if total <= 0:
        return fraction, result
    y = values[valid]
    f = y / total
    jacobian = (np.eye(len(y)) - f[:, None]) / total
    cov = covariance[np.ix_(valid, valid)]
    fraction[valid] = f
    result[np.ix_(valid, valid)] = jacobian @ cov @ jacobian.T
    return fraction, result


def covariance_rows(stats, cross_section_factor, edges, identifiers):
    factor = np.broadcast_to(np.asarray(cross_section_factor, float), stats["corrected"].shape)
    rows = []
    for component in ["mc", "projected", "correction_mc", "selected"]:
        cov = stats["cov_" + component]
        _, shape_cov = normalized_shape(stats["corrected"], cov)
        for i in range(len(factor)):
            for j in range(len(factor)):
                rows.append(dict(identifiers, component=component, bin_i=i, bin_j=j,
                                 bin_i_low=edges[i], bin_i_high=edges[i+1],
                                 bin_j_low=edges[j], bin_j_high=edges[j+1],
                                 covariance_corrected_visible=cov[i, j],
                                 covariance_cross_section=cov[i, j]*factor[i]*factor[j],
                                 covariance_shape_fraction=shape_cov[i, j],
                                 covariance_efficiency_mc=stats["efficiency_mc_cov"][i, j],
                                 covariance_truth_reco_mc=stats["truth_reco_mc_cov"][i, j]))
    return rows


def write_covariance(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
