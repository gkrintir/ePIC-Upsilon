"""Statistical models for fixed-response unfolding of weighted MC templates.

MC mode retains the historical clipped Gaussian sumw2 toys and central68%
half-width. Projected mode treats scaled reconstructed means as Poisson event
counts, reruns the unfolding for every toy, and reports the full covariance.
The response and prior stay fixed in both modes. Response-MC, backgrounds,
normalization systematics and other nuisance parameters are not included.
"""

from dataclasses import dataclass
from pathlib import Path
import sys
from statistics import NormalDist

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projected_intervals import (
    DEFAULT_COVERAGE, add_interval_arguments, draw_absolute_intervals, sampling_interval,
)


@dataclass
class ToyStatistics:
    error: np.ndarray
    covariance: np.ndarray
    standard_deviation: np.ndarray
    q16: np.ndarray
    q50: np.ndarray
    q84: np.ndarray
    mean: np.ndarray
    toys: np.ndarray
    interval_low: np.ndarray
    interval_high: np.ndarray


def add_statistical_arguments(parser):
    parser.add_argument("--stat-mode", choices=["mc", "projected"], default="mc",
                        help="mc: weighted-template MC precision; projected: Poisson data at the reconstructed level, with fixed MC response")
    parser.add_argument("--projection-lumi-invfb", type=float, default=None,
                        help="Target luminosity for projected mode; defaults to --lumi-invfb")
    add_interval_arguments(parser)
    parser.epilog = ("For unfolded spectra, --count-interval median_garwood selects central sampling bands "
                     "from complete reconstructed-Poisson toys, not Garwood confidence limits. "
                     "The nominal expectation can lie outside a sampling band; covariance is unchanged.")


def configure_statistics(args, parser):
    if not np.isfinite(args.tuple_lumi_invfb) or args.tuple_lumi_invfb <= 0:
        parser.error("--tuple-lumi-invfb must be finite and positive")
    if args.n_toys < 0:
        parser.error("--n-toys cannot be negative")
    if args.n_iterations < 1:
        parser.error("--n-iterations must be positive")
    if not 0 < getattr(args, "interval_coverage", DEFAULT_COVERAGE) < 1:
        parser.error("--interval-coverage must be strictly between 0 and 1")
    if args.projection_lumi_invfb is None:
        args.projection_lumi_invfb = args.lumi_invfb
    if not np.isfinite(args.projection_lumi_invfb) or args.projection_lumi_invfb <= 0:
        parser.error("--projection-lumi-invfb must be finite and positive")
    if args.stat_mode == "projected" and args.n_toys < 2:
        parser.error("Projected unfolding requires at least two toys for a covariance estimate")
    args.source_test_lumi_invfb = args.tuple_lumi_invfb * (1.0 - args.train_fraction)
    args.projection_scale = (args.projection_lumi_invfb / args.source_test_lumi_invfb
                             if args.stat_mode == "projected" else 1.0)
    # Legacy MC conversion used --lumi-invfb rather than changing stored weights.
    # Keep it numerically unchanged; projected mode uses the actual target L.
    args.yield_lumi_invfb = (args.projection_lumi_invfb if args.stat_mode == "projected"
                            else args.lumi_invfb * (1.0 - args.train_fraction))


def require_projectable_weights(weights, args):
    if args.stat_mode == "projected":
        weights = np.asarray(weights, dtype=float)
        if not np.all(np.isfinite(weights)) or np.any(weights < 0):
            raise ValueError("Projected Poisson statistics require finite nonnegative event weights; signed-weight samples are unsupported")


def check_response_support(response, expected_truth, expected_reco, args):
    """Do not present unidentifiable projected bins as zero with zero error."""
    response = np.asarray(response, dtype=float)
    missing_truth = np.flatnonzero((np.asarray(expected_truth) > 0) & (response.sum(axis=1) <= 0))
    missing_reco = np.flatnonzero((np.asarray(expected_reco) > 0) & (response.sum(axis=0) <= 0))
    args.unsupported_truth_bins = ",".join(map(str, missing_truth))
    args.unsupported_reco_bins = ",".join(map(str, missing_reco))
    args.zero_efficiency_truth_bins = ",".join(map(str, np.flatnonzero(response.sum(axis=1) <= 0)))
    if args.stat_mode == "projected":
        if not np.all(np.isfinite(response)) or np.any(response < 0):
            raise ValueError("Projected unfolding requires a finite nonnegative response")
        if len(missing_truth) or len(missing_reco):
            raise ValueError("Projected unfolding has unsupported bins: positive truth with zero efficiency "
                             f"[{args.unsupported_truth_bins}]; positive reco means without response support "
                             f"[{args.unsupported_reco_bins}]. Merge bins or improve the response MC.")


def sample_unfolded_statistics(means, mc_errors, unfold, n_toys, seed, mode, coverage=DEFAULT_COVERAGE):
    """Fluctuate reconstructed bins, then unfold each complete toy spectrum.

    `unfold` is a callback with a fixed training response and prior. Its return
    value is one unfolded spectrum; covariance therefore contains migrations
    and nonlinear iterative-unfolding correlations, conditional on that model.
    """
    means = np.asarray(means, dtype=float)
    mc_errors = np.asarray(mc_errors, dtype=float)
    central = np.asarray(unfold(means), dtype=float)
    if mode not in ("mc", "projected"):
        raise ValueError("Unknown statistical mode")
    if mode == "projected":
        if not np.all(np.isfinite(means)) or np.any(means < 0):
            raise ValueError("Projected reconstructed means must be finite and nonnegative")
        if n_toys < 2:
            raise ValueError("Projected covariance requires at least two toys")
    if not np.all(np.isfinite(mc_errors)) or np.any(mc_errors < 0):
        raise ValueError("MC template errors must be finite and nonnegative")
    if n_toys <= 0:
        zeros = np.zeros_like(central)
        return ToyStatistics(zeros, np.zeros((len(central), len(central))), zeros,
                             central.copy(), central.copy(), central.copy(), central.copy(),
                             np.empty((0, len(central))), central.copy(), central.copy())
    rng = np.random.default_rng(seed)
    toys = []
    for _ in range(n_toys):
        if mode == "projected":
            varied = rng.poisson(means).astype(float)
        else:
            varied = np.clip(rng.normal(means, mc_errors), 0.0, None)
        toys.append(np.asarray(unfold(varied), dtype=float))
    toys = np.vstack(toys)
    if not np.all(np.isfinite(toys)):
        raise ValueError("Unfolding produced nonfinite toy results")
    covariance = (np.atleast_2d(np.cov(toys, rowvar=False, ddof=1))
                  if n_toys > 1 else np.zeros((len(central), len(central))))
    covariance = 0.5 * (covariance + covariance.T)
    standard_deviation = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    q16, q50, q84 = np.quantile(toys, [0.16, 0.5, 0.84], axis=0)
    interval_low, interval_high = sampling_interval(toys, coverage=coverage, axis=0)
    error = standard_deviation if mode == "projected" else 0.5 * (q84 - q16)
    return ToyStatistics(error, covariance, standard_deviation, q16, q50, q84, toys.mean(axis=0), toys,
                         interval_low, interval_high)


def identifiable_intervals(stats, efficiency, mode):
    """Zero-efficiency bins have no inferred constraint, including empty bins."""
    low, high = stats.interval_low.copy(), stats.interval_high.copy()
    if mode == "projected":
        low[np.asarray(efficiency) <= 0] = np.nan
        high[np.asarray(efficiency) <= 0] = np.nan
    return low, high


def gaussian_interval_scale(coverage=DEFAULT_COVERAGE):
    return NormalDist().inv_cdf((1.0 + coverage) / 2.0)


def summary_interval_bounds(summary, value, error, lower, upper):
    mode = str(summary["stat_mode"].iloc[0]) if "stat_mode" in summary else "mc"
    selection = (str(summary["count_interval"].iloc[0]) if "count_interval" in summary
                 else "standard_deviation")
    if mode == "projected" and selection != "standard_deviation":
        lo, hi = summary[lower].to_numpy(float).copy(), summary[upper].to_numpy(float).copy()
    else:
        scale = (gaussian_interval_scale(float(summary["interval_coverage"].iloc[0]))
                 if mode == "projected" and "interval_coverage" in summary else 1.0)
        lo = (summary[value] - scale * summary[error]).to_numpy(float)
        hi = (summary[value] + scale * summary[error]).to_numpy(float)
    if mode == "projected" and "unfolding_identifiable" in summary:
        mask = ~summary["unfolding_identifiable"].to_numpy(bool)
        lo[mask], hi[mask] = np.nan, np.nan
    return lo, hi


def draw_summary_intervals(ax, summary, x, xerr, value, error, lower, upper, **kwargs):
    mode = str(summary["stat_mode"].iloc[0]) if "stat_mode" in summary else "mc"
    selection = (str(summary["count_interval"].iloc[0]) if "count_interval" in summary
                 else "standard_deviation")
    y = summary[value].to_numpy(float).copy()
    if mode == "projected" and "unfolding_identifiable" in summary:
        y[~summary["unfolding_identifiable"].to_numpy(bool)] = np.nan
    if mode == "projected" and selection != "standard_deviation":
        lo, hi = summary_interval_bounds(summary, value, error, lower, upper)
        return draw_absolute_intervals(ax, x, y, lo, hi, xerr=xerr, **kwargs)
    scale = (gaussian_interval_scale(float(summary["interval_coverage"].iloc[0]))
             if mode == "projected" and "interval_coverage" in summary else 1.0)
    return ax.errorbar(x, y, xerr=xerr, yerr=scale * summary[error], **kwargs)


def save_unfolded_toys(path, stats, bin_edges, args):
    """Persist correlated spectra for nonlinear downstream transformations."""
    np.savez_compressed(path, unfolded_visible=stats.toys,
                        bin_edges=np.asarray(bin_edges, dtype=float),
                        stat_mode=np.asarray(args.stat_mode),
                        yield_lumi_invfb=np.asarray(args.yield_lumi_invfb),
                        count_interval=np.asarray(getattr(args, "count_interval", "median_garwood")),
                        interval_coverage=np.asarray(getattr(args, "interval_coverage", DEFAULT_COVERAGE)))


def covariance_to_csv(path, covariance, labels):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(covariance, index=labels, columns=labels).to_csv(path, index_label="truth_bin")


def scale_covariance(covariance, factors):
    factors = np.asarray(factors, dtype=float)
    return np.asarray(covariance, dtype=float) * np.outer(factors, factors)


def statistical_metadata(args):
    projected = args.stat_mode == "projected"
    return {
        "stat_mode": args.stat_mode,
        "count_interval": getattr(args, "count_interval", "median_garwood"),
        "interval_coverage": getattr(args, "interval_coverage", DEFAULT_COVERAGE),
        "unfolded_interval_kind": "Pointwise central sampling quantiles of complete unfolded toys; not confidence intervals",
        "unfolded_display_interval": ("Gaussian symmetric band from covariance SD" if getattr(args, "count_interval", "median_garwood") == "standard_deviation"
                                       else "Pointwise central sampling quantiles; not Garwood confidence limits") if projected else "Legacy MC central68% half-width",
        "unfolded_interval_center": "Nominal unfolded expectation; may lie outside the central sampling interval",
        "unfolded_zero_mean_policy": "Zero assumed Poisson mean gives zero toys; this is not an observed-zero upper confidence limit",
        "normalized_shape_interval_conditioning": "Unit-area toy transformations exclude zero-total toys and are conditional on positive total yield",
        "projection_lumi_invfb": args.projection_lumi_invfb if projected else np.nan,
        "source_tuple_lumi_invfb": args.tuple_lumi_invfb,
        "source_test_lumi_invfb": args.source_test_lumi_invfb,
        "yield_lumi_invfb": args.yield_lumi_invfb,
        "projection_scale": args.projection_scale,
        "yield_scope": ("Expected full projected-luminosity yields extrapolated from held-out MC"
                        if projected else "Held-out test subset; yields are not scaled to the full tuple luminosity"),
        "uncertainty_convention": ("Reconstructed-bin Poisson toys propagated through fixed-response unfolding; main error=sqrt(covariance diagonal)"
                                   if projected else "Clipped Gaussian weighted-MC sumw2 toys; legacy main error=half central68% width"),
        "uncertainty_model": "projected_reco_poisson_fixed_response" if projected else "weighted_MC_sumw2",
        "covariance_valid": args.n_toys > 1,
        "covariance_toys": args.n_toys,
        "unsupported_truth_bins": getattr(args, "unsupported_truth_bins", ""),
        "unsupported_reco_bins": getattr(args, "unsupported_reco_bins", ""),
        "zero_efficiency_truth_bins": getattr(args, "zero_efficiency_truth_bins", ""),
        "zero_efficiency_policy": "Unconstrained truth bins; projected interval endpoints are NaN and markers are omitted",
        "response_uncertainty_included": False,
        "prior_uncertainty_included": False,
        "template_mc_uncertainty_combined": False,
        "template_mc_uncertainty_output": "Separate Gaussian-toy covariance from held-out reco sumw2; excludes training-response and prior uncertainty",
        "truth_reference_uncertainty": "Scaled held-out MC sumw2 precision, not projected-data Poisson uncertainty",
        "closure_error_convention": "Numerator-only error with fixed truth denominator; no reco-truth covariance or denominator uncertainty propagated",
        "excluded_uncertainties": "Response/prior MC statistics, backgrounds, luminosity, branching fraction, flux and other systematics",
    }


def mode_annotation(mode, count_interval="standard_deviation", coverage=DEFAULT_COVERAGE):
    if mode == "projected":
        if count_interval != "standard_deviation":
            return f"Signal-only; {100 * coverage:.2f}% central sampling band, not a confidence interval; fixed response"
        return f"Signal-only projected statistics; Gaussian {100 * coverage:.2f}% band; fixed response"
    return "Signal-only MC; statistical bands reflect MC precision"
