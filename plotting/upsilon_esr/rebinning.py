"""Exact adjacent-bin aggregation of the saved W/y estimators.

This is a display merge, not a refit with three truth parameters or a new
optimization. It preserves the source validation status. Sampling quantiles
are computed from sums of the original joint Poisson toy estimates.
"""
from copy import deepcopy
import json
import warnings

import numpy as np
import pandas as pd


def requested_edges(config, state, observable):
    return config.rebin_edges.get(state, {}).get(observable)


def merge_matrix(old_edges, new_edges):
    old = np.asarray(old_edges, float)
    if old.ndim != 1 or len(old) < 2 or not np.isfinite(old).all() or np.any(np.diff(old) <= 0):
        raise ValueError("Saved edges must be finite and strictly increasing.")
    if isinstance(new_edges, str):
        if new_edges != "keep_outer_bins":
            raise ValueError(f"Unknown display-rebin rule: {new_edges!r}")
        # Preserve the previous first/interior/last grouping using actual
        # saved boundaries. This is a display choice, not an optimization.
        if len(old) > 4:
            new = old[[0, 1, -2, -1]]
        else:
            new = old.copy()
            if len(old) < 4:
                warnings.warn(f"Only {len(old)-1} saved bins are available; keeping them unchanged. "
                              "Three bins would require splitting a saved bin.", RuntimeWarning)
    else:
        new = np.asarray(new_edges, float)
    if new.ndim != 1 or len(new) < 2 or not np.isfinite(new).all():
        raise ValueError("Rebin edges must be a finite one-dimensional list.")
    if np.any(np.diff(old) <= 0) or np.any(np.diff(new) <= 0):
        raise ValueError("Rebin edges must be strictly increasing.")
    indices = []
    for edge in new:
        matches = np.flatnonzero(np.isclose(old, edge, rtol=0, atol=1e-9))
        if len(matches) != 1:
            raise ValueError(f"Rebin edge {edge:g} is not a saved edge. Existing bins cannot be split. "
                             f"Saved edges: {old.tolist()}")
        indices.append(int(matches[0]))
    if indices[0] != 0 or indices[-1] != len(old)-1:
        raise ValueError("Rebinning must preserve the full saved range.")
    matrix = np.zeros((len(new)-1, len(old)-1))
    for i, (lo, hi) in enumerate(zip(indices[:-1], indices[1:])):
        matrix[i, lo:hi] = 1.
    return matrix, old[indices]


def maybe_rebin(config, state, observable, frame, metadata, entry):
    edges = requested_edges(config, state, observable)
    if edges is None:
        return frame, metadata
    if state != "1S" or observable not in ("w", "y"):
        raise ValueError("This display-rebin option supports 1S W and rapidity only.")
    old_edges = np.r_[frame.x_low.to_numpy(float), frame.x_high.iloc[-1]]
    group, new_edges = merge_matrix(old_edges, edges)
    if group.shape[0] == group.shape[1]:
        return frame, metadata
    from .rebin_toys import replay_joint_toys
    values, replay = replay_joint_toys(entry, cache_dir=config.output_dir / "rebin_cache")
    merged_toys = values @ group.T
    n = len(new_edges)-1
    settings = entry["inference"]["settings"]
    alpha = 1. - float(settings["coverage"])
    low, high = np.quantile(merged_toys, [alpha/2., 1.-alpha/2.], axis=0)
    original_scale = entry["cross_section_scaling"]
    scale = deepcopy(original_scale)
    base = 1./(1000. * float(scale["luminosity_invfb"]) * float(scale["branching_fraction_ee"]))
    scale["bin_edges"] = new_edges.tolist()
    if observable == "w":
        flux = group @ np.asarray(scale["bin_flux"], float)
        factors = base / flux
        scale["bin_flux"] = flux.tolist()
    else:
        factors = base / np.diff(new_edges)
        scale["alternate_factors"] = (factors / float(entry["flux_integral"])).tolist()
    scale["factors"] = factors.tolist()
    old_factors = np.asarray(original_scale["factors"], float)
    display_factors = (np.asarray(scale["alternate_factors"], float) if observable == "y" else factors)
    old_display_factors = (np.asarray(original_scale["alternate_factors"], float) if observable == "y" else old_factors)
    transform = display_factors[:, None] * group / old_display_factors[None, :]
    old_covariance = np.asarray(metadata["covariance"], float)
    covariance = transform @ old_covariance @ transform.T
    toy_covariance = np.atleast_2d(np.cov(merged_toys, rowvar=False, ddof=1))
    count_covariance = group @ np.asarray(entry["inference"]["coarse"]["covariance"], float) @ group.T
    if not np.allclose(toy_covariance, count_covariance, rtol=2e-5, atol=1e-7):
        raise ValueError("Merged toy covariance disagrees with the propagated saved covariance.")
    output = pd.DataFrame(dict(bin=np.arange(n), x_low=new_edges[:-1], x_high=new_edges[1:]))
    for key in ("state", "observable", "status", "overall_status", "stage_status", "failed_checks",
                "publication_ready", "quantity", "quantity_units", "x_units", "interval_interpretation",
                "proton_selection_mode", "proton_selection_contract",
                "alternate_quantity", "alternate_units", "display_quantity", "target_lumi_invfb"):
        if key in frame:
            distinct = frame[key].astype(str).unique()
            if len(distinct) != 1:
                raise ValueError(f"Cannot merge inconsistent {key} metadata.")
            output[key] = frame[key].iloc[0]
    for key in ("fitted_truth", "reference_truth", "toy_mean"):
        output[key] = group @ frame[key].to_numpy(float)
    output["toy_lower"], output["toy_upper"] = low, high
    output["sampling_sd"] = np.sqrt(np.maximum(np.diag(count_covariance), 0.))
    output["conversion_factor"] = factors
    if observable == "y":
        output["alternate_conversion_factor"] = display_factors
    for key in ("fitted_truth", "reference_truth", "toy_mean", "toy_lower", "toy_upper", "sampling_sd"):
        output[key + "_cross_section"] = output[key] * factors
        if observable == "y":
            output[key + "_alternate_cross_section"] = output[key] * display_factors
    output["mean"] = output.fitted_truth * display_factors
    output["lower"] = low * display_factors
    output["upper"] = high * display_factors
    output["mean_count"] = group @ frame.mean_count.to_numpy(float)
    output["x_center"] = (new_edges[:-1] + new_edges[1:]) / 2.
    output["bin_width"] = np.diff(new_edges)
    output["valid"] = True
    output["rebinned"] = True
    output["source_bins"] = [json.dumps(np.flatnonzero(row).tolist()) for row in group]
    meta = deepcopy(metadata)
    meta["cross_section_scaling"] = scale
    meta["covariance"] = covariance.tolist()
    meta["source_covariance_metadata"] = meta.pop("covariance_metadata")
    meta["covariance_metadata"] = {
        "overall_status": metadata["overall_status"], "publication_ready": False,
        "truth_edges": new_edges.tolist(), "conversion_factors": factors.tolist(),
        "interpretation": "Saved full sampling covariance propagated by H C H^T; not an interval",
        "truth_yield_covariance": count_covariance.tolist(),
        "primary_cross_section_covariance": (count_covariance*factors[:,None]*factors[None,:]).tolist()}
    meta["rebinning"] = {
        "edge_selection": edges if isinstance(edges, str) else "explicit_saved_edges",
        "kind": "Adjacent-bin aggregation of original estimator; no three-parameter refit",
        "source_edges": old_edges.tolist(), "display_edges": new_edges.tolist(),
        "yield_merge_matrix": group.tolist(), "cross_section_merge_matrix": transform.tolist(),
        "uncertainties": "Central quantiles of summed joint Poisson toy estimates; not sums of marginal endpoints",
        "replay": replay, "validation_status_changed": False,
        "additional_validation": "No new binning optimization, coverage or stress validation is claimed"}
    return output, meta
