"""Training-only bin selection and migration responses for projected spectra.

The response includes failed reconstruction. Reconstructed events whose true
variable lies outside the fit interval are a separate feed-in template. No
observed/held-out counts are consulted when choosing edges.
"""

from itertools import combinations
import math

import numpy as np


DEFAULTS = {
    "fit_range": [0.0, 1.0],
    "candidate_edges": [0.0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .6, .7, .8, 1.0],
    "reco_control_max": 5.0,
    "reco_control_min": None,
    "fine_truth_step": .02,
    "min_bins": 3,
    "max_bins": 5,
    "min_expected_reco": 5.0,
    "min_purity": .4,
    "min_selected_stability": .4,
    "min_effective_mc": 100.0,
    "max_condition_number": 20.0,
    "max_candidates": 100000,
}

OBSERVABLES = {
    "pt2": ("truth_pt2", "reco_pt2", "true_pT**2", "reco_pT**2", "GeV^2", "GeV^-2"),
    "abs_t_via_pt2": ("truth_abs_t", "reco_pt2", "true_t (positive convention)", "reco_pT**2", "GeV^2", "GeV^-2"),
    "w": ("truth_w", "reco_w", "true_W", "reco_W", "GeV", "dimensionless"),
    "y": ("truth_y", "reco_y", "true_y (Upsilon rapidity)", "reco_y (Upsilon rapidity)", "dimensionless", "dimensionless"),
}


def settings(config, observable):
    """Resolve documented global, optimization, then observable overrides."""
    if observable not in OBSERVABLES:
        raise ValueError("Unknown observable: " + str(observable))
    result = dict(DEFAULTS)
    result.update({k: config[k] for k in DEFAULTS if k in config})
    result.update(config.get("optimization", {}))
    observations = config.get("observables", {})
    if isinstance(observations, dict):
        result.update(observations.get(observable, {}))
    lo, hi = map(float, result["fit_range"])
    if not np.isfinite([lo, hi]).all() or (observable != "y" and lo < 0) or hi <= lo:
        raise ValueError("fit_range must be finite and increasing; only Upsilon rapidity may be negative")
    if result["reco_control_max"] is not None and (
            not np.isfinite(result["reco_control_max"]) or result["reco_control_max"] <= hi):
        raise ValueError("reco_control_max must exceed the fit-range upper edge")
    if result["reco_control_min"] is not None and (
            not np.isfinite(result["reco_control_min"]) or result["reco_control_min"] >= lo):
        raise ValueError("reco_control_min must be below the fit-range lower edge")
    if observable != "y" and result["reco_control_min"] is not None and result["reco_control_min"] < 0:
        raise ValueError("Only rapidity may have a negative reconstructed control boundary")
    for key in ("min_bins", "max_bins", "max_candidates"):
        value = result[key]
        if isinstance(value, bool) or int(value) != value or value < 1:
            raise ValueError(key + " must be a positive integer")
        result[key] = int(value)
    if observable in ("pt2", "abs_t_via_pt2") and result["min_bins"] < 3:
        raise ValueError("Slope bin optimization requires at least three bins")
    if observable in ("w", "y") and result["min_bins"] < 2:
        raise ValueError("Spectrum bin optimization requires at least two bins")
    if result["max_bins"] < result["min_bins"]:
        raise ValueError("max_bins must be at least min_bins")
    for key in ("min_expected_reco", "min_effective_mc", "fine_truth_step", "max_condition_number"):
        if not np.isfinite(result[key]) or result[key] <= 0:
            raise ValueError(key + " must be positive and finite")
    for key in ("min_purity", "min_selected_stability"):
        if not np.isfinite(result[key]) or not 0 <= result[key] <= 1:
            raise ValueError(key + " must be between zero and one")
    result["fit_range"] = [lo, hi]
    return result


def reco_edges_for_fit(edges, cfg):
    """Attach optional lower/upper control intervals without assuming either."""
    parts = []
    if cfg.get("reco_control_min") is not None:
        parts.append(np.array([float(cfg["reco_control_min"])]))
    parts.append(np.asarray(edges, dtype=float))
    if cfg.get("reco_control_max") is not None:
        parts.append(np.array([float(cfg["reco_control_max"])]))
    return _edges(np.concatenate(parts))


def fit_reco_indices(reco_edges, fit_range):
    """Indices of reconstructed bins inside the extraction interval."""
    reco_edges = np.asarray(reco_edges, dtype=float)
    lo, hi = fit_range
    return np.flatnonzero((reco_edges[:-1] >= lo - 1e-12)
                          & (reco_edges[1:] <= hi + 1e-12))


def _edges(values, fit_range=None):
    edges = np.asarray(values, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError("Edges must be a finite, strictly increasing vector")
    if fit_range is not None and not np.allclose(edges[[0, -1]], fit_range, rtol=0, atol=1e-12):
        raise ValueError("Edges must span the entire configured fit_range")
    return edges


def _arrays(sample, observable, split, tilt, fit_range, shape_mode="combined"):
    truth_key, reco_key, truth_label, reco_label, units, tilt_units = OBSERVABLES[observable]
    x = np.asarray(sample[truth_key], dtype=float)
    y = np.asarray(sample[reco_key], dtype=float)
    w = np.asarray(sample["weight"], dtype=float)
    truth = np.asarray(sample["truth_selected"])
    reco = np.asarray(sample.get("reco_selected_" + observable, sample["reco_selected"]))
    train = np.asarray(sample["is_train"])
    arrays = (x, y, w, truth, reco, train)
    if any(a.ndim != 1 for a in arrays) or len({len(a) for a in arrays}) != 1:
        raise ValueError("Sample arrays must be equally sized one-dimensional arrays")
    if not len(w):
        raise ValueError("Cannot build a response from an empty sample")
    if any(a.dtype.kind != "b" for a in (truth, reco, train)):
        raise ValueError("Selection and split arrays must be boolean")
    if not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError("Source projection weights must be finite and positive")
    invalid_truth = ~np.isfinite(x) | ((x == -999.0) if observable == "y" else (x < 0))
    if np.any(truth & invalid_truth):
        raise ValueError("Selected truth values must be finite and valid for " + observable)
    if np.any(reco & ~truth):
        raise ValueError("Current differential contract requires reco_selected within truth_selected")
    invalid_reco = ~np.isfinite(y) | ((y == -999.0) if observable == "y" else (y < 0))
    invalid_reco_selected = int(np.sum(reco & invalid_reco))
    if observable == "y":
        # Rapidity has an independent producer validity requirement. It must
        # neither reject negative physical values nor affect other observables.
        reco = reco & ~invalid_reco
    elif np.any(reco & invalid_reco):
        raise ValueError("Selected reco values must be finite and nonnegative")
    fraction = float(sample["train_fraction"])
    if not np.isfinite(fraction) or not 0 < fraction < 1:
        raise ValueError("train_fraction must be strictly between zero and one")
    if split == "train":
        chosen, exposure_fraction = train, fraction
    elif split in ("test", "validation", "holdout"):
        chosen, exposure_fraction = ~train, 1.0 - fraction
    elif split == "all":
        chosen, exposure_fraction = np.ones(len(w), dtype=bool), 1.0
    else:
        raise ValueError("split must be train, test, validation, holdout, or all")
    if not np.isfinite(tilt):
        raise ValueError("shape_tilt must be finite")
    if shape_mode not in ("combined", "shape_only"):
        raise ValueError("shape_mode must be 'shape_only' or 'combined'")
    truth = truth & chosen
    reco = reco & chosen
    weights = w / exposure_fraction
    before = float(weights[truth].sum())
    in_range = truth & (x >= fit_range[0]) & (x <= fit_range[1])
    in_range_before = float(weights[in_range].sum())
    outside_before = float(weights[truth & ~in_range].sum())
    normalization = 1.0
    if tilt != 0:
        target = in_range if shape_mode == "shape_only" else truth
        target_before = in_range_before if shape_mode == "shape_only" else before
        if not target.any() or target_before <= 0:
            raise ValueError("Cannot apply a shape tilt to an empty truth selection")
        # Combined mode deliberately changes the outside-range yield through
        # normalization of the full selected truth population. Shape-only mode
        # modifies and normalizes in-range events only, leaving tails unchanged.
        tilt_coordinate = np.clip(x[target], *fit_range)
        if observable in ("w", "y"):
            tilt_coordinate = (tilt_coordinate - fit_range[0]) / (fit_range[1] - fit_range[0])
        log_multiplier = -float(tilt) * tilt_coordinate
        offset = float(np.max(log_multiplier))
        multiplier = np.exp(log_multiplier - offset)
        denominator = float(np.sum(weights[target] * multiplier))
        if not np.isfinite(denominator) or denominator <= 0:
            raise ValueError("Shape tilt has no finite weighted support")
        normalization = target_before / denominator
        weights = weights.copy()
        weights[target] *= multiplier * normalization
    metadata = {
        "observable": observable,
        "truth_variable": truth_label,
        "reco_variable": reco_label,
        "observable_units": units,
        "split": split,
        "split_exposure_fraction": exposure_fraction,
        "source_lumi_invfb": sample.get("source_lumi"),
        "target_lumi_invfb": sample.get("target_lumi"),
        "shape_tilt": float(tilt),
        "shape_tilt_units": tilt_units,
        "shape_mode": shape_mode,
        "shape_tilt_formula": ("exp(-shape_tilt * clip((truth_variable-fit_min)/(fit_max-fit_min), 0, 1)), followed by normalization preservation"
            if observable in ("w", "y") else "exp(-shape_tilt * clip(truth_variable, fit_min, fit_max)), followed by normalization preservation"),
        "shape_tilt_fit_range": list(map(float, fit_range)),
        "shape_tilt_normalization_scope": ("In-range selected truth only; outside-range truth and feed-in unchanged"
            if shape_mode == "shape_only" else "All selected truth events in this split, including those outside the variable fit range"),
        "truth_selected_expected_before_tilt": before,
        "truth_selected_expected_after_tilt": float(weights[truth].sum()),
        "in_range_truth_expected_before_tilt": in_range_before,
        "in_range_truth_expected_after_tilt": float(weights[in_range].sum()),
        "outside_range_truth_expected_before_tilt": outside_before,
        "outside_range_truth_expected_after_tilt": float(weights[truth & ~in_range].sum()),
        "invalid_reco_observable_selected_rows": invalid_reco_selected,
        "raw_categories_valid_for_multinomial_response_bootstrap": bool(tilt == 0 and np.allclose(w, w[0], rtol=1e-8, atol=0)),
        "selection_scope": "Current differential truth mass, W and Q2 fiducial selection; variable-range feed-in is included",
        "proton_selection": sample.get("metadata", {}).get("proton_selection", {"mode": "unrecorded"}),
        "proton_selection_cutflow_full_sample": sample.get("metadata", {}).get("proton_selection_cutflow", {}),
    }
    if observable in ("pt2", "abs_t_via_pt2"):
        metadata["shape_tilt_GeVminus2"] = float(tilt)
    return x, y, weights, truth, reco, metadata


def _bin_indices(values, edges, fit_hi=None):
    """NumPy-style bins, except the fit upper boundary stays in the fit bin."""
    result = np.searchsorted(edges, values, side="right") - 1
    result[values == edges[-1]] = len(edges) - 2
    if fit_hi is not None:
        match = np.flatnonzero(np.isclose(edges, fit_hi, rtol=0, atol=1e-12))
        if len(match) != 1:
            raise ValueError("fit boundary missing or duplicated in reco edges")
        result[values == fit_hi] = int(match[0]) - 1
    result[(~np.isfinite(values)) | (values < edges[0]) | (values > edges[-1])] = -1
    return result


def _hist(index, weights, count):
    selected = (index >= 0) & (index < count)
    return np.bincount(index[selected], weights=weights[selected], minlength=count).astype(float)


def _matrix(reco_index, truth_index, weights, n_reco, n_truth):
    selected = (reco_index >= 0) & (reco_index < n_reco) & (truth_index >= 0) & (truth_index < n_truth)
    flat = reco_index[selected] * n_truth + truth_index[selected]
    return np.bincount(flat, weights=weights[selected], minlength=n_reco * n_truth).reshape(n_reco, n_truth).astype(float)


def _histograms(sample, observable, truth_edges, reco_edges, split, tilt, shape_mode="combined"):
    x, y, w, truth, reco, metadata = _arrays(sample, observable, split, tilt, truth_edges[[0, -1]], shape_mode)
    nt, nr = len(truth_edges) - 1, len(reco_edges) - 1
    ti = _bin_indices(x, truth_edges)
    ri = _bin_indices(y, reco_edges, truth_edges[-1])
    ti[~truth] = -1
    ri[~reco] = -1
    truth_count = _hist(ti, np.ones(len(w)), nt)
    truth_sum = _hist(ti, w, nt)
    truth_sumw2 = _hist(ti, w * w, nt)
    matched_sum = _matrix(ri, ti, w, nr, nt)
    matched_sumw2 = _matrix(ri, ti, w * w, nr, nt)
    matched_count = _matrix(ri, ti, np.ones(len(w)), nr, nt)
    reco_sum = _hist(ri, w, nr)
    outside_truth = (ti < 0) & reco
    feed_in = _hist(np.where(outside_truth, ri, -1), w, nr)
    selected_ti = np.where(reco, ti, -1)
    selected_truth_sum = _hist(selected_ti, w, nt)
    selected_truth_sumw2 = _hist(selected_ti, w * w, nt)
    # Invalid/nonfinite reconstructed entries have already become inefficiency.
    # All other selected reco events outside the modeled range are also misses.
    missed_count = truth_count - matched_count.sum(axis=0)
    if np.any(missed_count < -1e-9):
        raise RuntimeError("Response categories exceed their truth denominator")
    metadata.update({
        "reco_boundary_convention": "Left-closed bins; fit upper edge belongs to final fit bin; control interval is open at that edge",
        "truth_boundary_convention": "Left-closed bins, final upper edge included",
        "reco_expected_in_model_range": float(reco_sum.sum()),
        "reco_expected_outside_model_range": float(w[reco & (ri < 0)].sum()),
        "feed_in_expected": float(feed_in.sum()),
        "feed_in_interpretation": "Fixed known MC expectation from truth-variable values outside fit range within the same truth mass/W/Q2 fiducial selection",
        "generated_truth_expected_in_fit_range": float(truth_sum.sum()),
        "response_orientation": "Rows reconstructed bins, columns truth bins",
        "raw_category_orientation": "Rows reconstructed bins followed by missed category, columns truth bins",
        "fit_reco_indices": fit_reco_indices(reco_edges, truth_edges[[0, -1]]).tolist(),
        "control_reco_indices": np.setdiff1d(np.arange(nr), fit_reco_indices(reco_edges, truth_edges[[0, -1]])).tolist(),
    })
    return {
        "truth_edges": np.asarray(truth_edges), "reco_edges": np.asarray(reco_edges),
        "truth_expected": truth_sum, "truth_sumw2": truth_sumw2,
        "truth_raw_count": truth_count, "matched_expected": matched_sum,
        "matched_sumw2": matched_sumw2,
        "reco_expected": reco_sum, "feed_in": feed_in,
        "selected_truth_expected": selected_truth_sum,
        "selected_truth_sumw2": selected_truth_sumw2,
        "raw_categories": np.vstack((matched_count, missed_count)).astype(np.int64),
        "metadata": metadata,
    }


def _divide(num, den):
    num, den = np.broadcast_arrays(num, den)
    return np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0)


def _finish(hist, coarse):
    truth = hist["truth_expected"]
    matched = hist["matched_expected"]
    response = _divide(matched, truth[np.newaxis, :])
    signal = matched.sum(axis=1)
    result = dict(hist, response=response, signal_reco_expected=signal)
    singular = np.linalg.svd(response, compute_uv=False)
    tolerance = (max(response.shape) * np.finfo(float).eps * singular[0]) if len(singular) else 0
    rank = int(np.sum(singular > tolerance))
    condition = float(singular[0] / singular[-1]) if rank == response.shape[1] and singular[-1] > 0 else float("inf")
    result["response_rank"] = rank
    result["condition_number"] = condition
    result["metadata"] = dict(hist["metadata"],
        missing_truth_columns=np.flatnonzero(truth <= 0).tolist(),
        all_missed_truth_columns=np.flatnonzero((truth > 0) & (matched.sum(axis=0) <= 0)).tolist(),
        min_raw_generated=int(hist["truth_raw_count"].min()) if len(truth) else 0,
        raw_generated_per_truth_cell=hist["truth_raw_count"].astype(int).tolist(),
        response_finite_mc_bootstrap_convention="Per-truth-cell multinomial resampling including missed events, fixed generated count; allowed only for uniform un-tilted source weights",
    )
    result["status"] = "MISSING_TRUTH_SUPPORT" if np.any(truth <= 0) else "READY"
    if coarse:
        n = len(truth)
        fit_indices = np.asarray(hist["metadata"]["fit_reco_indices"], dtype=int)
        if len(fit_indices) != n:
            raise ValueError("Coarse response requires matching truth and reconstructed fit bins")
        diagonal = np.diag(matched[fit_indices])
        generated_mc = _divide(truth * truth, hist["truth_sumw2"])
        selected_mc = _divide(hist["selected_truth_expected"] ** 2, hist["selected_truth_sumw2"])
        result["metrics"] = {
            "purity": _divide(diagonal, hist["reco_expected"][fit_indices]),
            "selected_stability": _divide(diagonal, hist["selected_truth_expected"]),
            "diagonal_acceptance": _divide(diagonal, truth),
            "efficiency_selected": _divide(hist["selected_truth_expected"], truth),
            "efficiency_modeled": _divide(matched.sum(axis=0), truth),
            "truth_effective_mc": generated_mc,
            "selected_effective_mc": selected_mc,
            "expected_reco_signal": signal[fit_indices],
            "feed_in_fraction": _divide(hist["feed_in"][fit_indices], hist["reco_expected"][fit_indices]),
        }
        result["metadata"]["metric_definitions"] = {
            "purity": "Diagonal matched count / all selected reco in that bin, including outside-fit truth feed-in",
            "selected_stability": "Diagonal matched count / all selected reconstruction from that truth bin, including reco outside modeled range",
            "diagonal_acceptance": "Diagonal matched count / all generated events in that truth bin",
            "efficiency_selected": "All selected reconstruction from truth bin / generated truth bin",
            "efficiency_modeled": "Selected reco in modeled range from truth bin / generated truth bin",
            "effective_mc": "(sum weight)^2 / sum(weight^2); invariant to split-exposure/luminosity scaling",
        }
    return result


def evaluate_edges(sample, observable, edges, config, split="train", shape_tilt=0.0, shape_mode="combined"):
    """Evaluate fixed edges without selecting or changing them."""
    cfg = settings(config, observable)
    edges = _edges(edges, cfg["fit_range"])
    reco = reco_edges_for_fit(edges, cfg)
    return _finish(_histograms(sample, observable, edges, reco, split, shape_tilt, shape_mode), coarse=True)


def _aggregate(hist, positions):
    """Contiguous aggregation of the precomputed candidate-grid histogram."""
    positions = np.asarray(positions, dtype=int)
    n_fine = len(hist["truth_expected"])
    intervals = list(zip(positions[:-1], positions[1:]))
    fine_fit_rows = np.asarray(hist["metadata"]["fit_reco_indices"], dtype=int)
    offset = int(fine_fit_rows[0])
    n_reco = len(hist["reco_expected"])
    rows = [(j, j + 1) for j in range(offset)]
    rows += [(a + offset, b + offset) for a, b in intervals]
    rows += [(j, j + 1) for j in range(offset + n_fine, n_reco)]
    reconstructed_edges = [float(hist["reco_edges"][rows[0][0]])]
    reconstructed_edges += [float(hist["reco_edges"][b]) for _, b in rows]
    result = {"truth_edges": hist["truth_edges"][positions],
              "reco_edges": np.asarray(reconstructed_edges),
              "metadata": dict(hist["metadata"])}
    fit_indices = fit_reco_indices(result["reco_edges"], result["truth_edges"][[0, -1]])
    result["metadata"].update(fit_reco_indices=fit_indices.tolist(),
        control_reco_indices=np.setdiff1d(np.arange(len(rows)), fit_indices).tolist())
    for key in ("truth_expected", "truth_sumw2", "truth_raw_count", "selected_truth_expected", "selected_truth_sumw2"):
        result[key] = np.array([hist[key][a:b].sum() for a, b in intervals])
    for key in ("reco_expected", "feed_in"):
        result[key] = np.array([hist[key][a:b].sum() for a, b in rows])
    for key in ("matched_expected", "matched_sumw2"):
        result[key] = np.array([[hist[key][a:b, c:d].sum() for c, d in intervals] for a, b in rows])
    categories = hist["raw_categories"]
    merged = np.array([[categories[a:b, c:d].sum() for c, d in intervals] for a, b in rows])
    missed = np.array([categories[-1, a:b].sum() for a, b in intervals])
    result["raw_categories"] = np.vstack((merged, missed)).astype(np.int64)
    return _finish(result, coarse=True)


def candidate_summary(evaluation, config, observable):
    """Serializable optimizer criteria, also useful for external common edges."""
    cfg = settings(config, observable)
    metrics = evaluation["metrics"]
    counts = metrics["expected_reco_signal"]
    row = {
        "edges": evaluation["truth_edges"].tolist(),
        "n_bins": int(len(counts)),
        "min_expected_reco": float(np.min(counts)),
        "min_purity": float(np.min(metrics["purity"])),
        "min_selected_stability": float(np.min(metrics["selected_stability"])),
        "min_diagonal_acceptance": float(np.min(metrics["diagonal_acceptance"])),
        "min_truth_effective_mc": float(np.min(metrics["truth_effective_mc"])),
        "min_selected_effective_mc": float(np.min(metrics["selected_effective_mc"])),
        "condition_number": float(evaluation["condition_number"]),
        "response_rank": int(evaluation["response_rank"]),
        "count_balance": float(np.min(counts) / np.max(counts)) if np.max(counts) > 0 else 0.0,
        "max_feed_in_fraction": float(np.max(metrics["feed_in_fraction"])),
    }
    criteria = {
        "minimum_bin_count": row["n_bins"] >= cfg["min_bins"],
        "maximum_bin_count": row["n_bins"] <= cfg["max_bins"],
        "expected_signal_reco_count": row["min_expected_reco"] >= cfg["min_expected_reco"],
        "purity_including_feed_in": row["min_purity"] >= cfg["min_purity"],
        "selected_migration_stability": row["min_selected_stability"] >= cfg["min_selected_stability"],
        "generated_mc_support": row["min_truth_effective_mc"] >= cfg["min_effective_mc"],
        "selected_mc_support": row["min_selected_effective_mc"] >= cfg["min_effective_mc"],
        "full_column_rank": row["response_rank"] == row["n_bins"],
        "conditioning": np.isfinite(row["condition_number"]) and row["condition_number"] <= cfg["max_condition_number"],
    }
    row["passes"] = bool(all(criteria.values()))
    row["failure_reasons"] = [name for name, passes in criteria.items() if not passes]
    # Nonfinite values are represented explicitly, not as nonstandard JSON Infinity.
    if not np.isfinite(row["condition_number"]):
        row["condition_number"] = None
    return row


def optimize(sample, observable, config):
    """Select bins using training MC only, never a held-out response/count."""
    cfg = settings(config, observable)
    grid = _edges(cfg["candidate_edges"], cfg["fit_range"])
    max_bins = min(cfg["max_bins"], len(grid) - 1)
    n_candidates = sum(math.comb(len(grid) - 2, n - 1) for n in range(cfg["min_bins"], max_bins + 1))
    if n_candidates > cfg["max_candidates"]:
        raise ValueError("Candidate scan exceeds max_candidates; reduce the configured grid")
    hist = _histograms(sample, observable, grid, reco_edges_for_fit(grid, cfg), "train", 0.0)
    candidates = []
    for n_bins in range(cfg["min_bins"], max_bins + 1):
        for interior in combinations(range(1, len(grid) - 1), n_bins - 1):
            evaluation = _aggregate(hist, (0,) + interior + (len(grid) - 1,))
            candidates.append(candidate_summary(evaluation, config, observable))
    passing = [row for row in candidates if row["passes"]]
    passing.sort(key=lambda row: (-row["n_bins"], row["condition_number"], -row["count_balance"], tuple(row["edges"])))
    selected = passing[0] if passing else None
    return {
        "status": "SELECTED_TRAINING_ONLY" if selected else (
            "NO_FEASIBLE_SLOPE_BINNING" if observable in ("pt2", "abs_t_via_pt2") else "NO_FEASIBLE_SPECTRUM_BINNING"),
        "observable": observable,
        "selected_edges": selected["edges"] if selected else None,
        "selected_candidate": selected,
        "candidate_grid": grid.tolist(),
        "candidates": candidates,
        "n_candidates": len(candidates),
        "n_passing": len(passing),
        "thresholds": {key: cfg[key] for key in DEFAULTS if key not in ("candidate_edges", "fine_truth_step")},
        "selection_rule": "Maximize passing bin count, then minimize response condition number, then maximize count balance; lexicographic edges break exact ties",
        "data_used": "Training split only; bin edges must be frozen before held-out and shape-stress validation",
        "scope": "Design thresholds are not a low-count confidence-coverage guarantee",
    }


def build_model(sample, observable, edges, config, split="train", shape_tilt=0.0, shape_mode="combined"):
    """Fine-truth forward-folding model and its coarse extraction companion."""
    cfg = settings(config, observable)
    edges = _edges(edges, cfg["fit_range"])
    lo, hi = cfg["fit_range"]
    # Subdivide each optimized interval separately. Merging two floating-point
    # grids via np.unique can create near-zero cells at 0.3 versus 0.30000000004.
    fine = np.concatenate([
        np.linspace(a, b, int(math.ceil((b - a) / cfg["fine_truth_step"])) + 1)[:-1]
        for a, b in zip(edges[:-1], edges[1:])
    ] + [np.array([hi])])
    reco = reco_edges_for_fit(edges, cfg)
    result = _finish(_histograms(sample, observable, fine, reco, split, shape_tilt, shape_mode), coarse=False)
    result["coarse"] = evaluate_edges(sample, observable, edges, config, split, shape_tilt, shape_mode)
    result["metadata"].update({
        "optimized_edges": edges.tolist(),
        "max_fine_truth_width": float(np.max(np.diff(fine))),
        "fine_response_assumption": "Within each fine truth cell, the detector response is averaged using its MC population; check finer-cell and shape-stress sensitivity",
    })
    return result
