"""Read saved v4 inference; never reproject counts or refit corrected points.

The existing renderers keep their typography and labels. Their input adapter
retains the source case status, sampling covariance and interval interpretation.
Rapidity deliberately uses the saved flux-normalized companion for continuity
with notebook 03; the primary ep spectrum is retained in its points CSV too.
"""
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .common import read_json, read_csv, STATES


class SavedSlopeUnavailable(ValueError):
    """The producer reports no identifiable slope or no usable interval."""


class SavedCaseUnavailable(ValueError):
    """A consistent, completed producer decision has no numerical extraction."""

    def __init__(self, message, metadata, sources, entry):
        super().__init__(message)
        self.metadata, self.sources, self.entry = metadata, sources, entry


def _unavailable_case(config, state, observable, source_config, entry,
                      statuses, table_path, paths, proton_policy):
    """Recognize only the explicit no-feasible-slope outcome, never missing data."""
    status = "NO_FEASIBLE_SLOPE_BINNING"
    optimization = entry.get("optimization") or {}
    inference = entry.get("inference") or {}
    if (observable != "pt2" or optimization.get("status") != status
            or entry.get("status", status) != status):
        raise ValueError(f"{state}/{observable}: inconsistent no-feasible-binning status.")
    if (optimization.get("selected_edges") not in (None, [])
            or inference.get("coarse") or inference.get("slope")
            or inference.get("status") not in (None, "NOT_RUN", status)
            or entry.get("cross_section_scaling")):
        raise ValueError(f"{state}/{observable}: unavailable case contains numerical extraction metadata.")
    for column in ("status", "overall_status"):
        if set(statuses[column]) != {status}:
            raise ValueError(f"{state}/{observable}: result and exported overall statuses disagree.")
    if set(statuses.stage_status) != {"NOT_RUN"}:
        raise ValueError(f"{state}/{observable}: unavailable extraction must have stage_status NOT_RUN.")
    for column, allowed in (("optimization_status", {status}), ("inference_status", {"NOT_RUN"})):
        if column in statuses and set(statuses[column]) != allowed:
            raise ValueError(f"{state}/{observable}: inconsistent {column}.")
    if (entry.get("publication_ready") not in (None, False)
            or statuses.publication_ready.astype(str).str.lower().isin(["true", "1"]).any()):
        raise ValueError("Unavailable extraction unexpectedly claims publication readiness.")
    table = read_csv(table_path, required=("state", "observable", "status", "overall_status", "stage_status"))
    if ((table.state == state) & (table.observable == observable)).any():
        raise ValueError(f"{state}/{observable}: unavailable case contains numerical CSV rows.")
    if proton_policy is not None:
        from .proton_inputs import require_policy
        if "proton_selection" in statuses:
            recorded = require_policy(json.loads(statuses.iloc[0].proton_selection), paths[2])
            if recorded != proton_policy:
                raise ValueError(f"{state}/{observable}: status-table proton-selection metadata disagree.")
        if "proton_selection_mode" in statuses and set(statuses.proton_selection_mode) != {"required"}:
            raise ValueError(f"{state}/{observable}: status-table proton-selection mode disagrees.")
    state_cfg = entry["state_config"]
    metadata = dict(state=state, observable=observable, available=False,
        source_lumi=float(source_config["target_lumi_invfb"]),
        q2_max=state_cfg["q2_max"], q2_cut_events=state_cfg["q2_max"], q2_max_flux=state_cfg["q2_max"],
        category=state_cfg["category"], w_range=source_config["w_range"],
        w_min=source_config["w_range"][0], w_max=source_config["w_range"][1],
        source=str(paths[1]), overall_status=status, stage_status="NOT_RUN", publication_ready=False,
        optimization_status=status, failed_checks=entry.get("failed_checks") or [],
        reason="No candidate binning met the configured slope-study requirements; no spectrum or slope was extracted.",
        yield_luminosity_rescaled=False, bins=0)
    if proton_policy is not None:
        metadata.update(proton_selection=proton_policy, proton_selection_mode=proton_policy["mode"])
    raise SavedCaseUnavailable(f"{state}/{observable}: {metadata['reason']}", metadata,
                               [_fingerprint(p) for p in paths], entry)


def _fingerprint(path):
    return {"path": str(path), "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def _same(actual, expected, what):
    if np.shape(actual) != np.shape(expected) or not np.allclose(actual, expected, rtol=2e-8, atol=1e-11):
        raise ValueError(f"Reoptimization {what} disagree; copy one complete run, without mixing files.")


def _case_rows(path, state, observable):
    table = read_csv(path, required=("state", "observable", "status", "overall_status", "stage_status"))
    table = table.loc[(table.state == state) & (table.observable == observable)].copy()
    if table.empty:
        raise ValueError(f"{path}: missing {state}/{observable}; no fallback to legacy results.")
    return table


def load_case(config, state, observable):
    run = Path(config.reoptimization_run)
    base = run / "artifacts/reoptimization"
    config_path = run / "effective_config.json"
    result_path = base / state / observable / "result.json"
    status_path = base / "validation_status.csv"
    table_path = base / "differential_results.csv"
    source_config, entry = read_json(config_path), read_json(result_path)
    proton_policy = None
    if config.require_proton_tag:
        from .proton_inputs import reoptimization_policy
        proton_policy = reoptimization_policy(source_config, entry, config_path, result_path)
    expected_energy = f"{source_config['electron_energy']:g}x{source_config['proton_energy']:g}"
    if config.energy != expected_energy:
        raise ValueError(f"Requested energy {config.energy} differs from saved {expected_energy}.")
    saved_lumi = float(source_config["target_lumi_invfb"])
    if not np.isclose(config.lumi_invfb, saved_lumi):
        raise ValueError(f"Saved inference is at {saved_lumi:g} fb^-1. Rerun reoptimization for another luminosity; do not rescale its intervals.")
    _same(config.coverage, source_config["inference"]["coverage"], "interval coverage")
    if entry.get("state") != state or entry.get("observable") != observable or not entry.get("complete"):
        raise ValueError(f"{result_path}: wrong case or incomplete inference.")
    state_cfg = entry["state_config"]
    _same(state_cfg["q2_max"], config.q2_max[state], "Q2 selection")
    if state_cfg["category"] != config.categories[state]:
        raise ValueError(f"{state}: configured category disagrees with saved production metadata.")
    status = str(entry["overall_status"])
    statuses = _case_rows(status_path, state, observable)
    if len(statuses) != 1:
        raise ValueError(f"{status_path}: expected one final status for {state}/{observable}.")
    if status == "NO_FEASIBLE_SLOPE_BINNING":
        _unavailable_case(config, state, observable, source_config, entry, statuses,
                          table_path, (config_path, result_path, status_path, table_path), proton_policy)
    table = _case_rows(table_path, state, observable).sort_values("bin").reset_index(drop=True)
    for source in (statuses, table):
        if set(source.status) != {status} or set(source.overall_status) != {status}:
            raise ValueError(f"{state}/{observable}: result and exported overall statuses disagree.")
        if source.publication_ready.astype(str).str.lower().isin(["true", "1"]).any():
            raise ValueError("Conditional reoptimization output unexpectedly claims publication readiness.")
    if status != "PASS_CONDITIONAL" and not config.allow_diagnostic_results:
        raise ValueError(f"{state}/{observable}: {status}; diagnostic rendering is disabled. "
                         "Inspect validation_status.csv and failed checks. To intentionally render "
                         "a diagnostic, set allow_diagnostic_results=true in the plotting profile "
                         "or pass --allow-diagnostics; its validation status will remain unchanged.")
    if status != "PASS_CONDITIONAL":
        warnings.warn(f"{state}/{observable}: {status}; rendering saved diagnostic only. Failed checks: {entry.get('failed_checks', [])}", RuntimeWarning)
    n = len(table)
    _same(table.bin.to_numpy(), np.arange(n), "bin indices")
    edges = np.asarray(entry["optimization"]["selected_edges"], float)
    _same(table.x_low.to_numpy(float), edges[:-1], "lower edges")
    _same(table.x_high.to_numpy(float), edges[1:], "upper edges")
    if not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError("Invalid saved bin edges.")
    scale = entry["cross_section_scaling"]
    _same(scale["luminosity_invfb"], saved_lumi, "conversion luminosity")
    _same(scale["branching_fraction_ee"], state_cfg["br_ee"], "branching fraction")
    factors = np.asarray(scale["factors"], float)
    _same(table.conversion_factor.to_numpy(float), factors, "conversion factors")
    base_factor = 1. / (1000. * saved_lumi * state_cfg["br_ee"])
    expected = {"w": "sigma_gamma_p_w", "pt2": "dsigma_gamma_p_dpt2", "y": "dsigma_ep_dy"}[observable]
    if scale["quantity"] != expected or set(table.quantity) != {expected}:
        raise ValueError("Unexpected cross-section convention.")
    units = "pb/GeV^2" if observable == "pt2" else "pb"
    if scale["units"] != units or set(table.quantity_units) != {units}:
        raise ValueError("Unexpected cross-section units.")
    if observable == "w":
        _same(factors, base_factor / np.asarray(scale["bin_flux"]), "W flux normalization (no Delta W)")
    elif observable == "y":
        _same(factors, base_factor / np.diff(edges), "primary ep rapidity normalization")
    else:
        _same(factors, base_factor / float(entry["flux_integral"]) / np.diff(edges), "differential normalization")
    extraction = entry["inference"]["coarse"]
    _same(extraction["truth_edges"], edges, "inference edges")
    for key in ("fitted_truth", "reference_truth", "toy_lower", "toy_upper", "sampling_sd"):
        _same(table[key].to_numpy(float), extraction[key], key)
        _same(table[key + "_cross_section"].to_numpy(float), table[key].to_numpy(float) * factors, key + " cross section")
    suffix, display_factors = "_cross_section", factors
    covariance_name = f"{state}_{observable}_cross_section_covariance.csv"
    if observable == "y":
        display_factors = np.asarray(scale["alternate_factors"], float)
        _same(display_factors, factors / float(entry["flux_integral"]), "flux-normalized rapidity")
        _same(table.alternate_conversion_factor.to_numpy(float), display_factors, "alternate conversion factors")
        if scale["alternate_quantity"] != "flux_normalized_ep_rapidity":
            raise ValueError("Rapidity companion is not the expected flux-normalized convention.")
        suffix = "_alternate_cross_section"
        covariance_name = f"{state}_{observable}_alternate_cross_section_covariance.csv"
        for key in ("fitted_truth", "reference_truth", "toy_lower", "toy_upper", "sampling_sd"):
            _same(table[key + suffix].to_numpy(float), table[key].to_numpy(float) * display_factors, key + " alternate")
    covariance_path = base / covariance_name
    covariance = np.loadtxt(covariance_path, delimiter=",", ndmin=2)
    count_covariance = np.asarray(extraction["covariance"], float)
    _same(covariance, count_covariance * display_factors[:, None] * display_factors[None, :], "full covariance")
    covariance_meta_path = base / f"{state}_{observable}_covariance_metadata.json"
    cov_meta = read_json(covariance_meta_path)
    if cov_meta["overall_status"] != status:
        raise ValueError("Covariance overall status disagrees with case status.")
    projected = table.copy()
    projected["mean"] = table["fitted_truth" + suffix]
    projected["lower"] = table["toy_lower" + suffix]
    projected["upper"] = table["toy_upper" + suffix]
    projected["x_center"] = (table.x_low + table.x_high) / 2.
    projected["bin_width"] = table.x_high - table.x_low
    projected["display_quantity"] = scale["alternate_quantity"] if observable == "y" else scale["quantity"]
    projected["target_lumi_invfb"] = saved_lumi
    projected["valid"] = np.isfinite(projected[["mean", "lower", "upper"]]).all(axis=1)
    if not projected.valid.all() or (projected.lower < 0).any() or (projected.upper < projected.lower).any():
        raise ValueError("Saved spectrum contains invalid sampling quantiles.")
    # Reco intervals may include extra control bins. Use the producer's map;
    # these means are diagnostics only and are never used to remake error bars.
    coarse = entry["validation_model"]["coarse"]
    indices = np.asarray(coarse["metadata"]["fit_reco_indices"], int)
    if indices.shape != (n,) or len(set(indices)) != n:
        raise ValueError("Invalid reconstructed-bin mapping.")
    reco_edges = np.asarray(coarse["reco_edges"], float)
    if np.any(indices < 0) or np.any(indices + 1 >= len(reco_edges)):
        raise ValueError("Reconstructed-bin index outside edge range.")
    _same(reco_edges[indices], edges[:-1], "reconstructed lower edges")
    _same(reco_edges[indices + 1], edges[1:], "reconstructed upper edges")
    projected["mean_count"] = np.asarray(coarse["reco_expected"], float)[indices]
    metadata = dict(state=state, observable=observable, available=True,
        diagnostic_rendering_explicit=bool(config.allow_diagnostic_results), q2_max=state_cfg["q2_max"],
        q2_cut_events=state_cfg["q2_max"], q2_max_flux=state_cfg["q2_max"],
        category=state_cfg["category"], w_range=source_config["w_range"],
        w_min=source_config["w_range"][0], w_max=source_config["w_range"][1],
        source_lumi=saved_lumi, source=str(result_path), overall_status=status,
        stage_status=extraction["status"], publication_ready=False,
        failed_checks=entry.get("failed_checks", []),
        uncertainty_scope=entry.get("scope", source_config.get("uncertainty_scope")),
        interval_interpretation="Central toy sampling quantiles, not confidence limits",
        cross_section_scaling=scale, displayed_quantity=projected.display_quantity.iloc[0],
        covariance=covariance.tolist(), covariance_metadata=cov_meta,
        yield_luminosity_rescaled=False, bin_center_convention="Arithmetic midpoint; value is a bin average")
    if proton_policy is not None:
        metadata.update(proton_selection=proton_policy, proton_selection_mode=proton_policy["mode"])
        projected["proton_selection_mode"] = proton_policy["mode"]
        projected["proton_selection_contract"] = proton_policy["contract"]
    sources = [_fingerprint(p) for p in (config_path, result_path, status_path, table_path, covariance_path, covariance_meta_path)]
    return projected, metadata, sources, entry


def differential_inputs(config, state, variable):
    observable = "pt2" if variable == "pT2" else "y"
    frame, meta, sources, entry = load_case(config, state, observable)
    from .rebinning import maybe_rebin
    frame, meta = maybe_rebin(config, state, observable, frame, meta, entry)
    return frame, meta, sources


def w_inputs(config, state):
    projected, meta, sources, entry = load_case(config, state, "w")
    from .rebinning import maybe_rebin
    projected, meta = maybe_rebin(config, state, "w", projected, meta, entry)
    frame = projected.rename(columns={"x_low": "W_low", "x_high": "W_high", "x_center": "W_center"}).copy()
    frame["sigma_gamma_p_pb"] = projected["mean"]
    frame["sigma_gamma_p_truth_closure_pb"] = projected["reference_truth_cross_section"]
    frame["plotted"] = projected["valid"]
    check = {"status": "passed", "scope": "Saved conversion factors and CSV/JSON agreement; not a validation of inference or production rate",
             "source_lumi_invfb": config.lumi_invfb,
             "additional_branching_or_luminosity_factor_applied": False}
    return frame, projected, meta, sources, check


def saved_slope(config, state, fit_range):
    from .differential import exponential_bin_average
    frame, meta, sources, entry = load_case(config, state, "pt2")
    edges = np.r_[frame.x_low.to_numpy(float), frame.x_high.iloc[-1]]
    model_edges = np.asarray(entry["model"]["truth_edges"], float)
    _same(list(fit_range), [model_edges[0], model_edges[-1]], "saved slope fit range")
    slope = entry["inference"]["slope"]
    if slope.get("b") is None or slope.get("amplitude") is None:
        raise SavedSlopeUnavailable(f"Saved slope unavailable: {slope.get('status')}")
    b = float(slope["b"])
    low, high = slope.get("interval_low"), slope.get("interval_high")
    if low is None:
        raise SavedSlopeUnavailable("Saved slope interval is missing.")
    if not np.isfinite([b, low]).all() or low < 0 or low > b or (high is not None and (not np.isfinite(high) or high < b)):
        raise ValueError("Saved slope interval is inconsistent.")
    amplitude_factor = np.asarray(entry["cross_section_scaling"]["factors"]) * np.diff(edges)
    _same(amplitude_factor, np.full(len(frame), entry["cross_section_factor"]), "slope amplitude normalization")
    amplitude = float(slope["amplitude"]) * amplitude_factor[0]
    fit = dict(slope, state=state, status=meta["overall_status"], overall_status=meta["overall_status"],
               stage_status=slope["status"], publication_ready=False, failed_checks=meta["failed_checks"],
               b_slope_GeVminus2=b, b_interval_low=low, b_interval_high=high,
               b_slope_err_low=b-low, b_slope_err_high=None if high is None else high-b,
               amplitude_pb_per_GeV2=amplitude, lumi_invfb=config.lumi_invfb,
               requested_fit_range_GeV2=list(fit_range), used_fit_bins_GeV2=np.column_stack([edges[:-1], edges[1:]]).tolist(),
               fit_statistic="Saved fine-response forward-folded Poisson fit; no notebook refit",
               source_files=sources, proton_selection=meta.get("proton_selection"))
    predictions = pd.DataFrame(dict(x_low=edges[:-1], x_high=edges[1:], included_in_fit=True,
        fit_bin_average_pb_per_GeV2=exponential_bin_average(edges[:-1], edges[1:], b, amplitude)))
    profile = pd.DataFrame(slope["profile"])
    for output in (predictions, profile):
        output["overall_status"] = meta["overall_status"]
        output["status"] = meta["overall_status"]
        output["stage_status"] = slope["status"]
        output["publication_ready"] = False
        if config.require_proton_tag:
            output["proton_selection_mode"] = meta["proton_selection_mode"]
    return fit, predictions, profile


def preflight_reoptimization(config, section):
    if section not in ("w", "differential"):
        raise ValueError("Reoptimization inputs apply only to notebooks 02/03. Keep notebook 01 on its tagged full-framework runs.")
    rows = []
    for state in STATES:
        for observable in (("w",) if section == "w" else ("pt2", "y")):
            try:
                frame, meta, _, _ = load_case(config, state, observable)
            except SavedCaseUnavailable as error:
                meta = error.metadata
                rows.append({key: meta[key] for key in ("state", "observable", "source_lumi", "q2_max", "category", "overall_status", "stage_status", "publication_ready")})
                rows[-1].update(available=False, bins=0, display_bins=0, reason=str(error))
                if config.require_proton_tag:
                    rows[-1]["proton_selection_mode"] = meta["proton_selection_mode"]
                warnings.warn(str(error), RuntimeWarning)
                continue
            rows.append({key: meta[key] for key in ("state", "observable", "source_lumi", "q2_max", "category", "overall_status", "stage_status", "publication_ready")})
            rows[-1]["available"] = True
            if config.require_proton_tag:
                rows[-1]["proton_selection_mode"] = meta["proton_selection_mode"]
            rows[-1]["bins"] = len(frame)
            from .rebinning import requested_edges, merge_matrix
            edges = requested_edges(config, state, observable)
            if edges is not None:
                old_edges = np.r_[frame.x_low.to_numpy(float), frame.x_high.iloc[-1]]
                matrix, _ = merge_matrix(old_edges, edges)
                rows[-1]["display_bins"] = matrix.shape[0]
            else:
                rows[-1]["display_bins"] = len(frame)
    return pd.DataFrame(rows)
