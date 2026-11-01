"""Figure 3: differential spectra and finite-bin projected slope forecasts.

No ROOT dependency. Inputs are the framework's bin-by-bin extraction CSVs.
The fit is conditional on fixed correction factors and a signal-only template.
Its profile-deviance intervals are asymptotic forecasts, not exact low-count CIs.
"""

from pathlib import Path
import hashlib
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import minimize_scalar, brentq
from scipy.special import logsumexp, xlogy
from scipy.stats import chi2

from .common import (
    COLORS, MARKERS, FONT_SIZES, annotate, draw_intervals, project_cross_section,
    read_csv, resolve_artifacts, sample_metadata, save_figure, set_style,
    state_label, write_json,
)


SLOPE_REFERENCES = [
    dict(experiment="ZEUS", state="1S", observable="abs_t", slope=4.3,
         stat_low=1.3, stat_high=2.0, syst_low=0.6, syst_high=0.5,
         x_min=0.0, x_max=5.0, w_min=60.0, w_max=220.0, q2_max=1.0,
         kinematics="60 < W < 220 GeV; Q2 < 1 GeV2; |t| < 5 GeV2",
         source="https://arxiv.org/pdf/1111.2133", doi="10.1016/j.physletb.2012.01.009",
         provenance="Sections 4 and 6-8; the reconstructed distribution is fitted in four |t| bins up to 5 GeV2."),
    dict(experiment="CMS", state="nS", observable="pT2", slope=6.0,
         stat_low=2.1, stat_high=2.1, syst_low=0.3, syst_high=0.3,
         x_min=0.01, x_max=1.0, w_min=91.0, w_max=826.0, q2_max=None,
         kinematics="91 < W < 826 GeV; |y| < 2.2; 0.01 < pT2 < 1 GeV2; combined nS",
         source="https://arxiv.org/pdf/1809.11080", doi="10.1140/epjc/s10052-019-6774-8",
         provenance="Section 7.1; slope of branching-weighted combined Upsilon(nS) pPb spectrum, not state-resolved gamma-p slope.",
         erratum="https://link.springer.com/content/pdf/10.1140/epjc/s10052-022-10276-2.pdf",
         erratum_note="The 2022 erratum corrects an author name only; the quoted slope is unchanged."),
]


def exponential_bin_average(low, high, slope, amplitude=1.0):
    """Average of A exp(-b x) across each bin (A has density units)."""
    low, high = np.broadcast_arrays(np.asarray(low, float), np.asarray(high, float))
    width = high - low
    if np.any(width <= 0) or not np.isfinite(slope) or slope < 0:
        raise ValueError("Exponential average needs positive bin widths and nonnegative slope")
    if abs(slope) < 1e-10:
        return np.full_like(low, float(amplitude))
    return amplitude * np.exp(-slope * low) * (-np.expm1(-slope * width)) / (slope * width)


def fit_projected_exponential(frame, fit_range=(0.0, 1.0), delta_deviance=1.0):
    """Profile Poisson deviance of target-luminosity expected reconstructed means.

    At each b the density normalization A is profiled analytically. The model
    is bin-integrated before conversion back to reconstructed means, using the
    fixed per-bin cross-section/count factors. Fractional values here are Asimov
    expected means, never rounded observations. The displayed Garwood endpoints
    are not used as likelihood variances.
    """
    lo = frame["x_low"].to_numpy(float)
    hi = frame["x_high"].to_numpy(float)
    counts = frame["mean_count"].to_numpy(float)
    factors = frame["count_to_cross_section"].to_numpy(float)
    mask = frame["valid"].to_numpy(bool).copy()
    low_fit, high_fit = map(float, fit_range)
    if not np.isfinite([low_fit, high_fit]).all() or low_fit < 0 or high_fit <= low_fit:
        raise ValueError("fit_range must be finite, ordered and nonnegative")
    if not np.isfinite(delta_deviance) or delta_deviance <= 0:
        raise ValueError("delta_deviance must be positive")
    geometric = (lo >= low_fit - 1e-10) & (hi <= high_fit + 1e-10)
    mask &= geometric & np.isfinite(counts) & (counts >= 0)
    mask &= np.isfinite(factors) & (factors > 0) & (hi > lo)
    if mask.sum() < 3:
        raise ValueError("Slope fit needs at least three supported complete bins in fit_range; widen the range or regenerate finer supported bins")
    low, high, n, c = lo[mask], hi[mask], counts[mask], factors[mask]
    n_total = float(n.sum())
    if n_total <= 0:
        raise ValueError("Slope is not identifiable with zero expected reconstructed yield")

    def log_shape(slope):
        width = high - low
        if slope < 1e-10:
            return -np.log(c)
        return (-slope * low + np.log(-np.expm1(-slope * width))
                - np.log(slope) - np.log(width) - np.log(c))

    def profile(slope, return_parameters=False):
        log_h = log_shape(slope)
        log_sum_h = logsumexp(log_h)
        log_mu = np.log(n_total) + log_h - log_sum_h
        # Stable deviance; sum(mu)=sum(n) by the profiled normalization.
        deviance = float(2.0 * np.sum(xlogy(n, n) - n * log_mu))
        if return_parameters:
            return deviance, float(np.log(n_total) - log_sum_h), np.exp(log_mu)
        return deviance

    upper_search = 50.0
    result = minimize_scalar(profile, bounds=(0.0, upper_search), method="bounded",
                             options={"xatol": 1e-9})
    while result.x > 0.95 * upper_search and upper_search < 1000:
        upper_search *= 2
        result = minimize_scalar(profile, bounds=(0.0, upper_search), method="bounded",
                                 options={"xatol": 1e-9})
    if not result.success:
        raise RuntimeError("Projected exponential fit failed: " + str(result.message))
    slope = float(result.x) if profile(result.x) < profile(0.0) else 0.0
    d_min, log_amplitude, mu = profile(slope, True)
    threshold = d_min + delta_deviance
    at_lower_boundary = profile(0.0) <= threshold
    lower = 0.0 if at_lower_boundary else float(brentq(lambda b: profile(b) - threshold, 0.0, slope))
    upper = max(slope + 1.0, 2 * slope)
    while profile(upper) < threshold and upper < 10000:
        upper *= 2.0
    upper = (float(brentq(lambda b: profile(b) - threshold, slope, upper))
             if profile(upper) >= threshold else None)
    amplitude = float(np.exp(log_amplitude))
    predictions = exponential_bin_average(lo, hi, slope, amplitude)
    profile_max = max(slope + 3.0, upper * 1.2 if upper is not None else slope * 2 + 3)
    grid = np.linspace(0.0, profile_max, 241)
    summary = {
        "b_slope_GeVminus2": slope, "b_interval_low": lower,
        "b_interval_high": upper, "b_slope_err_low": slope - lower,
        "b_slope_err_high": None if upper is None else upper - slope,
        "amplitude_pb_per_GeV2": amplitude, "deviance": max(d_min, 0.0),
        "ndof": int(mask.sum() - 2), "expected_reconstructed_count_in_fit": n_total,
        "requested_fit_range_GeV2": [low_fit, high_fit],
        "used_fit_bins_GeV2": np.column_stack([low, high]).tolist(),
        "included_bin_indices": np.flatnonzero(mask).tolist(),
        "excluded_bin_indices": np.flatnonzero(~mask).tolist(),
        "partly_overlapping_bins_excluded": np.flatnonzero(
            ~geometric & (hi > low_fit) & (lo < high_fit)).tolist(),
        "delta_deviance": float(delta_deviance),
        "asymptotic_nominal_coverage": float(chi2.cdf(delta_deviance, 1)),
        "lower_interval_at_physical_boundary": bool(at_lower_boundary),
        "upper_interval_unbounded": upper is None,
        "fit_statistic": "Poisson deviance on expected reconstructed means; finite-bin exponential; fixed correction factors",
        "interval_interpretation": "Conditional Asimov profile-deviance forecast with a one-parameter chi-square threshold; not an exact low-count confidence interval.",
        "limitations": ["No background, mass-component separation uncertainty, efficiency-MC uncertainty or systematics",
                        "Bin-by-bin correction; no response-matrix forward folding",
                        "A pT2 slope is not automatically the proton-vertex |t| slope"],
    }
    bin_table = pd.DataFrame({"x_low": lo, "x_high": hi,
                              "fit_bin_average_pb_per_GeV2": predictions,
                              "included_in_fit": mask})
    profile_table = pd.DataFrame({"b_GeVminus2": grid,
                                  "delta_deviance": [max(profile(b) - d_min, 0) for b in grid]})
    return summary, bin_table, profile_table


def _bool(series):
    return series.astype(str).str.lower().isin(["true", "1", "1.0"]).to_numpy()


def _hash(path):
    return {"path": str(path), "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def _differential_inputs(config, state, variable):
    if config.reoptimization_run is not None:
        from .reoptimization_inputs import differential_inputs
        return differential_inputs(config, state, variable)
    artifacts = resolve_artifacts(config.runs[state])
    stem = f"{config.energy}_{state.lower()}"
    path = artifacts / "differential" / f"{stem}_{state}_{variable}.csv"
    frame = read_csv(path, required=("x_low", "x_high", "x_center", "bin_width",
                                     "n_reco_visible", "dsigma_dx_pb", "passes_min_eff",
                                     "output_convention"), state=state)
    frame = frame.sort_values("x_low").reset_index(drop=True)
    if set(frame["output_convention"].astype(str)) != {"inclusive_upsilon"}:
        raise ValueError(f"{path}: expected inclusive_upsilon cross sections; rerun with the desired convention instead of silently changing branching fractions")
    low, high = frame.x_low.to_numpy(float), frame.x_high.to_numpy(float)
    if np.any(high <= low) or np.any(low[1:] < high[:-1] - 1e-10):
        raise ValueError(f"{path}: overlapping or invalid bins")
    if not np.allclose(high - low, frame.bin_width.to_numpy(float)):
        raise ValueError(f"{path}: bin_width disagrees with edges")
    meta_path = artifacts / "differential" / f"{stem}_metadata.csv"
    meta = read_csv(meta_path, required=("q2_max_flux", "q2_cut_events", "w_min", "w_max"), state=state)
    if len(meta) != 1:
        raise ValueError(f"{meta_path}: expected one metadata row for {state}")
    metadata = meta.iloc[0].to_dict()
    if not np.isclose(metadata["q2_max_flux"], metadata["q2_cut_events"]):
        raise ValueError(f"{meta_path}: flux and event Q2 limits differ; cannot publish this normalization")
    sample = sample_metadata(config, state)
    if not np.isclose(float(metadata["q2_cut_events"]), float(sample["q2_max"])):
        raise ValueError(f"{meta_path}: Q2 cut disagrees with the sample performance audit")
    if not np.allclose([metadata["w_min"], metadata["w_max"]], sample["w_range"]):
        raise ValueError(f"{meta_path}: W interval disagrees with the sample performance audit")
    scale_column = "reco_count_to_dsigma_dx_pb" if "reco_count_to_dsigma_dx_pb" in frame else None
    projected = project_cross_section(frame, config.lumi_invfb, "dsigma_dx_pb",
                                     count_column="n_reco_visible", scale_column=scale_column,
                                     coverage=config.coverage)
    projected["valid"] &= _bool(frame["passes_min_eff"])
    projected["x_low"], projected["x_high"] = low, high
    return projected, dict(metadata, category=sample["category"], input_validation=sample), [_hash(path), _hash(meta_path), _hash(sample["source"])]


def _slope_inset(ax, fits, references):
    # A shorter legend places its Fit entry higher. Shift the smaller inset
    # right when states are unavailable so its y label cannot overlap that entry.
    available_states = [state for state in ("1S", "2S", "3S")
                        if state in fits and fits[state].get("b_interval_high") is not None]
    bounds = [0.62, 0.51, 0.35, 0.39] if len(available_states) < 3 else [0.53, 0.51, 0.44, 0.39]
    inset = ax.inset_axes(bounds)
    inset.set_facecolor("white")
    labels = []
    for i, state in enumerate(available_states):
        fit = fits[state]
        #labels.append(state)
        labels.append("ePIC" if state == "1S" else state)
        inset.errorbar(i, fit["b_slope_GeVminus2"],
                       yerr=[[fit["b_slope_err_low"]], [fit["b_slope_err_high"]]],
                       fmt=MARKERS[state], color=COLORS[state], ms=4, capsize=2, elinewidth=1.0)
    for j, ref in enumerate(references, len(available_states)):
        low = np.hypot(float(ref["stat_low"]), float(ref["syst_low"]))
        high = np.hypot(float(ref["stat_high"]), float(ref["syst_high"]))
        inset.errorbar(j, float(ref["slope"]), yerr=[[low], [high]],
                       fmt="s" if ref["experiment"] == "ZEUS" else "D",
                       color="#62676c", mfc="white", ms=4, capsize=2, elinewidth=1.0)
        labels.append(ref["experiment"])
    inset.set_xticks(range(len(labels)), labels, fontsize=FONT_SIZES["inset_tick"])
    inset.set_ylabel(r"$b$ [GeV$^{-2}$]", fontsize=FONT_SIZES["inset_label"], labelpad=3)
    inset.set_ylim(bottom=0)
    inset.set_xlim(-0.6, len(labels)-0.4)
    inset.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=FONT_SIZES["inset_tick"])
    top = max([8.6] + [fit["b_interval_high"] * 1.25 for fit in fits.values()
                      if fit.get("b_interval_high") is not None])
    inset.set_ylim(0, top)
    return inset


def plot_differential(config, fit_range=(0.0, 1.0), reference_csv=None):
    """Create the note's combined pT2/rapidity figure and auditable fit tables.

    Returns (figure, audit). ``reference_csv=None`` uses the verified ZEUS/CMS
    slope values embedded above; a supplied CSV must use the same numeric keys.
    Panel (a) displays complete bins within ``fit_range``. All source bins are
    retained in the points CSV, with an explicit ``displayed`` flag.
    """
    set_style()
    config.output_dir.mkdir(parents=True, exist_ok=True)
    references = (SLOPE_REFERENCES if reference_csv is None else
                  read_csv(Path(reference_csv), required=("experiment", "slope", "stat_low", "stat_high", "syst_low", "syst_high", "kinematics", "source")).to_dict("records"))
    fit_range = tuple(map(float, fit_range))
    if len(fit_range) != 2 or not np.isfinite(fit_range).all() or not 0 <= fit_range[0] < fit_range[1]:
        raise ValueError("fit_range must be finite, ordered and nonnegative")
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 8.0))
    fig.subplots_adjust(left=0.16, right=0.97, top=0.81, bottom=0.10, hspace=0.40)
    audit = {"figure": "upsilon_es_differential", "energy": config.energy,
             "target_lumi_invfb": config.lumi_invfb,
             "uncertainty": "Signal-only projected counts, fixed bin-by-bin corrections; median expected Garwood display intervals",
             "source_files": [], "samples": {}, "slope_fits": {},
             "pT2_display_range_GeV2": list(fit_range),
             "pT2_display_policy": "Only complete bins within the requested fit interval are shown; all source bins remain in the points CSV.",
             "slope_references": references,
             "warnings": ["State-resolved signal-only errors do not include simultaneous mass-fit correlations or continuum backgrounds."]}
    fits, tables, all_points, availability = {}, {}, [], []
    for state in ("1S", "2S", "3S"):
        tables[state] = {}
        for variable, ax in zip(("pT2", "y"), axes):
            try:
                frame, metadata, sources = _differential_inputs(config, state, variable)
            except ValueError as error:
                if config.reoptimization_run is None:
                    raise
                from .reoptimization_inputs import SavedCaseUnavailable
                if not isinstance(error, SavedCaseUnavailable):
                    raise
                metadata, sources = dict(error.metadata), error.sources
                audit["source_files"].extend(sources)
                audit["samples"][state] = metadata
                audit.setdefault("validation", {}).setdefault(state, {})[variable] = metadata
                availability.append(dict(state=state, variable=variable,
                    observable="pt2" if variable == "pT2" else "y", available=False,
                    bins=0, displayed_bins=0, status=metadata["overall_status"],
                    overall_status=metadata["overall_status"], stage_status="NOT_RUN",
                    publication_ready=False, failed_checks=metadata.get("failed_checks", []),
                    reason=str(error)))
                if variable == "pT2":
                    audit["slope_fits"][state] = dict(state=state,
                        status=metadata["overall_status"], overall_status=metadata["overall_status"],
                        stage_status="NOT_RUN", publication_ready=False,
                        failed_checks=metadata.get("failed_checks", []), reason=str(error),
                        proton_selection=metadata.get("proton_selection"),
                        source_files=sources)
                audit["warnings"].append(str(error))
                warnings.warn(str(error), RuntimeWarning)
                continue
            audit["source_files"].extend(sources)
            audit["samples"][state] = metadata
            if config.reoptimization_run is not None:
                audit.setdefault("validation", {}).setdefault(state, {})[variable] = metadata
                if metadata["overall_status"] != "PASS_CONDITIONAL":
                    audit["warnings"].append(f"{state}/{variable}: {metadata['overall_status']}; saved diagnostic, not a validated measurement.")
            tables[state][variable] = frame
            mask = frame["valid"].to_numpy(bool) & np.isfinite(frame["mean"])
            if variable == "pT2":
                mask &= frame["mean"].to_numpy(float) > 0
                mask &= ((frame["x_low"].to_numpy(float) >= fit_range[0] - 1e-10)
                         & (frame["x_high"].to_numpy(float) <= fit_range[1] + 1e-10))
            frame["displayed"] = mask
            availability.append(dict(state=state, variable=variable,
                observable="pt2" if variable == "pT2" else "y", available=True,
                bins=len(frame), displayed_bins=int(mask.sum()),
                status=metadata.get("overall_status", "LEGACY_INPUT"),
                overall_status=metadata.get("overall_status", "LEGACY_INPUT"),
                stage_status=metadata.get("stage_status"), publication_ready=False,
                failed_checks=metadata.get("failed_checks", []), reason=""))
            x = frame["x_center"].to_numpy(float)
            half = 0.5 * frame["bin_width"].to_numpy(float)
            label = state_label(state)
            draw_intervals(ax, x[mask], frame["mean"].to_numpy(float)[mask],
                           frame["lower"].to_numpy(float)[mask], frame["upper"].to_numpy(float)[mask],
                           xerr=half[mask], fmt=MARKERS[state], color=COLORS[state],
                           label=label, ms=4.5, capsize=2, elinewidth=1.0, lw=1.0)
            rows = frame.copy()
            rows["state"], rows["variable"] = state, variable
            rows["displayed"] = mask
            all_points.append(rows)
            if variable == "pT2":
                try:
                    if config.reoptimization_run is not None:
                        from .reoptimization_inputs import saved_slope
                        fit, predictions, profile = saved_slope(config, state, fit_range)
                    else:
                        fit, predictions, profile = fit_projected_exponential(
                            frame, fit_range, delta_deviance=chi2.ppf(config.coverage, 1))
                    fit.update(state=state, lumi_invfb=config.lumi_invfb)
                    fits[state] = fit
                    audit["slope_fits"][state] = fit
                    predictions.insert(0, "state", state)
                    predictions.to_csv(config.output_dir / f"upsilon_es_{state.lower()}_pt2_fit_bins.csv", index=False)
                    profile.to_csv(config.output_dir / f"upsilon_es_{state.lower()}_pt2_fit_profile.csv", index=False)
                    inside = predictions["included_in_fit"].to_numpy(bool)
                    xx = np.linspace(fit_range[0], fit_range[1], 240)
                    ax.plot(xx, fit["amplitude_pb_per_GeV2"] * np.exp(-fit["b_slope_GeVminus2"]*xx),
                            color=COLORS[state], ls="--", lw=1.2, alpha=0.9, zorder=1)
                    for row in predictions.loc[inside].itertuples():
                        ax.plot([row.x_low, row.x_high], [row.fit_bin_average_pb_per_GeV2]*2,
                                color=COLORS[state], lw=1.6, alpha=0.9, zorder=1)
                except ValueError as error:
                    if config.reoptimization_run is not None:
                        from .reoptimization_inputs import SavedSlopeUnavailable
                        if not isinstance(error, SavedSlopeUnavailable):
                            raise
                    audit["slope_fits"][state] = {"status": "unavailable", "reason": str(error)}
                    if config.reoptimization_run is not None:
                        audit["slope_fits"][state].update(
                            state=state, status=metadata["overall_status"],
                            overall_status=metadata["overall_status"], stage_status="UNAVAILABLE",
                            publication_ready=False, failed_checks=metadata["failed_checks"])
                    audit["warnings"].append(f"{state} slope unavailable: {error}")
                    warnings.warn(f"{state}: {error}", RuntimeWarning)
        audit["samples"][state]["unsupported_bins"] = {
            variable: np.flatnonzero(~frame["valid"].to_numpy(bool)).tolist()
            for variable, frame in tables[state].items()}

    old_text_count = len(axes[0].texts)
    annotate(axes[0], config, panel="a")
    header_texts = axes[0].texts[old_text_count:]
    # Publication header above the panel; slopes are now inside the panel.
    for text, position in zip(header_texts, ((.035, 1.43), (.035, 1.27), (.025, .95))):
        text.set_position(position)
    header_texts[-1].set_ha("left")
    # Main labels follow the underlying truth-bin convention of the extraction.
    axes[0].set_xlabel(r"$p_{T,\Upsilon}^{2}$ [GeV$^2$]")
    axes[0].set_ylabel(r"$d\sigma_{\gamma p}/dp_{T,\Upsilon}^{2}$ [pb/GeV$^2$]")
    axes[0].set_yscale("log")
    axes[0].set_xlim(*fit_range)
    pt2_tables = {state: pair["pT2"] for state, pair in tables.items() if "pT2" in pair}
    y_tables = {state: pair["y"] for state, pair in tables.items() if "y" in pair}
    positive = (np.concatenate([frame.loc[frame["displayed"], "mean"].to_numpy(float)
                               for frame in pt2_tables.values()]) if pt2_tables else np.array([]))
    if len(positive):
        positive_lowers = np.concatenate([frame.loc[frame["displayed"] & (frame["lower"] > 0), "lower"].to_numpy(float) for frame in pt2_tables.values()])
        uppers = np.concatenate([frame.loc[frame["displayed"] & np.isfinite(frame["upper"]), "upper"].to_numpy(float) for frame in pt2_tables.values()])
        floor = min(positive.min()/5, positive_lowers.min()/1.5 if len(positive_lowers) else np.inf)
        floor = max(floor, np.finfo(float).tiny)
        ceiling = max(uppers.max() if len(uppers) else positive.max(), positive.max()) * 3.0
        # Reserve space for the inset and the category labels below its axes.
        inset_left = fit_range[0] + 0.53 * (fit_range[1] - fit_range[0])
        right_uppers = np.concatenate([
            f.loc[f["displayed"] & (f.x_high >= inset_left)
                  & np.isfinite(f["upper"]), "upper"].to_numpy(float)
            for f in pt2_tables.values()])
        if len(right_uppers):
            ceiling = max(ceiling, floor * (max(right_uppers.max() * 1.25 / floor, 1.0)) ** (1 / 0.43))
        # Reserve the upper-left block for the relocated legend as well as
        # the upper-right slope inset. Keep labels clear of the data/fit.
        legend_right = fit_range[0] + 0.51 * (fit_range[1] - fit_range[0])
        left_uppers = np.concatenate([
            f.loc[f["displayed"] & (f.x_low <= legend_right)
                  & np.isfinite(f["upper"]), "upper"].to_numpy(float)
            for f in pt2_tables.values()])
        fit_tops = [f["amplitude_pb_per_GeV2"] * np.exp(-f["b_slope_GeVminus2"] * fit_range[0])
                    for f in fits.values()]
        legend_under = max(list(left_uppers) + fit_tops + [floor])
        ceiling = max(ceiling, floor * max(legend_under * 1.15 / floor, 1.) ** (1 / 0.40))
        axes[0].set_ylim(floor, ceiling)
        for state, points in pt2_tables.items():
            at_zero = points["displayed"] & (points["lower"] <= 0)
            axes[0].plot(points.loc[at_zero, "x_center"], np.full(at_zero.sum(), floor*1.08),
                         "v", color=COLORS[state], ms=3, clip_on=False)
        audit["zero_lower_endpoint_display"] = "Downward triangles at the log-axis floor mark limits that extend to zero."
    handles, labels = axes[0].get_legend_handles_labels()
    if fits:
        handles.append(Line2D([], [], color="#555555", lw=1.4, ls="--"))
        labels.append(rf"Fit: ${fit_range[0]:g}<p_T^2<{fit_range[1]:g}$ GeV$^2$")
    if handles:
        axes[0].legend(handles, labels, loc="upper left", bbox_to_anchor=(0.02, 0.88),
                       frameon=False, fontsize=FONT_SIZES["legend"], ncol=1, labelspacing=.25,
                       columnspacing=1.0, handlelength=1.5)
    if any(fit.get("b_interval_high") is not None for fit in fits.values()):
        _slope_inset(axes[0], fits, references)
    if not pt2_tables:
        axes[0].set_ylim(1., 10.)
        axes[0].text(.5, .5, "No available $p_T^2$ spectrum", transform=axes[0].transAxes,
                     ha="center", va="center", fontsize=FONT_SIZES["legend"])
    axes[1].set_xlabel(r"$y_{\Upsilon}$ (laboratory)")
    axes[1].set_ylabel(r"$d\sigma_{\gamma p}/dy_{\Upsilon}$ [pb]")
    axes[1].text(0.025, 0.95, "(b)", transform=axes[1].transAxes, va="top", fontsize=FONT_SIZES["panel"])
    if y_tables:
        axes[1].set_xlim(min(f.x_low.min() for f in y_tables.values()), max(f.x_high.max() for f in y_tables.values()))
    y_uppers = (np.concatenate([f.loc[f["displayed"] & np.isfinite(f["upper"]), "upper"].to_numpy(float)
                                for f in y_tables.values()]) if y_tables else np.array([]))
    axes[1].set_ylim(0, max(y_uppers.max() * 1.35, 1.0) if len(y_uppers) else 1.0)
    if y_tables:
        axes[1].legend(loc="upper right", frameon=False, fontsize=FONT_SIZES["legend"], labelspacing=.25)
    else:
        axes[1].text(.5, .5, "No available rapidity spectrum", transform=axes[1].transAxes,
                     ha="center", va="center", fontsize=FONT_SIZES["legend"])
    q2_text = "; ".join(f"{state}: {m['q2_cut_events']:g}" for state, m in audit["samples"].items())
    ranges = {(float(m["w_min"]), float(m["w_max"])) for m in audit["samples"].values()}
    if len(ranges) == 1:
        wl, wh = next(iter(ranges))
        selection_text = rf"${wl:g}<W_{{\gamma p}}<{wh:g}$ GeV; $Q^2_{{\max}}$ [GeV$^2$]: " + q2_text
    else:
        selection_text = "State-specific W ranges: " + "; ".join(f"{s} {m['w_min']:g}-{m['w_max']:g}" for s,m in audit["samples"].items()) + " GeV"
        audit["warnings"].append("The states do not share a common W interval.")
    categories = "; ".join(f"{s}: {m['category']}" for s, m in audit["samples"].items())
    reference_notes = (r"ZEUS: $\Upsilon(1S)$, $|t|<5$ GeV$^2$, $60<W<220$ GeV, $Q^2<1$ GeV$^2$." + "\n"
                       r"CMS: $\Upsilon(nS)$, $0.01<p_T^2<1$ GeV$^2$, $|y|<2.2$, $91<W<826$ GeV.")
    if reference_csv is not None:
        reference_notes = "\n".join(f"{r['experiment']}: {r['kinematics']}" for r in references)
    footer = (selection_text + "\n" + categories + "\n"
              "Signal-only statistics; fixed bin-by-bin corrections. Slopes: conditional profile forecasts.\n"
              + reference_notes + "\n"
              "Reference errors: stat. + syst. in quadrature; different regions. Downward marks: zero limits.")
    # Caption information remains auditable without a block of text on the PDF.
    audit["caption_notes"] = footer.splitlines()
    audit["publication_footer_displayed"] = False
    if len({float(m["q2_cut_events"]) for m in audit["samples"].values()}) > 1:
        audit["warnings"].append("Different Q2 regions across states; state-to-state ratios are not matched-acceptance measurements.")
    if config.reoptimization_run is not None:
        audit.update(input_run=str(config.reoptimization_run), publication_ready=False,
                     uncertainty="Saved central toy sampling quantiles, not confidence limits; full covariance retained in validation metadata",
                     slope_source="Saved forward-folded fine-response Poisson profiles; no notebook refit",
                     rapidity_convention="Saved flux_normalized_ep_rapidity companion; primary ep values retained in points CSV",
                     yield_luminosity_rescaled=False)
        audit["caption_notes"] = [selection_text, categories,
                                  "Signal-only fixed-response and fixed-feed-in projection; response-MC variations remain separate.",
                                  "Displayed bars: central toy sampling quantiles. Slope inset: saved conditional profile endpoints.",
                                  "Final validation status applies to each state/observable; NOT_VALIDATED results are diagnostics.",
                                  "Rapidity retains (1/Phi_total) d sigma_ep/dy, a beam/selection-dependent flux-normalized convention.",
                                  reference_notes]
    points = (pd.concat(all_points, ignore_index=True) if all_points else
              pd.DataFrame(columns=["state", "variable", "x_low", "x_high", "mean", "lower", "upper", "displayed"]))
    points.to_csv(config.output_dir / "upsilon_es_differential_points.csv", index=False)
    audit["case_availability"] = availability
    availability_path = config.output_dir / "upsilon_es_differential_availability.csv"
    pd.DataFrame(availability).to_csv(availability_path, index=False)
    audit["availability_file"] = str(availability_path)
    slope_rows = audit["slope_fits"].values() if config.reoptimization_run is not None else fits.values()
    pd.DataFrame(slope_rows).to_csv(config.output_dir / "upsilon_es_slopes.csv", index=False)
    audit["output_files"] = save_figure(fig, config, "upsilon_es_differential")
    write_json(config.output_dir / "upsilon_es_differential_audit.json", audit)
    return fig, audit
