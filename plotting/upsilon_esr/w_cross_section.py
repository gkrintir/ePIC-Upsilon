"""Note W figure from saved bin-by-bin cross sections, without renormalizing them.

The default reference curve is the *same simulated sample's* saved truth cross
section. It is deliberately not labelled as an independent generator calculation.
Supply ``model_csv`` only for a separately evaluated, flux-averaged input model.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .common import (
    COLORS, MARKERS, FONT_SIZES, annotate, draw_intervals, project_cross_section,
    read_csv, read_json, resolve_artifacts, sample_metadata, save_figure,
    set_style, source_luminosity, state_label, write_json,
)


ZEUS_SOURCE = {
    "collaboration": "ZEUS",
    "citation": "S. Chekanov et al., Phys. Lett. B 680 (2009) 4-12, Table 1",
    "doi": "10.1016/j.physletb.2009.07.066",
    "url": "https://arxiv.org/pdf/0903.4205",
    "table": "Table 1, PDF page 19 (zero-based page 18)",
    "q2_max_GeV2": 1.0,
    "state": "1S",
    "interpretation": (
        "Published cross sections at reference W0=100 and 180 GeV, not bin "
        "midpoints. Statistical and asymmetric systematic errors are added "
        "in quadrature on each side. The paper excludes its 2.6% luminosity "
        "uncertainty; no modern branching-fraction rescaling is applied."
    ),
}


def zeus_reference():
    """Published Table 1 numbers; retain asymmetric errors and W0 values."""
    result = pd.DataFrame({
        "state": ["1S", "1S"], "W_low": [60.0, 130.0],
        "W_high": [130.0, 220.0], "W_ref": [100.0, 180.0],
        "sigma_pb": [160.0, 321.0], "stat_pb": [51.0, 88.0],
        "syst_low_pb": [21.0, 114.0], "syst_high_pb": [48.0, 46.0],
    })
    result["err_low_pb"] = np.hypot(result.stat_pb, result.syst_low_pb)
    result["err_high_pb"] = np.hypot(result.stat_pb, result.syst_high_pb)
    return result


def _fingerprint(path):
    path = Path(path)
    return {"path": str(path.resolve()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _bool(values):
    return pd.Series(values).astype(str).str.lower().isin(["true", "1", "1.0"]).to_numpy()


def _reference(reference_csv):
    if reference_csv is None:
        return zeus_reference(), ZEUS_SOURCE
    path = Path(reference_csv)
    frame = read_csv(path, required=("state", "W_low", "W_high", "W_ref",
                                    "sigma_pb", "err_low_pb", "err_high_pb"), state="1S")
    sidecar = path.with_suffix(".json")
    if not sidecar.exists():
        raise ValueError(f"Custom reference requires a provenance sidecar: {sidecar}")
    metadata = read_json(sidecar)
    if not metadata.get("citation") or not metadata.get("url"):
        raise ValueError("Custom reference metadata requires citation and url.")
    return frame, metadata


def _model(model_csv, state, frame, config, metadata):
    """Require exact matching bins/phase space, rather than interpolating a curve."""
    path = Path(model_csv)
    model = read_csv(path, required=("state", "energy", "W_low", "W_high",
                                    "sigma_pb", "q2_max_GeV2"), state=state)
    model = model.loc[model.energy.astype(str) == str(config.energy)].copy()
    sidecar = path.with_suffix(".json")
    if not sidecar.exists():
        raise ValueError(f"Independent model requires provenance sidecar: {sidecar}")
    provenance = read_json(sidecar)
    for key in ("label", "generator_revision", "source", "averaging"):
        if not provenance.get(key):
            raise ValueError(f"Independent model metadata missing {key} in {sidecar}")
    if provenance["averaging"] != "same_flux_and_bins":
        raise ValueError("Model sidecar averaging must be 'same_flux_and_bins'.")
    model = model.sort_values(["W_low", "W_high"]).reset_index(drop=True)
    if len(model) != len(frame) or not np.allclose(
            model[["W_low", "W_high"]], frame[["W_low", "W_high"]],
            rtol=1e-10, atol=1e-10):
        raise ValueError(f"{state}: independent model and framework W bins differ.")
    if not np.allclose(model.q2_max_GeV2, metadata["q2_max"], rtol=1e-10, atol=0):
        raise ValueError(f"{state}: independent model and sample Q2 ranges differ.")
    if not np.all(np.isfinite(model.sigma_pb) & (model.sigma_pb >= 0)):
        raise ValueError("Independent model has nonfinite or negative cross sections.")
    provenance = {**provenance, "file": _fingerprint(path),
                  "sidecar": _fingerprint(sidecar)}
    return model.sigma_pb.to_numpy(float), str(provenance["label"]), provenance


def _check_normalization(frame, normalization, state):
    """Audit saved inclusive values using their own luminosity, flux and BR.

    No factor from this check is applied to the plotted cross sections.
    The notebook consumes the current framework's visible-decay convention.
    """
    if len(normalization) != 1:
        raise ValueError(f"{state}: expected exactly one normalization row.")
    br = float(normalization.iloc[0]["br_ee"])
    if not 0 < br < 1:
        raise ValueError(f"{state}: invalid saved branching fraction {br}.")
    lumi = source_luminosity(frame)
    flux = frame.flux_int.to_numpy(float)
    truth = frame.n_true_weighted_visible.to_numpy(float)
    saved = frame.sigma_gamma_p_truth_closure_pb.to_numpy(float)
    if np.any(~np.isfinite(flux)) or np.any(flux <= 0):
        raise ValueError(f"{state}: saved flux integral must be positive and finite.")
    calculated = truth / (1000.0 * lumi * br * flux)
    if not np.allclose(calculated, saved, rtol=1e-5, atol=1e-10, equal_nan=True):
        raise ValueError(
            f"{state}: saved truth cross section is inconsistent with "
            "Ntruth/(1000 * L * B_ee * flux). Check the source convention; "
            "the notebook will not repair or renormalize it silently."
        )
    return {"status": "passed", "saved_br_ee": br,
            "source_lumi_invfb": lumi,
            "test": "saved truth sigma = Ntruth / (1000 L B_ee flux)",
            "additional_branching_or_luminosity_factor_applied": False,
            "normalization_row": normalization.iloc[0].to_dict()}


def _steps(ax, frame, values, **kwargs):
    """Draw bin averages as disconnected horizontal segments, without interpolation."""
    lo, hi = frame.W_low.to_numpy(float), frame.W_high.to_numpy(float)
    values = np.asarray(values, float)
    good = np.isfinite(values) & (values > 0)
    # A NaN separator preserves gaps and non-contiguous accepted W regions.
    xx = np.column_stack([lo[good], hi[good], np.full(good.sum(), np.nan)]).ravel()
    yy = np.column_stack([values[good], values[good], np.full(good.sum(), np.nan)]).ravel()
    return ax.plot(xx, yy, **kwargs)[0]


def plot_w_cross_section(config, reference_csv=None, model_csv=None):
    """Create ``upsilon_es_w_cross_section.pdf`` and companion audit/points CSV.

    Parameters
    ----------
    config : upsilon_esr.common.Config
        Explicit per-state run directories, energy and target luminosity.
    reference_csv : path, optional
        Defaults to the published ZEUS Table 1 values embedded above. A custom
        file requires a JSON provenance sidecar (same stem).
    model_csv : path, optional
        An independently calculated generator prediction, already averaged with
        the same flux and W bins. The saved MC truth baseline is retained in numerical audits only; no MC
        curves or entries are drawn in any panel.
    """
    set_style()
    states = ["1S", "2S", "3S"]
    reference, reference_info = _reference(reference_csv)
    sources, data, metadata, records = [], {}, {}, []
    audit = {
        "figure": "upsilon_es_w_cross_section", "energy": config.energy,
        "target_luminosity_invfb": float(config.lumi_invfb),
        "method": "Saved bin-by-bin cross sections; fixed corrections",
        "projection": "Ensemble-median expected Garwood limits on reconstructed counts",
        "central_cross_sections_changed": False,
        "background_systematics_or_efficiency_MC_included": False,
        "same_MC_closure_is_independent_validation": False,
        "reference": reference_info, "states": {},
    }
    for state in states:
        if config.reoptimization_run is not None:
            from .reoptimization_inputs import w_inputs
            frame, projection, metadata[state], case_sources, normalization_check = w_inputs(config, state)
            sources.extend(case_sources)
        else:
            artifacts = resolve_artifacts(config.runs[state])
            path = artifacts / "cross_sections" / f"{config.energy}_sigma_gamma_p_{state.lower()}.csv"
            normalization_path = path.with_name(path.stem + "_normalization.csv")
            frame = read_csv(path, required=(
                "W_low", "W_high", "W_center", "sigma_gamma_p_pb",
                "n_reco_weighted_visible", "n_true_weighted_visible", "flux_int",
                "sigma_gamma_p_truth_closure_pb", "reco_count_to_sigma_pb", "plotted",
            ), state=state).sort_values(["W_low", "W_high"]).reset_index(drop=True)
            normalization = read_csv(normalization_path, required=("br_ee",), state=state)
            metadata[state] = sample_metadata(config, state)
            normalization_check = _check_normalization(frame, normalization, state)
            projection = project_cross_section(
                frame, config.lumi_invfb, "sigma_gamma_p_pb",
                count_column="n_reco_weighted_visible", scale_column="reco_count_to_sigma_pb",
                coverage=config.coverage)
        valid = np.asarray(projection["valid"], bool) & _bool(frame.plotted)
        valid &= np.isfinite(frame.sigma_gamma_p_pb.to_numpy(float))
        valid &= frame.sigma_gamma_p_pb.to_numpy(float) > 0
        if not valid.any():
            raise ValueError(f"{state}: no supported W bins remain for the figure.")
        if model_csv is None:
            model = frame.sigma_gamma_p_truth_closure_pb.to_numpy(float)
            model_label = "MC"
            model_info = {"kind": "held_out_truth" if config.reoptimization_run is not None else "same_sample_truth", "column": "sigma_gamma_p_truth_closure_pb",
                          "independent_generator_evaluation": False}
        else:
            model, model_label, model_info = _model(model_csv, state, frame, config, metadata[state])
        if config.reoptimization_run is None:
            sources.extend([_fingerprint(path), _fingerprint(normalization_path)])
        audit["states"][state] = {
            "sample_metadata": metadata[state], "normalization_check": normalization_check,
            "model": model_info, "displayed_row_indices": np.flatnonzero(valid).tolist(),
            "omitted_row_indices": np.flatnonzero(~valid).tolist(),
            "omission_reason": "Upstream plotted=false, unsupported correction, or nonpositive cross section",
        }
        data[state] = (frame, projection, valid, model, model_label)
        for i, row in frame.iterrows():
            records.append({
                "state": state, "energy": config.energy, "W_low": row.W_low,
                "W_high": row.W_high, "W_center": row.W_center,
                "sigma_gamma_p_pb": row.sigma_gamma_p_pb,
                "expected_reconstructed_count": float(projection["mean_count"].iloc[i]),
                "interval_low_pb": float(projection["lower"].iloc[i]),
                "interval_high_pb": float(projection["upper"].iloc[i]),
                "target_lumi_invfb": config.lumi_invfb,
                "q2_max_GeV2": metadata[state]["q2_max"],
                "baseline_sigma_pb": float(model[i]), "displayed": bool(valid[i]),
                "rebinned": "rebinning" in metadata[state],
                "source_bins": projection["source_bins"].iloc[i] if "source_bins" in projection else str([i]),
                **({"overall_status": metadata[state]["overall_status"],
                    "status": metadata[state]["overall_status"],
                    "stage_status": metadata[state]["stage_status"], "publication_ready": False,
                    "interval_interpretation": metadata[state]["interval_interpretation"]}
                   if config.reoptimization_run is not None else {}),
            })

    fig = plt.figure(figsize=(7.1, 7.7), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.12, wspace=0.13)
    axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    for index, (ax, state) in enumerate(zip(axes, states)):
        frame, projected, valid, baseline, baseline_label = data[state]
        center = frame.W_center.to_numpy(float)
        xerr = np.vstack([center - frame.W_low, frame.W_high - center])
        reco_handle = draw_intervals(
            ax, center[valid], projected["mean"].to_numpy(float)[valid],
            projected["lower"].to_numpy(float)[valid],
            projected["upper"].to_numpy(float)[valid], xerr=xerr[:, valid],
            color=COLORS[state], fmt=MARKERS[state], ms=4.1, capsize=2,
            elinewidth=1.0, label="RECO", zorder=6,
        )
        upper = projected["upper"].to_numpy(float)[valid]
        positive = projected["mean"].to_numpy(float)[valid]
        ymax = float(np.nanmax(upper))
        xmax = float(frame.W_high.max())
        if state == "1S":
            yerr = np.vstack([reference.err_low_pb, reference.err_high_pb])
            reference_handle = ax.errorbar(reference.W_ref, reference.sigma_pb,
                        xerr=np.vstack([reference.W_ref - reference.W_low,
                                        reference.W_high - reference.W_ref]),
                        yerr=yerr, fmt="s", color="#555555", mfc="white", mew=1.1,
                        ms=5.1, capsize=2.5, elinewidth=1.1,
                        label="ZEUS", zorder=4)
            ymax = max(ymax, float(np.max(reference.sigma_pb + reference.err_high_pb)))
            xmax = max(xmax, float(reference.W_high.max()))
            annotate(ax, config)
            ax.text(-0.075, 1.03, "(a)", transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=FONT_SIZES["panel"])
            ax.legend([reco_handle, reference_handle],
                      ["RECO", "ZEUS"],
                      loc="lower right", fontsize=FONT_SIZES["legend"], frameon=False,
                      labelspacing=.3, handlelength=1.5)
        else:
            ax.text(0.04, 0.95, "(" + "abc"[index] + ")", transform=ax.transAxes,
                    ha="left", va="top", fontsize=FONT_SIZES["panel"])
        ax.text(0.97, 0.96, state_label(state), transform=ax.transAxes,
                ha="right", va="top", fontsize=FONT_SIZES["axis"], color=COLORS[state])
        if state == "1S":
            phase = rf"$Q^2<{metadata[state]['q2_max']:g}$ GeV$^2$"
            ax.text(0.97, 0.82, phase, transform=ax.transAxes,
                    ha="right", va="top", fontsize=FONT_SIZES["annotation"])
        ax.set_yscale("log")
        ax.set_xlim(max(0, float(frame.W_low.min()) - 2), xmax + 5)
        lower = projected["lower"].to_numpy(float)[valid]
        positive_limits = lower[np.isfinite(lower) & (lower > 0)]
        ymin = float(np.min(positive)) / 5.0
        if len(positive_limits):
            ymin = min(ymin, float(np.min(positive_limits)) / 1.4)
        ymin = max(np.finfo(float).tiny, ymin)
        # Reserve the upper portion for readable campaign/selection labels.
        ceiling = max(ymax * 2.4, ymin * (max(ymax * 1.10 / ymin, 1.0)) ** (1 / 0.59))
        ax.set_ylim(ymin, ceiling)
        # A finite log axis cannot display zero. Down arrows explicitly mark
        # zero lower confidence endpoints instead of implying a positive limit.
        for coordinate in center[valid][lower <= 0]:
            ax.annotate("", xy=(coordinate, ymin * 1.025), xytext=(coordinate, ymin * 1.5),
                        arrowprops={"arrowstyle": "-|>", "color": COLORS[state],
                                    "lw": .8, "mutation_scale": 6})
        ax.set_xlabel(r"$W_{\gamma p}$ [GeV]", loc="right")
        ax.set_ylabel(r"$\sigma_{\gamma p}$ [pb]", loc="top")
        ax.tick_params(which="both", direction="in", top=True, right=True)
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5 if index == 0 else 4))

    audit["source_files"] = sources
    audit["display_annotations"] = {
        "electron_id_labels_shown": False,
        "q2_labels_shown_for_states": ["1S"],
        "legend_order": ["RECO", "ZEUS"],
        "MC_curves_displayed": False,
        "RECO_drawn_above_ZEUS": True,
        "RECO": "Reconstructed signal-only projection with the saved bin-by-bin corrections; not uncorrected reconstructed counts.",
        "baseline": ("MC denotes the saved truth cross section from the same simulated sample."
                     if model_csv is None else "The independent model retains its supplied provenance label."),
        "ZEUS": "Published ZEUS 2009 reference; statistical and systematic errors combined as recorded in reference metadata.",
        "selection_metadata": "All states' Q2 and electron-ID settings remain in states[].sample_metadata.",
    }
    audit["notes"] = [
        "Cross sections and fluxes are read unchanged. Luminosity only changes the expected selected counts and projected statistical limits.",
        "Horizontal bars show W bin widths; ZEUS markers use published W0 reference energies.",
        "On the logarithmic y axes, downward arrows at the lower boundary denote zero lower count limits.",
        "The same-sample truth baseline is an algebraic closure reference, not an independent generator validation.",
        "Bin-by-bin corrections do not provide a full migration unfolding or a simultaneous three-state mass-fit covariance.",
    ]
    if config.reoptimization_run is not None:
        audit.update(method="Saved nonnegative response-matrix spectrum inference",
                     projection="Saved central toy sampling quantiles, not confidence limits; no luminosity rescaling",
                     input_run=str(config.reoptimization_run), publication_ready=False)
        audit["display_annotations"]["RECO"] = "Saved response-corrected expected spectrum; diagnostic if overall status is NOT_VALIDATED."
        audit["display_annotations"]["baseline"] = "MC denotes held-out truth; not an independent generator calculation."
        audit["notes"] = ["Saved bin-integrated photon flux conversions are used once; no Delta W division.",
                          "Horizontal bars show bin widths; point coordinates are arithmetic midpoints.",
                          "Sampling quantiles and covariance are read unchanged at the saved luminosity.",
                          "Final case statuses and failed checks are retained in each state's metadata and points CSV.",
                          "NOT_VALIDATED cases are diagnostics, not validated measurements."]
    audit["display_annotations"]["baseline_display"] = "Hidden in every panel; retained in numerical companion only."
    audit["rebinned_states"] = [s for s, m in metadata.items() if "rebinning" in m]
    if audit["rebinned_states"]:
        audit["central_cross_sections_changed"] = True
        audit["notes"].append("Requested adjacent bins are merged with integrated photon-flux weights; sampling limits come from replayed joint toys. Overall status is unchanged.")
    output = Path(config.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(output / "upsilon_es_w_cross_section_points.csv", index=False)
    reference.to_csv(output / "zeus_2009_w_reference.csv", index=False)
    write_json(output / "zeus_2009_w_reference.json", reference_info)
    paths = save_figure(fig, config, "upsilon_es_w_cross_section")
    audit["outputs"] = paths
    write_json(output / "upsilon_es_w_cross_section_audit.json", audit)
    return fig, audit
