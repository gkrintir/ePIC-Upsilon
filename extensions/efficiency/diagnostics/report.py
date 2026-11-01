"""Render descriptive efficiency diagnostics from counted MC CSV tables.

This module never changes selections, fits an efficiency, or rescales counts.
The curves are individual finite-MC Wilson intervals, not projected-data errors.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


COLORS = {"all": "#333333", "low": "#0072B2", "high": "#D55E00"}
STATE_COLORS = {"1S": "#0072B2", "2S": "#D55E00", "3S": "#009E73"}
VARIABLE_LABELS = {
    "true_pT": r"True $p_{T,\Upsilon}$ [GeV]",
    "true_W": r"True $W_{\gamma p}$ [GeV]",
    "true_y": r"True $y_{\Upsilon}$",
    "true_Q2": r"True $Q^2$ [GeV$^2$]",
}
STAGES = ("reco_mass", "reco_W", "sigma_epz", "proton_tag")
STAGE_COLORS = {"reco_mass": "#0072B2", "reco_W": "#D55E00", "sigma_epz": "#009E73", "proton_tag": "#CC79A7"}
MARKERS = {"1S": "o", "2S": "s", "3S": "^"}


def _read(path):
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _selection(frame, **values):
    result = frame
    for key, value in values.items():
        if key not in result:
            return result.iloc[0:0]
        result = result.loc[result[key] == value]
    return result


def _style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 11,
        "axes.labelsize": 12, "axes.titlesize": 11,
        "legend.fontsize": 9, "xtick.labelsize": 10,
        "ytick.labelsize": 10, "axes.linewidth": 1.0,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "savefig.dpi": 170, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def _q2_label(band, cap, split):
    if band == "all":
        return r"$0\leq Q^2\leq %g$" % cap
    if band == "low":
        return r"$0\leq Q^2\leq %g$" % min(cap, split)
    return r"$%g<Q^2\leq %g$" % (split, cap)


def _short_sample(sample):
    state = sample["state"]
    mode = str(sample.get("id_mode", "unspecified ID")).replace("_", " ")
    return r"$\Upsilon(%s)$, %s" % (state, mode)


def _efficiency(ax, frame, label, color, marker="o", sparse_threshold=50,
                xlo="low", xhi="high"):
    if frame.empty:
        return False
    frame = frame.sort_values(xlo)
    good = (frame.denominator.to_numpy(float) > 0)
    good &= np.isfinite(frame.efficiency.to_numpy(float))
    frame = frame.loc[good]
    if frame.empty:
        return False
    x = (frame[xlo].to_numpy(float) + frame[xhi].to_numpy(float)) / 2
    dx = (frame[xhi].to_numpy(float) - frame[xlo].to_numpy(float)) / 2
    y = frame.efficiency.to_numpy(float) * 100
    lo = frame.wilson_low.to_numpy(float) * 100
    hi = frame.wilson_high.to_numpy(float) * 100
    # Wilson intervals include the binomial MLE, including k=0 or k=n.
    if np.any(~np.isfinite(lo) | ~np.isfinite(hi) | (lo > y + 1e-10) | (hi < y - 1e-10)):
        raise ValueError("Malformed finite-MC efficiency intervals in report inputs")
    sparse = frame.denominator.to_numpy(float) < sparse_threshold
    ax.errorbar(x, y, xerr=dx, yerr=np.maximum(0, [y-lo, hi-y]),
                fmt="none", ecolor=color, elinewidth=.9, capsize=1.7, zorder=2)
    ax.plot(x, y, ls="none", marker=marker, ms=4.5, color=color,
            markerfacecolor="white", markeredgewidth=1.1, label=label, zorder=3)
    if np.any(~sparse):
        ax.plot(x[~sparse], y[~sparse], ls="none", marker=marker, ms=4.5,
                color=color, markeredgewidth=1.1, label="_nolegend_", zorder=4)
    return True


def _finish_axes(axes, ymin=0, ymax=None):
    for ax in np.asarray(axes, dtype=object).flat:
        ax.grid(axis="y", alpha=.15, lw=.6)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.set_ylim(bottom=ymin)
        if ymax is not None:
            ax.set_ylim(top=ymax)


def _save(fig, output, stem, title, config, footer):
    settings = config.get("settings", {})
    beam = config.get("energy_label", "%gx%g GeV" % (
        settings.get("electron_energy", 9), settings.get("proton_energy", 275)))
    partial = "  |  PARTIAL RUN" if config.get("execution_status_before_render") == "PARTIAL" else ""
    fig.suptitle("ePIC simulation diagnostics  |  %s%s\n%s" % (beam, partial, title + (" | proton tagged" if settings.get("proton_selection", "required") == "required" else " | proton inclusive")),
                 fontsize=13, fontweight="bold", y=.995)
    fig.text(.5, .014, footer, ha="center", va="bottom", fontsize=8, color="#444444")
    if config.get("synthetic", False):
        fig.text(.5, .045, "SYNTHETIC TEST INPUTS — NOT PHYSICS RESULTS",
                 ha="center", va="bottom", color="#9e1b1b", weight="bold", fontsize=10)
    fig.tight_layout(rect=(.01, .075 if config.get("synthetic", False) else .05, .99, .94))
    paths = []
    for suffix in ("pdf", "png"):
        target = output / "plots" / (stem + "." + suffix)
        fig.savefig(target, bbox_inches="tight")
        paths.append(str(target.relative_to(output)))
    plt.close(fig)
    return paths


def _q2_split(eff, sample, settings, output, config):
    split = settings["q2_split"]
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.9))
    for ax, variable in zip(axes, ("true_pT", "true_W", "true_y")):
        for band in ("all", "low", "high"):
            f = _selection(eff, sample=sample["key"], q2_band=band, variable=variable)
            _efficiency(ax, f, _q2_label(band, sample["q2_max"], split), COLORS[band],
                        sparse_threshold=settings["min_denominator"])
        ax.set_xlabel(VARIABLE_LABELS[variable])
        ax.set_ylabel("Selection efficiency [%]")
    axes[0].legend(frameon=False)
    _finish_axes(axes)
    return _save(fig, output, "1s_newid_q2_split", _short_sample(sample) + r": $Q^2$ split [GeV$^2$]", config,
                 "Individual finite-MC Wilson intervals. Hollow markers: fewer than %d eligible rows; no yield or luminosity rescaling."
                 % settings["min_denominator"])


def _id_comparison(eff, first, second, settings, output, config):
    split = settings["q2_split"]
    cases = [("low", min(first["q2_max"], split), min(second["q2_max"], split))]
    if np.isclose(first["q2_max"], second["q2_max"], rtol=0, atol=1e-12):
        if first["q2_max"] > split:
            cases.append(("all", first["q2_max"], second["q2_max"]))
    cases = [x for x in cases if np.isclose(x[1], x[2], rtol=0, atol=1e-12)]
    if not cases:
        return []
    fig, axes = plt.subplots(len(cases), 3, figsize=(12.6, 3.7*len(cases)), squeeze=False)
    for row, (band, cap, _) in enumerate(cases):
        for ax, variable in zip(axes[row], ("true_pT", "true_W", "true_y")):
            for sample, color, marker in ((first, "#0072B2", "o"), (second, "#D55E00", "s")):
                _efficiency(ax, _selection(eff, sample=sample["key"], q2_band=band, variable=variable),
                            str(sample["id_mode"]).replace("_", " "), color, marker,
                            settings["min_denominator"])
            ax.set_xlabel(VARIABLE_LABELS[variable])
            ax.set_ylabel("Selection efficiency [%]")
            ax.set_title(_q2_label("all", cap, split) + r" GeV$^2$")
        axes[row, 0].legend(frameon=False)
    _finish_axes(axes)
    return _save(fig, output, "1s_matched_q2_id_comparison", r"$\Upsilon(1S)$ ID comparison at matched cuts", config,
                 "Individual Wilson intervals; event pairing is not assumed. Differences also require matched production/eligible-row populations.")


def _state_comparison(eff, samples, settings, output, config):
    chosen = []
    for state in ("1S", "2S", "3S"):
        choices = [s for s in samples if s["state"] == state]
        if choices:
            chosen.append(next((s for s in choices if "new" in s.get("id_mode", "").lower()), choices[0]))
    if len(chosen) < 2:
        return []
    low_caps = [min(s["q2_max"], settings["q2_split"]) for s in chosen]
    rows = ["all"] + (["low"] if np.allclose(low_caps, low_caps[0], rtol=0, atol=1e-12) else [])
    fig, axes = plt.subplots(len(rows), 3, figsize=(12.5, 3.7*len(rows)), squeeze=False)
    for row, band in enumerate(rows):
        for ax, variable in zip(axes[row], ("true_pT", "true_W", "true_y")):
            for sample in chosen:
                label = _short_sample(sample)
                if band == "all":
                    label += r", $Q^2\leq%g$" % sample["q2_max"]
                _efficiency(ax, _selection(eff, sample=sample["key"], q2_band=band, variable=variable),
                            label, STATE_COLORS[sample["state"]], MARKERS[sample["state"]],
                            settings["min_denominator"])
            ax.set_xlabel(VARIABLE_LABELS[variable])
            ax.set_ylabel("Selection efficiency [%]")
            ax.set_title("Nominal sample selections" if band == "all" else r"Common $Q^2\leq%g$ GeV$^2$" % low_caps[0])
        axes[row, 0].legend(frameon=False, fontsize=8)
    _finish_axes(axes)
    return _save(fig, output, "states_nominal_and_common_q2", "State comparisons with selections retained", config,
                 "Q² limits are in GeV². Matching Q² alone does not match electron ID or W/rapidity populations. Hollow markers indicate limited eligible MC.")


def _cutflow(cutflow, nminus1, sample, settings, output, config):
    source = _selection(cutflow, sample=sample["key"])
    if source.empty:
        return []
    bands = [b for b in ("all", "low", "high") if b in set(source.q2_band)]
    fig, axes = plt.subplots(2, len(bands), figsize=(4.2*len(bands), 7.1), squeeze=False)
    stage_names = {"reco_mass": "Mass 8–12 GeV", "reco_W": "+ reconstructed W", "sigma_epz": r"+ $\Sigma(E-p_z)$", "proton_tag": "+ proton policy"}
    isolated_names = {"reco_mass": "Mass cut", "reco_W": "W cut", "sigma_epz": r"$\Sigma(E-p_z)$ cut", "proton_tag": "Proton policy"}
    for j, band in enumerate(bands):
        for stage in STAGES:
            f = _selection(source, q2_band=band, stage=stage).copy()
            if not f.empty:
                f["denominator"] = f.baseline
                _efficiency(axes[0, j], f, stage_names[stage], STAGE_COLORS[stage],
                            sparse_threshold=settings["min_denominator"], xlo="pt_low", xhi="pt_high")
            nm = _selection(nminus1, sample=sample["key"], q2_band=band, cut=stage)
            baseline = _selection(source, q2_band=band, stage="truth")
            if not nm.empty and not baseline.empty:
                m = nm.merge(baseline[["pt_low", "pt_high", "baseline"]], on=["pt_low", "pt_high"], validate="one_to_one")
                m = m.loc[m.baseline > 0].sort_values("pt_low")
                # This is a descriptive count fraction, not a fitted estimate or significance.
                axes[1, j].plot((m.pt_low + m.pt_high)/2, 100*m.failed_without_cut/m.baseline,
                                marker="o", ms=3.5, lw=1, color=STAGE_COLORS[stage], label=isolated_names[stage])
        axes[0, j].set_title(_q2_label(band, sample["q2_max"], settings["q2_split"]) + r" GeV$^2$")
        axes[0, j].set_ylabel("Cumulative efficiency [%]")
        axes[1, j].set_ylabel("Fail this cut only / eligible [%]")
        axes[1, j].set_xlabel(VARIABLE_LABELS["true_pT"])
    axes[0, 0].legend(frameon=False)
    axes[1, 0].legend(frameon=False)
    _finish_axes(axes)
    return _save(fig, output, sample["key"] + "_selection_cutflow", _short_sample(sample) + ": stored-tuple selections", config,
                 "Top: cumulative cuts with Wilson intervals. Bottom: descriptive count fractions passing all other final cuts; overlapping failures are excluded.")


def _grid(frame, column):
    xedges = np.unique(np.r_[frame.x_low.to_numpy(float), frame.x_high.to_numpy(float)])
    yedges = np.unique(np.r_[frame.y_low.to_numpy(float), frame.y_high.to_numpy(float)])
    z = np.full((len(yedges)-1, len(xedges)-1), np.nan)
    xi = np.searchsorted(xedges, frame.x_low.to_numpy(float))
    yi = np.searchsorted(yedges, frame.y_low.to_numpy(float))
    z[yi, xi] = frame[column].to_numpy(float)
    return xedges, yedges, z


def _maps(maps, sample, settings, output, config):
    source = _selection(maps, sample=sample["key"], q2_band="all")
    if source.empty:
        return []
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.3))
    for j, variable in enumerate(("true_W", "true_y", "true_Q2")):
        f = _selection(source, y_variable=variable)
        if f.empty:
            continue
        xe, ye, value = _grid(f, "efficiency")
        _, _, den = _grid(f, "denominator")
        value = np.ma.masked_where((den <= 0) | ~np.isfinite(value), value*100)
        counts = np.ma.masked_where((den <= 0) | ~np.isfinite(den), den)
        im = axes[0, j].pcolormesh(xe, ye, value, norm=Normalize(0, 100), cmap="viridis", shading="flat")
        fig.colorbar(im, ax=axes[0, j], label="Selection efficiency [%]", fraction=.046, pad=.03)
        maximum = max(1.01, float(np.nanmax(den)))
        im = axes[1, j].pcolormesh(xe, ye, counts, norm=LogNorm(1, maximum), cmap="magma", shading="flat")
        fig.colorbar(im, ax=axes[1, j], label="Eligible MC rows", fraction=.046, pad=.03)
        sparse = (den > 0) & (den < settings["min_denominator"])
        xx, yy = np.meshgrid((xe[:-1]+xe[1:])/2, (ye[:-1]+ye[1:])/2)
        axes[0, j].scatter(xx[sparse], yy[sparse], marker="x", s=9, color="white", linewidth=.5)
        for ax in axes[:, j]:
            ax.set_xlabel(VARIABLE_LABELS["true_pT"])
            ax.set_ylabel(VARIABLE_LABELS[variable])
            if variable == "true_Q2":
                positive = ye[ye > 0]
                ax.set_yscale("symlog", linthresh=float(positive[0]) if len(positive) else 1e-6)
                ax.axhline(settings["q2_split"], color="#009E73", lw=1.2, ls="--")
    return _save(fig, output, sample["key"] + "_pt_correlations", _short_sample(sample) + ": efficiency and population maps", config,
                 "Nominal Q² range. White crosses: fewer than %d eligible rows. Blank: no eligible rows. Green line marks the Q² split; all band maps are saved in CSV."
                 % settings["min_denominator"])


def _composition(composition, sample, settings, output, config):
    source = _selection(composition, sample=sample["key"], q2_band="all")
    if source.empty:
        return []
    regions = sorted(set(zip(source.pt_low, source.pt_high)))
    # Show low pT, the dip region, and the final tail region without assuming exact edges.
    if len(regions) > 3:
        regions = [regions[0], regions[len(regions)//2], regions[-1]]
    fig, axes = plt.subplots(len(regions), 3, figsize=(12.5, 2.9*len(regions)), squeeze=False)
    for i, (low, high) in enumerate(regions):
        for j, variable in enumerate(("true_W", "true_y", "true_Q2")):
            f = _selection(source, pt_low=low, pt_high=high, variable=variable).sort_values("low")
            ax = axes[i, j]
            if not f.empty:
                edges = np.r_[f.low.to_numpy(float), f.high.iloc[-1]]
                ax.stairs(f.truth_fraction.to_numpy(float), edges, color="#444444", lw=1.6, label="Eligible")
                ax.stairs(f.selected_fraction.to_numpy(float), edges, color="#0072B2", lw=1.6, label="Selected")
                finite = f[["truth_fraction", "selected_fraction"]].to_numpy(float)
                finite = finite[np.isfinite(finite)]
                ax.set_ylim(0, max(.01, 1.4*float(finite.max())) if finite.size else 1)
                ax.text(.97, .94, "Eligible %d\nSelected %d" % (f.denominator.sum(), f.selected.sum()),
                        transform=ax.transAxes, ha="right", va="top", fontsize=8)
                if variable == "true_Q2":
                    positive = edges[edges > 0]
                    ax.set_xscale("symlog", linthresh=float(positive[0]) if len(positive) else 1e-6)
                    ax.axvline(settings["q2_split"], color="#009E73", ls="--", lw=1)
            ax.set_xlabel(VARIABLE_LABELS[variable])
            ax.set_ylabel("Fraction in displayed bin")
            last_region = (low, high) == max(set(zip(source.pt_low, source.pt_high)))
            relation = r"\leq" if last_region else "<"
            ax.set_title(r"$%g\leq p_{T,\Upsilon}%s%g$ GeV" % (low, relation, high), fontsize=10)
    axes[0, 0].legend(frameon=False, loc="upper left")
    _finish_axes(axes)
    return _save(fig, output, sample["key"] + "_tail_composition", _short_sample(sample) + ": distributions within pT regions", config,
                 "Eligible and selected distributions each normalized within the displayed variable range; these are bin fractions, not densities. See audited under/overflows.")


def _tail_counts(eff, samples, settings, output, config):
    if not samples:
        return []
    fig, axes = plt.subplots(2, len(samples), figsize=(6.2*len(samples), 7.2), squeeze=False)
    for j, sample in enumerate(samples):
        for band in ("all", "low", "high"):
            f = _selection(eff, sample=sample["key"], q2_band=band, variable="true_pT")
            f = f.loc[f.low >= 1.6-1e-10].sort_values("low")
            _efficiency(axes[0, j], f, _q2_label(band, sample["q2_max"], settings["q2_split"]),
                        COLORS[band], sparse_threshold=settings["min_denominator"])
            if not f.empty:
                x = (f.low+f.high)/2
                axes[1, j].plot(x, f.denominator, color=COLORS[band], marker="o", ms=4, lw=1.2)
                axes[1, j].plot(x, f.selected, color=COLORS[band], marker="s", ms=3.5, lw=1.2, ls="--")
        axes[0, j].set_title(_short_sample(sample))
        axes[0, j].set_ylabel("Selection efficiency [%]")
        axes[0, j].legend(frameon=False)
        axes[1, j].set_xlabel(VARIABLE_LABELS["true_pT"])
        axes[1, j].set_ylabel("Raw MC rows")
        axes[1, j].set_yscale("symlog", linthresh=1)
        axes[1, j].legend(handles=[Line2D([], [], color="black", marker="o", label="Eligible"),
                                          Line2D([], [], color="black", marker="s", ls="--", label="Selected")], frameon=False)
    _finish_axes(axes)
    return _save(fig, output, "1s_high_pt_counts", "High-pT efficiencies and the underlying MC counts", config,
                 "Colors identify Q² ranges in both rows. Counts are unweighted simulated rows; the count axis is linear near zero and logarithmic above one.")


def _mixture(mixture, sample, output, config):
    source = _selection(mixture, sample=sample["key"], variable="true_pT").sort_values("low") if not mixture.empty else mixture
    if source.empty:
        return []
    f = source.loc[source.n_all > 0]
    if f.empty:
        return []
    fig, axes = plt.subplots(1, 2, figsize=(12.3, 4.2))
    x = (f.low + f.high)/2
    axes[0].plot(x, 100*f.high_truth_fraction, color="#444444", marker="o", ms=4, lw=1.2, label="Among eligible rows")
    axes[0].plot(x, 100*f.high_selected_fraction, color="#D55E00", marker="s", ms=4, lw=1.2, label="Among selected rows")
    axes[0].set_ylabel("Higher-Q² fraction [%]")
    axes[0].legend(frameon=False)
    axes[1].plot(x, 100*f.eff_all, color="#333333", marker="o", ms=4, lw=1.4, label="Total efficiency")
    axes[1].plot(x, 100*f.low_contribution, color=COLORS["low"], marker="s", ms=3.5, lw=1.2, ls="--", label="Low-Q² contribution")
    axes[1].plot(x, 100*f.high_contribution, color=COLORS["high"], marker="^", ms=3.5, lw=1.2, ls="--", label="Higher-Q² contribution")
    axes[1].set_ylabel("Contribution / all eligible rows [%]")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.set_xlabel(VARIABLE_LABELS["true_pT"])
    _finish_axes(axes)
    axes[0].set_ylim(0, 100)
    return _save(fig, output, "1s_newid_q2_mixture", _short_sample(sample) + ": exact Q²-mixture decomposition", config,
                 "Descriptive fractions from raw counts. Low + higher contributions equal total efficiency; conditional intervals are in the Q²-split plot.")


def _mass_failures(reasons, sample, settings, output, config):
    source = _selection(reasons, sample=sample["key"], cut="reco_mass")
    if source.empty:
        return []
    bands = [b for b in ("all", "low", "high") if b in set(source.q2_band)]
    fig, axes = plt.subplots(1, len(bands), figsize=(4.3*len(bands), 4.1), squeeze=False)
    categories = (("Unavailable/sentinel mass", ("nonfinite", "minus_999_sentinel"), "#CC79A7", "o"),
                  ("Finite mass outside range", ("finite_lower_rejection", "finite_upper_rejection"), "#D55E00", "s"))
    for ax, band in zip(axes.flat, bands):
        f = _selection(source, q2_band=band)
        groups = f.groupby(["pt_low", "pt_high"])
        if (groups.denominator.nunique() > 1).any():
            raise ValueError("Mass-failure categories do not share a common eligible denominator")
        baseline = groups.denominator.first()
        for label, included, color, marker in categories:
            counts = f.loc[f.reason.isin(included)].groupby(["pt_low", "pt_high"])["count"].sum()
            merged = pd.concat([baseline.rename("denominator"), counts.rename("count")], axis=1).reset_index()
            merged = merged.loc[merged.denominator > 0].sort_values("pt_low")
            ax.plot((merged.pt_low+merged.pt_high)/2, 100*merged["count"]/merged.denominator,
                    color=color, marker=marker, ms=4, lw=1.2, label=label)
        ax.set_title(_q2_label(band, sample["q2_max"], settings["q2_split"]) + r" GeV$^2$")
        ax.set_xlabel(VARIABLE_LABELS["true_pT"])
        ax.set_ylabel("Mass-failure category / eligible [%]")
    axes[0, 0].legend(frameon=False, fontsize=8.4)
    _finish_axes(axes, ymax=100)
    return _save(fig, output, "1s_newid_mass_failure_modes", _short_sample(sample) + ": distinguish unavailable mass from finite out-of-range mass", config,
                 "Descriptive count fractions before the other reconstructed cuts. Unavailable = nonfinite or exactly −999; finite out-of-range values do not establish physical mass tails.")


def _summary(output, config, samples, audits, eff):
    settings = config.get("settings", {})
    execution = config.get("execution_status_before_render", "UNKNOWN")
    lines = ["# Efficiency diagnostic summary", "", "Execution status before rendering: **%s**." % execution,
             "Physics interpretation status: **GENERATED_DIAGNOSTICS_NOT_CAUSAL_VALIDATION**", "",
             "The outputs count eligible stored tuple rows. They do not establish the complete generated-event denominator, independently verify the producer ID decisions, or validate full experimental efficiencies.", "",
             "## Requested inputs and reference checks", "",
             "Every requested sample is listed, including failed inputs. A partial execution is not a complete comparison; a reference MISMATCH requires inspection before using the corresponding curves.", "",
             "| Requested sample | Processing | Saved efficiency reference | Detail |",
             "|---|---|---|---|"]
    for sample in samples:
        audit = audits.get(sample["key"], {})
        reference = audit.get("reference_efficiency_check", {})
        ref_status = reference.get("status", "NOT_RUN" if audit.get("status") == "FAILED" else "NOT_CONFIGURED")
        detail = str(audit.get("error", reference.get("reason", ""))).replace("|", "\\|").replace("\n", " ")
        lines.append("| %s | %s | %s | %s |" % (sample.get("label", sample["key"]), audit.get("status", "UNKNOWN"), ref_status, detail))
    lines += ["", "Input/reference evidence: [sample_audits.json](sample_audits.json). Matched-input checks: [input_comparability.json](input_comparability.json). These checks are not proof of event pairing or identical producer configurations.", "",
             "## Selection and uncertainty contract", "",
             "- Truth: finite state mass within the configured window, nonnegative true Q² within the sample/band limit, and true W in the configured range.",
             "- Selected: truth selection plus finite reconstructed mass in 8–12 GeV, reconstructed W in the configured range, and |Σ(E−pz)−2Ee| ≤ 0.2Ee.",
             "- Efficiencies and Wilson intervals use unweighted MC counts. No luminosity, cross-section, branching-fraction, flux or event-weight rescaling is applied.",
             "- ID comparisons use identical selected Q² and W cuts, but differences in production or stored-row populations may remain. Event pairing is not assumed; no paired significance or ratio uncertainty is reported.",
             "- Cutflow stages are analysis selections on stored branches. They do not isolate tracking, cluster matching, tagger or producer-ID failures individually.", "",
             "- A failure of the 8–12 GeV reconstructed-mass requirement includes unavailable/sentinel values. It must not be interpreted directly as a physical mass-window loss. `1s_newid_mass_failure_modes` separates nonfinite/exactly −999 values from finite out-of-range values; even finite out-of-range values do not by themselves establish physical mass tails. All cut-specific categories are in `tables/failure_reasons.csv`.", "",
             "## Counts in the displayed pT range", "",
             "The following totals include only finite, in-range pT values. Underflows, overflows and invalid entries remain in the per-sample audit.", "",
             "| Sample | Q² band | Eligible rows | Selected rows | Selected/eligible |",
             "|---|---|---:|---:|---:|"]
    for sample in samples:
        for band in ("all", "low", "high"):
            f = _selection(eff, sample=sample["key"], q2_band=band, variable="true_pT")
            if f.empty:
                continue
            n, k = int(f.denominator.sum()), int(f.selected.sum())
            value = "undefined" if n == 0 else "%.3f%%" % (100*k/n)
            lines.append("| %s | %s | %d | %d | %s |" % (sample.get("label", sample["key"]), band, n, k, value))
    lines += ["", "## Interpretation sequence", "",
              "1. Inspect `1s_newid_q2_split`: determine whether the high-pT recovery occurs within the low-Q² population, the added higher-Q² population, or their changing mixture.",
              "   `1s_newid_q2_mixture` displays the higher-Q² eligible/selected fractions and exact contribution identity. Its full numerical decomposition is saved in `q2_mixture.csv`; these descriptive fractions do not add independent evidence beyond the original counts.",
              "2. Inspect `1s_matched_q2_id_comparison`: compare the two 1S reconstructions at identical Q² and W cuts. Check truth-distribution and stored-row differences before attributing every change to ID.",
              "3. Inspect `1s_high_pt_counts` and each sample's selection cutflow: distinguish a stable pattern from limited MC support and identify which stored analysis requirement rejects events.",
              "4. Inspect the pT correlation maps, tail-composition figures and `strata.csv`: test whether W/rapidity composition changes explain the one-dimensional trend. Maps for all/low/high bands are retained in `maps.csv`.",
              "5. A causal ID or detector explanation requires matched production and, where necessary, producer-level association/cutflow information. These diagnostics do not automatically declare a cause.", "",
              "## Reporting limits", "",
              "- Hollow markers and crossed map cells denote fewer than %d eligible rows. This is a display flag, not a universal statistical-validity criterion." % settings.get("min_denominator", 50),
              "- Zero selected events with a positive denominator have a nonzero Wilson upper bound. No eligible rows give undefined efficiency and remain blank.",
              "- Nominal and low-Q² samples overlap; all = low + high in counts. Curves are not statistically independent.",
              "- Composition plots show separately normalized bin fractions. They are not probability densities when bin widths differ.",
              "- A numerical trend seen here is a diagnostic observation, not proof of a detector-performance mechanism.", ""]
    if config.get("synthetic", False):
        lines.insert(5, "**SYNTHETIC TEST INPUTS — NOT PHYSICS RESULTS.**\n")
    path = output / "summary.md"
    path.write_text("\n".join(lines))
    return str(path.relative_to(output))


def render_report(output_dir, config, samples, sample_audits):
    """Create PNG/PDF diagnostics and a Markdown summary from saved tables."""
    output = Path(output_dir)
    (output / "plots").mkdir(parents=True, exist_ok=True)
    settings = dict(config.get("settings", {}))
    settings.setdefault("q2_split", .01)
    settings.setdefault("min_denominator", 50)
    _style()
    tables = {name: _read(output/"tables"/(name+".csv"))
              for name in ("efficiency", "cutflow", "nminus1", "maps", "composition", "q2_mixture", "failure_reasons")}
    eff = tables["efficiency"]
    if eff.empty:
        raise ValueError("No efficiency table rows are available for plotting")
    requested_samples = list(samples)
    present = set(eff["sample"])
    samples = [s for s in requested_samples if s["key"] in present]
    one_s = [s for s in samples if s["state"] == "1S"]
    new = next((s for s in one_s if "new" in str(s.get("id_mode", "")).lower() or "newid" in s["key"].lower()), None)
    standard = next((s for s in one_s if s is not new), None)
    paths = []
    if new is not None:
        paths += _q2_split(eff, new, settings, output, config)
        paths += _mixture(tables["q2_mixture"], new, output, config)
        paths += _mass_failures(tables["failure_reasons"], new, settings, output, config)
    if new is not None and standard is not None:
        paths += _id_comparison(eff, new, standard, settings, output, config)
    paths += _state_comparison(eff, samples, settings, output, config)
    for sample in one_s:
        paths += _cutflow(tables["cutflow"], tables["nminus1"], sample, settings, output, config)
        paths += _maps(tables["maps"], sample, settings, output, config)
        paths += _composition(tables["composition"], sample, settings, output, config)
    paths += _tail_counts(eff, one_s, settings, output, config)
    paths.append(_summary(output, config, requested_samples, sample_audits, eff))
    return paths
