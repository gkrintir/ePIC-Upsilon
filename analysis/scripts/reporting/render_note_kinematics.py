#!/usr/bin/env python3
"""Render the reviewed note's split kinematic pages from prepared single-state Upsilon tuples.

Single-sample projected pages retain notebook 03's broad mass/W>=0 selection.
Two-sample projected pages and t-method pages additionally apply exclusivity
and each sample's requested reconstructed W range.
Stored electron-ID and cut-flow histograms are producer counts, not tuple yields.
MC mode reads weights as supplied. Projected mode scales weighted 1D yields
by the target/source luminosity ratio, never by their post-selection sum.
"""
import argparse
import csv
import json
import sys
import warnings
from pathlib import Path

# Shared explicit reconstructed-proton selection; truth denominators are unchanged.
from pathlib import Path as _ProtonPath
import sys as _proton_sys
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from note_style import set_style, finish
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projected_intervals import (DEFAULT_COVERAGE, add_interval_arguments,
                                 expected_count_intervals, draw_absolute_intervals)

COLORS = ["#0072B2", "#E69F00"]
STYLES = ["-", "--"]
MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
MASS_1S = MASSES["1S"]  # Legacy import compatibility.


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_sample(path, label, ee, wmin, wmax, args):
    state = getattr(args, "state", "1S")
    mass = MASSES[state]
    with uproot.open(path) as source:
        tree = source[args.tree]
        data = {name: tree[name].array(library="np").astype(float) for name in tree.keys()}
        stored = {}
        for name in ["h_CutFlow", "h_EleID_EpRatio", "h_EleID_IsoFraction"]:
            if name in source:
                h = source[name]
                stored[name] = {"values": h.values(), "edges": h.axis().edges(),
                                "labels": h.axis().labels(), "flow": h.values(flow=True),
                                "variance": h.variances()}
    required = ["weight", "weight_proj", "true_mass", "true_Q2", "reco_mass", "reco_W", "sigma_epz"]
    required += proton_branches()
    missing = [name for name in required if name not in data]
    if missing:
        raise ValueError(f"{path}: missing branches {missing}")
    weight = data["weight"]
    if args.stat_mode == "projected":
        validate_projection_weights(weight, data["weight_proj"], path)
    truth = np.isfinite(weight) & (weight >= 0) & np.isfinite(data["true_mass"])
    truth &= np.abs(data["true_mass"] - mass) <= .2
    truth &= np.isfinite(data["true_Q2"]) & (data["true_Q2"] >= 0) & (data["true_Q2"] <= args.q2_cut)
    projected = truth & np.isfinite(data["reco_mass"]) & (data["reco_mass"] >= 8) & (data["reco_mass"] <= 12)
    projected &= np.isfinite(data["reco_W"]) & (data["reco_W"] >= 0)
    projected &= proton_tag_mask(data)
    selected = projected & np.isfinite(data["sigma_epz"]) & (np.abs(data["sigma_epz"]-2*ee) <= .2*ee)
    selected &= (data["reco_W"] >= wmin) & (data["reco_W"] <= wmax)
    return {"path": str(Path(path).resolve()), "label": label, "ee": ee, "wmin": wmin,
            "wmax": wmax, "state": state, "mass": mass, "data": data, "stored": stored,
            "truth": truth, "projected": projected, "selected": selected}


def edges(lo, hi, n=50):
    return np.linspace(lo, hi, n+1)


def panel(branch, xlabel, lo, hi, n=50, log=False, weight="weight"):
    return {"branch": branch, "xlabel": xlabel, "edges": edges(lo, hi, n),
            "log": log, "weight": weight}


def validate_projection_weights(weight, corrected_weight, context="input"):
    """Require a positive underlying event intensity for projected data errors."""
    if np.any(~np.isfinite(weight)) or np.any(weight < 0):
        raise ValueError(f"{context}: projected statistics require finite nonnegative weight; signed templates need a separate model")
    if np.any(~np.isfinite(corrected_weight)) or np.any(corrected_weight < 0):
        raise ValueError(f"{context}: projected statistics require finite nonnegative weight_proj")
    if np.any((corrected_weight > 0) & (weight <= 0)):
        raise ValueError(f"{context}: nonzero weight_proj requires positive underlying weight")


def histogram_statistics(values, weight, base_weight, bins, stat_mode="mc",
                         luminosity_ratio=1.0, corrected=False, raw_count=False,
                         count_interval="median_garwood", interval_coverage=DEFAULT_COVERAGE):
    """Single-fill histograms with separate source-MC and projected variances.

    Each MC row represents a projected Poisson intensity r*w. An absolute
    corrected weight w_corr supplies an event correction c=w_corr/w, hence
    E[Y]=r*sum(w_corr), Var_data(Y)=r*sum(w*c**2). The finite-MC variance of
    the scaled template remains r**2*sum(w_corr**2), and is not added here.
    """
    sums = np.histogram(values, bins=bins, weights=weight)[0]
    sumw2 = np.histogram(values, bins=bins, weights=weight**2)[0]
    counts = np.histogram(values, bins=bins)[0]
    scale = luminosity_ratio if stat_mode == "projected" and not raw_count else 1.0
    projected_variance = None
    if raw_count:
        variance = sumw2
        model = "raw_mc_count"
    elif stat_mode == "projected":
        validate_projection_weights(base_weight, weight if corrected else np.zeros_like(weight))
        if corrected:
            contribution = np.divide(weight**2, base_weight,
                                     out=np.zeros_like(weight), where=base_weight > 0)
            projected_variance = scale*np.histogram(values, bins=bins, weights=contribution)[0]
            model = "projected_corrected_poisson"
        else:
            projected_variance = scale*sums
            model = "projected_poisson"
        variance = projected_variance
    else:
        variance = sumw2
        model = "finite_mc_sumw2"
    result = {"source_sum_weight": sums, "source_sum_weight2": sumw2, "entries": counts,
              "yield": scale*sums, "variance": variance,
              "finite_mc_variance": scale**2*sumw2,
              "projected_data_variance": projected_variance,
              "luminosity_scale": scale, "uncertainty_model": model}
    result.update(histogram_intervals(values, weight, base_weight, bins, result,
                                     stat_mode, corrected, raw_count,
                                     count_interval, interval_coverage))
    return result


def histogram_intervals(values, weight, base_weight, bins, stats, stat_mode,
                        corrected, raw_count, method, coverage):
    """Keep confidence-bound forecasts separate from means and variances.

    Garwood endpoints are evaluated at the median integer Poisson observation,
    not at the fractional expected count. A corrected yield has a scaled
    Poisson interval only when its event correction is constant within the bin.
    Variable corrections retain their exact compound-Poisson variance but have
    no count confidence interval in this renderer.
    """
    mean = stats["yield"]
    sd = np.sqrt(stats["variance"])
    n = len(mean)
    result = {"covariance_sd": sd, "interval_lower": mean-sd, "interval_upper": mean+sd,
              "interval_reference_count": np.full(n, np.nan),
              "interval_method": np.full(n, "mc_standard_deviation", dtype=object),
              "interval_coverage": None, "interval_correction_factor": np.full(n, np.nan),
              "interval_underlying_mean": np.full(n, np.nan), "interval_warnings": []}
    if stat_mode != "projected" or raw_count:
        return result
    result["interval_coverage"] = coverage
    if not corrected:
        interval = expected_count_intervals(mean, coverage=coverage, method=method)
        result.update(interval_lower=interval["lower"], interval_upper=interval["upper"],
                      interval_reference_count=interval["reference_count"],
                      interval_method=np.full(n, interval["method"], dtype=object),
                      interval_correction_factor=np.ones(n), interval_underlying_mean=mean.copy())
        return result
    if method == "standard_deviation":
        unit = expected_count_intervals(np.ones(n), coverage=coverage, method=method)
        halfwidth = (unit["upper"]-1)*sd
        result.update(interval_lower=mean-halfwidth, interval_upper=mean+halfwidth,
                      interval_method=np.full(n, "compound_poisson_standard_deviation", dtype=object))
        return result

    # Ignore zero-intensity MC rows when deciding whether the event correction
    # is constant. The last histogram edge is included, matching np.histogram.
    positive = base_weight > 0
    correction = np.divide(weight, base_weight, out=np.full_like(weight, np.nan), where=positive)
    indices = np.searchsorted(bins, values, side="right")-1
    indices = np.where(values == bins[-1], n-1, indices)
    factors = np.full(n, np.nan)
    reasons = np.full(n, "unavailable_variable_correction", dtype=object)
    for j in range(n):
        local = correction[positive & (indices == j)]
        if len(local) and np.allclose(local, local[0], rtol=1e-10, atol=1e-12):
            factors[j] = local[0]
        elif not len(local):
            # No corrected events determine a correction in this bin. Do not
            # infer it from other bins, even if their factors happen to agree.
            reasons[j] = "unavailable_empty_bin_correction"
    available = np.isfinite(factors)
    underlying = stats["luminosity_scale"]*np.histogram(values, bins=bins, weights=base_weight)[0]
    lower, upper, reference = (np.full(n, np.nan) for _ in range(3))
    if np.any(available):
        interval = expected_count_intervals(underlying[available], coverage=coverage, method=method)
        lower[available] = factors[available]*interval["lower"]
        upper[available] = factors[available]*interval["upper"]
        reference[available] = interval["reference_count"]
        reasons[available] = str(interval["method"])+"_scaled_constant_correction"
    result.update(interval_lower=lower, interval_upper=upper,
                  interval_reference_count=reference, interval_method=reasons,
                  interval_correction_factor=factors, interval_underlying_mean=underlying)
    if np.any(~available):
        result["interval_warnings"].append(
            "Count confidence bounds unavailable in zero-based bins "
            + ",".join(map(str, np.flatnonzero(~available)))
            + ": variable or unknown empty-bin correction; the compound-Poisson variance is retained. "
              "Use standard_deviation for a variance-scale display, not a Poisson confidence interval.")
    return result


def arrays(sample, spec, mask_name):
    data = sample["data"]
    branch = spec["branch"]
    if branch == "reco_pT2":
        values = data["reco_pT"]**2
    else:
        values = data[branch]
    weight = np.ones(len(values), dtype=float) if spec["weight"] == "mc_count" else data[spec["weight"]]
    mask = sample[mask_name] & np.isfinite(values) & (values > -900) & np.isfinite(weight) & (weight >= 0)
    # Projective corrected yields retain their supplied absolute weights and subset.
    if spec["weight"] == "weight_proj":
        mask &= weight > 0
    return values[mask], weight[mask], data["weight"][mask]


def draw_pages(samples, specs, stem, title, subtitle, mask_name, outdir, rows, audit,
               density=False, stat_mode="mc", source_lumi_invfb=10.0,
               projection_lumi_invfb=None, count_interval="median_garwood",
               interval_coverage=DEFAULT_COVERAGE):
    target_lumi = source_lumi_invfb if projection_lumi_invfb is None else projection_lumi_invfb
    luminosity_ratio = target_lumi/source_lumi_invfb
    nrows = 1 if len(specs) <= 2 else 2
    fig, axs = plt.subplots(nrows, 2, figsize=(12.4, 5.6 if nrows == 1 else 9.8), squeeze=False)
    for ax, spec in zip(axs.flat, specs):
        positive = []
        if spec["log"]:
            ax.set_yscale("log")
        for i, sample in enumerate(samples):
            if spec["branch"] != "reco_pT2" and spec["branch"] not in sample["data"]:
                raise ValueError(f"{sample['path']}: missing {spec['branch']}")
            values, weight, base_weight = arrays(sample, spec, mask_name)
            bins = spec["edges"]
            raw_count = spec["weight"] == "mc_count"
            stats = histogram_statistics(values, weight, base_weight, bins, stat_mode,
                                         luminosity_ratio, spec["weight"] == "weight_proj", raw_count,
                                         count_interval, interval_coverage)
            sums, sumw2, counts = stats["source_sum_weight"], stats["source_sum_weight2"], stats["entries"]
            display_divisor = np.diff(bins) if density else np.ones(len(sums))
            display = stats["yield"]/display_divisor
            error = np.sqrt(stats["variance"])/display_divisor
            lower = stats["interval_lower"]/display_divisor
            upper = stats["interval_upper"]/display_divisor
            ax.stairs(display, bins, color=COLORS[i], linestyle=STYLES[i], lw=1.7)
            centers = .5*(bins[1:]+bins[:-1])
            if stat_mode == "projected" and not raw_count:
                draw_absolute_intervals(ax, centers, display, lower, upper,
                                        fmt="o", markersize=2.2, color=COLORS[i], alpha=.7, linewidth=.7)
                positive.extend(upper[np.isfinite(upper) & (upper > 0)])
            else:
                good = display > 0 if spec["log"] else np.ones(len(display), dtype=bool)
                ax.errorbar(centers[good], display[good], yerr=error[good], fmt="none", color=COLORS[i], lw=.7, alpha=.7)
            positive.extend(display[display > 0])
            for warning in stats["interval_warnings"]:
                warnings.warn(f"{stem}/{sample['label']}/{spec['branch']}: {warning}", RuntimeWarning)
            if stats["interval_warnings"]:
                ax.text(.98, .98, "Some count bounds unavailable; see audit", transform=ax.transAxes,
                        ha="right", va="top", fontsize=9, color=COLORS[i])
            for j in range(len(sums)):
                rows.append({"page": stem, "sample": sample["label"], "observable": spec["branch"],
                             "weight_branch": spec["weight"], "bin_low": bins[j], "bin_high": bins[j+1],
                             "entries": int(counts[j]), "sum_weight": sums[j], "sum_weight2": sumw2[j],
                             "stat_mode": "mc" if raw_count else stat_mode,
                             "uncertainty_model": stats["uncertainty_model"],
                             "source_luminosity_invfb": "" if raw_count else source_lumi_invfb,
                             "display_luminosity_invfb": "" if raw_count else (target_lumi if stat_mode == "projected" else source_lumi_invfb),
                             "luminosity_scale": stats["luminosity_scale"], "display_bin_width": display_divisor[j],
                             "yield": stats["yield"][j], "yield_variance": stats["variance"][j],
                             "finite_mc_variance": stats["finite_mc_variance"][j],
                             "projected_data_variance": "" if stats["projected_data_variance"] is None else stats["projected_data_variance"][j],
                             "yield_covariance_sd": stats["covariance_sd"][j],
                             "yield_interval_lower": stats["interval_lower"][j], "yield_interval_upper": stats["interval_upper"][j],
                             "interval_method": stats["interval_method"][j],
                             "interval_coverage": "" if stats["interval_coverage"] is None else stats["interval_coverage"],
                             "interval_reference_count": stats["interval_reference_count"][j],
                             "interval_underlying_mean": stats["interval_underlying_mean"][j],
                             "interval_correction_factor": stats["interval_correction_factor"][j],
                             "display_value": display[j], "display_error": error[j],
                             "display_covariance_sd": error[j],
                             "display_interval_lower": lower[j], "display_interval_upper": upper[j]})
            audit.append({"page": stem, "sample": sample["label"], "observable": spec["branch"],
                          "mask": mask_name, "weight_branch": spec["weight"], "valid_entries": len(values),
                          "sum_weight": float(weight.sum()), "underflow_entries": int(np.sum(values < bins[0])),
                          "overflow_entries": int(np.sum(values > bins[-1])),
                          "underflow_weight": float(weight[values < bins[0]].sum()),
                          "overflow_weight": float(weight[values > bins[-1]].sum()),
                          "window_sum_weight": float(sums.sum()),
                          "window_display_yield": float(stats["yield"].sum()),
                          "window_display_yield_variance": float(stats["variance"].sum()),
                          "luminosity_scale": stats["luminosity_scale"],
                          "uncertainty_model": stats["uncertainty_model"],
                          "interval_methods": sorted(set(stats["interval_method"])),
                          "interval_coverage": stats["interval_coverage"],
                          "interval_unavailable_bins": np.flatnonzero(~np.isfinite(stats["interval_lower"])).tolist(),
                          "interval_warnings": stats["interval_warnings"],
                          "fills_per_event": 1})
        ax.set_xlabel(spec["xlabel"])
        if spec["weight"] == "mc_count":
            ylabel = "Generated MC events / bin"
        elif spec["weight"] == "weight_proj":
            ylabel = r"$dN_{\rm corrected}/d|t|$ [GeV$^{-2}$]" if density else "Corrected visible yield / bin"
        else:
            ylabel = r"$dN_{\rm visible}/d|t|$ [GeV$^{-2}$]" if density else "Expected visible events / bin"
        ax.set_ylabel(ylabel)
        ax.set_xlim(spec["edges"][0], spec["edges"][-1])
        if spec["log"]:
            ax.set_yscale("log")
            if positive:
                ax.set_ylim(max(min(positive)*.5, 1e-6), max(positive)*2)
        else:
            ax.set_ylim(bottom=0)
        ax.grid(alpha=.18)
    for ax in list(axs.flat)[len(specs):]:
        ax.set_visible(False)
    handles = [Line2D([0], [0], color=COLORS[i], ls=STYLES[i], lw=2) for i in range(len(samples))]
    finish(fig, outdir, stem, title, subtitle, handles=handles,
           labels=[s["label"] + rf", $\Upsilon({s.get('state', '1S')})$" for s in samples], legend_ncol=2)


def correlations(samples, args, rows):
    specs = [("true_decay_eMinus_pT", "true_ePlus_pT", edges(0, 10, 50), edges(0, 10, 50),
              r"Generated decay $e^-$ $p_T$ [GeV]", r"Generated decay $e^+$ $p_T$ [GeV]"),
             ("true_ePlus_pT", "true_ePlus_eta", edges(0, 10, 50), edges(-4, 4, 50),
              r"Generated decay $e^{\pm}$ $p_T$ [GeV]", r"Generated decay $e^{\pm}$ $\eta$"),
             ("true_W", "true_t", edges(0, max(s["wmax"] for s in samples), 60), edges(0, 1.5, 50),
              r"Generated $W_{\gamma p}$ [GeV]", r"Generated $|t|$ [GeV$^2$]")]
    fig, axs = plt.subplots(2, 2, figsize=(12.4, 9.8))
    for ax, (xb, yb, xe, ye, xl, yl) in zip(axs.flat, specs):
        for i, sample in enumerate(samples):
            d = sample["data"]
            x, y, w = d[xb], d[yb], np.ones(len(d[xb]), dtype=float)
            mask = sample["truth"] & np.isfinite(x) & np.isfinite(y) & (x > -900) & (y > -900)
            x, y, w = x[mask], y[mask], w[mask]
            if xb == "true_ePlus_pT":
                x2, y2 = d["true_decay_eMinus_pT"], d["true_decay_eMinus_eta"]
                m2 = sample["truth"] & np.isfinite(x2) & np.isfinite(y2) & (x2 > -900) & (y2 > -900)
                x, y, w = np.r_[x, x2[m2]], np.r_[y, y2[m2]], np.r_[w, np.ones(int(m2.sum()), dtype=float)]
            h, _, _ = np.histogram2d(x, y, bins=[xe, ye], weights=w)
            if h.max() > 0:
                ax.contour(.5*(xe[1:]+xe[:-1]), .5*(ye[1:]+ye[:-1]), h.T,
                           levels=h.max()*np.array([.1,.3,.5,.7,.9]), colors=[COLORS[i]],
                           linestyles=STYLES[i], linewidths=1.2)
            for j, k in np.ndindex(h.shape):
                rows.append({"page": "generator_kinematic_correlations", "sample": sample["label"],
                             "x_branch": xb, "y_branch": yb, "x_low": xe[j], "x_high": xe[j+1],
                             "y_low": ye[k], "y_high": ye[k+1], "sum_weight": h[j,k],
                             "stat_mode": "mc_diagnostic", "normalization": "raw_mc_entries", "luminosity_invfb": ""})
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(alpha=.18)
    axs.flat[-1].set_visible(False)
    finish(fig, args.out_dir, "generator_kinematic_correlations", "Generated kinematic correlations",
           rf"$\Upsilon({getattr(args, 'state', '1S')})$ | raw MC | contours at 10, 30, 50, 70 and 90% of each histogram maximum",
           handles=[Line2D([0],[0], color=COLORS[i], ls=STYLES[i]) for i in range(len(samples))],
           labels=[s["label"] for s in samples])


def stored_pages(samples, args, rows, audit):
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.6))
    for ax, name, xlabel, log in zip(axs, ["h_EleID_EpRatio", "h_EleID_IsoFraction"],
                                    [r"Electron ID $E_{\rm cal}/p_{\rm track}$", r"Electron ID $E_{\rm cand}/E_{\rm cone}$"], [True, False]):
        found = False
        for i, sample in enumerate(samples):
            if name not in sample["stored"]:
                audit.append({"sample":sample["label"],"missing_stored_histogram":name});continue
            h = sample["stored"][name]; found=True
            ax.stairs(h["values"], h["edges"], color=COLORS[i], ls=STYLES[i])
            for j, value in enumerate(h["values"]):
                rows.append({"sample":sample["label"], "histogram":name, "bin":j+1,
                             "bin_low":h["edges"][j],"bin_high":h["edges"][j+1], "count":value,
                             "definition":"stored producer entries", "loss_from_previous_percent":""})
        ax.set_xlabel(xlabel);ax.set_ylabel("Stored candidate entries");ax.grid(alpha=.18)
        if log and found:ax.set_yscale("log")
        else:ax.set_ylim(bottom=0)
        if not found:ax.text(.5,.5,"Stored histogram unavailable",ha="center",transform=ax.transAxes)
    handles=[Line2D([0],[0],color=COLORS[i],ls=STYLES[i]) for i in range(len(samples))]
    finish(fig,args.out_dir,"electron_identification","Electron-identification distributions",
           "Stored producer histograms | raw entries, before this task's tuple selections",handles=handles,labels=[s["label"] for s in samples])

    definitions = [
        "All processed events",
        "Generated Upsilon decay found",
        "Both generated leptons have |eta| < 3.5",
        "At least one identified central positron",
        "Identified central positron and electron available",
        "Second central electron or eligible tagger available",
        "Valid lepton assignment and mass in [8, 12] GeV",
        "Producer |sum(E - pz) - 2 Ee| <= 0.2 Ee",
        "Reconstructed photon-proton W >= 0",
        "RP or OMD proton candidate",
    ]
    fig = plt.figure(figsize=(14, 8.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1.2], height_ratios=[3, 1], hspace=.08, wspace=.20)
    ax = fig.add_subplot(gs[0, 0])
    lossax = fig.add_subplot(gs[1, 0], sharex=ax)
    textax = fig.add_subplot(gs[:, 1])
    textax.set_axis_off()
    found = False
    for i, sample in enumerate(samples):
        if "h_CutFlow" not in sample["stored"]:
            audit.append({"sample": sample["label"], "missing_stored_histogram": "h_CutFlow"})
            continue
        histogram = sample["stored"]["h_CutFlow"]
        counts = np.asarray(histogram["values"], dtype=float)
        if len(counts) != 10:
            raise ValueError(f"{sample['path']}: expected ten producer-counter stages, found {len(counts)}")
        found = True
        losses = np.full(10, np.nan)
        losses[1:] = np.divide(100 * (counts[:-1] - counts[1:]), counts[:-1],
                               out=np.full(9, np.nan), where=counts[:-1] > 0)
        # Reconstruction starts from all processed events, independently of truth stages 2 and 3.
        losses[3] = np.nan
        ax.stairs(counts[:3], np.arange(.5, 4), baseline=None, color=COLORS[i], ls=STYLES[i])
        ax.stairs(counts[3:], np.arange(3.5, 11), baseline=None, color=COLORS[i], ls=STYLES[i])
        lossax.stairs(losses, np.arange(.5, 11), baseline=None, color=COLORS[i], ls=STYLES[i])
        for j, count in enumerate(counts):
            rows.append({"sample": sample["label"], "histogram": "h_CutFlow", "bin": j + 1,
                         "bin_low": j + .5, "bin_high": j + 1.5, "count": count,
                         "definition": definitions[j],
                         "loss_from_previous_percent": float(losses[j]) if np.isfinite(losses[j]) else ""})
        audit.append({"sample": sample["label"], "histogram": "h_CutFlow",
                      "stored_labels": histogram["labels"],
                      "nested_stage_sequences": [[1, 2, 3], [1, 4, 5, 6, 7, 8, 9, 10]],
                      "undefined_adjacent_loss_stages": (np.flatnonzero(~np.isfinite(losses)) + 1).tolist()})
    ax.set_ylabel("Simulated events")
    ax.grid(alpha=.18)
    lossax.grid(alpha=.18)
    if found:
        ax.set_yscale("log")
    ax.axvline(3.5, color="0.55", ls=":", lw=1)
    ax.tick_params(labelbottom=False)
    lossax.set_xticks(np.arange(1, 11))
    lossax.set_xlabel("Producer stage")
    lossax.set_ylabel("Step loss [%]")
    lossax.set_ylim(-2, 102)
    lossax.axvspan(3.5, 4.5, color="0.95", zorder=0)
    lossax.text(4, 65, "Separate\nchain", ha="center", va="center", fontsize=10)
    ax.set_xlim(.5, 10.5)
    textax.text(0, 1, "Producer counter definitions", fontweight="bold", fontsize=14,
                va="top", transform=textax.transAxes)
    for j, definition in enumerate(definitions):
        textax.text(0, .92 - j * .086, f"{j + 1}. {definition}", fontsize=10,
                    va="top", transform=textax.transAxes, wrap=True)
    finish(fig, args.out_dir, "event_selection_cutflow_numbered", r"$\Upsilon\to e^+e^-$ producer event counters",
           "Raw MC counts | reconstruction starts from all processed events\nStep losses compare adjacent counters within the truth and reconstruction chains",
           handles=handles, labels=[sample["label"] for sample in samples])


def tracking_page(samples,args,rows):
    fig,axs=plt.subplots(1,len(samples),figsize=(8 if len(samples)==1 else 14,7),squeeze=False)
    for ax,sample in zip(axs.flat,samples):
        d=sample["data"];x=d["reco_nTracks_Central"];y=d["reco_nTracks_Tagger"];w=d["weight"]
        mask=sample["selected"]&np.isfinite(x)&np.isfinite(y)&(x>=0)&(y>=0)
        h,xe,ye=np.histogram2d(x[mask],y[mask],bins=[np.arange(-.5,15.5),np.arange(-.5,8.5)],weights=w[mask])
        positive=h[h>0]
        if len(positive):
            mesh=ax.pcolormesh(xe,ye,np.ma.masked_less_equal(h.T,0),norm=LogNorm(vmin=positive.min(),vmax=max(positive.max(),positive.min()*1.1)),cmap="viridis")
            fig.colorbar(mesh,ax=ax,label="Weighted MC yield at source luminosity")
        for j,k in np.ndindex(h.shape):rows.append({"page":"central_tagger_correlation","sample":sample["label"],"x_branch":"reco_nTracks_Central","y_branch":"reco_nTracks_Tagger","x_low":xe[j],"x_high":xe[j+1],"y_low":ye[k],"y_high":ye[k+1],"sum_weight":h[j,k],"stat_mode":"mc_diagnostic","normalization":"source_weight","luminosity_invfb":args.lumi_invfb})
        ax.set_title(sample["label"]);ax.set_xlabel("Central reconstructed tracks");ax.set_ylabel("Tagger reconstructed tracks")
        ax.set_xticks(np.arange(0,15,2));ax.set_yticks(np.arange(0,8));ax.grid(alpha=.15)
    finish(fig,args.out_dir,"central_tagger_correlation","Central and low-Q² tagger track multiplicity",
           f"MC diagnostic at source luminosity {args.lumi_invfb:g} fb⁻¹\nSelected {getattr(args, 'state', '1S')} events | mass 8–12 GeV, beam-specific W range and exclusivity selection")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state", choices=list(MASSES), default="1S")
    ap.add_argument("--input",required=True);ap.add_argument("--tree",default="UpsilonTree")
    ap.add_argument("--energy-label",default="9x130 GeV");ap.add_argument("--electron-energy",type=float,default=9)
    ap.add_argument("--w-min",type=float,default=20);ap.add_argument("--w-max",type=float,default=70)
    ap.add_argument("--q2-cut",type=float,default=1);ap.add_argument("--lumi-invfb",type=float,default=10)
    ap.add_argument("--stat-mode", choices=["mc", "projected"], default="mc",
                    help="Weighted 1D errors: finite-MC sumw2 or projected data counting statistics")
    ap.add_argument("--projection-lumi-invfb", type=float,
                    help="Target luminosity for projected weighted 1D pages; defaults to source --lumi-invfb")
    add_interval_arguments(ap)
    ap.add_argument("--second-input");ap.add_argument("--second-energy-label",default="9x275 GeV")
    ap.add_argument("--second-electron-energy",type=float,default=9);ap.add_argument("--second-w-min",type=float,default=20)
    ap.add_argument("--second-w-max",type=float,default=100);ap.add_argument("--out-dir",required=True)
    args=ap.parse_args();args.out_dir=Path(args.out_dir);args.out_dir.mkdir(parents=True,exist_ok=True)
    if args.q2_cut<0:ap.error("--q2-cut must be nonnegative")
    if not np.isfinite(args.lumi_invfb) or args.lumi_invfb <= 0:ap.error("--lumi-invfb must be finite and positive")
    if args.projection_lumi_invfb is None:args.projection_lumi_invfb=args.lumi_invfb
    if not np.isfinite(args.projection_lumi_invfb) or args.projection_lumi_invfb <= 0:ap.error("--projection-lumi-invfb must be finite and positive")
    if not np.isfinite(args.interval_coverage) or not 0 < args.interval_coverage < 1:ap.error("--interval-coverage must be finite and between zero and one")
    if args.stat_mode == "mc" and args.projection_lumi_invfb != args.lumi_invfb:
        ap.error("--projection-lumi-invfb differing from source luminosity requires --stat-mode projected")
    samples=[load_sample(args.input,args.energy_label,args.electron_energy,args.w_min,args.w_max,args)]
    if args.second_input:samples.append(load_sample(args.second_input,args.second_energy_label,args.second_electron_energy,args.second_w_min,args.second_w_max,args))
    set_style();rows=[];rows2d=[];storedrows=[];audit=[]
    maxw=max(s["wmax"] for s in samples);maxee=max(s["ee"] for s in samples)
    display_lumi=args.projection_lumi_invfb if args.stat_mode == "projected" else args.lumi_invfb
    error_label=(f"{100*args.interval_coverage:.1f}% median count bounds" if args.count_interval == "median_garwood"
                 else "projected covariance bands") if args.stat_mode == "projected" else "finite-MC errors"
    subtitle=rf"$\mathcal{{L}}_{{ep}}={display_lumi:g}$ fb$^{{-1}}$ | $0\leq Q^2_{{\rm true}}\leq {args.q2_cut:g}$ GeV$^2$ | {error_label}"
    stats_options={"stat_mode":args.stat_mode,"source_lumi_invfb":args.lumi_invfb,
                   "projection_lumi_invfb":args.projection_lumi_invfb,
                   "count_interval":args.count_interval,"interval_coverage":args.interval_coverage}
    generator=[
      ("generator_upsilon_kinematics","Generator-level Upsilon kinematics",[
        panel("true_mass",r"Generated $m_{ee}$ [GeV]",8,12,80,True),panel("true_pT",r"Generated $\Upsilon$ $p_T$ [GeV]",0,2),
        panel("true_y",r"Generated $\Upsilon$ rapidity $y$",-4,4),panel("true_W",r"Generated $W_{\gamma p}$ [GeV]",0,maxw)]),
      ("generator_production_scattered_electron","Generated production and scattered electron",[
        panel("true_Q2",r"Generated $Q^2$ [GeV$^2$]",0,args.q2_cut,60,True),panel("true_t",r"Generated $|t|$ [GeV$^2$]",0,1.5,50,True),
        panel("true_scat_eMinus_E",r"Generated scattered $e^-$ energy [GeV]",0,maxee*1.15),panel("true_scat_eMinus_eta",r"Generated scattered $e^-$ $\eta$",-12,4,64,True)]),
      ("generator_decay_leptons","Generated decay-lepton kinematics",[
        panel("true_ePlus_pT",r"Generated decay $e^+$ $p_T$ [GeV]",0,10),panel("true_ePlus_eta",r"Generated decay $e^+$ $\eta$",-4,4),
        panel("true_decay_eMinus_pT",r"Generated decay $e^-$ $p_T$ [GeV]",0,10),panel("true_decay_eMinus_eta",r"Generated decay $e^-$ $\eta$",-4,4)])]
    projected=[
      (f"projected_{args.state.lower()}_upsilon_kinematics","Projected reconstructed Upsilon kinematics",[
        panel("reco_mass",r"Reconstructed $m_{ee}$ [GeV]",8,12,80,True),panel("reco_pT",r"Reconstructed $\Upsilon$ $p_T$ [GeV]",0,2),
        panel("reco_y",r"Reconstructed $\Upsilon$ rapidity $y$",-4,4),panel("reco_W",r"Reconstructed $W_{\gamma p}$ [GeV]",0,maxw)]),
      (f"projected_{args.state.lower()}_production_exclusivity","Projected production and exclusivity",[
        panel("reco_Q2",r"Reconstructed $Q^2$ [GeV$^2$]",0,args.q2_cut,60,True),panel("reco_t_proj",r"Projective $|t|$ [GeV$^2$]",0,1.5,50,True,"weight_proj"),
        panel("sigma_epz",r"Exclusivity $\Sigma(E-p_z)$ [GeV]",0,maxee*2.6,65)]),
      (f"projected_{args.state.lower()}_decay_leptons","Projected decay-lepton kinematics",[
        panel("reco_ePlus_pT",r"Reconstructed decay $e^+$ $p_T$ [GeV]",0,10),panel("reco_ePlus_eta",r"Reconstructed decay $e^+$ $\eta$",-4,4),
        panel("reco_decay_eMinus_pT",r"Reconstructed decay $e^-$ $p_T$ [GeV]",0,10),panel("reco_decay_eMinus_eta",r"Reconstructed decay $e^-$ $\eta$",-4,4)]),
      (f"projected_{args.state.lower()}_scattered_electron_tracking","Projected scattered electron and tracking",[
        panel("reco_scat_eMinus_E",r"Reconstructed scattered $e^-$ energy [GeV]",0,maxee*1.15),panel("reco_scat_eMinus_eta",r"Reconstructed scattered $e^-$ $\eta$",-12,4,64,True),
        panel("reco_nTracks_Central","Central reconstructed tracks",-.5,14.5,15,True),panel("reco_nTracks_Tagger","Tagger reconstructed tracks",-.5,7.5,8,True)])]
    gen_subtitle = rf"$\Upsilon({args.state})\to e^+e^-$ | generated MC entries | $0\leq Q^2_{{\rm true}}\leq {args.q2_cut:g}$ GeV$^2$"
    for stem,title,specs in generator:
        for spec in specs: spec["weight"] = "mc_count"
        draw_pages(samples,specs,stem,title,gen_subtitle,"truth",args.out_dir,rows,audit)
    correlations(samples,args,rows2d)
    projected_mask = "selected" if args.second_input else "projected"
    projected_subtitle = ("\nMass 8–12 GeV; beam-specific reconstructed W range and exclusivity selection"
                          if args.second_input else
                          "\nBroad reco selection: mass 8–12 GeV and W ≥ 0; exclusivity not applied")
    for stem,title,specs in projected:
        draw_pages(samples,specs,stem,title,subtitle+projected_subtitle,projected_mask,args.out_dir,rows,audit,**stats_options)
    methods=[("t_estimators_electron_hybrid","Reconstructed electron and hybrid t estimators",[
              panel("reco_t_eX",r"$|t_{eX}|$ [GeV$^2$]",0,5,50,True),panel("reco_t_BABE",r"$|t_{\rm BABE}|$ [GeV$^2$]",0,5,50,True),panel("reco_t_eXBABE",r"$|t_{eX,BABE}|$ [GeV$^2$]",0,5,50,True)]),
             ("t_estimators_projective_transverse","Projective and transverse momentum transfer",[
              panel("reco_t_proj",r"Projective $|t|$ [GeV$^2$]",0,5,50,True,"weight_proj"),panel("reco_pT2",r"$p^2_{T,\Upsilon}$ [GeV$^2$]",0,5,50,True)])]
    for stem,title,specs in methods:draw_pages(samples,specs,stem,title,subtitle+"\nBeam-specific reconstructed W range and exclusivity selection","selected",args.out_dir,rows,audit,True,**stats_options)
    stored_pages(samples,args,storedrows,audit);tracking_page(samples,args,rows2d)
    write_csv(args.out_dir/"kinematic_bins.csv",rows);write_csv(args.out_dir/"kinematic_2d_bins.csv",rows2d);write_csv(args.out_dir/"stored_histogram_bins.csv",storedrows)
    selection=[]
    for sample in samples:
        d=sample["data"]
        selection.append({"input":sample["path"],"label":sample["label"],"entries":len(d["weight"]),"input_sum_weight":float(d["weight"].sum()),
                          "input_sum_weight_proj":float(d["weight_proj"].sum()),"electron_energy":sample["ee"],"w_min":sample["wmin"],"w_max":sample["wmax"],
                          "masks":{k:{"entries":int(sample[k].sum()),"sum_weight":float(d["weight"][sample[k]].sum())} for k in ["truth","projected","selected"]}})
    info={"proton_selection": proton_selection_metadata(), "normalization":"Generator pages use unit weights (raw MC entries). Weighted 1D pages use prepared input weights, scaled only by target/source luminosity in projected mode. No branching-ratio correction or post-selection normalization.",
          "stat_mode":args.stat_mode,"source_luminosity_invfb":args.lumi_invfb,
          "projection_luminosity_invfb":args.projection_lumi_invfb,
          "weighted_1d_luminosity_scale":display_lumi/args.lumi_invfb,
          "luminosity_label_invfb":display_lumi,"luminosity_note":"Input weights must already carry the source luminosity normalization; 2D diagnostics retain source normalization.",
          "state":args.state,"q2_cut":args.q2_cut,"truth_definition":f"|true_mass - {MASSES[args.state]:.5f}| <= 0.2, finite nonnegative weight and 0 <= true_Q2 <= q2_cut",
          "projected_definition":"truth plus 8 <= reco_mass <= 12 and reco_W >= 0, with the declared proton selection",
          "projected_pages_mask":projected_mask,
          "projected_pages_selection_reason":("Two-energy note comparison adds beam-specific reconstructed W and exclusivity selections."
                                                if args.second_input else "Single-sample projected-state note pages use the broad reconstructed selection."),
          "selected_definition":"projected plus |sigma_epz - 2 Ee| <= 0.2 Ee and beam-specific reconstructed W range",
          "projective_definition":"Stored absolute weight_proj used for positive finite weight_proj, with target/source luminosity scaling only in projected mode. Projected data variance is r*sum(weight_proj**2/weight); no additional eXBABE cut is imposed by notebook03 t-page selection.",
          "stored_histograms":"Raw producer counts. The truth sequence is 1-2-3 and reconstruction independently follows 1-4-5-6-7-8-9-10. Adjacent loss is undefined at stages 1 and 4 and for zero preceding counts. Stage 6 requires an additional central electron or eligible tagger candidate. Stage 9 requires reconstructed W >= 0. State/Q2/W/proton selections and weight normalization are applied separately to tuple masks.",
          "errors":("Projected weighted 1D errors: sqrt(r*sum(weight)) for ordinary yields; sqrt(r*sum(weight_proj**2/weight)) for corrected yields. Finite-MC variances are retained separately, not added."
                    if args.stat_mode == "projected" else "sqrt(sum of squared event weights), finite-Monte-Carlo errors; not projected data counting errors"),
          "count_interval":args.count_interval,"interval_coverage":args.interval_coverage,
          "interval_display":"Projected median_garwood endpoints are ensemble-median expected confidence bounds, evaluated at the median integer Poisson count. Expected-mean points are unchanged, not rounded. For constant within-bin projective corrections, underlying count limits are multiplied by that correction. Variable or unknown empty-bin corrections have unavailable bounds and audited warnings; no correction is guessed. Standard_deviation displays normal-quantile-scaled covariance SD, not count confidence limits. MC bars retain their original SD convention.",
          "interval_columns":"CSV yield/display_interval_lower and upper are absolute displayed endpoints. interval_reference_count is the underlying integer Poisson median. display_error is retained as a compatibility alias for display_covariance_sd, not an interval half-width.",
          "error_scope":"Signal-only counting statistics with fixed normalization and correction factors. Empty expected bins have zero model variance but ordinary median-Garwood forecasts have a nonzero upper count bound. Confidence-bound display never changes the covariance. No background, correction-factor, or other systematic uncertainty is included.",
          "bin_width_rule":"Display density, standard deviation and interval endpoints divide by bin width; density variance divides by bin width squared. CSV sum_weight/sum_weight2 retain the source weights; yield/variance use the selected display luminosity.",
          "covariance_scope":"Each weighted 1D panel has one fill per event and disjoint bins. No joint covariance between different observables or samples is supplied.",
          "diagnostic_pages":"Generator pages/correlations retain raw MC entries. Stored electron-ID/cutflow retain producer counts. The track correlation retains source weighted MC yields. No projected uncertainties or target luminosity scaling are applied to these pages, including the combined two-lepton correlation.",
          "correlations":"Contours at fixed fractions of each histogram peak, not probability-content contours; combined decay +/- panel counts both leptons",
          "samples":selection,"histogram_ranges":audit}
    (args.out_dir/"kinematic_selection_audit.json").write_text(json.dumps(info,indent=2))
    print(f"Wrote reviewed-note kinematic pages to {args.out_dir}")

if __name__=="__main__":main()
