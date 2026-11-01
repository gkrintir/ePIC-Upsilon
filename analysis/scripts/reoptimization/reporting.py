"""Numerical tables and publication-style diagnostics for reoptimization.

The plots deliberately distinguish detector-response diagnostics, conditional
profile limits and sampling bands.  A successful renderer never upgrades a
failed or unavailable inference to a validated measurement.
"""

from pathlib import Path
import csv
import json
import math
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {"1S": "#0072B2", "2S": "#D55E00", "3S": "#009E73"}
OBSERVABLES = {
    "pt2": (r"$p_{T,\Upsilon}^2$", r"$b_{p_T^2}$", "pT2"),
    "abs_t_via_pt2": (r"$|t|$", r"$b_{|t|}$", "absolute t via reconstructed pT2"),
    "w": (r"$W$", "", "photon-proton center-of-mass energy"),
    "y": (r"$y_\Upsilon$", "", "Upsilon rapidity"),
}
SCOPE = ("Truth-fiducial signal-only projection with fixed response and fixed outside-variable-range "
         "feed-in. No mass-component extraction, background or systematic "
         "uncertainty is included. Response-MC variations are reported separately.")


def _clean(value):
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (np.ndarray, tuple, list)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def _write_json(path, data):
    Path(path).write_text(json.dumps(_clean(data), indent=2, allow_nan=False) + "\n")


def _write_csv(path, rows, fields=None):
    common = ["state", "observable", "status", "overall_status", "stage_status",
              "failed_checks", "publication_ready"]
    fields = list(dict.fromkeys(list(fields or []) + common + [k for row in rows for k in row]))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            clean = _clean(row)
            writer.writerow({k: (json.dumps(v, separators=(",", ":"))
                                    if isinstance(v, (dict, list)) else v)
                             for k, v in clean.items()})


def _number(value, digits=3):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "not available"
    return f"{value:.{digits}g}" if np.isfinite(value) else "not available"


def _tex(value):
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
                    "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
                    "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(char, char) for char in str(value))


def _style():
    return {
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
        "font.size": 12, "axes.labelsize": 15, "axes.titlesize": 13,
        "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11.5,
        "axes.linewidth": 1.0, "axes.grid": False, "legend.frameon": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "lines.linewidth": 1.7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 200,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "text.usetex": False,
    }


def _label(state, observable, config, subtitle):
    energy = f"{float(config.get('electron_energy', 9)):g}"
    energy += r"\times" + f"{float(config.get('proton_energy', 275)):g}"
    lumi = float(config.get("target_lumi_invfb", 2.5))
    obs = OBSERVABLES.get(observable, (observable, "b", observable))[0]
    prefix = "SYNTHETIC SOFTWARE TEST\n" if config.get("synthetic", False) else ""
    return (prefix + rf"$\Upsilon({state})$, ${energy}$ GeV, "
            rf"$\mathcal{{L}}={lumi:g}$ fb$^{{-1}}$" + "\n" +
            f"{subtitle}  ·  {obs}")


def _overall_status(entry):
    inference=entry.get("inference") or {}
    return entry.get("overall_status", entry.get("status",
        inference.get("overall_status", "NOT_VALIDATED")))


def _inference_notice(entry):
    # Always display the actual final case outcome, including on training-only
    # diagnostics. A passing component must not hide a failed overall result.
    return "Overall status: " + str(_overall_status(entry))


def _save(fig, base, captions, caption):
    paths = []
    for suffix in ("pdf", "png"):
        path = base.with_suffix("." + suffix)
        fig.savefig(path, bbox_inches="tight", pad_inches=0.12)
        paths.append(str(path))
    plt.close(fig)
    captions[base.name] = caption
    return paths


def _array(mapping, key, shape=None):
    if key not in mapping or mapping[key] is None:
        return None
    value = np.asarray(mapping[key], dtype=float)
    if shape is not None and value.shape != shape:
        raise ValueError(f"{key} shape {value.shape} does not match {shape}")
    return value


def _coarse(entry):
    model = entry.get("model") or {}
    return model.get("coarse") or entry.get("coarse") or {}


def _safe_stem(state, observable):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{state}_{observable}")


def _bins_text(edges):
    return ", ".join(f"{float(x):g}" for x in edges)


def _unit(observable, tex=False):
    return {"w": "GeV", "y": ""}.get(observable, r"GeV$^2$" if tex else "GeV²")


def _axis(observable, reconstructed=False):
    key = "pt2" if reconstructed and observable == "abs_t_via_pt2" else observable
    label = OBSERVABLES.get(key, (key,))[0]
    unit = _unit(observable, True)
    return label + (f" [{unit}]" if unit else "")


def _cross_section(entry, observable, edges):
    """Resolve an already-complete count-to-quantity conversion exactly once."""
    n = len(edges) - 1
    supplied = entry.get("cross_section_scaling")
    if supplied is not None:
        scale = _array(supplied, "factors", (n,))
        if scale is None:
            raise ValueError("cross_section_scaling must contain factors")
        labels = {"w": r"$\langle\sigma_{\gamma p}\rangle_{\Phi}$ [pb]",
                  "y": r"$d\sigma_{ep}/dy_\Upsilon$ [pb]",
                  "pt2": r"$d\sigma_{\gamma p}/dp_T^2$ [pb/GeV$^2$]",
                  "abs_t_via_pt2": r"$d\sigma_{\gamma p}/d|t|$ [pb/GeV$^2$]"}
        meta = dict(supplied)
        meta.setdefault("ylabel", labels[observable])
    else:
        factor = entry.get("cross_section_factor")
        if observable in ("w", "y") and factor is not None:
            raise ValueError("W and rapidity require explicit observable cross_section_scaling")
        scale = float(factor) / np.diff(edges) if factor is not None else np.ones(n)
        tex = "p_T^2" if observable == "pt2" else "|t|"
        meta = {"quantity": "gamma_p_differential" if factor is not None else "truth_yield",
                "units": "pb/GeV^2" if factor is not None else "events",
                "ylabel": rf"$d\sigma_{{\gamma p}}/d{tex}$ [pb/GeV$^2$]" if factor is not None else "Expected truth-bin events"}
    if not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ValueError("A finite positive cross-section conversion is required")
    return scale, meta


def _fit_reco_indices(coarse):
    n = len(coarse["truth_edges"]) - 1
    indices = np.asarray((coarse.get("metadata") or {}).get("fit_reco_indices", np.arange(n)), int)
    if indices.shape != (n,) or len(set(indices)) != n:
        raise ValueError("Invalid fit reconstructed-bin index mapping")
    return indices


def _failed_checks(entry):
    recorded = entry.get("failed_checks")
    if recorded is not None:
        return list(dict.fromkeys([str(x) for x in recorded] +
                                  [str(x) for x in entry.get("errors", [])]))
    failed = []
    for prefix, gates in (("inference", (entry.get("inference") or {}).get("gates", {})),
                          ("additional", entry.get("additional_validation_gates", {}))):
        if isinstance(gates, dict):
            failed.extend(prefix + "." + key for key, value in gates.items() if value is False)
    failed.extend(str(x) for x in entry.get("errors", []))
    return failed


def _plot_response(base, state, observable, entry, config, captions):
    coarse = _coarse(entry)
    if not coarse:
        return []
    matrix = _array(coarse, "response")
    truth_edges = _array(coarse, "truth_edges")
    reco_edges = _array(coarse, "reco_edges")
    if matrix is None or truth_edges is None or reco_edges is None:
        return []
    n = len(truth_edges) - 1
    if matrix.shape != (len(reco_edges) - 1, n):
        raise ValueError("Response and bin edges have incompatible shapes")
    fig, ax = plt.subplots(figsize=(7.4, 6.3), layout="constrained")
    img = ax.imshow(matrix, origin="lower", aspect="auto", cmap="Blues", vmin=0,
                    vmax=max(float(np.nanmax(matrix)), .01))
    ax.set_xticks(np.arange(n), [f"{a:g}–{b:g}" for a, b in zip(truth_edges[:-1], truth_edges[1:])])
    ax.set_yticks(np.arange(len(reco_edges)-1),
                  [f"{a:g}–{b:g}" for a, b in zip(reco_edges[:-1], reco_edges[1:])])
    for r in range(matrix.shape[0]):
        for t in range(matrix.shape[1]):
            v = matrix[r, t]
            if np.isfinite(v):
                ax.text(t, r, f"{v:.3f}", ha="center", va="center", fontsize=11,
                        color="white" if v > .6*img.norm.vmax else "#16283e")
    ax.set_xlabel("Truth " + _axis(observable))
    ax.set_ylabel("Reconstructed " + _axis(observable, True))
    ax.set_title(_label(state, observable, config, "Training response\n"+_inference_notice(entry)), pad=12)
    fig.colorbar(img, ax=ax, label="Selected / generated truth-bin events", shrink=.85)
    return _save(fig, base.with_name(base.name + "_response"), captions,
                 _inference_notice(entry)+". "+
                 "Training-only response. Any configured control intervals are identified by the recorded bin-index mapping. "
                 "Each truth column is divided by all generated truth-bin events, including "
                 "missed reconstruction. Column sums are modeled efficiencies, not unity. " + SCOPE)


def _plot_metrics(base, state, observable, entry, config, captions):
    coarse = _coarse(entry)
    if not coarse:
        return []
    edges = _array(coarse, "truth_edges")
    if edges is None:
        return []
    n = len(edges)-1
    metrics = coarse.get("metrics") or {}
    reco = _array(coarse, "reco_expected")
    feed = _array(coarse, "feed_in")
    signal = _array(metrics, "expected_reco_signal")
    if reco is None or feed is None or signal is None:
        return []
    if len(reco) < n or len(feed) < n or len(signal) != n:
        raise ValueError("Count/metric shapes disagree with selected edges")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.1), layout="constrained")
    fit_indices = _fit_reco_indices(coarse)
    x = np.arange(n)
    color = COLORS.get(state, "#0072B2")
    axes[0].bar(x, signal, .65, color=color, alpha=.85, label="Selected signal")
    axes[0].bar(x, feed[fit_indices], .65, bottom=signal, color="#c7c7c7", hatch="//",
                edgecolor="#777777", linewidth=.7, label="Outside-range feed-in")
    threshold = (entry.get("optimization") or {}).get("thresholds", {}).get("min_expected_reco")
    if threshold is not None:
        axes[0].axhline(float(threshold), color="#555555", linestyle="--", linewidth=1.2,
                       label=f"Design threshold {float(threshold):g}")
    axes[0].set_ylabel("Expected reconstructed events")
    axes[0].set_ylim(0, max(float(np.max(reco[fit_indices])), float(threshold or 0), 1)*1.35)
    axes[0].legend(loc="upper right", fontsize=10.5)
    styles = (("purity", "Purity", "o", "#0072B2"),
              ("selected_stability", "Selected stability", "s", "#D55E00"),
              ("diagonal_acceptance", "Diagonal acceptance", "^", "#009E73"),
              ("efficiency_modeled", "Modeled efficiency", "D", "#777777"))
    for key, title, marker, c in styles:
        values = _array(metrics, key)
        if values is not None:
            if len(values) != n:
                raise ValueError(f"Metric {key} does not match binning")
            axes[1].plot(x, values*100, marker=marker, color=c, label=title, markersize=5)
    axes[1].set_ylabel("Fraction [%]")
    axes[1].set_ylim(0, 115)
    axes[1].legend(loc="upper right", fontsize=10.5)
    for panel, ax in zip(("a", "b"), axes):
        ax.set_xticks(x, [f"{a:g}–{b:g}" for a,b in zip(edges[:-1],edges[1:])])
        ax.set_xlabel("Reconstructed " + _axis(observable, True) if panel == "a"
                      else "Shared truth/reco edges" + (f" [{_unit(observable, True)}]" if _unit(observable) else ""))
        ax.text(.015, .97, f"({panel})", transform=ax.transAxes, va="top", fontsize=13)
    fig.suptitle(_label(state, observable, config, "Training binning diagnostics\n"+_inference_notice(entry)), fontsize=13)
    return _save(fig, base.with_name(base.name + "_binning"), captions,
                 _inference_notice(entry)+". "+
                 "Expected training-split yields rescaled to the full requested projection. "
                 "Purity includes feed-in. Selected stability divides by selected events "
                 "from the truth bin; diagonal acceptance divides by all generated events. "
                 "The count panel uses reconstructed-variable fit bins; the metrics use the shared numeric truth/reco boundaries. "
                 "The count threshold is a design criterion, not a confidence-coverage theorem.")


def _draw_absolute_band(ax, x, central, lower, upper, color, label, xerr):
    # A sampling-quantile interval need not contain the Asimov central value.
    ax.errorbar(x, central, xerr=xerr, fmt="o", color=color, ms=4.5, label=label)
    valid = np.isfinite(lower) & np.isfinite(upper)
    ax.vlines(x[valid], lower[valid], upper[valid], color=color, linewidth=1.5)
    ax.plot(x[valid], lower[valid], "_", color=color, ms=6)
    ax.plot(x[valid], upper[valid], "_", color=color, ms=6)


def _plot_coarse(base, state, observable, entry, config, captions):
    inference = entry.get("inference") or {}
    coarse = inference.get("coarse") or {}
    if not coarse or str(coarse.get("status", "")).startswith(("FAILED", "BLOCKED", "UNSUPPORTED")):
        return []
    edges = _array(coarse, "truth_edges")
    fitted = _array(coarse, "fitted_truth")
    reference = _array(coarse, "reference_truth")
    if edges is None or fitted is None or reference is None:
        return []
    n = len(edges)-1
    if fitted.shape != (n,) or reference.shape != (n,):
        raise ValueError("Coarse extraction shape disagrees with truth edges")
    lower = _array(coarse, "toy_lower", (n,))
    upper = _array(coarse, "toy_upper", (n,))
    width = np.diff(edges)
    scale, conversion = _cross_section(entry, observable, edges)
    centers = (edges[:-1]+edges[1:])/2
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 7.2), sharex=True,
                             gridspec_kw={"height_ratios": [2.2, 1]}, layout="constrained")
    color = COLORS.get(state, "#0072B2")
    axes[0].stairs(reference*scale, edges, color="#222222", label="Held-out MC truth")
    if lower is not None and upper is not None:
        _draw_absolute_band(axes[0], centers, fitted*scale, lower*scale, upper*scale,
                            color, "Poisson extraction; toy sampling band", width/2)
    else:
        axes[0].errorbar(centers, fitted*scale, xerr=width/2, fmt="o", color=color,
                         label="Poisson extraction; no interval available")
    axes[0].set_ylabel(conversion["ylabel"])
    heights = np.r_[reference*scale, fitted*scale,
                    upper*scale if upper is not None else []]
    finite_heights = heights[np.isfinite(heights)]
    # Reserve space for the two-line legend above the sampling bands.
    axes[0].set_ylim(0, 1.30*max(float(finite_heights.max()), 1e-12)
                     if len(finite_heights) else 1.)
    axes[0].legend(loc="upper right", fontsize=10)
    axes[0].set_title(_label(state, observable, config,
                            "Held-out extraction\n"+_inference_notice(entry)), pad=12)
    axes[0].text(.025, .97, "(a)", transform=axes[0].transAxes, va="top", fontsize=13)
    valid = reference > 0
    ratio = np.divide(fitted, reference, out=np.full(n, np.nan), where=valid)
    if lower is not None and upper is not None:
        _draw_absolute_band(axes[1], centers, ratio,
                            np.divide(lower, reference, out=np.full(n,np.nan), where=valid),
                            np.divide(upper, reference, out=np.full(n,np.nan), where=valid),
                            color, None, width/2)
    else:
        axes[1].errorbar(centers, ratio, xerr=width/2, fmt="o", color=color)
    axes[1].axhline(1, color="#555555", linestyle="--", linewidth=1)
    axes[1].set_ylabel("Extraction / truth")
    axes[1].set_ylim(bottom=0)
    axes[1].set_xlabel(_axis(observable))
    axes[1].text(.025,.97,"(b)",transform=axes[1].transAxes,va="top",fontsize=13)
    return _save(fig, base.with_name(base.name + "_differential"), captions,
                 _inference_notice(entry)+". The central result fits independent held-out reconstructed expected means "
                 "with the frozen training response. Vertical bars are central toy sampling "
                 "quantiles, not confidence limits. Truth in the ratio is fixed. " + SCOPE)


def _plot_slope(base, state, observable, entry, config, captions):
    inference = entry.get("inference") or {}
    slope = inference.get("slope") or {}
    profile = slope.get("profile") or []
    if not profile:
        return []
    if isinstance(profile, dict):
        x = np.asarray(profile["b"], float)
        delta = np.asarray(profile["delta_deviance"], float)
    else:
        x = np.asarray([row["b"] for row in profile],float)
        delta = np.asarray([row["delta_deviance"] for row in profile],float)
    fig, ax = plt.subplots(figsize=(7.4,5.5),layout="constrained")
    ax.plot(x, delta, color=COLORS.get(state,"#0072B2"), label="Held-out expected means")
    threshold = float(slope.get("delta_deviance_threshold",slope.get("delta_deviance",1.0)))
    ax.axhline(threshold, color="#555555", linestyle="--", linewidth=1.1,
               label=rf"$\Delta D={threshold:g}$, conditional profile")
    b = slope.get("b")
    if b is not None and np.isfinite(float(b)):
        ax.axvline(float(b), color=COLORS.get(state,"#0072B2"), linestyle=":", linewidth=1)
    finite = delta[np.isfinite(delta)]
    ax.set_ylim(0,min(max(float(np.max(finite)) if len(finite) else 4,2),8))
    ax.set_xlabel(OBSERVABLES.get(observable,("", "b"))[1]+r" [GeV$^{-2}$]")
    ax.set_ylabel(r"Profile $\Delta D$")
    ax.set_title(_label(state,observable,config,"Forward-folded slope diagnostic\n"+_inference_notice(entry)),pad=12)
    ax.legend(loc="upper right",fontsize=10.5)
    return _save(fig,base.with_name(base.name+"_slope_profile"),captions,
                 _inference_notice(entry)+". Finite-bin exponential prediction forward-folded through the frozen response. "
                 "Amplitude is profiled in reconstructed Poisson likelihood. The displayed "
                 "profile threshold is asymptotic and conditional. Its adequacy must be assessed "
                 "with the separately recorded independent toy coverage; this is not an exact "
                 "low-count confidence construction. " + SCOPE)


def _plot_coverage(base,state,observable,entry,config,captions):
    rows=(entry.get("inference") or {}).get("slope_toys") or []
    valid=[row for row in rows if row.get("coverage") is not None and row.get("input_b") is not None]
    if not valid:
        return []
    x=np.asarray([row["input_b"] for row in valid],float)
    coverage=np.asarray([row["coverage"] for row in valid],float)
    err=np.asarray([row.get("coverage_mc_se",np.nan) for row in valid],float)
    fig,axes=plt.subplots(1,2,figsize=(12,5.1),layout="constrained")
    color=COLORS.get(state,"#0072B2")
    axes[0].errorbar(x,100*coverage,yerr=100*err,fmt="o",color=color,capsize=3)
    nominal=float(config.get("inference",{}).get("coverage",0.6826894921370859))
    axes[0].axhline(100*nominal,color="#555555",linestyle="--",label=f"Nominal {100*nominal:.1f}%")
    axes[0].set_ylabel("Empirical coverage [%]")
    axes[0].set_ylim(0,105)
    axes[0].legend(loc="lower right",fontsize=10)
    bias=np.asarray([row.get("bias_over_sd",np.nan) for row in valid],float)
    axes[1].plot(x,bias,"o",color=color)
    axes[1].axhline(0,color="#555555",linestyle="--",linewidth=1)
    axes[1].set_ylabel(r"Slope bias / toy standard deviation")
    for panel,ax in zip(("a","b"),axes):
        ax.set_xlabel("Injected "+OBSERVABLES.get(observable,("","b"))[1]+r" [GeV$^{-2}$]")
        ax.text(.02,.97,f"({panel})",transform=ax.transAxes,va="top",fontsize=13)
    fig.suptitle(_label(state,observable,config,"Independent Poisson-toy validation\n"+_inference_notice(entry)),fontsize=13)
    return _save(fig,base.with_name(base.name+"_slope_coverage"),captions,
                 _inference_notice(entry)+". Empirical coverage at the tested injected slopes only. Error bars quantify "
                 "finite-toy binomial simulation uncertainty. This scan does not establish "
                 "uniform coverage at every slope or under omitted nuisance effects. "+SCOPE)


def _stress_description(case, observable):
    if case.get("kind") in ("tail_only", "feed_in_normalization_sensitivity"):
        return f"Outside-range rate ×{float(case.get('scale', 1)):g}"
    tilt = case.get("tilt_dimensionless", case.get("tilt", case.get("tilt_GeVminus2", 0)))
    if case.get("tilt_dimensionless") is not None or observable in ("w", "y"):
        return rf"Normalized-coordinate tilt ${float(tilt):+g}$"
    return rf"Truth-shape tilt ${float(tilt):+g}$ GeV$^{{-2}}$"


def _plot_stress(base,state,observable,entry,config,captions):
    groups = {}
    for row in entry.get("stress_tests", []):
        if row.get("coarse_closure", {}).get("truth") is not None:
            groups.setdefault(row.get("kind", "unknown"), []).append(row)
    coarse = _coarse(entry)
    if not groups or "truth_edges" not in coarse:
        return []
    edges = np.asarray(coarse["truth_edges"], float)
    n = len(edges)-1
    centers = (edges[:-1]+edges[1:])/2
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    titles = {"shape_only": "Shape-only stress with fixed outside-range truth",
              "tail_only": "Outside-range rate stress with fixed in-range truth",
              "combined": "Combined shape and tail stress",
              "truth_shape": "Legacy combined shape and tail stress",
              "feed_in_normalization_sensitivity": "Outside-range rate sensitivity"}
    output = []
    for kind, cases in groups.items():
        fig, ax = plt.subplots(figsize=(7.4, 5.5), layout="constrained")
        for i, case in enumerate(cases):
            closure = case["coarse_closure"]
            fitted = _array(closure, "truth", (n,))
            reference = _array(closure, "reference_truth", (n,))
            if fitted is None or reference is None:
                continue
            ratio = np.divide(fitted, reference, out=np.full(n, np.nan), where=reference>0)
            ax.errorbar(centers, ratio, xerr=np.diff(edges)/2, fmt="o", ms=4.5,
                        color=palette[i % len(palette)], label=_stress_description(case, observable))
        ax.axhline(1, color="#555555", linestyle="--", linewidth=1)
        ax.set_ylabel("Extraction / stressed truth")
        ax.set_xlabel(_axis(observable))
        title = titles.get(kind, kind)
        ax.set_title(_label(state, observable, config, title+"\n"+_inference_notice(entry)), pad=12)
        ax.legend(loc="best", fontsize=10.5)
        suffix = "shape_stress" if kind == "truth_shape" else kind + "_stress"
        caption = (_inference_notice(entry)+". "+title+". Central held-out closure with the training "
                   "response and extraction feed-in fixed. Horizontal bars show bin widths; no vertical "
                   "confidence interval is assigned. ")
        if kind == "shape_only":
            caption += "In-range truth normalization is preserved and outside-range truth weights are unchanged. "
        elif kind in ("combined", "truth_shape"):
            caption += ("Clipped truth-shape reweighting preserves all selected truth normalization and "
                        "also changes outside-range feed-in. This does not isolate in-range response sensitivity. ")
        elif kind in ("tail_only", "feed_in_normalization_sensitivity"):
            caption += "Only the outside-range expected contribution changes; in-range truth stays fixed. "
        output.extend(_save(fig, base.with_name(base.name+"_"+suffix), captions, caption+SCOPE))
    return output


def _plot_bin_coverage(base, state, observable, entry, config, captions):
    rows = (entry.get("inference") or {}).get("bin_coverage") or []
    valid = [r for r in rows if r.get("coverage") is not None and r.get("bin") is not None]
    if not valid:
        return []
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.1), layout="constrained")
    groups = {}
    for row in valid:
        groups.setdefault(row.get("injection", "nominal"), []).append(row)
    for name, points in groups.items():
        points = sorted(points, key=lambda p: p["bin"])
        x = np.array([r["bin"] for r in points]) + 1
        axes[0].errorbar(x, [100*r["coverage"] for r in points],
                        yerr=[100*r.get("coverage_mc_se", 0) for r in points],
                        fmt="o-", capsize=3, label=str(name))
        axes[1].plot(x, [r.get("bias_over_sd", np.nan) for r in points], "o-", label=str(name))
    nominal = float(config.get("inference", {}).get("coverage", .6826894921370859))
    axes[0].axhline(100*nominal, color="#555555", linestyle="--", label=f"Nominal {100*nominal:.1f}%")
    axes[0].set_ylabel("Binwise profile coverage [%]")
    axes[0].set_ylim(0, 105)
    axes[1].axhline(0, color="#555555", linestyle="--")
    axes[1].set_ylabel("Yield bias / toy standard deviation")
    for ax in axes:
        ax.set_xlabel("Truth bin number")
        ax.set_xticks(sorted({int(row["bin"])+1 for row in valid}))
        ax.legend(fontsize=9)
    fig.suptitle(_label(state, observable, config, "Independent Poisson-toy interval checks\n"+_inference_notice(entry)))
    return _save(fig, base.with_name(base.name+"_bin_coverage"), captions,
                 _inference_notice(entry)+". Binwise likelihood-ratio interval coverage at the recorded injected spectra. "
                 "Errors quantify finite-toy simulation uncertainty. This is not simultaneous confidence coverage "
                 "or an exact calibrated confidence belt. "+SCOPE)


def _metric_rows(state,observable,entry):
    coarse=_coarse(entry)
    if not coarse or "truth_edges" not in coarse:
        return []
    edges=np.asarray(coarse["truth_edges"],float)
    n=len(edges)-1
    metrics=coarse.get("metrics") or {}
    truth=_array(coarse,"truth_expected",(n,))
    rows=[]
    for i in range(n):
        row=dict(state=state,observable=observable,bin=i,x_low=edges[i],x_high=edges[i+1])
        if truth is not None:
            row["expected_truth_count"]=truth[i]
        for key,values in metrics.items():
            arr=np.asarray(values)
            if arr.shape==(n,):
                row[key]=arr[i]
        rows.append(row)
    return rows


def write_report(out_dir, results, config):
    """Write standalone tables, captions and diagnostic figures.

    ``results[state][observable]`` must contain an optimization dict, and may
    contain model/validation_model/inference dictionaries. Missing inference
    stays missing. No optimizer or statistical calculation is performed here.
    """
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    plot_dir=out/"plots"
    plot_dir.mkdir(exist_ok=True)
    selected=[]; candidates=[]; metrics=[]; statuses=[]; slopes=[]; coverages=[]
    spectra=[]; profiles=[]; response_mc=[]; adequacies=[]
    bin_intervals=[]; bin_coverages=[]; covariance_metadata=[]; entry_lookup={}
    stresses=[]; stress_bins=[]; fine_checks=[]; selections=[]
    figures=[]; captions={}; sections=[]
    with plt.rc_context(_style()):
        for state,state_results in results.items():
            if not isinstance(state_results,dict):
                continue
            for observable,entry in state_results.items():
                if not isinstance(entry,dict) or "optimization" not in entry:
                    continue
                entry_lookup[(state,observable)] = entry
                optimization=entry.get("optimization") or {}
                inference=entry.get("inference") or {}
                status=_overall_status(entry)
                statuses.append(dict(state=state,observable=observable,status=status,
                                     optimization_status=optimization.get("status","NOT_RUN"),
                                     inference_status=inference.get("status","NOT_RUN"),
                                     component_validation=entry.get("component_validation", {}),
                                     publication_ready=False,
                                     reasons=entry.get("errors",inference.get("gates",[])),
                                     additional_validation_gates=entry.get("additional_validation_gates",{})))
                edges=optimization.get("selected_edges")
                for candidate in optimization.get("candidates",[]):
                    candidates.append(dict(state=state,observable=observable,**candidate))
                section=[f"### {state}, {OBSERVABLES.get(observable,('', '',observable))[2]}",
                         "",f"Status: **{status}**.",""]
                if edges is None:
                    section.extend(["No qualifying binning was selected. No differential result or slope "
                                    "is certified for this state/observable.",""])
                    sections.extend(section)
                    continue
                if len(edges)<2 or np.any(np.diff(np.asarray(edges,float))<=0):
                    raise ValueError("Selected bin edges must be finite and strictly ordered")
                if not np.isfinite(np.asarray(edges,float)).all():
                    raise ValueError("Selected bin edges must be finite")
                selected.append(dict(state=state,observable=observable,status=status,
                                     edges=edges,n_bins=len(edges)-1,
                                     thresholds=optimization.get("thresholds",{})))
                metrics.extend(_metric_rows(state,observable,entry))
                section.extend([f"Selected edges ({_unit(observable) or 'dimensionless'}): {_bins_text(edges)}.",""])
                slope=inference.get("slope") or {}
                if slope:
                    row=dict(state=state,observable=observable,status=slope.get("status","NOT_REPORTED"),
                             overall_status=status)
                    for key in ("b","interval_low","interval_high","amplitude","expected_signal_reco",
                                "deviance","ndof","interval_interpretation","at_lower_boundary",
                                "lower_at_boundary","at_search_boundary","upper_unbounded",
                                "upper_limit_not_found_within_search","upper_search_limited",
                                "profile_set_components_on_display_grid","delta_deviance_threshold","nominal_coverage"):
                        if key in slope:
                            row[key]=slope[key]
                    slopes.append(row)
                    section.extend([f"Forward-folded slope diagnostic: {_number(slope.get('b'))} GeV⁻².",
                                    f"Conditional profile endpoints: [{_number(slope.get('interval_low'))}, "
                                    f"{_number(slope.get('interval_high'))}] GeV⁻².",
                                    "The endpoints are not automatically a validated experimental uncertainty.",""])
                for toy in inference.get("slope_toys",[]):
                    coverages.append(dict(state=state,observable=observable,**toy))
                for point in slope.get("profile",[]):
                    if isinstance(point,dict):
                        profiles.append(dict(state=state,observable=observable,**point))
                if inference.get("response_mc"):
                    response_mc.append(dict(state=state,observable=observable,**inference["response_mc"]))
                adequacy=inference.get("model_adequacy") or {}
                if adequacy:
                    adequacies.append(dict(state=state,observable=observable,**adequacy))
                    section.extend(["Expected-template exponential-model adequacy: "
                                    f"{adequacy.get('status','NOT_RUN')}. Refitted-toy tail fraction "
                                    f"{_number(adequacy.get('expected_template_discrepancy_tail_probability'))}. "
                                    "This compares expected-template discrepancy with projected statistical "
                                    "fluctuations. It is not an observed-data goodness-of-fit p value.",""])
                extraction=inference.get("coarse") or {}
                if extraction.get("truth_edges") is not None and extraction.get("fitted_truth") is not None:
                    xedges=np.asarray(extraction["truth_edges"],float)
                    nbins=len(xedges)-1
                    conversion, conversion_meta = _cross_section(entry, observable, xedges)
                    alternate = conversion_meta.get("alternate_factors")
                    if alternate is not None:
                        alternate = np.asarray(alternate, float)
                        if alternate.shape != (nbins,) or not np.isfinite(alternate).all() or np.any(alternate <= 0):
                            raise ValueError("Invalid alternate cross-section conversion")
                    for i in range(nbins):
                        row=dict(state=state,observable=observable,bin=i,
                                 x_low=xedges[i],x_high=xedges[i+1],
                                 status=extraction.get("status","NOT_REPORTED"),
                                 quantity=conversion_meta.get("quantity"),
                                 quantity_units=conversion_meta.get("units"),
                                 x_units=_unit(observable), conversion_factor=conversion[i],
                                 interval_interpretation="Central toy sampling quantiles, not confidence limits")
                        for key in ("fitted_truth","reference_truth","toy_mean","toy_lower","toy_upper","sampling_sd","bias_over_sd","boundary_fractions"):
                            if extraction.get(key) is not None:
                                arr=_array(extraction,key,(nbins,))
                                row[key]=arr[i]
                        for key in ("fitted_truth","reference_truth","toy_mean","toy_lower","toy_upper","sampling_sd"):
                            if key in row:
                                row[key+"_cross_section"] = row[key]*conversion[i]
                                if observable in ("pt2", "abs_t_via_pt2") and conversion_meta.get("units") != "events":
                                    row[key+"_pb_per_GeV2"] = row[key]*conversion[i]
                                if alternate is not None:
                                    row[key+"_alternate_cross_section"] = row[key]*alternate[i]
                        if alternate is not None:
                            row.update(alternate_quantity=conversion_meta.get("alternate_quantity"),
                                       alternate_units=conversion_meta.get("alternate_units"),
                                       alternate_conversion_factor=alternate[i])
                        spectra.append(row)
                    for interval in inference.get("bin_intervals", []):
                        row=dict(state=state, observable=observable, **interval)
                        i=int(row["bin"])
                        row.update(quantity=conversion_meta.get("quantity"), quantity_units=conversion_meta.get("units"))
                        for key in ("truth_yield", "interval_low", "interval_high"):
                            if row.get(key) is not None:
                                row[key+"_cross_section"] = row[key]*conversion[i]
                        bin_intervals.append(row)
                    for coverage in inference.get("bin_coverage", []):
                        bin_coverages.append(dict(state=state, observable=observable, **coverage))
                    if extraction.get("covariance") is not None:
                        covariance=_array(extraction,"covariance",(nbins,nbins))
                        stem=_safe_stem(state, observable)
                        files = {"truth_yield": stem+"_truth_yield_covariance.csv",
                                 "cross_section": stem+"_cross_section_covariance.csv"}
                        np.savetxt(out/files["truth_yield"],covariance,delimiter=",")
                        np.savetxt(out/files["cross_section"],covariance*conversion[:,None]*conversion[None,:],delimiter=",")
                        if alternate is not None:
                            files["alternate_cross_section"] = stem+"_alternate_cross_section_covariance.csv"
                            np.savetxt(out/files["alternate_cross_section"],covariance*alternate[:,None]*alternate[None,:],delimiter=",")
                        calibration=inference.get("response_mc") or {}
                        if calibration.get("covariance") is not None:
                            mc_covariance=_array(calibration,"covariance",(nbins,nbins))
                            files["response_mc_truth_yield"] = stem+"_response_mc_truth_yield_covariance.csv"
                            files["response_mc_cross_section"] = stem+"_response_mc_cross_section_covariance.csv"
                            np.savetxt(out/files["response_mc_truth_yield"],mc_covariance,delimiter=",")
                            np.savetxt(out/files["response_mc_cross_section"],mc_covariance*conversion[:,None]*conversion[None,:],delimiter=",")
                        covariance_metadata.append(dict(state=state, observable=observable,
                            status=extraction.get("status", "NOT_REPORTED"), files=files,
                            response_mc_stage_status=calibration.get("status","NOT_RUN"),
                            quantity=conversion_meta.get("quantity"), quantity_units=conversion_meta.get("units"),
                            truth_edges=xedges, conversion_factors=conversion,
                            alternate_quantity=conversion_meta.get("alternate_quantity"),
                            alternate_units=conversion_meta.get("alternate_units"),
                            alternate_conversion_factors=alternate,
                            interpretation="Sampling covariance of the nonnegative truth-yield estimator; not an interval"))
                meta=entry.get("state_config") or {}
                if meta:
                    sigma=meta.get("sigma")
                    mass_halfwidth=(float(sigma)*float(config.get("reco_mass_window_sigma",3))
                                    if sigma is not None else None)
                    wrange=config.get("w_range",[20,100])
                    selections.append(dict(state=state,observable=observable,q2_max=meta.get("q2_max"),
                                      truth_w_min=wrange[0],truth_w_max=wrange[1],
                                      reco_w_min=wrange[0],reco_w_max=wrange[1],
                                      truth_mass_halfwidth_GeV=config.get("truth_mass_window",.2),
                                      reco_mass_halfwidth_GeV=mass_halfwidth,
                                      fixed_selection_sigma_GeV=sigma,
                                      sigma_epz_center_GeV=config.get("sigma_epz_center",18),
                                      sigma_epz_halfwidth_GeV=config.get("sigma_epz_halfwidth",1.8),
                                      category=meta.get("category"),q2_note="Truth-fiducial selection",
                                      visible_xsec_fb=meta.get("visible_xsec_fb"),br_ee=meta.get("br_ee"),
                                      target_lumi_invfb=config.get("target_lumi_invfb",2.5)))
                    section.extend([f"Configured Q² maximum: {_number(meta.get('q2_max'))} GeV². "
                                    f"Visible cross section: {_number(meta.get('visible_xsec_fb'),6)} fb. "
                                    f"Dielectron branching fraction: {_number(meta.get('br_ee'),6)}. "
                                    f"Reconstructed mass half-window: {_number(mass_halfwidth)} GeV, "
                                    "using the fixed selection width rather than a fitted resolution.",""])
                for case in entry.get("stress_tests",[]):
                    fit=case.get("fit") or {}
                    row=dict(state=state,observable=observable,kind=case.get("kind"),
                             tilt_GeVminus2=case.get("tilt_GeVminus2"),
                             tilt_dimensionless=case.get("tilt_dimensionless"),tilt=case.get("tilt"),
                             tilt_units=case.get("tilt_units"),scale=case.get("scale"),
                             status=fit.get("status",(case.get("coarse_closure") or {}).get("status","NOT_RUN")),b=fit.get("b"),
                             closure_passes=case.get("closure_passes"),
                             shift_over_nominal_halfwidth=case.get("shift_over_nominal_halfwidth"),
                             sensitivity_flag=case.get("sensitivity_flag"),
                             normalization=case.get("normalization"))
                    stresses.append(row)
                    closure=case.get("coarse_closure") or {}
                    if closure.get("truth") is not None:
                        for i,truth in enumerate(closure["truth"]):
                            rec=dict(state=state,observable=observable,kind=case.get("kind"),
                                     tilt_GeVminus2=case.get("tilt_GeVminus2"),
                                     tilt_dimensionless=case.get("tilt_dimensionless"),tilt=case.get("tilt"),
                                     tilt_units=case.get("tilt_units"),scale=case.get("scale"),bin=i,fitted_truth=truth,
                                     status=closure.get("status"))
                            for key in ("reference_truth","relative_bias","bias_over_nominal_sampling_sd"):
                                if closure.get(key) is not None:
                                    rec[key]=closure[key][i]
                            stress_bins.append(rec)
                fine=entry.get("fine_cell_check") or {}
                if fine:
                    fit=fine.get("fit") or {}
                    fine_checks.append(dict(state=state,observable=observable,
                        step_GeV2=fine.get("step_GeV2"),b=fit.get("b"),fit_status=fit.get("status"),
                        shift_over_nominal_halfwidth=fine.get("shift_over_nominal_halfwidth"),
                        threshold=fine.get("threshold"),passes=fine.get("passes")))
                base=plot_dir/_safe_stem(state,observable)
                old_caption_keys=set(captions)
                figures.extend(_plot_response(base,state,observable,entry,config,captions))
                figures.extend(_plot_metrics(base,state,observable,entry,config,captions))
                figures.extend(_plot_coarse(base,state,observable,entry,config,captions))
                figures.extend(_plot_slope(base,state,observable,entry,config,captions))
                figures.extend(_plot_coverage(base,state,observable,entry,config,captions))
                figures.extend(_plot_stress(base,state,observable,entry,config,captions))
                figures.extend(_plot_bin_coverage(base,state,observable,entry,config,captions))
                scaling=entry.get("cross_section_scaling") or {}
                if observable == "y" and scaling.get("alternate_factors") is not None:
                    alternate_entry=dict(entry)
                    alternate_entry["cross_section_scaling"]=dict(scaling,
                        factors=scaling["alternate_factors"], quantity=scaling.get("alternate_quantity"),
                        units=scaling.get("alternate_units"),
                        ylabel=r"$(1/\Phi)\,d\sigma_{ep}/dy_\Upsilon$ [pb]")
                    alternate_base=base.with_name(base.name+"_flux_normalized")
                    figures.extend(_plot_coarse(alternate_base,state,observable,alternate_entry,config,captions))
                    captions[alternate_base.name+"_differential"] += (
                        " This companion spectrum divides the ep rapidity spectrum by the total selected W photon flux. "
                        "It is not a model-independent photon-proton cross section at a rapidity value.")
                if meta:
                    for key in set(captions)-old_caption_keys:
                        captions[key]+= (f" State-specific truth Q² maximum {_number(meta.get('q2_max'))} GeV²; "
                                         f"truth and reco W range {_bins_text(config.get('w_range',[20,100]))} GeV. "
                                         f"Fixed reconstructed mass half-window {_number(mass_halfwidth)} GeV.")
                sections.extend(section)
    # Every numerical table distinguishes the final validation status from the
    # component status. The generic status column is deliberately conservative
    # for old notebook consumers that read only that column.
    tables=[(selected,"selected"),(candidates,"candidate"),(metrics,"response"),(statuses,"inference"),
            (slopes,"slope"),(coverages,"slope_toys"),(spectra,"coarse"),(profiles,"slope"),
            (response_mc,"response_mc"),(adequacies,"model_adequacy"),(stresses,"stress"),
            (stress_bins,"stress"),(fine_checks,"fine"),(selections,"configuration"),
            (bin_intervals,"bin_intervals"),(bin_coverages,"bin_coverage"),(covariance_metadata,"coarse")]
    for rows, component in tables:
        for row in rows:
            entry=entry_lookup[(row["state"],row["observable"])]
            inference=entry.get("inference") or {}
            overall=_overall_status(entry)
            if component == "candidate":
                stage="PASS_DESIGN" if row.get("passes") else "FAIL_DESIGN"
            elif component == "selected":
                stage=(entry.get("optimization") or {}).get("status","NOT_RUN")
            elif component == "inference":
                stage=inference.get("status","NOT_RUN")
            elif component == "configuration":
                stage="CONFIGURED"
            else:
                current=inference.get(component) or {}
                default=current.get("status","RECORDED") if isinstance(current,dict) else "RECORDED"
                stage=row.get("stage_status",row.get("fit_status",row.get("status",default)))
            selection_meta=entry.get("selection_metadata") or {}
            proton_policy=entry.get("proton_selection", selection_meta.get("proton_selection", {"mode":"unrecorded"}))
            row.update(proton_selection=proton_policy,
                       proton_selection_mode=proton_policy.get("mode", "unrecorded"))
            if component == "configuration":
                row["proton_selection_cutflow"]=selection_meta.get("proton_selection_cutflow", {})
            row.update(status=overall,overall_status=overall,stage_status=stage,
                       failed_checks=_failed_checks(entry),publication_ready=False)
    for row in covariance_metadata:
        _write_json(out/(_safe_stem(row["state"],row["observable"])+"_covariance_metadata.json"),row)
    _write_csv(out/"bin_intervals.csv",bin_intervals,fields=["state","observable","bin","status","overall_status","stage_status"])
    _write_csv(out/"bin_coverage.csv",bin_coverages,fields=["state","observable","bin","status","overall_status","stage_status"])
    _write_csv(out/"covariance_metadata.csv",covariance_metadata,fields=["state","observable","status","overall_status","stage_status"])
    _write_csv(out/"selected_bins.csv",selected,fields=list(dict.fromkeys(
        ["state","observable","status","edges","n_bins","thresholds"])))
    _write_csv(out/"selected_bin_metrics.csv",metrics,
               fields=None if metrics else ["state","observable","bin","x_low","x_high"])
    _write_csv(out/"all_candidates.csv",candidates,
               fields=None if candidates else ["state","observable","edges","passes"])
    _write_csv(out/"validation_status.csv",statuses,
               fields=None if statuses else ["state","observable","status","publication_ready"])
    _write_csv(out/"slope_summary.csv",slopes,
               fields=None if slopes else ["state","observable","status","b","interval_low","interval_high"])
    _write_csv(out/"slope_coverage.csv",coverages,
               fields=None if coverages else ["state","observable","input_b","coverage"])
    _write_csv(out/"differential_results.csv",spectra,
               fields=None if spectra else ["state","observable","bin","fitted_truth","interval_interpretation"])
    _write_csv(out/"slope_profiles.csv",profiles,
               fields=None if profiles else ["state","observable","b","delta_deviance"])
    _write_csv(out/"response_mc_summary.csv",response_mc,
               fields=None if response_mc else ["state","observable","n_bootstrap","b_std"])
    _write_csv(out/"model_adequacy.csv",adequacies,
               fields=None if adequacies else ["state","observable","status","expected_template_discrepancy_tail_probability"])
    _write_csv(out/"stress_summary.csv",stresses,
               fields=None if stresses else ["state","observable","kind","status","b"])
    _write_csv(out/"shape_stress_closure.csv",stress_bins,
               fields=None if stress_bins else ["state","observable","tilt_GeVminus2","bin","relative_bias"])
    _write_csv(out/"fine_cell_checks.csv",fine_checks,
               fields=None if fine_checks else ["state","observable","step_GeV2","b","passes"])
    _write_csv(out/"selection_settings.csv",selections,
               fields=None if selections else ["state","observable","q2_max","fixed_selection_sigma_GeV"])
    if config.get("synthetic",False):
        captions={key:"SYNTHETIC SOFTWARE TEST. Not a physics projection. "+value
                  for key,value in captions.items()}
    _write_json(out/"figure_captions.json",captions)
    lumi=float(config.get("target_lumi_invfb",2.5))
    summary=["# Updated-yield binning and response validation", "",
             f"Projection luminosity: **{lumi:g} fb⁻¹** at "
             f"{config.get('electron_energy',9)}×{config.get('proton_energy',275)} GeV.","",
             SCOPE,"",
             "Binning is selected with training events and frozen before independent held-out "
             "validation. Design thresholds do not prove confidence-interval coverage. "
             "Missing or failed inference remains explicitly unverified.","",
             "W reports a photon-flux-weighted bin-average cross section in pb, with no extra W-bin-width division. "
             "Rapidity reports inclusive ep dσ/dy; its separately labeled companion divides by total photon flux. "
             "W and rapidity do not use an exponential slope fit. Binwise conditional profile coverage is tested separately.","",
             "The pt2 study targets true pT². The abs_t_via_pt2 study uses reconstructed "
             "pT² as an estimator of true |t|. Their slopes have different definitions.","",
             "Mass-window widths are fixed selection parameters, not the fitted mass resolutions. "
             "State-dependent Q² and ID selections must remain in the appendix description.","",
             "No result in this report is automatically marked ready for a full experimental "
             "measurement. Conditional validation status and omitted effects must be stated.","",
             *sections,"## Figure captions","",]
    if config.get("synthetic",False):
        summary[2:2]=["**SYNTHETIC SOFTWARE TEST. These are software-validation results, not physics projections.**",""]
    for stem,caption in captions.items():
        summary.extend([f"- `{stem}`: {caption}"])
    (out/"summary.md").write_text("\n".join(summary)+"\n")
    synthetic_caption=(r"\textbf{SYNTHETIC SOFTWARE TEST.} Not a physics projection. "
                       if config.get("synthetic",False) else "")
    tex=["% Generated diagnostics. Review scope and statuses before inclusion in the note."]
    if config.get("synthetic",False):
        tex.extend([r"\noindent\textbf{SYNTHETIC SOFTWARE TEST.} These tables test software "
                    r"and must not be presented as physics projections.",""])
    tex.extend([
         r"\begin{table}[htbp]",r"\centering",r"\small",
         r"\begin{tabular}{lllp{0.30\linewidth}}",r"\toprule",
         r"State & Observable & Edges (observable units) & Status \\",r"\midrule"])
    for row in statuses:
        match=next((s for s in selected if s["state"]==row["state"] and s["observable"]==row["observable"]),None)
        edges_text=_bins_text(match["edges"]) if match else "no qualifying solution"
        obs_tex=OBSERVABLES.get(row["observable"],(row["observable"],))[0]
        tex.append(f"{_tex(row['state'])} & {obs_tex} & {_tex(edges_text)} & {_tex(row['status'])} \\\\")
    tex.extend([r"\bottomrule",r"\end{tabular}",
                r"\caption{"+synthetic_caption+r"Training-selected binning for the updated $9\times275$ GeV projection at $"
                +f"{lumi:g}"+r"\,\mathrm{fb}^{-1}$. Status refers to conditional signal-only "
                r"validation, not complete experimental coverage. Failed or unavailable solutions "
                r"must not be quoted as validated cross sections.}",r"\end{table}",""])
    if metrics:
        tex.extend([r"\begin{table}[htbp]",r"\centering",r"\small",
                    r"\begin{tabular}{lllrrrr}",r"\toprule",
                    r"State & Observable & Bin (observable units) & $\mu_{\rm reco}$ & Purity & Stability & Diag. acc. \\",
                    r"\midrule"])
        for row in metrics:
            obs_tex=OBSERVABLES.get(row["observable"],(row["observable"],))[0]
            vals=[_tex(row["state"]),obs_tex,
                  f"{row['x_low']:g}--{row['x_high']:g}",
                  _number(row.get("expected_reco_signal"),3),
                  _number(row.get("purity"),3),
                  _number(row.get("selected_stability"),3),
                  _number(row.get("diagonal_acceptance"),3)]
            tex.append(" & ".join(vals)+r" \\")
        tex.extend([r"\bottomrule",r"\end{tabular}",
                    r"\caption{"+synthetic_caption+r"Training-bin diagnostics at the requested luminosity. The expected "
                    r"signal count excludes outside-range feed-in; purity includes it in the "
                    r"reconstructed denominator. Stability is conditional on selection, whereas "
                    r"diagonal acceptance uses the generated truth denominator.}",r"\end{table}",""])
    if slopes:
        tex.extend([r"\begin{table}[htbp]",r"\centering",r"\small",
                    r"\begin{tabular}{llllp{0.25\linewidth}}",r"\toprule",
                    r"State & Observable & $b$ [GeV$^{-2}$] & Profile endpoints & Overall status \\",r"\midrule"])
        for row in slopes:
            obs_tex=OBSERVABLES.get(row["observable"],(row["observable"],))[0]
            tex.append(_tex(row["state"])+" & "+obs_tex+" & "+_number(row.get("b"),4)+" & ["+
                       _number(row.get("interval_low"),4)+", "+_number(row.get("interval_high"),4)+"] & "+
                       _tex(row.get("overall_status","NOT_VALIDATED"))+r" \\")
        tex.extend([r"\bottomrule",r"\end{tabular}",
                    r"\caption{"+synthetic_caption+r"Slope diagnostics from the frozen-response Poisson fit to held-out "
                    r"expected means. These are conditional profile limits, not certified "
                    r"experimental confidence intervals. Quote them only with the accompanying "
                    r"coverage, bias and boundary diagnostics.}",r"\end{table}",""])
    tex.extend(["% Response definitions used by the numerical outputs.",
                r"The response is $R_{ri}=N^{\mathrm{selected}}_{ri}/N_i^{\mathrm{generated}}$, "
                r"including missed reconstruction in the truth denominator. Purity includes "
                r"migration from outside the fitted truth range. Selected migration stability "
                r"and diagonal acceptance are reported separately because they have different denominators.","",
                r"The forward-folded slope fit uses reconstructed Poisson counts and finite-bin "
                r"integrals of the exponential model. Its profile limits are conditional forecasts; "
                r"independent toy coverage is reported at the tested slopes. Bands on extracted "
                r"differential spectra are central toy sampling quantiles, not confidence limits.","",
                r"These studies retain truth-fiducial mass, $W$ and $Q^2$ selections. "
                r"Outside-variable-range feed-in is modeled within that truth fiducial region. "
                r"Mass-component extraction, backgrounds and systematic uncertainties are not "
                r"included. Finite response-MC variations are reported separately and are not "
                r"silently added to the projected-data interval."])
    (out/"appendix_tables.tex").write_text("\n".join(tex)+"\n")
    manifest={"files":[str(out/name) for name in ("summary.md","appendix_tables.tex",
               "selected_bins.csv","selected_bin_metrics.csv","all_candidates.csv",
               "validation_status.csv","slope_summary.csv","slope_coverage.csv","differential_results.csv",
               "slope_profiles.csv","response_mc_summary.csv","model_adequacy.csv","stress_summary.csv",
               "shape_stress_closure.csv","fine_cell_checks.csv","selection_settings.csv","figure_captions.json",
               "bin_intervals.csv","bin_coverage.csv","covariance_metadata.csv")],
              "figures":figures,"publication_ready":False,"scope":SCOPE,
              "status_contract":"status and overall_status are final case status; stage_status is component status",
              "cases":statuses,
              "synthetic":bool(config.get("synthetic",False))}
    _write_json(out/"report_manifest.json",manifest)
    return manifest
