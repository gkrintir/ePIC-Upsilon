#!/usr/bin/env python3
"""Render single-sample t-estimator benchmark pages in the reviewed-note style.

Computations come from validation/benchmark_t_estimators.py without formula
changes. These are finite-range response and closure diagnostics, not a repair
of any t estimator or a new unfolded measurement. The projected estimator keeps
its own selection and absolute corrected weights; statistical response metrics
continue to use ordinary event weights, exactly as in the source implementation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

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
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import benchmark_t_estimators as benchmark
from note_style import BLUE, set_style, finish

METHOD_KEYS = ["eX", "BABE", "eXBABE", "proj", "pT2"]


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def subtitle(args, additional=""):
    label = args.energy_label.replace("x", r"$\times$")
    result = (rf"$\Upsilon({getattr(args, 'state', '1S')})$  |  {label}  |  ${args.w_min:g}\leq W_{{\rm true,reco}}\leq {args.w_max:g}$ GeV"
              rf"  |  $0\leq Q^2_{{\rm true}}\leq {args.q2_cut:g}$ GeV$^2$" + "\n" +
              rf"${args.t_edges[0]:g}\leq |t|_{{\rm true,reco}}\leq {args.t_edges[-1]:g}$ GeV$^2$"
              rf"  |  $\mathcal{{L}}={args.lumi_invfb:g}$ fb$^{{-1}}$; normalized input weights")
    if additional:
        result += "  |  " + additional
    return result


def render_metrics(results, args, output):
    fields = [("average_efficiency", "Average efficiency"),
              ("min_purity_stability", "Minimum purity/stability"),
              ("median_delta_t", r"Median $\Delta t$ [GeV$^2$]"),
              ("core_sigma_t", r"Core $\sigma(\Delta t)$ [GeV$^2$]")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    x = np.arange(len(METHOD_KEYS))
    for ax, (field, ylabel) in zip(axes.flat, fields):
        values = np.array([results[key]["summary"][field] for key in METHOD_KEYS], dtype=float)
        valid = np.isfinite(values)
        ax.bar(x[valid], values[valid], width=0.65, color=BLUE, alpha=0.85)
        ax.set_xticks(x, [benchmark.METHODS[key]["label"] for key in METHOD_KEYS])
        ax.set_ylabel(ylabel)
        ax.axhline(0, color="gray", lw=0.7)
        ax.grid(axis="x", visible=False)
        if field in ["average_efficiency", "min_purity_stability"]:
            ax.set_ylim(0, max(0.05, float(np.nanmax(values)) * 1.25) if np.any(valid) else 1)
        if not np.all(valid):
            for pos in x[~valid]:
                ax.text(pos, 0.05, "no entries", ha="center", va="bottom", rotation=90, fontsize=8,
                        transform=ax.get_xaxis_transform())
    finish(fig, output, "t_methods_metrics", r"Reconstructed $|t|$ estimator response metrics", subtitle(args))


def render_residuals(results, args, output):
    fig, axes = plt.subplots(3, 2, figsize=(13, 12))
    rows, display_audit = [], {}
    for ax, method in zip(axes.flat, METHOD_KEYS):
        payload = results[method]["residual"]
        delta, weights = payload["delta"], payload["weights"]
        ax.set_title(benchmark.METHODS[method]["label"])
        ax.set_xlabel(r"$|t|_{\rm reco}-|t|_{\rm true}$ [GeV$^2$]")
        ax.set_ylabel("Projected events / bin")
        ax.axvline(0, color="gray", ls="--", lw=1)
        if len(delta):
            lo, hi = np.nanpercentile(delta, [1, 99])
            if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
                lo, hi = -2.0, 2.0
            pad = 0.1 * (hi - lo)
            edges = np.linspace(lo - pad, hi + pad, 60)
            counts, _ = np.histogram(delta, edges)
            weighted, _ = np.histogram(delta, edges, weights=weights)
            weighted2, _ = np.histogram(delta, edges, weights=weights**2)
            ax.stairs(weighted, edges, color=BLUE, lw=1.7)
            ax.axvline(payload["median"], color="red", lw=1.3)
            ax.set_yscale("log")
            positive = weighted[weighted > 0]
            if len(positive):
                ax.set_ylim(float(positive.min()) * 0.5, float(positive.max()) * 4)
            ax.set_xlim(edges[0], edges[-1])
            ax.text(0.98, 0.95, f"Median = {payload['median']:.3g}\nCore width = {payload['core']:.3g} GeV$^2$",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9)
            below, above = delta < edges[0], delta > edges[-1]
            display_audit[method] = {"display_range_GeV2": [float(edges[0]), float(edges[-1])],
                                    "range_rule": "1st to 99th percentiles of selected residuals plus 10% padding (source convention)",
                                    "underflow_mc_entries": int(below.sum()), "overflow_mc_entries": int(above.sum()),
                                    "underflow_sum_weight": float(weights[below].sum()), "overflow_sum_weight": float(weights[above].sum()),
                                    "selected_sum_weight": float(weights.sum())}
            for i in range(len(weighted)):
                rows.append({"method": method, "low_GeV2": float(edges[i]), "high_GeV2": float(edges[i+1]),
                             "mc_count": int(counts[i]), "projected_events": float(weighted[i]), "finite_mc_variance": float(weighted2[i])})
        else:
            ax.text(0.5, 0.5, "No selected entries", transform=ax.transAxes, ha="center")
            display_audit[method] = {"selected_mc_entries": 0}
    axes.flat[-1].set_axis_off()
    handles = [Line2D([], [], color=BLUE, lw=1.7), Line2D([], [], color="gray", ls="--"), Line2D([], [], color="red")]
    finish(fig, output, "t_methods_residuals", r"Reconstructed $|t|$ estimator residuals", subtitle(args),
           handles=handles, labels=["Selected signal MC", "Zero residual", "Median residual"], legend_ncol=3)
    benchmark.write_csv(output / "t_methods_residual_bins.csv", rows,
                        fieldnames=["method", "low_GeV2", "high_GeV2", "mc_count", "projected_events", "finite_mc_variance"])
    return display_audit


def render_closure(results, args, output):
    fig, axes = plt.subplots(3, 2, figsize=(13, 12))
    for ax, method in zip(axes.flat, METHOD_KEYS):
        rows = results[method]["bins"]
        center = np.asarray([row["t_center"] for row in rows])
        half_width = np.asarray([(row["t_high"] - row["t_low"]) / 2 for row in rows])
        closure = np.asarray([row["closure"] for row in rows])
        valid = np.isfinite(closure)
        ax.axhline(1, color="gray", ls="--", lw=1)
        ax.errorbar(center[valid], closure[valid], xerr=half_width[valid], fmt="o", color=BLUE, ms=4, lw=1.3)
        ax.set_title(benchmark.METHODS[method]["label"])
        ax.set_xlabel(r"True $|t|$ [GeV$^2$]")
        ax.set_ylabel("Unfolded / truth")
        # Keep the original 0--2 framing when it contains all points; expose any excursion.
        max_value = float(np.max(closure[valid])) if np.any(valid) else 2
        min_value = float(np.min(closure[valid])) if np.any(valid) else 0
        ax.set_ylim(min(0, min_value * 1.1), max(2, max_value * 1.1))
        ax.set_xlim(args.t_edges[0], args.t_edges[-1])
        if not np.any(valid):
            ax.text(0.5, 0.5, "No finite closure points", ha="center", transform=ax.transAxes)
    axes.flat[-1].set_axis_off()
    handles = [Line2D([], [], color=BLUE, marker="o", ls="none"), Line2D([], [], color="gray", ls="--")]
    finish(fig, output, "t_methods_closure", r"Independent-sample $|t|$ unfolding closure",
           subtitle(args, f"{args.n_iterations} Bayesian iterations"), handles=handles,
           labels=["Closure ratio; horizontal bars show bin widths", "Unity"], legend_ncol=2)


def main():
    # The preserved source adds hash((sample_name, method)) to its RNG seed.
    # Fix hash randomization before loading data so --seed is reproducible.
    if os.environ.get("PYTHONHASHSEED") in (None, "random"):
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = "0"
        os.execve(sys.executable, [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]], env)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", choices=list(benchmark.UPSILON), default="1S")
    parser.add_argument("--input", required=True)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--energy-label", required=True)
    parser.add_argument("--electron-energy", type=float, default=9)
    parser.add_argument("--w-min", type=float, default=20)
    parser.add_argument("--w-max", type=float, default=70)
    parser.add_argument("--q2-cut", type=float, default=1)
    parser.add_argument("--lumi-invfb", type=float, default=10)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--t-edges", default="0,0.2,0.5,1,5")
    parser.add_argument("--train-fraction", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--n-iterations", type=int, default=4)
    args = parser.parse_args()
    args.t_edges = benchmark.parse_edges(args.t_edges)
    if len(args.t_edges) < 2 or np.any(~np.isfinite(args.t_edges)) or np.any(np.diff(args.t_edges) <= 0) or args.t_edges[0] < 0:
        parser.error("t-edges must be finite, strictly increasing, and nonnegative")
    if not 0 < args.train_fraction < 1 or args.n_iterations < 1 or args.seed < 0:
        parser.error("Require 0 < train-fraction < 1, n-iterations >= 1, and seed >= 0")
    if args.w_min >= args.w_max or args.electron_energy <= 0 or args.q2_cut <= 0 or args.lumi_invfb <= 0:
        parser.error("Require an ordered W interval and positive energy, Q2 limit, and luminosity")
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    set_style()
    legacy_args = SimpleNamespace(state=args.state, tree=args.tree, truth_mass_window=0.2, q2_cut=args.q2_cut,
                                  state_norm_mode="raw", lumi_invfb=args.lumi_invfb, t_edges=args.t_edges,
                                  sigma_epz_halfwidth=None, sigma_epz_halfwidth_fraction=0.2,
                                  reco_mass_min=8.0, reco_mass_max=12.0, seed=args.seed,
                                  train_fraction=args.train_fraction, n_iterations=args.n_iterations)
    data, truth_state, source_meta = benchmark.load_sample(args.input, args.energy_label, args.electron_energy, None, legacy_args)
    if np.any(data["weight_visible"] < 0):
        raise RuntimeError("Statistical response metrics require nonnegative ordinary event weights")
    with uproot.open(args.input) as source:
        original_entries = source[args.tree].num_entries
    results, summary_rows, per_bin_rows, selection_details = {}, [], [], {}
    for method in METHOD_KEYS:
        summary, bins, residual = benchmark.analyze_method(data, truth_state, args.energy_label,
                            args.electron_energy, method, legacy_args, args.w_min, args.w_max)
        results[method] = {"summary": summary, "bins": bins, "residual": residual}
        summary_rows.append(summary)
        per_bin_rows.extend(bins)
        reco_t = benchmark.make_reco_t(data, method)
        truth_range = (truth_state & (data["true_t"] >= args.t_edges[0]) & (data["true_t"] <= args.t_edges[-1])
                       & (data["true_W"] >= args.w_min) & (data["true_W"] <= args.w_max))
        baseline_reco = (truth_range & (data["reco_mass"] >= 8) & (data["reco_mass"] <= 12)
                         & (data["reco_W"] >= args.w_min) & (data["reco_W"] <= args.w_max)
                         & (np.abs(data["sigma_epz"] - 2 * args.electron_energy) <= 0.2 * args.electron_energy))
        baseline_reco &= proton_tag_mask(data)
        method_baseline = baseline_reco.copy()
        if method == "proj":
            method_baseline &= (data["weight_proj"] > 0) & (data["reco_t_eXBABE"] > 0)
        selection_details[method] = {"baseline_selected_before_method_requirements": int(baseline_reco.sum()),
                                     "after_method_requirements_before_reco_t_range": int(method_baseline.sum()),
                                     "reco_t_below_range": int(np.sum(method_baseline & (reco_t < args.t_edges[0]))),
                                     "reco_t_above_range": int(np.sum(method_baseline & (reco_t > args.t_edges[-1]))),
                                     "true_t_equal_last_edge": int(np.sum(truth_range & (data["true_t"] == args.t_edges[-1]))),
                                     "reco_t_equal_last_edge": int(np.sum(method_baseline & (reco_t == args.t_edges[-1]))),
                                     "effective_rng_seed": args.seed + hash((args.energy_label, method)) % 100000}
    benchmark.write_csv(output / "t_methods_summary.csv", summary_rows)
    benchmark.write_csv(output / "t_methods_per_bin.csv", per_bin_rows)
    render_metrics(results, args, output)
    residual_audit = render_residuals(results, args, output)
    render_closure(results, args, output)
    source_path = Path(benchmark.__file__).resolve()
    audit = {"proton_selection": proton_selection_metadata(),
             "input": str(Path(args.input).resolve()), "tree": args.tree, "state": args.state, "energy_label": args.energy_label,
             "electron_energy_GeV": args.electron_energy, "input_tuple_entries": original_entries,
             "source_load_metadata": source_meta, "source_implementation": str(source_path),
             "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
             "source_functions": ["load_sample", "analyze_method"], "computation_changes": "Explicit declared proton selection in reconstructed numerator; truth denominator unchanged",
             "normalization_mode": "raw: normalized input weight and absolute weight_proj used without further scaling",
             "normalization_caveat": "Luminosity is a label only; input event weights must have that normalization",
             "luminosity_label_invfb": args.lumi_invfb, "true_mass_window_GeV": [benchmark.UPSILON[args.state]['mass'] - 0.2, benchmark.UPSILON[args.state]['mass'] + 0.2],
             "true_Q2_GeV2": [0, args.q2_cut], "true_and_reco_W_GeV": [args.w_min, args.w_max],
             "true_and_reco_t_edges_GeV2": args.t_edges, "reco_mass_GeV": [8, 12],
             "sigma_epz_GeV": [1.8 * args.electron_energy, 2.2 * args.electron_energy],
             "finite_entry_rule": "Source removes rows with any nonfinite required branch before all method comparisons",
             "finite_rows_removed": original_entries - len(data["weight"]),
             "projected_method_requirement": "weight_proj > 0 and reco_t_eXBABE > 0, in addition to all other cuts",
             "statistical_weights": "ordinary visible weight for response, residuals, and closure, including projected method",
             "projected_corrected_yield": "sum of absolute weight_proj, retained separately as sum_reco_corrected",
             "purity_stability_caveat": "Minimum includes empty bins as zero, exactly as implemented in source",
             "closure_caveat": "Independent randomized train/test split; no closure-ratio uncertainty is computed by the source",
             "interpretation": "Finite-range diagnostic; not an estimator correction, physical validation, or new measured spectrum",
             "train_fraction": args.train_fraction, "n_iterations": args.n_iterations, "seed": args.seed,
             "python_hash_seed": os.environ.get("PYTHONHASHSEED"), "methods": selection_details,
             "residual_display": residual_audit}
    (output / "note_t_studies_audit.json").write_text(json.dumps(json_safe(audit), indent=2, allow_nan=False) + "\n")
    print(json.dumps(json_safe({"output_dir": str(output), "summary": summary_rows}), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
