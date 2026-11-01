#!/usr/bin/env python3
"""
Study W_gamma_p reconstruction resolution for exclusive Upsilon tuples.

Uses true_W and reco_W in UpsilonTree. The resolution is computed for reconstructed
selected events only, because events with no reco candidate have no reco_W resolution.
Inefficiency is handled separately by response/unfolding.

Outputs:
  <out-prefix>_W_resolution_summary.csv
  <out-prefix>_W_resolution_vs_W.png
  <out-prefix>_W_relative_resolution_vs_W.png
  <out-prefix>_deltaW_distributions.png
  <out-prefix>_response_counts.png

Example:
  python3 study_W_resolution_swan.py /eos/user/g/gkrintir/merged_elecs_test_tuples_fix.root \
    --state 1S \
    --electron-energy 18 \
    --w-edges "37.1429,45.7143,54.2857,62.8571,71.4286,80,88.5714,97.1429,105.714,114.286" \
    --q2-cut 0.001 \
    --out-prefix Wres_1S_18x275_stableW
"""

import argparse
import os
import math
from pathlib import Path

import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass

UPSILON = {
    "1S": {"mass": 9.46030, "sigma": 0.053},
    "2S": {"mass": 10.02326, "sigma": 0.056},
    "3S": {"mass": 10.35520, "sigma": 0.059},
}


def parse_edges(s):
    return np.array([float(x) for x in s.split(",") if x.strip()], dtype=float)


def branch_exists(path, tree_name, branch):
    f = ROOT.TFile.Open(os.path.expanduser(path))
    if not f or f.IsZombie():
        return False
    t = f.Get(tree_name)
    ok = bool(t and t.GetBranch(branch))
    f.Close()
    return ok


def weighted_quantile(values, quantiles, weights=None):
    values = np.asarray(values, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)

    if weights is None:
        return np.quantile(values, quantiles)

    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[mask]
    weights = weights[mask]

    if len(values) == 0:
        return np.full_like(quantiles, np.nan, dtype=float)

    order = np.argsort(values)
    values = values[order]
    weights = weights[order]

    cdf = np.cumsum(weights)
    cdf /= cdf[-1]

    return np.interp(quantiles, cdf, values)


def hist2d_weighted(x, y, w, xbins, ybins):
    H, _, _ = np.histogram2d(x, y, bins=[xbins, ybins], weights=w)
    return H


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--state", choices=["1S", "2S", "3S"], default="1S")
    p.add_argument("--electron-energy", type=float, default=18.0)
    p.add_argument("--w-edges", required=True)
    p.add_argument("--q2-cut", type=float, default=None)
    p.add_argument("--mass-window-sigma", type=float, default=3.0)
    p.add_argument("--true-mass-window", type=float, default=0.2)
    p.add_argument("--sigma-epz-center", type=float, default=None)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    p.add_argument("--out-prefix", default="W_resolution")
    args = p.parse_args()

    infile = os.path.expanduser(args.input)
    edges = parse_edges(args.w_edges)
    props = UPSILON[args.state]

    if args.sigma_epz_center is None:
        args.sigma_epz_center = 2.0 * args.electron_energy
    if args.sigma_epz_halfwidth is None:
        args.sigma_epz_halfwidth = 0.2 * args.electron_energy

    has_true_q2 = branch_exists(infile, args.tree, "true_Q2")

    branches = [
        "weight", "true_mass", "true_W", "reco_mass", "reco_W", "sigma_epz"
    ]
    if has_true_q2:
        branches.append("true_Q2")

    rdf = ROOT.RDataFrame(args.tree, infile)
    arr = rdf.AsNumpy(columns=branches)
    df = pd.DataFrame({k: np.asarray(v, dtype=float) for k, v in arr.items()})
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    mass = props["mass"]
    truth_state = (
        (df["true_mass"] >= mass - args.true_mass_window) &
        (df["true_mass"] <= mass + args.true_mass_window)
    )
    if args.q2_cut is not None and has_true_q2:
        truth_state = truth_state & (df["true_Q2"] >= 0.0) & (df["true_Q2"] <= args.q2_cut)

    reco_mwin = args.mass_window_sigma * props["sigma"]
    selected = (
        truth_state &
        (df["true_W"] >= edges[0]) & (df["true_W"] <= edges[-1]) &
        (df["reco_mass"] >= mass - reco_mwin) &
        (df["reco_mass"] <= mass + reco_mwin) &
        (np.abs(df["sigma_epz"] - args.sigma_epz_center) <= args.sigma_epz_halfwidth) &
        (df["reco_W"] >= edges[0]) & (df["reco_W"] <= edges[-1])
    )

    sel = df[selected].copy()
    if len(sel) == 0:
        raise RuntimeError("No selected events for W resolution study")

    sel["delta_W"] = sel["reco_W"] - sel["true_W"]
    sel["rel_delta_W"] = sel["delta_W"] / sel["true_W"]

    rows = []
    for i in range(len(edges) - 1):
        low, high = edges[i], edges[i + 1]
        b = sel[(sel["true_W"] >= low) & (sel["true_W"] < high)].copy()
        if len(b) == 0:
            rows.append({
                "bin": i, "W_low": low, "W_high": high, "W_center": 0.5 * (low + high),
                "n_selected": 0, "sum_weight": 0.0,
                "median_delta_W": np.nan, "core_sigma_W": np.nan,
                "mean_delta_W": np.nan, "rms_delta_W": np.nan,
                "median_rel_delta_W": np.nan, "core_sigma_rel": np.nan,
                "bin_width": high - low, "bin_width_over_2core_sigma": np.nan,
            })
            continue

        w = b["weight"].to_numpy()
        d = b["delta_W"].to_numpy()
        rd = b["rel_delta_W"].to_numpy()

        q16, q50, q84 = weighted_quantile(d, [0.16, 0.50, 0.84], w)
        rq16, rq50, rq84 = weighted_quantile(rd, [0.16, 0.50, 0.84], w)

        core = 0.5 * (q84 - q16)
        rcore = 0.5 * (rq84 - rq16)
        width = high - low

        rows.append({
            "bin": i,
            "W_low": low,
            "W_high": high,
            "W_center": 0.5 * (low + high),
            "n_selected": len(b),
            "sum_weight": float(w.sum()),
            "median_delta_W": q50,
            "core_sigma_W": core,
            "q16_delta_W": q16,
            "q84_delta_W": q84,
            "mean_delta_W": float(np.average(d, weights=w)) if w.sum() > 0 else np.nan,
            "rms_delta_W": float(np.sqrt(np.average((d - np.average(d, weights=w))**2, weights=w))) if w.sum() > 0 else np.nan,
            "median_rel_delta_W": rq50,
            "core_sigma_rel": rcore,
            "q16_rel_delta_W": rq16,
            "q84_rel_delta_W": rq84,
            "bin_width": width,
            "bin_width_over_2core_sigma": width / (2.0 * core) if core > 0 else np.nan,
            "abs_bias_over_core": abs(q50) / core if core > 0 else np.nan,
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(f"{args.out_prefix}_W_resolution_summary.csv", index=False)

    # Response counts for the same selected events
    H = hist2d_weighted(sel["true_W"], sel["reco_W"], sel["weight"], edges, edges)
    pd.DataFrame(
        H,
        index=[f"{edges[i]:g}-{edges[i+1]:g}" for i in range(len(edges)-1)],
        columns=[f"{edges[i]:g}-{edges[i+1]:g}" for i in range(len(edges)-1)]
    ).to_csv(f"{args.out_prefix}_response_counts_selected.csv")

    # Plot resolution vs W
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(
        summary["W_center"], summary["median_delta_W"],
        yerr=summary["core_sigma_W"], fmt="o", color="royalblue", capsize=3,
        label=r"median $\Delta W$ with 16--84\% core"
    )
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel(r"$W_{\gamma p}^{\rm true}$ [GeV]")
    ax.set_ylabel(r"$W_{\rm reco}-W_{\rm true}$ [GeV]")
    ax.set_title(rf"{args.state} W resolution")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_W_resolution_vs_W.png", dpi=250)
    plt.close(fig)

    # Relative resolution
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(
        summary["W_center"], summary["median_rel_delta_W"],
        yerr=summary["core_sigma_rel"], fmt="o", color="purple", capsize=3,
        label=r"median relative residual with 16--84\% core"
    )
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel(r"$W_{\gamma p}^{\rm true}$ [GeV]")
    ax.set_ylabel(r"$(W_{\rm reco}-W_{\rm true})/W_{\rm true}$")
    ax.set_title(rf"{args.state} relative W resolution")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_W_relative_resolution_vs_W.png", dpi=250)
    plt.close(fig)

    # Delta distributions in bins
    n_bins = len(edges) - 1
    ncols = 3
    nrows = int(np.ceil(n_bins / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.6 * nrows), squeeze=False)
    for i in range(n_bins):
        ax = axes[i // ncols][i % ncols]
        low, high = edges[i], edges[i + 1]
        b = sel[(sel["true_W"] >= low) & (sel["true_W"] < high)]
        if len(b):
            ax.hist(b["delta_W"], bins=40, weights=b["weight"], histtype="step", color="black")
        ax.axvline(0, color="gray", linestyle="--")
        ax.set_title(f"{low:g} < Wtrue < {high:g}")
        ax.set_xlabel(r"$W_{\rm reco}-W_{\rm true}$ [GeV]")
        ax.set_ylabel("weighted events")
    for j in range(n_bins, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_deltaW_distributions.png", dpi=220)
    plt.close(fig)

    # Response heatmap
    fig, ax = plt.subplots(figsize=(7, 6))
    mesh = ax.pcolormesh(edges, edges, H)
    fig.colorbar(mesh, ax=ax, label="weighted selected events")
    ax.set_xlabel(r"$W_{\rm reco}$ [GeV]")
    ax.set_ylabel(r"$W_{\rm true}$ [GeV]")
    ax.set_title("Selected-event response counts")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_response_counts.png", dpi=250)
    plt.close(fig)

    print("=== W resolution summary ===")
    print(summary.to_string(index=False))
    print(f"\nWrote {args.out_prefix}_W_resolution_summary.csv")
    print(f"Wrote {args.out_prefix}_response_counts_selected.csv")
    print(f"Wrote {args.out_prefix}_W_resolution_vs_W.png")
    print(f"Wrote {args.out_prefix}_W_relative_resolution_vs_W.png")
    print(f"Wrote {args.out_prefix}_deltaW_distributions.png")
    print(f"Wrote {args.out_prefix}_response_counts.png")


if __name__ == "__main__":
    main()
