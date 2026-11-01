#!/usr/bin/env python3
"""
Make a compact set of diagnostic/comparison plots from the final Upsilon sigma(gamma p)
CSV outputs.

Inputs must be supplied explicitly. This avoids silently selecting stale files
whose filenames and recorded phase space disagree.
  - flux_ratio_scan_W20to26p5_Q2lt0p001.csv
  - compare_10x100_18x275_truth_W20to26p5_Q2lt0p001_stateNorm.csv
  - sigma_gp_family_18x275_fixed_allbins.csv   optional
  - overlay_sigma_gp_1S_10x100_18x275.csv      optional

Outputs:
  - upsilon_selected_comparisons.pdf
  - upsilon_selected_overlay_1S.png
  - upsilon_flux_validation.png
  - upsilon_truth_low_high_ratio.png
  - upsilon_18x275_family_efficiency.png

Run:
  python3 plot_selected_upsilon_comparisons.py --indir . --out-prefix upsilon_selected
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


FINAL_LOW_1S = ""
FINAL_HIGH_FAMILY = ""
FLUX_VALIDATION = "flux_ratio_scan_W20to26p5_Q2lt0p001.csv"
TRUTH_COMPARISON = "compare_10x100_18x275_truth_W20to26p5_Q2lt0p001_stateNorm.csv"
HIGH_ALLBINS = "sigma_gp_family_18x275_fixed_allbins.csv"
OVERLAY_1S = "overlay_sigma_gp_1S_10x100_18x275.csv"


STATE_STYLE = {
    "1S": {"color": "royalblue", "marker": "o", "label": r"$\Upsilon(1S)$"},
    "2S": {"color": "crimson", "marker": "s", "label": r"$\Upsilon(2S)$"},
    "3S": {"color": "forestgreen", "marker": "^", "label": r"$\Upsilon(3S)$"},
}


def read_csv(indir, name, required=True):
    path = Path(indir) / name
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        print(f"[skip] {path} not found")
        return None
    return pd.read_csv(path)


def finite_positive(df, col):
    return np.isfinite(df[col]) & (df[col] > 0)


def add_epic_label(ax, text_lines=None):
    lines = ["ePIC Performance", r"exclusive $\Upsilon$ photoproduction"]
    if text_lines:
        lines.extend(text_lines)
    y = 0.94
    for i, line in enumerate(lines):
        ax.text(
            0.04, y - i * 0.06, line,
            transform=ax.transAxes,
            fontsize=12 if i == 0 else 11,
            fontweight="bold" if i == 0 else None,
            ha="left", va="top",
        )


def plot_final_overlay(low, high, min_eff, out_png=None):
    fig, ax = plt.subplots(figsize=(9.5, 7.5))

    # HERA exclusive reference
    hera_w = np.array([95.0, 175.0])
    hera_w_err = np.array([35.0, 45.0])
    hera_sig = np.array([160.0, 321.0])
    hera_err = np.array([70.0, 144.0])
    ax.errorbar(
        hera_w, hera_sig, xerr=hera_w_err, yerr=hera_err,
        fmt="^", color="gray", markersize=8, linewidth=1.6,
        label=r"HERA exclusive $\Upsilon(1S)$",
        zorder=1,
    )

    # 10x100 1S
    ldf = low.copy()
    ldf = ldf[(ldf["state"] == "1S") & finite_positive(ldf, "sigma_gamma_p_pb")]
    ldf = ldf[ldf["eff_reco_over_true"] > min_eff]
    if len(ldf):
        ax.errorbar(
            ldf["W_center"], ldf["sigma_gamma_p_pb"], yerr=ldf["sigma_gamma_p_err_pb"],
            fmt="o", markerfacecolor="white", markeredgecolor="royalblue",
            markeredgewidth=1.8, color="royalblue", capsize=3,
            label=r"10$\times$100 GeV $\Upsilon(1S)$",
            zorder=3,
        )

    # 18x275 family
    hdf = high.copy()
    hdf = hdf[finite_positive(hdf, "sigma_gamma_p_pb")]
    hdf = hdf[hdf["eff_reco_over_true"] > min_eff]
    for state, sdf in hdf.groupby("state"):
        st = STATE_STYLE.get(state, {"color": "black", "marker": "o", "label": state})
        ax.errorbar(
            sdf["W_center"], sdf["sigma_gamma_p_pb"], yerr=sdf["sigma_gamma_p_err_pb"],
            fmt=st["marker"], color=st["color"], capsize=3,
            label=rf"18$\times$275 GeV {st['label']}",
            zorder=4,
        )

    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    ax.set_ylabel(r"$\sigma(\gamma p\rightarrow \Upsilon(nS)p)$ [pb]")
    ax.set_yscale("log")
    ax.set_xlim(18, 225)
    ax.set_ylim(0.5, 2000)
    add_epic_label(ax, [r"$Q^2<10^{-3}$ GeV$^2$", rf"$\epsilon_{{reco}}>{min_eff:g}$"])
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=300)
    return fig


def plot_flux_validation(flux, out_png=None):
    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    x = flux["W_center"] if "W_center" in flux.columns else 0.5 * (flux["W_low"] + flux["W_high"])
    y = flux["internal_over_analytic_ratio"]

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.5)
    ax.plot(x, y, "o-", color="black", label="eSTARlight ratio / analytic WW ratio")

    for _, r in flux.iterrows():
        xc = 0.5 * (r["W_low"] + r["W_high"])
        ax.text(xc, r["internal_over_analytic_ratio"] + 0.006, f"{r['internal_over_analytic_ratio']:.3f}",
                ha="center", va="bottom", fontsize=9)

    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    ax.set_ylabel(r"$R_{\rm eSTARlight}/R_{\rm WW}$")
    ax.set_ylim(0.94, 1.06)
    ax.set_title("Photon-flux validation, 10×100 over 18×275")
    add_epic_label(ax, [r"$Q^2_{\max}=10^{-3}$ GeV$^2$"])
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=300)
    return fig


def plot_truth_ratio(comp, out_png=None):
    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    x = 0.5 * (comp["W_low"] + comp["W_high"])
    ratio = comp["ratio_low_over_high"]

    ax.axhline(1.0, color="gray", linestyle="--")
    ax.plot(x, ratio, "o-", color="purple", label=r"tuple truth-level 10$\times$100 / 18$\times$275")

    for _, r in comp.iterrows():
        xc = 0.5 * (r["W_low"] + r["W_high"])
        ax.text(xc, r["ratio_low_over_high"] + 0.08, f"{r['ratio_low_over_high']:.2f}",
                ha="center", va="bottom", fontsize=9)

    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    ax.set_ylabel(r"$\sigma_{\gamma p}^{10\times100}/\sigma_{\gamma p}^{18\times275}$")
    ax.set_ylim(0, max(3.0, ratio.max() * 1.25))
    ax.set_title("Truth-level low/high comparison after normalization fixes")
    add_epic_label(ax, [r"$Q^2<10^{-3}$ GeV$^2$", "state-scoped normalization"])
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=300)
    return fig


def plot_efficiency_vs_w(high, out_png=None):
    fig, ax = plt.subplots(figsize=(9, 6.5))

    for state, sdf in high.groupby("state"):
        st = STATE_STYLE.get(state, {"color": "black", "marker": "o", "label": state})
        sdf = sdf.sort_values("W_center")
        ax.plot(
            sdf["W_center"], sdf["eff_reco_over_true"],
            marker=st["marker"], color=st["color"], label=st["label"]
        )

    ax.axhline(0.05, color="gray", linestyle="--", label=r"$\epsilon=0.05$")
    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    ax.set_ylabel(r"$\epsilon_{\rm reco}=N_{\rm reco}/N_{\rm truth}$")
    ax.set_ylim(0, 0.8)
    ax.set_title(r"18$\times$275 reconstruction efficiency by state")
    add_epic_label(ax, [r"$Q^2<10^{-3}$ GeV$^2$"])
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=300)
    return fig



def plot_closure(low, high, min_eff, out_png=None):
    fig, ax = plt.subplots(figsize=(9, 6.5))

    frames = []
    if "closure_unfold_over_truth" in low.columns:
        tmp = low.copy()
        tmp["sample_label"] = r"10$\times$100 GeV"
        tmp["energy_label"] = "low"
        frames.append(tmp)
    if "closure_unfold_over_truth" in high.columns:
        tmp = high.copy()
        tmp["sample_label"] = r"18$\times$275 GeV"
        tmp["energy_label"] = "high"
        frames.append(tmp)

    if not frames:
        ax.text(0.5, 0.5, "No closure_unfold_over_truth column found", ha="center", va="center")
        return fig

    df = pd.concat(frames, ignore_index=True)
    df = df[np.isfinite(df["closure_unfold_over_truth"])]
    df = df[df["eff_reco_over_true"] > min_eff]

    markers = {"low": "o", "high": "s"}

    for (state, energy), sdf in df.groupby(["state", "energy_label"]):
        st = STATE_STYLE.get(state, {"color": "black", "marker": "o", "label": state})
        label_state = st["label"]
        label_energy = sdf["sample_label"].iloc[0]
        marker = markers.get(energy, "o")
        face = "white" if energy == "low" else st["color"]

        ax.plot(
            sdf["W_center"],
            sdf["closure_unfold_over_truth"],
            marker=marker,
            linestyle="-",
            color=st["color"],
            markerfacecolor=face,
            markeredgecolor=st["color"],
            markeredgewidth=1.6,
            label=f"{label_energy} {label_state}",
        )

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.5)
    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    ax.set_ylabel(r"$\sigma_{\rm unfolded}/\sigma_{\rm truth}$")
    ax.set_ylim(0.85, 1.15)
    ax.set_title("Bin-by-bin unfolding closure")
    add_epic_label(ax, [r"$Q^2<10^{-3}$ GeV$^2$", rf"$\epsilon_{{reco}}>{min_eff:g}$"])
    ax.legend(frameon=False, loc="best", fontsize=9)
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=300)
    return fig


def print_summary(low, high, flux, comp, min_eff, low_name, high_name, flux_name, comparison_name):
    print("\n=== Recommended final inputs ===")
    print(f"10x100 final 1S: {low_name}")
    print(f"18x275 final family: {high_name}")
    print(f"Flux validation: {flux_name}")
    print(f"Truth comparison: {comparison_name}")

    print("\n=== 10x100 plotted bins ===")
    l = low[(low["state"] == "1S") & finite_positive(low, "sigma_gamma_p_pb") & (low["eff_reco_over_true"] > min_eff)]
    print(l[["W_low", "W_high", "sigma_gamma_p_pb", "sigma_gamma_p_err_pb", "eff_reco_over_true"]].to_string(index=False))

    print("\n=== 18x275 plotted bins by state ===")
    h = high[finite_positive(high, "sigma_gamma_p_pb") & (high["eff_reco_over_true"] > min_eff)]
    print(h[["state", "W_low", "W_high", "sigma_gamma_p_pb", "sigma_gamma_p_err_pb", "eff_reco_over_true"]].to_string(index=False))

    print("\n=== Flux closure ===")
    print(flux[["W_low", "W_high", "internal_over_analytic_ratio"]].to_string(index=False))

    print("\n=== Truth low/high ratio ===")
    print(comp[["W_low", "W_high", "ratio_low_over_high"]].to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=".")
    ap.add_argument("--out-prefix", default="upsilon_selected")
    ap.add_argument("--min-eff", type=float, default=0.05)
    ap.add_argument("--low-csv", required=True)
    ap.add_argument("--high-csv", required=True)
    ap.add_argument("--flux-csv", default=FLUX_VALIDATION)
    ap.add_argument("--truth-compare-csv", default=TRUTH_COMPARISON)
    args = ap.parse_args()

    indir = Path(args.indir)

    low = read_csv(indir, args.low_csv)
    high = read_csv(indir, args.high_csv)
    flux = read_csv(indir, args.flux_csv)
    comp = read_csv(indir, args.truth_compare_csv)

    print_summary(
        low, high, flux, comp, args.min_eff,
        args.low_csv, args.high_csv, args.flux_csv, args.truth_compare_csv,
    )

    pdf_path = f"{args.out_prefix}_comparisons.pdf"

    figs = []
    figs.append(plot_final_overlay(low, high, args.min_eff, f"{args.out_prefix}_overlay_1S_family.png"))
    figs.append(plot_flux_validation(flux, f"{args.out_prefix}_flux_validation.png"))
    figs.append(plot_truth_ratio(comp, f"{args.out_prefix}_truth_low_high_ratio.png"))
    figs.append(plot_efficiency_vs_w(high, f"{args.out_prefix}_18x275_family_efficiency.png"))
    figs.append(plot_closure(low, high, args.min_eff, f"{args.out_prefix}_closure.png"))

    with PdfPages(pdf_path) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)

    print(f"\nSaved {pdf_path}")
    print(f"Saved {args.out_prefix}_overlay_1S_family.png")
    print(f"Saved {args.out_prefix}_flux_validation.png")
    print(f"Saved {args.out_prefix}_truth_low_high_ratio.png")
    print(f"Saved {args.out_prefix}_18x275_family_efficiency.png")
    print(f"Saved {args.out_prefix}_closure.png")


if __name__ == "__main__":
    main()
