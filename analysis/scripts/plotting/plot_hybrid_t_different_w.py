#!/usr/bin/env python3
"""
Build a hybrid |t| distribution using the low-energy sample for selected low-|t| bins
and the high-energy sample for the remaining bins.

Inputs are CSVs produced by unfold_t_distribution_swan.py:
  <prefix>_t_unfold_summary.csv

Default behavior:
  - take bins with t_high <= --switch-t from the low-energy CSV
  - take bins with t_low >= --switch-t from the high-energy CSV
  - if a bin crosses the boundary, keep the lower-energy version only when its center < switch-t

This is meant for a detector-feasibility / coverage plot. If the two samples cover
different W ranges, label the result as hybrid-energy and do not present it as a
single physical dN/d|t| spectrum without a common W reweighting.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


def read_csv(path):
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def load_zeus(path):
    if not path:
        return None
    return read_csv(path)


def zeus_centers(zeus):
    if {"t_low", "t_high"}.issubset(zeus.columns):
        return 0.5 * (zeus["t_low"].to_numpy() + zeus["t_high"].to_numpy())
    if "t_center" in zeus.columns:
        return zeus["t_center"].to_numpy()
    return zeus["t_GeV2"].to_numpy()

def zeus_xerr(zeus):
    if {"t_low", "t_high"}.issubset(zeus.columns):
        return 0.5 * (zeus["t_high"].to_numpy() - zeus["t_low"].to_numpy())
    return None


def scale_model_to_zeus(x, y, zeus):
    scale_num = 0.0
    scale_den = 0.0
    zx = zeus_centers(zeus)
    for idx_row, r in zeus.iterrows():
        idx = int(np.argmin(np.abs(x - zx[idx_row])))
        m = y[idx]
        if not np.isfinite(m) or m <= 0:
            continue
        err = 0.0
        if "err_low" in zeus.columns and "err_high" in zeus.columns:
            err = 0.5 * (float(r["err_low"]) + float(r["err_high"]))
        if err <= 0:
            err = max(float(r["dNdt_GeVminus2"]) * 0.2, 1.0)
        w = 1.0 / (err * err)
        scale_num += float(r["dNdt_GeVminus2"]) * m * w
        scale_den += m * m * w
    return scale_num / scale_den if scale_den > 0 else 1.0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--low-csv", required=True, help="Low-energy t unfolding summary CSV")
    p.add_argument("--high-csv", required=True, help="High-energy t unfolding summary CSV")
    p.add_argument("--switch-t", type=float, default=0.2, help="Use low-energy bins below this |t| edge")
    p.add_argument("--low-label", default=r"10$\times$100 GeV")
    p.add_argument("--high-label", default=r"18$\times$275 GeV")
    p.add_argument("--zeus-csv", default="")
    p.add_argument("--zeus-mode", choices=["scale_to_zeus", "none"], default="scale_to_zeus")
    p.add_argument("--out-prefix", default="hybrid_t_low_high")
    p.add_argument("--y-min", type=float, default=0.1)
    p.add_argument("--y-max", type=float, default=1000.0)
    args = p.parse_args()

    low = read_csv(args.low_csv).copy()
    high = read_csv(args.high_csv).copy()

    low["sample"] = "low"
    low["sample_label"] = args.low_label
    high["sample"] = "high"
    high["sample_label"] = args.high_label

    # Pick low-energy for bins below switch-t, high-energy for bins above.
    low_sel = low[low["t_high"] <= args.switch_t].copy()
    high_sel = high[high["t_low"] >= args.switch_t].copy()

    if low_sel.empty:
        # Fallback: take the first low bin.
        low_sel = low.sort_values("t_center").head(1).copy()
        print("WARNING: no low-energy bin has t_high <= switch-t; using the first low-energy bin.")

    hybrid = pd.concat([low_sel, high_sel], ignore_index=True).sort_values("t_center")
    hybrid.to_csv(args.out_prefix + "_hybrid.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.5, 6.8))

    # Plot low and high source points faintly for diagnostics.
    for df, marker, color, label in [
        (low, "o", "tab:blue", args.low_label + " all bins"),
        (high, "s", "tab:orange", args.high_label + " all bins"),
    ]:
        x = df["t_center"]
        xerr = 0.5 * (df["t_high"] - df["t_low"])
        y = df["dNdt_unfolded_visible"]
        yerr = df["dNdt_unfolded_visible_err"]
        m = np.isfinite(y) & (y > 0)
        ax.errorbar(x[m], y[m], xerr=xerr[m], yerr=yerr[m], fmt=marker, color=color, alpha=0.25, capsize=2, label=label)

    # Hybrid points.
    x = hybrid["t_center"]
    xerr = 0.5 * (hybrid["t_high"] - hybrid["t_low"])
    y = hybrid["dNdt_unfolded_visible"]
    yerr = hybrid["dNdt_unfolded_visible_err"]

    zeus = load_zeus(args.zeus_csv)
    if zeus is not None and args.zeus_mode == "scale_to_zeus":
        scale = scale_model_to_zeus(x.to_numpy(), y.to_numpy(), zeus)
        y_plot = y * scale
        yerr_plot = yerr * scale
        ylabel = r"scaled $dN/d|t|$ [arb.]"
        hybrid_label = rf"hybrid ePIC shape $\times {scale:.2g}$"
    else:
        y_plot = y
        yerr_plot = yerr
        ylabel = r"$dN_{ee}/d|t|$ [GeV$^{-2}$]"
        hybrid_label = "hybrid ePIC"

    m = np.isfinite(y_plot) & (y_plot > 0)
    ax.errorbar(x[m], y_plot[m], xerr=xerr[m], yerr=yerr_plot[m], fmt="D-", color="black", capsize=3, label=hybrid_label, zorder=5)

    if zeus is not None:
        zx = zeus_centers(zeus)
        zxerr = zeus_xerr(zeus)
        zy = zeus["dNdt_GeVminus2"].to_numpy()
        zlo = zeus["err_low"].to_numpy() if "err_low" in zeus.columns else np.zeros(len(zeus))
        zhi = zeus["err_high"].to_numpy() if "err_high" in zeus.columns else np.zeros(len(zeus))
        if {"t_low", "t_high"}.issubset(zeus.columns):
            step_edges = np.r_[zeus["t_low"].to_numpy(), zeus["t_high"].to_numpy()[-1]]
            step_vals = np.r_[zy, zy[-1]]
            ax.step(step_edges, step_vals, where="post", color="black", linestyle=":", linewidth=2.2,
                    label=r"ZEUS exclusive $\Upsilon(1S)$ component", zorder=4)
            ax.errorbar(zx, zy, xerr=zxerr, yerr=np.vstack([zlo, zhi]),
                        fmt="none", ecolor="black", capsize=2, zorder=4)
        else:
            ax.errorbar(zx, zy, yerr=np.vstack([zlo, zhi]),
                        fmt="ko", markerfacecolor="white", capsize=3,
                        label=r"ZEUS exclusive $\Upsilon(1S)$ component", zorder=4)

    for _, r in hybrid.iterrows():
        ax.text(r["t_center"], y_plot.loc[r.name] * 1.15 if r.name in y_plot.index else r["dNdt_unfolded_visible"] * 1.15,
                r["sample_label"], fontsize=8, ha="center", va="bottom", rotation=20)

    ax.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.set_xlim(min(low["t_low"].min(), high["t_low"].min()), max(low["t_high"].max(), high["t_high"].max()))
    ax.set_ylim(args.y_min, args.y_max)
    ax.set_title(r"Hybrid low/high-energy $|t|$ unfolding comparison")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()

    fig.savefig(args.out_prefix + "_hybrid_overlay.png", dpi=250)
    fig.savefig(args.out_prefix + "_hybrid_overlay.pdf")

    print("=== Hybrid bins ===")
    print(hybrid[[
        "sample_label", "t_low", "t_high", "t_center",
        "efficiency", "purity", "stability",
        "closure_unfold_over_truth",
        "dNdt_unfolded_visible", "dNdt_unfolded_visible_err"
    ]].to_string(index=False))
    print(f"\nWrote {args.out_prefix}_hybrid.csv")
    print(f"Wrote {args.out_prefix}_hybrid_overlay.png")
    print(f"Wrote {args.out_prefix}_hybrid_overlay.pdf")


if __name__ == "__main__":
    main()
