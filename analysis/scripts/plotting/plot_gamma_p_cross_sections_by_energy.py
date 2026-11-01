#!/usr/bin/env python3
"""
Overlay two ePIC sigma(gamma p -> Upsilon(nS) p) extractions.

Input CSVs should be produced by
  extract_sigma_gamma_p_family_fixed_swan_v2.py

The script expects columns such as
  state
  W_center
  W_low
  W_high
  sigma_gamma_p_pb
  sigma_gamma_p_err_pb
  eff_reco_over_true
  plotted

It can overlay 10x100 and 18x275 outputs, with optional HERA reference points.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


STATE_STYLE = {
    "1S": {"color": "royalblue", "label": r"$\Upsilon(1S)$"},
    "2S": {"color": "crimson", "label": r"$\Upsilon(2S)$"},
    "3S": {"color": "forestgreen", "label": r"$\Upsilon(3S)$"},
}


def parse_states(s):
    return [x.strip() for x in s.split(",") if x.strip()]


def load_extraction_csv(path, sample_label, energy_label, min_eff, states, require_plotted):
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    needed = {
        "state",
        "W_center",
        "W_low",
        "W_high",
        "sigma_gamma_p_pb",
        "sigma_gamma_p_err_pb",
        "eff_reco_over_true",
    }
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns {sorted(missing)}")

    df = df.copy()
    df["sample_label"] = sample_label
    df["energy_label"] = energy_label
    df["source_csv"] = path

    if states:
        df = df[df["state"].isin(states)]

    df = df[np.isfinite(df["sigma_gamma_p_pb"])]
    df = df[df["sigma_gamma_p_pb"] > 0]
    df = df[df["eff_reco_over_true"] > min_eff]

    if require_plotted and "plotted" in df.columns:
        # Accept both boolean and string representations.
        plotted = df["plotted"]
        if plotted.dtype == object:
            plotted = plotted.astype(str).str.lower().isin(["true", "1", "yes"])
        df = df[plotted]

    return df


def add_hera(ax):
    hera_w = np.array([95.0, 175.0])
    hera_w_err = np.array([35.0, 45.0])
    hera_sig = np.array([160.0, 321.0])
    hera_err = np.array([70.0, 144.0])

    ax.errorbar(
        hera_w,
        hera_sig,
        xerr=hera_w_err,
        yerr=hera_err,
        fmt="^",
        color="gray",
        label=r"HERA exclusive $\gamma p \rightarrow \Upsilon(1S)p$",
        markersize=8,
        linewidth=1.8,
        capsize=0,
        zorder=1,
    )


def apply_x_offsets(df, offset_fraction):
    if offset_fraction == 0:
        df["W_plot"] = df["W_center"]
        return df

    df = df.copy()
    widths = df["W_high"] - df["W_low"]
    # Shift low energy slightly left and high energy slightly right.
    direction = np.where(df["sample_label"].astype(str).str.contains("10", case=False, regex=False), -1.0, 1.0)
    df["W_plot"] = df["W_center"] + direction * offset_fraction * widths
    return df


def main():
    parser = argparse.ArgumentParser(description="Overlay two ePIC Upsilon gamma-p cross section extractions.")
    parser.add_argument("--low-csv", required=True, help="CSV from 10x100 extraction")
    parser.add_argument("--high-csv", required=True, help="CSV from 18x275 extraction")
    parser.add_argument("--low-label", default=r"ePIC 10$\times$100 GeV")
    parser.add_argument("--high-label", default=r"ePIC 18$\times$275 GeV")
    parser.add_argument("--states", default="1S", help="Comma-separated list such as 1S or 1S,2S,3S")
    parser.add_argument("--min-eff", type=float, default=0.05)
    parser.add_argument("--require-plotted", action="store_true", help="Only use rows marked plotted=True in the extraction CSV")
    parser.add_argument("--offset-fraction", type=float, default=0.06, help="Horizontal offset as fraction of W-bin width")
    parser.add_argument("--out-prefix", default="sigma_gp_overlay_two_energies")
    parser.add_argument("--title", default="")
    parser.add_argument("--q2-label", default=r"$Q^2 < 0.001$ GeV$^2$")
    parser.add_argument("--lumi-label", default=r"$L_{ep}=10$ fb$^{-1}$")
    parser.add_argument("--show-hera", action="store_true", default=True)
    parser.add_argument("--no-hera", action="store_false", dest="show_hera")
    parser.add_argument("--x-min", type=float, default=None)
    parser.add_argument("--x-max", type=float, default=None)
    parser.add_argument("--y-min", type=float, default=0.5)
    parser.add_argument("--y-max", type=float, default=2000.0)
    args = parser.parse_args()

    states = parse_states(args.states)

    low = load_extraction_csv(
        args.low_csv,
        sample_label=args.low_label,
        energy_label="low",
        min_eff=args.min_eff,
        states=states,
        require_plotted=args.require_plotted,
    )
    high = load_extraction_csv(
        args.high_csv,
        sample_label=args.high_label,
        energy_label="high",
        min_eff=args.min_eff,
        states=states,
        require_plotted=args.require_plotted,
    )

    combined = pd.concat([low, high], ignore_index=True)
    combined = apply_x_offsets(combined, args.offset_fraction)

    out_csv = os.path.expanduser(args.out_prefix + ".csv")
    out_png = os.path.expanduser(args.out_prefix + ".png")
    out_pdf = os.path.expanduser(args.out_prefix + ".pdf")
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)

    combined.to_csv(out_csv, index=False)

    fig, ax = plt.subplots(figsize=(9.5, 7.5))

    if args.show_hera:
        add_hera(ax)

    markers = {
        "low": "o",
        "high": "s",
    }

    for state in states:
        sdf = combined[combined["state"] == state]
        if sdf.empty:
            continue

        style = STATE_STYLE.get(state, {"color": "black", "label": state})

        for energy_label, edf in sdf.groupby("energy_label"):
            edf = edf.sort_values("W_center")
            sample_label = edf["sample_label"].iloc[0]
            marker = markers.get(energy_label, "o")

            ax.errorbar(
                edf["W_plot"],
                edf["sigma_gamma_p_pb"],
                yerr=edf["sigma_gamma_p_err_pb"],
                fmt=marker,
                color=style["color"],
                markerfacecolor=style["color"] if energy_label == "high" else "white",
                markeredgecolor=style["color"],
                markeredgewidth=1.8,
                label=f"{sample_label} {style['label']}",
                markersize=7.5,
                linewidth=1.6,
                elinewidth=1.3,
                capsize=3,
                zorder=3,
            )

    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]", fontsize=14)
    ax.set_ylabel(r"$\sigma(\gamma p \rightarrow \Upsilon(nS)p)$ [pb]", fontsize=14)
    ax.set_yscale("log")

    if args.x_min is None:
        x_min = max(0.1, combined["W_low"].min() * 0.92) if len(combined) else 10
    else:
        x_min = args.x_min

    if args.x_max is None:
        x_max = combined["W_high"].max() * 1.05 if len(combined) else 150
        if args.show_hera:
            x_max = max(x_max, 225)
    else:
        x_max = args.x_max

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(args.y_min, args.y_max)

    ax.text(0.04, 0.94, "ePIC Performance", transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="top")
    ax.text(0.04, 0.88, r"$e+p$ exclusive $\Upsilon$ photoproduction", transform=ax.transAxes, fontsize=11, ha="left", va="top")
    ax.text(0.04, 0.82, args.q2_label, transform=ax.transAxes, fontsize=11, ha="left", va="top")
    ax.text(0.04, 0.76, args.lumi_label, transform=ax.transAxes, fontsize=11, ha="left", va="top")

    if args.title:
        ax.set_title(args.title, fontsize=14)

    ax.legend(loc="lower right", fontsize=9.5, frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)

    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")
    print(f"Saved {out_csv}")
    print("")
    print("Rows used in overlay")
    print(combined[[
        "sample_label",
        "state",
        "W_low",
        "W_high",
        "W_center",
        "eff_reco_over_true",
        "sigma_gamma_p_pb",
        "sigma_gamma_p_err_pb",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
