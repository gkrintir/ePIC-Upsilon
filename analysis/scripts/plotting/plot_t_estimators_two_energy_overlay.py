#!/usr/bin/env python3
"""
Overlay reconstructed |t| estimators for Υ(1S), comparing 10x100 and 18x275
in each panel with the correct energy-specific selections.
"""

import argparse
import os

import ROOT
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

hep.style.use(hep.style.ROOT)

TREE_NAME = "UpsilonTree"
EXPECTED_LUMI_INVFB = 10.0

SAMPLES = {
    "low": {
        "file": "~/merged_elecs_test_10t100_10x100.root",
        "label": "10×100 GeV",
        "campaign": "26.04.1 Campaign",
        "electron_energy": 10.0,
        "visible_xsec_fb": 104.793,
        "color": "dodgerblue",
        "w_min": 20.0,
        "w_max": 26.5,
    },
    "high": {
        "file": "~/merged_elecs_test_tuples_fix2.root",
        "label": "18×275 GeV",
        "campaign": "26.03.1 Campaign",
        "electron_energy": 18.0,
        "visible_xsec_fb": 327.584,
        "color": "darkorange",
        "w_min": 37.0,
        "w_max": 114.0,
    },
}

STATE_NORM_MODE = "target_visible_xsec"   # "target_visible_xsec" or "raw"
APPLY_SIGMA_EPZ_CUT = True
SIGMA_EPZ_HALF_WIDTH_FRACTION = 0.20

Y1S_MASS = 9.46030
Y1S_TRUTH_WINDOW = 0.20

T_BINS = np.linspace(0.0, 5.0, 51)

METHODS = [
    {
        "key": "eX",
        "branch": "reco_t_eX",
        "weight": "weight_visible",
        "label": r"$t_{eX}$",
        "xlabel": r"Reco $|t_{eX}|$ [GeV$^2$]",
    },
    {
        "key": "BABE",
        "branch": "reco_t_BABE",
        "weight": "weight_visible",
        "label": r"$t_{\rm BABE}$",
        "xlabel": r"Reco $|t_{\rm BABE}|$ [GeV$^2$]",
    },
    {
        "key": "eXBABE",
        "branch": "reco_t_eXBABE",
        "weight": "weight_visible",
        "label": r"$t_{eX{\rm BABE}}$",
        "xlabel": r"Reco $|t_{eX{\rm BABE}}|$ [GeV$^2$]",
    },
    {
        "key": "proj",
        "branch": "reco_t_proj",
        "weight": "weight_proj_visible",
        "label": r"$t_{\rm proj}$",
        "xlabel": r"Projective $|t|$ [GeV$^2$]",
    },
    {
        "key": "pT2",
        "branch": "reco_pT2",
        "weight": "weight_visible",
        "label": r"$p_{T,\Upsilon}^{2}$",
        "xlabel": r"Reco $p_{T,\Upsilon}^{2}$ [GeV$^2$]",
    },
]


def make_hist(data, bins, weights):
    counts, _ = np.histogram(data, bins=bins, weights=weights)
    variances, _ = np.histogram(data, bins=bins, weights=weights * weights)
    return counts, np.sqrt(variances)


def truth_1s_expr():
    return f"true_mass >= {Y1S_MASS - Y1S_TRUTH_WINDOW:.12g} && true_mass <= {Y1S_MASS + Y1S_TRUTH_WINDOW:.12g}"


def sigma_epz_info(sample_key):
    cfg = SAMPLES[sample_key]
    center = 2.0 * cfg["electron_energy"]
    half_width = SIGMA_EPZ_HALF_WIDTH_FRACTION * cfg["electron_energy"]
    return center, half_width


def sigma_epz_cut_expr(sample_key):
    center, half_width = sigma_epz_info(sample_key)
    return f"abs(sigma_epz - {center:.12g}) <= {half_width:.12g}"


def add_visible_weight_columns(df, sample_key):
    cfg = SAMPLES[sample_key]
    scale = 1.0

    if STATE_NORM_MODE == "target_visible_xsec":
        sumw = df.Filter(truth_1s_expr()).Sum("weight").GetValue()
        target_sumw = cfg["visible_xsec_fb"] * EXPECTED_LUMI_INVFB
        scale = target_sumw / sumw if sumw > 0 else 1.0
        raw_xsec = sumw / EXPECTED_LUMI_INVFB if EXPECTED_LUMI_INVFB > 0 else 0.0
        print(f"{cfg['label']} {cfg['campaign']}: raw visible 1S xsec = {raw_xsec:.3f} fb, "
              f"target = {cfg['visible_xsec_fb']:.3f} fb, scale = {scale:.6g}")
    elif STATE_NORM_MODE != "raw":
        raise ValueError(f"Unknown STATE_NORM_MODE: {STATE_NORM_MODE}")

    return (
        df.Define(f"state_norm_{sample_key}", f"{scale}")
          .Define("weight_visible", f"weight * state_norm_{sample_key}")
          .Define("weight_proj_visible", f"weight_proj * state_norm_{sample_key}")
          .Define("reco_pT2", "reco_pT * reco_pT")
    )


def base_selection(sample_key):
    cfg = SAMPLES[sample_key]
    cuts = [
        "reco_mass >= 8.0 && reco_mass <= 12.0",
        f"reco_W >= {cfg['w_min']:.12g} && reco_W <= {cfg['w_max']:.12g}",
        f"({truth_1s_expr()})",
    ]
    if APPLY_SIGMA_EPZ_CUT:
        cuts.append(sigma_epz_cut_expr(sample_key))
    return " && ".join(cuts)


def load_sample_arrays(sample_key, branches):
    cfg = SAMPLES[sample_key]
    df = ROOT.RDataFrame(TREE_NAME, cfg["file"])
    df = add_visible_weight_columns(df, sample_key)
    df_sel = df.Filter(base_selection(sample_key))
    return df_sel.AsNumpy(branches)


def draw_common_labels(ax):
    ax.text(0.04, 0.94, "ePIC Performance", transform=ax.transAxes,
            fontsize=11, fontweight="bold", ha="left", va="top")
    ax.text(0.04, 0.88, r"$e+p\ \Upsilon(1S)\to e^+e^-$", transform=ax.transAxes,
            fontsize=10, ha="left", va="top")
    ax.text(0.04, 0.82, r"Visible $e^+e^-$ yields, $L=10~\mathrm{fb}^{-1}$", transform=ax.transAxes,
            fontsize=10, ha="left", va="top")

    cL, hwL = sigma_epz_info("low")
    cH, hwH = sigma_epz_info("high")

    lines = [
        rf"10×100, {SAMPLES['low']['campaign']}, $W_{{\gamma p}}\in[{SAMPLES['low']['w_min']:g},{SAMPLES['low']['w_max']:g}]$ GeV",
        rf"$|\Sigma(E-p_z)-{cL:g}|<{hwL:g}$ GeV",
        rf"18×275, {SAMPLES['high']['campaign']}, $W_{{\gamma p}}\in[{SAMPLES['high']['w_min']:g},{SAMPLES['high']['w_max']:g}]$ GeV",
        rf"$|\Sigma(E-p_z)-{cH:g}|<{hwH:g}$ GeV",
    ]
    y = 0.76
    for line in lines:
        ax.text(0.04, y, line, transform=ax.transAxes, fontsize=9.5, ha="left", va="top")
        y -= 0.06


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay reconstructed Upsilon(1S) |t| estimators at two beam energies."
    )
    parser.add_argument("--low-file", required=True)
    parser.add_argument("--high-file", required=True)
    parser.add_argument("--tree", default=TREE_NAME)
    parser.add_argument("--low-visible-xsec-fb", type=float, default=104.793)
    parser.add_argument("--high-visible-xsec-fb", type=float, default=327.584)
    parser.add_argument("--lumi-invfb", type=float, default=10.0)
    parser.add_argument("--state-norm-mode", choices=["raw", "target_visible_xsec"], default="target_visible_xsec")
    parser.add_argument("--no-sigma-epz", action="store_true")
    parser.add_argument("--out-prefix", default="upsilon_1s_t_estimators_two_energy_overlay")
    return parser.parse_args()


def main():
    global TREE_NAME, EXPECTED_LUMI_INVFB, STATE_NORM_MODE, APPLY_SIGMA_EPZ_CUT

    args = parse_args()
    TREE_NAME = args.tree
    EXPECTED_LUMI_INVFB = args.lumi_invfb
    STATE_NORM_MODE = args.state_norm_mode
    APPLY_SIGMA_EPZ_CUT = not args.no_sigma_epz
    SAMPLES["low"]["file"] = os.path.expanduser(args.low_file)
    SAMPLES["high"]["file"] = os.path.expanduser(args.high_file)
    SAMPLES["low"]["visible_xsec_fb"] = args.low_visible_xsec_fb
    SAMPLES["high"]["visible_xsec_fb"] = args.high_visible_xsec_fb
    print("\n=== Υ(1S) overlay of reconstructed |t| estimators in each panel ===")
    print(f"STATE_NORM_MODE = {STATE_NORM_MODE}")
    print(f"APPLY_SIGMA_EPZ_CUT = {APPLY_SIGMA_EPZ_CUT}\n")

    branches = [
        "reco_t_eX", "reco_t_BABE", "reco_t_eXBABE", "reco_t_proj",
        "reco_pT", "reco_pT2", "weight_visible", "weight_proj_visible"
    ]
    sample_data = {k: load_sample_arrays(k, branches) for k in ["low", "high"]}

    fig, axs = plt.subplots(3, 2, figsize=(15, 16), sharex=True)
    axs = axs.flatten()

    summary_rows = []

    for i, method_cfg in enumerate(METHODS):
        ax = axs[i]
        max_density = 0.0

        for sample_key in ["low", "high"]:
            cfg = SAMPLES[sample_key]
            data = sample_data[sample_key]

            arr = data[method_cfg["branch"]]
            wgt = data[method_cfg["weight"]]

            valid = np.isfinite(arr) & np.isfinite(wgt) & (wgt > 0) & (arr >= T_BINS[0]) & (arr <= T_BINS[-1])

            if np.count_nonzero(valid) == 0:
                continue

            counts, errs = make_hist(arr[valid], T_BINS, wgt[valid])
            density = counts / np.diff(T_BINS)
            density_err = errs / np.diff(T_BINS)

            hep.histplot(density, bins=T_BINS, yerr=density_err, histtype="step",
                         color=cfg["color"], ax=ax, linewidth=2.0,
                         label=rf"{cfg['label']}")
            hep.histplot(density, bins=T_BINS, histtype="fill",
                         color=cfg["color"], alpha=0.12, ax=ax)

            max_density = max(max_density, float(np.max(density)) if len(density) else 0.0)

            summary_rows.append({
                "method": method_cfg["key"],
                "sample": sample_key,
                "label": cfg["label"],
                "campaign": cfg["campaign"],
                "sum_weight": float(np.sum(wgt[valid])),
                "n_entries": int(np.count_nonzero(valid)),
                "w_min": cfg["w_min"],
                "w_max": cfg["w_max"],
            })

        ax.set_title(method_cfg["label"], fontsize=16, pad=8)
        ax.set_xlabel(method_cfg["xlabel"], fontsize=13)
        ax.set_ylabel(r"$dN/d|t|$ [GeV$^{-2}$]", fontsize=13)
        ax.set_yscale("log")
        ax.set_xlim(T_BINS[0], T_BINS[-1])
        ax.set_ylim(0.1, max_density * 100.0 if max_density > 0 else 10.0)
        draw_common_labels(ax)
        ax.legend(frameon=False, fontsize=10, loc="upper right")

    # Hide last empty pad
    if len(axs) > len(METHODS):
        axs[-1].axis("off")

    fig.suptitle(r"Reconstructed $|t|$ estimators for $\Upsilon(1S)$, overlaid by beam energy",
                 fontsize=20, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])

    out_png = args.out_prefix + ".png"
    fig.savefig(out_png, dpi=300)
    print(f"Saved {out_png}")

    import csv
    out_csv = out_png.replace(".png", "_summary.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Saved {out_csv}")


if __name__ == "__main__":
    main()
