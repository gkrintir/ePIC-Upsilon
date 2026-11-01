#!/usr/bin/env python3
"""
Side-by-side reconstructed |t| distributions for two beam energies.

This is intentionally simpler than the response study. It plots the reconstructed
distributions passing the selected reco-level criteria, one row per |t| method
and two columns for the two beam energies.

No BR correction is applied. Weights are visible e+e- weights, with optional
state normalization to the supplied visible cross sections.
"""

import argparse
import csv
import os
import tempfile
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join(tempfile.gettempdir(), f"{os.environ.get('USER', 'user')}_mplconfig")
)

import uproot
import numpy as np
import matplotlib

matplotlib.use("Agg")
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

METHODS = {
    "eX": {"branch": "reco_t_eX", "label": r"$t_{eX}$"},
    "BABE": {"branch": "reco_t_BABE", "label": r"$t_{\rm BABE}$"},
    "eXBABE": {"branch": "reco_t_eXBABE", "label": r"$t_{eX{\rm BABE}}$"},
    "proj": {"branch": "reco_t_proj", "label": r"$t_{\rm proj}$"},
    "pT2": {"branch": "reco_pT", "label": r"$p_{T,\Upsilon}^{2}$"},
}


def parse_edges(s):
    return np.asarray([float(x) for x in s.split(",") if x.strip()], dtype=float)


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def reco_t(data, method):
    if method == "pT2":
        return np.asarray(data["reco_pT"], dtype=float) ** 2
    return np.asarray(data[METHODS[method]["branch"]], dtype=float)


def load_sample(path, sample_name, Ee, xsec_fb, args):
    props = UPSILON[args.state]
    cols = [
        "weight", "weight_proj",
        "true_mass",
        "reco_mass", "reco_W", "sigma_epz", "reco_pT",
        "reco_t_eX", "reco_t_BABE", "reco_t_eXBABE", "reco_t_proj",
    ]
    cols += proton_branches()
    if args.q2_cut is not None:
        cols.append("true_Q2")
    with uproot.open(os.path.expanduser(path)) as root_file:
        tree = root_file[args.tree]
        arr = {name: tree[name].array(library="np") for name in cols}
    data = {k: np.asarray(v, dtype=float) for k, v in arr.items()}

    # Validate every producer tag before the legacy finite-value filter.
    proton_tag_mask(data)
    finite = np.ones(len(data["weight"]), dtype=bool)
    for name, v in data.items():
        if name not in proton_branches():
            finite &= np.isfinite(v)
    data = {k: v[finite] for k, v in data.items()}

    truth_state = (
        (data["true_mass"] >= props["mass"] - args.truth_mass_window)
        & (data["true_mass"] <= props["mass"] + args.truth_mass_window)
    )

    if args.q2_cut is not None:
        truth_state &= (data["true_Q2"] >= 0.0) & (data["true_Q2"] <= args.q2_cut)

    if args.state_norm_mode == "target_visible_xsec":
        sumw = float(np.sum(data["weight"][truth_state]))
        scale = (xsec_fb * args.lumi_invfb) / sumw if sumw > 0 else 1.0
    else:
        sumw = float(np.sum(data["weight"][truth_state]))
        scale = 1.0

    data["weight_visible"] = data["weight"] * scale
    data["weight_proj_visible"] = data["weight_proj"] * scale

    meta = {
        "proton_selection": proton_selection_metadata()["mode"],
        "proton_selection_contract": proton_selection_metadata()["contract"],
        "proton_selection_predicate": proton_selection_metadata()["predicate"],
        "sample": sample_name,
        "raw_visible_xsec_fb": sumw / args.lumi_invfb if args.lumi_invfb > 0 else np.nan,
        "target_visible_xsec_fb": xsec_fb,
        "state_weight_scale": scale,
    }
    return data, truth_state, meta


def selection_mask(data, truth_state, method, Ee, w_min, w_max, args):
    sigma_half = args.sigma_epz_halfwidth
    if sigma_half is None:
        sigma_half = args.sigma_epz_halfwidth_fraction * Ee
    sigma_center = 2.0 * Ee

    rt = reco_t(data, method)

    mask = (
        truth_state
        & proton_tag_mask(data)
        & (data["reco_mass"] >= args.reco_mass_min)
        & (data["reco_mass"] <= args.reco_mass_max)
        & (np.abs(data["sigma_epz"] - sigma_center) <= sigma_half)
        & (data["reco_W"] >= w_min)
        & (data["reco_W"] <= w_max)
        & (rt >= args.t_edges[0])
        & (rt <= args.t_edges[-1])
    )

    if method == "proj":
        mask &= (data["weight_proj"] > 0) & (data["reco_t_eXBABE"] > 0)

    return mask, rt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--low-file", required=True)
    ap.add_argument("--high-file", required=True)
    ap.add_argument("--low-label", default="10x100")
    ap.add_argument("--high-label", default="18x275")
    ap.add_argument("--tree", default="UpsilonTree")
    ap.add_argument("--state", choices=list(UPSILON.keys()), default="1S")
    ap.add_argument("--methods", default="eX,BABE,eXBABE,proj,pT2")
    ap.add_argument("--low-visible-xsec-fb", type=float, default=104.793)
    ap.add_argument("--high-visible-xsec-fb", type=float, default=327.584)
    ap.add_argument("--lumi-invfb", type=float, default=10.0)
    ap.add_argument("--state-norm-mode", choices=["raw", "target_visible_xsec"], default="target_visible_xsec")
    ap.add_argument("--truth-mass-window", type=float, default=0.20)
    ap.add_argument("--q2-cut", type=float, default=None,
                    help="Optional true Q² upper cut in GeV², applied to truth and reconstructed selections")
    ap.add_argument("--reco-mass-min", type=float, default=8.0)
    ap.add_argument("--reco-mass-max", type=float, default=12.0)
    ap.add_argument("--low-electron-energy", type=float, default=10.0)
    ap.add_argument("--high-electron-energy", type=float, default=18.0)
    ap.add_argument("--low-w-min", type=float, default=20.0)
    ap.add_argument("--low-w-max", type=float, default=26.5)
    ap.add_argument("--high-w-min", type=float, default=37.0)
    ap.add_argument("--high-w-max", type=float, default=114.0)
    ap.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    ap.add_argument("--sigma-epz-halfwidth-fraction", type=float, default=0.2)
    ap.add_argument("--t-edges", default="0,0.5,1,5")
    ap.add_argument("--shape-normalize", action="store_true", help="Normalize each histogram to unit area before plotting")
    ap.add_argument("--out-prefix", default="reco_t_distributions_two_energies")
    args = ap.parse_args()

    args.methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    args.t_edges = parse_edges(args.t_edges)

    low_data, low_truth, low_meta = load_sample(args.low_file, args.low_label, args.low_electron_energy, args.low_visible_xsec_fb, args)
    high_data, high_truth, high_meta = load_sample(args.high_file, args.high_label, args.high_electron_energy, args.high_visible_xsec_fb, args)

    print("=== Normalization ===")
    print(low_meta)
    print(high_meta)

    fig, axes = plt.subplots(len(args.methods), 2, figsize=(12, 3.0 * len(args.methods)), squeeze=False)
    summary_rows = []

    for i, method in enumerate(args.methods):
        for j, sample in enumerate(["low", "high"]):
            ax = axes[i][j]

            if sample == "low":
                data, truth, Ee = low_data, low_truth, args.low_electron_energy
                w_min, w_max = args.low_w_min, args.low_w_max
                label = args.low_label
                color = "dodgerblue"
            else:
                data, truth, Ee = high_data, high_truth, args.high_electron_energy
                w_min, w_max = args.high_w_min, args.high_w_max
                label = args.high_label
                color = "darkorange"

            mask, rt = selection_mask(data, truth, method, Ee, w_min, w_max, args)
            weights = data["weight_proj_visible"] if method == "proj" else data["weight_visible"]

            vals = rt[mask]
            wgts = weights[mask]

            counts, errs = hist_weighted(vals, wgts, args.t_edges)

            if args.shape_normalize:
                total = np.sum(counts)
                if total > 0:
                    # Convert bin yields to a unit-area density, including unequal widths.
                    normalization = total * np.diff(args.t_edges)
                    counts = counts / normalization
                    errs = errs / normalization
                ylabel = r"Unit-normalized $dN/d|t|$"
            else:
                widths = np.diff(args.t_edges)
                counts = counts / widths
                errs = errs / widths
                ylabel = r"$dN/d|t|$ [GeV$^{-2}$]"

            hep.histplot(counts, bins=args.t_edges, yerr=errs, histtype="step", color=color, ax=ax, linewidth=2)
            hep.histplot(counts, bins=args.t_edges, histtype="fill", color=color, alpha=0.15, ax=ax)

            ax.set_title(f"{label}, {METHODS[method]['label']}, W=[{w_min:g},{w_max:g}] GeV")
            ax.set_xlabel(r"$|t|_{\rm reco}$ [GeV$^2$]")
            ax.set_ylabel(ylabel)
            ax.set_yscale("log")
            ax.set_ylim(bottom=1e-3 if args.shape_normalize else 0.1)
            ax.legend([METHODS[method]["label"]], frameon=False, loc="upper right")

            summary_rows.append({
                "proton_selection": proton_selection_metadata()["mode"],
                "proton_selection_contract": proton_selection_metadata()["contract"],
                "sample": sample,
                "method": method,
                "w_min": w_min,
                "w_max": w_max,
                "sum_weight_visible_selected": float(np.sum(wgts)),
                "n_entries_selected": int(np.sum(mask)),
            })

    fig.suptitle(r"Reconstructed $|t|$ distributions by method and beam energy", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.out_prefix + ".png", dpi=250)
    fig.savefig(args.out_prefix + ".pdf")
    write_csv(args.out_prefix + "_summary.csv", summary_rows)
    write_csv(args.out_prefix + "_normalization.csv", [low_meta, high_meta])

    print(f"Wrote {args.out_prefix}.png")
    print(f"Wrote {args.out_prefix}.pdf")
    print(f"Wrote {args.out_prefix}_summary.csv")


if __name__ == "__main__":
    main()
