#!/usr/bin/env python3
"""
Make a polished 2D correlation plot between central-track and low-Q2 tagger-track
multiplicities.

Default use case
----------------
For the reconstructed Upsilon tuple, this plots

    reco_nTracks_Central  vs  reco_nTracks_Tagger

after a standard exclusive selection.

The plot includes:
  - weighted 2D occupancy map
  - optional integer-bin annotations
  - weighted means and RMS values
  - weighted Pearson correlation
  - dominant class, e.g. Ncentral=3, Ntagger=1
  - ePIC-style label block

Example
-------
python3 plot_central_tagger_correlation.py \
  /eos/user/g/gkrintir/merged_elecs_test_tuples_fix2.root \
  --electron-energy 18 \
  --energy-label "18x275 GeV" \
  --campaign "26.03.1 Campaign" \
  --state 1S \
  --w-min 20 \
  --w-max 140 \
  --out-prefix central_tagger_corr_18x275_1S

For 10x100:
python3 plot_central_tagger_correlation.py \
  /eos/user/g/gkrintir/merged_elecs_test_10t100.root \
  --electron-energy 10 \
  --energy-label "10x100 GeV" \
  --campaign "26.04.1 Campaign" \
  --state 1S \
  --w-min 20 \
  --w-max 26.5 \
  --out-prefix central_tagger_corr_10x100_1S
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
from matplotlib.colors import LogNorm

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


STATES = {
    "all": None,
    "1S": {"mass": 9.46030, "window": 0.20},
    "2S": {"mass": 10.02326, "window": 0.20},
    "3S": {"mass": 10.35520, "window": 0.20},
}


def weighted_mean(values, weights):
    s = np.sum(weights)
    return float(np.sum(values * weights) / s) if s > 0 else np.nan


def weighted_rms(values, weights):
    m = weighted_mean(values, weights)
    s = np.sum(weights)
    return float(np.sqrt(np.sum(weights * (values - m) ** 2) / s)) if s > 0 else np.nan


def weighted_corr(x, y, weights):
    wx = weighted_mean(x, weights)
    wy = weighted_mean(y, weights)
    s = np.sum(weights)
    if s <= 0:
        return np.nan
    cov = np.sum(weights * (x - wx) * (y - wy)) / s
    sx = weighted_rms(x, weights)
    sy = weighted_rms(y, weights)
    if sx <= 0 or sy <= 0:
        return np.nan
    return float(cov / (sx * sy))


def build_cut(args):
    cuts = []

    if args.q2_cut is not None:
        cuts.append(f"true_Q2 >= 0 && true_Q2 <= {args.q2_cut:.12g}")

    if args.state != "all":
        st = STATES[args.state]
        cuts.append(
            f"true_mass >= {st['mass'] - st['window']:.12g} && "
            f"true_mass <= {st['mass'] + st['window']:.12g}"
        )

    if args.apply_reco_selection:
        cuts.extend([
            f"reco_mass >= {args.mass_min:.12g}",
            f"reco_mass <= {args.mass_max:.12g}",
            f"reco_W >= {args.w_min:.12g}",
            f"reco_W <= {args.w_max:.12g}",
        ])

        if args.apply_sigma_epz:
            center = args.sigma_epz_center
            if center is None:
                center = 2.0 * args.electron_energy
            half = args.sigma_epz_halfwidth
            if half is None:
                half = 0.2 * args.electron_energy
            cuts.append(f"abs(sigma_epz - {center:.12g}) <= {half:.12g}")

    return " && ".join(cuts) if cuts else "1"


def selection_mask(data, args):
    """NumPy equivalent of the selection displayed by build_cut()."""
    mask = np.ones(len(data[args.weight_branch]), dtype=bool)
    if args.q2_cut is not None:
        mask &= (data["true_Q2"] >= 0) & (data["true_Q2"] <= args.q2_cut)
    if args.state != "all":
        st = STATES[args.state]
        mask &= (data["true_mass"] >= st["mass"] - st["window"])
        mask &= (data["true_mass"] <= st["mass"] + st["window"])
    if args.apply_reco_selection:
        mask &= (data["reco_mass"] >= args.mass_min) & (data["reco_mass"] <= args.mass_max)
        mask &= (data["reco_W"] >= args.w_min) & (data["reco_W"] <= args.w_max)
        if args.apply_sigma_epz:
            center = args.sigma_epz_center
            if center is None:
                center = 2.0 * args.electron_energy
            half = args.sigma_epz_halfwidth
            if half is None:
                half = 0.2 * args.electron_energy
            mask &= np.abs(data["sigma_epz"] - center) <= half
    return mask


def draw_epic_labels(ax, args, n_entries, sumw, corr, dominant):
    center = args.sigma_epz_center if args.sigma_epz_center is not None else 2.0 * args.electron_energy
    half = args.sigma_epz_halfwidth if args.sigma_epz_halfwidth is not None else 0.2 * args.electron_energy

    state_label = r"$\Upsilon(nS)$" if args.state == "all" else rf"$\Upsilon({args.state})$"

    ax.text(0.035, 0.965, "ePIC Performance", transform=ax.transAxes,
            fontsize=12, fontweight="bold", ha="left", va="top")
    ax.text(0.035, 0.905, args.campaign, transform=ax.transAxes,
            fontsize=10.5, ha="left", va="top")
    ax.text(0.035, 0.85, rf"e+p {state_label}$\to e^+e^-$, {args.energy_label}",
            transform=ax.transAxes, fontsize=10.5, ha="left", va="top")
    ax.text(0.035, 0.795, rf"$W_{{\gamma p}}\in[{args.w_min:g},{args.w_max:g}]$ GeV",
            transform=ax.transAxes, fontsize=10.5, ha="left", va="top")
    ax.text(0.035, 0.74, rf"$|\Sigma(E-p_z)-{center:g}|<{half:g}$ GeV",
            transform=ax.transAxes, fontsize=10.5, ha="left", va="top")
    ax.text(0.035, 0.685, rf"$\rho_w={corr:.3f}$, dominant class: {dominant}",
            transform=ax.transAxes, fontsize=10.5, ha="left", va="top")
    ax.text(0.035, 0.63, rf"entries={n_entries:,}, $\sum w={sumw:.1f}$",
            transform=ax.transAxes, fontsize=10.5, ha="left", va="top")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input", help="Input ROOT file with UpsilonTree")
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--out-prefix", default="central_tagger_correlation")

    p.add_argument("--state", choices=list(STATES.keys()), default="1S")
    p.add_argument("--electron-energy", type=float, default=18.0)
    p.add_argument("--q2-cut", type=float, default=None,
                   help="Optional event-level true_Q2 cut; requires the true_Q2 branch.")
    p.add_argument("--energy-label", default="18x275 GeV")
    p.add_argument("--campaign", default="26.03.1 Campaign")

    p.add_argument("--central-branch", default="reco_nTracks_Central")
    p.add_argument("--tagger-branch", default="reco_nTracks_Tagger")
    p.add_argument("--weight-branch", default="weight")

    p.add_argument("--apply-reco-selection", action="store_true", default=True)
    p.add_argument("--no-reco-selection", dest="apply_reco_selection", action="store_false")
    p.add_argument("--apply-sigma-epz", action="store_true", default=True)
    p.add_argument("--no-sigma-epz", dest="apply_sigma_epz", action="store_false")

    p.add_argument("--mass-min", type=float, default=8.0)
    p.add_argument("--mass-max", type=float, default=12.0)
    p.add_argument("--w-min", type=float, default=20.0)
    p.add_argument("--w-max", type=float, default=140.0)
    p.add_argument("--sigma-epz-center", type=float, default=None)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)

    p.add_argument("--max-central", type=int, default=25)
    p.add_argument("--max-tagger", type=int, default=12)
    p.add_argument("--min-bin-label", type=float, default=0.0,
                   help="Annotate bins with weighted counts above this value. Use 0 to disable.")
    p.add_argument("--normalize", choices=["none", "total"], default="none")
    args = p.parse_args()

    filename = os.path.expanduser(args.input)
    cut = build_cut(args)

    columns = [
        args.central_branch,
        args.tagger_branch,
        args.weight_branch,
        "true_mass",
        "reco_mass",
        "reco_W",
        "sigma_epz",
    ]
    columns += proton_branches()
    # unique preserving order
    columns = list(dict.fromkeys(columns))
    if args.q2_cut is not None:
        columns.append("true_Q2")

    # Check every producer tag before the existing backend applies selection.
    # Keep the archived ROOT and later uproot loaders otherwise unchanged.
    if proton_branches():
        if "uproot" in globals():
            with uproot.open(filename) as tag_source:
                tag_tree = tag_source[args.tree]
                tag_arrays = {name: tag_tree[name].array(library="np") for name in proton_branches()}
        elif "ROOT" in globals():
            tag_arrays = ROOT.RDataFrame(args.tree, filename).AsNumpy(columns=proton_branches())
        else:
            raise RuntimeError("Unsupported track-correlation input backend")
        proton_tag_mask(tag_arrays)

    print("=== Central/tagger correlation ===")
    print(f"Input: {filename}")
    print(f"Cut: {cut}")
    print(f"Branches: {args.central_branch}, {args.tagger_branch}, weight={args.weight_branch}")

    with uproot.open(filename) as source:
        tree = source[args.tree]
        arr = {name: tree[name].array(library="np") for name in columns}
    selected = selection_mask(arr, args)
    arr = {name: values[selected] for name, values in arr.items()}

    selected_proton = proton_tag_mask(arr)
    x = np.asarray(arr[args.central_branch], dtype=float)
    y = np.asarray(arr[args.tagger_branch], dtype=float)
    w = np.asarray(arr[args.weight_branch], dtype=float)

    valid = (
        selected_proton & np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
        & (w > 0)
        & (x >= 0) & (x <= args.max_central)
        & (y >= 0) & (y <= args.max_tagger)
    )
    x = x[valid]
    y = y[valid]
    w = w[valid]

    if len(x) == 0:
        raise RuntimeError("No events passed the selection and plotting range.")

    xbins = np.arange(-0.5, args.max_central + 1.5, 1.0)
    ybins = np.arange(-0.5, args.max_tagger + 1.5, 1.0)

    H, _, _ = np.histogram2d(x, y, bins=[xbins, ybins], weights=w)
    # H is [xbin, ybin], transpose for pcolormesh y vs x.
    Hplot = H.T

    H_for_plot = Hplot.copy()
    if args.normalize == "total" and np.sum(H_for_plot) > 0:
        H_for_plot = H_for_plot / np.sum(H_for_plot)

    positive = H_for_plot[H_for_plot > 0]
    norm = LogNorm(vmin=max(np.min(positive), 1e-6), vmax=np.max(positive)) if len(positive) else None

    max_idx = np.unravel_index(np.argmax(H), H.shape)
    dominant = rf"$N_{{central}}={int(max_idx[0])}$, $N_{{tagger}}={int(max_idx[1])}$"

    corr = weighted_corr(x, y, w)

    fig, ax = plt.subplots(figsize=(10.5, 7.4))
    mesh = ax.pcolormesh(xbins, ybins, H_for_plot, cmap="viridis", norm=norm)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.015)
    cbar.set_label("Weighted events" if args.normalize == "none" else "Fraction of selected yield")

    # Optional bin labels.
    if args.min_bin_label > 0:
        for ix in range(H_for_plot.shape[1]):
            for iy in range(H_for_plot.shape[0]):
                val = H_for_plot[iy, ix]
                if val >= args.min_bin_label:
                    ax.text(ix, iy, f"{val:.0f}", ha="center", va="center",
                            fontsize=8, color="white" if val < np.max(H_for_plot) * 0.5 else "black")

    ax.set_xlabel(r"$N_{\mathrm{central}}$", fontsize=16)
    ax.set_ylabel(r"$N_{\mathrm{tagger}}$", fontsize=16)
    ax.set_xlim(-0.5, args.max_central + 0.5)
    ax.set_ylim(-0.5, args.max_tagger + 0.5)
    ax.set_xticks(np.arange(0, args.max_central + 1, 2))
    ax.set_yticks(np.arange(0, args.max_tagger + 1, 1))
    ax.set_title(r"Central vs low-$Q^2$ tagger track multiplicity", fontsize=18, pad=10)

    draw_epic_labels(ax, args, len(x), float(np.sum(w)), corr, dominant)

    fig.tight_layout()
    out_png = args.out_prefix + ".png"
    out_pdf = args.out_prefix + ".pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")

    # Write compact summary and bin map.
    summary_csv = args.out_prefix + "_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=[
            "proton_selection", "entries", "sum_weight", "mean_central", "rms_central",
            "mean_tagger", "rms_tagger", "weighted_corr",
            "dominant_ncentral", "dominant_ntagger", "dominant_weight"
        ])
        wr.writeheader()
        wr.writerow({
            "proton_selection": proton_selection_metadata()["mode"],
            "entries": len(x),
            "sum_weight": float(np.sum(w)),
            "mean_central": weighted_mean(x, w),
            "rms_central": weighted_rms(x, w),
            "mean_tagger": weighted_mean(y, w),
            "rms_tagger": weighted_rms(y, w),
            "weighted_corr": corr,
            "dominant_ncentral": int(max_idx[0]),
            "dominant_ntagger": int(max_idx[1]),
            "dominant_weight": float(H[max_idx]),
        })
    print(f"Wrote {summary_csv}")

    map_csv = args.out_prefix + "_map.csv"
    with open(map_csv, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["ncentral", "ntagger", "weight"])
        for ix in range(H.shape[0]):
            for iy in range(H.shape[1]):
                if H[ix, iy] > 0:
                    wr.writerow([ix, iy, H[ix, iy]])
    print(f"Wrote {map_csv}")


if __name__ == "__main__":
    main()
