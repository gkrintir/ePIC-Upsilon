#!/usr/bin/env python3
"""
Compare all reconstructed |t| estimators side by side for two beam energies.

Outputs
-------
<out-prefix>_summary.csv
<out-prefix>_per_bin.csv
<out-prefix>_residuals.png
<out-prefix>_metrics.png
<out-prefix>_closure.png

The script avoids pandas on purpose, for SWAN environments where a local pandas
module can shadow the real package.

Methods compared
----------------
eX       reco_t_eX
BABE     reco_t_BABE
eXBABE   reco_t_eXBABE
proj     reco_t_proj, with optional projected corrected weights
pT2      reco_pT^2

Default physics convention
--------------------------
Weights are visible e+e- weights.  No BR correction is applied.
Optional state normalization rescales the selected truth state so sum(weight)/L
matches the visible e+e- cross section supplied for each beam energy.
"""

import argparse
import csv
import os
import tempfile
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument,
                                branching_fractions_from_args,
                                branching_fraction_metadata)

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
    "1S": {"mass": 9.46030, "sigma": 0.053, "br": 0.0234},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br": 0.0191},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br": 0.0218},
}

METHODS = {
    "eX": {
        "branch": "reco_t_eX",
        "label": r"$t_{eX}$",
        "needs": ["reco_t_eX"],
        "extra_cut": None,
    },
    "BABE": {
        "branch": "reco_t_BABE",
        "label": r"$t_{\rm BABE}$",
        "needs": ["reco_t_BABE"],
        "extra_cut": None,
    },
    "eXBABE": {
        "branch": "reco_t_eXBABE",
        "label": r"$t_{eX{\rm BABE}}$",
        "needs": ["reco_t_eXBABE"],
        "extra_cut": None,
    },
    "proj": {
        "branch": "reco_t_proj",
        "label": r"$t_{\rm proj}$",
        "needs": ["reco_t_proj", "reco_t_eXBABE", "weight_proj"],
        "extra_cut": "proj",
    },
    "pT2": {
        "branch": "reco_pT",
        "label": r"$p_{T,\Upsilon}^{2}$",
        "needs": ["reco_pT"],
        "extra_cut": None,
    },
}


def parse_edges(s):
    return np.asarray([float(x) for x in s.split(",") if x.strip()], dtype=float)


def write_csv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def weighted_quantile(values, quantiles, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[ok]
    weights = weights[ok]
    if len(values) == 0:
        return np.full(len(quantiles), np.nan)
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return np.interp(quantiles, cdf, values)


def make_reco_t(data, method):
    if method == "pT2":
        return np.asarray(data["reco_pT"], dtype=float) ** 2
    return np.asarray(data[METHODS[method]["branch"]], dtype=float)


def response_matrix(true_t, reco_t, weights, edges, reco_selected):
    nb = len(edges) - 1
    truth_hist, truth_err = hist_weighted(true_t, weights, edges)

    H = np.zeros((nb, nb), dtype=float)
    ib = np.digitize(true_t, edges) - 1
    jb = np.digitize(reco_t, edges) - 1
    ok = reco_selected & (ib >= 0) & (ib < nb) & (jb >= 0) & (jb < nb)

    for i, j, w in zip(ib[ok], jb[ok], weights[ok]):
        H[i, j] += w

    A = np.divide(H, truth_hist[:, None], out=np.zeros_like(H), where=truth_hist[:, None] > 0)
    eff = A.sum(axis=1)
    return H, truth_hist, truth_err, A, eff


def purity_stability(H):
    diag = np.diag(H)
    reco_sum = H.sum(axis=0)
    truth_sum = H.sum(axis=1)
    purity = np.divide(diag, reco_sum, out=np.zeros_like(diag), where=reco_sum > 0)
    stability = np.divide(diag, truth_sum, out=np.zeros_like(diag), where=truth_sum > 0)
    return purity, stability


def iterative_bayes(data_reco, A, prior_counts, n_iterations):
    data = np.asarray(data_reco, dtype=float)
    n_truth, n_reco = A.shape
    eff = A.sum(axis=1)

    prior = np.asarray(prior_counts, dtype=float)
    if prior.sum() <= 0:
        prior = np.ones(n_truth)
    prior = prior / prior.sum()

    unfolded = np.zeros(n_truth, dtype=float)
    for _ in range(int(n_iterations)):
        denom = A.T @ prior
        M = np.zeros((n_truth, n_reco), dtype=float)
        good = denom > 0
        M[:, good] = A[:, good] * prior[:, None] / denom[good][None, :]
        unfolded = M @ data
        unfolded = np.divide(unfolded, eff, out=np.zeros_like(unfolded), where=eff > 0)
        if unfolded.sum() > 0:
            prior = unfolded / unfolded.sum()
    return unfolded


def load_sample(path, sample_name, electron_energy, visible_xsec_fb, args):
    props = dict(UPSILON[args.state], br=branching_fractions_from_args(args)[args.state])

    base_columns = [
        "weight", "weight_proj",
        "true_mass", "true_t", "true_W", "true_Q2",
        "reco_mass", "reco_W", "sigma_epz",
        "reco_pT",
        "reco_t_eX", "reco_t_BABE", "reco_t_eXBABE", "reco_t_proj",
    ]

    base_columns += proton_branches()
    with uproot.open(os.path.expanduser(path)) as root_file:
        tree = root_file[args.tree]
        arrs = {name: tree[name].array(library="np") for name in base_columns}
    data = {k: np.asarray(v, dtype=float) for k, v in arrs.items()}

    # Validate the full producer contract before any existing finite-value filter.
    proton_tag_mask(data)
    finite = np.ones(len(data["weight"]), dtype=bool)
    for name, v in data.items():
        if name not in proton_branches():
            finite &= np.isfinite(v)
    data = {k: v[finite] for k, v in data.items()}

    mass = props["mass"]
    truth_state = (
        (data["true_mass"] >= mass - args.truth_mass_window)
        & (data["true_mass"] <= mass + args.truth_mass_window)
    )

    if args.q2_cut is not None:
        truth_state &= (data["true_Q2"] >= 0.0) & (data["true_Q2"] <= args.q2_cut)

    if args.state_norm_mode == "target_visible_xsec":
        sumw = float(np.sum(data["weight"][truth_state]))
        target_sumw = visible_xsec_fb * args.lumi_invfb
        scale = target_sumw / sumw if sumw > 0 else 1.0
    elif args.state_norm_mode == "raw":
        sumw = float(np.sum(data["weight"][truth_state]))
        scale = 1.0
    else:
        raise ValueError(f"Unknown state norm mode {args.state_norm_mode}")

    data["weight_visible"] = data["weight"] * scale
    data["weight_proj_visible"] = data["weight_proj"] * scale

    meta = {
        "proton_selection": proton_selection_metadata(),
        **branching_fraction_metadata(args),
        "branching_fraction_correction_applied": False,
        "sample": sample_name,
        "raw_visible_xsec_fb": sumw / args.lumi_invfb if args.lumi_invfb > 0 else np.nan,
        "target_visible_xsec_fb": visible_xsec_fb,
        "state_weight_scale": scale,
        "entries_after_finite": len(data["weight"]),
        "entries_truth_state": int(np.sum(truth_state)),
    }

    return data, truth_state, meta


def analyze_method(data, truth_state, sample_name, electron_energy, method, args, w_min, w_max):
    props = dict(UPSILON[args.state], br=branching_fractions_from_args(args)[args.state])
    reco_t = make_reco_t(data, method)
    true_t = data["true_t"]

    sigma_center = 2.0 * electron_energy
    sigma_half = args.sigma_epz_halfwidth if args.sigma_epz_halfwidth is not None else args.sigma_epz_halfwidth_fraction * electron_energy

    truth_mask = (
        truth_state
        & (true_t >= args.t_edges[0])
        & (true_t <= args.t_edges[-1])
        & (data["true_W"] >= w_min)
        & (data["true_W"] <= w_max)
    )

    reco_mask = (
        truth_mask
        & proton_tag_mask(data)
        & (data["reco_mass"] >= args.reco_mass_min)
        & (data["reco_mass"] <= args.reco_mass_max)
        & (np.abs(data["sigma_epz"] - sigma_center) <= sigma_half)
        & (data["reco_W"] >= w_min)
        & (data["reco_W"] <= w_max)
        & (reco_t >= args.t_edges[0])
        & (reco_t <= args.t_edges[-1])
    )

    if method == "proj":
        reco_mask &= (data["weight_proj"] > 0) & (data["reco_t_eXBABE"] > 0)

    # Statistical response metrics use visible event weights.  The projected corrected
    # yield is additionally stored for the proj method.
    w_stat = data["weight_visible"]
    w_corr = data["weight_proj_visible"] if method == "proj" else data["weight_visible"]

    H, truth_hist, truth_err, A, eff = response_matrix(true_t[truth_mask], reco_t[truth_mask], w_stat[truth_mask], args.t_edges, reco_mask[truth_mask])
    purity, stability = purity_stability(H)

    reco_hist, reco_err = hist_weighted(reco_t[reco_mask], w_stat[reco_mask], args.t_edges)

    # Independent closure split.
    rng = np.random.default_rng(args.seed + hash((sample_name, method)) % 100000)
    selected_indices = np.where(truth_mask)[0]
    train_bool_local = rng.random(len(selected_indices)) < args.train_fraction
    train_idx = selected_indices[train_bool_local]
    test_idx = selected_indices[~train_bool_local]

    reco_sel_train = reco_mask[train_idx]
    reco_sel_test = reco_mask[test_idx]

    H_train, truth_train, truth_train_err, A_train, eff_train = response_matrix(
        true_t[train_idx], reco_t[train_idx], w_stat[train_idx], args.t_edges, reco_sel_train
    )

    reco_test, reco_test_err = hist_weighted(reco_t[test_idx][reco_sel_test], w_stat[test_idx][reco_sel_test], args.t_edges)
    truth_test, truth_test_err = hist_weighted(true_t[test_idx], w_stat[test_idx], args.t_edges)

    unfolded = iterative_bayes(reco_test, A_train, truth_train, args.n_iterations)
    closure = np.divide(unfolded, truth_test, out=np.full_like(unfolded, np.nan), where=truth_test > 0)

    # Residuals on selected reconstructed events.
    sel = reco_mask
    delta = reco_t[sel] - true_t[sel]
    rel = np.divide(delta, true_t[sel], out=np.full_like(delta, np.nan), where=true_t[sel] > 0)
    ww = w_stat[sel]

    q16, q50, q84 = weighted_quantile(delta, [0.16, 0.5, 0.84], ww)
    rq16, rq50, rq84 = weighted_quantile(rel, [0.16, 0.5, 0.84], ww)
    core = 0.5 * (q84 - q16)
    core_rel = 0.5 * (rq84 - rq16)

    min_purity = float(np.nanmin(purity)) if len(purity) else np.nan
    min_stability = float(np.nanmin(stability)) if len(stability) else np.nan
    min_ps = float(np.nanmin(np.minimum(purity, stability))) if len(purity) else np.nan
    closure_abs_dev = np.nanmax(np.abs(closure - 1.0)) if np.any(np.isfinite(closure)) else np.nan

    summary = {
        "proton_selection": proton_selection_metadata()["mode"],
        "br_ee": props["br"],
        "branching_fraction_convention": branching_fraction_metadata(args)["branching_fraction_convention"],
        "branching_fraction_correction_applied": False,
        "sample": sample_name,
        "method": method,
        "method_label": METHODS[method]["label"],
        "w_min": w_min,
        "w_max": w_max,
        "sum_truth_visible": float(np.sum(w_stat[truth_mask])),
        "sum_reco_visible": float(np.sum(w_stat[reco_mask])),
        "sum_reco_corrected": float(np.sum(w_corr[reco_mask])),
        "average_efficiency": float(np.sum(w_stat[reco_mask]) / np.sum(w_stat[truth_mask])) if np.sum(w_stat[truth_mask]) > 0 else np.nan,
        "min_purity": min_purity,
        "min_stability": min_stability,
        "min_purity_stability": min_ps,
        "median_delta_t": float(q50),
        "core_sigma_t": float(core),
        "median_rel_delta_t": float(rq50),
        "core_sigma_rel": float(core_rel),
        "closure_max_abs_deviation": float(closure_abs_dev),
        "n_selected_reco_entries": int(np.sum(reco_mask)),
        "n_truth_entries": int(np.sum(truth_mask)),
    }

    per_bin_rows = []
    for i in range(len(args.t_edges) - 1):
        per_bin_rows.append({
            "proton_selection": proton_selection_metadata()["mode"],
            "sample": sample_name,
            "method": method,
            "t_low": args.t_edges[i],
            "t_high": args.t_edges[i+1],
            "t_center": 0.5 * (args.t_edges[i] + args.t_edges[i+1]),
            "truth_visible": truth_hist[i],
            "reco_visible": reco_hist[i],
            "efficiency": eff[i],
            "purity": purity[i],
            "stability": stability[i],
            "closure": closure[i],
        })

    residual_payload = {
        "delta": delta,
        "rel": rel,
        "weights": ww,
        "median": q50,
        "core": core,
    }

    return summary, per_bin_rows, residual_payload


def plot_residuals(results, args):
    methods = args.methods
    samples = ["low", "high"]
    fig, axes = plt.subplots(len(methods), 2, figsize=(12, 3.0 * len(methods)), squeeze=False)

    for i, method in enumerate(methods):
        for j, sample in enumerate(samples):
            ax = axes[i][j]
            payload = results[(sample, method)]["residual"]
            delta = payload["delta"]
            weights = payload["weights"]

            if len(delta) > 0:
                lo = np.nanpercentile(delta, 1)
                hi = np.nanpercentile(delta, 99)
                if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
                    lo, hi = -2, 2
                pad = 0.1 * (hi - lo)
                bins = np.linspace(lo - pad, hi + pad, 60)
                ax.hist(delta, bins=bins, weights=weights, histtype="step", color="black")
                ax.axvline(0, color="gray", linestyle="--")
                ax.axvline(payload["median"], color="tab:red", linestyle="-")
                ax.text(
                    0.04, 0.92,
                    rf"median={payload['median']:.3g}" + "\n" + rf"core={payload['core']:.3g}",
                    transform=ax.transAxes, ha="left", va="top", fontsize=9
                )

            label = args.low_label if sample == "low" else args.high_label
            ax.set_title(f"{label}, {METHODS[method]['label']}")
            ax.set_xlabel(r"$|t|_{\rm reco}-|t|_{\rm true}$ [GeV$^2$]")
            ax.set_ylabel("Weighted events")
            ax.set_yscale("log")
            ax.set_ylim(bottom=0.5)

    fig.tight_layout()
    fig.savefig(args.out_prefix + "_residuals.png", dpi=250)
    plt.close(fig)


def plot_metrics(summary_rows, args):
    methods = args.methods
    sample_names = ["low", "high"]
    x = np.arange(len(methods))
    width = 0.36

    metrics = [
        ("average_efficiency", "Average efficiency"),
        ("min_purity_stability", "Minimum purity/stability"),
        ("median_delta_t", r"Median $\Delta t$ [GeV$^2$]"),
        ("core_sigma_t", r"Core $\sigma(\Delta t)$ [GeV$^2$]"),
    ]

    rowmap = {(r["sample"], r["method"]): r for r in summary_rows}

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes = axes.flatten()

    for ax, (metric, ylabel) in zip(axes, metrics):
        low_vals = [rowmap[("low", m)][metric] for m in methods]
        high_vals = [rowmap[("high", m)][metric] for m in methods]

        ax.bar(x - width/2, low_vals, width, label=args.low_label)
        ax.bar(x + width/2, high_vals, width, label=args.high_label)

        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=25)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(frameon=False)

    fig.suptitle(r"Comparison of reconstructed $|t|$ methods in two beam-energy samples", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out_prefix + "_metrics.png", dpi=250)
    plt.close(fig)


def plot_closure(per_bin_rows, args):
    methods = args.methods
    samples = ["low", "high"]
    fig, axes = plt.subplots(len(methods), 2, figsize=(12, 3.0 * len(methods)), squeeze=False)

    for i, method in enumerate(methods):
        for j, sample in enumerate(samples):
            ax = axes[i][j]
            rows = [r for r in per_bin_rows if r["sample"] == sample and r["method"] == method]
            x = np.array([r["t_center"] for r in rows])
            xerr = 0.5 * np.array([r["t_high"] - r["t_low"] for r in rows])
            y = np.array([r["closure"] for r in rows], dtype=float)

            ax.axhline(1.0, color="gray", linestyle="--")
            ax.errorbar(x, y, xerr=xerr, fmt="o", color="black")
            ax.set_ylim(0.0, 2.0)
            label = args.low_label if sample == "low" else args.high_label
            ax.set_title(f"{label}, {METHODS[method]['label']}")
            ax.set_xlabel(r"$|t|_{\rm true}$ bin [GeV$^2$]")
            ax.set_ylabel("Unfolded / truth")

    fig.tight_layout()
    fig.savefig(args.out_prefix + "_closure.png", dpi=250)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    add_branching_fraction_argument(ap)
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
    ap.add_argument("--reco-mass-min", type=float, default=8.0)
    ap.add_argument("--reco-mass-max", type=float, default=12.0)
    ap.add_argument("--q2-cut", type=float, default=0.001)
    ap.add_argument("--low-electron-energy", type=float, default=10.0)
    ap.add_argument("--high-electron-energy", type=float, default=18.0)
    ap.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    ap.add_argument("--sigma-epz-halfwidth-fraction", type=float, default=0.2)
    # Separate W ranges for the two samples.
    # If --w-min/--w-max are provided, they are used as a legacy common range.
    ap.add_argument("--w-min", type=float, default=None, help="Legacy common W min for both samples")
    ap.add_argument("--w-max", type=float, default=None, help="Legacy common W max for both samples")
    ap.add_argument("--low-w-min", type=float, default=20.0)
    ap.add_argument("--low-w-max", type=float, default=26.5)
    ap.add_argument("--high-w-min", type=float, default=37.0)
    ap.add_argument("--high-w-max", type=float, default=114.0)
    ap.add_argument("--t-edges", default="0,0.5,1,5")
    ap.add_argument("--train-fraction", type=float, default=0.5)
    ap.add_argument("--n-iterations", type=int, default=4)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out-prefix", default="compare_t_methods_two_energies")
    args = ap.parse_args()
    branching_fractions_from_args(args)

    args.methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    args.t_edges = parse_edges(args.t_edges)

    if args.w_min is not None:
        args.low_w_min = args.w_min
        args.high_w_min = args.w_min
    if args.w_max is not None:
        args.low_w_max = args.w_max
        args.high_w_max = args.w_max

    for m in args.methods:
        if m not in METHODS:
            raise ValueError(f"Unknown method {m}. Available methods are {list(METHODS)}")

    low_data, low_truth, low_meta = load_sample(args.low_file, "low", args.low_electron_energy, args.low_visible_xsec_fb, args)
    high_data, high_truth, high_meta = load_sample(args.high_file, "high", args.high_electron_energy, args.high_visible_xsec_fb, args)

    print("=== Sample normalization ===")
    for m in [low_meta, high_meta]:
        print(m)

    summary_rows = []
    per_bin_rows = []
    results = {}

    for sample, data, truth, Ee, w_min, w_max in [
        ("low", low_data, low_truth, args.low_electron_energy, args.low_w_min, args.low_w_max),
        ("high", high_data, high_truth, args.high_electron_energy, args.high_w_min, args.high_w_max),
    ]:
        for method in args.methods:
            print(f"Processing {sample}, method {method}, W in [{w_min}, {w_max}]")
            summary, bins, residual = analyze_method(data, truth, sample, Ee, method, args, w_min, w_max)
            summary_rows.append(summary)
            per_bin_rows.extend(bins)
            results[(sample, method)] = {"summary": summary, "bins": bins, "residual": residual}

    write_csv(args.out_prefix + "_summary.csv", summary_rows)
    write_csv(args.out_prefix + "_per_bin.csv", per_bin_rows)

    plot_residuals(results, args)
    plot_metrics(summary_rows, args)
    plot_closure(per_bin_rows, args)

    print("\n=== Summary ===")
    for r in summary_rows:
        print(
            f"{r['sample']:>4s} {r['method']:>6s} "
            f"eff={r['average_efficiency']:.3f} "
            f"minPS={r['min_purity_stability']:.3f} "
            f"median_dt={r['median_delta_t']:.3f} "
            f"core={r['core_sigma_t']:.3f} "
            f"closureMaxDev={r['closure_max_abs_deviation']:.3f}"
        )

    print(f"\nWrote {args.out_prefix}_summary.csv")
    print(f"Wrote {args.out_prefix}_per_bin.csv")
    print(f"Wrote {args.out_prefix}_residuals.png")
    print(f"Wrote {args.out_prefix}_metrics.png")
    print(f"Wrote {args.out_prefix}_closure.png")


if __name__ == "__main__":
    main()
