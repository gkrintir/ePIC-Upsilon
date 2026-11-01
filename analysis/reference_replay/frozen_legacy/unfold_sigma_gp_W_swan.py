#!/usr/bin/env python3
"""
Response-matrix unfolding for exclusive Upsilon sigma(gamma p -> Upsilon p) versus W_gamma_p.

This script uses a flattened UpsilonTree with branches:
  weight, true_mass, true_W, true_Q2, reco_mass, reco_W, sigma_epz

Background is set to zero by default. The unfolding is signal-only.

Workflow:
  1. Select a Upsilon state using truth mass.
  2. Renormalize event weights to the visible e+e- eSTARlight cross section for that state.
  3. Split events into response-training and pseudo-data closure samples.
  4. Build a response matrix A_{ij} = P(reco bin j | truth bin i), including inefficiency.
  5. Unfold the pseudo-data reconstructed spectrum with iterative Bayesian unfolding.
  6. Convert unfolded visible ep yields to sigma(gamma p) by dividing by L * photon flux * BR_ee.
  7. Write ROOT-like CSV outputs and plots.

Recommended 18x275 example:
  python3 unfold_sigma_gp_W_swan.py /eos/user/g/gkrintir/merged_elecs_test_tuples_fix.root \
    --state 1S \
    --electron-energy 18 \
    --proton-energy 275 \
    --w-edges "20,37.1429,45.7143,54.2857,62.8571,71.4286,80,88.5714,97.1429,105.714,114.286,122.857,140" \
    --q2-max 0.001 \
    --q2-cut 0.001 \
    --visible-xsec-fb 327.584 \
    --n-iterations 4 \
    --out-prefix unfold_1S_18x275

Recommended 10x100 example:
  python3 unfold_sigma_gp_W_swan.py /eos/user/g/gkrintir/merged_elecs_test_10t100.root \
    --state 1S \
    --electron-energy 10 \
    --proton-energy 100 \
    --w-edges "20,22.1667,24.3333,26.5" \
    --q2-max 0.001 \
    --q2-cut 0.001 \
    --visible-xsec-fb 104.793 \
    --n-iterations 4 \
    --out-prefix unfold_1S_10x100
"""

import argparse
import math
import os
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


ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813

UPSILON = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0234},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218},
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


def photon_flux_dW(W, e_e, e_p, q2_max):
    s_minus_mp2 = 4.0 * e_e * e_p
    y = (W * W - M_P * M_P) / s_minus_mp2
    if y <= 0.0 or y >= 1.0:
        return 0.0

    q2_min = (M_E * M_E * y * y) / (1.0 - y)
    if q2_min <= 0.0 or q2_min >= q2_max:
        return 0.0

    log_term = math.log(q2_max / q2_min)
    dndy = (ALPHA_EM / (2.0 * math.pi)) * (
        ((1.0 + (1.0 - y) ** 2) / y) * log_term
        - (2.0 * (1.0 - y) / y) * (1.0 - q2_min / q2_max)
    )
    return dndy * (2.0 * W / s_minus_mp2)


def integrated_flux(w_low, w_high, e_e, e_p, q2_max, n_steps=500):
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    return float(np.trapz(vals, grid))


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def weighted_response(true_w, reco_w, weights, truth_bins, reco_bins, reco_selected):
    n_true_bins = len(truth_bins) - 1
    n_reco_bins = len(reco_bins) - 1

    truth_hist, truth_err = hist_weighted(true_w, weights, truth_bins)

    H = np.zeros((n_true_bins, n_reco_bins), dtype=float)
    H2 = np.zeros_like(H)

    tbin = np.digitize(true_w, truth_bins) - 1
    rbin = np.digitize(reco_w, reco_bins) - 1

    valid = (
        reco_selected
        & (tbin >= 0) & (tbin < n_true_bins)
        & (rbin >= 0) & (rbin < n_reco_bins)
    )

    for i, j, w in zip(tbin[valid], rbin[valid], weights[valid]):
        H[i, j] += w
        H2[i, j] += w * w

    A = np.divide(H, truth_hist[:, None], out=np.zeros_like(H), where=truth_hist[:, None] > 0)
    eff = A.sum(axis=1)

    return H, np.sqrt(H2), truth_hist, truth_err, A, eff


def iterative_bayes_unfold(data_reco, A, prior_counts, n_iterations=4, background=None):
    """
    Iterative Bayesian unfolding.

    A has shape [truth bin, reco bin] and contains P(reco | truth), including inefficiency.
    data_reco is the reco spectrum after background subtraction.
    background is set to zero if not provided.
    """
    A = np.asarray(A, dtype=float)
    data_reco = np.asarray(data_reco, dtype=float)

    if background is None:
        background = np.zeros_like(data_reco)
    data = data_reco - background
    data = np.clip(data, 0.0, None)

    n_truth, n_reco = A.shape
    eff = A.sum(axis=1)

    prior = np.asarray(prior_counts, dtype=float)
    if prior.sum() <= 0:
        prior = np.ones(n_truth, dtype=float)
    prior = prior / prior.sum()

    unfolded = np.zeros(n_truth, dtype=float)

    for _ in range(int(n_iterations)):
        denom = A.T @ prior
        M = np.zeros((n_truth, n_reco), dtype=float)

        valid_reco = denom > 0
        M[:, valid_reco] = (A[:, valid_reco] * prior[:, None]) / denom[valid_reco][None, :]

        unfolded = M @ data
        unfolded = np.divide(unfolded, eff, out=np.zeros_like(unfolded), where=eff > 0)

        if unfolded.sum() > 0:
            prior = unfolded / unfolded.sum()

    return unfolded


def response_purity_stability(H):
    diag = np.diag(H)
    reco_sum = H.sum(axis=0)
    truth_sum = H.sum(axis=1)

    purity = np.divide(diag, reco_sum, out=np.zeros_like(diag), where=reco_sum > 0)
    stability = np.divide(diag, truth_sum, out=np.zeros_like(diag), where=truth_sum > 0)
    return purity, stability


def make_selection_arrays(df, props, args, has_true_q2):
    mass = props["mass"]

    truth_state = (
        (df["true_mass"] >= mass - args.true_mass_window)
        & (df["true_mass"] <= mass + args.true_mass_window)
    )

    if args.q2_cut is not None and has_true_q2:
        truth_state = truth_state & (df["true_Q2"] >= 0.0) & (df["true_Q2"] <= args.q2_cut)

    reco_mass_window = args.mass_window_sigma * props["sigma"]
    reco_selected = (
        truth_state
        & (df["reco_mass"] >= mass - reco_mass_window)
        & (df["reco_mass"] <= mass + reco_mass_window)
        & (np.abs(df["sigma_epz"] - args.sigma_epz_center) <= args.sigma_epz_halfwidth)
        & (df["reco_W"] >= args.w_edges[0])
        & (df["reco_W"] <= args.w_edges[-1])
    )

    truth_in_range = (
        truth_state
        & (df["true_W"] >= args.w_edges[0])
        & (df["true_W"] <= args.w_edges[-1])
    )

    return truth_state, truth_in_range, reco_selected


def split_mask(n, train_fraction, seed):
    rng = np.random.default_rng(seed)
    train = rng.random(n) < train_fraction
    test = ~train
    return train, test


def convert_to_sigma(unfolded_visible, unfolded_err_visible, args, props, w_edges):
    flux = np.array([
        integrated_flux(w_edges[i], w_edges[i + 1], args.electron_energy, args.proton_energy, args.q2_max)
        for i in range(len(w_edges) - 1)
    ])

    br_corr = 1.0 / props["br_ee"] if args.weight_convention == "visible_decay" else 1.0
    denom = args.lumi_invfb * 1000.0 * flux

    sigma = np.divide(unfolded_visible, denom, out=np.full_like(unfolded_visible, np.nan), where=denom > 0) * br_corr
    sigma_err = np.divide(unfolded_err_visible, denom, out=np.full_like(unfolded_err_visible, np.nan), where=denom > 0) * br_corr

    return sigma, sigma_err, flux


def toy_uncertainty(data_reco, data_reco_err, A, prior_counts, args):
    if args.n_toys <= 0:
        return np.zeros(A.shape[0], dtype=float)

    rng = np.random.default_rng(args.seed + 999)
    toys = []
    for _ in range(args.n_toys):
        varied = rng.normal(data_reco, data_reco_err)
        varied = np.clip(varied, 0.0, None)
        u = iterative_bayes_unfold(varied, A, prior_counts, args.n_iterations)
        toys.append(u)

    if not toys:
        return np.zeros(A.shape[0], dtype=float)

    toys = np.vstack(toys)
    return 0.5 * (np.nanquantile(toys, 0.84, axis=0) - np.nanquantile(toys, 0.16, axis=0))


def plot_matrix(matrix, x_edges, y_edges, outpath, title, xlabel, ylabel, log=False):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    m = np.asarray(matrix)
    if log:
        import matplotlib.colors as colors
        positive = m[m > 0]
        vmin = positive.min() if len(positive) else 1e-12
        mesh = ax.pcolormesh(x_edges, y_edges, m, norm=colors.LogNorm(vmin=vmin, vmax=max(m.max(), vmin)))
    else:
        mesh = ax.pcolormesh(x_edges, y_edges, m)
    fig.colorbar(mesh, ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=250)
    plt.close(fig)


def plot_spectra(summary, outpath):
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})

    x = summary["W_center"]
    xerr = 0.5 * (summary["W_high"] - summary["W_low"])

    ax.errorbar(
        x, summary["sigma_unfolded_pb"], xerr=xerr, yerr=summary["sigma_unfolded_err_pb"],
        fmt="o", color="royalblue", label="unfolded closure"
    )
    ax.errorbar(
        x, summary["sigma_truth_test_pb"], xerr=xerr, yerr=summary["sigma_truth_test_err_pb"],
        fmt="s", color="black", markerfacecolor="white", label="truth test"
    )

    ax.set_yscale("log")
    ax.set_ylabel(r"$\sigma(\gamma p\rightarrow\Upsilon p)$ [pb]")
    ax.legend(frameon=False)

    rax.axhline(1.0, color="gray", linestyle="--")
    rax.errorbar(
        x, summary["closure_unfolded_over_truth"],
        xerr=xerr, yerr=summary["closure_ratio_err"],
        fmt="o", color="purple"
    )
    rax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    rax.set_ylabel("unfold/truth")
    rax.set_ylim(0.5, 1.5)

    fig.tight_layout()
    fig.savefig(outpath, dpi=250)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--out-prefix", default="unfold_sigma_gp_W")
    parser.add_argument("--state", choices=["1S", "2S", "3S"], default="1S")

    parser.add_argument("--electron-energy", type=float, required=True)
    parser.add_argument("--proton-energy", type=float, required=True)
    parser.add_argument("--q2-max", type=float, default=0.001)
    parser.add_argument("--q2-cut", type=float, default=None)

    parser.add_argument("--w-edges", required=True)
    parser.add_argument("--lumi-invfb", type=float, default=10.0)
    parser.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    parser.add_argument("--visible-xsec-fb", type=float, required=True)
    parser.add_argument("--weight-convention", choices=["visible_decay", "inclusive_upsilon"], default="visible_decay")

    parser.add_argument("--sigma-epz-center", type=float, default=None)
    parser.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    parser.add_argument("--mass-window-sigma", type=float, default=3.0)
    parser.add_argument("--true-mass-window", type=float, default=0.2)

    parser.add_argument("--train-fraction", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--n-iterations", type=int, default=4)
    parser.add_argument("--n-toys", type=int, default=300)
    args = parser.parse_args()

    args.w_edges = parse_edges(args.w_edges)
    if args.sigma_epz_center is None:
        args.sigma_epz_center = 2.0 * args.electron_energy
    if args.sigma_epz_halfwidth is None:
        args.sigma_epz_halfwidth = 0.2 * args.electron_energy

    infile = os.path.expanduser(args.input)
    outdir = Path(os.path.dirname(os.path.abspath(args.out_prefix)) or ".")
    outdir.mkdir(parents=True, exist_ok=True)

    props = UPSILON[args.state]
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

    truth_state, truth_in_range, reco_selected = make_selection_arrays(df, props, args, has_true_q2)

    sumw_state_raw = df.loc[truth_state, "weight"].sum()
    inferred_visible_xsec_fb = sumw_state_raw / args.tuple_lumi_invfb
    weight_scale = args.visible_xsec_fb / inferred_visible_xsec_fb

    df["w_visible"] = df["weight"] * weight_scale

    selected_df = df[truth_in_range].copy()
    selected_reco = reco_selected[truth_in_range].to_numpy()
    selected_df = selected_df.reset_index(drop=True)

    train, test = split_mask(len(selected_df), args.train_fraction, args.seed)

    truth_train = selected_df.loc[train, "true_W"].to_numpy()
    reco_train = selected_df.loc[train, "reco_W"].to_numpy()
    w_train = selected_df.loc[train, "w_visible"].to_numpy()
    sel_train = selected_reco[train]

    truth_test = selected_df.loc[test, "true_W"].to_numpy()
    reco_test = selected_df.loc[test, "reco_W"].to_numpy()
    w_test = selected_df.loc[test, "w_visible"].to_numpy()
    sel_test = selected_reco[test]

    H, Herr, truth_train_hist, truth_train_err, A, eff = weighted_response(
        truth_train, reco_train, w_train, args.w_edges, args.w_edges, sel_train
    )

    data_reco, data_reco_err = hist_weighted(reco_test[sel_test], w_test[sel_test], args.w_edges)
    truth_test_hist, truth_test_err = hist_weighted(truth_test, w_test, args.w_edges)

    unfolded = iterative_bayes_unfold(data_reco, A, truth_train_hist, args.n_iterations)
    unfolded_toy_err = toy_uncertainty(data_reco, data_reco_err, A, truth_train_hist, args)

    sigma_unfold, sigma_unfold_err, flux = convert_to_sigma(unfolded, unfolded_toy_err, args, props, args.w_edges)
    sigma_truth, sigma_truth_err, _ = convert_to_sigma(truth_test_hist, truth_test_err, args, props, args.w_edges)

    purity, stability = response_purity_stability(H)

    W_low = args.w_edges[:-1]
    W_high = args.w_edges[1:]
    W_center = 0.5 * (W_low + W_high)

    closure = np.divide(sigma_unfold, sigma_truth, out=np.full_like(sigma_unfold, np.nan), where=sigma_truth > 0)
    closure_err = np.divide(sigma_unfold_err, sigma_truth, out=np.full_like(sigma_unfold_err, np.nan), where=sigma_truth > 0)

    summary = pd.DataFrame({
        "state": args.state,
        "W_low": W_low,
        "W_high": W_high,
        "W_center": W_center,
        "flux_int": flux,
        "efficiency": eff,
        "purity": purity,
        "stability": stability,
        "truth_train_visible": truth_train_hist,
        "truth_test_visible": truth_test_hist,
        "reco_test_visible": data_reco,
        "unfolded_visible": unfolded,
        "unfolded_visible_err_toy": unfolded_toy_err,
        "sigma_unfolded_pb": sigma_unfold,
        "sigma_unfolded_err_pb": sigma_unfold_err,
        "sigma_truth_test_pb": sigma_truth,
        "sigma_truth_test_err_pb": sigma_truth_err,
        "closure_unfolded_over_truth": closure,
        "closure_ratio_err": closure_err,
    })

    norm = pd.DataFrame([{
        "state": args.state,
        "sumw_state_raw": sumw_state_raw,
        "inferred_visible_xsec_fb": inferred_visible_xsec_fb,
        "target_visible_xsec_fb": args.visible_xsec_fb,
        "weight_scale": weight_scale,
        "train_fraction": args.train_fraction,
        "n_iterations": args.n_iterations,
        "n_toys": args.n_toys,
        "q2_max_flux": args.q2_max,
        "q2_cut": args.q2_cut if args.q2_cut is not None else np.nan,
        "background_assumption": 0.0,
    }])

    prefix = Path(args.out_prefix)
    summary.to_csv(str(prefix) + "_unfolded_summary.csv", index=False)
    norm.to_csv(str(prefix) + "_normalization.csv", index=False)
    pd.DataFrame(H, index=[f"{a:g}-{b:g}" for a, b in zip(W_low, W_high)],
                 columns=[f"{a:g}-{b:g}" for a, b in zip(W_low, W_high)]).to_csv(str(prefix) + "_response_counts.csv")
    pd.DataFrame(A, index=[f"{a:g}-{b:g}" for a, b in zip(W_low, W_high)],
                 columns=[f"{a:g}-{b:g}" for a, b in zip(W_low, W_high)]).to_csv(str(prefix) + "_response_prob.csv")

    plot_matrix(
        H, args.w_edges, args.w_edges, str(prefix) + "_response_counts.png",
        "Response counts", r"$W_{\rm reco}$ [GeV]", r"$W_{\rm truth}$ [GeV]", log=True
    )
    plot_matrix(
        A, args.w_edges, args.w_edges, str(prefix) + "_response_probability.png",
        r"Response probability $P(W_{\rm reco}|W_{\rm truth})$",
        r"$W_{\rm reco}$ [GeV]", r"$W_{\rm truth}$ [GeV]", log=False
    )
    plot_spectra(summary, str(prefix) + "_closure.png")

    print("=== Normalization ===")
    print(norm.to_string(index=False))
    print("\n=== Unfolding summary ===")
    print(summary.to_string(index=False))
    print(f"\nWrote {prefix}_unfolded_summary.csv")
    print(f"Wrote {prefix}_normalization.csv")
    print(f"Wrote {prefix}_response_counts.csv")
    print(f"Wrote {prefix}_response_prob.csv")
    print(f"Wrote {prefix}_response_counts.png")
    print(f"Wrote {prefix}_response_probability.png")
    print(f"Wrote {prefix}_closure.png")


if __name__ == "__main__":
    main()
