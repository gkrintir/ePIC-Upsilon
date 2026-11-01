#!/usr/bin/env python3
# Response-matrix unfolding for dN/d|t| in exclusive Upsilon production.
# Background is set to zero. ZEUS overlay is the black exclusive-Upsilon reference distribution.

import argparse
from pathlib import Path as _BranchingPath
import sys as _branching_sys
_branching_sys.path.insert(0, str(_BranchingPath(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument, branching_fractions_from_args,
                                branching_fraction_metadata)
import os
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from statistical_modes import (
    add_statistical_arguments, configure_statistics, require_projectable_weights, check_response_support,
    sample_unfolded_statistics, covariance_to_csv, scale_covariance, save_unfolded_toys,
    statistical_metadata, mode_annotation,
    identifiable_intervals, draw_summary_intervals, summary_interval_bounds,
    draw_absolute_intervals, sampling_interval, gaussian_interval_scale,
)

import uproot
# Shared explicit reconstructed-proton selection; truth denominators are unchanged.
from pathlib import Path as _ProtonPath
import sys as _proton_sys
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass

UPSILON = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0234, "color": "royalblue"},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191, "color": "crimson"},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218, "color": "forestgreen"},
}

METHODS = {
    "eX": "reco_t_eX",
    "BABE": "reco_t_BABE",
    "eXBABE": "reco_t_eXBABE",
    "proj": "reco_t_proj",
    "pT2": "reco_pT",
}


def parse_edges(s):
    return np.asarray([float(x) for x in s.split(",") if x.strip()], dtype=float)


def branch_exists(path, tree_name, branch):
    with uproot.open(os.path.expanduser(path)) as source:
        return tree_name in source and branch in source[tree_name].keys()


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def response_matrix(true_x, reco_x, weights, edges, reco_selected):
    nb = len(edges) - 1
    truth_hist, truth_err = hist_weighted(true_x, weights, edges)

    H = np.zeros((nb, nb), dtype=float)
    H2 = np.zeros_like(H)

    ib = np.digitize(true_x, edges) - 1
    jb = np.digitize(reco_x, edges) - 1
    # Match the final-bin right-edge inclusion of np.histogram.
    ib[np.asarray(true_x) == edges[-1]] = nb - 1
    jb[np.asarray(reco_x) == edges[-1]] = nb - 1
    valid = reco_selected & (ib >= 0) & (ib < nb) & (jb >= 0) & (jb < nb)

    for i, j, w in zip(ib[valid], jb[valid], weights[valid]):
        H[i, j] += w
        H2[i, j] += w*w

    A = np.divide(H, truth_hist[:, None], out=np.zeros_like(H), where=truth_hist[:, None] > 0)
    eff = A.sum(axis=1)
    return H, np.sqrt(H2), truth_hist, truth_err, A, eff


def purity_stability(H):
    diag = np.diag(H)
    reco_sum = H.sum(axis=0)
    truth_sum = H.sum(axis=1)

    purity = np.divide(diag, reco_sum, out=np.zeros_like(diag), where=reco_sum > 0)
    stability = np.divide(diag, truth_sum, out=np.zeros_like(diag), where=truth_sum > 0)
    return purity, stability


def iterative_bayes(data_reco, A, prior_counts, n_iterations=4):
    A = np.asarray(A, dtype=float)
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
        valid = denom > 0
        M[:, valid] = A[:, valid] * prior[:, None] / denom[valid][None, :]
        unfolded = M @ data
        unfolded = np.divide(unfolded, eff, out=np.zeros_like(unfolded), where=eff > 0)
        if unfolded.sum() > 0:
            prior = unfolded / unfolded.sum()

    return unfolded


def toy_uncertainty(data_reco, data_err, A, prior_counts, args):
    return sample_unfolded_statistics(
        data_reco, data_err,
        lambda varied: iterative_bayes(varied, A, prior_counts, args.n_iterations),
        args.n_toys, args.seed + 1000, getattr(args, "stat_mode", "mc"),
    ).error


def weighted_quantile(values, q, weights):
    values = np.asarray(values)
    weights = np.asarray(weights)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[mask]
    weights = weights[mask]
    if len(values) == 0:
        return np.full(len(q), np.nan)
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    cdf /= cdf[-1]
    return np.interp(q, cdf, values)


def make_reco_t(df, method):
    if method == "pT2":
        return np.asarray(df["reco_pT"], dtype=float)**2
    return np.asarray(df[METHODS[method]], dtype=float)


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

def scale_model_to_zeus(model_x, model_y, zeus):
    # Match model bin centers to ZEUS histogram bin centers.  If ZEUS has no
    # uncertainties, use a 20% relative weighting for the scale fit.
    xs = np.asarray(model_x)
    ys = np.asarray(model_y)
    zx = zeus_centers(zeus)
    scale_num = 0.0
    scale_den = 0.0
    for idx_row, r in zeus.iterrows():
        idx = int(np.argmin(np.abs(xs - zx[idx_row])))
        m = ys[idx]
        if not np.isfinite(m) or m <= 0:
            continue
        err = 0.0
        if "err_low" in zeus.columns and "err_high" in zeus.columns:
            err = 0.5 * (float(r["err_low"]) + float(r["err_high"]))
        if err <= 0:
            err = max(float(r["dNdt_GeVminus2"]) * 0.2, 1.0)
        w = 1.0 / (err*err)
        scale_num += float(r["dNdt_GeVminus2"]) * m * w
        scale_den += m * m * w
    return scale_num / scale_den if scale_den > 0 else 1.0


def plot_response(H, edges, out_prefix, state=None):
    fig, ax = plt.subplots(figsize=(7,6))
    mesh = ax.pcolormesh(edges, edges, H)
    fig.colorbar(mesh, ax=ax, label="weighted selected events")
    ax.set_xlabel(r"$|t|_{\rm reco}$ [GeV$^2$]")
    ax.set_ylabel(r"$|t|_{\rm true}$ [GeV$^2$]")
    ax.set_title((rf"$\Upsilon({state})$: " if state else "") + "Selected-event response counts")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_response_counts.png", dpi=250)
    plt.close(fig)


def plot_closure(summary, out_prefix):
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8.5,8), sharex=True, gridspec_kw={"height_ratios":[3,1]})
    x = summary["t_center"]
    xerr = 0.5 * (summary["t_high"] - summary["t_low"])

    mode = str(summary["stat_mode"].iloc[0]) if "stat_mode" in summary else "mc"
    unfolded_label = "unfolded projection (data stat.)" if mode == "projected" else "unfolded test subset (MC stat.)"
    draw_summary_intervals(ax, summary, x, xerr, "dNdt_unfolded_visible", "dNdt_unfolded_visible_err",
                           "dNdt_unfolded_visible_interval_low", "dNdt_unfolded_visible_interval_high",
                           fmt="o", label=unfolded_label)
    ax.errorbar(x, summary["dNdt_truth_test_visible"], xerr=xerr, yerr=summary["dNdt_truth_test_visible_err"], fmt="s", markerfacecolor="white", color="black", label="truth reference (MC stat.)")
    ax.set_yscale("log")
    ax.set_ylabel(r"$dN_{ee}/d|t|$ [GeV$^{-2}$]")
    title = (f"Expected yields at {summary['yield_lumi_invfb'].iloc[0]:g} fb$^{{-1}}$"
             if mode == "projected" else "Held-out test subset yields")
    if "state" in summary:
        title = rf"$\Upsilon({summary['state'].iloc[0]})$: " + title
    ax.set_title(title, fontsize=14)
    ax.legend(frameon=False, fontsize=11)

    rax.axhline(1.0, color="gray", linestyle="--")
    draw_summary_intervals(rax, summary, x, xerr, "closure_unfold_over_truth", "closure_err",
                           "closure_interval_low", "closure_interval_high", fmt="o", color="purple")
    rax.set_xlabel(r"$|t|$ [GeV$^2$]")
    rax.set_ylabel("unfold/truth")
    lower, upper = summary_interval_bounds(summary, "closure_unfold_over_truth", "closure_err",
                                           "closure_interval_low", "closure_interval_high")
    finite = np.isfinite(lower) & np.isfinite(upper)
    low = min(0.3, float(np.min(lower[finite]))) if np.any(finite) else 0.3
    high = max(1.7, float(np.max(upper[finite]))) if np.any(finite) else 1.7
    nominal = summary["closure_unfold_over_truth"].to_numpy(float)
    if np.any(np.isfinite(nominal)):
        low, high = min(low, float(np.nanmin(nominal))), max(high, float(np.nanmax(nominal)))
    margin = 0.04 * (high - low)
    rax.set_ylim(low - margin, high + margin)
    fig.text(0.5, 0.012, mode_annotation(mode, summary["count_interval"].iloc[0], summary["interval_coverage"].iloc[0]),
             ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(f"{out_prefix}_closure.png", dpi=250)
    plt.close(fig)


def plot_zeus_overlay(summary, zeus, args, out_prefix, stats=None):
    if args.state != "1S":
        raise ValueError("The ZEUS reference measures Upsilon(1S), not the requested excited state")
    fig, ax = plt.subplots(figsize=(8.5,6.8))

    x = summary["t_center"].to_numpy()
    xerr = 0.5 * (summary["t_high"] - summary["t_low"]).to_numpy()
    y = summary["dNdt_unfolded_visible"].to_numpy()
    yerr = summary["dNdt_unfolded_visible_err"].to_numpy()
    mode = getattr(args, "stat_mode", "mc")
    toy_density = (stats.toys / summary["bin_width"].to_numpy()[None, :]
                   if stats is not None and len(stats.toys) > 1 else None)
    transformed = toy_density
    identifiable = (summary["unfolding_identifiable"].to_numpy(bool)
                    if "unfolding_identifiable" in summary else np.ones(len(summary), dtype=bool))
    if mode == "projected" and args.zeus_mode != "absolute" and not np.all(identifiable):
        raise ValueError("A normalized ZEUS shape comparison requires identifiable truth bins over its full support")

    if args.zeus_mode == "scale_to_zeus":
        scale = scale_model_to_zeus(x, y, zeus)
        y_plot = scale * y
        yerr_plot = scale * yerr
        if mode == "projected" and toy_density is not None:
            transformed = np.asarray([scale_model_to_zeus(x, toy, zeus) * toy for toy in toy_density])
            yerr_plot = np.std(transformed, axis=0, ddof=1)
        label = rf"ePIC unfolded $Υ({args.state})$ shape $\times {scale:.2g}$"
        ylabel = r"$dN/d|t|$ [GeV$^{-2}$]"
    elif args.zeus_mode == "unit_area":
        widths = (summary["t_high"] - summary["t_low"]).to_numpy()
        area = np.nansum(y * widths)
        y_plot = y / area if area > 0 else y
        yerr_plot = yerr / area if area > 0 else yerr
        if mode == "projected" and toy_density is not None:
            toy_area = np.sum(toy_density * widths[None, :], axis=1)
            valid_toys = toy_area > 0
            if np.sum(valid_toys) < 2:
                raise ValueError("Unit-area projected shape requires at least two positive-yield toys")
            transformed = toy_density[valid_toys] / toy_area[valid_toys, None]
            yerr_plot = np.std(transformed, axis=0, ddof=1)
        # Normalize ZEUS too
        # Use approximate bin widths inferred from ZEUS binning 0-0.2,0.2-0.5,0.5-1,1-5
        label = "ePIC unfolded unit-area shape"
        ylabel = r"unit-normalized shape"
    else:
        y_plot = y
        yerr_plot = yerr
        label = (r"ePIC projected visible $dN/d|t|$ (data stat.)" if mode == "projected"
                 else r"ePIC test-subset visible $dN/d|t|$ (MC stat.)")
        ylabel = r"$dN_{ee}/d|t|$ [GeV$^{-2}$]"

    if mode == "projected":
        y_plot = np.asarray(y_plot, dtype=float).copy()
        y_plot[~identifiable] = np.nan
    if mode == "projected" and args.count_interval != "standard_deviation":
        lower, upper = sampling_interval(transformed, coverage=args.interval_coverage, axis=0)
        lower[~identifiable], upper[~identifiable] = np.nan, np.nan
        draw_absolute_intervals(ax, x, y_plot, lower, upper, xerr=xerr, fmt="o-", color="royalblue", capsize=3, label=label)
    else:
        if mode == "projected":
            yerr_plot = gaussian_interval_scale(args.interval_coverage) * yerr_plot
        ax.errorbar(x, y_plot, xerr=xerr, yerr=yerr_plot, fmt="o-", color="royalblue", capsize=3, label=label)

    z = zeus.copy()
    zx = zeus_centers(z)
    zxerr = zeus_xerr(z)
    if args.zeus_mode == "unit_area":
        if {"t_low", "t_high"}.issubset(z.columns):
            widths = (z["t_high"] - z["t_low"]).to_numpy()
        else:
            widths = np.array([0.2, 0.3, 0.5, 4.0])[:len(z)]
        area = np.sum(z["dNdt_GeVminus2"].to_numpy() * widths)
        z_y = z["dNdt_GeVminus2"].to_numpy() / area
        z_elo = z["err_low"].to_numpy() / area if "err_low" in z.columns else np.zeros(len(z))
        z_ehi = z["err_high"].to_numpy() / area if "err_high" in z.columns else np.zeros(len(z))
    else:
        z_y = z["dNdt_GeVminus2"].to_numpy()
        z_elo = z["err_low"].to_numpy() if "err_low" in z.columns else np.zeros(len(z))
        z_ehi = z["err_high"].to_numpy() if "err_high" in z.columns else np.zeros(len(z))

    # Draw the ZEUS exclusive component as a step histogram when bin edges are available.
    if {"t_low", "t_high"}.issubset(z.columns):
        step_edges = np.r_[z["t_low"].to_numpy(), z["t_high"].to_numpy()[-1]]
        step_vals = np.r_[z_y, z_y[-1]]
        ax.step(step_edges, step_vals, where="post", color="black", linestyle=":", linewidth=2.2,
                label=r"ZEUS exclusive $\Upsilon(1S)$ component")
        ax.errorbar(zx, z_y, xerr=zxerr, yerr=np.vstack([z_elo, z_ehi]),
                    fmt="none", ecolor="black", capsize=2)
    else:
        ax.errorbar(zx, z_y, yerr=np.vstack([z_elo, z_ehi]),
                    fmt="ko", capsize=3, label=r"ZEUS exclusive $\Upsilon(1S)$ component")

    ax.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.set_xlim(args.t_edges[0], args.t_edges[-1])
    ax.set_ylim(args.y_min, args.y_max)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(r"$|t|$ unfolding and ZEUS exclusive $Υ$ comparison")
    fig.text(0.5, 0.012, mode_annotation(mode, args.count_interval, args.interval_coverage),
             ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(f"{out_prefix}_zeus_overlay.png", dpi=250)
    plt.close(fig)


def plot_resolution(resolution, out_prefix, state=None):
    fig, ax = plt.subplots(figsize=(8,6))
    ax.errorbar(resolution["t_center"], resolution["median_delta_t"], yerr=resolution["core_sigma_t"], fmt="o", capsize=3)
    ax.axhline(0, color="gray", linestyle="--")
    ax.set_xlabel(r"$|t|_{\rm true}$ [GeV$^2$]")
    ax.set_ylabel(r"$|t|_{\rm reco}-|t|_{\rm true}$ [GeV$^2$]")
    ax.set_title((rf"$\Upsilon({state})$: " if state else "") + r"$|t|$ resolution")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}_t_resolution.png", dpi=250)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--out-prefix", default="unfold_t")
    p.add_argument("--state", choices=["1S","2S","3S"], default="1S")
    p.add_argument("--method", choices=list(METHODS.keys()), default="eX")
    p.add_argument("--t-edges", required=True)

    p.add_argument("--electron-energy", type=float, default=18.0)
    p.add_argument("--w-min", type=float, default=None)
    p.add_argument("--w-max", type=float, default=None)
    p.add_argument("--q2-cut", type=float, default=None)
    p.add_argument("--visible-xsec-fb", type=float, default=None)
    p.add_argument("--state-norm-mode", choices=["raw", "target_visible_xsec"],
                   default="target_visible_xsec",
                   help="raw preserves stored weights after analysis cuts; target_visible_xsec rescales the selected truth phase space.")
    p.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    p.add_argument("--lumi-invfb", type=float, default=10.0)
    p.add_argument("--weight-convention", choices=["visible_decay","inclusive_upsilon"], default="visible_decay")

    p.add_argument("--sigma-epz-center", type=float, default=None)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    p.add_argument("--mass-window-sigma", type=float, default=3.0)
    p.add_argument("--true-mass-window", type=float, default=0.2)

    p.add_argument("--train-fraction", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--n-iterations", type=int, default=4)
    p.add_argument("--n-toys", type=int, default=300)

    p.add_argument("--zeus-csv", default="")
    p.add_argument("--zeus-mode", choices=["scale_to_zeus","unit_area","absolute"], default="scale_to_zeus")
    p.add_argument("--y-min", type=float, default=0.1)
    p.add_argument("--y-max", type=float, default=1000.0)

    add_statistical_arguments(p)
    add_branching_fraction_argument(p)
    args = p.parse_args()
    if not 0.0 < args.train_fraction < 1.0:
        p.error("--train-fraction must be strictly between 0 and 1")
    if not np.isfinite(args.lumi_invfb) or args.lumi_invfb <= 0.0:
        p.error("--lumi-invfb must be finite and positive")
    configure_statistics(args, p)
    args.t_edges = parse_edges(args.t_edges)

    if args.sigma_epz_center is None:
        args.sigma_epz_center = 2.0 * args.electron_energy
    if args.sigma_epz_halfwidth is None:
        args.sigma_epz_halfwidth = 0.2 * args.electron_energy

    infile = os.path.expanduser(args.input)
    props = dict(UPSILON[args.state], br_ee=branching_fractions_from_args(args)[args.state])
    has_true_q2 = branch_exists(infile, args.tree, "true_Q2")

    branches = ["weight","true_mass","true_t","true_W","reco_mass","reco_W","sigma_epz","reco_pT"]
    if METHODS[args.method] not in branches and args.method != "pT2":
        branches.append(METHODS[args.method])
    if has_true_q2:
        branches.append("true_Q2")

    finite_columns = list(branches)
    branches += proton_branches()
    with uproot.open(infile) as source:
        tree = source[args.tree]
        arr = {name: tree[name].array(library="np") for name in dict.fromkeys(branches)}
    df = pd.DataFrame({k: np.asarray(v, dtype=float) for k,v in arr.items()})
    require_projectable_weights(df["weight"].to_numpy(), args)
    proton_tag_mask(arr)  # Validate tag schema without removing truth failures.
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=finite_columns)

    mass = props["mass"]
    truth_state = (df["true_mass"] >= mass - args.true_mass_window) & (df["true_mass"] <= mass + args.true_mass_window)
    if args.q2_cut is not None and has_true_q2:
        truth_state = truth_state & (df["true_Q2"] >= 0.0) & (df["true_Q2"] <= args.q2_cut)
    if args.w_min is not None:
        truth_state = truth_state & (df["true_W"] >= args.w_min)
    if args.w_max is not None:
        truth_state = truth_state & (df["true_W"] <= args.w_max)

    sumw_state_raw = df.loc[truth_state, "weight"].sum()
    if args.tuple_lumi_invfb <= 0 or not np.isfinite(sumw_state_raw) or sumw_state_raw <= 0:
        raise ValueError("A positive tuple luminosity and positive selected sum of weights are required")
    inferred_xsec = sumw_state_raw / args.tuple_lumi_invfb
    if args.state_norm_mode == "target_visible_xsec":
        if args.visible_xsec_fb is None or not np.isfinite(args.visible_xsec_fb) or args.visible_xsec_fb <= 0:
            raise ValueError("target_visible_xsec requires a positive --visible-xsec-fb")
        scale = args.visible_xsec_fb / inferred_xsec
    else:
        scale = 1.0
    df["w_visible"] = df["weight"] * scale

    # Truth phase space in t.
    truth_in_t = truth_state & (df["true_t"] >= args.t_edges[0]) & (df["true_t"] <= args.t_edges[-1])

    reco_t = make_reco_t(df, args.method)
    tagged_proton = proton_tag_mask(df)
    reco_mass_window = args.mass_window_sigma * props["sigma"]
    reco_selected = (
        truth_in_t & tagged_proton &
        (df["reco_mass"] >= mass - reco_mass_window) &
        (df["reco_mass"] <= mass + reco_mass_window) &
        (np.abs(df["sigma_epz"] - args.sigma_epz_center) <= args.sigma_epz_halfwidth) &
        (reco_t >= args.t_edges[0]) & (reco_t <= args.t_edges[-1])
    )
    if args.w_min is not None:
        reco_selected = reco_selected & (df["reco_W"] >= args.w_min)
    if args.w_max is not None:
        reco_selected = reco_selected & (df["reco_W"] <= args.w_max)

    sub = df[truth_in_t].copy().reset_index(drop=True)
    reco_t_sub = reco_t[truth_in_t.to_numpy()]
    reco_sel_sub = reco_selected[truth_in_t].to_numpy()

    rng = np.random.default_rng(args.seed)
    train = rng.random(len(sub)) < args.train_fraction
    test = ~train
    if not np.any(train) or not np.any(test):
        raise ValueError("Both training and held-out test samples must contain events")

    true_train = sub.loc[train, "true_t"].to_numpy()
    reco_train = reco_t_sub[train]
    w_train = sub.loc[train, "w_visible"].to_numpy()
    sel_train = reco_sel_sub[train]

    true_test = sub.loc[test, "true_t"].to_numpy()
    reco_test = reco_t_sub[test]
    w_test = sub.loc[test, "w_visible"].to_numpy()
    sel_test = reco_sel_sub[test]

    H, Herr, truth_train_hist, truth_train_err, A, eff = response_matrix(true_train, reco_train, w_train, args.t_edges, sel_train)

    data_reco, data_reco_err = hist_weighted(reco_test[sel_test], w_test[sel_test], args.t_edges)
    truth_test_hist, truth_test_err = hist_weighted(true_test, w_test, args.t_edges)
    source_data_reco = data_reco.copy()
    source_truth_test = truth_test_hist.copy()
    data_reco *= args.projection_scale
    data_reco_err *= args.projection_scale
    truth_test_hist *= args.projection_scale
    truth_test_err *= args.projection_scale
    check_response_support(A, truth_test_hist, data_reco, args)

    unfold = lambda varied: iterative_bayes(varied, A, truth_train_hist, args.n_iterations)
    unfolded = unfold(data_reco)
    stats = sample_unfolded_statistics(data_reco, data_reco_err, unfold, args.n_toys, args.seed + 1000, args.stat_mode,
                                      coverage=args.interval_coverage)
    unfolded_err = stats.error
    template_stats = (sample_unfolded_statistics(data_reco, data_reco_err, unfold, args.n_toys, args.seed + 101000, "mc")
                      if args.stat_mode == "projected" else stats)

    widths = np.diff(args.t_edges)
    centers = 0.5 * (args.t_edges[:-1] + args.t_edges[1:])

    br_corr = 1.0 / props["br_ee"] if args.weight_convention == "visible_decay" else 1.0

    dNdt = unfolded / widths
    dNdt_err = unfolded_err / widths
    truth_dNdt = truth_test_hist / widths
    truth_dNdt_err = truth_test_err / widths

    # MC retains held-out yields; projected spectra represent the full target L.
    test_lumi_invfb = args.yield_lumi_invfb
    dsigma_ep_dt_fb = unfolded / test_lumi_invfb / widths * br_corr
    dsigma_ep_dt_err_fb = unfolded_err / test_lumi_invfb / widths * br_corr
    interval_low, interval_high = identifiable_intervals(stats, eff, args.stat_mode)
    dsigma_interval_low = interval_low / test_lumi_invfb / widths * br_corr
    dsigma_interval_high = interval_high / test_lumi_invfb / widths * br_corr

    purity, stability = purity_stability(H)

    closure = np.divide(unfolded, truth_test_hist, out=np.full_like(unfolded, np.nan), where=truth_test_hist > 0)
    closure_err = np.divide(unfolded_err, truth_test_hist, out=np.full_like(unfolded_err, np.nan), where=truth_test_hist > 0)
    closure_interval_low = np.divide(interval_low, truth_test_hist, out=np.full_like(interval_low, np.nan), where=truth_test_hist > 0)
    closure_interval_high = np.divide(interval_high, truth_test_hist, out=np.full_like(interval_high, np.nan), where=truth_test_hist > 0)

    # Resolution summary from selected events only.
    res_rows = []
    for i in range(len(widths)):
        b = (sub["true_t"].to_numpy() >= args.t_edges[i]) & (sub["true_t"].to_numpy() < args.t_edges[i+1]) & reco_sel_sub
        if not np.any(b):
            res_rows.append({"bin": i, "t_low": args.t_edges[i], "t_high": args.t_edges[i+1], "t_center": centers[i],
                             "median_delta_t": np.nan, "core_sigma_t": np.nan, "median_rel_delta_t": np.nan, "core_sigma_rel": np.nan})
            continue
        delta = reco_t_sub[b] - sub["true_t"].to_numpy()[b]
        rel = delta / sub["true_t"].to_numpy()[b]
        ww = sub["w_visible"].to_numpy()[b]
        q16, q50, q84 = weighted_quantile(delta, [0.16,0.5,0.84], ww)
        rq16, rq50, rq84 = weighted_quantile(rel, [0.16,0.5,0.84], ww)
        res_rows.append({
            "bin": i, "t_low": args.t_edges[i], "t_high": args.t_edges[i+1], "t_center": centers[i],
            "median_delta_t": q50, "core_sigma_t": 0.5*(q84-q16),
            "median_rel_delta_t": rq50, "core_sigma_rel": 0.5*(rq84-rq16),
            "n_selected": int(np.sum(b)), "sumw_selected": float(ww.sum())
        })

    resolution = pd.DataFrame(res_rows)

    summary = pd.DataFrame({
        "state": args.state,
        "method": args.method,
        "stat_mode": args.stat_mode,
        "count_interval": args.count_interval,
        "interval_coverage": args.interval_coverage,
        "unfolding_identifiable": eff > 0,
        "yield_lumi_invfb": args.yield_lumi_invfb,
        "t_low": args.t_edges[:-1],
        "t_high": args.t_edges[1:],
        "t_center": centers,
        "bin_width": widths,
        "efficiency": eff,
        "purity": purity,
        "stability": stability,
        "truth_test_visible": truth_test_hist,
        "reco_test_visible": data_reco,
        "source_truth_test_visible": source_truth_test,
        "source_reco_test_visible": source_data_reco,
        "reco_template_mc_err": data_reco_err,
        "truth_reference_mc_err": truth_test_err,
        "unfolded_visible": unfolded,
        "unfolded_visible_err": unfolded_err,
        "unfolded_visible_cov_std": stats.standard_deviation,
        "unfolded_visible_interval_low": interval_low,
        "unfolded_visible_interval_high": interval_high,
        "unfolded_visible_toy_q16": stats.q16,
        "unfolded_visible_toy_q50": stats.q50,
        "unfolded_visible_toy_q84": stats.q84,
        "unfolded_visible_toy_mean": stats.mean,
        "unfolded_template_mc_err": template_stats.standard_deviation,
        "dNdt_unfolded_visible": dNdt,
        "dNdt_unfolded_visible_err": dNdt_err,
        "dNdt_unfolded_visible_interval_low": interval_low / widths,
        "dNdt_unfolded_visible_interval_high": interval_high / widths,
        "dNdt_truth_test_visible": truth_dNdt,
        "dNdt_truth_test_visible_err": truth_dNdt_err,
        "dsigma_ep_dt_inclusive_fb_per_GeV2": dsigma_ep_dt_fb,
        "dsigma_ep_dt_inclusive_err_fb_per_GeV2": dsigma_ep_dt_err_fb,
        "dsigma_ep_dt_inclusive_interval_low_fb_per_GeV2": dsigma_interval_low,
        "dsigma_ep_dt_inclusive_interval_high_fb_per_GeV2": dsigma_interval_high,
        "closure_unfold_over_truth": closure,
        "closure_err": closure_err,
        "closure_interval_low": closure_interval_low,
        "closure_interval_high": closure_interval_high,
    })

    norm = pd.DataFrame([{
        "state": args.state,
        "method": args.method,
        "sumw_state_raw": sumw_state_raw,
        "inferred_visible_xsec_fb": inferred_xsec,
        "target_visible_xsec_fb": args.visible_xsec_fb,
        "weight_scale": scale,
        "state_norm_mode": args.state_norm_mode,
        "train_fraction": args.train_fraction,
        "test_fraction": 1.0 - args.train_fraction,
        "effective_test_lumi_invfb": args.source_test_lumi_invfb,
        **statistical_metadata(args),
        "br_ee": props["br_ee"],
        **branching_fraction_metadata(args),
        "n_iterations": args.n_iterations,
        "n_toys": args.n_toys,
        "background_assumption": 0.0,
        "proton_selection": proton_selection_metadata(),
    }])

    prefix = args.out_prefix
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(f"{prefix}_t_unfold_summary.csv", index=False)
    norm.to_csv(f"{prefix}_normalization.csv", index=False)
    labels = [f"{lo:g}-{hi:g}" for lo, hi in zip(args.t_edges[:-1], args.t_edges[1:])]
    sigma_factors = br_corr / args.yield_lumi_invfb / widths
    covariance_to_csv(f"{prefix}_unfolded_covariance.csv", stats.covariance, labels)
    covariance_to_csv(f"{prefix}_density_covariance.csv", scale_covariance(stats.covariance, 1.0 / widths), labels)
    covariance_to_csv(f"{prefix}_cross_section_covariance.csv", scale_covariance(stats.covariance, sigma_factors), labels)
    covariance_to_csv(f"{prefix}_template_mc_covariance.csv", template_stats.covariance, labels)
    covariance_to_csv(f"{prefix}_template_mc_cross_section_covariance.csv", scale_covariance(template_stats.covariance, sigma_factors), labels)
    save_unfolded_toys(f"{prefix}_unfolded_toys.npz", stats, args.t_edges, args)
    resolution.to_csv(f"{prefix}_t_resolution_summary.csv", index=False)
    pd.DataFrame(H).to_csv(f"{prefix}_response_counts.csv", index=False)
    pd.DataFrame(A).to_csv(f"{prefix}_response_prob.csv", index=False)

    plot_response(H, args.t_edges, prefix, state=args.state)
    plot_closure(summary, prefix)
    plot_resolution(resolution, prefix, state=args.state)

    use_zeus_reference = bool(args.zeus_csv and args.state == "1S")
    if args.zeus_csv and not use_zeus_reference:
        print(f"Skipping ZEUS Upsilon(1S) reference for Upsilon({args.state}); no ZEUS overlay produced.")
    if use_zeus_reference:
        zeus = pd.read_csv(os.path.expanduser(args.zeus_csv))
        plot_zeus_overlay(summary, zeus, args, prefix, stats)

    print("=== Normalization ===")
    print(norm.to_string(index=False))
    print("\n=== t unfolding summary ===")
    print(summary.to_string(index=False))
    print("\n=== t resolution summary ===")
    print(resolution.to_string(index=False))
    print(f"\nWrote {prefix}_t_unfold_summary.csv")
    print(f"Wrote {prefix}_t_resolution_summary.csv")
    print(f"Wrote {prefix}_closure.png")
    print(f"Wrote {prefix}_response_counts.png")
    if use_zeus_reference:
        print(f"Wrote {prefix}_zeus_overlay.png")


if __name__ == "__main__":
    main()
