#!/usr/bin/env python3
# Hybrid |t| unfolding with low/high beam-energy samples reweighted to a common W_gamma_p spectrum.
# This version deliberately does NOT import pandas.  It is meant for SWAN setups where a
# local or broken pandas module shadows the real package.

import argparse
from pathlib import Path as _BranchingPath
import sys as _branching_sys
_branching_sys.path.insert(0, str(_BranchingPath(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument, branching_fractions_from_args,
                                branching_fraction_metadata)
import csv
import os
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.join(tempfile.gettempdir(), f"{os.environ.get('USER', 'user')}_mplconfig")
)

import ROOT
# Shared explicit reconstructed-proton selection; truth denominators are unchanged.
from pathlib import Path as _ProtonPath
import sys as _proton_sys
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


UPSILON = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0234},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218},
}


class InputValidationError(ValueError):
    """Raised when an input sample cannot support the requested measurement."""

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
    f = ROOT.TFile.Open(os.path.expanduser(path))
    if not f or f.IsZombie():
        return False
    t = f.Get(tree_name)
    ok = bool(t and t.GetBranch(branch))
    f.Close()
    return ok


def write_csv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def write_matrix_csv(path, matrix, row_edges=None, col_edges=None):
    matrix = np.asarray(matrix, dtype=float)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        if col_edges is not None:
            header = ["truth/reco"] + [f"{col_edges[i]:g}-{col_edges[i+1]:g}" for i in range(len(col_edges)-1)]
        else:
            header = ["row"] + [f"col{j}" for j in range(matrix.shape[1])]
        writer.writerow(header)
        for i in range(matrix.shape[0]):
            if row_edges is not None:
                label = f"{row_edges[i]:g}-{row_edges[i+1]:g}"
            else:
                label = f"row{i}"
            writer.writerow([label] + list(matrix[i]))


def read_zeus_csv(path):
    if not path:
        return None
    path = os.path.expanduser(path)
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            out = {}
            for k, v in r.items():
                try:
                    out[k] = float(v)
                except Exception:
                    out[k] = v
            rows.append(out)
    return rows


def zeus_centers(rows):
    if not rows:
        return np.array([])
    if "t_low" in rows[0] and "t_high" in rows[0]:
        return np.asarray([0.5 * (r["t_low"] + r["t_high"]) for r in rows], dtype=float)
    if "t_center" in rows[0]:
        return np.asarray([r["t_center"] for r in rows], dtype=float)
    return np.asarray([r["t_GeV2"] for r in rows], dtype=float)


def zeus_xerr(rows):
    if not rows:
        return None
    if "t_low" in rows[0] and "t_high" in rows[0]:
        return np.asarray([0.5 * (r["t_high"] - r["t_low"]) for r in rows], dtype=float)
    return None


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def make_reco_t(data, method):
    if method == "pT2":
        return np.asarray(data["reco_pT"], dtype=float) ** 2
    return np.asarray(data[METHODS[method]], dtype=float)


def finite_mask(data, columns):
    mask = np.ones(len(next(iter(data.values()))), dtype=bool)
    for c in columns:
        mask &= np.isfinite(data[c])
    return mask


def response_matrix(true_x, reco_x, weights, edges, reco_selected):
    nb = len(edges) - 1
    truth_hist, truth_err = hist_weighted(true_x, weights, edges)

    H = np.zeros((nb, nb), dtype=float)
    ib = np.digitize(true_x, edges) - 1
    jb = np.digitize(reco_x, edges) - 1
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


def toy_uncertainty(data_reco, data_err, A, prior_counts, n_iterations, n_toys, seed):
    if n_toys <= 0:
        return np.zeros(A.shape[0])
    rng = np.random.default_rng(seed)
    toys = []
    for _ in range(n_toys):
        varied = np.clip(rng.normal(data_reco, data_err), 0.0, None)
        toys.append(iterative_bayes(varied, A, prior_counts, n_iterations))
    toys = np.vstack(toys)
    return 0.5 * (np.nanquantile(toys, 0.84, axis=0) - np.nanquantile(toys, 0.16, axis=0))


def weighted_quantile(values, quantiles, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values = values[mask]
    weights = weights[mask]
    if len(values) == 0:
        return np.full(len(quantiles), np.nan)
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return np.interp(quantiles, cdf, values)


def load_sample(path, label, visible_xsec_fb, args, electron_energy):
    props = dict(UPSILON[args.state], br_ee=branching_fractions_from_args(args)[args.state])
    has_q2 = branch_exists(path, args.tree, "true_Q2")

    columns = ["weight", "true_mass", "true_t", "true_W", "reco_mass", "reco_W", "sigma_epz", "reco_pT"]
    if args.method != "pT2":
        columns.append(METHODS[args.method])
    if has_q2:
        columns.append("true_Q2")

    finite_columns = list(columns)
    columns += proton_branches()
    rdf = ROOT.RDataFrame(args.tree, os.path.expanduser(path))
    arrs = rdf.AsNumpy(columns=list(dict.fromkeys(columns)))
    data = {k: np.asarray(v, dtype=float) for k, v in arrs.items()}

    proton_tag_mask(data)  # Validate all tags before any legacy finite-row filtering.
    mask = finite_mask(data, finite_columns)
    data = {k: v[mask] for k, v in data.items()}

    mass = props["mass"]
    state_mask = (data["true_mass"] >= mass - args.true_mass_window) & (data["true_mass"] <= mass + args.true_mass_window)
    if args.q2_cut is not None and has_q2:
        state_mask &= (data["true_Q2"] >= 0.0) & (data["true_Q2"] <= args.q2_cut)

    if args.tuple_lumi_invfb <= 0:
        raise InputValidationError("--tuple-lumi-invfb must be positive")
    if visible_xsec_fb <= 0:
        raise InputValidationError(f"{label}: target visible cross section must be positive")

    sumw_state_raw = float(np.sum(data["weight"][state_mask]))
    if sumw_state_raw <= 0:
        raise InputValidationError(
            f"{label}: no positive-weight Upsilon({args.state}) truth sample is available; "
            "refusing to generate zero-valued unfolding outputs"
        )

    inferred_visible_xsec_fb = sumw_state_raw / args.tuple_lumi_invfb
    scale = visible_xsec_fb / inferred_visible_xsec_fb

    data["w_visible"] = data["weight"] * scale
    data["reco_t_use"] = make_reco_t(data, args.method)

    truth_common = (
        state_mask
        & (data["true_W"] >= args.w_edges[0])
        & (data["true_W"] <= args.w_edges[-1])
        & (data["true_t"] >= args.t_edges[0])
        & (data["true_t"] <= args.t_edges[-1])
    )

    if not np.any(truth_common):
        raise InputValidationError(
            f"{label}: no Upsilon({args.state}) events enter the requested common-W and |t| phase space"
        )

    sigma_center = 2.0 * electron_energy
    sigma_half = args.sigma_epz_halfwidth if args.sigma_epz_halfwidth is not None else 0.2 * electron_energy
    reco_mass_window = args.mass_window_sigma * props["sigma"]

    reco_selected_all = (
        truth_common & proton_tag_mask(data)
        & (data["reco_mass"] >= mass - reco_mass_window)
        & (data["reco_mass"] <= mass + reco_mass_window)
        & (np.abs(data["sigma_epz"] - sigma_center) <= sigma_half)
        & (data["reco_W"] >= args.w_edges[0])
        & (data["reco_W"] <= args.w_edges[-1])
        & (data["reco_t_use"] >= args.t_edges[0])
        & (data["reco_t_use"] <= args.t_edges[-1])
    )

    if not np.any(reco_selected_all):
        raise InputValidationError(
            f"{label}: the requested common-W selection has zero reconstructed events; "
            "the response matrix would have zero efficiency in every bin"
        )

    sample = {k: v[truth_common] for k, v in data.items()}
    reco_selected = reco_selected_all[truth_common]

    meta = {
        "proton_selection": proton_selection_metadata(),
        "sample": label,
        "br_ee": props["br_ee"],
        "branching_conversion_applied": False,
        **branching_fraction_metadata(args),
        "entries_all_after_finite": len(data["weight"]),
        "entries_truth_common": len(sample["weight"]),
        "sumw_state_raw": sumw_state_raw,
        "inferred_visible_xsec_fb": inferred_visible_xsec_fb,
        "target_visible_xsec_fb": visible_xsec_fb,
        "state_weight_scale": scale,
        "has_true_Q2": int(has_q2),
    }
    return sample, reco_selected, meta


def build_common_w_weights(low, high, args):
    low_hist, _ = hist_weighted(low["true_W"], low["w_visible"], args.w_edges)
    high_hist, _ = hist_weighted(high["true_W"], high["w_visible"], args.w_edges)

    if args.target_w_mode == "low":
        target = low_hist.copy()
    elif args.target_w_mode == "high":
        target = high_hist.copy()
    elif args.target_w_mode == "flat":
        target = np.ones_like(low_hist) * 0.5 * (low_hist.sum() + high_hist.sum()) / len(low_hist)
    else:
        target = 0.5 * (low_hist + high_hist)

    def factors(sample, hist):
        ib = np.digitize(sample["true_W"], args.w_edges) - 1
        out = np.zeros(len(sample["true_W"]), dtype=float)
        ok = (ib >= 0) & (ib < len(hist)) & (hist[ib] > 0)
        out[ok] = target[ib[ok]] / hist[ib[ok]]
        return out

    rows = []
    low_fac_bins = np.divide(target, low_hist, out=np.zeros_like(target), where=low_hist > 0)
    high_fac_bins = np.divide(target, high_hist, out=np.zeros_like(target), where=high_hist > 0)
    for i in range(len(target)):
        rows.append({
            "W_low": args.w_edges[i],
            "W_high": args.w_edges[i+1],
            "low_visible_yield_before": low_hist[i],
            "high_visible_yield_before": high_hist[i],
            "target_common_visible_yield": target[i],
            "low_to_target_factor": low_fac_bins[i],
            "high_to_target_factor": high_fac_bins[i],
        })
    return factors(low, low_hist), factors(high, high_hist), rows


def unfold_sample(sample, reco_selected, args, seed):
    rng = np.random.default_rng(seed)
    train = rng.random(len(sample["true_t"])) < args.train_fraction
    test = ~train

    H, truth_train, truth_train_err, A, eff = response_matrix(
        sample["true_t"][train],
        sample["reco_t_use"][train],
        sample["w_common"][train],
        args.t_edges,
        reco_selected[train],
    )

    reco_test, reco_test_err = hist_weighted(
        sample["reco_t_use"][test & reco_selected],
        sample["w_common"][test & reco_selected],
        args.t_edges,
    )
    truth_test, truth_test_err = hist_weighted(sample["true_t"][test], sample["w_common"][test], args.t_edges)

    unfolded = iterative_bayes(reco_test, A, truth_train, args.n_iterations)
    unfolded_err = toy_uncertainty(reco_test, reco_test_err, A, truth_train, args.n_iterations, args.n_toys, seed + 1000)

    purity, stability = purity_stability(H)

    widths = np.diff(args.t_edges)
    centers = 0.5 * (args.t_edges[:-1] + args.t_edges[1:])
    closure = np.divide(unfolded, truth_test, out=np.full_like(unfolded, np.nan), where=truth_test > 0)
    closure_err = np.divide(unfolded_err, truth_test, out=np.full_like(unfolded_err, np.nan), where=truth_test > 0)

    rows = []
    for i in range(len(widths)):
        rows.append({
            "t_low": args.t_edges[i],
            "t_high": args.t_edges[i+1],
            "t_center": centers[i],
            "bin_width": widths[i],
            "efficiency": eff[i],
            "purity": purity[i],
            "stability": stability[i],
            "truth_test_visible_commonW": truth_test[i],
            "reco_test_visible_commonW": reco_test[i],
            "unfolded_visible_commonW": unfolded[i],
            "unfolded_visible_commonW_err": unfolded_err[i],
            "dNdt_unfolded_visible_commonW": unfolded[i] / widths[i],
            "dNdt_unfolded_visible_commonW_err": unfolded_err[i] / widths[i],
            "closure_unfold_over_truth": closure[i],
            "closure_err": closure_err[i],
        })

    # Resolution from selected events.
    res_rows = []
    for i in range(len(widths)):
        b = (sample["true_t"] >= args.t_edges[i]) & (sample["true_t"] < args.t_edges[i+1]) & reco_selected
        if np.any(b):
            delta = sample["reco_t_use"][b] - sample["true_t"][b]
            rel = delta / sample["true_t"][b]
            ww = sample["w_common"][b]
            q16, q50, q84 = weighted_quantile(delta, [0.16, 0.5, 0.84], ww)
            rq16, rq50, rq84 = weighted_quantile(rel, [0.16, 0.5, 0.84], ww)
            res_rows.append({
                "t_low": args.t_edges[i],
                "t_high": args.t_edges[i+1],
                "t_center": centers[i],
                "median_delta_t": q50,
                "core_sigma_t": 0.5 * (q84 - q16),
                "median_rel_delta_t": rq50,
                "core_sigma_rel": 0.5 * (rq84 - rq16),
                "n_selected": int(np.sum(b)),
                "sumw_selected": float(np.sum(ww)),
            })
        else:
            res_rows.append({
                "t_low": args.t_edges[i],
                "t_high": args.t_edges[i+1],
                "t_center": centers[i],
                "median_delta_t": np.nan,
                "core_sigma_t": np.nan,
                "median_rel_delta_t": np.nan,
                "core_sigma_rel": np.nan,
                "n_selected": 0,
                "sumw_selected": 0.0,
            })

    return rows, res_rows, H, A


def scale_model_to_zeus(x, y, rows):
    if not rows:
        return 1.0
    if "t_low" in rows[0] and "t_high" in rows[0]:
        zx = np.asarray([0.5 * (r["t_low"] + r["t_high"]) for r in rows])
    elif "t_center" in rows[0]:
        zx = np.asarray([r["t_center"] for r in rows])
    else:
        zx = np.asarray([r["t_GeV2"] for r in rows])

    num = 0.0
    den = 0.0
    for k, r in enumerate(rows):
        idx = int(np.argmin(np.abs(x - zx[k])))
        m = y[idx]
        if not np.isfinite(m) or m <= 0:
            continue
        err = 0.5 * (float(r.get("err_low", 0.0)) + float(r.get("err_high", 0.0)))
        if err <= 0:
            err = max(float(r["dNdt_GeVminus2"]) * 0.2, 1.0)
        ww = 1.0 / (err * err)
        num += float(r["dNdt_GeVminus2"]) * m * ww
        den += m * m * ww
    return num / den if den > 0 else 1.0


def plot_hybrid(hybrid_rows, low_rows, high_rows, zeus_rows, args):
    if zeus_rows and args.state != "1S":
        raise ValueError("The ZEUS reference measures Upsilon(1S), not the requested excited state")
    fig, ax = plt.subplots(figsize=(8.5, 6.8))

    for rows, label, color, marker in [
        (low_rows, args.low_label + " all bins", "tab:blue", "o"),
        (high_rows, args.high_label + " all bins", "tab:orange", "s"),
    ]:
        x = np.asarray([r["t_center"] for r in rows])
        xerr = 0.5 * np.asarray([r["t_high"] - r["t_low"] for r in rows])
        y = np.asarray([r["dNdt_unfolded_visible_commonW"] for r in rows])
        yerr = np.asarray([r["dNdt_unfolded_visible_commonW_err"] for r in rows])
        m = np.isfinite(y) & (y > 0)
        ax.errorbar(x[m], y[m], xerr=xerr[m], yerr=yerr[m], fmt=marker, color=color, alpha=0.25, capsize=2, label=label)

    x = np.asarray([r["t_center"] for r in hybrid_rows])
    xerr = 0.5 * np.asarray([r["t_high"] - r["t_low"] for r in hybrid_rows])
    y = np.asarray([r["dNdt_unfolded_visible_commonW"] for r in hybrid_rows])
    yerr = np.asarray([r["dNdt_unfolded_visible_commonW_err"] for r in hybrid_rows])

    if zeus_rows and args.zeus_mode == "scale_to_zeus":
        scale = scale_model_to_zeus(x, y, zeus_rows)
        y_plot = y * scale
        yerr_plot = yerr * scale
        ylabel = r"scaled $dN/d|t|$ [arb.]"
        label = rf"hybrid common-$W$ shape $\times {scale:.2g}$"
    else:
        y_plot = y
        yerr_plot = yerr
        ylabel = r"$dN_{ee}/d|t|$ [GeV$^{-2}$]"
        label = r"hybrid common-$W$"

    m = np.isfinite(y_plot) & (y_plot > 0)
    ax.errorbar(x[m], y_plot[m], xerr=xerr[m], yerr=yerr_plot[m], fmt="D-", color="black", capsize=3, label=label, zorder=5)

    if zeus_rows:
        if "t_low" in zeus_rows[0] and "t_high" in zeus_rows[0]:
            step_edges = np.asarray([r["t_low"] for r in zeus_rows] + [zeus_rows[-1]["t_high"]])
            vals = np.asarray([r["dNdt_GeVminus2"] for r in zeus_rows])
            step_vals = np.r_[vals, vals[-1]]
            ax.step(step_edges, step_vals, where="post", color="black", linestyle=":", linewidth=2.2,
                    label=r"ZEUS exclusive $\Upsilon(1S)$ component")
        else:
            zx = np.asarray([r.get("t_GeV2", r.get("t_center")) for r in zeus_rows])
            zy = np.asarray([r["dNdt_GeVminus2"] for r in zeus_rows])
            ax.plot(zx, zy, "k:", label=r"ZEUS exclusive $\Upsilon(1S)$ component")

    ax.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.set_xlim(args.t_edges[0], args.t_edges[-1])
    ax.set_ylim(args.y_min, args.y_max)
    ax.set_title(rf"$\Upsilon({args.state})$: hybrid $|t|$ unfolding after common-$W_{{\gamma p}}$ reweighting")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out_prefix + "_hybrid_commonW_overlay.png", dpi=250)
    fig.savefig(args.out_prefix + "_hybrid_commonW_overlay.pdf")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--low-file", required=True)
    p.add_argument("--high-file", required=True)
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--state", default="1S", choices=list(UPSILON.keys()))
    p.add_argument("--method", default="pT2", choices=list(METHODS.keys()))
    p.add_argument("--low-visible-xsec-fb", type=float, required=True)
    p.add_argument("--high-visible-xsec-fb", type=float, required=True)
    p.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    p.add_argument("--t-edges", required=True)
    p.add_argument("--switch-t", type=float, default=0.5)
    p.add_argument("--w-edges", required=True)
    p.add_argument("--target-w-mode", choices=["average", "low", "high", "flat"], default="average")
    p.add_argument("--q2-cut", type=float, default=None)
    p.add_argument("--electron-energy-low", type=float, default=10.0)
    p.add_argument("--electron-energy-high", type=float, default=18.0)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    p.add_argument("--mass-window-sigma", type=float, default=3.0)
    p.add_argument("--true-mass-window", type=float, default=0.2)
    p.add_argument("--train-fraction", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--n-iterations", type=int, default=4)
    p.add_argument("--n-toys", type=int, default=300)
    p.add_argument("--low-label", default=r"10$\times$100 GeV")
    p.add_argument("--high-label", default=r"18$\times$275 GeV")
    p.add_argument("--zeus-csv", default="")
    p.add_argument("--zeus-mode", choices=["scale_to_zeus", "none"], default="scale_to_zeus")
    p.add_argument("--out-prefix", default="hybrid_t_commonW")
    p.add_argument("--y-min", type=float, default=0.1)
    p.add_argument("--y-max", type=float, default=1000.0)
    add_branching_fraction_argument(p)
    args = p.parse_args()

    args.t_edges = parse_edges(args.t_edges)
    args.w_edges = parse_edges(args.w_edges)

    low, low_sel, low_meta = load_sample(args.low_file, "low", args.low_visible_xsec_fb, args, args.electron_energy_low)
    high, high_sel, high_meta = load_sample(args.high_file, "high", args.high_visible_xsec_fb, args, args.electron_energy_high)

    low_factor, high_factor, w_rows = build_common_w_weights(low, high, args)
    low["w_common"] = low["w_visible"] * low_factor
    high["w_common"] = high["w_visible"] * high_factor

    low_rows, low_res, low_H, low_A = unfold_sample(low, low_sel, args, args.seed + 10)
    high_rows, high_res, high_H, high_A = unfold_sample(high, high_sel, args, args.seed + 20)

    for r in low_rows:
        r["sample"] = "low"
        r["sample_label"] = args.low_label
    for r in high_rows:
        r["sample"] = "high"
        r["sample_label"] = args.high_label

    low_pick = [r for r in low_rows if r["t_high"] <= args.switch_t]
    if not low_pick and low_rows:
        low_pick = [sorted(low_rows, key=lambda x: x["t_center"])[0]]
    high_pick = [r for r in high_rows if r["t_low"] >= args.switch_t]
    hybrid_rows = sorted(low_pick + high_pick, key=lambda x: x["t_center"])

    meta_rows = [low_meta, high_meta]
    for r in meta_rows:
        r["target_w_mode"] = args.target_w_mode
        r["common_W_low"] = args.w_edges[0]
        r["common_W_high"] = args.w_edges[-1]

    prefix = args.out_prefix
    write_csv(prefix + "_normalization_meta.csv", meta_rows)
    write_csv(prefix + "_W_reweighting_table.csv", w_rows)
    write_csv(prefix + "_low_t_unfold_commonW.csv", low_rows)
    write_csv(prefix + "_high_t_unfold_commonW.csv", high_rows)
    write_csv(prefix + "_low_t_resolution_commonW.csv", low_res)
    write_csv(prefix + "_high_t_resolution_commonW.csv", high_res)
    write_csv(prefix + "_hybrid_commonW.csv", hybrid_rows)

    write_matrix_csv(prefix + "_low_response_counts_commonW.csv", low_H, args.t_edges, args.t_edges)
    write_matrix_csv(prefix + "_high_response_counts_commonW.csv", high_H, args.t_edges, args.t_edges)
    write_matrix_csv(prefix + "_low_response_prob_commonW.csv", low_A, args.t_edges, args.t_edges)
    write_matrix_csv(prefix + "_high_response_prob_commonW.csv", high_A, args.t_edges, args.t_edges)

    use_zeus_reference = bool(args.zeus_csv and args.state == "1S")
    if args.zeus_csv and not use_zeus_reference:
        print(f"Skipping ZEUS Upsilon(1S) reference for Upsilon({args.state}).")
    zeus_rows = read_zeus_csv(args.zeus_csv) if use_zeus_reference else None
    plot_hybrid(hybrid_rows, low_rows, high_rows, zeus_rows, args)

    print("=== W reweighting table ===")
    for r in w_rows:
        print(r)

    print("\n=== Hybrid common-W bins ===")
    for r in hybrid_rows:
        print({k: r[k] for k in [
            "sample_label", "t_low", "t_high", "t_center",
            "efficiency", "purity", "stability", "closure_unfold_over_truth",
            "dNdt_unfolded_visible_commonW", "dNdt_unfolded_visible_commonW_err"
        ]})

    print(f"\nWrote {prefix}_hybrid_commonW.csv")
    print(f"Wrote {prefix}_hybrid_commonW_overlay.png")
    print(f"Wrote {prefix}_W_reweighting_table.csv")


if __name__ == "__main__":
    main()
