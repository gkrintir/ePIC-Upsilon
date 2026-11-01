#!/usr/bin/env python3
"""
Compare generator-level extracted sigma(gamma p -> Upsilon(1S) p)
from two eSTARlight/ePIC UpsilonTree samples.

Purpose:
  Diagnose discontinuities between lower- and higher-sqrt(s_ep) samples.
  Uses truth-level W only, so detector/reconstruction does not enter.

Key convention:
  Provide the visible e+e- generated-sample cross section from the eSTARlight log
  via --xsec-visible-fb. The script rescales the merged tuple weights so that
  sum(weight)/tuple_lumi equals that visible cross section, then unfolds BR_ee
  to obtain inclusive Upsilon(1S) gamma-p cross sections.

Example:
  python3 compare_sigma_gp_two_energies_swan.py \
    --low-file /eos/user/g/gkrintir/merged_elecs_test_10t100.root \
    --low-ee 10 --low-ep 100 --low-xsec-visible-fb 104.793 \
    --high-file /eos/user/g/gkrintir/merged_elecs_test_tuples.root \
    --high-ee 18 --high-ep 275 --high-xsec-visible-fb 327.584 \
    --w-min 10 --w-max 63 --n-w-bins 12 \
    --q2-max-flux 1.0 \
    --q2-cut 1.0 \
    --out-prefix compare_10x100_18x275
"""

import argparse
import math
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import ROOT
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument,
                                branching_fractions_from_args,
                                branching_fraction_metadata)

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass

ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813
M_UPSILON_1S = 9.46030
BR_EE_1S = 0.0234


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
        ((1.0 + (1.0 - y)**2) / y) * log_term
        - (2.0 * (1.0 - y) / y) * (1.0 - q2_min / q2_max)
    )
    return dndy * (2.0 * W / s_minus_mp2)


def integrated_flux(w_low, w_high, e_e, e_p, q2_max, n_steps=500):
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    return float(np.trapezoid(vals, grid))


def branch_exists(tree, branch):
    return tree.GetBranch(branch) is not None


def load_truth_sample(path, tree_name, target_visible_xsec_fb, tuple_lumi_invfb, mass_half_width, q2_cut,
                      branching_fraction=BR_EE_1S):
    path = os.path.expanduser(path)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {path}")
    t = f.Get(tree_name)
    if not t:
        raise RuntimeError(f"Could not find tree {tree_name} in {path}")

    has_true_q2 = branch_exists(t, "true_Q2")
    f.Close()

    cols = ["true_W", "true_mass", "weight"]
    if has_true_q2:
        cols.append("true_Q2")

    rdf = ROOT.RDataFrame(tree_name, path)
    arr = rdf.AsNumpy(columns=cols)

    df = pd.DataFrame({k: np.asarray(v, dtype=float) for k, v in arr.items()})
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    sumw_all = df["weight"].sum()

    # Normalize to the raw weight sum for the requested 1S truth-mass subset,
    # not to the whole file. This matters for merged 1S+2S+3S tuples.
    state_sel_for_norm = (df["true_mass"] > M_UPSILON_1S - mass_half_width) & (df["true_mass"] < M_UPSILON_1S + mass_half_width)
    sumw_state_raw = df.loc[state_sel_for_norm, "weight"].sum()
    inferred_visible_xsec_fb = sumw_state_raw / tuple_lumi_invfb
    scale = target_visible_xsec_fb / inferred_visible_xsec_fb if inferred_visible_xsec_fb > 0 else np.nan

    # Visible-normalized, then BR-unfolded to inclusive Upsilon.
    df["w_visible_corr"] = df["weight"] * scale
    df["w_inclusive_corr"] = df["w_visible_corr"] / branching_fraction

    sel = state_sel_for_norm.copy()
    if q2_cut is not None and has_true_q2:
        sel = sel & (df["true_Q2"] >= 0.0) & (df["true_Q2"] <= q2_cut)

    out = df[sel].copy()
    meta = {
        "path": path,
        "entries_all": len(df),
        "entries_selected": len(out),
        "has_true_Q2": has_true_q2,
        "sumw_all_raw": float(sumw_all),
        "sumw_state_raw": float(sumw_state_raw),
        "tuple_lumi_invfb": tuple_lumi_invfb,
        "inferred_visible_xsec_fb_raw_state": float(inferred_visible_xsec_fb),
        "target_visible_xsec_fb": float(target_visible_xsec_fb),
        "weight_scale_to_visible": float(scale),
        "br_ee": branching_fraction,
        "target_inclusive_xsec_pb": float(target_visible_xsec_fb / branching_fraction / 1000.0),
    }
    return out, meta


def extract_sigma(df, e_e, e_p, w_edges, q2_max_flux, lumi_invfb):
    w = df["true_W"].to_numpy()
    weights = df["w_inclusive_corr"].to_numpy()

    n, _ = np.histogram(w, bins=w_edges, weights=weights)
    n2, _ = np.histogram(w, bins=w_edges, weights=weights * weights)
    err = np.sqrt(n2)

    flux = np.array([integrated_flux(w_edges[i], w_edges[i+1], e_e, e_p, q2_max_flux) for i in range(len(w_edges)-1)])
    denom = lumi_invfb * 1000.0 * flux

    sigma = np.divide(n, denom, out=np.full_like(n, np.nan, dtype=float), where=denom > 0)
    sigma_err = np.divide(err, denom, out=np.full_like(err, np.nan, dtype=float), where=denom > 0)

    return pd.DataFrame({
        "W_low": w_edges[:-1],
        "W_high": w_edges[1:],
        "W_center": 0.5 * (w_edges[:-1] + w_edges[1:]),
        "flux_int": flux,
        "yield_inclusive_weighted": n,
        "yield_inclusive_weighted_err": err,
        "sigma_gp_pb": sigma,
        "sigma_gp_err_pb": sigma_err,
    })


def main():
    p = argparse.ArgumentParser()
    add_branching_fraction_argument(p)
    p.add_argument("--low-file", required=True)
    p.add_argument("--high-file", required=True)
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--low-ee", type=float, required=True)
    p.add_argument("--low-ep", type=float, required=True)
    p.add_argument("--high-ee", type=float, required=True)
    p.add_argument("--high-ep", type=float, required=True)
    p.add_argument("--low-xsec-visible-fb", type=float, required=True)
    p.add_argument("--high-xsec-visible-fb", type=float, required=True)
    p.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    p.add_argument("--lumi-invfb", type=float, default=10.0)
    p.add_argument("--w-min", type=float, default=10.0)
    p.add_argument("--w-max", type=float, default=63.0)
    p.add_argument("--n-w-bins", type=int, default=12)
    p.add_argument("--q2-max-flux", type=float, default=1.0)
    p.add_argument("--q2-cut", type=float, default=None)
    p.add_argument("--mass-half-width", type=float, default=0.20)
    p.add_argument("--out-prefix", default="compare_sigma_gp_two_energies")
    p.add_argument("--min-yield-for-ratio", type=float, default=1e-9,
                   help="Only plot low/high ratios when both weighted yields exceed this value.")
    args = p.parse_args()
    br_ee = branching_fractions_from_args(args)["1S"]

    w_edges = np.linspace(args.w_min, args.w_max, args.n_w_bins + 1)

    low_df, low_meta = load_truth_sample(
        args.low_file, args.tree, args.low_xsec_visible_fb,
        args.tuple_lumi_invfb, args.mass_half_width, args.q2_cut, branching_fraction=br_ee
    )
    high_df, high_meta = load_truth_sample(
        args.high_file, args.tree, args.high_xsec_visible_fb,
        args.tuple_lumi_invfb, args.mass_half_width, args.q2_cut, branching_fraction=br_ee
    )

    low_meta.update(branching_fraction_metadata(args))
    high_meta.update(branching_fraction_metadata(args))
    low = extract_sigma(low_df, args.low_ee, args.low_ep, w_edges, args.q2_max_flux, args.lumi_invfb)
    high = extract_sigma(high_df, args.high_ee, args.high_ep, w_edges, args.q2_max_flux, args.lumi_invfb)

    merged = low.add_prefix("low_").join(high.add_prefix("high_"))
    merged["W_low"] = low["W_low"]
    merged["W_high"] = low["W_high"]
    merged["W_center"] = low["W_center"]
    with np.errstate(divide="ignore", invalid="ignore"):
        merged["ratio_low_over_high"] = merged["low_sigma_gp_pb"] / merged["high_sigma_gp_pb"]

    # A ratio is meaningful only when both samples have support in the bin.
    good_ratio = (
        np.isfinite(merged["ratio_low_over_high"]) &
        (merged["low_yield_inclusive_weighted"] > args.min_yield_for_ratio) &
        (merged["high_yield_inclusive_weighted"] > args.min_yield_for_ratio)
    )
    merged["ratio_low_over_high_valid"] = good_ratio

    csv = args.out_prefix + ".csv"
    meta_csv = args.out_prefix + "_meta.csv"
    png = args.out_prefix + ".png"

    merged.to_csv(csv, index=False)
    pd.DataFrame([{"sample": "low", **low_meta}, {"sample": "high", **high_meta}]).to_csv(meta_csv, index=False)

    print("=== Normalization metadata ===")
    print(pd.DataFrame([{"sample": "low", **low_meta}, {"sample": "high", **high_meta}]).to_string(index=False))
    print("\n=== Comparison ===")
    cols = ["W_low", "W_high", "low_sigma_gp_pb", "high_sigma_gp_pb", "ratio_low_over_high",
            "low_yield_inclusive_weighted", "high_yield_inclusive_weighted"]
    print(merged[cols].to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})

    lm = np.isfinite(low["sigma_gp_pb"]) & (low["sigma_gp_pb"] > 0)
    hm = np.isfinite(high["sigma_gp_pb"]) & (high["sigma_gp_pb"] > 0)

    ax.errorbar(low["W_center"][lm], low["sigma_gp_pb"][lm], yerr=low["sigma_gp_err_pb"][lm],
                fmt="o", label=f"low {args.low_ee:g}x{args.low_ep:g}")
    ax.errorbar(high["W_center"][hm], high["sigma_gp_pb"][hm], yerr=high["sigma_gp_err_pb"][hm],
                fmt="s", label=f"high {args.high_ee:g}x{args.high_ep:g}")
    ax.set_yscale("log")
    ax.set_ylabel(r"$\sigma_{\gamma p}$ [pb]")
    ax.legend(frameon=False)

    rm = merged["ratio_low_over_high_valid"].to_numpy(dtype=bool)
    rax.axhline(1.0, color="gray", linestyle="--")
    if rm.any():
        rax.plot(merged["W_center"][rm], merged["ratio_low_over_high"][rm], "ko-")
        ymax = max(3.0, float(np.nanmax(merged.loc[rm, "ratio_low_over_high"])) * 1.2)
    else:
        ymax = 3.0
    rax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    rax.set_ylabel("low/high")
    rax.set_ylim(0, ymax)

    fig.tight_layout()
    fig.savefig(png, dpi=250)
    print(f"\nWrote {csv}")
    print(f"Wrote {meta_csv}")
    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
