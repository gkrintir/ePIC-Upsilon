#!/usr/bin/env python3
"""
Extract sigma(gamma p -> Upsilon(1S) p) versus W with robust tuple-weight renormalization.

This version addresses the common ePIC/SWAN pitfall where separate output tuples are produced
per shard/run with the full generator cross section, and then merged. In that case the merged
sum of weights is too large by approximately the number of shards.

Use:
  --renormalize-sumw-to-xsec 119.685

This rescales all tuple weights so that sum(weight)/tuple_lumi_invfb equals 119.685 fb.
"""

import argparse
from pathlib import Path as _BranchingPath
import sys as _branching_sys
_branching_sys.path.insert(0, str(_BranchingPath(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument, branching_fractions_from_args,
                                branching_fraction_metadata)
import math
import os
import sys

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
M_UPSILON_1S = 9.46030
SIGMA_MASS_1S = 0.053
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
        ((1.0 + (1.0 - y) ** 2) / y) * log_term
        - (2.0 * (1.0 - y) / y) * (1.0 - q2_min / q2_max)
    )
    return dndy * (2.0 * W / s_minus_mp2)


def integrated_flux(w_low, w_high, e_e, e_p, q2_max, n_steps=500):
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    return float(np.trapezoid(vals, grid))


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def load_theory_csv(path):
    if not path:
        return None
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    raw = pd.read_csv(path)
    low_cols = [c for c in ["muR_min", "muF_min", "PDFerr_min"] if c in raw.columns]
    high_cols = [c for c in ["muR_max", "muF_max", "PDFerr_max"] if c in raw.columns]

    if "s_gp" not in raw.columns or not low_cols or not high_cols:
        raise ValueError("Theory CSV must contain s_gp and min/max columns")

    low_nb = raw[low_cols].min(axis=1).to_numpy()
    high_nb = raw[high_cols].max(axis=1).to_numpy()
    mid_nb = 0.5 * (low_nb + high_nb)

    return pd.DataFrame({
        "W": raw["s_gp"].to_numpy(dtype=float),
        "low_pb": 1000.0 * low_nb,
        "mid_pb": 1000.0 * mid_nb,
        "high_pb": 1000.0 * high_nb,
    }).sort_values("W")


def theory_interp(theory, W):
    if theory is None:
        return np.full_like(W, np.nan, dtype=float), np.full_like(W, np.nan, dtype=float), np.full_like(W, np.nan, dtype=float)
    Wtab = theory["W"].to_numpy()
    return (
        np.interp(W, Wtab, theory["low_pb"].to_numpy(), left=np.nan, right=np.nan),
        np.interp(W, Wtab, theory["mid_pb"].to_numpy(), left=np.nan, right=np.nan),
        np.interp(W, Wtab, theory["high_pb"].to_numpy(), left=np.nan, right=np.nan),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--out-prefix", default="sigma_gp_1S_10x100_fixed")
    p.add_argument("--theory-csv", default="")
    p.add_argument("--electron-energy", type=float, default=10.0)
    p.add_argument("--proton-energy", type=float, default=100.0)
    p.add_argument("--q2-max", type=float, default=0.01)
    p.add_argument("--w-min", type=float, default=10.0)
    p.add_argument("--w-max", type=float, default=50.0)
    p.add_argument("--n-w-bins", type=int, default=14)
    p.add_argument("--lumi-invfb", type=float, default=10.0)
    p.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    p.add_argument("--expected-generator-xsec-fb", type=float, default=119.685)
    p.add_argument("--renormalize-sumw-to-xsec", type=float, default=0.0,
                   help="If >0, multiply all weights so sum(weight)/tuple_lumi_invfb equals this fb value.")
    p.add_argument("--renormalize-sumw-scope", choices=["all", "truth_1s"], default="truth_1s",
                   help="Use all tuple entries or only true-mass 1S entries to determine the weight scale. Use truth_1s for merged 1S+2S+3S tuples.")
    p.add_argument("--weight-convention", choices=["inclusive_upsilon", "visible_decay"], default="inclusive_upsilon")
    p.add_argument("--mass-window-sigma", type=float, default=3.0)
    p.add_argument("--true-mass-window", type=float, default=0.2)
    p.add_argument("--sigma-epz-center", type=float, default=None)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    add_branching_fraction_argument(p)
    args = p.parse_args()

    infile = os.path.expanduser(args.input)
    df = ROOT.RDataFrame(args.tree, infile)

    # Decide the scope used to determine the normalization scale.
    # For merged 1S+2S+3S tuples, normalizing with the sum of all weights would
    # under-normalize the selected 1S subset. Use the 1S truth-mass window instead.
    mass_window_norm = args.true_mass_window
    true_mass_min_norm = M_UPSILON_1S - mass_window_norm
    true_mass_max_norm = M_UPSILON_1S + mass_window_norm

    sumw_all = float(df.Sum("weight").GetValue())
    if args.renormalize_sumw_scope == "truth_1s":
        norm_filter = f"true_mass >= {true_mass_min_norm} && true_mass <= {true_mass_max_norm}"
        sumw_norm = float(df.Filter(norm_filter).Sum("weight").GetValue())
        norm_scope_label = f"truth_1s [{true_mass_min_norm:.3f},{true_mass_max_norm:.3f}]"
    else:
        sumw_norm = sumw_all
        norm_scope_label = "all"

    inferred_xsec = sumw_norm / args.tuple_lumi_invfb

    weight_scale = 1.0
    if args.renormalize_sumw_to_xsec > 0:
        weight_scale = args.renormalize_sumw_to_xsec / inferred_xsec

    print("=== Normalization ===")
    print(f"sum weight all                       = {sumw_all:.8g}")
    print(f"normalization scope                  = {norm_scope_label}")
    print(f"sum weight in norm scope             = {sumw_norm:.8g}")
    print(f"inferred xsec in norm scope          = {inferred_xsec:.8g} fb")
    print(f"expected generator xsec              = {args.expected_generator_xsec_fb:.8g} fb")
    print(f"renormalization target               = {args.renormalize_sumw_to_xsec:.8g} fb")
    print(f"weight scale applied                 = {weight_scale:.8g}")

    if args.renormalize_sumw_to_xsec <= 0 and args.expected_generator_xsec_fb > 0:
        ratio = inferred_xsec / args.expected_generator_xsec_fb
        if abs(ratio - 1.0) > 0.05:
            print(f"WARNING: inferred/expected = {ratio:.4g}. Use --renormalize-sumw-to-xsec {args.expected_generator_xsec_fb:g} if this is unintended.")

    sigma_epz_center = 2.0 * args.electron_energy if args.sigma_epz_center is None else args.sigma_epz_center
    sigma_epz_halfwidth = 0.2 * args.electron_energy if args.sigma_epz_halfwidth is None else args.sigma_epz_halfwidth

    w_bins = np.linspace(args.w_min, args.w_max, args.n_w_bins + 1)
    w_centers = 0.5 * (w_bins[:-1] + w_bins[1:])

    mass_window = args.mass_window_sigma * SIGMA_MASS_1S
    reco_mass_min = M_UPSILON_1S - mass_window
    reco_mass_max = M_UPSILON_1S + mass_window
    true_mass_min = M_UPSILON_1S - args.true_mass_window
    true_mass_max = M_UPSILON_1S + args.true_mass_window

    reco_filter = (
        f"reco_mass >= {reco_mass_min} && reco_mass <= {reco_mass_max}"
        f" && abs(sigma_epz - {sigma_epz_center}) <= {sigma_epz_halfwidth}"
        f" && reco_W >= {args.w_min} && reco_W <= {args.w_max}"
    )
    true_filter = (
        f"true_mass >= {true_mass_min} && true_mass <= {true_mass_max}"
        f" && true_W >= {args.w_min} && true_W <= {args.w_max}"
    )

    reco = df.Filter(reco_filter).AsNumpy(["reco_W", "weight"])
    truth = df.Filter(true_filter).AsNumpy(["true_W", "weight"])

    reco_w = np.asarray(reco["weight"], dtype=float) * weight_scale
    true_w = np.asarray(truth["weight"], dtype=float) * weight_scale

    n_reco, err_reco = hist_weighted(np.asarray(reco["reco_W"], dtype=float), reco_w, w_bins)
    n_true, err_true = hist_weighted(np.asarray(truth["true_W"], dtype=float), true_w, w_bins)

    eff = np.divide(n_reco, n_true, out=np.zeros_like(n_reco), where=n_true > 0)

    flux_int = np.array([integrated_flux(w_bins[i], w_bins[i+1], args.electron_energy, args.proton_energy, args.q2_max)
                         for i in range(len(w_centers))])

    br_corr = 1.0 / branching_fractions_from_args(args)["1S"] if args.weight_convention == "visible_decay" else 1.0

    denom = args.lumi_invfb * 1000.0 * flux_int
    sigma_truth = np.divide(n_true, denom, out=np.full_like(n_true, np.nan), where=denom > 0) * br_corr
    sigma_reco = np.divide(n_reco, denom * eff, out=np.full_like(n_reco, np.nan), where=(denom > 0) & (eff > 0)) * br_corr
    err_sigma = np.divide(err_reco, denom * eff, out=np.full_like(err_reco, np.nan), where=(denom > 0) & (eff > 0)) * br_corr

    theory = load_theory_csv(args.theory_csv) if args.theory_csv else None
    th_low, th_mid, th_high = theory_interp(theory, w_centers)

    out = pd.DataFrame({
        "W_low": w_bins[:-1],
        "W_high": w_bins[1:],
        "W_center": w_centers,
        "flux_int": flux_int,
        "n_true_weighted": n_true,
        "n_reco_weighted": n_reco,
        "eff_reco_over_true": eff,
        "sigma_gp_from_truth_pb": sigma_truth,
        "sigma_gp_from_reco_pb": sigma_reco,
        "sigma_gp_stat_err_pb": err_sigma,
        "theory_low_pb": th_low,
        "theory_mid_pb": th_mid,
        "theory_high_pb": th_high,
        "ratio_to_theory_mid": sigma_truth / th_mid,
    })

    csv_path = os.path.expanduser(args.out_prefix + ".csv")
    png_path = os.path.expanduser(args.out_prefix + ".png")
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    out["br_ee"] = branching_fractions_from_args(args)["1S"]
    out["branching_fraction_convention"] = branching_fraction_metadata(args)["branching_fraction_convention"]
    out["branching_fraction_uncertainty_included"] = False
    out.to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(9, 7))

    if theory is not None:
        ax.fill_between(theory["W"], theory["low_pb"], theory["high_pb"], color="tab:red", alpha=0.18, linewidth=0, label="CT18NLO Y(1S)")
        ax.plot(theory["W"], theory["mid_pb"], color="tab:red", linewidth=2)

    m = np.isfinite(sigma_reco) & (sigma_reco > 0)
    ax.errorbar(w_centers[m], sigma_reco[m], yerr=err_sigma[m], fmt="o", color="royalblue", label="ePIC projected Y(1S)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(args.w_min, max(args.w_max, 60))
    ax.set_ylim(0.01, 200)
    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]")
    ax.set_ylabel(r"$\sigma(\gamma p\rightarrow \Upsilon(1S)p)$ [pb]")
    ax.text(0.04, 0.94, "ePIC Performance", transform=ax.transAxes, fontweight="bold", va="top")
    ax.text(0.04, 0.88, rf"$e+p$, {args.electron_energy:g}$\times${args.proton_energy:g} GeV", transform=ax.transAxes, va="top")
    ax.text(0.04, 0.82, rf"$Q^2 < {args.q2_max:g}$ GeV$^2$", transform=ax.transAxes, va="top")
    ax.text(0.04, 0.76, rf"$L_{{ep}} = {args.lumi_invfb:g}$ fb$^{{-1}}$", transform=ax.transAxes, va="top")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(png_path, dpi=250)

    print("\n=== Extracted bins ===")
    print(out[["W_low", "W_high", "n_true_weighted", "n_reco_weighted", "eff_reco_over_true", "sigma_gp_from_truth_pb", "theory_mid_pb", "ratio_to_theory_mid"]].to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print(f"\nWrote {csv_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
