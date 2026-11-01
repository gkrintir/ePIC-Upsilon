#!/usr/bin/env python3
"""
Debug and extract sigma(gamma p -> Upsilon(1S) p) versus W for an ePIC UpsilonTree.

This is intentionally verbose. It is meant to diagnose normalization mismatches by printing:
  - sum of tuple weights
  - selected and truth-weighted yields
  - flux integral per W bin
  - extracted gamma-p cross section per W bin
  - CT18NLO gamma-p theory in the same bins
  - CT18NLO ep prediction obtained by folding the CT18 gamma-p curve with the same photon flux

Assumptions:
  - The tuple weight is an ep expected-yield weight for tuple_lumi_invfb.
  - If the tuple was produced with weight = sigma_ep_fb * tuple_lumi_invfb / Ngen,
    set --tuple-lumi-invfb to that value.
  - For Y(1S), use --weight-convention inclusive_upsilon unless the generator xsec is definitely sigma x BR(Y->ee).
"""

import argparse
import math
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument,
                                branching_fractions_from_args,
                                branching_fraction_metadata)

ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813
BR_EE_1S = 0.0234
M_UPSILON_1S = 9.46030
SIGMA_MASS_1S = 0.053


def photon_flux_dW(W, e_e, e_p, q2_max):
    """Equivalent photon flux dN_gamma/dW."""
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
    dydw = 2.0 * W / s_minus_mp2
    return dndy * dydw


def integrated_flux(w_low, w_high, e_e, e_p, q2_max, n_steps=500):
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    return float(np.trapezoid(vals, grid))


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def load_theory_gamma_pb(path):
    """Load CT18NLO gamma-p table. Input columns are assumed nb. Output pb."""
    if path is None or path == "":
        return None

    path = os.path.expanduser(path)
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    raw = pd.read_csv(path)
    if "s_gp" not in raw.columns:
        raise ValueError(f"Theory file {path} does not contain s_gp column")

    # Envelope. The CSV has no explicit central column in the files used here.
    low_cols = [c for c in ["muR_min", "muF_min", "PDFerr_min"] if c in raw.columns]
    high_cols = [c for c in ["muR_max", "muF_max", "PDFerr_max"] if c in raw.columns]
    if len(low_cols) == 0 or len(high_cols) == 0:
        raise ValueError("Could not identify min/max columns in theory CSV")

    low_nb = raw[low_cols].min(axis=1).to_numpy()
    high_nb = raw[high_cols].max(axis=1).to_numpy()
    mid_nb = 0.5 * (low_nb + high_nb)

    out = pd.DataFrame({
        "W": raw["s_gp"].to_numpy(dtype=float),
        "sigma_low_pb": 1000.0 * low_nb,
        "sigma_mid_pb": 1000.0 * mid_nb,
        "sigma_high_pb": 1000.0 * high_nb,
    })
    out = out.sort_values("W")
    return out


def theory_bin_average_and_ep(theory, w_bins, e_e, e_p, q2_max):
    """Return CT18 bin-averaged gamma-p pb and flux-folded ep cross section in fb."""
    if theory is None:
        return None

    Wtab = theory["W"].to_numpy()
    sig_mid = theory["sigma_mid_pb"].to_numpy()
    sig_low = theory["sigma_low_pb"].to_numpy()
    sig_high = theory["sigma_high_pb"].to_numpy()

    rows = []
    for lo, hi in zip(w_bins[:-1], w_bins[1:]):
        grid = np.linspace(lo, hi, 500)
        flux = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
        smid = np.interp(grid, Wtab, sig_mid, left=np.nan, right=np.nan)
        slow = np.interp(grid, Wtab, sig_low, left=np.nan, right=np.nan)
        shigh = np.interp(grid, Wtab, sig_high, left=np.nan, right=np.nan)

        valid = np.isfinite(smid) & np.isfinite(flux)
        if not np.any(valid):
            rows.append({
                "theory_gamma_mid_pb": np.nan,
                "theory_gamma_low_pb": np.nan,
                "theory_gamma_high_pb": np.nan,
                "theory_ep_mid_fb": np.nan,
                "theory_flux_int": np.nan,
            })
            continue

        f_int = float(np.trapezoid(flux[valid], grid[valid]))
        if f_int <= 0:
            rows.append({
                "theory_gamma_mid_pb": np.nan,
                "theory_gamma_low_pb": np.nan,
                "theory_gamma_high_pb": np.nan,
                "theory_ep_mid_fb": np.nan,
                "theory_flux_int": f_int,
            })
            continue

        # Flux-weighted gamma-p bin average.
        gmid = float(np.trapezoid(flux[valid] * smid[valid], grid[valid]) / f_int)
        glow = float(np.trapezoid(flux[valid] * slow[valid], grid[valid]) / f_int)
        ghigh = float(np.trapezoid(flux[valid] * shigh[valid], grid[valid]) / f_int)

        # ep cross section pb = integral flux(W) sigma_gamma(W) dW.
        ep_pb = float(np.trapezoid(flux[valid] * smid[valid], grid[valid]))
        ep_fb = 1000.0 * ep_pb

        rows.append({
            "theory_gamma_mid_pb": gmid,
            "theory_gamma_low_pb": glow,
            "theory_gamma_high_pb": ghigh,
            "theory_ep_mid_fb": ep_fb,
            "theory_flux_int": f_int,
        })

    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    add_branching_fraction_argument(p)
    p.add_argument("input", help="Input ROOT file with UpsilonTree")
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--out-csv", default="sigma_gp_1S_debug_table.csv")
    p.add_argument("--theory-csv", default="", help="CT18NLO gamma-p CSV, e.g. Fig5b_top_plot_CT18NLO.csv")
    p.add_argument("--electron-energy", type=float, default=10.0)
    p.add_argument("--proton-energy", type=float, default=100.0)
    p.add_argument("--q2-max", type=float, default=0.01)
    p.add_argument("--w-min", type=float, default=10.0)
    p.add_argument("--w-max", type=float, default=50.0)
    p.add_argument("--n-w-bins", type=int, default=14)
    p.add_argument("--lumi-invfb", type=float, default=10.0, help="Luminosity to use for extracted cross sections")
    p.add_argument("--tuple-lumi-invfb", type=float, default=10.0, help="Luminosity used when writing the tuple weights")
    p.add_argument("--expected-generator-xsec-fb", type=float, default=119.685, help="Expected total ep generator xsec used in tuple weight, for sanity checks")
    p.add_argument("--weight-convention", choices=["inclusive_upsilon", "visible_decay"], default="inclusive_upsilon")
    p.add_argument("--mass-window-sigma", type=float, default=3.0)
    p.add_argument("--true-mass-window", type=float, default=0.2)
    p.add_argument("--sigma-epz-center", type=float, default=None)
    p.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    p.add_argument("--no-selection", action="store_true", help="Skip reco mass/Sigma/W selection and just inspect true W")
    args = p.parse_args()
    br_ee = branching_fractions_from_args(args)["1S"]

    infile = os.path.expanduser(args.input)
    if not os.path.exists(infile):
        print(f"ERROR: input file not found: {infile}", file=sys.stderr)
        sys.exit(1)

    sigma_epz_center = 2.0 * args.electron_energy if args.sigma_epz_center is None else args.sigma_epz_center
    sigma_epz_halfwidth = 0.2 * args.electron_energy if args.sigma_epz_halfwidth is None else args.sigma_epz_halfwidth

    df = ROOT.RDataFrame(args.tree, infile)

    # Basic checks
    n_entries = int(df.Count().GetValue())
    sum_weight_all = float(df.Sum("weight").GetValue())
    inferred_xsec_fb = sum_weight_all / args.tuple_lumi_invfb if args.tuple_lumi_invfb > 0 else np.nan

    print("\n=== Tuple normalization sanity check ===")
    print(f"entries                     = {n_entries}")
    print(f"sum weight all              = {sum_weight_all:.6g}")
    print(f"tuple_lumi_invfb            = {args.tuple_lumi_invfb:.6g}")
    print(f"inferred xsec from weights  = {inferred_xsec_fb:.6g} fb")
    print(f"expected generator xsec     = {args.expected_generator_xsec_fb:.6g} fb")
    if args.expected_generator_xsec_fb > 0:
        print(f"ratio inferred/expected     = {inferred_xsec_fb / args.expected_generator_xsec_fb:.6g}")

    if abs(inferred_xsec_fb / args.expected_generator_xsec_fb - 1.0) > 0.05:
        print("WARNING: sum(weight)/tuple_lumi_invfb is not close to the expected generator xsec.")
        print("         Check whether the tuple is already selected, merged, or using a different xsec/lumi.")

    w_bins = np.linspace(args.w_min, args.w_max, args.n_w_bins + 1)
    w_centers = 0.5 * (w_bins[:-1] + w_bins[1:])

    mass_window = args.mass_window_sigma * SIGMA_MASS_1S
    reco_mass_min = M_UPSILON_1S - mass_window
    reco_mass_max = M_UPSILON_1S + mass_window
    true_mass_min = M_UPSILON_1S - args.true_mass_window
    true_mass_max = M_UPSILON_1S + args.true_mass_window

    if args.no_selection:
        reco_filter = f"true_mass >= {true_mass_min} && true_mass <= {true_mass_max} && true_W >= {args.w_min} && true_W <= {args.w_max}"
        reco_cols = ["true_W", "weight"]
        reco_w_branch = "true_W"
    else:
        reco_filter = (
            f"reco_mass >= {reco_mass_min} && reco_mass <= {reco_mass_max}"
            f" && abs(sigma_epz - {sigma_epz_center}) <= {sigma_epz_halfwidth}"
            f" && reco_W >= {args.w_min} && reco_W <= {args.w_max}"
        )
        reco_cols = ["reco_W", "weight"]
        reco_w_branch = "reco_W"

    true_filter = (
        f"true_mass >= {true_mass_min} && true_mass <= {true_mass_max}"
        f" && true_W >= {args.w_min} && true_W <= {args.w_max}"
    )

    reco = df.Filter(reco_filter).AsNumpy(reco_cols)
    truth = df.Filter(true_filter).AsNumpy(["true_W", "weight"])

    n_reco, err_reco = hist_weighted(reco[reco_w_branch], reco["weight"], w_bins)
    n_true, err_true = hist_weighted(truth["true_W"], truth["weight"], w_bins)

    eff = np.divide(n_reco, n_true, out=np.zeros_like(n_reco), where=n_true > 0.0)

    flux_int = np.array([
        integrated_flux(w_bins[i], w_bins[i + 1], args.electron_energy, args.proton_energy, args.q2_max)
        for i in range(len(w_centers))
    ])

    # Extracted gamma-p cross section. Algebraically n_reco/eff = n_true if eff is nonzero.
    br_corr = 1.0 / br_ee if args.weight_convention == "visible_decay" else 1.0
    denom = args.lumi_invfb * 1000.0 * flux_int
    sigma_from_truth_pb = np.divide(n_true, denom, out=np.full_like(n_true, np.nan), where=denom > 0) * br_corr
    sigma_from_reco_pb = np.divide(n_reco, args.lumi_invfb * 1000.0 * flux_int * eff,
                                   out=np.full_like(n_reco, np.nan),
                                   where=(flux_int > 0) & (eff > 0)) * br_corr

    # Statistical error from reco numerator, with same bin-by-bin correction.
    err_sigma_pb = np.divide(err_reco, args.lumi_invfb * 1000.0 * flux_int * eff,
                             out=np.full_like(err_reco, np.nan),
                             where=(flux_int > 0) & (eff > 0)) * br_corr

    theory = load_theory_gamma_pb(args.theory_csv) if args.theory_csv else None
    theory_bins = theory_bin_average_and_ep(theory, w_bins, args.electron_energy, args.proton_energy, args.q2_max)
    if theory_bins is None:
        theory_bins = pd.DataFrame({})

    out = pd.DataFrame({
        "br_ee": br_ee,
        "branching_fraction_convention": branching_fraction_metadata(args)["branching_fraction_convention"],
        "branching_fraction_correction_applied": args.weight_convention == "visible_decay",
        "branching_fraction_uncertainty_included": False,
        "tuple_weight_convention": args.weight_convention,
        "W_low": w_bins[:-1],
        "W_high": w_bins[1:],
        "W_center": w_centers,
        "flux_int": flux_int,
        "n_true_weighted": n_true,
        "n_reco_weighted": n_reco,
        "eff_reco_over_true": eff,
        "sigma_gp_from_truth_pb": sigma_from_truth_pb,
        "sigma_gp_from_reco_pb": sigma_from_reco_pb,
        "sigma_gp_stat_err_pb": err_sigma_pb,
        "ep_xsec_in_bin_inclusive_upsilon_fb": n_true * br_corr / args.lumi_invfb,
    })

    if len(theory_bins) == len(out):
        out = pd.concat([out, theory_bins], axis=1)
        out["ratio_data_to_theory_mid"] = out["sigma_gp_from_truth_pb"] / out["theory_gamma_mid_pb"]
        out["ep_xsec_in_bin_from_tuple_fb"] = out["n_true_weighted"] / args.lumi_invfb
        # Preserve the archived diagnostic exactly, but expose its convention.
        # For visible-decay tuples, only the new inclusive ratio is comparable
        # directly with the inclusive CT18 photon-flux-folded prediction.
        out["ratio_ep_bin_tuple_to_theory"] = out["ep_xsec_in_bin_from_tuple_fb"] / out["theory_ep_mid_fb"]
        out["ratio_ep_bin_tuple_to_theory_is_legacy"] = True
        out["ratio_ep_bin_tuple_to_theory_convention"] = args.weight_convention + "_over_inclusive_upsilon"
        out["ratio_ep_bin_inclusive_to_theory"] = out["ep_xsec_in_bin_inclusive_upsilon_fb"] / out["theory_ep_mid_fb"]
        out["theory_ep_convention"] = "inclusive_upsilon"

    out_path = os.path.expanduser(args.out_csv)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    out.to_csv(out_path, index=False)

    print("\n=== Per-bin debug table ===")
    cols = [
        "W_low", "W_high", "flux_int", "n_true_weighted", "n_reco_weighted", "eff_reco_over_true",
        "sigma_gp_from_truth_pb"
    ]
    if "theory_gamma_mid_pb" in out.columns:
        cols += ["theory_gamma_mid_pb", "ratio_data_to_theory_mid", "ep_xsec_in_bin_from_tuple_fb",
                 "ep_xsec_in_bin_inclusive_upsilon_fb", "theory_ep_mid_fb",
                 "ratio_ep_bin_inclusive_to_theory"]
    print(out[cols].to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    print(f"\nWrote {out_path}")

    if "theory_ep_mid_fb" in out.columns:
        ep_tuple_fb = float(np.nansum(out["ep_xsec_in_bin_from_tuple_fb"]))
        ep_inclusive_fb = float(np.nansum(out["ep_xsec_in_bin_inclusive_upsilon_fb"]))
        ep_theory_fb = float(np.nansum(out["theory_ep_mid_fb"]))
        print("\n=== Integrated over requested W range ===")
        print(f"tuple ep xsec ({args.weight_convention}) = {ep_tuple_fb:.6g} fb")
        print(f"tuple inclusive ep xsec       = {ep_inclusive_fb:.6g} fb")
        print(f"CT18 inclusive ep xsec        = {ep_theory_fb:.6g} fb")
        print(f"historical tuple/theory ratio = {ep_tuple_fb / ep_theory_fb if ep_theory_fb > 0 else np.nan:.6g} ({args.weight_convention}/inclusive_upsilon)")
        print(f"inclusive tuple/theory ratio  = {ep_inclusive_fb / ep_theory_fb if ep_theory_fb > 0 else np.nan:.6g}")
        print("Use the inclusive tuple/theory ratio for a common decay convention before unfolding to gamma-p.")


if __name__ == "__main__":
    main()
