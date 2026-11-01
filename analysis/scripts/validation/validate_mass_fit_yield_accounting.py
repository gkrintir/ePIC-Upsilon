#!/usr/bin/env python3
"""
Cross-check the Υ signal yields used in the mass fit and in the sigma(gamma p) extraction.

Purpose
-------
Clarify the difference between:
  1. reconstructed visible e+e- yields in the mass fit,
  2. reconstructed yields with BR unfolded,
  3. efficiency-corrected visible yields used as the sigma(gamma p) numerator,
  4. efficiency-corrected inclusive Υ yields after BR unfolding.

The script computes all of these from the same tuple.

No pandas dependency.
"""

import argparse
import os
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata
import uproot
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument,
                                branching_fractions_from_args,
                                branching_fraction_metadata)

STATES = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br": 0.0234, "visible_xsec_fb": 327.584},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br": 0.0191, "visible_xsec_fb": 84.780},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br": 0.0218, "visible_xsec_fb": 68.929},
}


def parse_xsecs(s):
    vals = {}
    if not s:
        return vals
    for item in s.split(","):
        if not item.strip():
            continue
        key, val = item.split(":")
        vals[key.strip()] = float(val)
    return vals


def hist_weighted(values, weights, bins):
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def weighted_sum(w, mask):
    return float(np.sum(w[mask]))


def main():
    ap = argparse.ArgumentParser()
    add_branching_fraction_argument(ap)
    ap.add_argument("input", help="Input ROOT file with UpsilonTree")
    ap.add_argument("--tree", default="UpsilonTree")
    ap.add_argument("--states", default="1S,2S,3S",
                    help="Comma-separated truth states to report, e.g. 1S")
    ap.add_argument("--lumi-invfb", type=float, default=10.0)
    ap.add_argument("--electron-energy", type=float, default=18.0)
    ap.add_argument("--sigma-epz-center", type=float, default=None)
    ap.add_argument("--sigma-epz-halfwidth", type=float, default=None)

    # Mass-fit selection
    ap.add_argument("--massfit-mass-min", type=float, default=8.0)
    ap.add_argument("--massfit-mass-max", type=float, default=12.0)
    ap.add_argument("--massfit-w-min", type=float, default=0.0)

    # Cross-section selection
    ap.add_argument("--w-min", type=float, default=20.0)
    ap.add_argument("--w-max", type=float, default=140.0)
    ap.add_argument("--n-w-bins", type=int, default=14)
    ap.add_argument("--truth-mass-window", type=float, default=0.20)
    ap.add_argument("--q2-cut", type=float, default=None,
                    help="Optional true Q² upper cut in GeV², applied to truth and reconstructed selections")
    ap.add_argument("--reco-mass-window-sigma", type=float, default=3.0)

    # Optional state-scoped normalization to expected visible xsecs.
    ap.add_argument(
        "--state-normalize",
        action="store_true",
        help="Also compute columns after scaling each truth state so sum(weight)/L equals the visible xsec.",
    )
    ap.add_argument(
        "--visible-xsecs-fb",
        default="1S:327.584,2S:84.780,3S:68.929",
        help="Comma-separated visible xsecs, e.g. 1S:327.584,2S:84.780,3S:68.929",
    )
    ap.add_argument("--out-csv", default="yield_crosscheck_massfit_vs_sigma.csv")
    args = ap.parse_args()
    branching_fractions = branching_fractions_from_args(args)

    selected_states = [state.strip() for state in args.states.split(",") if state.strip()]
    if not selected_states or any(state not in STATES for state in selected_states):
        ap.error("--states must contain one or more of 1S,2S,3S")
    selected_states = list(dict.fromkeys(selected_states))

    if args.sigma_epz_center is None:
        args.sigma_epz_center = 2.0 * args.electron_energy
    if args.sigma_epz_halfwidth is None:
        args.sigma_epz_halfwidth = 0.2 * args.electron_energy

    visible_xsecs = parse_xsecs(args.visible_xsecs_fb)

    filename = os.path.expanduser(args.input)
    cols = ["weight", "true_mass", "true_W", "reco_mass", "reco_W", "sigma_epz"]
    if args.q2_cut is not None:
        cols.append("true_Q2")
    cols += proton_branches()
    with uproot.open(filename) as root_file:
        tree = root_file[args.tree]
        arr = {name: tree[name].array(library="np") for name in cols}
    data = {k: np.asarray(v, dtype=float) for k, v in arr.items()}

    # Validate the full producer contract before any existing finite-value filter.
    proton_tag_mask(data)
    finite = np.ones(len(data["weight"]), dtype=bool)
    for name, v in data.items():
        if name not in proton_branches():
            finite &= np.isfinite(v)
    data = {k: v[finite] for k, v in data.items()}

    w_bins = np.linspace(args.w_min, args.w_max, args.n_w_bins + 1)

    rows = []
    q2_mask = np.ones(len(data["weight"]), dtype=bool)
    if args.q2_cut is not None:
        q2_mask &= (data["true_Q2"] >= 0.0) & (data["true_Q2"] <= args.q2_cut)

    common_massfit = (
        q2_mask
        & proton_tag_mask(data)
        & (data["reco_mass"] >= args.massfit_mass_min)
        & (data["reco_mass"] <= args.massfit_mass_max)
        & (np.abs(data["sigma_epz"] - args.sigma_epz_center) <= args.sigma_epz_halfwidth)
        & (data["reco_W"] >= args.massfit_w_min)
    )

    for name in selected_states:
        props = dict(STATES[name], br=branching_fractions[name],
                     visible_xsec_fb=visible_xsecs.get(name, STATES[name]["visible_xsec_fb"]))
        mass = props["mass"]
        br = props["br"]
        visible_xsec = props["visible_xsec_fb"]

        truth_state_all = (
            (data["true_mass"] >= mass - args.truth_mass_window)
            & (data["true_mass"] <= mass + args.truth_mass_window)
        )

        sumw_state_raw_all = weighted_sum(data["weight"], truth_state_all)
        raw_visible_xsec_fb = sumw_state_raw_all / args.lumi_invfb if args.lumi_invfb > 0 else np.nan
        state_scale = visible_xsec / raw_visible_xsec_fb if raw_visible_xsec_fb > 0 else np.nan

        # Raw and optionally normalized weights.
        weights_modes = [("raw", data["weight"])]
        if args.state_normalize:
            weights_modes.append(("stateNorm", data["weight"] * state_scale))

        reco_state_window = (
            (data["reco_mass"] >= mass - args.reco_mass_window_sigma * props["sigma"])
            & (data["reco_mass"] <= mass + args.reco_mass_window_sigma * props["sigma"])
        )

        truth_w_range = (
            truth_state_all
            & q2_mask
            & (data["true_W"] >= args.w_min)
            & (data["true_W"] <= args.w_max)
        )

        reco_sigma_range = (
            reco_state_window
            & proton_tag_mask(data)
            & q2_mask
            & (np.abs(data["sigma_epz"] - args.sigma_epz_center) <= args.sigma_epz_halfwidth)
            & (data["reco_W"] >= args.w_min)
            & (data["reco_W"] <= args.w_max)
        )

        massfit_truth_labeled = common_massfit & truth_state_all

        for mode, w in weights_modes:
            massfit_visible = weighted_sum(w, massfit_truth_labeled)
            massfit_inclusive = massfit_visible / br

            n_reco, err_reco = hist_weighted(data["reco_W"][reco_sigma_range], w[reco_sigma_range], w_bins)
            n_true, err_true = hist_weighted(data["true_W"][truth_w_range], w[truth_w_range], w_bins)

            eff = np.divide(n_reco, n_true, out=np.zeros_like(n_reco), where=n_true > 0)

            # Usable bin-by-bin unfolded numerator only where reco and eff are nonzero.
            good = (n_reco > 0) & (eff > 0)
            unfolded_visible_bins = np.zeros_like(n_reco)
            unfolded_visible_bins[good] = n_reco[good] / eff[good]

            reco_visible_sum = float(np.sum(n_reco))
            truth_visible_sum = float(np.sum(n_true))
            unfolded_visible_sum = float(np.sum(unfolded_visible_bins))

            rows.append({
                "proton_selection": proton_selection_metadata()["mode"],
                "state": name,
                "mode": mode,
                "br_ee": br,
                "branching_fraction_convention": branching_fraction_metadata(args)["branching_fraction_convention"],
                "branching_fraction_uncertainty_included": False,
                "visible_xsec_target_fb": visible_xsec,
                "raw_visible_xsec_from_truth_state_fb": raw_visible_xsec_fb,
                "state_scale_applied": state_scale if mode == "stateNorm" else 1.0,

                # What mass fit should compare to.
                "massfit_truth_labeled_reco_visible_yield": massfit_visible,
                "massfit_truth_labeled_reco_inclusive_yield": massfit_inclusive,

                # What W cross-section extraction uses.
                "sigma_reco_visible_yield_in_W_range": reco_visible_sum,
                "sigma_truth_visible_yield_in_W_range": truth_visible_sum,
                "sigma_unfolded_visible_yield_sum_good_bins": unfolded_visible_sum,

                "sigma_reco_inclusive_yield_in_W_range": reco_visible_sum / br,
                "sigma_truth_inclusive_yield_in_W_range": truth_visible_sum / br,
                "sigma_unfolded_inclusive_yield_sum_good_bins": unfolded_visible_sum / br,

                "average_eff_in_W_range": reco_visible_sum / truth_visible_sum if truth_visible_sum > 0 else np.nan,
                "n_good_W_bins": int(np.sum(good)),
                "n_total_W_bins": len(w_bins) - 1,
            })

    fieldnames = list(rows[0].keys())
    with open(args.out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print("=== Yield cross-check ===")
    for r in rows:
        print("")
        print(f"{r['state']} [{r['mode']}]")
        print(f"  raw visible xsec from truth state: {r['raw_visible_xsec_from_truth_state_fb']:.3f} fb")
        print(f"  state scale applied: {r['state_scale_applied']:.6g}")
        print(f"  mass-fit comparable reco visible yield: {r['massfit_truth_labeled_reco_visible_yield']:.3f}")
        print(f"  mass-fit comparable reco inclusive yield: {r['massfit_truth_labeled_reco_inclusive_yield']:.3f}")
        print(f"  sigma reco visible yield in W range: {r['sigma_reco_visible_yield_in_W_range']:.3f}")
        print(f"  sigma truth/unfolded visible yield in W range: {r['sigma_truth_visible_yield_in_W_range']:.3f} / {r['sigma_unfolded_visible_yield_sum_good_bins']:.3f}")
        print(f"  sigma truth/unfolded inclusive yield in W range: {r['sigma_truth_inclusive_yield_in_W_range']:.3f} / {r['sigma_unfolded_inclusive_yield_sum_good_bins']:.3f}")
        print(f"  average efficiency in W range: {r['average_eff_in_W_range']:.4f}")
        print(f"  good W bins: {r['n_good_W_bins']}/{r['n_total_W_bins']}")

    print(f"\nWrote {args.out_csv}")
    print("\nInterpretation:")
    print("  * Mass-fit visible mode should compare to massfit_truth_labeled_reco_visible_yield.")
    print("  * Mass-fit inclusive mode should compare to massfit_truth_labeled_reco_inclusive_yield.")
    print("  * sigma(gamma p) extraction uses the efficiency-corrected numerator:")
    print("      visible:   sigma_unfolded_visible_yield_sum_good_bins")
    print("      inclusive: sigma_unfolded_inclusive_yield_sum_good_bins")
    print("    This is larger than the mass-fit reco yield by roughly 1/efficiency.")


if __name__ == "__main__":
    main()
