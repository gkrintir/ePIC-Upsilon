#!/usr/bin/env python3
"""
Extract exclusive sigma(gamma p -> Upsilon(nS) p) versus W_gamma_p from a flattened UpsilonTree.

Main fixes relative to the earlier script
  - no notebook argument override
  - visible e+e- cross-section normalization per state
  - state-scoped normalization, so merged 1S+2S+3S tuples are handled correctly
  - optional Q2 truth cut
  - q2-max default set by user, not hidden
  - CT18NLO overlay is opt-in and labelled as inelastic, not as the exclusive prediction
  - writes a CSV debug table with yields, efficiencies, fluxes, bin-by-bin unfolded cross sections, closure, and uncertainties

Recommended convention
  Store or renormalize tuple weights to the visible e+e- cross section.
  Use --weight-convention visible_decay to unfold the branching ratio for the inclusive
  gamma-p exclusive production cross section.

Example for 18x275 Upsilon(1S)
  python3 extract_sigma_gamma_p_family_fixed_swan.py merged_elecs_test_tuples_fix.root \
    --states 1S \
    --electron-energy 18 \
    --proton-energy 275 \
    --w-min 20 \
    --w-max 140 \
    --q2-max 0.001 \
    --lumi-invfb 10 \
    --tuple-lumi-invfb 10 \
    --visible-xsec-1s-fb 327.584 \
    --weight-convention visible_decay \
    --out-prefix sigma_gp_1S_18x275_fixed

Example for 10x100 Upsilon(1S)
  python3 extract_sigma_gamma_p_family_fixed_swan.py merged_elecs_test_10t100.root \
    --states 1S \
    --electron-energy 10 \
    --proton-energy 100 \
    --w-min 20 \
    --w-max 50 \
    --q2-max 0.001 \
    --lumi-invfb 10 \
    --tuple-lumi-invfb 10 \
    --visible-xsec-1s-fb 104.793 \
    --weight-convention visible_decay \
    --out-prefix sigma_gp_1S_10x100_fixed
"""

import argparse
import math
import os
import sys
from typing import Dict, List, Tuple

import uproot
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

sys.path.insert(0, os.path.dirname(__file__))
from statistical_components import (configure_statistics, metadata as stat_metadata, statistical_label,
                                    binned_components, covariance_rows, write_covariance,
                                    add_interval_arguments, binned_intervals, draw_interval_bars)

ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813

UPSILON_STATES = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0234, "color": "royalblue"},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191, "color": "crimson"},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218, "color": "forestgreen"},
}


def branch_exists(root_file: str, tree_name: str, branch_name: str) -> bool:
    with uproot.open(os.path.expanduser(root_file)) as source:
        return tree_name in source and branch_name in source[tree_name].keys()


def photon_flux_dW(W: float, e_e: float, e_p: float, q2_max: float) -> float:
    """Equivalent photon flux dN_gamma/dW for ep photoproduction."""
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


def integrated_flux(w_low: float, w_high: float, e_e: float, e_p: float, q2_max: float, n_steps: int = 500) -> float:
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:
        trapezoid = np.trapz
    return float(trapezoid(vals, grid))


def hist_weighted(values: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    y, _ = np.histogram(values, bins=bins, weights=weights)
    e2, _ = np.histogram(values, bins=bins, weights=weights * weights)
    return y.astype(float), np.sqrt(e2.astype(float))


def parse_states(s: str) -> List[str]:
    out = [x.strip() for x in s.split(",") if x.strip()]
    bad = [x for x in out if x not in UPSILON_STATES]
    if bad:
        raise ValueError(f"Unknown states {bad}. Available states are {list(UPSILON_STATES)}")
    return out


def load_theory_gamma_pb(path: str):
    """Load the uploaded CT18NLO inelastic gamma-p curve. This is not the exclusive prediction."""
    if not path:
        return None
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        print(f"Theory CSV {path} not found. Skipping.")
        return None

    data = np.genfromtxt(path, delimiter=",", names=True)
    required = ["s_gp", "muR_min", "muF_min", "PDFerr_min", "muR_max", "muF_max", "PDFerr_max"]
    if not all(name in data.dtype.names for name in required):
        print(f"Theory CSV {path} does not contain the expected CT18 columns. Skipping.")
        return None
    return data


def plot_ct18_inelastic(ax, data):
    w = data["s_gp"]
    low_nb = np.minimum.reduce([data["muR_min"], data["muF_min"], data["PDFerr_min"]])
    high_nb = np.maximum.reduce([data["muR_max"], data["muF_max"], data["PDFerr_max"]])
    mid_nb = 0.5 * (low_nb + high_nb)

    ax.fill_between(
        w, 1000.0 * low_nb, 1000.0 * high_nb,
        color="tab:red", alpha=0.14, linewidth=0,
        label=r"CT18NLO inelastic $\Upsilon(1S)+X$",
        zorder=1,
    )
    ax.plot(w, 1000.0 * mid_nb, color="tab:red", linewidth=2.0, zorder=2)


def state_visible_xsecs(args) -> Dict[str, float]:
    return {
        "1S": args.visible_xsec_1s_fb,
        "2S": args.visible_xsec_2s_fb,
        "3S": args.visible_xsec_3s_fb,
    }


def make_mask_for_state(data, props, args, has_true_q2, reco=False):
    if reco:
        mass_window = args.mass_window_sigma * props["sigma"]
        mass_min = props["mass"] - mass_window
        mass_max = props["mass"] + mass_window
        mask = (
            (data["reco_mass"] >= mass_min) & (data["reco_mass"] <= mass_max)
            & (np.abs(data["sigma_epz"] - args.sigma_epz_center) <= args.sigma_epz_halfwidth)
            & (data["reco_W"] >= args.w_min) & (data["reco_W"] <= args.w_max)
        )
    else:
        mass_min = props["mass"] - args.true_mass_window
        mass_max = props["mass"] + args.true_mass_window
        mask = (
            (data["true_mass"] >= mass_min) & (data["true_mass"] <= mass_max)
            & (data["true_W"] >= args.w_min) & (data["true_W"] <= args.w_max)
        )

    if args.q2_cut is not None and has_true_q2:
        mask &= (data["true_Q2"] >= 0) & (data["true_Q2"] <= args.q2_cut)

    return mask


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input ROOT file containing UpsilonTree")
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--out-prefix", default="sigma_gp_upsilon_fixed")
    parser.add_argument("--states", default="1S", help="Comma-separated list, e.g. 1S or 1S,2S,3S")

    parser.add_argument("--lumi-invfb", type=float, default=10.0)
    parser.add_argument("--tuple-lumi-invfb", type=float, default=10.0)
    parser.add_argument("--state-norm-mode", choices=["raw", "target_visible_xsec"],
                        default="target_visible_xsec",
                        help="raw preserves stored weights after all analysis cuts; target_visible_xsec rescales the selected normalization scope.")

    parser.add_argument("--electron-energy", type=float, required=True)
    parser.add_argument("--proton-energy", type=float, required=True)
    parser.add_argument("--q2-max", type=float, default=0.001)
    parser.add_argument("--q2-cut", type=float, default=None, help="Optional event-level true_Q2 cut. Leave unset if the tuple already has the right Q2 range.")

    parser.add_argument("--w-min", type=float, default=20.0)
    parser.add_argument("--w-max", type=float, default=140.0)
    parser.add_argument("--n-w-bins", type=int, default=14)

    parser.add_argument("--sigma-epz-center", type=float, default=None)
    parser.add_argument("--sigma-epz-halfwidth", type=float, default=None)
    parser.add_argument("--mass-window-sigma", type=float, default=3.0)
    parser.add_argument("--true-mass-window", type=float, default=0.2)

    parser.add_argument("--visible-xsec-1s-fb", type=float, default=327.584)
    parser.add_argument("--visible-xsec-2s-fb", type=float, default=84.780)
    parser.add_argument("--visible-xsec-3s-fb", type=float, default=68.929)

    parser.add_argument(
        "--weight-convention",
        choices=["visible_decay", "inclusive_upsilon"],
        default="visible_decay",
        help="Use visible_decay when the normalization cross sections are for Upsilon -> ee. Use inclusive_upsilon only if xsecs are already inclusive.",
    )
    parser.add_argument(
        "--normalization-scope",
        choices=["truth_state", "all"],
        default="truth_state",
        help="Use truth_state for merged 1S+2S+3S tuples. all is only for pure single-state files.",
    )

    parser.add_argument("--min-eff", type=float, default=0.0, help="Only plot bins with reco/true efficiency above this threshold.")
    parser.add_argument("--include-eff-mc-unc", action="store_true", help="Legacy alias for --include-efficiency-mc-uncertainty 1; projected mode only.")
    parser.add_argument("--show-ct18-inelastic", action="store_true", help="Overlay CT18NLO inelastic curve as a reference only.")
    parser.add_argument("--theory-csv", default="/eos/user/g/gkrintir/Fig5b_top_plot_CT18NLO.csv")
    parser.add_argument("--no-hera", action="store_true")

    parser.add_argument("--stat-mode", choices=["mc", "projected"], default="mc")
    parser.add_argument("--projection-lumi-invfb", type=float, default=None,
                        help="Target luminosity, defaults to --lumi-invfb. Projected yields scale from --tuple-lumi-invfb; MC ignores this target.")
    parser.add_argument("--include-efficiency-mc-uncertainty", type=int, choices=[0, 1], default=0)
    add_interval_arguments(parser)
    args = configure_statistics(parser.parse_args())

    filename = os.path.expanduser(args.input)
    states = parse_states(args.states)
    visible_xsec = state_visible_xsecs(args)

    if args.sigma_epz_center is None:
        args.sigma_epz_center = 2.0 * args.electron_energy
    if args.sigma_epz_halfwidth is None:
        args.sigma_epz_halfwidth = 0.2 * args.electron_energy

    out_png = os.path.expanduser(args.out_prefix + ".png")
    out_csv = os.path.expanduser(args.out_prefix + ".csv")
    out_norm_csv = os.path.expanduser(args.out_prefix + "_normalization.csv")
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)

    has_true_q2 = branch_exists(filename, args.tree, "true_Q2")
    if args.q2_cut is not None and not has_true_q2:
        print("Requested --q2-cut but true_Q2 branch does not exist. Ignoring Q2 event cut.")

    columns = ["weight", "true_mass", "true_W", "reco_mass", "reco_W", "sigma_epz"]
    if has_true_q2:
        columns.append("true_Q2")
    with uproot.open(filename) as source:
        tree = source[args.tree]
        data = {name: tree[name].array(library="np") for name in columns}
    w_bins = np.linspace(args.w_min, args.w_max, args.n_w_bins + 1)
    w_centers = 0.5 * (w_bins[:-1] + w_bins[1:])

    theory_data = load_theory_gamma_pb(args.theory_csv) if args.show_ct18_inelastic else None

    fig, ax = plt.subplots(figsize=(9.5, 7.5))

    if theory_data is not None:
        plot_ct18_inelastic(ax, theory_data)

    if not args.no_hera:
        hera_w = np.array([95.0, 175.0])
        hera_w_err = np.array([35.0, 45.0])
        hera_sig = np.array([160.0, 321.0])
        hera_err = np.array([70.0, 144.0])
        ax.errorbar(
            hera_w, hera_sig, xerr=hera_w_err, yerr=hera_err,
            fmt="^", color="gray",
            label=r"HERA exclusive $\gamma p \rightarrow \Upsilon(1S)p$",
            markersize=8, linewidth=1.8,
        )

    all_rows = []
    norm_rows = []

    for name in states:
        props = UPSILON_STATES[name]
        target_visible_fb = visible_xsec[name]

        if args.state_norm_mode == "target_visible_xsec" and target_visible_fb <= 0:
            print(f"Skipping {name} because visible cross section is nonpositive")
            continue

        # Determine state-scoped normalization.
        if args.normalization_scope == "truth_state":
            norm_mask = (
                (data["true_mass"] >= props["mass"] - args.true_mass_window)
                & (data["true_mass"] <= props["mass"] + args.true_mass_window)
            )
            if args.q2_cut is not None and has_true_q2:
                norm_mask &= (data["true_Q2"] >= 0) & (data["true_Q2"] <= args.q2_cut)
            sumw_norm = float(np.sum(data["weight"][norm_mask], dtype=float))
            norm_scope_label = "truth_state"
        else:
            sumw_norm = float(np.sum(data["weight"], dtype=float))
            norm_scope_label = "all"

        inferred_xsec_fb = sumw_norm / args.tuple_lumi_invfb if args.tuple_lumi_invfb > 0 else np.nan
        if args.tuple_lumi_invfb <= 0:
            raise ValueError("--tuple-lumi-invfb must be positive")
        if not np.isfinite(inferred_xsec_fb) or inferred_xsec_fb <= 0:
            raise ValueError(
                f"No positive-weight normalization sample is available for Upsilon({name}); "
                "refusing to emit zero or NaN cross sections"
            )
        weight_scale = target_visible_fb / inferred_xsec_fb if args.state_norm_mode == "target_visible_xsec" else 1.0

        norm_rows.append({
            "state": name,
            "normalization_scope": norm_scope_label,
            "state_norm_mode": args.state_norm_mode,
            **stat_metadata(args),
            "target_visible_xsec_fb": target_visible_fb,
            "sum_weight_norm_scope_raw": sumw_norm,
            "inferred_xsec_norm_scope_fb": inferred_xsec_fb,
            "weight_scale": weight_scale,
            "br_ee": props["br_ee"],
            "target_inclusive_xsec_pb_if_visible": target_visible_fb / props["br_ee"] / 1000.0,
        })

        true_mask = make_mask_for_state(data, props, args, has_true_q2, reco=False)
        reco_mask = make_mask_for_state(data, props, args, has_true_q2, reco=True)
        truth = {name: data[name][true_mask] for name in ["true_W", "weight"]}
        reco = {name: data[name][reco_mask] for name in ["reco_W", "weight"]}

        if len(truth["true_W"]) == 0:
            print(f"Upsilon({name}) skipped because truth selection is empty")
            continue

        stats = binned_components(data["true_W"], data["reco_W"],
                                  np.asarray(data["weight"], float) * weight_scale,
                                  true_mask, reco_mask, w_bins, args)
        intervals = binned_intervals(stats, args)
        n_true, n_reco = stats["truth"], stats["reco"]
        err_true, err_reco = np.sqrt(stats["truth_mc_var"]), np.sqrt(stats["reco_mc_var"])
        eff = stats["efficiency"]
        eff_err = np.sqrt(np.maximum(np.diag(stats["efficiency_mc_cov"]), 0))
        br_corr = 1.0 / props["br_ee"] if args.weight_convention == "visible_decay" else 1.0
        fluxes = np.asarray([integrated_flux(w_bins[i], w_bins[i+1], args.electron_energy,
                                           args.proton_energy, args.q2_max)
                             for i in range(len(w_centers))])
        denom = args.lumi_invfb * 1000.0 * fluxes
        factor = np.divide(br_corr, denom, out=np.full_like(denom, np.nan), where=denom > 0)
        n_unfolded = stats["corrected"]
        sigma_unfold = n_unfolded * factor
        sigma_truth_closure = n_true * factor
        err_sigma_mc = np.sqrt(np.maximum(np.diag(stats["cov_mc"]), 0)) * factor
        err_sigma_projected = np.sqrt(np.maximum(np.diag(stats["cov_projected"]), 0)) * factor
        err_sigma_eff_mc = np.sqrt(np.maximum(np.diag(stats["cov_correction_mc"]), 0)) * factor
        err_sigma_total = stats["error"] * factor
        # This legacy column remains the reconstructed-count term with the
        # correction held fixed. It is not the same-event MC primary error.
        reco_stat_error = np.sqrt(n_reco) if args.stat_mode == "projected" else err_reco
        err_n_unfolded_reco_stat = np.divide(reco_stat_error, eff,
                                            out=np.full_like(eff, np.nan), where=stats["valid"])
        err_sigma_reco_stat = err_n_unfolded_reco_stat * factor
        write_covariance(args.out_prefix + "_" + name + "_statistical_covariance.csv",
                         covariance_rows(stats, factor, w_bins, {"state": name, **stat_metadata(args)}))

        for i in range(len(w_centers)):
            all_rows.append({
                "state": name,
                **stat_metadata(args),
                "W_low": w_bins[i],
                "W_high": w_bins[i + 1],
                "W_center": w_centers[i],
                "flux_int": fluxes[i],
                "n_true_weighted_visible": n_true[i],
                "n_true_weighted_visible_err": err_true[i],
                "n_reco_weighted_visible": n_reco[i],
                "n_reco_weighted_visible_err": err_reco[i],
                "eff_reco_over_true": eff[i],
                "eff_err_mc_approx": eff_err[i],
                "eff_err_mc_delta": eff_err[i],
                "efficiency_uncertainty_method": "same-event weighted ratio delta covariance; migrations retained",
                "n_unfolded_visible": n_unfolded[i],
                "n_unfolded_visible_err_reco_stat": err_n_unfolded_reco_stat[i],
                "sigma_gamma_p_pb": sigma_unfold[i],
                "sigma_gamma_p_err_pb": err_sigma_total[i],
                "reco_count_to_sigma_pb": factor[i]/eff[i] if stats["valid"][i] else np.nan,
                "sigma_gamma_p_interval_low_pb": intervals["lower"][i]*factor[i],
                "sigma_gamma_p_interval_high_pb": intervals["upper"][i]*factor[i],
                "sigma_gamma_p_err_low_pb": intervals["error_low"][i]*factor[i],
                "sigma_gamma_p_err_high_pb": intervals["error_high"][i]*factor[i],
                "sigma_gamma_p_count_interval_low_pb": intervals["count_lower"][i]*factor[i],
                "sigma_gamma_p_count_interval_high_pb": intervals["count_upper"][i]*factor[i],
                "n_reco_count_interval_low": intervals["reco_lower"][i],
                "n_reco_count_interval_high": intervals["reco_upper"][i],
                "n_reco_interval_reference_count": intervals["reference_count"][i],
                "sigma_gamma_p_err_mc_pb": err_sigma_mc[i],
                "sigma_gamma_p_err_projected_reco_stat_pb": err_sigma_projected[i],
                "n_unfolded_visible_err": stats["error"][i],
                "sigma_gamma_p_err_reco_stat_pb": err_sigma_reco_stat[i],
                "sigma_gamma_p_err_eff_mc_pb": err_sigma_eff_mc[i],
                "sigma_gamma_p_truth_closure_pb": sigma_truth_closure[i],
                "closure_unfold_over_truth": (
                    sigma_unfold[i] / sigma_truth_closure[i]
                    if np.isfinite(sigma_unfold[i]) and np.isfinite(sigma_truth_closure[i]) and sigma_truth_closure[i] > 0
                    else np.nan
                ),
                "plotted": bool(np.isfinite(sigma_unfold[i]) and sigma_unfold[i] > 0.0 and eff[i] > args.min_eff),
            })

        mask = np.isfinite(sigma_unfold) & (sigma_unfold > 0.0) & (eff > args.min_eff)
        draw_interval_bars(
            ax, w_centers[mask], sigma_unfold[mask],
            (intervals["lower"]*factor)[mask], (intervals["upper"]*factor)[mask],
            fmt="o", color=props["color"],
            label=rf"ePIC $\Upsilon({name})$, $\epsilon>{args.min_eff:g}$",
            markersize=7, linewidth=1.8, elinewidth=1.4, capsize=3,
        )

    pd.DataFrame(all_rows).to_csv(out_csv, index=False)
    pd.DataFrame(norm_rows).to_csv(out_norm_csv, index=False)

    sqrts = math.sqrt(M_P * M_P + 4.0 * args.electron_energy * args.proton_energy)
    ax.set_xlabel(r"$W_{\gamma p}$ [GeV]", fontsize=14)
    ax.set_ylabel(r"$\sigma(\gamma p \rightarrow \Upsilon(nS)p)$ [pb]", fontsize=14)
    ax.set_yscale("log")
    #ax.set_xscale("log")
    ax.set_xlim(args.w_min, args.w_max)
    ax.set_ylim(0.5, 2000.0)

    ax.text(0.04, 0.94, "ePIC Performance", transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="top")
    ax.text(0.04, 0.88, rf"$e+p$, {args.electron_energy:.0f}$\times${args.proton_energy:.0f} GeV", transform=ax.transAxes, fontsize=11, ha="left", va="top")
    ax.text(0.04, 0.82, rf"$\sqrt{{s}} \simeq {sqrts:.0f}$ GeV, $Q^2 < {args.q2_max:g}$ GeV$^2$", transform=ax.transAxes, fontsize=11, ha="left", va="top")
    ax.text(0.04, 0.76, rf"$\mathcal{{L}}_{{ep}} = {args.lumi_invfb:g}$ fb$^{{-1}}$", transform=ax.transAxes, fontsize=11, ha="left", va="top")

    ax.legend(loc="lower right", fontsize=10, frameon=False)
    fig.text(0.5, 0.012, statistical_label(args),
             ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(out_png, dpi=300)

    print("=== Normalization summary ===")
    print(pd.DataFrame(norm_rows).to_string(index=False))
    print(f"\nSaved {out_png}")
    print(f"Saved {out_csv}")
    print(f"Saved {out_norm_csv}")


if __name__ == "__main__":
    main()
