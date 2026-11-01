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
  - writes a CSV debug table with yields, efficiencies, fluxes, and extracted cross sections

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
from typing import Dict, List, Tuple

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

UPSILON_STATES = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0234, "color": "royalblue"},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191, "color": "crimson"},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218, "color": "forestgreen"},
}


def branch_exists(root_file: str, tree_name: str, branch_name: str) -> bool:
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        return False
    t = f.Get(tree_name)
    ok = bool(t and t.GetBranch(branch_name))
    f.Close()
    return ok


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
    return float(np.trapz(vals, grid))


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


def make_filter_for_state(props, args, has_true_q2, reco=False):
    if reco:
        mass_window = args.mass_window_sigma * props["sigma"]
        mass_min = props["mass"] - mass_window
        mass_max = props["mass"] + mass_window
        expr = (
            f"reco_mass >= {mass_min:.12g} && reco_mass <= {mass_max:.12g}"
            f" && abs(sigma_epz - {args.sigma_epz_center:.12g}) <= {args.sigma_epz_halfwidth:.12g}"
            f" && reco_W >= {args.w_min:.12g} && reco_W <= {args.w_max:.12g}"
        )
    else:
        mass_min = props["mass"] - args.true_mass_window
        mass_max = props["mass"] + args.true_mass_window
        expr = (
            f"true_mass >= {mass_min:.12g} && true_mass <= {mass_max:.12g}"
            f" && true_W >= {args.w_min:.12g} && true_W <= {args.w_max:.12g}"
        )

    if args.q2_cut is not None and has_true_q2:
        expr += f" && true_Q2 >= 0 && true_Q2 <= {args.q2_cut:.12g}"

    return expr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input ROOT file containing UpsilonTree")
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--out-prefix", default="sigma_gp_upsilon_fixed")
    parser.add_argument("--states", default="1S", help="Comma-separated list, e.g. 1S or 1S,2S,3S")

    parser.add_argument("--lumi-invfb", type=float, default=10.0)
    parser.add_argument("--tuple-lumi-invfb", type=float, default=10.0)

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
    parser.add_argument("--show-ct18-inelastic", action="store_true", help="Overlay CT18NLO inelastic curve as a reference only.")
    parser.add_argument("--theory-csv", default="/eos/user/g/gkrintir/Fig5b_top_plot_CT18NLO.csv")
    parser.add_argument("--no-hera", action="store_true")

    args = parser.parse_args()

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

    df = ROOT.RDataFrame(args.tree, filename)
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

        if target_visible_fb <= 0:
            print(f"Skipping {name} because visible cross section is nonpositive")
            continue

        # Determine state-scoped normalization.
        if args.normalization_scope == "truth_state":
            norm_expr = (
                f"true_mass >= {props['mass'] - args.true_mass_window:.12g}"
                f" && true_mass <= {props['mass'] + args.true_mass_window:.12g}"
            )
            if args.q2_cut is not None and has_true_q2:
                norm_expr += f" && true_Q2 >= 0 && true_Q2 <= {args.q2_cut:.12g}"
            sumw_norm = float(df.Filter(norm_expr).Sum("weight").GetValue())
            norm_scope_label = "truth_state"
        else:
            sumw_norm = float(df.Sum("weight").GetValue())
            norm_scope_label = "all"

        inferred_xsec_fb = sumw_norm / args.tuple_lumi_invfb if args.tuple_lumi_invfb > 0 else np.nan
        weight_scale = target_visible_fb / inferred_xsec_fb if inferred_xsec_fb > 0 else np.nan

        norm_rows.append({
            "state": name,
            "normalization_scope": norm_scope_label,
            "target_visible_xsec_fb": target_visible_fb,
            "sum_weight_norm_scope_raw": sumw_norm,
            "inferred_xsec_norm_scope_fb": inferred_xsec_fb,
            "weight_scale": weight_scale,
            "br_ee": props["br_ee"],
            "target_inclusive_xsec_pb_if_visible": target_visible_fb / props["br_ee"] / 1000.0,
        })

        true_filter = make_filter_for_state(props, args, has_true_q2, reco=False)
        reco_filter = make_filter_for_state(props, args, has_true_q2, reco=True)

        truth = df.Filter(true_filter).AsNumpy(["true_W", "weight"])
        reco = df.Filter(reco_filter).AsNumpy(["reco_W", "weight"])

        if len(truth["true_W"]) == 0:
            print(f"Upsilon({name}) skipped because truth selection is empty")
            continue

        true_w = np.asarray(truth["weight"], dtype=float) * weight_scale
        true_W = np.asarray(truth["true_W"], dtype=float)

        if len(reco["reco_W"]) > 0:
            reco_w = np.asarray(reco["weight"], dtype=float) * weight_scale
            reco_W = np.asarray(reco["reco_W"], dtype=float)
        else:
            reco_w = np.array([], dtype=float)
            reco_W = np.array([], dtype=float)

        n_true, err_true = hist_weighted(true_W, true_w, w_bins)
        n_reco, err_reco = hist_weighted(reco_W, reco_w, w_bins)

        eff = np.divide(n_reco, n_true, out=np.zeros_like(n_true), where=n_true > 0.0)

        br_corr = 1.0
        if args.weight_convention == "visible_decay":
            br_corr = 1.0 / props["br_ee"]

        sigma = np.full_like(n_true, np.nan, dtype=float)
        err_sigma = np.full_like(n_true, np.nan, dtype=float)
        fluxes = np.zeros_like(n_true, dtype=float)

        for i in range(len(w_centers)):
            flux_int = integrated_flux(
                w_bins[i], w_bins[i + 1],
                args.electron_energy, args.proton_energy, args.q2_max
            )
            fluxes[i] = flux_int

            denom_truth = args.lumi_invfb * 1000.0 * flux_int
            if n_true[i] > 0.0 and denom_truth > 0.0:
                # This is equivalent to n_reco / efficiency when efficiency is nonzero,
                # but remains defined in bins with no reconstructed events.
                sigma[i] = n_true[i] / denom_truth * br_corr
                err_sigma[i] = err_true[i] / denom_truth * br_corr

        for i in range(len(w_centers)):
            all_rows.append({
                "state": name,
                "W_low": w_bins[i],
                "W_high": w_bins[i + 1],
                "W_center": w_centers[i],
                "flux_int": fluxes[i],
                "n_true_weighted_visible": n_true[i],
                "n_reco_weighted_visible": n_reco[i],
                "eff_reco_over_true": eff[i],
                "sigma_gamma_p_pb": sigma[i],
                "sigma_gamma_p_err_pb": err_sigma[i],
                "plotted": bool(np.isfinite(sigma[i]) and sigma[i] > 0.0 and eff[i] > args.min_eff),
            })

        mask = np.isfinite(sigma) & (sigma > 0.0) & (eff > args.min_eff)
        ax.errorbar(
            w_centers[mask], sigma[mask], yerr=err_sigma[mask],
            fmt="o", color=props["color"],
            label=rf"ePIC $\Upsilon({name})$, $\epsilon>{args.min_eff:g}$",
            markersize=7, linewidth=1.8,
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
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)

    print("=== Normalization summary ===")
    print(pd.DataFrame(norm_rows).to_string(index=False))
    print(f"\nSaved {out_png}")
    print(f"Saved {out_csv}")
    print(f"Saved {out_norm_csv}")


if __name__ == "__main__":
    main()
