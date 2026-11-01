#!/usr/bin/env python3
"""Render current single-state Upsilon samples with the analysis note's mass/efficiency layouts.

Mass pages are new signal-only DSCB fits to MC counts, subsequently scaled by
constant normalized tuple weights. Projected mode additionally reports a separate
signal-only Poisson Fisher forecast; efficiency errors remain finite-MC errors.
The figures do not reproduce the historical
three-state plus background model. Efficiency curves retain notebook 02 cell 3
definitions: acceptance A/G, conditional reconstruction (A and R)/A, and total
R/G. The nested (A and R)/G quantity is recorded separately in CSV/audit output.
Eligible stored tuple rows provide the denominator. No additional decay branching
fraction or cross-section scaling is applied here.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

# Shared explicit reconstructed-proton selection; truth denominators are unchanged.
from pathlib import Path as _ProtonPath
import sys as _proton_sys
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from plot_latest_sample_performance import (fit_mass, crystal_ball_integral, wilson,
                                            statistical_configuration, mass_projection_forecast,
                                            mass_interval_label, projected_yield_interval)
from projected_intervals import (add_interval_arguments, expected_count_intervals,
                                 draw_absolute_intervals)
from note_style import set_style, finish, SIGNAL_BLUE

MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
MASS = MASSES["1S"]  # Legacy import compatibility.
COLORS = ["dodgerblue", "darkorange"]
REQUIRED = ["weight", "true_mass", "true_Q2", "true_W", "true_y", "true_pT", "true_t",
            "true_ePlus_eta", "true_decay_eMinus_eta", "reco_mass", "reco_y", "reco_W", "sigma_epz"]


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_sample(path, tree_name, energy_label, electron_energy, w_min, w_max, q2_cut, state="1S"):
    mass = MASSES[state]
    with uproot.open(path) as source:
        tree = source[tree_name]
        required = REQUIRED + proton_branches()
        missing = [name for name in required if name not in tree.keys()]
        if missing:
            raise RuntimeError("Missing required tuple branches: " + ", ".join(missing))
        data = {name: tree[name].array(library="np") for name in required}
    weights = np.asarray(data["weight"], float)
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise RuntimeError("Normalized tuple weights must be finite and nonnegative")
    truth = (np.isfinite(data["true_mass"]) & (np.abs(data["true_mass"] - mass) <= 0.2)
             & np.isfinite(data["true_Q2"]) & (data["true_Q2"] >= 0) & (data["true_Q2"] <= q2_cut)
             & np.isfinite(data["true_W"]) & (data["true_W"] >= w_min) & (data["true_W"] <= w_max))
    accepted = (truth & np.isfinite(data["true_ePlus_eta"]) & np.isfinite(data["true_decay_eMinus_eta"])
                & (np.abs(data["true_ePlus_eta"]) < 3.5) & (np.abs(data["true_decay_eMinus_eta"]) < 3.5))
    reco = (truth & np.isfinite(data["reco_mass"]) & (data["reco_mass"] >= 8) & (data["reco_mass"] <= 12)
            & np.isfinite(data["reco_W"]) & (data["reco_W"] >= w_min) & (data["reco_W"] <= w_max)
            & np.isfinite(data["sigma_epz"]) & (np.abs(data["sigma_epz"] - 2 * electron_energy) <= 0.2 * electron_energy))
    reco &= proton_tag_mask(data)
    if not np.any(truth):
        raise RuntimeError(f"No eligible {state} truth rows for {energy_label}")
    audit = {"proton_selection": proton_selection_metadata(),
             "truth_rows_with_proton_tag": int(np.sum(truth & proton_tag_mask(data))),
             "input": str(Path(path).resolve()), "tree": tree_name, "state": state, "energy_label": energy_label,
             "electron_energy_GeV": electron_energy, "tuple_rows": len(weights), "eligible_truth_rows": int(truth.sum()),
             "truth_mass_GeV": [mass - 0.2, mass + 0.2], "true_Q2_GeV2": [0, q2_cut],
             "true_W_GeV": [w_min, w_max], "reco_W_GeV": [w_min, w_max], "reco_mass_GeV": [8, 12],
             "sigma_epz_GeV": [1.8 * electron_energy, 2.2 * electron_energy],
             "decay_lepton_acceptance": "both abs(true decay-lepton eta) < 3.5",
             "truth_in_acceptance_rows": int(accepted.sum()), "selected_reco_rows": int(reco.sum()),
             "selected_reco_in_acceptance_rows": int((accepted & reco).sum()),
             "selected_reco_outside_acceptance_rows": int((reco & ~accepted).sum()),
             "efficiency_denominator": "Eligible stored UpsilonTree rows, not an independently validated generated-event denominator",
             "efficiency_definitions": {"acceptance": "N(A)/N(G)", "reconstruction": "N(A and R)/N(A)",
                                        "total": "N(R)/N(G)", "nested_total_audit_only": "N(A and R)/N(G)"},
             "efficiency_source_convention": "notebooks/02_efficiency_and_two_energy_studies.ipynb cell 3",
             "acceptance_times_reconstruction_identity": "A * epsilon_reco = nested_total; total additionally includes selected events outside truth lepton acceptance",
             "total_minus_acceptance_times_reconstruction": float((reco & ~accepted).sum() / truth.sum()),
             "total_minus_acceptance_times_reconstruction_percentage_points": float(100 * (reco & ~accepted).sum() / truth.sum()),
             "mass_selection": "R (does not additionally require truth lepton acceptance)",
             "normalization": "weight branch used directly; no additional luminosity or branching-fraction scaling"}
    return {"data": data, "truth": truth, "accepted": accepted, "reco": reco,
            "energy_label": energy_label, "electron_energy": electron_energy,
            "w_min": w_min, "w_max": w_max, "q2_cut": q2_cut, "state": state, "audit": audit}


def smooth_model(mass_points, pars, yield_in_range, bin_width):
    """Expected counts per fixed-width bin from the smooth, range-normalized DSCB."""
    mean, sigma, al, nl, ar, nr = pars
    z = (mass_points - mean) / sigma
    density = np.empty_like(z)
    left = z < -al
    right = z > ar
    core = ~(left | right)
    density[core] = np.exp(-0.5 * z[core]**2)
    density[left] = np.exp(nl * np.log(nl / al) - al**2 / 2 - nl * np.log(nl / al - al - z[left]))
    density[right] = np.exp(nr * np.log(nr / ar) - ar**2 / 2 - nr * np.log(nr / ar - ar + z[right]))
    integral = np.diff(crystal_ball_integral((np.array([8.0, 12.0]) - mean) / sigma, al, nl, ar, nr))[0] * sigma
    return density * yield_in_range * bin_width / integral


def clean_energy_label(label):
    return label.replace("x", r"$\times$")


def render_mass_pages(sample, output, stat):
    lumi = stat["active_luminosity_invfb"]
    data = sample["data"]
    settings = [("mass_fit_inclusive", "Inclusive reconstructed rapidity", None),
                ("mass_fit_rapidity_minus1_to_0", r"$-1 \leq y_{\rm reco} < 0$", (-1, 0)),
                ("mass_fit_rapidity_0_to_1", r"$0 \leq y_{\rm reco} < 1$", (0, 1))]
    mass_edges = np.linspace(8, 12, 201)
    centers = (mass_edges[:-1] + mass_edges[1:]) / 2
    fit_parameters, bins_rows, fit_summaries = [], [], {}
    for stem, region_label, region in settings:
        selected = sample["reco"].copy()
        if region is not None:
            selected &= np.isfinite(data["reco_y"]) & (data["reco_y"] >= region[0]) & (data["reco_y"] < region[1])
        mass = data["reco_mass"][selected]
        weights = data["weight"][selected] * stat["yield_weight_scale"]
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_xlabel(r"Reco $M_{e^+e^-}$ [GeV/$c^2$]")
        ax.set_ylabel("Projected events / 20 MeV")
        ax.set_xlim(8, 12)
        ax.grid(alpha=0.16)
        projected, _ = np.histogram(mass, mass_edges, weights=weights)
        sumw2, _ = np.histogram(mass, mass_edges, weights=weights**2)
        poisson_variance = projected.copy()
        shown_variance = poisson_variance if stat["stat_mode"] == "projected" else sumw2
        interval_settings = stat["count_interval"]
        intervals = (expected_count_intervals(projected, coverage=interval_settings["coverage"], method=interval_settings["method"])
                     if stat["stat_mode"] == "projected" else None)
        yield_interval = (projected_yield_interval(weights.sum(), interval_settings["coverage"], interval_settings["method"])
                          if stat["stat_mode"] == "projected" else None)
        if intervals is not None:
            draw_absolute_intervals(ax, centers, projected, intervals["lower"], intervals["upper"],
                                    fmt="o", ms=3, color="black", elinewidth=0.8,
                                    label="Projected expected mean", zorder=5)
        else:
            ax.errorbar(centers, projected, yerr=np.sqrt(shown_variance), fmt="o", ms=3, color="black", elinewidth=0.8,
                        label="Signal simulation", zorder=5)
        parameter_note = "Insufficient selected MC events for a six-parameter fit"
        if len(mass) >= 100 and weights.sum() > 0:
            if not np.allclose(weights, weights[0], rtol=1e-8, atol=1e-14):
                raise RuntimeError("Mass fit-to-projection scaling requires constant normalized event weights")
            counts, fitted_counts, summary = fit_mass(mass, mass_edges, state=sample.get("state", "1S"))
            forecast = (mass_projection_forecast(mass_edges, summary, weights.sum(), stat["projection_shape_mode"])
                        if stat["stat_mode"] == "projected" else {"status": "not_requested", "valid_local_covariance": False})
            summary["projected_forecast"] = forecast
            summary["statistical_configuration"] = stat
            summary["projected_yield_interval"] = yield_interval
            params = summary["parameters"]
            pars = [params[key]["value"] for key in ["mean_GeV", "sigma_core_GeV", "alpha_left", "n_left", "alpha_right", "n_right"]]
            x = np.linspace(8, 12, 2401)
            y = smooth_model(x, pars, float(weights.sum()), 0.02)
            ax.fill_between(x, 0, y, color=SIGNAL_BLUE, alpha=0.28, label=rf"$\Upsilon({sample.get('state', '1S')})$ signal")
            ax.plot(x, y, color=SIGNAL_BLUE, ls="--", lw=1.5)
            ax.plot(x, y, color="red", lw=2, label="DSCB model")
            residual = mass - data["true_mass"][selected]
            quantiles = np.quantile(residual, [0.16, 0.5, 0.84])
            summary.update({"reco_rapidity_region": list(region) if region else None,
                            "selected_projected_yield": float(weights.sum()), "selected_mc_count": len(mass),
                            "mass_residual_quantiles_GeV": dict(zip(["q16", "median", "q84"], map(float, quantiles))),
                            "mass_residual_central68_halfwidth_GeV": float((quantiles[2] - quantiles[0]) / 2),
                            "fit_model_scope": f"One signal-only {sample.get('state', '1S')} DSCB, no fitted background or other-state components",
                            "plot_uncertainties": mass_interval_label(interval_settings["method"], interval_settings["coverage"]) if stat["stat_mode"] == "projected" else "sqrt(sum of squared event weights), finite MC uncertainty",
                            "fit_uncertainties": "MC calibration and projected data Fisher covariance are reported separately"})
            if stat["stat_mode"] == "projected" and forecast.get("valid_local_covariance"):
                error_mean = forecast["parameters"]["mean_GeV"]["projected_stat_error"]
                error_sigma = forecast["parameters"]["sigma_core_GeV"]["projected_stat_error"]
                yield_low = weights.sum() - yield_interval["lower"]
                yield_high = yield_interval["upper"] - weights.sum()
                parameter_note = (rf"$m={pars[0]:.4f}\pm{error_mean:.4f}$ GeV, "
                                  rf"$\sigma_{{\rm core}}={1000 * pars[1]:.1f}\pm{1000 * error_sigma:.1f}$ MeV, "
                                  rf"$N={weights.sum():.1f}^{{+{yield_high:.1f}}}_{{-{yield_low:.1f}}}$" + "\n" +
                                  f"Shape errors: local Fisher ({stat['projection_shape_mode']}); yield bounds: {interval_settings['method'].replace('_', ' ')}")
            else:
                parameter_note = (rf"$m={pars[0]:.4f}$ GeV, $\sigma_{{\rm core}}={1000 * pars[1]:.1f}$ MeV, "
                                  rf"$N_{{\rm projected}}={weights.sum():.1f}$  |  MC deviance/ndf = {summary['poisson_deviance']:.1f}/{summary['ndof']}")
                if stat["stat_mode"] == "projected":
                    parameter_note += "\nProjected covariance unavailable: " + forecast["status"]
            if summary["shape_mismatch_flag"]:
                ax.text(0.98, 0.55, "DSCB is a descriptive model\nTail mismatch remains", ha="right", va="top",
                        transform=ax.transAxes, fontsize=10, color="0.35")
            if forecast.get("low_expected_statistics_flag"):
                ax.text(0.98, 0.33, "Low expected yield\nLocal shape approximation may be inaccurate", ha="right", va="top",
                        transform=ax.transAxes, fontsize=9, color="0.4")
            for name, parameter in params.items():
                fit_parameters.append({"region": stem, "parameter": name, "value": parameter["value"],
                                       "mc_stat_error": parameter["mc_stat_error"],
                                       "projected_stat_error": forecast.get("parameters", {}).get(name, {}).get("projected_stat_error"),
                                       "projected_parameter_fixed": forecast.get("parameters", {}).get(name, {}).get("fixed"),
                                       "selected_stat_error": forecast.get("parameters", {}).get(name, {}).get("projected_stat_error") if stat["stat_mode"] == "projected" else parameter["mc_stat_error"],
                                       "at_bound": parameter["at_bound"], "stat_mode": stat["stat_mode"]})
        else:
            counts, _ = np.histogram(mass, mass_edges)
            fitted_counts = np.full_like(counts, np.nan, dtype=float)
            forecast = {"status": "insufficient_mc_calibration", "valid_local_covariance": False}
            summary = {"statistical_configuration": stat, "projected_forecast": forecast, "projected_yield_interval": yield_interval, "status": "insufficient_mc_statistics", "selected_mc_count": len(mass),
                       "selected_projected_yield": float(weights.sum()), "reco_rapidity_region": list(region) if region else None}
            ax.text(0.5, 0.5, "Insufficient selected events for a mass fit", transform=ax.transAxes, ha="center")
        fit_summaries[stem] = summary
        peak_with_error = float(np.max(intervals["upper"] if intervals is not None else projected + np.sqrt(shown_variance))) if len(projected) else 0
        ax.set_ylim(0, max(peak_with_error * 1.22, 1e-3))
        ax.legend(loc="upper right", frameon=False, fontsize=10)
        ax.text(0.98, 0.45, (mass_interval_label(interval_settings["method"], interval_settings["coverage"]) + "\nMC calibration remains separate") if stat["stat_mode"] == "projected" else "Error bars: finite MC statistics", ha="right", va="top", transform=ax.transAxes,
                fontsize=9, color="0.4")
        subtitle = (f"ePIC performance  |  {clean_energy_label(sample['energy_label'])}  |  "
                    rf"$\mathcal{{L}}={lumi:g}$ fb$^{{-1}}$  |  {region_label}" + "\n" + parameter_note)
        finish(fig, output, stem, rf"$\Upsilon({sample.get('state', '1S')})\rightarrow e^+e^-$ invariant-mass fit", subtitle)
        for i in range(len(counts)):
            bins_rows.append({"region": stem, "low_GeV": float(mass_edges[i]), "high_GeV": float(mass_edges[i+1]),
                              "mc_count": int(counts[i]), "model_mc_expectation": float(fitted_counts[i]) if np.isfinite(fitted_counts[i]) else None,
                              "projected_events": float(projected[i]), "projected_mc_variance": float(sumw2[i]),
                              "projected_data_variance": float(poisson_variance[i]), "displayed_variance": float(shown_variance[i]),
                              "projected_interval_low": float(intervals["lower"][i]) if intervals is not None else None,
                              "projected_interval_high": float(intervals["upper"][i]) if intervals is not None else None,
                              "interval_reference_count": float(intervals["reference_count"][i]) if intervals is not None and intervals["reference_count"] is not None and np.isfinite(intervals["reference_count"][i]) else None,
                              "count_interval_method": interval_settings["method"] if intervals is not None else "not_applied_mc",
                              "interval_coverage": interval_settings["coverage"] if intervals is not None else None,
                              "stat_mode": stat["stat_mode"], "mc_effective_entries": float(projected[i]**2/sumw2[i]) if sumw2[i] > 0 else None,
                              "model_asimov_expectation": float(forecast["asimov_expected_bin_counts"][i]) if "asimov_expected_bin_counts" in forecast else None})
    write_csv(output / "mass_fit_parameters.csv", fit_parameters)
    write_csv(output / "mass_fit_bins.csv", bins_rows)
    (output / "mass_fit_quality.json").write_text(json.dumps(fit_summaries, indent=2) + "\n")
    return fit_summaries


def render_efficiencies(samples, output):
    w_min = min(sample["w_min"] for sample in samples)
    w_max = max(sample["w_max"] for sample in samples)
    variables = [("true_pT", np.linspace(0, 2.5, 25), r"Generator-level $\Upsilon$ $p_T$ [GeV/$c$]"),
                 ("true_y", np.linspace(-1.5, 3, 24), r"Generator-level $\Upsilon$ rapidity $y_\Upsilon$"),
                 ("true_W", np.linspace(w_min, w_max, max(6, int(np.ceil((w_max - w_min) / 5))) + 1), r"Generator-level $W_{\gamma p}$ [GeV]"),
                 ("true_t", np.linspace(0, 1.5, 15), r"Generator-level $|t|$ [GeV$^2$]")]
    fig, axes = plt.subplots(2, 2, figsize=(14, 11.5))
    handles, labels, rows, integrated = [], [], [], []
    for index, sample in enumerate(samples):
        data, gen, accepted, reco = sample["data"], sample["truth"], sample["accepted"], sample["reco"]
        categories = [("acceptance", accepted, gen, "o", "-", "$A$"),
                      ("reconstruction", accepted & reco, accepted, "s", "--", r"$\epsilon_{\rm reco}$"),
                      ("total", reco, gen, "^", ":", "total"),
                      ("nested_total_audit_only", accepted & reco, gen, None, None, "nested total (audit only)")]
        sample_audit = sample["audit"]
        sample_audit["efficiency_variable_coverage"] = {}
        for kind, passed, denominator_mask, marker, linestyle, symbol in categories:
            n_total, n_pass = int(denominator_mask.sum()), int(passed.sum())
            mean, lo, hi = wilson(np.array([n_pass]), np.array([n_total]))
            integrated.append({"sample": sample["energy_label"], "efficiency_type": kind, "numerator": n_pass,
                               "denominator": n_total, "efficiency": float(mean[0]) if n_total else None,
                               "wilson_low": float(lo[0]) if n_total else None, "wilson_high": float(hi[0]) if n_total else None})
            label = f"{clean_energy_label(sample['energy_label'])} {symbol} [{100 * mean[0]:.1f}%]"
            for panel_index, (ax, (branch, edges, xlabel)) in enumerate(zip(axes.flat, variables)):
                finite = np.isfinite(data[branch]) & (data[branch] > -900)
                total, _ = np.histogram(data[branch][denominator_mask & finite], edges)
                numerator, _ = np.histogram(data[branch][passed & finite], edges)
                value, low, high = wilson(numerator, total)
                populated = total > 0
                centers = (edges[:-1] + edges[1:]) / 2
                if marker is not None:
                    line = ax.errorbar(centers[populated], value[populated],
                                      yerr=[np.maximum(value[populated] - low[populated], 0), np.maximum(high[populated] - value[populated], 0)],
                                      color=COLORS[index], marker=marker, ls=linestyle, ms=3.5, lw=1.3, capsize=2,
                                      label=label)
                    if panel_index == 0:
                        handles.append(line)
                        labels.append(label)
                for i in range(len(total)):
                    rows.append({"sample": sample["energy_label"], "variable": branch, "efficiency_type": kind,
                                 "low": float(edges[i]), "high": float(edges[i+1]), "numerator": int(numerator[i]),
                                 "denominator": int(total[i]), "efficiency": float(value[i]) if total[i] else None,
                                 "wilson_low": float(low[i]) if total[i] else None, "wilson_high": float(high[i]) if total[i] else None})
                sample_audit["efficiency_variable_coverage"].setdefault(branch, {})[kind] = {
                    "denominator_below_plot": int(np.sum(denominator_mask & finite & (data[branch] < edges[0]))),
                    "denominator_above_plot": int(np.sum(denominator_mask & finite & (data[branch] > edges[-1]))),
                    "denominator_invalid": int(np.sum(denominator_mask & ~finite))}
    for ax, (branch, edges, xlabel) in zip(axes.flat, variables):
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Efficiency")
        ax.set_ylim(-0.025, 1.45)
        ax.set_xlim(edges[0], edges[-1])
        ax.axhline(1, color="gray", ls="--", lw=1)
        ax.grid(alpha=0.16)
    cut_text = "; ".join(rf"{clean_energy_label(s['energy_label'])}: ${s['w_min']:g} \leq W \leq {s['w_max']:g}$ GeV" for s in samples)
    subtitle = (r"Decay-lepton acceptance $|\eta_{e^+}|,|\eta_{e^-}|<3.5$  |  Truth kinematics; eligible tuple rows" + "\n" +
                rf"$0\leq Q^2\leq {samples[0]['q2_cut']:g}$ GeV$^2$  |  " + cut_text)
    finish(fig, output, f"efficiency_{samples[0].get('state', '1S')}_truth_rapidity",
           rf"$\Upsilon({samples[0].get('state', '1S')})$ acceptance and reconstruction", subtitle,
           handles=handles, labels=labels, legend_ncol=3)
    write_csv(output / "efficiency_points.csv", rows)
    write_csv(output / "efficiency_integrated.csv", integrated)
    return integrated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", choices=list(MASSES), default="1S")
    parser.add_argument("--input", required=True)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--energy-label", required=True)
    parser.add_argument("--electron-energy", type=float, default=9)
    parser.add_argument("--w-min", type=float, default=20)
    parser.add_argument("--w-max", type=float, default=70)
    parser.add_argument("--q2-cut", type=float, default=1)
    parser.add_argument("--lumi-invfb", type=float, default=10, help="Luminosity already represented by input tuple weights")
    parser.add_argument("--stat-mode", choices=["mc", "projected"], default="mc")
    parser.add_argument("--projection-lumi-invfb", type=float)
    parser.add_argument("--projection-shape-mode", choices=["floating", "fixed_tails"], default="floating")
    add_interval_arguments(parser)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--second-input")
    parser.add_argument("--second-energy-label", default="9x275 GeV")
    parser.add_argument("--second-electron-energy", type=float, default=9)
    parser.add_argument("--second-w-min", type=float, default=20)
    parser.add_argument("--second-w-max", type=float, default=100)
    args = parser.parse_args()
    if args.w_min >= args.w_max or args.q2_cut <= 0 or args.lumi_invfb <= 0 or args.electron_energy <= 0:
        parser.error("Require an ordered W interval and positive Q2 limit, luminosity, and electron energy")
    if args.second_input and (args.second_w_min >= args.second_w_max or args.second_electron_energy <= 0):
        parser.error("Second sample requires an ordered W interval and positive electron energy")
    try:
        stat = statistical_configuration(args.stat_mode, args.lumi_invfb, args.projection_lumi_invfb, args.projection_shape_mode)
        expected_count_intervals(np.array([0.]), coverage=args.interval_coverage, method=args.count_interval)
    except ValueError as error:
        parser.error(str(error))
    output = Path(args.out_dir)
    stat["count_interval"] = {"method": args.count_interval, "coverage": args.interval_coverage,
                              "applied": args.stat_mode == "projected", "covariance_changed": False}
    output.mkdir(parents=True, exist_ok=True)
    set_style()
    sample = load_sample(args.input, args.tree, args.energy_label, args.electron_energy, args.w_min, args.w_max, args.q2_cut, args.state)
    samples = [sample]
    if args.second_input:
        samples.append(load_sample(args.second_input, args.tree, args.second_energy_label, args.second_electron_energy,
                                   args.second_w_min, args.second_w_max, args.q2_cut, args.state))
    fit_summaries = render_mass_pages(sample, output, stat)
    efficiency_summary = render_efficiencies(samples, output)
    audit = {"state": args.state, "purpose": f"Current {args.state} results rendered with reviewed analysis-note layouts",
             "statistical_configuration": stat,
             "normalization_caveat": "Input weights must represent --lumi-invfb; projected yields scale once by target/source luminosity",
             "empty_bin_caveat": "Zero empirical template counts do not establish a precisely known zero physical expectation; DSCB Asimov counts are separate",
             "interval_display_note": "Projected count-interval endpoints do not alter expected means, variances, or local Fisher covariance; zero-mean bins remain at zero with positive median-Garwood upper endpoints",
             "mass_fit_scope": f"First input only, new signal-only {args.state} fit; historical three-state plus exponential model is not reproduced",
             "efficiency_intervals": "68.27% Wilson intervals from unweighted MC counts",
             "samples": [s["audit"] for s in samples],
             "outputs": ["mass_fit_inclusive", "mass_fit_rapidity_minus1_to_0", "mass_fit_rapidity_0_to_1", f"efficiency_{args.state}_truth_rapidity"]}
    (output / "note_mass_efficiency_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"output_dir": str(output), "mass_regions": {key: {"selected_mc_count": value["selected_mc_count"],
           "selected_projected_yield": value["selected_projected_yield"], "converged": value.get("converged", False)} for key, value in fit_summaries.items()},
           "efficiencies": efficiency_summary}, indent=2))


if __name__ == "__main__":
    main()
