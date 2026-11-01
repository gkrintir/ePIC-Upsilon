#!/usr/bin/env python3
"""Audit projected counts and binning using saved CSV/JSON files, without ROOT.

This checks numerical bookkeeping. It does not certify that sparse bins support
a simultaneous mass fit or that different extraction selections are identical.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import poisson

from settings import load_profile, make_config, new_output_dir
from upsilon_esr.common import (
    boolean, expected_intervals, project_cross_section, read_csv, read_json,
    resolve_artifacts, sample_metadata, source_luminosity, write_json,
)
from upsilon_esr.w_cross_section import _check_normalization


def _check(report, state, observable, name, status, **evidence):
    report["checks"].append(dict(state=state, observable=observable,
                                 check=name, status=status, evidence=evidence))


def _error(report, state, observable, error):
    report["errors"].append(dict(state=state, observable=observable,
                                 error_type=type(error).__name__, message=str(error)))


def audit_bin_counts(frame, *, state, observable, target_lumi, coverage,
                     pt2_display_range=(0., 1.)):
    """Independent count audit; preserve bad rows as diagnostics, never repair.

    Count intervals are median-Garwood forecast endpoints for future integer
    counts. They are not an interval constructed from a fractional observation.
    """
    is_w = observable == "W"
    lo_key, hi_key, center_key = (("W_low", "W_high", "W_center") if is_w
                                  else ("x_low", "x_high", "x_center"))
    count_key = "n_reco_weighted_visible" if is_w else "n_reco_visible"
    value_key = "sigma_gamma_p_pb" if is_w else "dsigma_dx_pb"
    factor_key = "reco_count_to_sigma_pb" if is_w else "reco_count_to_dsigma_dx_pb"
    f = frame.sort_values(lo_key).reset_index(drop=True)
    source_lumi = source_luminosity(f)
    ratio = float(target_lumi / source_lumi)
    count = pd.to_numeric(f[count_key], errors="coerce").to_numpy(float)
    central = pd.to_numeric(f[value_key], errors="coerce").to_numpy(float)
    lo, hi, center = [pd.to_numeric(f[k], errors="coerce").to_numpy(float)
                      for k in (lo_key, hi_key, center_key)]
    explicit_factor = factor_key in f
    source_factor = (pd.to_numeric(f[factor_key], errors="coerce").to_numpy(float)
                     if explicit_factor else
                     np.divide(central, count, out=np.full(len(f), np.nan), where=count > 0))
    mu, factor = count * ratio, source_factor / ratio
    count_ok = np.isfinite(mu) & (mu >= 0)
    supported = (count_ok & np.isfinite(factor) & (factor > 0)
                 & np.isfinite(central) & (central >= 0))
    closure = supported & np.isclose(central, mu * factor, rtol=2e-6, atol=1e-9)
    upstream = np.ones(len(f), dtype=bool)
    flag_valid = np.ones(len(f), dtype=bool)
    for key in ("plotted", "passes_min_eff"):
        if key in f:
            flag_valid &= f[key].astype(str).str.lower().isin(("true", "false", "1", "0", "1.0", "0.0")).to_numpy()
            upstream &= boolean(f[key])
    valid = supported & upstream & closure & flag_valid
    displayed = valid.copy()
    if observable in ("W", "pT2"):
        displayed &= central > 0
    in_fit_range = ((lo >= pt2_display_range[0] - 1e-10) & (hi <= pt2_display_range[1] + 1e-10)
                    if observable == "pT2" else np.zeros(len(f), dtype=bool))
    if observable == "pT2":
        displayed &= in_fit_range
    median, lower, upper, p0, p_lt5 = [np.full(len(f), np.nan) for _ in range(5)]
    if count_ok.any():
        interval = expected_intervals(mu[count_ok], coverage)
        median[count_ok] = interval["reference_count"]
        lower[count_ok], upper[count_ok] = interval["lower"], interval["upper"]
        p0[count_ok], p_lt5[count_ok] = np.exp(-mu[count_ok]), poisson.cdf(4, mu[count_ok])
    out = pd.DataFrame(dict(
        state=state, observable=observable, bin_index=np.arange(len(f)),
        x_low=lo, x_high=hi, x_center=center, bin_width=hi-lo,
        source_lumi_invfb=source_lumi, target_lumi_invfb=target_lumi,
        source_expected_reco_count=count, expected_reco_count=mu,
        count_interval_reference=median, count_interval_low=lower,
        count_interval_high=upper, probability_zero=p0,
        probability_below_5=p_lt5, supported=supported,
        upstream_selected=upstream, upstream_flags_valid=flag_valid,
        valid=valid, displayed=displayed, used_by_slope_fit=valid & in_fit_range,
        count_factor_source="stored" if explicit_factor else "inferred_from_central_over_count",
        count_factor_closure=np.where(~supported, "UNSUPPORTED", np.where(closure, "PASS" if explicit_factor else "DERIVED_NOT_INDEPENDENT", "FAIL")),
        central_cross_section=central, count_to_cross_section=factor,
        reconstructed_central_cross_section=mu*factor,
        lower_cross_section=np.where(supported, lower*factor, np.nan),
        upper_cross_section=np.where(supported, upper*factor, np.nan),
        cross_section_units="pb" if is_w else "pb per unit of observable",
    ))
    return out


def _summarize_bins(frame, mass_yield):
    values = frame.expected_reco_count.to_numpy(float)
    finite = np.isfinite(values) & (values >= 0)
    summary = {"n_input_bins": len(frame), "count_rows_complete": bool(finite.all()),
               "cross_observable_selection_match": "NOT_VERIFIED"}
    for name, mask in (("full", np.ones(len(frame), bool)),
                       ("valid", frame.valid.to_numpy(bool)),
                       ("displayed", frame.displayed.to_numpy(bool))):
        selected = values[mask & finite]
        total = float(selected.sum())
        summary[f"{name}_expected_reco_sum"] = total
        summary[f"{name}_n_bins"] = int(mask.sum())
        summary[f"{name}_minimum_expected_count"] = float(selected.min()) if len(selected) else None
        for threshold in (1, 3, 5):
            summary[f"{name}_n_bins_mu_lt_{threshold}"] = int((selected < threshold).sum())
        summary[f"{name}_to_mass_expected_yield_ratio"] = (
            total / mass_yield if mass_yield is not None and mass_yield > 0 else None)
    summary["sum_interpretation"] = (
        "Counts are summed over the listed bins only; invalid counts are omitted from the sum and flagged. "
        "Ratios to the mass expected selected yield are diagnostics. A difference is not proof of a "
        "normalization error without independently matched selections, ranges and overflow accounting.")
    if len(frame) and frame.observable.iloc[0] == "pT2":
        selected = frame.used_by_slope_fit.to_numpy(bool) & finite
        fit_count = float(values[selected].sum())
        summary["slope_fit_count_diagnostic"] = dict(
            n_supported_complete_bins=int(selected.sum()),
            has_minimum_three_supported_bins=bool(selected.sum() >= 3),
            expected_reconstructed_count=fit_count,
            probability_zero_events=float(np.exp(-fit_count)),
            interval_coverage_verified=False,
            interpretation="Eligibility/count diagnostic only; neither exact interval coverage nor mass-fit feasibility is established.")
    return summary


def _audit_geometry(report, state, observable, frame):
    lo, hi, center = [frame[k].to_numpy(float) for k in ("x_low", "x_high", "x_center")]
    finite = np.isfinite(lo) & np.isfinite(hi) & np.isfinite(center)
    good = bool(finite.all() and np.all(hi > lo) and
                np.all(lo[1:] >= hi[:-1] - 1e-10) and
                np.all((center >= lo) & (center <= hi)))
    _check(report, state, observable, "finite_nonoverlapping_bins_and_centers", "PASS" if good else "FAIL")
    if observable != "W":
        midpoint = bool(np.allclose(center, (lo+hi)/2, rtol=1e-8, atol=1e-10))
        _check(report, state, observable, "center_matches_symmetric_horizontal_errors",
               "PASS" if midpoint else "FAIL")


def audit_projected_binning(config, out=None, *, pt2_display_range=(0., 1.), write_outputs=True):
    """Return and optionally save a bookkeeping report and per-bin count table.

    Errors in one observable are recorded while the remaining inputs are
    inspected. No binning or normalization is modified. No overall binning PASS
    is issued because background, state separation and migration are not tested.
    """
    fit_range = tuple(map(float, pt2_display_range))
    if len(fit_range) != 2 or not np.isfinite(fit_range).all() or not 0 <= fit_range[0] < fit_range[1]:
        raise ValueError("pt2_display_range must be a finite increasing nonnegative pair")
    report = dict(status="DIAGNOSTICS_ONLY", energy=config.energy,
                  target_lumi_invfb=float(config.lumi_invfb), coverage=float(config.coverage),
                  interval_method="median_garwood", pt2_display_range=list(fit_range),
                  cross_observable_selection_match="NOT_VERIFIED", checks=[], errors=[], states={},
                  limitations=[
                      "This is a count and algebra audit, not proof that the chosen bins support a measurement.",
                      "Signal-only reconstructed counts and fixed corrections omit background, mass separation, migration and calibration uncertainty.",
                      "Median-Garwood endpoints are expected endpoints for future integer observations, not fractional-data confidence intervals.",
                      "Different observables can select different events. No common-selection evidence is inferred from equal sums or matching Q2/W labels.",
                      "No automatic merging is performed; any physics rebinning must recompute corrections/flux and uncertainties.",
                  ])
    all_bins = []
    for state in ("1S", "2S", "3S"):
        state_info = dict(mass_expected_selected_yield=None, input_metadata={}, observables={})
        report["states"][state] = state_info
        try:
            base = resolve_artifacts(config.runs[state])
        except (OSError, ValueError) as error:
            _error(report, state, "all", error)
            continue
        performance, calibration = None, None
        try:
            path = base / "performance" / config.energy
            performance = read_json(path / "performance_audit.json")
            calibration = read_json(path / "mass_fit.json")
            state_info["input_metadata"]["performance_audit"] = performance
            state_info["input_metadata"]["mass_fit"] = calibration
            mass_lumi = source_luminosity(pd.DataFrame(), calibration)
            perf_lumi = source_luminosity(pd.DataFrame(), performance)
            same_lumi = bool(np.isclose(mass_lumi, perf_lumi, rtol=1e-8, atol=0))
            _check(report, state, "mass", "mass_fit_vs_performance_source_lumi",
                   "PASS" if same_lumi else "FAIL", mass_fit=mass_lumi, performance=perf_lumi)
            state_info["mass_source_lumi_invfb"] = mass_lumi
            state_info["performance_source_lumi_invfb"] = perf_lumi
            mass = read_csv(path / "mass_bins.csv", required=("low", "high", "projected_yield"))
            mass_low = mass.low.to_numpy(float)
            mass_high = mass.high.to_numpy(float)
            mass_geometry_ok = bool(len(mass) and np.isfinite(mass_low).all() and
                                    np.isfinite(mass_high).all() and np.all(mass_high > mass_low) and
                                    np.allclose(mass_low[1:], mass_high[:-1], rtol=1e-8, atol=1e-10))
            _check(report, state, "mass", "finite_positive_width_contiguous_mass_bins",
                   "PASS" if mass_geometry_ok else "FAIL")
            if not mass_geometry_ok:
                raise ValueError("Mass histogram edges must be finite, increasing and contiguous")
            vals = mass.projected_yield.to_numpy(float)
            if not np.isfinite(vals).all() or np.any(vals < 0):
                raise ValueError("Mass histogram has invalid expected yields")
            expected = float(calibration["selected_projected_yield"])
            if not np.isfinite(expected) or expected < 0:
                raise ValueError("Saved selected mass yield must be finite and nonnegative")
            _check(report, state, "mass", "mass_histogram_vs_selected_yield",
                   "PASS" if np.isclose(vals.sum(), expected, rtol=1e-7, atol=1e-10) else "FAIL",
                   histogram_sum=float(vals.sum()), saved_selected_yield=expected)
            state_info["mass_expected_selected_yield"] = expected * config.lumi_invfb / mass_lumi
            state_info["mass_histogram_expected_sum"] = float(vals.sum()) * config.lumi_invfb / mass_lumi
            state_info["mass_interval"] = [float(mass.low.min()), float(mass.high.max())]
            state_info["mass_yield_definition"] = "Saved selected expected yield, not a fitted component yield; scaled to target luminosity."
            state_info["mass_reference_metadata_consistent"] = same_lumi
            sample_metadata(config, state)
        except (OSError, ValueError, KeyError, TypeError) as error:
            _error(report, state, "mass", error)

        for observable in ("W", "pT2", "y"):
            try:
                if observable == "W":
                    file = base / "cross_sections" / f"{config.energy}_sigma_gamma_p_{state.lower()}.csv"
                    raw = read_csv(file, state=state, required=("W_low", "W_high", "W_center", "n_reco_weighted_visible", "sigma_gamma_p_pb"))
                    norm = read_csv(file.with_name(file.stem + "_normalization.csv"), state=state, required=("br_ee",))
                    state_info["input_metadata"]["w_normalization"] = norm.to_dict("records")
                    # Older W normalization tables do not store their event/flux
                    # cuts. Algebraic normalization closure cannot establish them.
                    phase_keys = ("q2_cut_events", "q2_max_flux", "w_min", "w_max")
                    if len(norm) == 1 and all(key in norm for key in phase_keys) and performance is not None:
                        n = norm.iloc[0]
                        truth = performance["truth_selection"]
                        matching = bool(np.isclose(n.q2_cut_events, n.q2_max_flux) and
                                        np.isclose(n.q2_cut_events, truth["Q2_GeV2"][1]) and
                                        np.allclose([n.w_min, n.w_max], truth["W_GeV"]))
                        _check(report, state, observable, "recorded_Q2_W_vs_performance",
                               "PASS" if matching else "FAIL")
                    else:
                        _check(report, state, observable, "recorded_Q2_W_vs_performance", "NOT_VERIFIED",
                               missing_normalization_columns=[key for key in phase_keys if key not in norm],
                               interpretation="Saved flux/count algebra does not independently verify the W event and flux selections.")
                    try:
                        _check_normalization(raw, norm, state)
                        _check(report, state, observable, "saved_truth_flux_normalization", "PASS")
                    except (ValueError, KeyError, TypeError) as error:
                        _error(report, state, observable, error)
                    central_key, count_key, factor_key = "sigma_gamma_p_pb", "n_reco_weighted_visible", "reco_count_to_sigma_pb"
                else:
                    stem = f"{config.energy}_{state.lower()}"
                    file = base / "differential" / f"{stem}_{state}_{observable}.csv"
                    raw = read_csv(file, state=state, required=("x_low", "x_high", "x_center", "bin_width", "n_reco_visible", "dsigma_dx_pb"))
                    meta = read_csv(base / "differential" / f"{stem}_metadata.csv", state=state)
                    state_info["input_metadata"]["differential_metadata"] = meta.to_dict("records")
                    if len(meta) != 1:
                        raise ValueError("Expected one differential metadata row")
                    m = meta.iloc[0]
                    if performance is not None:
                        truth = performance["truth_selection"]
                        matching = bool(np.isclose(m.q2_cut_events, m.q2_max_flux) and
                                        np.isclose(m.q2_cut_events, truth["Q2_GeV2"][1]) and
                                        np.allclose([m.w_min, m.w_max], truth["W_GeV"]))
                        _check(report, state, observable, "recorded_Q2_W_vs_performance",
                               "PASS" if matching else "FAIL")
                    widths_ok = bool(np.allclose(raw.x_high-raw.x_low, raw.bin_width))
                    _check(report, state, observable, "stored_bin_width", "PASS" if widths_ok else "FAIL")
                    central_key, count_key, factor_key = "dsigma_dx_pb", "n_reco_visible", "reco_count_to_dsigma_dx_pb"
                rows = audit_bin_counts(raw, state=state, observable=observable,
                                        target_lumi=config.lumi_invfb, coverage=config.coverage,
                                        pt2_display_range=fit_range)
                rows["energy"] = config.energy
                rows["source_file"] = str(file)
                all_bins.append(rows)
                _audit_geometry(report, state, observable, rows)
                count_ok = np.isfinite(rows.expected_reco_count) & (rows.expected_reco_count >= 0)
                _check(report, state, observable, "all_expected_reco_counts_finite_nonnegative",
                       "PASS" if count_ok.all() else "FAIL",
                       invalid_bin_indices=rows.loc[~count_ok, "bin_index"].tolist())
                active = rows.upstream_selected & rows.upstream_flags_valid
                central_ok = np.isfinite(rows.central_cross_section) & (rows.central_cross_section >= 0)
                factor_ok = np.isfinite(rows.count_to_cross_section) & (rows.count_to_cross_section > 0)
                active_ok = (~active) | (central_ok & factor_ok)
                _check(report, state, observable, "active_bins_have_valid_central_and_correction",
                       "PASS" if active_ok.all() else "FAIL",
                       invalid_active_bin_indices=rows.loc[~active_ok, "bin_index"].tolist(),
                       interpretation="Missing/invalid corrections are allowed only for explicitly masked bins; missing flags do not authorize masking.")
                for key in ("plotted", "passes_min_eff"):
                    if key in raw:
                        known = raw[key].astype(str).str.lower().isin(("true", "false", "1", "0", "1.0", "0.0"))
                        _check(report, state, observable, "strict_boolean_" + key,
                               "PASS" if known.all() else "FAIL",
                               invalid_source_rows=raw.index[~known].tolist(),
                               invalid_values=raw.loc[~known, key].astype(str).tolist())
                failed = rows.loc[rows.count_factor_closure == "FAIL", "bin_index"].tolist()
                n_tested = int(rows.supported.sum())
                _check(report, state, observable, "reco_count_factor_closure",
                       "FAIL" if failed else ("PASS" if factor_key in raw and n_tested else "WARN"), failed_bin_indices=failed,
                       n_tested_bins=n_tested,
                       independent_factor_present=bool(factor_key in raw),
                       interpretation="Independent stored count-factor closure" if factor_key in raw else "Factor inferred from central/count; this closure is tautological and not independent.",
                       unsupported_bin_indices=rows.loc[~rows.supported, "bin_index"].tolist())
                # Cross-check the independent calculation against the renderer.
                try:
                    projected = project_cross_section(raw.sort_values("W_low" if observable == "W" else "x_low").reset_index(drop=True),
                        config.lumi_invfb, central_key, count_column=count_key,
                        scale_column=factor_key, coverage=config.coverage)
                    use = rows.supported.to_numpy(bool)
                    matching = bool(np.allclose(rows.expected_reco_count, projected.mean_count, equal_nan=True) and
                                    np.allclose(rows.lower_cross_section[use], projected.lower[use], equal_nan=True) and
                                    np.allclose(rows.upper_cross_section[use], projected.upper[use], equal_nan=True))
                    _check(report, state, observable, "renderer_projection_agreement",
                           "WARN" if not use.any() else ("PASS" if matching else "FAIL"),
                           n_tested_supported_bins=int(use.sum()))
                except (ValueError, KeyError, TypeError) as error:
                    _error(report, state, observable, error)
                state_info["observables"][observable] = _summarize_bins(rows, state_info["mass_expected_selected_yield"])
            except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
                _error(report, state, observable, error)

    if report["errors"] or any(c["status"] == "FAIL" for c in report["checks"]):
        report["status"] = "INPUT_ISSUES"
    bins = pd.concat(all_bins, ignore_index=True) if all_bins else pd.DataFrame(columns=["state", "observable", "expected_reco_count"])
    report["bin_records"] = bins.to_dict("records")
    if write_outputs:
        directory = Path(out) if out is not None else Path(config.output_dir) / "binning_audit"
        directory.mkdir(parents=True, exist_ok=True)
        csv_path, json_path = directory / "projected_bins.csv", directory / "projected_binning_audit.json"
        bins.to_csv(csv_path, index=False)
        report["outputs"] = dict(bins_csv=str(csv_path), audit_json=str(json_path))
        write_json(json_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, help="Shared plot profile; no automatic run discovery")
    parser.add_argument("--out", type=Path, help="Fresh output directory; existing directories are refused")
    parser.add_argument("--pt2-range", nargs=2, type=float, default=(0., 1.), metavar=("LOW", "HIGH"))
    args = parser.parse_args()
    profile = load_profile(args.profile)
    if args.out is not None:
        out = args.out.expanduser().resolve()
        if out.exists():
            parser.error("Output directory exists; use a fresh --out")
    else:
        out = new_output_dir(profile, "binning_audit")
    cfg = make_config(profile=profile, output_dir=out)
    report = audit_projected_binning(cfg, out, pt2_display_range=args.pt2_range)
    print("Audit status:", report["status"], "(not a binning-validation verdict)")
    for state, result in report["states"].items():
        print(state, "mass expected selected yield:", result["mass_expected_selected_yield"])
        for observable, summary in result["observables"].items():
            print(" ", observable, "sum:", summary["full_expected_reco_sum"],
                  "displayed bins with mu<5:", summary["displayed_n_bins_mu_lt_5"],
                  "/", summary["displayed_n_bins"])
    print("Cross-observable selection match: NOT_VERIFIED")
    print("Outputs:", json.dumps(report.get("outputs", {})))
    if report["errors"]:
        print("Input errors:", json.dumps(report["errors"], indent=2))
    raise SystemExit(2 if report["status"] == "INPUT_ISSUES" else 0)


if __name__ == "__main__":
    main()
