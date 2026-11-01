#!/usr/bin/env python3
"""Independent ROOT normalization and EPA-flux audit; never rewrites inputs.

Uses only NumPy/uproot plus the standard library, without importing framework
extractors. Historical rates retain the archived convention; latest rates use
the reproduced July-card generator reference. Linking an EOS tuple to that
exact production configuration remains a separate provenance check.
Exit 1 means a numerical/input check failed;
unresolved physical provenance is separately recorded even when checks pass.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import (add_branching_fraction_argument,
                                branching_fractions_from_args,
                                branching_fraction_metadata, parse_branching_fractions)

ALPHA = 1.0 / 137.035999084
ME = 0.00051099895
MP = 0.9382720813
LUMI = 10.0
# The archived flux scan used rounded Lorentz factors, not exact nominal
# 10x100 and 18x275 energies. These inputs are recorded in the old
# config/pipeline.toml task validate_estarlight_flux.
REFERENCE_GAMMAS = {"low_electron": 19569., "low_proton": 106.6,
                    "high_electron": 35225., "high_proton": 293.2}
STATES = {
    "1S": {"mass": 9.46030, "br": 0.0234, "sigma": 0.053},
    "2S": {"mass": 10.02326, "br": 0.0191, "sigma": 0.056},
    "3S": {"mass": 10.35520, "br": 0.0218, "sigma": 0.059},
}
SAMPLES = {
    "old_low": {"ee": 10., "ep": 100., "mode": "state_sumw", "targets": {"1S": 104.793},
                "w_min": 20., "w_max": 26.5, "n_bins": 10,
                "rate_source": "Historical note configuration: 104.793 fb visible 1S at 10x100"},
    "old_high": {"ee": 18., "ep": 275., "mode": "state_sumw",
                 "targets": {"1S": 327.584, "2S": 84.780, "3S": 68.929},
                 "w_min": 37., "w_max": 114., "n_bins": 9,
                 "rate_source": "Historical 18x275 note normalization; report table labels these rates 10x275, appendix labels 18x275"},
    "new_low": {"ee": 9., "ep": 130., "mode": "all_generated_constant", "targets": {"1S": 119.685},
                "w_min": 20., "w_max": 70., "n_bins": 10,
                "branching_fractions": {"1S": 0.0238, "2S": 0.0191, "3S": 0.0218},
                "normalization_profile": "jul26-estarlight-939b11a",
                "rate_reference_reproduced": True,
                "rate_source": "Reproduced eSTARlight 939b11a with Jul26 Upsilon_Sim_Cards f62965a: 9x130, 1S, Q2 0-1 GeV^2, visible dielectron 119.685 fb, generator B_ee=0.0238"},
    "new_high": {"ee": 9., "ep": 275., "mode": "all_generated_constant", "targets": {"1S": 210.136},
                 "w_min": 20., "w_max": 100., "n_bins": 16,
                 "branching_fractions": {"1S": 0.0238, "2S": 0.0191, "3S": 0.0218},
                 "normalization_profile": "jul26-estarlight-939b11a",
                 "rate_reference_reproduced": True,
                 "rate_source": "Reproduced eSTARlight 939b11a with Jul26 Upsilon_Sim_Cards f62965a: 9x275, 1S, Q2 0-1 GeV^2, visible dielectron 210.136 fb, generator B_ee=0.0238"},
}
REQUIRED = ("weight", "true_mass", "true_Q2", "true_W", "reco_mass", "reco_W", "sigma_epz")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def state_mask(mass, state):
    # Normalization precedes all Q2, W, and reconstruction selections.
    m = STATES[state]["mass"]
    return np.isfinite(mass) & (mass >= m - 0.2) & (mass <= m + 0.2)


def normalized_weights(data, role, state, generated, branching_fractions=None):
    """Return weights and normalization ledger, with no kinematic preselection."""
    spec = SAMPLES[role]
    br = parse_branching_fractions(spec.get("branching_fractions") if branching_fractions is None else branching_fractions)[state]
    raw = np.asarray(data["weight"], dtype=float)
    if raw.size == 0 or not np.all(np.isfinite(raw) & (raw > 0)):
        raise ValueError("Positive, finite raw weights and a nonempty tuple are required")
    target = spec["targets"][state]
    mask = state_mask(data["true_mass"], state)
    if spec["mode"] == "state_sumw":
        raw_scope = float(np.sum(raw[mask], dtype=float))
        if raw_scope <= 0:
            raise ValueError(f"No positive raw normalization weight for {role} {state}")
        scale = target * LUMI / raw_scope
        weights = raw * scale
        normalized_scope = float(np.sum(weights[mask], dtype=float))
        scope = "truth state within +/-0.2 GeV, before Q2/W/reconstruction cuts"
        scope_entries = int(mask.sum())
        constant = None
    else:
        if not np.isfinite(generated) or not np.isclose(generated, raw.size, rtol=0, atol=1e-6):
            raise ValueError("h_EventCount must equal tuple entries for all-generated constant weighting")
        raw_scope = float(np.sum(raw, dtype=float))
        constant = target * LUMI / generated
        weights = np.full(raw.shape, constant, dtype=float)
        normalized_scope = float(np.sum(weights, dtype=float))
        scale = None  # Replacement is not a common scale if raw job weights vary.
        scope = "all tuple entries, assumed all unweighted generated events, before every selection"
        scope_entries = int(raw.size)
    ledger = dict(sample=role, state=state, electron_energy_GeV=spec["ee"], proton_energy_GeV=spec["ep"],
                  normalization_mode=spec["mode"], normalization_scope=scope, normalization_scope_entries=scope_entries,
                  visible_cross_section_assigned_fb=target, luminosity_invfb=LUMI,
                  raw_sum_weight_scope=raw_scope, raw_inferred_visible_cross_section_fb=raw_scope / LUMI,
                  weight_scale=scale, replacement_constant_weight=constant,
                  normalized_sum_weight_scope=normalized_scope, target_sum_weight_scope=target * LUMI,
                  normalized_sum_weight_truth_state=float(np.sum(weights[mask], dtype=float)),
                  truth_state_entries=int(mask.sum()), truth_state_fraction_of_all_entries=float(mask.mean()),
                  branching_fraction=br, branching_divisions_in_ep_to_gamma_p=1,
                  branching_fraction_convention="sample_default" if branching_fractions is None else "explicit_cli",
                  branching_fraction_uncertainty_included=False,
                  normalization_profile=spec.get("normalization_profile", "historical_note"),
                  rate_reference_reproduced=spec.get("rate_reference_reproduced", False),
                  assigned_inclusive_ep_cross_section_pb=target / (1000 * br),
                  rate_source=spec["rate_source"], physical_rate_provenance_verified=False)
    return weights, ledger


def photon_flux_dW(w, ee, ep, q2_max):
    """Transverse EPA flux, integrated over kinematic Q2_min through Q2_max.

    Uses y=(W^2-mp^2)/(4EeEp), the same low-Q2 approximation as the
    framework, independently implemented here. Values outside support vanish.
    """
    w = np.asarray(w, dtype=float)
    y = (w * w - MP * MP) / (4 * ee * ep)
    valid = np.isfinite(w) & (w > 0) & (y > 0) & (y < 1)
    yy = np.where(valid, y, 0.5)
    qmin = ME * ME * yy * yy / (1 - yy)
    valid &= qmin < q2_max
    dndy = ALPHA / (2 * math.pi) * (
        (1 + (1 - yy) ** 2) / yy * np.log(q2_max / qmin)
        - 2 * (1 - yy) / yy * (1 - qmin / q2_max))
    return np.where(valid, dndy * 2 * w / (4 * ee * ep), 0.)


def integrated_flux(low, high, ee, ep, q2_max):
    if not all(np.isfinite([low, high, ee, ep, q2_max])) or low >= high or min(ee, ep, q2_max) <= 0:
        raise ValueError("Finite increasing W limits, positive beam energies and positive Q2_max are required")
    # Split at the exact support endpoint of this approximate EPA integrand.
    ymax = 2 / (1 + math.sqrt(1 + 4 * ME * ME / q2_max))
    low = max(low, MP)
    high = min(high, math.sqrt(MP * MP + 4 * ee * ep * ymax))
    if high <= low:
        return 0.
    x, weights = np.polynomial.legendre.leggauss(96)
    w = (high + low) / 2 + (high - low) / 2 * x
    return float((high - low) / 2 * np.dot(weights, photon_flux_dW(w, ee, ep, q2_max)))


def check(checks, name, passed, **details):
    checks.append(dict(check=name, passed=bool(passed), **details))


def compare_flux_reference(path, checks):
    """Re-evaluate analytic columns; reuse, never regenerate, generator columns."""
    result = []
    with Path(path).open(newline="") as stream:
        refs = list(csv.DictReader(stream))
    if not refs:
        raise ValueError("Empty historical flux reference")
    for i, row in enumerate(refs):
        lo, hi = float(row["W_low"]), float(row["W_high"])
        low_ee, low_ep = REFERENCE_GAMMAS["low_electron"] * ME, REFERENCE_GAMMAS["low_proton"] * MP
        high_ee, high_ep = REFERENCE_GAMMAS["high_electron"] * ME, REFERENCE_GAMMAS["high_proton"] * MP
        low = integrated_flux(lo, hi, low_ee, low_ep, 0.001)
        high = integrated_flux(lo, hi, high_ee, high_ep, 0.001)
        low_old, high_old = float(row["analytic_flux_low"]), float(row["analytic_flux_high"])
        ratio = low / high
        internal = float(row["low_total_pb"]) / float(row["high_total_pb"])
        double = internal / ratio
        result.append(dict(W_low=lo, W_high=hi, Q2_max_GeV2=0.001,
                           actual_low_electron_energy_GeV=low_ee, actual_low_proton_energy_GeV=low_ep,
                           actual_high_electron_energy_GeV=high_ee, actual_high_proton_energy_GeV=high_ep,
                           recomputed_flux_10x100=low, recorded_flux_10x100=low_old,
                           relative_difference_flux_10x100=low / low_old - 1,
                           recomputed_flux_18x275=high, recorded_flux_18x275=high_old,
                           relative_difference_flux_18x275=high / high_old - 1,
                           recomputed_analytic_low_over_high=ratio,
                           recorded_generator_low_total_pb=float(row["low_total_pb"]),
                           recorded_generator_high_total_pb=float(row["high_total_pb"]),
                           recorded_generator_low_generated_pb=float(row["low_generated_pb"]),
                           recorded_generator_high_generated_pb=float(row["high_generated_pb"]),
                           recorded_generator_ratio_recomputed=internal,
                           generator_over_analytic_double_ratio=double,
                           recorded_double_ratio=float(row["internal_over_analytic_ratio"]),
                           generator_was_rerun=False))
        for name, actual, expected in [
            ("low_flux", low, low_old), ("high_flux", high, high_old),
            ("analytic_ratio", ratio, float(row["analytic_flux_ratio_low_over_high"])),
            ("recorded_generator_ratio", internal, float(row["internal_ratio_low_over_high"])),
            ("double_ratio", double, float(row["internal_over_analytic_ratio"]))]:
            check(checks, f"flux_reference_bin{i}_{name}", np.isclose(actual, expected, rtol=2e-7, atol=1e-14),
                  actual=actual, expected=expected, relative_tolerance=2e-7)
    return result


def range_summary(values, mask):
    selected = np.asarray(values)[mask & np.isfinite(values)]
    return dict(min=float(selected.min()) if len(selected) else None,
                max=float(selected.max()) if len(selected) else None)


def audit_sample(path, role, checks, branching_fractions=None):
    spec = SAMPLES[role]
    with uproot.open(path) as root:
        tree = root["UpsilonTree"]
        missing = sorted(set(REQUIRED) - set(tree.keys()))
        if missing:
            raise ValueError(f"Missing branches: {missing}")
        data = tree.arrays(list(REQUIRED), library="np")
        generated = float(np.sum(root["h_EventCount"].values(flow=False))) if "h_EventCount" in root else float("nan")
        entries = int(tree.num_entries)
    raw = np.asarray(data["weight"], float)
    finite_rows = np.logical_and.reduce([np.isfinite(data[k]) for k in REQUIRED])
    physical_truth = np.isfinite(data["true_mass"]) & (data["true_mass"] > 0) & np.isfinite(data["true_W"]) & (data["true_W"] > 0) & np.isfinite(data["true_Q2"]) & (data["true_Q2"] >= 0)
    masses = {state: int(state_mask(data["true_mass"], state).sum()) for state in STATES}
    weight_bad = int(np.sum(~np.isfinite(raw) | (raw <= 0)))
    check(checks, role + "_positive_finite_weights", weight_bad == 0, bad_rows=weight_bad)
    check(checks, role + "_event_count_matches_rows", np.isfinite(generated) and np.isclose(generated, entries, rtol=0, atol=1e-6),
          entries=entries, histogram_generated_events=generated)
    unique, counts = np.unique(raw[np.isfinite(raw)], return_counts=True)
    summary = dict(sample=role, input=str(Path(path).resolve()), sha256=digest(path), entries=entries,
                   h_EventCount=generated, raw_sum_weight=float(np.sum(raw)),
                   finite_all_required_branches_rows=int(finite_rows.sum()), physical_truth_rows=int(physical_truth.sum()),
                   nonfinite_by_branch={k: int(np.sum(~np.isfinite(data[k]))) for k in REQUIRED},
                   invalid_physical_truth_by_branch={
                       "true_mass_nonpositive_or_nonfinite": int(np.sum(~np.isfinite(data["true_mass"]) | (data["true_mass"] <= 0))),
                       "true_W_nonpositive_or_nonfinite": int(np.sum(~np.isfinite(data["true_W"]) | (data["true_W"] <= 0))),
                       "true_Q2_negative_or_nonfinite": int(np.sum(~np.isfinite(data["true_Q2"]) | (data["true_Q2"] < 0)))},
                   nonpositive_or_nonfinite_weight_rows=weight_bad, raw_weight_unique_count=len(unique),
                   raw_weight_values=unique.tolist() if len(unique) <= 100 else unique[:100].tolist(),
                   raw_weight_counts=counts.tolist() if len(unique) <= 100 else counts[:100].tolist(),
                   raw_weight_values_truncated=bool(len(unique) > 100),
                   mass_state_entries=masses, entries_outside_all_mass_windows=entries - sum(masses.values()),
                   sample_assumptions=spec,
                   caveat="Histogram/tuple agreement is internal consistency, not independent proof that no generated events were filtered upstream.")
    ledgers, coverage, support, bin_rows = [], [], [], []
    for state in spec["targets"]:
        weights, ledger = normalized_weights(data, role, state, generated, branching_fractions)
        ledgers.append(ledger)
        check(checks, f"{role}_{state}_sum_weights_equals_assigned_sigma_times_L",
              np.isclose(ledger["normalized_sum_weight_scope"], ledger["target_sum_weight_scope"], rtol=2e-12, atol=1e-10),
              actual=ledger["normalized_sum_weight_scope"], expected=ledger["target_sum_weight_scope"])
        sm = state_mask(data["true_mass"], state)
        valid_truth = sm & np.isfinite(data["true_W"]) & (data["true_W"] > 0) & np.isfinite(data["true_Q2"]) & (data["true_Q2"] >= 0)
        invalid_w = sm & (~np.isfinite(data["true_W"]) | (data["true_W"] <= 0))
        ledger.update(truth_state_invalid_W_entries=int(invalid_w.sum()),
                      truth_state_invalid_W_raw_weight_fraction=float(raw[invalid_w].sum() / raw[sm].sum()),
                      truth_state_invalid_W_normalized_weight=float(weights[invalid_w].sum()),
                      truth_state_valid_physical_coordinates_entries=int(valid_truth.sum()))
        q = data["true_Q2"]
        state_sum = float(weights[sm].sum())
        for qmax in (0.001, 1.0):
            qm = sm & np.isfinite(q) & (q >= 0) & (q <= qmax)
            coverage.append(dict(sample=role, state=state, Q2_max_GeV2=qmax,
                                 state_entries=int(sm.sum()), selected_Q2_entries=int(qm.sum()),
                                 entry_fraction_of_state=float(qm.sum() / sm.sum()),
                                 selected_Q2_weight=float(weights[qm].sum()),
                                 weight_fraction_of_normalized_state=float(weights[qm].sum() / state_sum),
                                 raw_Q2_selected_state_sum=float(raw[qm].sum()),
                                 hypothetical_scale_if_normalized_after_Q2=(spec["targets"][state] * LUMI / float(raw[qm].sum())) if float(raw[qm].sum()) > 0 else None,
                                 post_Q2_vs_precut_state_scale_ratio=(float(raw[sm].sum()) / float(raw[qm].sum())) if float(raw[qm].sum()) > 0 else None,
                                 state_Q2_negative_or_nonfinite_entries=int(np.sum(sm & (~np.isfinite(q) | (q < 0)))),
                                 state_Q2_above_cut_entries=int(np.sum(sm & np.isfinite(q) & (q > qmax)))))
            wm = valid_truth & (q <= qmax)
            r = range_summary(data["true_W"], wm)
            support.append(dict(sample=role, state=state, Q2_max_GeV2=qmax, truth_entries=int(wm.sum()),
                                W_min_GeV=r["min"], W_max_GeV=r["max"],
                                W_kinematic_endpoint_GeV=math.sqrt(MP * MP + 4 * spec["ee"] * spec["ep"]),
                                entries_below_Upsilon_proton_threshold=int(np.sum(wm & (data["true_W"] < STATES[state]["mass"] + MP))),
                                entries_above_approximate_beam_endpoint=int(np.sum(wm & (data["true_W"] > math.sqrt(MP * MP + 4 * spec["ee"] * spec["ep"]))))))
            for label, edges in [("note_configuration", np.linspace(spec["w_min"], spec["w_max"], spec["n_bins"] + 1)),
                                 ("common_overlap", np.array([20., 22.1667, 24.3333, 26.5]))]:
                # Reco selection mirrors the W extractor; it does not require a
                # truth-state tag. The normalization scope is nevertheless the
                # truth-state window above, and remains fixed for both Q2 cuts.
                p = dict(STATES[state], br=ledger["branching_fraction"])
                rm = (np.isfinite(data["reco_mass"]) & (np.abs(data["reco_mass"] - p["mass"]) <= 3 * p["sigma"])
                      & np.isfinite(data["sigma_epz"]) & (np.abs(data["sigma_epz"] - 2 * spec["ee"]) <= 0.2 * spec["ee"])
                      & np.isfinite(q) & (q >= 0) & (q <= qmax) & np.isfinite(data["reco_W"]))
                truth = np.histogram(data["true_W"][wm], bins=edges, weights=weights[wm])[0]
                reco = np.histogram(data["reco_W"][rm], bins=edges, weights=weights[rm])[0]
                for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                    phi = integrated_flux(lo, hi, spec["ee"], spec["ep"], qmax)
                    t, r = float(truth[i]), float(reco[i])
                    eff = r / t if t > 0 else None
                    supported = t > 0 and r > 0 and phi > 0
                    corrected = r / eff if supported else None
                    br = p["br"]
                    truth_sigma = t / (LUMI * 1000 * phi * br) if phi > 0 else None
                    sigma = corrected / (LUMI * 1000 * phi * br) if supported else None
                    row = dict(sample=role, state=state, binning=label, Q2_max_GeV2=qmax,
                               branching_fraction=br,
                               W_low=lo, W_high=hi, integrated_flux=phi, truth_weight=t, reco_weight=r,
                               reco_over_truth=eff, supported=bool(supported), corrected_visible_yield=corrected,
                               sigma_gamma_p_pb=sigma, sigma_gamma_p_truth_diagnostic_pb=truth_sigma,
                               inverse_conversion_visible_yield=sigma * LUMI * 1000 * phi * br if supported else None,
                               sigma_gamma_p_times_Bee_pb=sigma * br if supported else None,
                               raw_truth_state_fraction_in_bin=t / state_sum,
                               note="Unsupported corrected values are null, not zero; truth-only diagnostic is not a measurement.")
                    bin_rows.append(row)
                    if supported:
                        check(checks, f"{role}_{state}_{label}_q{qmax}_bin{i}_normalization_identity",
                              np.isclose(sigma, truth_sigma, rtol=2e-12, atol=1e-12)
                              and np.isclose(row["inverse_conversion_visible_yield"], t, rtol=2e-12, atol=1e-10),
                              corrected_gamma_p_pb=sigma, truth_gamma_p_pb=truth_sigma)
    return summary, ledgers, coverage, support, bin_rows


def plain(value):
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows([plain(row) for row in rows])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    br_argument = add_branching_fraction_argument(parser)
    br_argument.help = ("Override dielectron fractions for all requested samples. Without this option, "
                        "historical samples use 1S:0.0234 and July samples use 1S:0.0238; "
                        "2S:0.0191 and 3S:0.0218 in both conventions. No branching uncertainty is added.")
    parser.add_argument("--old-low", required=True, type=Path)
    parser.add_argument("--old-high", required=True, type=Path)
    parser.add_argument("--new-low", type=Path)
    parser.add_argument("--new-high", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--flux-reference", type=Path,
                        default=Path(__file__).resolve().parents[2] / "reference_data/flux_ratio_scan_W20to26p5_Q2lt0p001.csv")
    args = parser.parse_args(argv)
    branching_fractions = None if args.branching_fractions is None else branching_fractions_from_args(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    checks, samples, ledger, coverage, support, bins = [], [], [], [], [], []
    for role in SAMPLES:
        path = getattr(args, role)
        if path is None:
            continue
        try:
            one = audit_sample(path, role, checks, branching_fractions)
            samples.append(one[0])
            for target, rows in zip((ledger, coverage, support, bins), one[1:]):
                target.extend(rows)
        except Exception as error:
            check(checks, role + "_input_audit", False, input=str(path), error=f"{type(error).__name__}: {error}")
    flux = []
    try:
        flux = compare_flux_reference(args.flux_reference, checks)
    except Exception as error:
        check(checks, "historical_flux_reference", False, error=f"{type(error).__name__}: {error}")
    br_metadata = branching_fraction_metadata(args) if branching_fractions is not None else {
        "branching_fraction_convention": "per_sample_defaults",
        "branching_fractions_by_sample": {role: parse_branching_fractions(spec.get("branching_fractions")) for role, spec in SAMPLES.items()},
        "branching_fraction_uncertainty_included": False,
    }
    provenance = dict(**br_metadata,
                      created_utc=datetime.now(timezone.utc).isoformat(), validator_sha256=digest(__file__),
                      generator_was_rerun=False, historical_flux_reference=str(args.flux_reference.resolve()),
                      historical_flux_reference_sha256=digest(args.flux_reference) if args.flux_reference.exists() else None,
                      generator_reference_status="Stored historical generator cross sections only; independently recomputed analytic EPA columns and ratios.",
                      historical_flux_scan_beam_gammas=REFERENCE_GAMMAS,
                      historical_flux_beam_source="Old framework config/pipeline.toml, validate_estarlight_flux task; gamma-to-energy conversion uses this validator's electron/proton masses.",
                      assigned_rates=SAMPLES, input_files=[{k: s[k] for k in ("sample", "input", "sha256", "entries")} for s in samples],
                      units="Input visible cross sections: fb; L: inverse fb; output gamma-p cross sections: pb; W: GeV; Q2: GeV^2",
                      limitations=["Numerical closure does not independently establish generator cross-section provenance.",
                                   "The new constant weighting assumes an unweighted generated sample and no missing generated events.",
                                   "h_EventCount versus row count is an internal consistency check, not proof of complete generation coverage.",
                                   "Finite sentinels such as -999 are recorded separately from valid physical truth values.",
                                   "Recorded true_Q2 coverage concerns this tuple branch; it does not certify the generator's event definition.",
                                   "This audit normalizes before Q2 cuts. Historical extractors can include Q2 inside their truth_state normalization scope; q2_coverage.csv exposes that alternative scale for comparison.",
                                   "Same-sample R/(R/T)=T is algebraic closure, not independent detector-model validation.",
                                   "EPA neglects Q2 in the W-to-y relation and does not validate the generator's full finite-Q2 integration."])
    summary = dict(numerical_checks_passed=all(c["passed"] for c in checks), checks_total=len(checks),
                   failed_checks=[c for c in checks if not c["passed"]], samples=samples,
                   interpretation_flags=[dict(sample=r["sample"], state=r["state"],
                                              issue="Truth-state normalization includes rows with missing/nonphysical W; inspect input reconstruction before interpreting differential rates.",
                                              entries=r["truth_state_invalid_W_entries"],
                                              raw_state_weight_fraction=r["truth_state_invalid_W_raw_weight_fraction"])
                                         for r in ledger if r["truth_state_invalid_W_entries"] > 0],
                   physical_normalization_provenance_verified=False, generator_was_rerun=False)
    for name, rows in [("normalization_ledger", ledger), ("q2_coverage", coverage), ("w_support", support),
                       ("w_bin_checks", bins), ("flux_reference_check", flux), ("numeric_checks", checks)]:
        write_csv(args.out_dir / (name + ".csv"), rows)
    for name, obj in [("audit_summary", summary), ("provenance", provenance)]:
        (args.out_dir / (name + ".json")).write_text(json.dumps(plain(obj), indent=2, allow_nan=False) + "\n")
    print(f"{'PASS' if summary['numerical_checks_passed'] else 'FAIL'}: {len(checks)} numerical checks; physical rate provenance remains unverified")
    print(f"Audit outputs: {args.out_dir.resolve()}")
    return 0 if summary["numerical_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
