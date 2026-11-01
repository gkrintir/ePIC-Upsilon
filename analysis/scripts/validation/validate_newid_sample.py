#!/usr/bin/env python3
"""Validate an all-generated UpsilonTree and its electron-ID/fEMC configuration.

This reads flattened analysis tuples, not PODIO events. It never selects rows or
changes weights. Formula checks establish stored beam constants, not generator
provenance or the correctness of the externally supplied total cross section.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot

COUNTERS = ["ele_id_n" + name for name in ("Matched", "Missing", "Masked", "Ambiguous")]
BEAM_FIELDS = ["true_W", "true_Q2", "true_scat_eMinus_pT", "true_scat_eMinus_eta", "true_scat_eMinus_E"]
ID_FIELDS = ["ele_id_mode"] + COUNTERS
FEMC_MODES = {"full": 0, "keep_xpos": 3, "keep_yneg": 4}
REQUIRED = ID_FIELDS + ["femc_mode"] + BEAM_FIELDS


def validate(path, electron_energy, proton_energy, tree_name="UpsilonTree", hash_input=False,
             *, electron_id="reco", femc_scenario="full", allow_missing_legacy_id_metadata=False,
             allow_missing_full_femc_metadata=False):
    """Audit raw tuple metadata and beam formula without selecting or changing rows.

    Missing metadata allowances are explicit legacy-only exceptions. They record
    unverified fields and never override contradictory values present in a file.
    Half-plane fEMC modes always require their stored mode and legacy electron ID.
    """
    report = {"input": str(Path(path).expanduser().resolve()), "tree": tree_name,
              "status": "FAIL", "errors": [], "warnings": [],
              "requested_electron_id": electron_id, "requested_femc_scenario": femc_scenario,
              "metadata_validation": {},
              "normalization_scope": "All generated rows, before ID, reconstruction, mass, Q2 or W cuts",
              "beam_check_scope": "Stored truth W/Q2 formula only; not generator beam, rate or phase-space provenance"}
    errors = report["errors"]
    try:
        if electron_id not in ("legacy", "reco") or femc_scenario not in FEMC_MODES:
            raise ValueError("Unsupported electron ID or fEMC scenario")
        if electron_id == "reco" and femc_scenario != "full":
            raise ValueError("Half-plane fEMC scenarios require legacy electron ID in the supplied producer")
        if allow_missing_legacy_id_metadata and electron_id != "legacy":
            raise ValueError("Missing legacy-ID metadata allowance requires --electron-id legacy")
        if allow_missing_full_femc_metadata and (femc_scenario != "full" or electron_id != "legacy"):
            raise ValueError("Missing full-fEMC metadata allowance is only for a full legacy-ID baseline")
        if not all(np.isfinite(x) and x > 0 for x in (electron_energy, proton_energy)):
            raise ValueError("Beam energies must be finite and positive")
        with uproot.open(report["input"]) as root:
            if tree_name not in root:
                raise ValueError("Missing flattened analysis tree " + tree_name)
            tree = root[tree_name]
            if tree.classname != "TTree":
                raise ValueError("Expected a TTree compatible with prepare_latest_sample.py")
            rows = int(tree.num_entries)
            report["entries"] = rows
            if rows == 0:
                raise ValueError("Empty tuple")
            if "framework_normalization" in root:
                errors.append("Already normalized input; use the original merged file")
            if "h_EventCount" not in root:
                errors.append("Missing h_EventCount; generated-event normalization cannot be established")
            else:
                generated = float(np.sum(root["h_EventCount"].values(flow=False)))
                report["generated_events_h_EventCount"] = generated if np.isfinite(generated) else None
                if not np.isfinite(generated) or not np.isclose(generated, rows, rtol=0, atol=1e-6):
                    errors.append("h_EventCount differs from tuple rows; refusing selected-event normalization")
            available = set(tree.keys())
            required = set(BEAM_FIELDS)
            if not allow_missing_legacy_id_metadata:
                required.update(ID_FIELDS)
            if not allow_missing_full_femc_metadata:
                required.add("femc_mode")
            missing = sorted(required - available)
            if missing:
                raise ValueError("Missing required branches: " + ", ".join(missing))
            names = sorted((set(REQUIRED) | {"reco_nTracks_Central", "femc_nClustersRemoved"}) & available)
            arrays = {name: tree[name].array(library="np") for name in names}
        metadata = {}
        for name in ID_FIELDS + ["femc_mode", "femc_nClustersRemoved"]:
            if name not in arrays:
                continue
            values = arrays[name]
            if values.ndim != 1 or not np.all(np.isfinite(values) & (values == np.floor(values))):
                raise ValueError("Noninteger or nonscalar metadata in " + name)
            unique, counts = np.unique(values, return_counts=True)
            metadata[name] = {str(int(value)): int(count) for value, count in zip(unique, counts)}
        report["metadata_counts"] = metadata
        expected_id = 1 if electron_id == "reco" else 0
        missing_id = [name for name in ID_FIELDS if name not in arrays]
        id_verified = "ele_id_mode" in arrays and bool(np.all(arrays["ele_id_mode"] == expected_id))
        report["metadata_validation"]["electron_id"] = {
            "expected_mode": expected_id, "mode_verified": id_verified,
            "missing_branches": missing_id, "missing_metadata_allowed": allow_missing_legacy_id_metadata}
        if missing_id:
            report["warnings"].append("Legacy electron-ID metadata incomplete and unverified for missing branches: " +
                                      ", ".join(missing_id) + "; explicit allowance used, no metadata fabricated")
        if "ele_id_mode" in arrays and not id_verified:
            errors.append("ele_id_mode must be " + str(expected_id) + " on every row for " + electron_id + " ID")
        expected_femc = FEMC_MODES[femc_scenario]
        femc_verified = "femc_mode" in arrays and bool(np.all(arrays["femc_mode"] == expected_femc))
        report["metadata_validation"]["femc"] = {
            "expected_mode": expected_femc, "mode_verified": femc_verified,
            "missing_metadata_allowed": allow_missing_full_femc_metadata}
        if "femc_mode" not in arrays:
            report["warnings"].append("Full-fEMC mode unverified because femc_mode is absent; explicit legacy-baseline allowance used")
        elif not femc_verified:
            errors.append("femc_mode must be " + str(expected_femc) + " for scenario " + femc_scenario)
        if "femc_nClustersRemoved" in arrays:
            if np.any(arrays["femc_nClustersRemoved"] < 0):
                errors.append("femc_nClustersRemoved must be nonnegative")
            if femc_scenario == "full" and np.any(arrays["femc_nClustersRemoved"] != 0):
                errors.append("Full-fEMC sample contains removed clusters")
        if electron_id == "legacy":
            for name in COUNTERS:
                if name in arrays and np.any(arrays[name] != -1):
                    errors.append("Legacy ID counter " + name + " must remain -1 in the supplied producer")
        else:
            if any(np.any(arrays[name] < 0) for name in COUNTERS):
                errors.append("Reco ID counters must be nonnegative on every row; producer evaluates ID before every Fill")
            if np.any(arrays["ele_id_nMasked"] != 0):
                errors.append("Full-fEMC sample contains masked electron-ID candidates")
            if "reco_nTracks_Central" in arrays:
                counter_sum = sum(arrays[name] for name in COUNTERS)
                if not np.array_equal(counter_sum, arrays["reco_nTracks_Central"]):
                    errors.append("ID counter sum differs from reco_nTracks_Central")
            else:
                report["warnings"].append("No reco_nTracks_Central branch; counter partition check unavailable")
            if np.sum(arrays["ele_id_nMatched"]) == 0:
                report["warnings"].append("No matched ID candidates recorded; inspect reconstruction before interpreting efficiency")
        # Stable on-shell pz avoids loss of precision in ROOT eta at tiny pT.
        # The producer uses SetXYZM with m_e=0.000511 and a head-on massless beam electron.
        pt, eta, energy = (arrays[name] for name in BEAM_FIELDS[2:])
        usable = np.logical_and.reduce([np.isfinite(arrays[name]) & (arrays[name] != -999) for name in BEAM_FIELDS])
        usable &= (pt >= 0) & (energy > 0) & (energy * energy - pt * pt >= 0.000511 ** 2)
        report["beam_formula"] = {"electron_energy_GeV": electron_energy, "proton_momentum_GeV": proton_energy,
                                  "checked_rows": int(usable.sum()), "unavailable_rows": int((~usable).sum()),
                                  "absolute_tolerance": 0.001, "relative_tolerance": 0.001,
                                  "pz_method": "sign(eta)*sqrt(E^2-pT^2-m_e^2), m_e=0.000511 GeV"}
        if not usable.any():
            errors.append("No usable truth scattered electrons to validate stored beam constants")
        else:
            pt, eta, energy = pt[usable], eta[usable], energy[usable]
            pz = np.sign(eta) * np.sqrt(energy * energy - pt * pt - 0.000511 ** 2)
            q0, qz = electron_energy - energy, -electron_energy - pz
            q2 = pt * pt + qz * qz - q0 * q0
            w2 = (np.sqrt(proton_energy ** 2 + 0.938272 ** 2) + q0) ** 2 - pt * pt - (proton_energy + qz) ** 2
            for name, calculated in (("true_W", np.sign(w2) * np.sqrt(np.abs(w2))), ("true_Q2", q2)):
                stored = arrays[name][usable]
                mismatches = int(np.sum(~np.isclose(calculated, stored, atol=0.001, rtol=0.001)))
                report["beam_formula"][name] = {"mismatched_rows": mismatches,
                    "max_absolute_difference": float(np.max(np.abs(calculated - stored)))}
                if mismatches:
                    errors.append(name + " disagrees with requested beam constants on " + str(mismatches) + " rows")
        if hash_input:
            sha = hashlib.sha256()
            with open(report["input"], "rb") as stream:
                for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                    sha.update(chunk)
            report["sha256"] = sha.hexdigest()
    except Exception as exc:
        errors.append(type(exc).__name__ + ": " + str(exc))
    report["status"] = "PASS" if not errors else "FAIL"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--electron-energy", type=float, required=True)
    parser.add_argument("--proton-energy", type=float, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--hash-input", action="store_true", help="Optional extra input hash; workflow already hashes inputs")
    parser.add_argument("--electron-id", choices=("legacy", "reco"), default="reco")
    parser.add_argument("--femc-scenario", choices=tuple(FEMC_MODES), default="full")
    parser.add_argument("--allow-missing-legacy-id-metadata", action="store_true",
                        help="Explicitly accept missing legacy ID fields, recorded as unverified in the audit")
    parser.add_argument("--allow-missing-full-femc-metadata", action="store_true",
                        help="Explicitly accept absent femc_mode only for a full legacy-ID baseline")
    args = parser.parse_args()
    report = validate(args.input, args.electron_energy, args.proton_energy, args.tree, args.hash_input,
                      electron_id=args.electron_id, femc_scenario=args.femc_scenario,
                      allow_missing_legacy_id_metadata=args.allow_missing_legacy_id_metadata,
                      allow_missing_full_femc_metadata=args.allow_missing_full_femc_metadata)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(report["status"] + " " + args.electron_id + " ID / " + args.femc_scenario + " input: " + report["input"])
    for warning in report["warnings"]:
        print("WARNING: " + warning)
    for error in report["errors"]:
        print("ERROR: " + error)
    raise SystemExit(0 if report["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
