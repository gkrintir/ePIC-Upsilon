"""Locate a reconstructed pilot in a merged UpsilonTree and compare every branch."""
import argparse
import contextlib
import json
from pathlib import Path

import numpy as np
import uproot


# Beam-independent generated kinematics. The search examines every possible
# contiguous start in the merged file, without assuming filenames or basket order.
TRUTH = ("true_pT", "true_mass", "true_y", "true_scat_eMinus_E",
         "true_scat_eMinus_pT", "true_scat_eMinus_eta")
ID_FIELDS = ("ele_id_nMatched", "ele_id_nMissing", "ele_id_nMasked", "ele_id_nAmbiguous")
ATOL = 1e-12
RTOL = 1e-12


def is_close(a, b):
    return np.isclose(a, b, atol=ATOL, rtol=RTOL, equal_nan=True)


def locate(pilot, reference):
    kinematics = {k: pilot[k].array(library="np") for k in TRUTH}
    valid = np.flatnonzero(np.isfinite(kinematics["true_pT"]) & (kinematics["true_pT"] > 0))
    if not len(valid):
        raise ValueError("Pilot has no usable generated transverse momentum")
    row = int(valid[0])
    all_pt = reference["true_pT"].array(library="np")
    starts = np.flatnonzero(is_close(kinematics["true_pT"][row], all_pt)) - row
    starts = [int(i) for i in starts if 0 <= i <= reference.num_entries - pilot.num_entries]
    for key in TRUTH:
        starts = [i for i in starts if np.all(is_close(
            kinematics[key], reference[key].array(entry_start=i, entry_stop=i + pilot.num_entries, library="np")))]
    if len(starts) != 1:
        raise ValueError(f"Found {len(starts)} contiguous matches in six generated observables")
    return starts[0]


def field_differences(a, b):
    if a.ndim != 1 or b.ndim != 1 or a.dtype.kind not in "biuf" or b.dtype.kind not in "biuf":
        raise ValueError("Expected a flat numerical UpsilonTree branch")
    if a.dtype != b.dtype or a.shape != b.shape:
        return {"schema_mismatch": True, "pilot_dtype": str(a.dtype), "reference_dtype": str(b.dtype)}
    sentinel_changes = (a == -999) != (b == -999)
    unequal = a != b
    if np.issubdtype(a.dtype, np.inexact):
        unequal &= ~(np.isnan(a) & np.isnan(b))
        significant = ~is_close(a, b) | sentinel_changes
        finite = np.isfinite(a) & np.isfinite(b)
        abs_diff = float(np.max(np.abs(a[finite] - b[finite]))) if np.any(finite) else None
    else:
        significant = unequal
        abs_diff = int(np.max(np.abs(a.astype(np.int64) - b.astype(np.int64))))
    a_bytes = np.ascontiguousarray(a).view("u1").reshape(len(a), -1)
    b_bytes = np.ascontiguousarray(b).view("u1").reshape(len(b), -1)
    byte_diffs = a_bytes != b_bytes
    return {
        "bit_differences": int(np.count_nonzero(np.any(byte_diffs, axis=1))),
        "value_differences": int(np.count_nonzero(unequal)),
        "outside_tolerance": int(np.count_nonzero(significant)),
        "sentinel_mask_differences": int(np.count_nonzero(sentinel_changes)),
        "nan_mask_differences": int(np.count_nonzero(np.isnan(a) != np.isnan(b))),
        "max_abs_finite_difference": abs_diff,
    }


def check(pilot_path, reference_path):
    with contextlib.ExitStack() as stack:
        files = []
        for path in (pilot_path, reference_path):
            stream = stack.enter_context(path.open("rb"))
            files.append(stack.enter_context(uproot.open(stream)))
        required_branches = set(TRUTH) | set(ID_FIELDS) | {"reco_topology", "reco_tagged_proton_p"}
        for label, source in zip(("Pilot", "Reference"), files):
            missing_objects = {"UpsilonTree", "h_EventCount"} - set(source.keys(cycle=False))
            if missing_objects:
                raise ValueError(f"{label} file is missing: {', '.join(sorted(missing_objects))}")
            event_counts = source["h_EventCount"].values(flow=True)
            if (np.any(~np.isfinite(event_counts)) or np.any(event_counts < 0)
                    or event_counts[0] != 0 or event_counts[-1] != 0
                    or event_counts.sum() != source["UpsilonTree"].num_entries):
                raise ValueError(f"{label} h_EventCount disagrees with the generated row count")
            missing = required_branches - set(source["UpsilonTree"].keys())
            if missing:
                raise ValueError(f"{label} UpsilonTree is missing required branches: {', '.join(sorted(missing))}")
        pilot, reference = [f["UpsilonTree"] for f in files]
        start = locate(pilot, reference)
        stop = start + pilot.num_entries
        names, ref_names = set(pilot.keys()), set(reference.keys())
        pilot_values = pilot.arrays(sorted(names & ref_names), library="np")
        ref_values = reference.arrays(sorted(names & ref_names), entry_start=start, entry_stop=stop, library="np")
        fields = {k: field_differences(pilot_values[k], ref_values[k]) for k in sorted(names & ref_names)}
        topology = pilot_values["reco_topology"]
        return {
            "pilot_entries": int(pilot.num_entries),
            "reference_entries": int(reference.num_entries),
            "reference_entry_start": start, "reference_entry_stop_exclusive": stop,
            "ordered_truth_match_fields": list(TRUTH), "atol": ATOL, "rtol": RTOL,
            "producer_event_counts_histograms": [float(f["h_EventCount"].values().sum()) for f in files],
            "only_in_pilot": sorted(names - ref_names), "only_in_reference": sorted(ref_names - names),
            "common_branches": len(fields),
            "n_bit_identical_branches": sum(d.get("bit_differences") == 0 for d in fields.values()),
            "n_value_identical_branches": sum(d.get("value_differences") == 0 for d in fields.values()),
            "n_branches_within_tolerance": sum(d.get("outside_tolerance") == 0 for d in fields.values()),
            "pilot_electron_id_sums": {k: int(pilot_values[k].sum()) for k in ID_FIELDS},
            "reference_electron_id_sums": {k: int(ref_values[k].sum()) for k in ID_FIELDS},
            "pilot_topology_counts": {int(k): int(np.count_nonzero(topology == k)) for k in np.unique(topology)},
            "pilot_rows_with_proton_momentum": int(np.count_nonzero(pilot_values["reco_tagged_proton_p"] > 0)),
            "fields": fields,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--json", action="store_true", help="Print the per-branch report as JSON")
    args = parser.parse_args()
    try:
        result = check(args.pilot, args.reference)
    except ValueError as exc:
        parser.exit(1, f"{exc}\n")
    if args.json:
        print(json.dumps(result, indent=2, allow_nan=False))
    else:
        print(f"Matched {result['pilot_entries']} consecutive events at zero-based reference entries "
              f"[{result['reference_entry_start']}, {result['reference_entry_stop_exclusive']}).")
        print(f"{result['n_bit_identical_branches']}/{result['common_branches']} common branches "
              f"are bit-identical, {result['n_branches_within_tolerance']} agree within the stated tolerance.")
        for group in ("only_in_pilot", "only_in_reference"):
            if result[group]:
                print(f"{group}: {', '.join(result[group])}")
        for name, stats in result["fields"].items():
            if stats.get("bit_differences") != 0:
                print(f"{name}: {stats}")
        print("Electron-ID sums:", result["pilot_electron_id_sums"])
        print("Rows with a tagged proton:", result["pilot_rows_with_proton_momentum"])
    if (result["only_in_pilot"] or result["only_in_reference"]
            or result["n_branches_within_tolerance"] != result["common_branches"]):
        raise SystemExit(1)
