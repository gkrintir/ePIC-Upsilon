#!/usr/bin/env python3
"""Check a pure-state tuple and its generated denominator before normalization.

This catches mismatched input files/card ranges. It does not authenticate the
production generator or establish its absolute cross section.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import uproot

MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}


def inspect_sample(path, tree_name, state, q2_max, mass_window=0.2,
                   min_state_fraction=0.95, max_q2_outlier_fraction=0.001):
    result = {"input": str(Path(path).expanduser().resolve()), "tree": tree_name,
              "state": state, "state_mass_GeV": MASSES[state], "q2_max_GeV2": q2_max,
              "truth_mass_window_GeV": mass_window,
              "min_state_fraction": min_state_fraction,
              "max_q2_outlier_fraction": max_q2_outlier_fraction,
              "scope": "Input state and generated-range sanity checks; production provenance unverified",
              "errors": [], "warnings": []}
    with uproot.open(path) as source:
        if tree_name not in source:
            result["errors"].append("Missing tree " + tree_name)
            result["status"] = "FAIL"
            return result
        tree = source[tree_name]
        entries = int(tree.num_entries)
        result["entries"] = entries
        if entries == 0:
            result["errors"].append("Generated tuple is empty")
        if "h_EventCount" not in source:
            result["errors"].append("Missing h_EventCount generated denominator")
        else:
            histogram = source["h_EventCount"]
            regular = np.asarray(histogram.values(flow=False), dtype=float)
            flow = np.asarray(histogram.values(flow=True), dtype=float)
            generated = float(regular.sum())
            result["generated_events_h_EventCount"] = generated if math.isfinite(generated) else None
            if not np.all(np.isfinite(flow)) or not math.isclose(generated, entries, rel_tol=0, abs_tol=1e-6):
                result["errors"].append("h_EventCount does not match all generated tuple rows")
            if flow.ndim != 1 or flow.size != regular.size + 2 or np.any(flow[[0, -1]] != 0):
                result["errors"].append("h_EventCount has nonzero or unrecognized flow bins")
        missing = sorted({"true_mass", "true_Q2"} - set(tree.keys()))
        if missing:
            result["errors"].append("Missing branches: " + ", ".join(missing))
        elif entries:
            state_counts = {key: 0 for key in MASSES}
            invalid_mass = 0
            q2_invalid = 0
            q2_above = 0
            q2_observed_max = None
            for arrays in tree.iterate(["true_mass", "true_Q2"], step_size="32 MB", library="np"):
                mass = np.asarray(arrays["true_mass"], dtype=float)
                q2 = np.asarray(arrays["true_Q2"], dtype=float)
                invalid_mass += int((~np.isfinite(mass)).sum())
                for key, value in MASSES.items():
                    state_counts[key] += int((np.isfinite(mass) & (np.abs(mass - value) <= mass_window)).sum())
                finite = np.isfinite(q2)
                q2_invalid += int((~finite | (q2 < 0)).sum())
                q2_above += int((finite & (q2 > q2_max * (1 + 1e-6))).sum())
                if finite.any():
                    maximum = float(q2[finite].max())
                    q2_observed_max = maximum if q2_observed_max is None else max(q2_observed_max, maximum)
            state_fraction = state_counts[state] / entries
            outlier_fraction = (q2_invalid + q2_above) / entries
            result.update({"counts_within_each_state_window": state_counts,
                           "selected_state_fraction": state_fraction,
                           "nonfinite_truth_mass_count": invalid_mass,
                           "q2_invalid_count": q2_invalid, "q2_above_card_count": q2_above,
                           "q2_outlier_fraction": outlier_fraction, "q2_observed_max_GeV2": q2_observed_max})
            if state_fraction < min_state_fraction:
                result["errors"].append("Tuple is inconsistent with the selected pure state; check the input path")
            if outlier_fraction > max_q2_outlier_fraction:
                result["errors"].append("Truth Q2 population is inconsistent with the assigned generation range; check the production card before assigning its rate")
            elif q2_invalid + q2_above:
                result["warnings"].append("Small truth-Q2 outlier population retained in the generated denominator and excluded by downstream event cuts")
    result["status"] = "FAIL" if result["errors"] else "WARN" if result["warnings"] else "PASS"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--state", required=True, choices=sorted(MASSES))
    parser.add_argument("--q2-max", type=float, required=True)
    parser.add_argument("--mass-window", type=float, default=0.2)
    parser.add_argument("--min-state-fraction", type=float, default=0.95)
    parser.add_argument("--max-q2-outlier-fraction", type=float, default=0.001)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    if not math.isfinite(args.q2_max) or args.q2_max <= 0:
        parser.error("--q2-max must be positive and finite")
    if not math.isfinite(args.mass_window) or args.mass_window <= 0:
        parser.error("--mass-window must be positive and finite")
    if not 0 < args.min_state_fraction <= 1 or not 0 <= args.max_q2_outlier_fraction < 1:
        parser.error("Invalid state/Q2 population thresholds")
    result = inspect_sample(args.input, args.tree, args.state, args.q2_max, args.mass_window,
                            args.min_state_fraction, args.max_q2_outlier_fraction)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(result["status"] + " " + args.state + " input: " + args.input)
    for error in result["errors"]:
        print("  " + error)
    return 2 if result["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
