#!/usr/bin/env python3
"""Normalize a merged, unweighted-generator tuple once at generated-event level.

The supported producer assigns sigma_visible * luminosity / N_job to every
event in each input job. Merging such jobs repeats the full projected yield.
This program creates a separate ROOT file with sigma_visible * luminosity /
N_generated per row and preserves weight_proj / weight event by event. It
requires h_EventCount to confirm that the tuple contains all generated rows.
No kinematic, truth-state, reconstruction, or Q2 selection enters this step.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import positive_branching_fraction

import numpy as np
import uproot


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--cross-section-fb", required=True, type=float)
    parser.add_argument("--lumi-invfb", required=True, type=float)
    parser.add_argument("--normalization-source", default="User-supplied visible dielectron cross section")
    parser.add_argument("--generator-branching-fraction", type=positive_branching_fraction, default=None,
                        help="Generator decay fraction for provenance only; supplied visible cross section is not rescaled")
    parser.add_argument("--normalization-profile", default=None,
                        help="Name of the verified normalization profile, recorded as provenance")
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    args = parser.parse_args()
    if not np.isfinite(args.cross_section_fb) or args.cross_section_fb <= 0:
        parser.error("--cross-section-fb must be finite and positive")
    if not np.isfinite(args.lumi_invfb) or args.lumi_invfb <= 0:
        parser.error("--lumi-invfb must be finite and positive")
    source = args.input.expanduser().resolve()
    output = args.out_root.expanduser().resolve()
    audit_path = args.out_json.expanduser().resolve()
    if source == output:
        parser.error("The output must be a separate file; the source is never modified")
    output.parent.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".partial")
    with uproot.open(source) as root_file:
        if "framework_normalization" in root_file:
            raise ValueError("Input already carries framework_normalization metadata; use the original merged file")
        tree = root_file[args.tree]
        rows = int(tree.num_entries)
        if rows == 0 or "h_EventCount" not in root_file:
            raise ValueError("A nonempty tuple and h_EventCount are required to establish the generated denominator")
        generated = float(np.sum(root_file["h_EventCount"].values(flow=False)))
        if not np.isfinite(generated) or not np.isclose(generated, rows, rtol=0, atol=1e-6):
            raise ValueError("h_EventCount does not equal tuple rows; generated-event coverage needs a separate audit")
        names = list(tree.keys())
        old_weight = tree["weight"].array(library="np")
        if not np.all(np.isfinite(old_weight) & (old_weight > 0)):
            raise ValueError("The supplied global-weight producer requires finite, positive ordinary weights")
        projected = tree["weight_proj"].array(library="np") if "weight_proj" in names else None
        if projected is not None and not np.all(np.isfinite(projected) & (projected >= 0)):
            raise ValueError("Projected weights must be finite and nonnegative")
        total_expected = args.cross_section_fb * args.lumi_invfb
        corrected_weight = total_expected / generated
        raw_sum = float(np.sum(old_weight))
        raw_unique, raw_counts = np.unique(old_weight, return_counts=True)
        factors = projected / old_weight if projected is not None else np.array([], dtype=float)
        factor_unique, factor_counts = np.unique(factors, return_counts=True)
        report = {
            "input": str(source), "sha256": digest(source), "tree": args.tree,
            "entries": rows, "generated_events_h_EventCount": generated,
            "normalization_source": args.normalization_source,
            "normalization_profile": args.normalization_profile,
            "generator_branching_fraction_ee": args.generator_branching_fraction,
            "generator_branching_fraction_applied_in_preparation": False,
            "cross_section_visible_ee_fb": args.cross_section_fb,
            "luminosity_invfb": args.lumi_invfb,
            "weight_convention": "visible dielectron expected events; no branching fraction division here",
            "producer_assumption": "Unweighted generated events; each job assigned a constant sigma*L/N_job weight",
            "normalization_scope": "All generated events, before truth-state, Q2, W, and reconstruction cuts",
            "raw_sum_weight": raw_sum,
            "raw_weight_values": raw_unique.tolist(), "raw_weight_counts": raw_counts.tolist(),
            "target_sum_weight": total_expected,
            "corrected_event_weight": corrected_weight,
            "raw_sum_over_target": raw_sum / total_expected,
            "projective_factors": factor_unique.tolist(), "projective_factor_counts": factor_counts.tolist(),
            "output": str(output),
            "histogram_policy": "Stored histograms copied unchanged as raw producer diagnostics; not projected yields",
        }
        types = {name: tree[name].interpretation.numpy_dtype for name in names}
        try:
            with uproot.recreate(temporary) as normalized:
                normalized.mktree(args.tree, types, title=tree.title)
                written_sum = 0.0
                written_projected_sum = 0.0
                written_rows = 0
                for arrays in tree.iterate(names, step_size="64 MB", library="np"):
                    original = arrays["weight"]
                    replacement = np.full(original.shape, corrected_weight, dtype=np.float64)
                    if "weight_proj" in arrays:
                        arrays["weight_proj"] = (arrays["weight_proj"] / original) * replacement
                        written_projected_sum += float(np.sum(arrays["weight_proj"]))
                    arrays["weight"] = replacement
                    written_rows += len(original)
                    written_sum += float(np.sum(replacement))
                    normalized[args.tree].extend(arrays)
                for name, class_name in root_file.classnames(recursive=False, cycle=False).items():
                    if name == args.tree:
                        continue
                    if class_name.startswith(("TH1", "TH2", "TH3")) or class_name == "TObjString":
                        normalized[name] = root_file[name]
                    else:
                        report.setdefault("uncopied_objects", []).append({"name": name, "class": class_name})
                if written_rows != rows or not np.isclose(written_sum, total_expected, rtol=1e-12, atol=1e-10):
                    raise RuntimeError("Output normalization failed the row-count or total-weight check")
                report["written_sum_weight"] = written_sum
                report["written_sum_weight_proj"] = written_projected_sum
                normalized["framework_normalization"] = json.dumps(report, sort_keys=True)
            os.replace(temporary, output)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    audit_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Prepared {rows:,} events: sum(weight) {raw_sum:.10g} -> {written_sum:.10g}")
    print(f"Visible sigma = {args.cross_section_fb:g} fb; L = {args.lumi_invfb:g} fb^-1; per-event weight = {corrected_weight:.12g}")
    print(f"Normalization source: {args.normalization_source}")
    print(f"ROOT output: {output}")
    print(f"Audit: {audit_path}")


if __name__ == "__main__":
    main()
