#!/usr/bin/env python3
"""Render existing numerical reoptimization results into a new report directory.

No ROOT files are opened, no fit is repeated, and no validation gate is changed.
The destination must be new and outside the original run directory.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

FRAMEWORK = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FRAMEWORK / "scripts"))
sys.path.insert(0, str(FRAMEWORK / "workflow"))
from proton_category_isolation import enforce_category_workspace
from reoptimization.reporting import write_report


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_tagged_saved_results(results, config):
    """Do not relabel historical selected results as this tagged category."""
    policy = config.get("proton_selection")
    if not isinstance(policy, dict) or policy.get("mode") != "required":
        raise ValueError("Category report refresh requires saved required proton-selection provenance; rerun the numerical workflow")
    for state, observables in results.items():
        if not isinstance(observables, dict):
            raise ValueError("Invalid saved reoptimization results for " + str(state))
        for observable, entry in observables.items():
            label = str(state) + "/" + str(observable)
            if not isinstance(entry, dict):
                raise ValueError("Invalid saved result for " + label)
            # Input/optimization errors contain no selected numerical results.
            failed_without_results = (entry.get("status") in ("INPUT_ERROR", "OPTIMIZATION_ERROR")
                                      and not any(entry.get(k) for k in ("model", "validation_model", "inference", "selection_metadata")))
            if failed_without_results:
                continue
            if entry.get("proton_selection") != policy or entry.get("selection_metadata", {}).get("proton_selection") != policy:
                raise ValueError("Missing or inconsistent tagged selection provenance in " + label)
            for key in ("model", "validation_model"):
                model = entry.get(key)
                if model is None:
                    if entry.get("inference"):
                        raise ValueError("Tagged report is missing response provenance in " + label + "/" + key)
                    continue
                if not isinstance(model, dict) or model.get("metadata", {}).get("proton_selection") != policy:
                    raise ValueError("Missing or inconsistent tagged response provenance in " + label + "/" + key)
                coarse = model.get("coarse")
                if not isinstance(coarse, dict) or coarse.get("metadata", {}).get("proton_selection") != policy:
                    raise ValueError("Missing or inconsistent tagged coarse-response provenance in " + label + "/" + key)


def refresh(input_run, output_dir):
    category = enforce_category_workspace(FRAMEWORK, output_dir)
    if category is not None:
        enforce_category_workspace(FRAMEWORK, input_run)
    source = Path(input_run).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    if destination == source or source in destination.parents:
        raise ValueError("Output must be outside the original run directory")
    if destination.exists():
        raise FileExistsError("Output directory already exists; choose a fresh report directory")
    results_path = source / "artifacts/reoptimization/all_results.json"
    config_path = source / "effective_config.json"
    results = json.loads(results_path.read_text())
    config = json.loads(config_path.read_text())
    if not isinstance(results, dict) or not isinstance(config, dict):
        raise ValueError("Expected object-valued saved results and effective configuration")
    if category is not None:
        _require_tagged_saved_results(results, config)
    source_paths = [results_path, config_path]
    original_manifest = source / "run_manifest.json"
    if original_manifest.is_file():
        source_paths.append(original_manifest)
    source_hashes = {str(path): sha(path) for path in source_paths}
    destination.mkdir(parents=True, exist_ok=False)
    report = write_report(destination, results, config)
    if any(sha(Path(path)) != digest for path, digest in source_hashes.items()):
        raise RuntimeError("Source changed while rendering; this report must not be used")
    record = {
        "operation": "report_refresh_only",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_run": str(source),
        "source_sha256": source_hashes,
        "output_dir": str(destination),
        "reporting_code_sha256": sha(FRAMEWORK / "scripts/reoptimization/reporting.py"),
        "numerical_results_recomputed": False,
        "validation_statuses_changed": False,
        "publication_ready": False,
        "interpretation": "Reporting repair only. Original fit results and validation outcomes are unchanged. Legacy truth_shape stress is labeled as a combined shape-and-feed-in test.",
        "report": report,
    }
    (destination / "refresh_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-run", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    result = refresh(args.input_run, args.output_dir)
    print("Report refreshed without rerunning inference:", result["output_dir"])
    print("status/overall_status retain final case outcomes; stage_status records component outcomes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
