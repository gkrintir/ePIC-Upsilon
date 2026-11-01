#!/usr/bin/env python3
"""Prevent statistics runs from reusing artifacts with another statistical contract.

The legacy runner caches by output existence. This always-run prerequisite
guards configuration and analysis source changes without modifying that runner.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

if __package__:
    from .normalization_contract import validate_flat_normalization, task_definitions_sha256, configuration_json_sha256
else:
    # Support the existing importlib-based guard checks as well as the CLI.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from normalization_contract import validate_flat_normalization, task_definitions_sha256, configuration_json_sha256


def build_contract(settings, project_root, task_definitions=None):
    values = {}
    for setting in settings:
        key, sep, value = setting.partition("=")
        if not sep or not key or key in values:
            raise ValueError("Each --setting must have a unique KEY=VALUE")
        values[key] = value
    validate_flat_normalization(values)
    mode = values.get("analysis.stat_mode")
    if mode not in ("mc", "projected"):
        raise ValueError("analysis.stat_mode must be mc or projected")
    for key in ("analysis.lumi_invfb", "analysis.tuple_lumi_invfb", "analysis.projection_lumi_invfb"):
        value = float(values[key])
        if not math.isfinite(value) or value <= 0:
            raise ValueError(key + " must be positive and finite")
    if float(values["analysis.lumi_invfb"]) != float(values["analysis.tuple_lumi_invfb"]):
        raise ValueError("Keep lumi_invfb equal to tuple_lumi_invfb in this pipeline; set projection_lumi_invfb for a data forecast")
    if mode == "mc" and int(values["analysis.include_efficiency_mc_uncertainty"]):
        raise ValueError("Efficiency MC cannot be added again to same-sample MC closure")
    if values.get("analysis.count_interval", "median_garwood") not in ("median_garwood", "standard_deviation"):
        raise ValueError("analysis.count_interval must be median_garwood or standard_deviation")
    interval_coverage = float(values.get("analysis.interval_coverage", math.erf(1 / math.sqrt(2))))
    if not math.isfinite(interval_coverage) or not 0 < interval_coverage < 1:
        raise ValueError("analysis.interval_coverage must be strictly between zero and one")
    sources = {}
    for folder in ("scripts", "workflow"):
        for path in sorted((project_root / folder).rglob("*.py")):
            sources[str(path.relative_to(project_root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    contract = {"schema_version": 1, "settings": values, "analysis_source_sha256": sources,
                "configuration_json_sha256": configuration_json_sha256(project_root)}
    if task_definitions is not None:
        contract["task_definitions_sha256"] = task_definitions_sha256(task_definitions)
    return contract


def check_contract(run_dir, contract):
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / "statistical_contract.json"
    if path.exists():
        if json.loads(path.read_text()) != contract:
            raise ValueError("Statistical settings or analysis code changed. Use a new --run-id to avoid cached results from another calculation")
        return "MATCH"
    if (run_dir / "run_manifest.json").exists() or any((run_dir / "artifacts").glob("**/*")):
        raise ValueError("This run contains results without a statistical contract. Use a new --run-id")
    with path.open("x") as handle:
        json.dump(contract, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return "CREATED"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--setting", action="append", default=[])
    args = parser.parse_args()
    contract = build_contract(args.setting, args.project_root)
    print("Statistical contract", check_contract(args.run_dir, contract))


if __name__ == "__main__":
    main()
