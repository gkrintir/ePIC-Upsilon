#!/usr/bin/env python3
"""Read-only ROOT diagnostics for the state/ID/Q2 efficiency comparison.

Outputs always go to a new run directory. No framework installation is required.
Compatible with Python 3.9; configuration uses JSON rather than tomllib.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REQUIRED = ("true_mass", "true_Q2", "true_W", "true_y", "true_pT",
            "reco_mass", "reco_W", "sigma_epz")
TRUTH_BRANCHES = REQUIRED[:5]


def required_branches(settings):
    from diagnostics.proton import PROTON_BRANCHES, policy
    return REQUIRED + (PROTON_BRANCHES if policy(settings) == "required" else ())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path, framework_root, sample_keys=None):
    config = json.loads(Path(path).read_text())
    if config.get("schema_version") != 1:
        raise ValueError("Require schema_version=1")
    settings = config["settings"]
    for field in ("electron_energy", "true_mass_window", "q2_split"):
        if not math.isfinite(settings[field]) or settings[field] <= 0:
            raise ValueError("Require finite positive settings." + field)
    if not 0 < settings["coverage"] < 1:
        raise ValueError("coverage must be between 0 and 1")
    if settings["min_denominator"] < 1:
        raise ValueError("min_denominator must be positive")
    for field in ("w_range", "pt_edges", "w_edges", "y_edges", "q2_edges",
                  "tail_edges", "strata_w_edges", "strata_y_edges"):
        values = settings[field]
        if len(values) < 2 or not all(math.isfinite(x) for x in values):
            raise ValueError("Invalid edges: " + field)
        if any(b <= a for a, b in zip(values[:-1], values[1:])):
            raise ValueError("Edges must increase: " + field)
    if len(settings["w_range"]) != 2:
        raise ValueError("w_range requires two endpoints")
    samples = config["samples"]
    keys = [sample["key"] for sample in samples]
    if len(set(keys)) != len(keys) or not samples:
        raise ValueError("Sample keys must be unique and nonempty")
    if sample_keys:
        unknown = set(sample_keys) - set(keys)
        if unknown:
            raise ValueError("Unknown samples: " + ", ".join(sorted(unknown)))
        samples = [sample for sample in samples if sample["key"] in sample_keys]
    for sample in samples:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", sample["key"]):
            raise ValueError("Invalid sample key")
        if sample["state"] not in ("1S", "2S", "3S"):
            raise ValueError("Unsupported state")
        if not math.isfinite(sample["q2_max"]) or sample["q2_max"] <= 0:
            raise ValueError("q2_max must be finite and positive")
        for field in ("path", "reference_efficiency_csv"):
            if sample.get(field):
                item = Path(sample[field]).expanduser()
                sample[field] = str((framework_root / item).resolve() if not item.is_absolute() else item)
    config["samples"] = samples
    config["tree"] = config.get("tree", "UpsilonTree")
    return config


def import_dependencies():
    try:
        import numpy as np
        import uproot
        return np, uproot
    except ImportError as error:
        raise RuntimeError("Missing dependency. In the environment used for this script, run "
                           "python3 -m pip install -r " + str(HERE / "requirements.txt")) from error


def inspect_sample(sample, tree_name, uproot, settings=None, validate_protons=False):
    settings = settings or {}
    required = required_branches(settings)
    path = Path(sample["path"])
    if not path.is_file():
        raise FileNotFoundError(str(path))
    stat = path.stat()
    with uproot.open(path) as source:
        if tree_name not in source:
            raise ValueError("Missing tree " + tree_name)
        tree = source[tree_name]
        if tree.num_entries == 0:
            raise ValueError("Input tree contains no rows")
        names = set(tree.keys())
        missing = set(required) - names
        if missing:
            raise ValueError("Missing branches: " + ", ".join(sorted(missing)))
        if validate_protons and settings.get("proton_selection", "required") == "required":
            from diagnostics.proton import PROTON_BRANCHES, proton_mask
            for arrays in tree.iterate(expressions=list(PROTON_BRANCHES), step_size="32 MB", library="np"):
                proton_mask(arrays, settings, len(arrays[PROTON_BRANCHES[0]]))
        return {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                "entries": int(tree.num_entries), "branches": sorted(names),
                "required_branches": list(required), "weight_branches_used": []}


def check_reference(rows, path):
    """Compare integer counts to an old performance CSV without changing it."""
    if not path:
        return {"status": "NOT_CONFIGURED"}
    if not Path(path).is_file():
        return {"status": "UNAVAILABLE", "path": path,
                "detail": "Diagnostic still runs; original numerical cross-check unavailable."}
    with Path(path).open(newline="") as stream:
        previous = list(csv.DictReader(stream))
    nominal = [row for row in rows if row["q2_band"] == "all"]
    errors = []
    compared = 0
    for old in previous:
        matches = [row for row in nominal if row["variable"] == old["variable"]
                   and math.isclose(row["low"], float(old["low"]), rel_tol=0, abs_tol=1e-10)
                   and math.isclose(row["high"], float(old["high"]), rel_tol=0, abs_tol=1e-10)]
        if len(matches) != 1:
            errors.append({"variable": old["variable"], "low": old["low"], "reason": "bin mismatch"})
            continue
        new = matches[0]
        for field in ("denominator", "selected"):
            if int(new[field]) != int(old[field]):
                errors.append({"variable": old["variable"], "low": old["low"], "field": field,
                               "saved": int(old[field]), "recomputed": int(new[field])})
        compared += 1
    if len(previous) != len(nominal):
        errors.append({"reason": "different number of bins", "saved": len(previous), "new": len(nominal)})
    return {"status": "MATCH" if not errors else "MISMATCH", "path": path,
            "sha256": sha256(path), "compared_bins": compared, "differences": errors,
            "interpretation": "Unweighted selected and eligible counts; mismatch can indicate different input or settings."}


def write_tables(output, tables):
    folder = output / "tables"
    folder.mkdir()
    for name, rows in tables.items():
        if not rows:
            continue
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (folder / (name + ".csv")).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def compare_truth(audits, samples):
    one_s = [sample for sample in samples if sample["state"] == "1S" and sample["key"] in audits
             and audits[sample["key"]].get("status") == "DONE"]
    comparisons = []
    for i, left in enumerate(one_s):
        for right in one_s[i + 1:]:
            a, b = audits[left["key"]], audits[right["key"]]
            same = a["input"]["entries"] == b["input"]["entries"] and a["truth_content_sha256"] == b["truth_content_sha256"]
            comparisons.append({"left": left["key"], "right": right["key"],
                                "same_ordered_truth_values": same,
                                "interpretation": ("All five stored truth columns agree entry-by-entry. " if same else
                                    "Ordered truth content differs; order or generated populations can differ. ") +
                                "No assumption of independent or paired ID uncertainties is made. "
                                "This does not authenticate that ID was the only production change."})
    return comparisons


def run(config, args, root):
    np, uproot = import_dependencies()
    from diagnostics.core import Accumulator
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_id or ""):
        raise ValueError("Provide --run-id containing letters, digits, _, . or -, starting with a letter/digit")
    run_dir = root / "results/generated/runs" / args.run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    output = run_dir / "artifacts/efficiency_diagnostics"
    output.mkdir(parents=True)
    config["framework_root"] = str(root)
    write_json(run_dir / "effective_config.json", config)
    manifest = {"run_id": args.run_id, "started_utc": datetime.now(timezone.utc).isoformat(),
                "status": "RUNNING", "samples": {}, "synthetic": bool(config.get("synthetic", False)),
                "source_hashes": {str(path.relative_to(HERE)): sha256(path)
                                  for path in sorted(HERE.rglob("*.py")) if "__pycache__" not in path.parts},
                "output": str(output), "interpretation": "Execution status is not a physics validation status."}
    write_json(run_dir / "run_manifest.json", manifest)
    tables, audits = {}, {}
    for sample in config["samples"]:
        key = sample["key"]
        print("LOAD", key, sample["path"], flush=True)
        try:
            metadata = inspect_sample(sample, config["tree"], uproot, config["settings"])
            if args.hash_inputs:
                print("HASH", key, flush=True)
                metadata["sha256"] = sha256(sample["path"])
            acc = Accumulator(sample, config["settings"])
            truth_hash = {branch: hashlib.sha256() for branch in TRUTH_BRANCHES}
            read = 0
            with uproot.open(sample["path"]) as source:
                for chunk in source[config["tree"]].iterate(expressions=list(required_branches(config["settings"])), step_size=args.step_size, library="np"):
                    arrays = ({name: chunk[name] for name in chunk.dtype.names}
                              if isinstance(chunk, np.ndarray) and chunk.dtype.names else dict(chunk))
                    acc.consume(arrays)
                    for branch in TRUTH_BRANCHES:
                        truth_hash[branch].update(np.asarray(arrays[branch], dtype="<f8").tobytes())
                    read += len(arrays["true_mass"])
                    print("READ", key, read, "/", metadata["entries"], flush=True)
            after = Path(sample["path"]).stat()
            if after.st_size != metadata["size_bytes"] or after.st_mtime_ns != metadata["mtime_ns"]:
                raise RuntimeError("Input changed during reading; do not interpret this run")
            if read != metadata["entries"]:
                raise RuntimeError("Read entry count does not match ROOT metadata")
            result = acc.finish()
            if result["audit"]["band_counts"]["all"]["eligible_truth"] == 0:
                raise ValueError("No eligible truth rows under the configured state/Q2/W selection")
            reference = check_reference(result["tables"]["inclusive_efficiency"], sample.get("reference_efficiency_csv"))
            reference["comparison_category"] = "inclusive_before_proton_selection"
            reference["validates_proton_tagged_numerator"] = False
            audit = dict(result["audit"], input=metadata, status="DONE", reference_efficiency_check=reference,
                         truth_content_sha256={name: digest.hexdigest() for name, digest in truth_hash.items()})
            audits[key] = audit
            for name, rows in result["tables"].items():
                tables.setdefault(name, []).extend(rows)
            manifest["samples"][key] = {"status": "DONE", "entries": read, "reference_check": reference["status"]}
            print("DONE", key, "reference", reference["status"], flush=True)
        except Exception as error:
            audits[key] = {"status": "FAILED", "error": str(error), "path": sample["path"]}
            manifest["samples"][key] = {"status": "FAILED", "error": str(error)}
            traceback.print_exc()
        write_json(output / "sample_audits.json", audits)
        write_json(run_dir / "run_manifest.json", manifest)
    if tables:
        write_tables(output, tables)
    comparisons = compare_truth(audits, config["samples"])
    write_json(output / "input_comparability.json", comparisons)
    failed = any(audit["status"] != "DONE" or audit.get("reference_efficiency_check", {}).get("status") == "MISMATCH"
                 for audit in audits.values())
    if tables:
        try:
            from diagnostics.report import render_report
            config["execution_status_before_render"] = "PARTIAL" if failed else "COMPLETE"
            manifest["rendered_files"] = render_report(output, config, config["samples"], audits)
        except Exception as error:
            failed = True
            manifest["render_error"] = str(error)
            traceback.print_exc()
    else:
        failed = True
    manifest["status"] = "PARTIAL" if failed else "COMPLETE"
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["causal_validation_status"] = "NOT_ESTABLISHED_BY_EXECUTION"
    write_json(run_dir / "run_manifest.json", manifest)
    write_json(output / "diagnostic_status.json", {
        "execution_status": manifest["status"], "causal_validation_status": manifest["causal_validation_status"],
        "reference_checks": {key: audit.get("reference_efficiency_check", {"status": audit["status"]}) for key, audit in audits.items()},
        "limitations": ["Efficiencies are conditional on eligible stored tuple rows.",
                        "No luminosity, cross section, branching fraction or weight branch enters these ratios.",
                        "Different ROOT productions need metadata checks before attributing a change solely to electron ID.",
                        "Tuple cutflows do not identify individual upstream tracking or cluster-association failures.",
                        "Displayed errors are individual finite-MC Wilson intervals; no ID-difference significance is inferred."]})
    print("Run status:", manifest["status"], flush=True)
    print("Outputs:", output, flush=True)
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "validate", "run"))
    parser.add_argument("--framework-root", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=HERE / "config_9x275.json")
    parser.add_argument("--run-id")
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--hash-inputs", action="store_true")
    parser.add_argument("--proton-selection", choices=("required", "inclusive"), default="required", help="Default requires an explicitly validated reconstructed proton; inclusive reproduces the old comparison category")
    parser.add_argument("--step-size", default="64 MB", help="Uproot chunk size, for example '64 MB'")
    args = parser.parse_args(argv)
    try:
        root = args.framework_root.expanduser().resolve()
        config = load_config(args.config, root, args.samples)
        config["settings"]["proton_selection"] = args.proton_selection
        if args.action == "run":
            if not args.run_id:
                parser.error("run requires --run-id; existing run directories are never overwritten")
            return run(config, args, root)
        uproot = import_dependencies()[1] if args.action == "validate" else None
        failed = False
        for sample in config["samples"]:
            try:
                if uproot:
                    info = inspect_sample(sample, config["tree"], uproot, config["settings"], validate_protons=True)
                    print("READY", sample["key"], info["entries"], "rows")
                else:
                    print("INPUT", sample["key"], sample["path"])
                print("     ", sample["state"], sample["id_mode"], "Q2 <=", sample["q2_max"], "GeV^2", "proton selection", args.proton_selection)
            except Exception as error:
                failed = True
                print("BLOCKED", sample["key"], str(error))
        print("Read-only inputs. Fresh run output: results/generated/runs/<run-id>/artifacts/efficiency_diagnostics")
        return 1 if failed else 0
    except Exception as error:
        print("ERROR:", str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
