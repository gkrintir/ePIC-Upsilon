#!/usr/bin/env python3
"""Four standalone truth-coordinate efficiencies after the final selection.

Run from lxplus or SWAN. A new run ID is mandatory and never overwritten.
This entry point does not execute the diagnostic dashboards or extraction FW.
"""
import argparse
import csv
import hashlib
import json
import re
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from run_efficiency_diagnostics import (
    HERE, REQUIRED, TRUTH_BRANCHES, import_dependencies, inspect_sample,
    load_config, required_branches, sha256, write_json,
)

MAP_NAMES = ("pt_q2", "pt_y", "xv_q2", "xb_q2")


def execute(config, args, root):
    np, uproot = import_dependencies()
    from diagnostics.final_maps import FinalMapAccumulator
    from diagnostics.final_maps_report import render_final_maps
    sample = config["samples"][0]
    # Validate configuration and schema before allocating a run directory.
    acc = FinalMapAccumulator(sample, config["settings"])
    metadata = inspect_sample(sample, config["tree"], uproot, config["settings"])
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_id or ""):
        raise ValueError("Provide a fresh --run-id using letters, digits, _, . or -")
    run_dir = root / "results/generated/runs" / args.run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    output = run_dir / "artifacts/final_efficiency_2d" / sample["key"]
    output.mkdir(parents=True)
    config["framework_root"] = str(root)
    config["workflow"] = "final_efficiency_2d"
    write_json(run_dir / "effective_config.json", config)
    manifest = {
        "run_id": args.run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "RUNNING", "sample": sample["key"],
        "synthetic": bool(config.get("synthetic", False)),
        "output": str(output), "input": metadata,
        "source_hashes": {str(p.relative_to(HERE)): sha256(p)
                          for p in sorted(HERE.rglob("*.py"))
                          if "__pycache__" not in p.parts},
        "maps": list(MAP_NAMES),
        "interpretation": "Execution status is not a physics-validation claim.",
    }
    write_json(run_dir / "run_manifest.json", manifest)
    try:
        if args.hash_inputs:
            print("HASH", sample["path"], flush=True)
            metadata["sha256"] = sha256(sample["path"])
        truth_hash = {name: hashlib.sha256() for name in TRUTH_BRANCHES}
        read = 0
        print("LOAD", sample["key"], sample["path"], flush=True)
        with uproot.open(sample["path"]) as source:
            for chunk in source[config["tree"]].iterate(
                    expressions=list(required_branches(config["settings"])), step_size=args.step_size, library="np"):
                arrays = ({name: chunk[name] for name in chunk.dtype.names}
                          if isinstance(chunk, np.ndarray) and chunk.dtype.names else dict(chunk))
                acc.consume(arrays)
                for name in TRUTH_BRANCHES:
                    truth_hash[name].update(np.asarray(arrays[name], dtype="<f8").tobytes())
                read += len(arrays["true_mass"])
                print("READ", read, "/", metadata["entries"], flush=True)
        after = Path(sample["path"]).stat()
        if after.st_size != metadata["size_bytes"] or after.st_mtime_ns != metadata["mtime_ns"]:
            raise RuntimeError("Input changed during reading; do not interpret this run")
        if read != metadata["entries"]:
            raise RuntimeError("Read entry count disagrees with ROOT metadata")
        result = acc.finish()
        audit = result["audit"]
        audit["input"] = metadata
        audit["truth_content_sha256"] = {name: h.hexdigest() for name, h in truth_hash.items()}
        audit["synthetic"] = manifest["synthetic"]
        audit["config"] = {"sample": sample, "settings": config["settings"]}
        write_json(output / "selection_audit.json", audit)
        table_dir = output / "tables"
        table_dir.mkdir()
        rows = result["rows"]
        if not rows:
            raise RuntimeError("No map rows were produced")
        with (table_dir / "final_efficiency_2d.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        if not any(r["denominator"] > 0 for r in rows):
            raise ValueError("No eligible events are inside any configured map; inspect selection_audit.json")
        plots = render_final_maps(rows, sample, config["settings"], output / "plots",
                                  synthetic=manifest["synthetic"])
        if len(plots) != 8:
            raise RuntimeError("Expected exactly four PDF and four PNG outputs")
        manifest["artifacts"] = {str(Path(p).relative_to(run_dir)): sha256(p) for p in plots}
        manifest["status"] = "COMPLETE"
        print("DONE", "four final-selection efficiency maps", flush=True)
    except Exception as error:
        manifest["status"] = "FAILED"
        manifest["error"] = str(error)
        traceback.print_exc()
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(run_dir / "run_manifest.json", manifest)
    print("Run status:", manifest["status"])
    print("Outputs:", output)
    return 0 if manifest["status"] == "COMPLETE" else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "validate", "run"))
    parser.add_argument("--framework-root", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=HERE / "config_9x275.json")
    parser.add_argument("--sample", default="1s_newid",
                        help="One sample per run, default 1s_newid; exactly four plots")
    parser.add_argument("--run-id")
    parser.add_argument("--hash-inputs", action="store_true")
    parser.add_argument("--proton-selection", choices=("required", "inclusive"), default="required", help="Explicit reconstructed-proton category; default required")
    parser.add_argument("--step-size", default="64 MB")
    args = parser.parse_args(argv)
    try:
        root = args.framework_root.expanduser().resolve()
        if not root.is_dir():
            raise ValueError("Framework root is not a directory: " + str(root))
        config = load_config(args.config.expanduser().resolve(), root, [args.sample])
        config["settings"]["proton_selection"] = args.proton_selection
        if args.action == "plan":
            sample = config["samples"][0]
            print(json.dumps({"sample": sample, "maps": list(MAP_NAMES),
                              "selection": "nominal final selection, no Q2 bands or comparison panels",
                              "proton_selection": args.proton_selection,
                              "output": "results/generated/runs/<run-id>/artifacts/final_efficiency_2d/"
                                        + sample["key"]}, indent=2))
            return 0
        if args.action == "validate":
            from diagnostics.final_maps import FinalMapAccumulator
            _, uproot = import_dependencies()
            sample = config["samples"][0]
            FinalMapAccumulator(sample, config["settings"])
            metadata = inspect_sample(sample, config["tree"], uproot, config["settings"], validate_protons=True)
            print("PASS", sample["key"], "entries", metadata["entries"])
            return 0
        return execute(config, args, root)
    except Exception as error:
        print("ERROR:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
