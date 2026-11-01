#!/usr/bin/env python3
"""Run all three pure-state workflows at 9x130 and 9x275 with separate caches."""
import argparse
import datetime as dt
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from proton_selection import proton_selection_metadata, guard_proton_run
from proton_category_isolation import enforce_category_workspace, safe_component

ROOT = Path(__file__).resolve().parents[1]
STATES = ("1S", "2S", "3S")


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be positive and finite")
    return number


def state_list(value):
    values = list(STATES) if value == "all" else [x.strip().upper() for x in value.split(",")]
    if not values or any(x not in STATES for x in values) or len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("use all or a comma-separated subset of 1S,2S,3S")
    return values


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "plan", "validate", "list"), nargs="?", default="run")
    parser.add_argument("--states", type=state_list, default=list(STATES))
    parser.add_argument("--energy", "--energies", choices=("all", "9x130", "9x275"), default="all")
    parser.add_argument("--batch-id", help="Stable batch prefix; state suffixes -1s, -2s and -3s are added")
    parser.add_argument("--stat-mode", choices=("mc", "projected"), default="projected")
    parser.add_argument("--projection-lumi", type=positive, default=10.)
    parser.add_argument("--tuple-lumi", type=positive, default=10.)
    parser.add_argument("--config", help="Optional TOML overrides applied to each selected state")
    parser.add_argument("--set", dest="set_values", action="append", default=[], help="Additional KEY=VALUE override")
    parser.add_argument("--targets", nargs="+", help="Explicit task/group targets instead of the default full state workflow")
    parser.add_argument("--hash-inputs", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-going", action="store_true", help="Accepted explicitly; independent tasks and states always continue")
    return parser


def command_for(args, state, batch_id):
    command = [sys.executable, str(ROOT / "workflow/run_state.py"), "--state", state]
    if args.config:
        command += ["--config", str(Path(args.config).expanduser().resolve())]
    overrides = ["analysis.stat_mode=" + args.stat_mode,
                 "analysis.projection_lumi_invfb=" + str(args.projection_lumi),
                 "analysis.tuple_lumi_invfb=" + str(args.tuple_lumi),
                 "analysis.lumi_invfb=" + str(args.tuple_lumi)] + args.set_values
    for value in overrides:
        command += ["--set", value]
    command += [args.command]
    if args.command != "list":
        targets = args.targets or ["state_all" if args.energy == "all" else "state_" + args.energy]
        command += targets + ["--run-id", batch_id + "-" + state.lower()]
    if args.command == "run":
        command += ["--keep-going"]
        if args.hash_inputs:
            command += ["--hash-inputs"]
        if args.force:
            command += ["--force"]
        if args.dry_run:
            command += ["--dry-run"]
    return command


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.targets and args.energy != "all":
        parser.error("Use either --targets or --energy; an explicit target may include either energy")
    if args.command != "run" and (args.force or args.dry_run):
        parser.error("--force and --dry-run apply only to run")
    batch_id = args.batch_id or ("jul26-allstates-" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,149}", batch_id):
        parser.error("--batch-id must contain letters, digits, '.', '_' or '-' and be at most 150 characters")
    manifest = {"proton_selection": proton_selection_metadata(), "batch_id": batch_id, "command": args.command, "states": args.states,
                "energy": args.energy, "stat_mode": args.stat_mode,
                "projection_lumi_invfb": args.projection_lumi, "tuple_lumi_invfb": args.tuple_lumi,
                "phase_space_note": "1S Q2max=1 GeV^2; 2S/3S Q2max=0.01 GeV^2. State yields have different generated Q2 ranges.",
                "input_defaults": {"1S": "hiAcc/output_newID/merged_Upsilon_final.root",
                                   "2S": "hiAcc/output/merged_Upsilon_final.root",
                                   "3S": "hiAcc/output/merged_Upsilon_final.root"},
                "input_provenance": "Actual inputs, including any explicit overrides, are recorded in each state run manifest.",
                "runs": [], "status": "RUNNING"}
    safe_component(batch_id, "batch ID")
    batch_dir = ROOT / "results/generated/batches" / batch_id
    enforce_category_workspace(ROOT, batch_dir)
    write_manifest = args.command == "run" and not args.dry_run
    if write_manifest:
        guard_proton_run({"run": {"dir": str(batch_dir)},
                          "proton_run_context": dict(vars(args))}, ROOT)
        batch_dir.mkdir(parents=True, exist_ok=True)
    for state in args.states:
        command = command_for(args, state, batch_id)
        print("\nUpsilon(%s): %s" % (state, shlex.join(command)), flush=True)
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        result = subprocess.run(command, cwd=ROOT)
        run_id = batch_id + "-" + state.lower()
        record = {"state": state, "run_id": run_id, "command": command,
                  "returncode": result.returncode, "started_utc": started,
                  "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                  "run_manifest": str(ROOT / "results/generated/runs" / run_id / "run_manifest.json"),
                  "artifacts": str(ROOT / "results/generated/runs" / run_id / "artifacts")}
        manifest["runs"].append(record)
        if write_manifest:
            (batch_dir / "batch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    failed = [r["state"] for r in manifest["runs"] if r["returncode"] != 0]
    manifest["status"] = "PARTIAL" if failed else "DRY_RUN" if args.dry_run else "COMPLETE"
    manifest["incomplete_states"] = failed
    if write_manifest:
        path = batch_dir / "batch_manifest.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n")
        print("Batch manifest: " + str(path))
    print("Batch " + manifest["status"] + (": " + ", ".join(failed) if failed else ""), flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
