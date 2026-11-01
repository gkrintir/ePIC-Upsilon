#!/usr/bin/env python3
"""Run independent fEMC note scenarios, continuing after missing/failed inputs.

The existing comparison launcher and physics normalization are unchanged.
Each child uses the statistics continuation entry point, which skips blocked
tasks while retaining the projected unfolding response-support guard.
"""

import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys


DEFAULT_BASE_DIR = "/eos/cms/store/group/phys_top/gkrintir/upsilon_data/Upsilon1S"
BEAMS = ("9x130", "9x275")
SCENARIO_DIRS = {"full": "output", "remove_outer": "output_remove_outer",
                 "remove_inner": "output_remove_inner"}
SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def utc_now():
    return dt.datetime.now(dt.timezone.utc)


def positive_luminosity(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError("Luminosity must be a finite positive number")
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("Luminosity must be a finite positive number")
    return number


def safe_batch_id(value):
    if len(value) > 64 or not SAFE_ID.fullmatch(value):
        raise argparse.ArgumentTypeError("Batch ID must start with a letter or digit and contain at most 64 letters, digits, '.', '_' or '-'")
    return value


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beam", choices=("all",) + BEAMS, default="all")
    parser.add_argument("--scenario", choices=("all",) + tuple(SCENARIO_DIRS), default="all")
    parser.add_argument("--stat-mode", choices=("mc", "projected"),
                        default=os.environ.get("STAT_MODE", "projected"))
    parser.add_argument("--projection-lumi", type=positive_luminosity,
                        default=os.environ.get("PROJECTION_LUMI", "10"),
                        help="Projected luminosity in inverse fb; source tuple luminosity is unchanged")
    parser.add_argument("--batch-id", type=safe_batch_id,
                        default=utc_now().strftime("%Y%m%dT%H%M%SZ"))
    parser.add_argument("--base-dir", type=Path, default=Path(DEFAULT_BASE_DIR))
    parser.add_argument("--pipeline", type=Path,
                        default=Path("config/pipeline_statistics.toml"),
                        help="Statistics pipeline, relative to the framework root unless absolute")
    parser.add_argument("--config", type=Path, default=None,
                        help="Optional local overrides, relative to the framework root unless absolute")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print every intended command; do not launch jobs or write a summary")
    return parser


def parse_args(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    # argparse does not check choices on environment-provided defaults.
    if args.stat_mode not in ("mc", "projected"):
        parser.error("STAT_MODE/--stat-mode must be mc or projected")
    return args


def luminosity_text(value):
    # Short round-trip representation; 10.0 becomes the familiar L10 tag.
    text = str(float(value))
    if text.endswith(".0"):
        text = text[:-2]
    return text


def run_id(args, beam, scenario):
    lumi_tag = luminosity_text(args.projection_lumi).replace("+", "")
    value = f"femc-{beam}-{scenario}-{args.stat_mode}-L{lumi_tag}-{args.batch_id}"
    if not SAFE_ID.fullmatch(value):
        raise ValueError("Generated run ID contains unsafe characters")
    return value


def input_path(args, beam, scenario):
    return (args.base_dir.expanduser().resolve() / beam / "hiAcc" /
            SCENARIO_DIRS[scenario] / "merged_Upsilon_final.root")


def pipeline_path(args, project_root):
    path = args.pipeline.expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def override_path(args, project_root):
    if args.config is None:
        return None
    path = args.config.expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def build_command(args, beam, scenario, project_root):
    project_root = Path(project_root).resolve()
    command = [sys.executable, str(project_root / "workflow/run_statistics_continue.py"),
               "--pipeline", str(pipeline_path(args, project_root))]
    config = override_path(args, project_root)
    if config is not None:
        command.extend(["--config", str(config)])
    return command + [
            "--set", f"inputs.root_{beam}={input_path(args, beam, scenario)}",
            "--set", f"analysis.stat_mode={args.stat_mode}",
            "--set", f"analysis.projection_lumi_invfb={luminosity_text(args.projection_lumi)}",
            "--set", "analysis.count_interval=median_garwood",
            "run", f"note_{beam}", "--run-id", run_id(args, beam, scenario),
            "--keep-going", "--hash-inputs"]


def write_summary(path, summary, exclusive=False):
    summary["updated_at"] = utc_now().isoformat()
    text = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if exclusive:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x", encoding="utf-8") as handle:
            handle.write(text)
    else:
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(path)


def run_batch(args, project_root=None, run_command=None):
    project_root = (Path(project_root) if project_root is not None else
                    Path(__file__).resolve().parents[1]).resolve()
    runner = subprocess.run if run_command is None else run_command
    beams = BEAMS if args.beam == "all" else (args.beam,)
    scenarios = tuple(SCENARIO_DIRS) if args.scenario == "all" else (args.scenario,)
    records = [{"beam": beam, "scenario": scenario,
                "input": str(input_path(args, beam, scenario)),
                "run_id": run_id(args, beam, scenario), "status": "pending",
                "run_directory": str(project_root / "results/generated/runs" / run_id(args, beam, scenario)),
                "command": build_command(args, beam, scenario, project_root)}
               for beam in beams for scenario in scenarios]
    if args.dry_run:
        for record in records:
            print(shlex.join(record["command"]), flush=True)
        return 0
    entrypoint = project_root / "workflow/run_statistics_continue.py"
    if not entrypoint.is_file():
        raise ValueError(f"Missing continuation entry point: {entrypoint}")
    if not pipeline_path(args, project_root).is_file():
        raise ValueError(f"Missing pipeline: {pipeline_path(args, project_root)}")
    config = override_path(args, project_root)
    if config is not None and not config.is_file():
        raise ValueError(f"Missing local overrides: {config}")
    path = project_root / "results/generated/femc_batches" / (args.batch_id + ".json")
    summary = {"schema_version": 1, "batch_id": args.batch_id,
               "started_at": utc_now().isoformat(), "status": "running",
               "stat_mode": args.stat_mode, "projection_lumi_invfb": args.projection_lumi,
               "pipeline": str(pipeline_path(args, project_root)),
               "config": str(config) if config is not None else None,
               "count_interval": "median_garwood", "runs": records,
               "normalization": "Inherited from the pipeline; no cross-section or source/tuple-luminosity overrides"}
    try:
        write_summary(path, summary, exclusive=True)
    except FileExistsError:
        raise ValueError(f"Batch summary already exists: {path}. Choose a fresh --batch-id")
    print(f"Batch summary: {path}", flush=True)
    active = None
    try:
        for record in records:
            active = record
            input_error = "Input ROOT file is missing or inaccessible"
            try:
                available = Path(record["input"]).is_file()
            except OSError as error:
                available = False
                input_error = str(error)
            if not available:
                record.update(status="missing_input", error=input_error)
                print(f"SKIP {record['beam']} {record['scenario']}: {record['input']}", flush=True)
                write_summary(path, summary)
                continue
            record.update(status="running", started_at=utc_now().isoformat())
            write_summary(path, summary)
            print(shlex.join(record["command"]), flush=True)
            try:
                result = runner(record["command"], cwd=str(project_root), check=False)
                record["returncode"] = int(result.returncode)
                if result.returncode in (130, -signal.SIGINT):
                    raise KeyboardInterrupt
                record["status"] = "completed" if result.returncode == 0 else "failed"
            except OSError as error:
                record.update(status="launch_failed", error=str(error))
            record["finished_at"] = utc_now().isoformat()
            write_summary(path, summary)
    except KeyboardInterrupt:
        if active is not None and active["status"] == "running":
            active.update(status="interrupted", finished_at=utc_now().isoformat())
        summary.update(status="interrupted", finished_at=utc_now().isoformat(), exit_code=130)
        write_summary(path, summary)
        raise
    failed = sum(record["status"] != "completed" for record in records)
    summary.update(status="completed" if failed == 0 else "partial_failure",
                   completed_runs=len(records) - failed, failed_or_missing_runs=failed,
                   finished_at=utc_now().isoformat(), exit_code=0 if failed == 0 else 1)
    write_summary(path, summary)
    print(f"Finished: {len(records) - failed}/{len(records)} runs completed; summary: {path}", flush=True)
    return summary["exit_code"]


def main(argv=None):
    args = parse_args(argv)
    try:
        return run_batch(args)
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("Interrupted; remaining scenarios were not started.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
