#!/usr/bin/env python3
"""Summarize a rematching batch's Condor processes and per-input worker receipts."""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import subprocess

from jobs import expected_output_path
from rematch_condor import default_batch, recorded_clusters, queue_rows


FILE_STATES = ("DONE", "FAILED", "started/partial", "untouched", "unreadable")
JOB_STATES = {1: "idle", 2: "running", 3: "removing", 4: "completed",
              5: "held", 6: "transferring", 7: "suspended"}


def file_state(directory, output=None, external_files=None):
    """Inspect one directory; leave ROOT I/O and checksum validation to retry and merge."""
    try:
        names = set(os.listdir(directory))
        if output is None:
            output = directory / "output.root"
        if external_files is not None and "output.root" in names:
            return "unreadable", 0, "unexpected output.root in the job control directory"
        has_tuple = output.name in (names if external_files is None else external_files)
        if "result.json" in names:
            result = json.loads((directory / "result.json").read_text())
            if not isinstance(result, dict):
                return "unreadable", 0, "malformed result.json"
            if result.get("status") == "FAILED":
                return "FAILED", 0, str(result.get("error", "worker reported failure"))
            rows = result.get("output_rows")
            if (result.get("status") == "DONE" and type(rows) is int and rows > 0
                    and result.get("input_rows") == rows and has_tuple
                    and result.get("job") == str(directory / "job.json")
                    and result.get("output") == str(output)):
                return "DONE", rows, ""
            return "unreadable", 0, f"unexpected receipt, row count or missing tuple {output}"
        if "job.json" not in names:
            return "unreadable", 0, "missing job.json"
        evidence = sorted(names.intersection({"execution.lock", "root.log", "reuse-pending.json"}))
        if has_tuple:
            evidence.append(output.name)
        if evidence:
            return "started/partial", 0, ", ".join(evidence)
        return "untouched", 0, ""
    except (OSError, ValueError) as error:
        return "unreadable", 0, str(error)


def scan_files(plan, limit, show_untouched=False):
    counts = Counter()
    total_rows = 0
    examples = {"FAILED": [], "unreadable": [], "started/partial": []}
    if show_untouched:
        examples["untouched"] = []
    # One EOS listing per report instead of opening or statting each tuple.
    external_files = set(os.listdir(plan["tuple_output_dir"])) if "tuple_output_dir" in plan else None
    for index, job in enumerate(plan["jobs"]):
        output = expected_output_path(plan, index)
        state, rows, reason = file_state(Path(job).parent, output, external_files)
        counts[state] += 1
        total_rows += rows
        if state in examples and len(examples[state]) < limit:
            examples[state].append((index, reason))
    return counts, total_rows, examples


def reference_rows(plan):
    reference = Path(__file__).resolve().parents[1] / "inputs" / "paper_inputs.csv"
    state = "Upsilon" + plan["state"]
    directory = "output_newID" if plan["state"] == "1S" else "output"
    name = f"{state}/{plan['beam']}/hiAcc/{directory}/merged_Upsilon_final.root"
    with reference.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["relative_path"] == name:
                return int(row["generated_rows"])
    return None


def latest_group_logs(batch, width, group, attempts):
    stem = f"group_{group:06d}" if width > 1 else f"{group:06d}"
    for attempt in attempts:
        prefix = attempt / "logs" / stem
        if any(prefix.with_suffix(ext).is_file() for ext in (".err", ".log", ".out")):
            return prefix
    return batch / "logs" / stem


def show_queue(batch, limit):
    try:
        clusters = recorded_clusters(batch)
        if clusters:
            print("Recorded cluster IDs:", ", ".join(map(str, clusters)))
        rows = queue_rows(batch, clusters)
        counts = Counter(JOB_STATES.get(int(row["JobStatus"]), "unknown") for row in rows)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"Condor queue unavailable: {error}")
        return None
    if counts:
        print("Current Condor processes:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
        held = (row for row in rows if int(row["JobStatus"]) == 5)
        for _, row in zip(range(limit), held):
            print(f"  held {row['ClusterId']}.{row['ProcId']}: {row.get('HoldReason') or 'reason unavailable'}")
    else:
        print("Current Condor processes: 0")
    return len(rows)


def report(batch, limit, query_queue=True, queue_only=False):
    batch = batch.resolve()
    plan = json.loads((batch / "plan.json").read_text())
    count = len(plan["jobs"])
    width = plan.get("files_per_job", 1)
    if (plan.get("work_dir") != str(batch) or not count or count != plan.get("job_count")
            or count != len(plan.get("inputs", [])) or type(width) is not int or width < 1):
        raise ValueError(f"Invalid plan: {batch / 'plan.json'}")
    print(f"{plan['state']} {plan['beam']}: {count:,} input files, up to {width} per Condor process")
    if "tuple_output_dir" in plan:
        print("Tuple directory:", plan["tuple_output_dir"])
    live = show_queue(batch, limit) if query_queue else None
    if queue_only:
        return
    counts, rows, examples = scan_files(plan, limit, show_untouched=(live == 0))
    done = counts["DONE"]
    print(f"Input files: {done:,}/{count:,} DONE ({100 * done / count:.1f}%), "
          + ", ".join(f"{name}={counts[name]:,}" for name in FILE_STATES[1:]))
    expected = reference_rows(plan)
    print(f"Reported generated rows: {rows:,}" +
          (f" (full-sample reference {expected:,})" if expected is not None else ""))
    if done == count:
        if expected is not None and rows != expected:
            print("All planned inputs have DONE receipts. Check the row total against the intended inventory.")
        else:
            print("All planned inputs have DONE receipts. The merge verifies the outputs and their checksums.")
    elif counts["started/partial"]:
        print("Started/partial includes inputs currently running and inputs interrupted before a final receipt.")
    if live == 0 and counts["untouched"]:
        print("An empty queue with untouched inputs can indicate that a launcher exited before starting the worker.")

    attempts = sorted((batch / "retries").glob("attempt-*"), reverse=True)
    for state, items in examples.items():
        for index, reason in items:
            directory = Path(plan["jobs"][index]).parent
            prefix = latest_group_logs(batch, width, index // width, attempts)
            print(f"  {state} input {index:06d}: {plan['inputs'][index]['path']}")
            if reason:
                print("    " + " ".join(reason.split())[:300])
            if (directory / "root.log").is_file():
                print("    ROOT:", directory / "root.log")
            print("    Condor stderr:", prefix.with_suffix(".err"))
            print("    Condor event log:", prefix.with_suffix(".log"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, default=default_batch(), help="Directory containing plan.json")
    parser.add_argument("--limit", type=int, default=2, help="Maximum examples per problem category (default: 2)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--queue-only", action="store_true", help="Show scheduler status without scanning per-file directories")
    mode.add_argument("--no-queue", action="store_true", help="Read worker receipts without querying Condor")
    args = parser.parse_args()
    if args.limit < 0:
        parser.error("--limit must be nonnegative")
    try:
        report(args.batch, args.limit, not args.no_queue, args.queue_only)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.exit(2, f"Cannot read batch status: {error}\n")


if __name__ == "__main__":
    main()
