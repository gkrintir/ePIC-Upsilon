#!/usr/bin/env python3
"""Retry unfinished files in a drained 1S 9x275 grouped Condor campaign."""
import argparse
from collections import Counter
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from jobs import REMATCH, root_command, sha256, submit_text
from rematch_condor import default_batch, queue_rows, recorded_clusters
from reuse_removed_outputs import checked_job, read_plan


IMAGE = "/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable"
WORKER_EVIDENCE = {"root.log", "output.root", "execution.lock", "result.json"}


def require_drained(batch):
    if (batch.parent / "recovery.lock").exists():
        raise ValueError("Removal recovery is in progress. Finish it before retrying this batch.")
    clusters = recorded_clusters(batch)
    if not clusters:
        raise ValueError("No confirmed grouped submission. Use the original grouped launcher first.")
    rows = queue_rows(batch, clusters)
    if rows:
        ids = ", ".join(str(i) for i in sorted({r["ClusterId"] for r in rows}))
        raise ValueError(f"{len(rows)} Condor processes are still in the queue (cluster {ids}). "
                         "Check the status/hold reason before retrying.")
    return clusters


def check_done(job, job_path, receipt):
    """The producing worker already checked ROOT counts; verify that its bytes are intact."""
    rows = receipt.get("output_rows")
    output = job_path.with_name("output.root")
    expected = dict(job=str(job_path), output=str(output), root_returncode=0)
    for key in ("state", "beam", "electron_id", "source_hashes", "environment_label", "normalization"):
        expected[key] = job[key]
    execution = dict(job)
    if "reused_from" in receipt:
        origin = receipt["reused_from"]
        if not isinstance(origin, dict) or not re.fullmatch(r"[0-9a-f]{64}", str(origin.get("result_sha256"))):
            raise ValueError(f"Invalid imported receipt provenance: {job_path}")
        old_job = Path(origin.get("job", ""))
        if (not old_job.is_absolute() or old_job.name != "job.json" or old_job.parent.name != job_path.parent.name
                or old_job.parent.parent.name != "jobs" or str(old_job.with_name("output.root")) != origin.get("output")):
            raise ValueError(f"Invalid imported receipt paths: {job_path}")
        execution.update(output=origin["output"], macro_dir=str(old_job.parents[2] / "snapshot/macros"))
    expected["command"] = root_command(execution)
    if (type(rows) is not int or rows <= 0 or type(receipt.get("input_rows")) is not int
            or receipt["input_rows"] != rows or any(receipt.get(k) != v for k, v in expected.items())
            or not output.is_file() or sha256(output) != receipt.get("output_sha256")):
        raise ValueError(f"DONE receipt or tuple failed validation; inspect {job_path.parent}. No automatic retry.")


def scan(plan):
    pending, archive, counts = [], {}, Counter()
    for index, filename in enumerate(plan["jobs"]):
        path = Path(filename)
        folder = path.parent
        if folder.resolve() != folder:
            raise ValueError(f"Unexpected job directory: {folder}")
        names = set(os.listdir(folder))
        if "reuse-pending.json" in names:
            raise ValueError(f"Unfinished output import at {folder}; finish removal recovery first.")
        evidence = sorted(names & WORKER_EVIDENCE)
        if any(not (folder / n).is_file() or (folder / n).is_symlink() for n in evidence):
            raise ValueError(f"Unexpected worker artifact in {folder}.")
        job = checked_job(plan, index)
        receipt = None
        if "result.json" in evidence:
            try:
                receipt = json.loads(path.with_name("result.json").read_bytes())
            except (json.JSONDecodeError, UnicodeDecodeError):
                pass   # Truncated evidence is archived intact before rerunning the input.
        if isinstance(receipt, dict) and receipt.get("status") == "DONE":
            check_done(job, path, receipt)
            counts["DONE"] += 1
            if "execution.lock" in names:
                archive[folder] = ["execution.lock"]
            continue
        if isinstance(receipt, dict) and receipt.get("status") not in (None, "FAILED"):
            raise ValueError(f"Unknown status in {folder}/result.json")
        if isinstance(receipt, dict) and receipt.get("status") == "FAILED":
            state = "FAILED"
        else:
            state = "partial" if evidence else "untouched"
        counts[state] += 1
        pending.append(filename)
        if evidence:
            archive[folder] = evidence
    return pending, archive, counts


def new_attempt(batch):
    base = batch / "retries"
    number = 1 + max((int(p.name[8:]) for p in base.glob("attempt-*") if p.name[8:].isdigit()), default=0)
    attempt = base / f"attempt-{number:04d}"
    attempt.mkdir()
    (attempt / "logs").mkdir()
    return attempt


def retry(batch, submit):
    plan = read_plan(batch)
    if (plan["state"], plan["beam"], plan["electron_id"], plan.get("accounting_group"),
            plan["environment_label"]) != ("1S", "9x275", REMATCH, "group_u_CMS.CAF.COMM", "eic_xl:26.07.1-stable"):
        raise ValueError("Expected a pinned 1S 9x275 rematching plan with the CMS CAF COMM account.")
    width = plan.get("files_per_job")
    if type(width) is not int or width < 1 or plan.get("condor_job_count") != (plan["job_count"] + width - 1) // width:
        raise ValueError("Expected a grouped plan with consistent Condor and input counts.")
    clusters = require_drained(batch)
    pending, archive, counts = scan(plan)
    groups = {int(Path(path).parent.name.removeprefix("job_")) // width for path in pending}
    print(f"Verified {counts['DONE']} DONE. Unfinished inputs: {counts['FAILED']} failed, "
          f"{counts['partial']} partial, {counts['untouched']} untouched.")
    print(f"Retry: {len(pending)} files across {len(groups)} Condor groups.", flush=True)
    if not pending or not submit:
        return

    require_drained(batch)
    attempt = new_attempt(batch)
    print(f"Retry record: {attempt}", flush=True)
    (attempt / "retry.sub").write_text(submit_text(plan, pending, log_dir=attempt / "logs"))
    (attempt / "file_jobs.txt").write_text("\n".join(pending) + "\n")
    (attempt / "prior-clusters.json").write_text(json.dumps(clusters) + "\n")
    for folder, names in archive.items():
        destination = attempt / "files" / folder.name
        destination.mkdir(parents=True)
        for name in names:
            (folder / name).rename(destination / name)

    require_drained(batch)
    command = ["condor_submit", "-maxjobs", str(len(groups)), "-append", "getenv = False", "-append",
               f'+SingularityImage = "{IMAGE}"', str(attempt / "retry.sub")]
    # The reservation precedes scheduler contact. An interrupted/ambiguous submit blocks any repetition.
    with (attempt / "submission.txt").open("x") as record:
        process = subprocess.run(command, stdout=record, stderr=subprocess.STDOUT)
    response = (attempt / "submission.txt").read_text()
    print(response, end="", flush=True)
    if process.returncode:
        raise RuntimeError(f"condor_submit exited with {process.returncode}. Reconcile the scheduler and {record.name} before retrying.")
    recorded_clusters(batch)  # require an identifiable cluster ID, also after an apparently successful call


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, default=default_batch(), help="Grouped batch directory containing plan.json")
    parser.add_argument("--submit", action="store_true", help="Archive failed/partial evidence and submit outstanding groups")
    args = parser.parse_args(argv)
    batch = args.batch.resolve()
    if not args.submit:
        retry(batch, False)
        return
    # Check the scheduler before creating an archive or persistent lock file.
    require_drained(batch)
    root = batch / "retries"
    root.mkdir(exist_ok=True)
    with (root / ".retry.lock").open("a") as lock:
        fcntl.lockf(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        retry(batch, True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
        sys.exit(f"ERROR: {error}")
