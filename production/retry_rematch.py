#!/usr/bin/env python3
"""Retry unfinished files in one drained grouped rematching campaign."""
import argparse
from collections import Counter
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

from jobs import DEFAULT_CHECK_WORKERS, PRODUCERS, REMATCH, add_check_workers, ordered_checks, submit_text
from rematch_condor import queue_rows, recorded_clusters
from reuse_removed_outputs import check_done, checked_job, read_plan


IMAGE = "/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable"
WORKER_EVIDENCE = {"root.log", "execution.lock", "result.json"}


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


def scan(plan, check_workers=DEFAULT_CHECK_WORKERS):
    pending, archive, counts = [], {}, Counter()
    count = len(plan["jobs"])
    if count >= 1000:
        print(f"Checking {count:,} input records and completed tuple checksums "
              f"with up to {check_workers} concurrent checks.", flush=True)
    external_files = set(os.listdir(plan["tuple_output_dir"])) if "tuple_output_dir" in plan else None

    def inspect(index):
        path = Path(plan["jobs"][index])
        folder = path.parent
        if folder.resolve() != folder:
            raise ValueError(f"Unexpected job directory: {folder}")
        names = set(os.listdir(folder))
        if "reuse-pending.json" in names:
            raise ValueError(f"Unfinished output import at {folder}; finish removal recovery first.")
        job = checked_job(plan, index)
        output = Path(job["output"])
        evidence = [folder / name for name in sorted(names & WORKER_EVIDENCE)]
        if external_files is not None and "output.root" in names:
            raise ValueError(f"Unexpected workspace tuple: {folder / 'output.root'}")
        if output.name in (names if external_files is None else external_files):
            evidence.append(output)
        if any(path.is_symlink() or not path.is_file() for path in evidence):
            raise ValueError(f"Unexpected worker artifact in {folder}.")
        receipt = None
        if "result.json" in names:
            try:
                receipt = json.loads(path.with_name("result.json").read_bytes())
            except (json.JSONDecodeError, UnicodeDecodeError):
                pass   # Truncated evidence is archived intact before rerunning the input.
        if isinstance(receipt, dict) and receipt.get("status") == "DONE":
            check_done(job, path, receipt)
            return "DONE", [folder / "execution.lock"] if "execution.lock" in names else []
        if isinstance(receipt, dict) and receipt.get("status") not in (None, "FAILED"):
            raise ValueError(f"Unknown status in {folder}/result.json")
        if isinstance(receipt, dict) and receipt.get("status") == "FAILED":
            state = "FAILED"
        else:
            state = "partial" if evidence else "untouched"
        return state, evidence

    for index, (state, evidence) in enumerate(ordered_checks(inspect, range(count), check_workers)):
        filename = plan["jobs"][index]
        counts[state] += 1
        if state != "DONE":
            pending.append(filename)
        if evidence:
            archive[Path(filename).parent] = evidence
        if (index + 1) % 1000 == 0:
            print(f"Checked {index + 1:,}/{count:,} inputs.", flush=True)
    return pending, archive, counts


def new_attempt(batch):
    base = batch / "retries"
    number = 1 + max((int(p.name[8:]) for p in base.glob("attempt-*") if p.name[8:].isdigit()), default=0)
    attempt = base / f"attempt-{number:04d}"
    attempt.mkdir()
    (attempt / "logs").mkdir()
    return attempt


def retry(batch, submit, check_workers=DEFAULT_CHECK_WORKERS):
    plan = read_plan(batch)
    if plan["state"] not in ("1S", "2S", "3S") or plan["beam"] not in PRODUCERS:
        raise ValueError("Unsupported Upsilon state or beam in the retry plan.")
    if (plan["electron_id"], plan.get("accounting_group"), plan["environment_label"]) != (
            REMATCH, "group_u_CMS.CAF.COMM", "eic_xl:26.07.1-stable"):
        raise ValueError("Expected the pinned rematching image and the CMS CAF COMM account.")
    width = plan.get("files_per_job")
    if type(width) is not int or width < 1 or plan.get("condor_job_count") != (plan["job_count"] + width - 1) // width:
        raise ValueError("Expected a grouped plan with consistent Condor and input counts.")
    clusters = require_drained(batch)
    pending, archive, counts = scan(plan, check_workers)
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
    if any(path.parent != folder for folder, paths in archive.items() for path in paths):
        # Reserve the EOS archive before moving evidence on either filesystem.
        base = Path(plan["tuple_output_dir"]) / ".retries"
        if base.is_symlink():
            raise ValueError(f"Unexpected symlink: {base}")
        base.mkdir(exist_ok=True)
        tuple_archive = base / attempt.name
        tuple_archive.mkdir()
        (attempt / "tuple-archive.txt").write_text(str(tuple_archive) + "\n")
    for folder, paths in archive.items():
        for path in paths:
            if path.parent == folder:
                destination = attempt / "files" / folder.name
                destination.mkdir(parents=True, exist_ok=True)
            else:
                destination = tuple_archive
            path.rename(destination / path.name)

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


def main(argv=None, *, default=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, default=default, required=default is None,
                        help="Grouped batch directory containing plan.json")
    parser.add_argument("--submit", action="store_true", help="Archive failed/partial evidence and submit outstanding groups")
    add_check_workers(parser)
    args = parser.parse_args(argv)
    if args.check_workers < 1:
        parser.error("--check-workers must be positive")
    batch = args.batch.resolve()
    if not args.submit:
        retry(batch, False, args.check_workers)
        return
    # Check the scheduler before creating an archive or persistent lock file.
    require_drained(batch)
    root = batch / "retries"
    root.mkdir(exist_ok=True)
    with (root / ".retry.lock").open("a") as lock:
        fcntl.lockf(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        retry(batch, True, args.check_workers)


def cli(default=None):
    try:
        main(default=default)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
        sys.exit(f"ERROR: {error}")


if __name__ == "__main__":
    cli()
