#!/usr/bin/env python3
"""Carry validated file outputs from a removed Condor batch into a fresh grouped plan."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

from jobs import root_command, sha256, submit_text


WORK_SETTINGS = ("producer_identity", "beam", "state", "data_version", "electron_id",
                 "femc", "eta_split", "source_hashes", "normalization",
                 "environment_label", "setup_script")


def read_json(path):
    return json.loads(path.read_bytes())


def replace_atomically(path, write):
    fd, filename = tempfile.mkstemp(dir=path.parent, prefix=".reuse-")
    temporary = Path(filename)
    try:
        with os.fdopen(fd, "wb") as stream:
            write(stream)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_plan(root):
    plan = read_json(root / "plan.json")
    count = len(plan["inputs"])
    if (plan["schema_version"] != 1 or plan["kind"] != "new_shared_filesystem_condor_plan"
            or plan["work_dir"] != str(root) or not count or count != plan["job_count"]
            or len(plan["jobs"]) != count
            or len({record["path"] for record in plan["inputs"]}) != count):
        raise ValueError(f"Invalid batch manifest: {root}")
    if sha256(root / "snapshot/jobs.py") != plan["worker_sha256"]:
        raise ValueError(f"Frozen worker checksum mismatch: {root}")
    for name, digest in plan["source_hashes"].items():
        if Path(name).name != name or sha256(root / "snapshot/macros" / name) != digest:
            raise ValueError(f"Frozen source mismatch: {root}, {name}")
    for i, p in enumerate(plan["jobs"]):
        if p != str(root / "jobs" / f"job_{i:06d}" / "job.json"):
            raise ValueError(f"Unexpected job path at index {i}: {p}")
    return plan


def checked_job(plan, index):
    path = Path(plan["jobs"][index])
    job = read_json(path)
    folder = path.parent
    if (folder.resolve() != folder or path.is_symlink() or job["job_id"] != index
            or job["job_dir"] != str(folder) or job["output"] != str(folder / "output.root")
            or job["macro_dir"] != str(Path(plan["work_dir"]) / "snapshot/macros")
            or job["inputs"] != [plan["inputs"][index]]
            or any(job.get(k) != v for k, v in plan.items() if k not in ("jobs", "inputs"))):
        raise ValueError(f"Job differs from the batch plan: {path}")
    return job


def carry_file(old, new, index):
    old_path, new_path = (Path(p["jobs"][index]) for p in (old, new))
    old_dir, new_dir = old_path.parent, new_path.parent
    old_names, new_names = set(os.listdir(old_dir)), set(os.listdir(new_dir))
    new_result = new_dir / "result.json"
    pending = new_dir / "reuse-pending.json"
    new_output = new_dir / "output.root"
    if new_names.intersection({"execution.lock", "root.log"}):
        raise ValueError(f"Grouped worker already started: {new_dir}")

    source = old_dir / "result.json"
    receipt = None
    if "execution.lock" not in old_names and "result.json" in old_names:
        if source.is_symlink():
            raise ValueError(f"Unexpected symlink: {source}")
        source_bytes = source.read_bytes()
        status = json.loads(source_bytes)
        if status["status"] == "DONE":
            receipt = status
        elif status["status"] != "FAILED":
            raise ValueError(f"Unrecognized worker status: {source}")
    if receipt is None:
        if new_names.intersection({"result.json", "output.root", "reuse-pending.json"}):
            raise ValueError(f"New attempt has evidence but no reusable source: {new_dir}")
        return False

    old_job = checked_job(old, index)
    checked_job(new, index)
    original = old_dir / "output.root"
    rows = receipt["input_rows"]
    if (receipt.get("job") != str(old_path) or receipt.get("output") != str(original)
            or type(rows) is not int or rows <= 0 or type(receipt["output_rows"]) is not int
            or rows != receipt["output_rows"] or receipt.get("root_returncode") != 0
            or receipt.get("command") != root_command(old_job)
            or any(receipt.get(k) != old_job[k] for k in
                   ("state", "beam", "electron_id", "source_hashes", "normalization", "environment_label"))
            or original.is_symlink() or not original.is_file() or original.stat().st_size == 0
            or sha256(original) != receipt["output_sha256"]):
        raise ValueError(f"DONE receipt or output failed validation: {source}")
    record = old["inputs"][index]
    raw = Path(record["path"])
    st = raw.stat()
    input_hash = record.get("sha256") or new["inputs"][index].get("sha256")
    if ((st.st_size, st.st_mtime_ns) != (record["size_bytes"], record["mtime_ns"])
            or (input_hash and sha256(raw) != input_hash)):
        raise ValueError(f"Original input changed: {raw}")

    imported = dict(receipt, job=str(new_path), output=str(new_output),
                    reused_from={"job": str(old_path), "output": str(original),
                                 "result_sha256": hashlib.sha256(source_bytes).hexdigest()})
    if any(path.is_symlink() for path in (new_result, new_output, pending)):
        raise ValueError(f"Unexpected symlink in grouped output: {new_dir}")
    if "result.json" in new_names:
        if read_json(new_result) != imported or sha256(new_output) != receipt["output_sha256"]:
            raise ValueError(f"Existing grouped result is not this verified import: {new_dir}")
        return True
    if "reuse-pending.json" in new_names:
        if read_json(pending) != imported:
            raise ValueError(f"Unrecognized import in progress: {pending}")
    else:
        if "output.root" in new_names:
            raise ValueError(f"Refusing to replace an existing grouped output: {new_output}")
        replace_atomically(pending, lambda f: f.write((json.dumps(imported, indent=2) + "\n").encode()))
    if not new_output.exists():
        def copy_output(f):
            with original.open("rb") as src:
                shutil.copyfileobj(src, f, length=1024 * 1024)
        replace_atomically(new_output, copy_output)
    if sha256(new_output) != receipt["output_sha256"]:
        raise ValueError(f"Copied tuple checksum mismatch: {new_output}")
    os.replace(pending, new_result)
    return True


def reuse(old_root, new_root):
    old_root, new_root = old_root.resolve(), new_root.resolve()
    if old_root == new_root or old_root in new_root.parents or new_root in old_root.parents:
        raise ValueError("Source and destination must be separate batch directories")
    old, new = read_plan(old_root), read_plan(new_root)
    if any(old[k] != new[k] for k in WORK_SETTINGS):
        raise ValueError("Incompatible source and destination production settings")
    if len(old["inputs"]) != len(new["inputs"]):
        raise ValueError("Different source and destination input counts")
    width = new.get("files_per_job")
    if (type(width) is not int or width < 1
            or new.get("condor_job_count") != (new["job_count"] + width - 1) // width):
        raise ValueError("Invalid destination grouping")
    for a, b in zip(old["inputs"], new["inputs"]):
        if (any(a[k] != b[k] for k in ("path", "size_bytes", "mtime_ns")) or
                (a.get("sha256") and b.get("sha256") and a["sha256"] != b["sha256"])):
            raise ValueError(f"Input inventories differ: {a['path']} vs {b['path']}")

    with (new_root / ".reuse.lock").open("a") as lock:
        fcntl.lockf(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        remaining = [path for i, path in enumerate(new["jobs"]) if not carry_file(old, new, i)]
        text = submit_text(new, remaining) if remaining else ""
        replace_atomically(new_root / "after-removal.sub", lambda f: f.write(text.encode()))
    total = len(new["jobs"])
    print(f"Reused {total - len(remaining)} verified outputs; {len(remaining)} of {total} inputs to run.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old_batch", type=Path)
    parser.add_argument("new_batch", type=Path)
    args = parser.parse_args()
    try:
        reuse(args.old_batch, args.new_batch)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.exit(2, f"ERROR: {error}\n")
