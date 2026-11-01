#!/usr/bin/env python3
"""Incremental BNL Upsilon downloads using native, parallel XRootD copy jobs.

Python 3.9+ and xrdfs/xrdcp are required. Run this file; do not source it.
Remote listings (including missing directories) expire after six hours. Use
--refresh to discover new production files immediately. The first local scan
stats every existing ROOT file once. Subsequent runs list names and reuse sizes
from the inventory. Use --rescan-local after editing/replacing files externally.
The fast inventory assumes immutable files and one writer per dataset. Separate
datasets may run concurrently; two runs of the same dataset are blocked,
including dry runs. A compatibility lock also blocks older versions that lock
the entire destination: let those finish or restart them with this version.
These locks coordinate this program, not unrelated download scripts.

The cap counts logical bytes in immediate .root files plus this program's staged
copies, not filesystem allocation, metadata, logs, or analysis subdirectories.
It assumes remote sizes remain unchanged between listing and copying. Final
sizes are checked; use --checksum adler32 for additional transfer verification.
Existing files with matching size are not thereby proven intact. Imported files
have unknown provenance; --repair-existing replaces ambiguous duplicates,
size mismatches, and files with known outdated source provenance.

Examples:
  python3 download_upsilon_fast.py --dataset Upsilon1S/9x275/hiAcc --dry-run
  python3 download_upsilon_fast.py --jobs 12
  python3 download_upsilon_fast.py --refresh --repair-existing --dry-run
  python3 download_upsilon_fast.py --refresh --repair-existing

Failed batches are discarded without publishing any of their files. A killed
process can leave a clearly named staging directory; the next run stops and
reports it rather than ignoring its bytes. No existing final file is deleted
before its replacement has transferred successfully.
"""

import argparse
import concurrent.futures
from contextlib import contextmanager
import errno
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import signal
import socket
import stat
import subprocess
import sys
import tempfile
import threading
import time

SERVER = "root://epicxrd1.sdcc.bnl.gov:1095"
PREFIX = "/eic/EPIC/RECO/26.07.1/epic_craterlake"
SUFFIX = "EXCLUSIVE/UPSILON_ABCONV/eSTARlight1.3.0-1.0"
STATES = ("Upsilon1S", "Upsilon2S", "Upsilon3S")
ENERGIES = ("9x275", "9x130")
CONFIGS = ("hiAcc", "hiDiv")
GIB = 1024 ** 3
STAGE_PREFIX = ".upsilon-transfer-part-"


def say(message):
    print(message, flush=True)


@contextmanager
def transfer_lock(path, mode, busy_message, owner=None):
    """Keep the lock inode in place; close the handle to release the lock."""
    with path.open("a+") as handle:
        try:
            fcntl.flock(handle, mode | fcntl.LOCK_NB)
        except OSError as error:
            detail = (f"{errno.errorcode.get(error.errno, 'UNKNOWN')} "
                      f"(errno {error.errno}): {error.strerror}")
            if error.errno in (errno.EAGAIN, errno.EWOULDBLOCK):
                recorded = ""
                if owner is not None:
                    try:
                        handle.seek(0)
                        holder = json.load(handle)
                        if isinstance(holder, dict):
                            recorded = (f" Recorded holder: host={holder.get('hostname', '?')}, "
                                        f"PID={holder.get('pid', '?')}, "
                                        f"started={holder.get('started', '?')}.")
                    except (OSError, ValueError):
                        pass
                raise RuntimeError(f"{busy_message}{recorded} Lock: {path}. {detail}") from error
            raise RuntimeError(f"Cannot acquire filesystem lock {path}: {detail}. "
                               "Check locking support and access to this destination.") from error
        if owner is not None:
            # Do not replace/unlink a lock file: that would change its inode.
            handle.seek(0)
            handle.truncate()
            json.dump(owner, handle)
            handle.write("\n")
            handle.flush()
        yield


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as out:
            json.dump(value, out, separators=(",", ":"))
            out.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read_json(path):
    with path.open() as handle:
        return json.load(handle)


def stop_child(child):
    if child.poll() is not None:
        return
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        child.wait(timeout=10)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL)
        child.wait()


def run_command(command, timeout=None, control=None):
    """Stop our child before cleaning any of its destinations on interruption."""
    child = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True,
                             start_new_session=True)
    try:
        if control is not None:
            with control["lock"]:
                if control["cancelled"]:
                    raise RuntimeError("Remote listing cancelled.")
                control["children"].add(child)
        stdout, stderr = child.communicate(timeout=timeout)
        return child.returncode, stdout, stderr
    except BaseException:
        stop_child(child)
        child.communicate()
        raise
    finally:
        if control is not None:
            with control["lock"]:
                control["children"].discard(child)


def parse_listing(text, directory):
    """Parse compact and owner/group variants of xrdfs ls -l.

    Compact: permissions date time bytes /path
    Extended: permissions owner group bytes date time /path
    """
    files = {}
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) == 5:
            mode, date, clock, size, path = fields
        elif len(fields) == 7:
            mode, _owner, _group, size, date, clock, path = fields
        else:
            mode = date = clock = size = path = ""
        if (not re.fullmatch(r"[-dlrwx]+", mode)
                or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date)
                or not re.fullmatch(r"\d{2}:\d{2}:\d{2}", clock)
                or not size.isdigit() or not path.startswith("/")):
            excerpt = line if len(line) <= 240 else line[:240] + "..."
            raise RuntimeError("Unrecognized xrdfs ls -l format at line "
                               f"{number} ({len(fields)} fields): {excerpt!r}")
        if mode.startswith("d") or not path.endswith(".root"):
            continue
        if not mode.startswith("-"):
            raise RuntimeError(f"Nonregular remote ROOT entry: {path}")
        if os.path.dirname(path).rstrip("/") != directory.rstrip("/"):
            raise RuntimeError(f"Unexpected path in remote listing: {path}")
        name = os.path.basename(path)
        record = {"name": name, "size": int(size), "source": SERVER + "/" + path,
                  "stamp": date + " " + clock}
        if record["size"] <= 0:
            raise RuntimeError(f"Remote ROOT file has zero bytes: {path}")
        if name in files and files[name] != record:
            raise RuntimeError(f"Conflicting entries in one listing: {path}")
        files[name] = record
    return list(files.values())


def remote_listing(directory, cache_dir, args, control=None):
    start = time.monotonic()
    key = hashlib.sha256((SERVER + directory).encode()).hexdigest()
    path = cache_dir / (key + ".json")
    if path.exists() and not args.refresh:
        cached = read_json(path)
        if (cached.get("server") == SERVER and cached.get("directory") == directory
                and 0 <= time.time() - cached["created"] < args.cache_hours * 3600):
            return cached["files"], "cache", time.monotonic() - start
    code, out, err = run_command(["xrdfs", SERVER, "ls", "-l", directory],
                                timeout=args.list_timeout, control=control)
    if code:
        if "[3011]" in err or "no such file or directory" in err.lower():
            files = []
            status = "missing"
        else:
            raise RuntimeError(f"xrdfs failed ({code}) for {directory}\n{err[-3000:]}")
    else:
        files = parse_listing(out, directory)
        status = "fresh"
    atomic_json(path, {"server": SERVER, "directory": directory,
                       "created": time.time(), "files": files})
    return files, status, time.monotonic() - start


def fetch_sources(dataset, cache_dir, args):
    state, energy, config = dataset.split("/")
    tail = f"{SUFFIX}/{state}/{energy}/q2_0to0.01/{config}"
    directories = [f"{PREFIX}/try3/{tail}", f"{PREFIX}/try2/{tail}", f"{PREFIX}/{tail}"]
    results = {}
    say("  Remote listings: try3, try2, base (three concurrent queries/cache reads)")
    control = {"lock": threading.Lock(), "children": set(), "cancelled": False}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        pending = {pool.submit(remote_listing, d, cache_dir, args, control): i
                   for i, d in enumerate(directories)}
        try:
            for future in concurrent.futures.as_completed(pending):
                index = pending[future]
                rows, status, elapsed = future.result()
                results[index] = rows
                say(f"    {('try3', 'try2', 'base')[index]}: {len(rows):,} files, "
                    f"{status}, {elapsed:.2f} s")
        except BaseException:
            # Cancel inside the context, before its shutdown waits for workers.
            with control["lock"]:
                control["cancelled"] = True
                children = list(control["children"])
            for child in children:
                stop_child(child)
            for future in pending:
                future.cancel()
            raise
    chosen, locations = {}, {}
    for index in range(3):
        for row in sorted(results[index], key=lambda value: value["name"]):
            name = row["name"]
            locations[name] = locations.get(name, 0) + 1
            chosen.setdefault(name, row)
    return chosen, {name for name, count in locations.items() if count > 1}


def local_inventory(directory, index_path, force):
    start = time.monotonic()
    previous = read_json(index_path).get("files", {}) if index_path.exists() else {}
    records, stated = {}, 0
    names = os.listdir(directory)
    leftovers = [name for name in names if name.startswith(STAGE_PREFIX)]
    if leftovers:
        raise RuntimeError("Unfinished staging directories require inspection after "
                           "confirming no old transfer is running: " + ", ".join(leftovers)
                           + ". Remove the abandoned staging copies before rerunning.")
    for name in names:
        if not name.endswith(".root"):
            continue
        old = previous.get(name)
        if old and not force:
            records[name] = old
            continue
        info = (directory / name).lstat()
        stated += 1
        if not stat.S_ISREG(info.st_mode):
            raise RuntimeError(f"Local .root entry is not a regular file: {directory / name}")
        unchanged = (old and old["size"] == info.st_size
                     and old.get("mtime_ns") == info.st_mtime_ns)
        records[name] = old if unchanged else {
            "size": info.st_size, "mtime_ns": info.st_mtime_ns,
            "source": None, "stamp": None}
    total = sum(record["size"] for record in records.values())
    say(f"  Local inventory: {len(records):,} files, {total / GIB:.3f} GiB, "
        f"{stated:,} stat calls, {time.monotonic() - start:.2f} s")
    atomic_json(index_path, {"files": records})
    return records, total


def make_plan(remote, duplicates, local, total, args):
    started = time.monotonic()
    repairs, missing, suspect = [], [], []
    for name, row in remote.items():
        old = local.get(name)
        if old is None:
            missing.append(dict(row, old_size=0))
            continue
        bad_size = old["size"] != row["size"]
        old_source = old.get("source")
        outdated = old_source and (old_source != row["source"]
                                   or old.get("stamp") != row["stamp"])
        ambiguous = name in duplicates and not old_source
        if bad_size or outdated or ambiguous:
            suspect.append(name)
            if args.repair_existing:
                repairs.append(dict(row, old_size=old["size"]))
    batches, batch = [], []
    current, staged, after, peak, blocked, selected = total, 0, total, total, 0, 0
    for job in repairs + missing:
        if args.max_files is not None and selected >= args.max_files:
            break
        if batch and (len(batch) >= args.batch_files or current + staged + job["size"] > args.cap_bytes):
            batches.append(batch)
            batch, staged, current = [], 0, after
        if current + staged + job["size"] > args.cap_bytes:
            blocked += 1
            continue
        batch.append(job)
        selected += 1
        staged += job["size"]
        after += job["size"] - job["old_size"]
        peak = max(peak, current + staged)
    if batch:
        batches.append(batch)
    say(f"  Compare/plan: {time.monotonic() - started:.3f} s; "
        f"{len(remote):,} unique remote names, {len(duplicates):,} duplicated across tries")
    say(f"  Missing: {len(missing):,}; existing files needing review/repair: {len(suspect):,}")
    if suspect and not args.repair_existing:
        say("  WARNING: existing suspect files are kept. Inspect with --repair-existing "
            "--dry-run, then use --repair-existing to replace them from the latest try.")
    say(f"  Selected: {selected:,} files in {len(batches):,} batches; "
        f"final {after / GIB:.3f} GiB; peak including staging {peak / GIB:.3f} GiB; "
        f"cap {args.cap_bytes / GIB:.3f} GiB")
    if blocked:
        say(f"  {blocked:,} candidates do not fit the cap, including temporary replacement space.")
    if total > args.cap_bytes:
        say("  WARNING: existing data already exceed the cap; no files are removed.")
    planned_names = {job["name"] for group in batches for job in group}
    deferred = [name for name in suspect if name not in planned_names]
    if deferred:
        reason = "cap or pilot limit" if args.repair_existing else "--repair-existing not enabled"
        say(f"  Repairs deferred: {len(deferred):,} ({reason}). These existing files remain unresolved.")
    return batches, suspect, deferred


def copy_batch(batch, directory, state_dir, local, args, batch_number):
    stage = Path(tempfile.mkdtemp(prefix=STAGE_PREFIX, dir=directory))
    start = time.monotonic()
    try:
        inputs = state_dir / "inputs.txt"
        inputs.write_text("".join(job["source"] + "\n" for job in batch))
        command = ["xrdcp", "--nopbar", "--parallel", str(args.jobs),
                   "--infiles", str(inputs)]
        if args.streams:
            command.extend(["--streams", str(args.streams)])
        if args.checksum:
            command.extend(["--cksum", args.checksum])
        command.append(str(stage) + "/")
        size = sum(job["size"] for job in batch)
        say(f"  Batch {batch_number}: {len(batch)} files, {size / GIB:.3f} GiB, "
            f"{args.jobs} parallel copies")
        code, out, err = run_command(command)
        if code:
            (state_dir / "last-xrdcp-error.log").write_text(out + "\n" + err)
            raise RuntimeError(f"xrdcp returned {code}. No files from this batch were "
                               f"published. Log: {state_dir / 'last-xrdcp-error.log'}\n{err[-2000:]}")
        prepared = {}
        for job in batch:
            name = job["name"]
            info = (stage / name).lstat()
            if not stat.S_ISREG(info.st_mode) or info.st_size != job["size"]:
                raise RuntimeError(f"Transferred size/type mismatch for {name}; refresh remote listings.")
            if name in local:
                old = (directory / name).lstat()
                if (old.st_size != local[name]["size"]
                        or old.st_mtime_ns != local[name]["mtime_ns"]):
                    raise RuntimeError(f"Existing file changed while copying: {name}. Rescan required.")
            elif os.path.lexists(directory / name):
                raise RuntimeError(f"Another writer created destination: {name}")
            prepared[name] = {"size": info.st_size, "mtime_ns": info.st_mtime_ns,
                              "source": job["source"], "stamp": job["stamp"]}
        # An interrupted publication must trigger an actual local rescan next time.
        marker = state_dir / "rescan-required"
        marker.write_text("Interrupted publication requires a fresh local inventory.\n")
        for job in batch:
            name = job["name"]
            os.replace(stage / name, directory / name)
            local[name] = prepared[name]
        atomic_json(state_dir / "local.json", {"files": local})
        marker.unlink()
        elapsed = time.monotonic() - start
        say(f"    Completed in {elapsed:.2f} s, "
            f"{size / max(elapsed, 1e-9) / 1024 ** 2:.2f} MiB/s including verification/publication")
    finally:
        shutil.rmtree(stage)


def parser():
    result = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    result.add_argument("--dest", type=Path, default=Path("./upsilon_data"))
    result.add_argument("--dataset", action="append", help="STATE/ENERGY/CONFIG; repeatable; default all 12")
    result.add_argument("--jobs", type=int, default=12)
    result.add_argument("--streams", type=int, default=0, help="data streams per file; 0 uses xrdcp default")
    result.add_argument("--batch-files", type=int, default=256)
    result.add_argument("--cap-gib", type=float, default=400)
    result.add_argument("--cap-bytes", type=int, help="override --cap-gib with an integer byte limit")
    result.add_argument("--max-files", type=int, help="limit selected transfers per dataset for a pilot")
    result.add_argument("--cache-hours", type=float, default=6)
    result.add_argument("--refresh", action="store_true", help="refresh all requested remote listings")
    result.add_argument("--rescan-local", action="store_true", help="stat all existing files and update the inventory")
    result.add_argument("--repair-existing", action="store_true", help="replace suspect/outdated files from latest try")
    result.add_argument("--dry-run", action="store_true", help="read/update metadata and print plans; do not copy ROOT files")
    result.add_argument("--checksum", help="pass checksum type to xrdcp, e.g. adler32 (extra remote metadata work)")
    result.add_argument("--list-timeout", type=float, default=300)
    return result


def main(argv=None):
    command_parser = parser()
    args = command_parser.parse_args(argv)
    if args.cap_bytes is None:
        args.cap_bytes = int(args.cap_gib * GIB)
    if (args.jobs < 1 or args.batch_files < 1 or args.cap_bytes < 0
            or not 0 <= args.streams <= 15 or args.cache_hours < 0
            or args.list_timeout <= 0 or (args.max_files is not None and args.max_files < 1)):
        command_parser.error("Invalid numeric option.")
    all_datasets = ["/".join(parts) for parts in itertools.product(STATES, ENERGIES, CONFIGS)]
    datasets = list(dict.fromkeys(args.dataset or all_datasets))
    if any(item not in all_datasets for item in datasets):
        command_parser.error("Dataset must be STATE/ENERGY/CONFIG from the configured production.")
    for executable in ("xrdfs",) if args.dry_run else ("xrdfs", "xrdcp"):
        if not shutil.which(executable):
            command_parser.error(f"{executable} is not in PATH; initialize your XRootD environment.")
    if not args.dry_run:
        _, out, err = run_command(["xrdcp", "--help"], timeout=30)
        if not all(flag in out + err for flag in ("--infiles", "--parallel")):
            command_parser.error("This xrdcp lacks --infiles/--parallel; use a compatible client.")
    args.dest = args.dest.expanduser().resolve()
    args.dest.mkdir(parents=True, exist_ok=True)
    # Shared in this version, exclusive in older versions. Keep this guard for
    # the entire run so old and new versions cannot overlap unsafely.
    with transfer_lock(
            args.dest / ".upsilon-transfer.lock", fcntl.LOCK_SH,
            "An older downloader holds the destination-wide exclusive lock. "
            "Let that run finish, or stop it cleanly and restart with this version."):
        root_state = args.dest / ".upsilon-transfer"
        failures, reviews = 0, 0
        for dataset in datasets:
            say(f"\nDataset {dataset}")
            try:
                directory = args.dest / dataset
                state_dir = root_state / dataset
                state_dir.mkdir(parents=True, exist_ok=True)
                owner = {"pid": os.getpid(), "hostname": socket.gethostname(),
                         "dataset": dataset,
                         "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
                with transfer_lock(
                        state_dir / "dataset.lock", fcntl.LOCK_EX,
                        f"Dataset {dataset} is already in use by another run "
                        "(including a dry run). Let it finish or stop it cleanly.", owner):
                    directory.mkdir(parents=True, exist_ok=True)
                    marker = state_dir / "rescan-required"
                    local, total = local_inventory(directory, state_dir / "local.json",
                                                    args.rescan_local or marker.exists())
                    if marker.exists():
                        marker.unlink()
                    if total >= args.cap_bytes and not args.repair_existing and not args.refresh:
                        say("  At/above cap; skipping remote queries. Use --repair-existing --dry-run "
                            "to audit existing files. Replacements require temporary headroom.")
                        continue
                    remote, duplicates = fetch_sources(dataset, root_state / "remote", args)
                    batches, suspect, deferred = make_plan(remote, duplicates, local, total, args)
                    plan_path = state_dir / "plan.json"
                    atomic_json(plan_path, {"dataset": dataset, "cap_bytes": args.cap_bytes,
                                           "selected": [job for batch in batches for job in batch],
                                           "existing_needing_review": suspect,
                                           "repairs_deferred": deferred})
                    say(f"  Plan details: {plan_path}")
                    reviews += bool(deferred)
                    if not args.dry_run:
                        for number, batch in enumerate(batches, 1):
                            copy_batch(batch, directory, state_dir, local, args, number)
            except (OSError, RuntimeError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
                failures += 1
                say(f"  ERROR: {error}")
        say(f"\n{'Dry run' if args.dry_run else 'Run'} finished: "
            f"{failures} dataset errors; {reviews} datasets have existing files awaiting review.")
        return 1 if failures else (2 if reviews else 0)


if __name__ == "__main__":
    # Turn SIGTERM into an exception so the copy child is stopped before cleanup.
    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        say("\nInterrupted. Completed files remain. Interrupted publication is rescanned next run.")
        sys.exit(130)
    except (RuntimeError, OSError) as error:
        say(f"ERROR: {error}")
        sys.exit(1)
