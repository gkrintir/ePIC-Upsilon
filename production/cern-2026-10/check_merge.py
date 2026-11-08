"""Compare a local CERN production archive with the reference 2S 9x130 tuple."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from zipfile import ZipFile

import numpy as np
import uproot


CAMPAIGN = Path(__file__).with_name("campaign.json")


def log_column(lines, prefix):
    return [line[len(prefix):] for line in lines if line.startswith(prefix)]


def logged_counts(archive, campaign):
    jobs = campaign["jobs"]
    list_bytes = archive.read("eic/list.txt")
    if hashlib.sha256(list_bytes).hexdigest() != campaign["input_list_sha256"]:
        raise ValueError("The archive's input list differs from the reference campaign")
    inputs = list_bytes.decode("ascii").splitlines()
    cluster_logs = [name for name in archive.namelist()
                    if re.fullmatch(r"eic/logs/upsilon_\d+\.log", name)]
    if len(cluster_logs) != 1:
        raise ValueError("Expected exactly one tuple-production Condor event log")
    event_log = archive.read(cluster_logs[0]).decode("utf-8")
    stem = cluster_logs[0][:-len(".log")]
    cluster_id = stem.rsplit("_", 1)[1]
    completed = re.findall(r"^005 \(" + cluster_id + r"\.(\d+)\.000\) [^\n]+\n"
                           r"\s*\(1\) Normal termination \(return value (\d+)\)", event_log, re.M)
    if sorted((int(proc), int(status)) for proc, status in completed) != [
            (proc, 0) for proc in range(jobs["processes"])]:
        raise ValueError("The event log does not confirm all processes completed normally")

    listed, outputs, counts = [], [], []
    for proc in range(jobs["processes"]):
        shard = archive.read(f"eic/list_part_{proc:03d}").decode("ascii").splitlines()
        lines = archive.read(f"{stem}_{proc}.out").decode("utf-8").splitlines()
        written = log_column(lines, "Flat Tree written to: ")
        expected_names = ["AnalysisOut_" + PurePosixPath(p).stem.replace(".", "_") + ".root" for p in shard]
        if (log_column(lines, "Processing: ") != shard or log_column(lines, "Output: ") != written
                or [PurePosixPath(p).name for p in written] != expected_names):
            raise ValueError(f"Input and completed output records disagree for process {proc}")
        sizes = log_column(lines, "fEMC clusters masked: ")
        mode = "full (mode 0), eta split = 1.68"
        electron_id = "reco (mode 1); 0.8 <= E/p <= 1.2; track-linked isolation >= 0.9 in R=0.4"
        if (len(sizes) != len(shard) or any(not row.endswith(f"; affected events: 0 / {jobs['events_per_logged_input']}")
                                            for row in sizes)
                or log_column(lines, "fEMC scenario: ") != [mode] * len(shard)
                or log_column(lines, "Electron ID: ") != [electron_id] * len(shard)):
            raise ValueError(f"Event counts or detector settings disagree for process {proc}")
        matches = log_column(lines, "Reco-ID match counts: ")
        if len(matches) != len(shard):
            raise ValueError(f"Missing per-input electron-ID counts for process {proc}")
        for row in matches:
            match = re.fullmatch(r"usable=(\d+), missing=(\d+), masked-only=(\d+), ambiguous=(\d+)", row)
            if not match:
                raise ValueError(f"Unrecognized electron-ID counts for process {proc}")
            counts.append(tuple(map(int, match.groups())))
        listed.extend(shard)
        outputs.extend(expected_names)
    if (listed != inputs or len(set(inputs)) != jobs["unique_reconstructed_inputs"]
            or len(set(outputs)) != len(inputs)):
        raise ValueError("The process lists do not exactly cover the reference list")
    return np.asarray(counts, dtype=np.int64)[np.argsort(outputs)]


def check(archive_path, merged_root):
    campaign = json.loads(CAMPAIGN.read_text())
    reference = campaign["merged_tuple_check"]
    with ZipFile(archive_path) as archive:
        logged = logged_counts(archive, campaign)
    with merged_root.open("rb") as stream:
        digest = hashlib.sha256()
        for chunk in iter(lambda: stream.read(4 * 1024 ** 2), b""):
            digest.update(chunk)
        if digest.hexdigest() != reference["sha256"]:
            raise ValueError("Merged file differs from the archived 2S 9x130 ROOT checksum")
        stream.seek(0)
        with uproot.open(stream) as root:
            tree = root["UpsilonTree"]
            if tree.num_entries != len(logged) * reference["block_size"]:
                raise ValueError("Unexpected number of event rows")
            block_sums = []
            for field in reference["signature_order"]:
                values = tree[field].array(library="np")
                block_sums.append(values.reshape(-1, reference["block_size"]).sum(axis=1))
    matched = int(np.count_nonzero(np.all(np.column_stack(block_sums) == logged, axis=1)))
    print(f"{matched}/{len(logged)} blocks match the logged counters in output filename order.")
    if matched != len(logged):
        raise ValueError("Electron-ID block comparison failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Local ZIP containing the original submission lists and logs")
    parser.add_argument("merged_root", type=Path, help="Reference merged 2S 9x130 ROOT file")
    options = parser.parse_args()
    check(options.archive, options.merged_root)
