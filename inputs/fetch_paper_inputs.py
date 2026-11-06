#!/usr/bin/env python3
"""Fetch or verify the six ROOT files for the Upsilon paper and two-energy comparison."""

import argparse
import csv
import hashlib
import http.client
import os
from pathlib import Path
import sys
from tempfile import NamedTemporaryFile
import time
from urllib.parse import quote
from urllib.request import Request, urlopen


MANIFEST = Path(__file__).with_name("paper_inputs.csv")
DAV_ROOT = "https://cernbox.cern.ch/remote.php/dav/public-files/UHPvC04xj29cxrp/"
CHUNK_SIZE = 4 * 1024 * 1024


def verify(path, size, checksum):
    if not path.is_file() or path.stat().st_size != size:
        raise ValueError(f"expected a regular file of {size} bytes")
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != checksum:
        raise ValueError("SHA-256 mismatch")


def read_chunk(url, start, end, size):
    expected_range = f"bytes {start}-{end}/{size}"
    length = end - start + 1
    request = Request(
        f"{url}?part={start}-{end}",
        headers={"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"},
    )
    for attempt in range(4):
        try:
            with urlopen(request, timeout=90) as response:
                if response.status != 206 or response.headers.get("Content-Range") != expected_range:
                    raise OSError(f"server did not return {expected_range}")
                chunk = response.read(length + 1)
                if len(chunk) != length:
                    raise OSError(f"received {len(chunk)} bytes for {expected_range}")
                return chunk
        except (OSError, http.client.HTTPException) as error:
            if attempt == 3:
                raise OSError(f"failed to read {expected_range}: {error}") from error
            time.sleep(2**attempt)


def fetch(relative_path, path, size, checksum):
    path.parent.mkdir(parents=True, exist_ok=True)
    url = DAV_ROOT + quote(relative_path, safe="/")
    digest = hashlib.sha256()
    with NamedTemporaryFile(dir=path.parent, prefix=".paper-input-", suffix=".part") as output:
        for start in range(0, size, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, size) - 1
            chunk = read_chunk(url, start, end, size)
            output.write(chunk)
            digest.update(chunk)
        output.flush()
        if os.fstat(output.fileno()).st_size != size or digest.hexdigest() != checksum:
            raise ValueError("download failed size or SHA-256 verification")
        try:
            # Link the completed temporary inode into place only if the target is free.
            os.link(output.name, path)
        except FileExistsError:
            verify(path, size, checksum)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dest", metavar="DEST", type=Path, help="base directory for the three Upsilon state trees")
    parser.add_argument("--check-only", action="store_true", help="verify every file locally without downloading")
    args = parser.parse_args()
    failed = False
    with MANIFEST.open(newline="", encoding="utf-8") as manifest:
        for row in csv.DictReader(manifest):
            relative_path = row["relative_path"]
            path = args.dest / relative_path
            size = int(row["size_bytes"])
            checksum = row["sha256"]
            try:
                if path.exists() or path.is_symlink():
                    verify(path, size, checksum)
                elif args.check_only:
                    raise FileNotFoundError("missing file")
                else:
                    print(f"Downloading {relative_path}", flush=True)
                    fetch(relative_path, path, size, checksum)
                print(f"Verified {relative_path}", flush=True)
            except (OSError, ValueError) as error:
                print(f"{path}: {error}", file=sys.stderr)
                failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
