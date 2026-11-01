#!/usr/bin/env python3
"""Execute selected plain-Python code cells from a legacy notebook.

This is a narrow compatibility bridge for the report-producing notebooks. It
does not emulate IPython magics. Paths are replaced in memory, and execution
occurs in a caller-supplied run directory so notebook defaults cannot create a
literal ``~`` directory inside the project.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys


def parse_replacement(raw: str) -> tuple[str, str]:
    if "=" not in raw:
        raise argparse.ArgumentTypeError("replacement must be OLD=NEW")
    old, new = raw.split("=", 1)
    return old, os.path.expandvars(os.path.expanduser(new))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebook")
    parser.add_argument("--cells", required=True, help="Comma-separated zero-based code-cell indices")
    parser.add_argument("--replace", action="append", default=[], type=parse_replacement)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()

    notebook_path = Path(args.notebook).expanduser().resolve()
    workdir = Path(args.workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    if "~" in workdir.parts:
        raise SystemExit(f"Refusing unresolved literal '~' output path: {workdir}")

    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    code_cells = [cell for cell in notebook.get("cells", []) if cell.get("cell_type") == "code"]
    selected = [int(item.strip()) for item in args.cells.split(",") if item.strip()]
    namespace = {"__name__": "__main__", "__file__": str(notebook_path)}
    previous_argv = sys.argv
    previous_cwd = Path.cwd()
    try:
        sys.argv = [str(notebook_path)]
        os.chdir(workdir)
        for index in selected:
            if index < 0 or index >= len(code_cells):
                raise IndexError(f"Notebook has {len(code_cells)} code cells, requested {index}")
            source = "".join(code_cells[index].get("source", []))
            for old, new in args.replace:
                source = source.replace(old, new)
            magic = [line for line in source.splitlines() if re.match(r"\s*(?:!|%[^%])", line)]
            if magic:
                raise RuntimeError(f"Cell {index} contains unsupported IPython magic: {magic[0]!r}")
            print(f"[notebook] executing {notebook_path.name} code cell {index}")
            exec(compile(source, f"{notebook_path}#code-cell-{index}", "exec"), namespace, namespace)
    finally:
        os.chdir(previous_cwd)
        sys.argv = previous_argv
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
