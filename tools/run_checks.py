#!/usr/bin/env python3
"""Run each recovered test suite in its own process and working directory."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               MPLBACKEND="Agg", YPSILON_TEST_FRAMEWORK=str(ROOT / "analysis"))
    suites = [
        ("repository hygiene", ROOT, [sys.executable, "-B", "tools/check_repository.py"]),
        ("repository interface", ROOT, [sys.executable, "-B", "-m", "pytest", "tests", "-q"]),
        ("analysis backend", ROOT / "analysis", [sys.executable, "-B", "-m", "pytest", "tests", "-q"]),
        ("plotting", ROOT / "plotting", [sys.executable, "-B", "-m", "pytest", "tests", "-q"]),
        ("efficiency", ROOT / "extensions/efficiency", [sys.executable, "-B", "-m", "pytest", "tests", "-q"]),
        ("9x130 extension", ROOT / "extensions/reoptimization_9x130", [sys.executable, "-B", "-m", "pytest", "tests", "-q"]),
    ]
    failures = []
    for label, directory, command in suites:
        print("\nCHECK " + label, flush=True)
        result = subprocess.run(command, cwd=directory, env=env)
        if result.returncode:
            failures.append(label)
    if failures:
        print("\nFAILED: " + ", ".join(failures), flush=True)
        return 1
    print("\nPASS all repository checks", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
