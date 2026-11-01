#!/usr/bin/env python3
"""Check the existing framework and install section rendering when needed."""
import ast
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = (
    "workflow/run_analysis.py", "workflow/check_statistical_run.py",
    "config/pipeline_statistics.toml", "scripts/validation/prepare_latest_sample.py",
    "scripts/reporting/render_note_physics.py",
)


def main():
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    if missing:
        print("Install into the existing framework after its statistics/note updates. Missing:")
        print("\n".join(missing))
        return 2
    target = ROOT / "scripts/reporting/render_note_physics.py"
    tree = ast.parse(target.read_text(encoding="utf-8"))
    section = any(isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                  and node.func.attr == "add_argument"
                  and any(isinstance(arg, ast.Constant) and arg.value == "--section" for arg in node.args)
                  for node in ast.walk(tree))
    if not section:
        result = subprocess.run([sys.executable, str(ROOT / "workflow/newid_renderer_update.py"),
                                 "--project-root", str(ROOT)], check=False)
        if result.returncode:
            return result.returncode
    else:
        print("Section rendering is already available.")
    result = subprocess.run([sys.executable, str(ROOT / "workflow/run_newid.py"), "list"], check=False)
    if result.returncode:
        return result.returncode
    print("Electron-ID and half-plane workflows ready. Use a fresh run ID after this update.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
