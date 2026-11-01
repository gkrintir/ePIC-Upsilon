#!/usr/bin/env python3
"""Install only the four renderer blocks needed for independent note tasks."""
import argparse
import ast
import datetime as dt
import json
import os
from pathlib import Path
import stat
import tempfile


def updated_source(source, replacements):
    tree = ast.parse(source)
    lines = source.splitlines(keepends=True)
    changes = []
    for replacement in replacements:
        nodes = tree.body
        scope = replacement["scope"]
        if scope != "module":
            parents = [node for node in nodes if isinstance(node, ast.ClassDef) and node.name == scope]
            if len(parents) != 1:
                raise ValueError("Unrecognized renderer class: " + scope)
            nodes = parents[0].body
        matches = [node for node in nodes if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                   and node.name == replacement["name"]]
        if len(matches) != 1:
            raise ValueError("Unrecognized renderer method: " + replacement["name"])
        node = matches[0]
        actual = "".join(lines[node.lineno - 1:node.end_lineno])
        # Compare without the final newline, which is outside the AST node.
        installed_blocks = [replacement["new"], *replacement.get("installed_variants", [])]
        if any(actual.rstrip("\r\n") == block.rstrip("\r\n") for block in installed_blocks):
            continue
        if actual.rstrip("\r\n") != replacement["old"].rstrip("\r\n"):
            raise ValueError("Renderer block was customized: " + scope + "." + node.name
                             + ". No changes were written. Retain this file for a targeted update.")
        ending = "\n" if actual.endswith("\n") else ""
        changes.append((node.lineno - 1, node.end_lineno, replacement["new"].rstrip("\r\n") + ending))
    for start, end, replacement in sorted(changes, reverse=True):
        lines[start:end] = [replacement]
    result = "".join(lines)
    ast.parse(result, feature_version=(3, 9))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--check", action="store_true", help="Check applicability without writing")
    args = parser.parse_args()
    target = args.project_root / "scripts/reporting/render_note_physics.py"
    patch_path = Path(__file__).with_name("continuation_renderer_patch.json")
    try:
        original = target.read_bytes()
        source = original.decode("utf-8")
        patch = json.loads(patch_path.read_text(encoding="utf-8"))
        updated = updated_source(source, patch["replacements"])
        if updated == source:
            print("Continuation renderer already installed.")
            return 0
        if args.check:
            print("Renderer update is applicable. No files changed.")
            return 0
        timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        backup = target.with_name(target.name + ".before_continuation_" + timestamp + ".bak")
        with backup.open("xb") as handle:
            handle.write(original)
        mode = stat.S_IMODE(target.stat().st_mode)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".continuation_", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(updated.encode("utf-8"))
                handle.flush()
                os.fsync(handle.fileno())
            temporary.chmod(mode)
            os.replace(str(temporary), str(target))
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()
        print("Updated " + str(target))
        print("Backup: " + str(backup))
        print("Use a new run ID: the statistical cache contract includes source-code hashes.")
        return 0
    except (OSError, ValueError, SyntaxError) as error:
        parser.exit(2, str(error) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
