#!/usr/bin/env python3
"""Configuration-driven runner for the Upsilon analysis stages.

The runner deliberately imports no analysis dependencies. Each ROOT-dependent
stage executes in its own subprocess, with logs and manifests isolated in its
run directory. A task may explicitly select another working directory.
"""

from __future__ import annotations

import argparse
import datetime as dt
import glob
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
from typing import Any


if __package__:
    from .normalization_contract import validate_flat_normalization, check_normalization_cache
else:
    # Support direct execution and importlib loading from an absolute file path.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from normalization_contract import validate_flat_normalization, check_normalization_cache


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PIPELINE = PROJECT_ROOT / "config" / "pipeline.toml"
PLACEHOLDER = re.compile(r"\{([A-Za-z0-9_.-]+)\}")


class ConfigurationError(RuntimeError):
    pass


def read_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def set_nested(config: dict[str, Any], dotted_key: str, value: Any) -> None:
    parts = dotted_key.split(".")
    cursor = config
    for part in parts[:-1]:
        child = cursor.setdefault(part, {})
        if not isinstance(child, dict):
            raise ConfigurationError(f"Cannot set {dotted_key}: {part} is not a table")
        cursor = child
    cursor[parts[-1]] = value


def parse_override(raw: str) -> tuple[str, Any]:
    if "=" not in raw:
        raise ConfigurationError(f"--set expects KEY=VALUE, got {raw!r}")
    key, value = raw.split("=", 1)
    try:
        parsed = tomllib.loads(f"value = {value}")["value"]
    except tomllib.TOMLDecodeError:
        parsed = value
    return key, parsed


def flatten(config: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in config.items():
        dotted = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            result.update(flatten(value, dotted))
        else:
            result[dotted] = value
    return result


def resolve_template(value: Any, values: dict[str, Any]) -> Any:
    if not isinstance(value, str):
        return value

    missing: set[str] = set()

    def replace(match: re.Match[str]) -> str:
        key = match.group(1)
        if key not in values:
            missing.add(key)
            return match.group(0)
        resolved = values[key]
        return "" if resolved is None else str(resolved)

    resolved = PLACEHOLDER.sub(replace, value)
    if missing:
        raise ConfigurationError(f"Unresolved placeholders in {value!r}: {sorted(missing)}")
    return os.path.expandvars(os.path.expanduser(resolved))


def resolve_path(value: str, project_root: Path = PROJECT_ROOT) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    if not path.is_absolute():
        path = project_root / path
    return path.resolve(strict=False)


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: Path, include_hash: bool = False) -> dict[str, Any]:
    stat = path.stat()
    result: dict[str, Any] = {
        "path": str(path),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if include_hash:
        result["sha256"] = sha256(path)
    return result


def load_configuration(args: argparse.Namespace, run_id: str) -> dict[str, Any]:
    pipeline_path = resolve_path(args.pipeline)
    config = read_toml(pipeline_path)
    if args.config:
        config = deep_merge(config, read_toml(resolve_path(args.config)))
    for raw in getattr(args, "set_values", []) or []:
        key, value = parse_override(raw)
        set_nested(config, key, value)

    run_dir = PROJECT_ROOT / "results" / "generated" / "runs" / run_id
    config.setdefault("paths", {})
    config["paths"]["project"] = str(PROJECT_ROOT)
    config["paths"]["reference"] = str(PROJECT_ROOT / "results" / "reference")
    config["paths"]["generated"] = str(PROJECT_ROOT / "results" / "generated")
    config.setdefault("run", {})
    config["run"]["id"] = run_id
    config["run"]["dir"] = str(run_dir)
    check_category_configuration(config)
    try:
        validate_flat_normalization(flatten(config))
    except ValueError as error:
        raise ConfigurationError(str(error)) from error
    return config


def task_index(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    tasks = config.get("task", [])
    if not isinstance(tasks, list):
        raise ConfigurationError("pipeline.toml must contain [[task]] entries")
    index: dict[str, dict[str, Any]] = {}
    for task in tasks:
        name = task.get("name")
        if not name or name in index:
            raise ConfigurationError(f"Task names must be nonempty and unique: {name!r}")
        index[name] = task
    return index


def expand_targets(config: dict[str, Any], requested: list[str]) -> list[str]:
    index = task_index(config)
    groups = config.get("groups", {})
    requested = requested or ["core"]
    expanded: list[str] = []
    for item in requested:
        if item in groups:
            values = groups[item]
            if not isinstance(values, list):
                raise ConfigurationError(f"Group {item} must be an array")
            expanded.extend(values)
        elif item in index:
            expanded.append(item)
        else:
            raise ConfigurationError(f"Unknown task or group: {item}")
    return list(dict.fromkeys(expanded))


def dependency_closure(config: dict[str, Any], targets: list[str]) -> list[str]:
    index = task_index(config)
    seen: set[str] = set()
    visiting: set[str] = set()
    ordered: list[str] = []

    def visit(name: str) -> None:
        if name in seen:
            return
        if name in visiting:
            raise ConfigurationError(f"Dependency cycle involving {name}")
        if name not in index:
            raise ConfigurationError(f"Unknown dependency {name}")
        visiting.add(name)
        for dep in index[name].get("deps", []):
            visit(dep)
        visiting.remove(name)
        seen.add(name)
        ordered.append(name)

    for target in targets:
        visit(target)
    return ordered


def resolved_task(task: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    values = flatten(config)
    resolved = dict(task)
    for key in ("script", "cwd"):
        if key in resolved:
            resolved[key] = resolve_template(resolved[key], values)
    for key in ("args", "requirements", "outputs", "fail_patterns"):
        resolved[key] = [resolve_template(item, values) for item in resolved.get(key, [])]
    return resolved


def check_category_configuration(config, tasks=None):
    """Check isolated-category paths before any runner output-directory creation."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from proton_category_isolation import enforce_category_configuration
    checked = config if tasks is None else dict(config, task=list(tasks))
    try:
        enforce_category_configuration(checked, PROJECT_ROOT, resolved_task)
    except ValueError as error:
        raise ConfigurationError(str(error)) from error


def requirement_result(requirement: str) -> tuple[bool, str]:
    kind, separator, payload = requirement.partition(":")
    if not separator:
        return False, f"malformed requirement {requirement!r}"
    payload = payload.strip()
    if kind == "file":
        if not payload:
            return False, "required file path is blank"
        path = resolve_path(payload)
        return path.is_file(), f"file {path}"
    if kind == "directory":
        if not payload:
            return False, "required directory path is blank"
        path = resolve_path(payload)
        return path.is_dir(), f"directory {path}"
    if kind == "module":
        return importlib.util.find_spec(payload) is not None, f"Python module {payload}"
    if kind in {"command", "executable"}:
        if not payload:
            return False, "required executable is blank"
        candidate = Path(payload)
        found = candidate.is_file() and os.access(candidate, os.X_OK)
        if not found:
            found = shutil.which(payload) is not None
        return found, f"executable {payload}"
    if kind == "value":
        return bool(payload), "nonempty configured value"
    return False, f"unknown requirement kind {kind!r}"


def validate_task(task: dict[str, Any], config: dict[str, Any]) -> list[str]:
    check_category_configuration(config, [task])
    errors: list[str] = []
    resolved = resolved_task(task, config)
    script = resolve_path(resolved["script"])
    if not script.is_file():
        errors.append(f"script missing: {script}")
    else:
        try:
            compile(script.read_text(encoding="utf-8"), str(script), "exec")
        except SyntaxError as exc:
            errors.append(f"script does not compile: line {exc.lineno}, {exc.msg}")

    if resolved.get("cwd"):
        working_directory = resolve_path(str(resolved["cwd"]))
        if not working_directory.is_dir():
            errors.append(f"working directory missing: {working_directory}")

    for requirement in resolved["requirements"]:
        ok, message = requirement_result(requirement)
        if not ok:
            errors.append(f"missing requirement: {message}")

    for output in resolved["outputs"]:
        if "~" in Path(output).parts:
            errors.append(f"output contains an unresolved literal '~': {output}")
        parent = Path(output.split("*", 1)[0]).parent
        # Validation must not create an uncontracted run directory.
        # Inspect the nearest existing ancestor; create parents at execution.
        probe = parent
        while not probe.exists() and probe != probe.parent:
            if probe.is_symlink():
                break
            probe = probe.parent
        if not probe.is_dir() or not os.access(probe, os.W_OK | os.X_OK):
            errors.append(f"output directory ancestor is not writable: {probe}")
    return errors


def configured_input_fingerprints(config: dict[str, Any], include_hash: bool) -> list[dict[str, Any]]:
    fingerprints: list[dict[str, Any]] = []
    for key, value in config.get("inputs", {}).items():
        if not key.endswith(("_root", "_csv", "_file", "_executable")) or not value:
            continue
        path = resolve_path(str(value))
        if path.is_file():
            item = file_fingerprint(path, include_hash=include_hash)
            item["config_key"] = f"inputs.{key}"
            fingerprints.append(item)
    return fingerprints


def verify_configured_checksums(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    checksums = config.get("checksums", {})
    for key, expected in checksums.items():
        if not expected:
            continue
        path_value = config.get("inputs", {}).get(key)
        if not path_value:
            errors.append(f"checksum configured for inputs.{key}, but its path is blank")
            continue
        path = resolve_path(str(path_value))
        if not path.is_file():
            errors.append(f"cannot verify checksum, file missing: {path}")
            continue
        observed = sha256(path)
        if observed.lower() != str(expected).lower():
            errors.append(f"checksum mismatch for inputs.{key}: expected {expected}, observed {observed}")
    return errors


def output_matches(pattern: str) -> list[Path]:
    return [Path(item) for item in glob.glob(pattern) if Path(item).is_file()]


def task_outputs_complete(task: dict[str, Any]) -> bool:
    if not task["outputs"]:
        return False
    return all(any(path.stat().st_size > 0 for path in output_matches(pattern)) for pattern in task["outputs"])


def command_for(task: dict[str, Any], config: dict[str, Any]) -> list[str]:
    python = str(config.get("executables", {}).get("python", sys.executable))
    script = str(resolve_path(task["script"]))
    return [python, script, *[str(arg) for arg in task["args"]]]


def run_task(
    task: dict[str, Any],
    config: dict[str, Any],
    *,
    force: bool,
    dry_run: bool,
    hash_inputs: bool,
) -> tuple[str, dict[str, Any]]:
    check_category_configuration(config, [task])
    resolved = resolved_task(task, config)
    command = command_for(resolved, config)
    task_dir = Path(config["run"]["dir"]) / "stages" / resolved["name"]
    if not dry_run:
        for output in resolved["outputs"]:
            parent = Path(output.split("*", 1)[0]).parent
            parent.mkdir(parents=True, exist_ok=True)
        task_dir.mkdir(parents=True, exist_ok=True)
    execution_cwd = resolve_path(str(resolved["cwd"])) if resolved.get("cwd") else task_dir

    if task_outputs_complete(resolved) and not force:
        print(f"CACHED  {resolved['name']}")
        return "CACHED", {
            "task": resolved["name"],
            "status": "CACHED",
            "command": command,
            "cwd": str(execution_cwd),
        }

    print(f"READY   {resolved['name']}")
    print("        " + shlex.join(command))
    if resolved.get("cwd"):
        print(f"        cwd: {execution_cwd}")
    if dry_run:
        return "READY", {
            "task": resolved["name"],
            "status": "READY",
            "command": command,
            "cwd": str(execution_cwd),
        }

    started = dt.datetime.now(dt.timezone.utc)
    log_path = task_dir / "stage.log"
    captured: list[str] = []
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=execution_cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log.write(line)
            captured.append(line)
        process.stdout.close()
        return_code = process.wait()
    finished = dt.datetime.now(dt.timezone.utc)
    combined = "".join(captured)

    status = "DONE"
    reasons: list[str] = []
    if return_code != 0:
        status = "FAILED"
        reasons.append(f"exit code {return_code}")
    for pattern in resolved["fail_patterns"]:
        if re.search(pattern, combined, flags=re.IGNORECASE):
            status = "INVALID_OUTPUT"
            reasons.append(f"log matched failure pattern {pattern!r}")
    if return_code == 0 and not task_outputs_complete(resolved):
        status = "INVALID_OUTPUT"
        reasons.append("one or more declared outputs are missing or empty")

    produced: list[dict[str, Any]] = []
    for pattern in resolved["outputs"]:
        for path in output_matches(pattern):
            produced.append(file_fingerprint(path, include_hash=True))

    record: dict[str, Any] = {
        "task": resolved["name"],
        "stage": resolved.get("stage", ""),
        "status": status,
        "reasons": reasons,
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "duration_seconds": (finished - started).total_seconds(),
        "command": command,
        "cwd": str(execution_cwd),
        "script": file_fingerprint(resolve_path(resolved["script"]), include_hash=True),
        "configured_inputs": configured_input_fingerprints(config, include_hash=hash_inputs),
        "outputs": produced,
        "log": str(log_path),
    }
    (task_dir / "manifest.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{status:<8}{resolved['name']}")
    return status, record


def write_run_manifest(config: dict[str, Any], records: list[dict[str, Any]]) -> Path:
    run_dir = Path(config["run"]["dir"])
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "run_id": config["run"]["id"],
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "project_root": str(PROJECT_ROOT),
        "configuration": config,
        "tasks": records,
    }
    path = run_dir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def cmd_list(args: argparse.Namespace) -> int:
    config = load_configuration(args, "listing")
    print(f"{'STAGE':<18} {'TASK':<42} {'MODE':<13} DESCRIPTION")
    for task in config.get("task", []):
        mode = "experimental" if task.get("experimental", False) else "default"
        print(f"{task.get('stage',''):<18} {task['name']:<42} {mode:<13} {task.get('description','')}")
    return 0


def cmd_plan(args: argparse.Namespace) -> int:
    config = load_configuration(args, args.run_id or "plan")
    targets = expand_targets(config, args.targets)
    order = dependency_closure(config, targets)
    index = task_index(config)
    for position, name in enumerate(order, start=1):
        task = resolved_task(index[name], config)
        print(f"{position:02d}. {name} [{task.get('stage','')}]")
        print("    " + shlex.join(command_for(task, config)))
        if task.get("cwd"):
            print(f"    cwd: {resolve_path(str(task['cwd']))}")
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    config = load_configuration(args, args.run_id or "validation")
    targets = expand_targets(config, args.targets or ["all"])
    order = dependency_closure(config, targets)
    index = task_index(config)
    failures = 0
    for checksum_error in verify_configured_checksums(config):
        print(f"INVALID  configuration: {checksum_error}")
        failures += 1
    for name in order:
        task = index[name]
        if task.get("experimental", False) and not args.include_experimental:
            print(f"SKIP     {name}: experimental")
            continue
        errors = validate_task(task, config)
        if errors:
            failures += 1
            print(f"BLOCKED  {name}")
            for error in errors:
                print(f"         {error}")
        else:
            print(f"READY    {name}")
    return 1 if failures else 0


def cmd_run(args: argparse.Namespace) -> int:
    run_id = args.run_id or dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    config = load_configuration(args, run_id)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from proton_selection import guard_proton_run
    try:
        guard_proton_run(config, PROJECT_ROOT, dry_run=args.dry_run)
    except ValueError as error:
        raise ConfigurationError(str(error)) from error
    if not args.dry_run:
        try:
            check_normalization_cache(config, PROJECT_ROOT)
        except ValueError as error:
            raise ConfigurationError(str(error)) from error
    targets = expand_targets(config, args.targets)
    order = dependency_closure(config, targets)
    index = task_index(config)

    checksum_errors = verify_configured_checksums(config)
    if checksum_errors:
        for error in checksum_errors:
            print(f"INVALID  configuration: {error}")
        return 2

    records: list[dict[str, Any]] = []
    statuses: dict[str, str] = {}
    failures = 0
    for name in order:
        task = index[name]
        if task.get("experimental", False) and not args.include_experimental:
            statuses[name] = "SKIPPED_EXPERIMENTAL"
            records.append({"task": name, "status": statuses[name]})
            print(f"SKIP     {name}: pass --include-experimental to run")
            continue
        blocked_deps = [dep for dep in task.get("deps", []) if statuses.get(dep) not in {"DONE", "CACHED", "READY"}]
        if blocked_deps:
            statuses[name] = "BLOCKED_DEPENDENCY"
            records.append({"task": name, "status": statuses[name], "blocked_by": blocked_deps})
            print(f"BLOCKED  {name}: dependencies {blocked_deps}")
            failures += 1
            if not args.keep_going:
                break
            continue
        errors = validate_task(task, config)
        if errors:
            statuses[name] = "BLOCKED_INPUT"
            records.append({"task": name, "status": statuses[name], "errors": errors})
            print(f"BLOCKED  {name}")
            for error in errors:
                print(f"         {error}")
            failures += 1
            if not args.keep_going:
                break
            continue
        status, record = run_task(
            task,
            config,
            force=args.force,
            dry_run=args.dry_run,
            hash_inputs=args.hash_inputs,
        )
        statuses[name] = status
        records.append(record)
        if status in {"FAILED", "INVALID_OUTPUT"}:
            failures += 1
            if not args.keep_going:
                break

    if args.dry_run:
        print("Dry run completed; no run artifacts or manifests were written.")
        return 1 if failures else 0
    manifest = write_run_manifest(config, records)
    print(f"Run manifest: {manifest}")
    return 1 if failures else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pipeline", default=str(DEFAULT_PIPELINE), help="Pipeline TOML file")
    parser.add_argument("--config", help="Local TOML overrides, normally config/local.toml")
    parser.add_argument("--set", dest="set_values", action="append", default=[], help="Override KEY=VALUE")
    sub = parser.add_subparsers(dest="command", required=True)

    list_parser = sub.add_parser("list", help="List every stage")
    list_parser.set_defaults(func=cmd_list)

    plan_parser = sub.add_parser("plan", help="Print the resolved DAG and commands")
    plan_parser.add_argument("targets", nargs="*", help="Task or group names")
    plan_parser.add_argument("--run-id")
    plan_parser.set_defaults(func=cmd_plan)

    validate_parser = sub.add_parser("validate", help="Check code, software, inputs, and output locations")
    validate_parser.add_argument("targets", nargs="*", help="Task or group names")
    validate_parser.add_argument("--run-id")
    validate_parser.add_argument("--include-experimental", action="store_true")
    validate_parser.set_defaults(func=cmd_validate)

    run_parser = sub.add_parser("run", help="Execute a task group in an isolated run directory")
    run_parser.add_argument("targets", nargs="*", help="Task or group names")
    run_parser.add_argument("--run-id")
    run_parser.add_argument("--force", action="store_true")
    run_parser.add_argument("--dry-run", action="store_true")
    run_parser.add_argument("--keep-going", action="store_true")
    run_parser.add_argument("--hash-inputs", action="store_true", help="SHA-256 input files in addition to stat metadata")
    run_parser.add_argument("--include-experimental", action="store_true")
    run_parser.set_defaults(func=cmd_run)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return int(args.func(args))
    except ConfigurationError as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
