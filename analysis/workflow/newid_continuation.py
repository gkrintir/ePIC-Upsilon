#!/usr/bin/env python3
"""Continue independent statistical analyses and render physics sections separately.

This opt-in entry point preserves the statistical preflight and the existing
runner. Its configuration transformation and runner hooks are scoped to main;
importing this module never changes the existing runner's functions.
"""
from contextlib import contextmanager
from copy import deepcopy
import ast
import datetime as dt
import json
from pathlib import Path
import re
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
    sys.modules["tomllib"] = tomllib

import run_analysis as runner
from check_statistical_run import build_contract, check_contract


SECTIONS = (
    "w_extraction", "differential", "momentum_shapes", "w_unfolding",
    "t_unfolding", "w_resolution", "t_resolution",
)
SECTION_DEPENDENCIES = {
    "w_extraction": "extract_{energy}_sigma_gamma_p",
    "differential": "extract_{energy}_differential_gamma_p",
    "momentum_shapes": "unfold_{energy}_t_spectrum",
    "w_unfolding": "unfold_{energy}_sigma_gamma_p_w",
    "t_unfolding": "unfold_{energy}_t_spectrum",
    "w_resolution": "validate_{energy}_w_resolution",
    "t_resolution": "unfold_{energy}_t_spectrum",
}
SECTION_STEMS = {
    "w_extraction": ("cross_section_w", "w_reconstruction_efficiency", "bin_by_bin_closure"),
    "differential": ("pt2_spectra", "rapidity_coverage", "b_slope_comparison"),
    "momentum_shapes": (
        "momentum_transfer_comparison", "impact_parameter_profiles", "momentum_transfer_shapes",
        "momentum_transfer_epic_shape_covariance", "momentum_transfer_zeus_shape_covariance",
        "momentum_transfer_epic_shape_toy_covariance",
        "impact_parameter_epic_covariance", "impact_parameter_zeus_covariance",
    ),
    "w_unfolding": ("w_unfolding_closure", "w_response_probability"),
    "t_unfolding": ("t_unfolding_closure", "t_response_probability"),
    "w_resolution": ("w_resolution", "w_relative_resolution"),
    "t_resolution": ("t_resolution",),
}
SUCCESS_STATUSES = {"DONE", "CACHED"}


def _option_value(arguments, option):
    value = None
    for i, argument in enumerate(arguments):
        if argument == option:
            if i+1 == len(arguments):
                raise runner.ConfigurationError("Missing value for " + option)
            value = arguments[i+1]
        elif isinstance(argument, str) and argument.startswith(option + "="):
            value = argument.split("=", 1)[1]
    return value


def _section_arguments(arguments, section):
    result = []
    skip = False
    for argument in arguments:
        if skip:
            skip = False
            continue
        if argument == "--section":
            skip = True
        elif not (isinstance(argument, str) and argument.startswith("--section=")):
            result.append(argument)
    return result + ["--section", section]


def _reference_requirement(requirement, reference):
    if not requirement.startswith("file:"):
        return False
    payload = requirement.partition(":")[2]
    if "{inputs.zeus_step_csv}" in payload:
        return True
    def normalized(value):
        return str(value).removeprefix("{paths.project}/").removeprefix("./")
    return reference is not None and normalized(payload) == normalized(reference)


def transform_configuration(config):
    """Return a copy with each monolithic physics task split into seven tasks.

    Existing arguments, execution settings and additional shared prerequisites
    are retained. Known physics prerequisites and output files are assigned to
    their consumer section. Unknown output ownership is an explicit error, not
    a reason to silently drop the user's declared output validation.
    """
    result = deepcopy(config)
    tasks = result.get("task", [])
    if not isinstance(tasks, list):
        raise runner.ConfigurationError("pipeline must contain [[task]] entries")
    existing_names = {task.get("name") for task in tasks}
    replacement = {}
    transformed = []
    output_owner = {stem: section for section, stems in SECTION_STEMS.items() for stem in stems}
    for original in tasks:
        name = original.get("name", "")
        match = re.fullmatch(r"render_(9x130|9x275)_note_physics", name)
        if match is None:
            transformed.append(original)
            continue
        energy = match.group(1)
        outdir = _option_value(original.get("args", []), "--out-dir")
        if not outdir:
            raise runner.ConfigurationError(name + " requires --out-dir to assign section audits")
        reference = _option_value(original.get("args", []), "--zeus-csv")
        names = [name + "_" + section for section in SECTIONS]
        if any(section_name in existing_names for section_name in names):
            raise runner.ConfigurationError(name + ": generated section task name already exists")
        replacement[name] = names
        known_deps = {pattern.format(energy=energy) for pattern in SECTION_DEPENDENCIES.values()}
        outputs = {section: [] for section in SECTIONS}
        for output in original.get("outputs", []):
            filename = Path(output).name
            if filename == "note_physics_audit.json":
                continue
            owner = output_owner.get(Path(filename).stem)
            if owner is None:
                raise runner.ConfigurationError(name + ": cannot assign output to a physics section: " + output)
            outputs[owner].append(output)
        for section, section_name in zip(SECTIONS, names):
            task = deepcopy(original)
            task["name"] = section_name
            task["description"] = original.get("description", "Render physics note") + " [" + section + "]"
            task["args"] = _section_arguments(original.get("args", []), section)
            required_dep = SECTION_DEPENDENCIES[section].format(energy=energy)
            task["deps"] = [dep for dep in original.get("deps", [])
                            if dep == required_dep or dep not in known_deps]
            if required_dep not in task["deps"]:
                if required_dep not in existing_names:
                    raise runner.ConfigurationError(section_name + ": required upstream task missing: " + required_dep)
                task["deps"].append(required_dep)
            task["requirements"] = [item for item in original.get("requirements", [])
                                    if section == "momentum_shapes" or not _reference_requirement(item, reference)]
            task["outputs"] = outputs[section] + [str(Path(outdir) / ("note_physics_" + section + "_audit.json"))]
            transformed.append(task)
    def expanded(items):
        values = [value for item in items for value in replacement.get(item, [item])]
        return list(dict.fromkeys(values))
    groups = result.setdefault("groups", {})
    for name, members in list(groups.items()):
        if not isinstance(members, list):
            raise runner.ConfigurationError("Group " + name + " must be an array")
        groups[name] = expanded(members)
    for name, names in replacement.items():
        if name in groups and groups[name] != names:
            raise runner.ConfigurationError("Cannot replace existing custom group alias " + name)
        groups[name] = names
    for task in transformed:
        task["deps"] = expanded(task.get("deps", []))
    result["task"] = transformed
    # The original runner validates dependency closure for the requested
    # targets. Do not let an unrelated, unrequested task block this transform.
    runner.task_index(result)
    return result


def prior_failure_reasons(task_name, run_dir):
    """Explain why existing output files are insufficient evidence for caching."""
    run_dir = Path(run_dir)
    reasons = []
    stage = run_dir / "stages" / task_name / "manifest.json"
    if stage.exists():
        try:
            status = json.loads(stage.read_text(encoding="utf-8")).get("status")
            if status != "DONE":
                reasons.append("previous stage status " + str(status))
        except (OSError, ValueError, AttributeError):
            reasons.append("previous stage manifest is unreadable or invalid")
    manifest = run_dir / "run_manifest.json"
    if manifest.exists():
        try:
            records = json.loads(manifest.read_text(encoding="utf-8"))["tasks"]
            if not isinstance(records, list):
                raise ValueError("tasks must be a list")
            for record in records:
                if record.get("task") == task_name and record.get("status") not in SUCCESS_STATUSES:
                    reasons.append("previous run status " + str(record.get("status")))
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            reasons.append("previous run manifest is unreadable or invalid")
    return reasons


def has_successful_execution_record(task_name, run_dir):
    """Output existence alone does not prove an interrupted task completed."""
    run_dir = Path(run_dir)
    stage = run_dir / "stages" / task_name / "manifest.json"
    if stage.exists():
        try:
            return json.loads(stage.read_text(encoding="utf-8")).get("status") == "DONE"
        except (OSError, ValueError, AttributeError):
            return False
    try:
        records = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))["tasks"]
        matching = [record for record in records if record.get("task") == task_name]
        return bool(matching) and matching[-1].get("status") in SUCCESS_STATUSES
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return False


def section_support_errors(task, config):
    """Give a clear install error instead of invoking an older renderer CLI."""
    if not re.fullmatch(r"render_(9x130|9x275)_note_physics_(" + "|".join(SECTIONS) + ")(?:_keep_(?:yneg|xpos))?", task.get("name", "")):
        return []
    resolved = runner.resolved_task(task, config)
    path = runner.resolve_path(resolved["script"])
    if path.name != "render_note_physics.py" or not path.is_file():
        return []
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError):
        return []  # Existing validation reports these errors.
    supported = any(isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument"
                    and any(isinstance(arg, ast.Constant) and arg.value == "--section" for arg in node.args)
                    for node in ast.walk(tree))
    if not supported:
        return ["renderer lacks --section support; run python3 workflow/install_newid_update.py"]
    return []


def summarize_records(records, dry_run=False):
    statuses = {}
    for record in records:
        statuses.setdefault(record.get("status", "UNKNOWN"), []).append(record.get("task", "unknown"))
    successful = SUCCESS_STATUSES | ({"READY"} if dry_run else set())
    incomplete = [record.get("task", "unknown") for record in records
                  if record.get("status") not in successful | {"SKIPPED_EXPERIMENTAL"}]
    status = "PARTIAL" if incomplete else "DRY_RUN" if dry_run else "COMPLETE"
    return {"status": status, "total_tasks": len(records), "incomplete_tasks": incomplete,
            "counts_by_status": {key: len(value) for key, value in statuses.items()},
            "tasks_by_status": statuses, "keep_going": True}


@contextmanager
def install_continuation(dry_run=False):
    """Temporarily adapt the existing runner; restore every hook on exit."""
    original_load = runner.load_configuration
    original_run_task = runner.run_task
    original_manifest = runner.write_run_manifest
    original_validate = runner.validate_task
    state = {"summary": None}
    def load_configuration(args, run_id):
        return transform_configuration(original_load(args, run_id))
    def run_task(task, config, *, force, dry_run, hash_inputs):
        reasons = prior_failure_reasons(task["name"], config["run"]["dir"])
        if not reasons and runner.task_outputs_complete(runner.resolved_task(task, config)):
            if not has_successful_execution_record(task["name"], config["run"]["dir"]):
                reasons.append("no successful execution record")
        if reasons and not force:
            print("RETRY   " + task["name"] + ": " + "; ".join(reasons), flush=True)
        started = dt.datetime.now(dt.timezone.utc)
        try:
            return original_run_task(task, config, force=force or bool(reasons),
                                     dry_run=dry_run, hash_inputs=hash_inputs)
        except OSError as error:
            # A missing interpreter, inaccessible working directory or failed
            # artifact read is a task failure. User interrupts and programming
            # exceptions retain their normal behavior.
            finished = dt.datetime.now(dt.timezone.utc)
            resolved = runner.resolved_task(task, config)
            stage_dir = Path(config["run"]["dir"]) / "stages" / task["name"]
            log = stage_dir / "stage.log"
            reason = type(error).__name__ + " during task execution: " + str(error)
            record = {"task": task["name"], "stage": task.get("stage", ""),
                      "status": "FAILED", "reasons": [reason],
                      "started_utc": started.isoformat(), "finished_utc": finished.isoformat(),
                      "duration_seconds": (finished-started).total_seconds(),
                      "command": runner.command_for(resolved, config), "log": str(log),
                      "error_type": type(error).__name__, "outputs": []}
            try:
                stage_dir.mkdir(parents=True, exist_ok=True)
                with log.open("a", encoding="utf-8") as stream:
                    stream.write("\nFAILED: " + reason + "\n")
                (stage_dir / "manifest.json").write_text(
                    json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            except OSError as write_error:
                record["stage_record_write_error"] = str(write_error)
            print("FAILED  " + task["name"] + ": " + reason, flush=True)
            return "FAILED", record
    def validate_task(task, config):
        return original_validate(task, config) + section_support_errors(task, config)
    def write_run_manifest(config, records):
        path = original_manifest(config, records)
        summary = summarize_records(records, dry_run=dry_run)
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["continuation_summary"] = summary
        manifest["run_status"] = summary["status"]
        path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        state["summary"] = summary
        counts = ", ".join(str(value) + " " + key for key, value in summary["counts_by_status"].items())
        print("Run summary: " + summary["status"] + " (" + counts + ")", flush=True)
        return path
    runner.load_configuration = load_configuration
    runner.run_task = run_task
    runner.write_run_manifest = write_run_manifest
    runner.validate_task = validate_task
    try:
        yield state
    finally:
        runner.load_configuration = original_load
        runner.run_task = original_run_task
        runner.write_run_manifest = original_manifest
        runner.validate_task = original_validate


def main(argv=None):
    parser = runner.build_parser()
    parser.description = __doc__
    parser.set_defaults(pipeline=str(Path(__file__).resolve().parents[1] / "config/pipeline_statistics.toml"))
    args = parser.parse_args(argv)
    try:
        with install_continuation(dry_run=getattr(args, "dry_run", False)) as state:
            if args.command == "run":
                args.keep_going = True
                args.run_id = args.run_id or dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                config = runner.load_configuration(args, args.run_id)
                settings = runner.flatten({key: value for key, value in config.items()
                                           if key not in ("groups", "task", "run", "paths")})
                settings = [key + "=" + str(value) for key, value in sorted(settings.items())]
                contract = build_contract(settings, runner.PROJECT_ROOT, config.get("task", []))
                if not args.dry_run:
                    status = check_contract(Path(config["run"]["dir"]), contract)
                    print("Statistical contract " + status, flush=True)
            result = int(args.func(args))
            if args.command == "run":
                if state["summary"] and state["summary"]["status"] == "PARTIAL":
                    return result or 1
                if result and state["summary"] is None:
                    print("Run summary: PARTIAL (runner stopped before writing a new task manifest)", flush=True)
            return result
    except (runner.ConfigurationError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
