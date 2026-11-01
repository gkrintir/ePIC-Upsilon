#!/usr/bin/env python3
"""Run the installed full statistics/note workflow on reconstructed-electron-ID samples.

Uses the installed pipeline and physics scripts. Adds newID input validation,
electron-ID comparisons, complete legacy-ID fEMC half-plane workflows, and
dependency-aware continuation without replacing the original runner or pipeline.
Precedence: installed pipeline, comparison defaults, newid.toml, optional
--config, then --set.
"""
from contextlib import contextmanager
from copy import copy, deepcopy
import datetime as dt
import json
from pathlib import Path
import re
import sys

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        raise SystemExit("Python <3.11 needs tomli: python3 -m pip install tomli")
    sys.modules["tomllib"] = tomllib

try:
    import newid_continuation as continuation
except ModuleNotFoundError as error:
    raise SystemExit("Missing framework dependency %s. Install this update in the existing framework with its statistics update." % error.name)

runner = continuation.runner
ROOT = Path(__file__).resolve().parents[1]
BEAMS = {"9x130": "low", "9x275": "high"}


def add_newid_tasks(config, comparison):
    result = deepcopy(config)
    tasks = result.setdefault("task", [])
    names = {task["name"] for task in tasks}
    for energy, slot in BEAMS.items():
        name = "validate_" + energy + "_newid"
        if name in names:
            raise runner.ConfigurationError("Task already exists: " + name)
        tasks.append({
            "name": name, "stage": "schema", "description": "Check new electron-ID mode, generated denominator and stored beam kinematics",
            "script": "scripts/validation/validate_newid_sample.py",
            "args": ["--input", "{inputs.root_" + energy + "}", "--tree", "{analysis.tree}",
                     "--electron-energy", "{beam." + slot + ".electron_energy}",
                     "--proton-energy", "{beam." + slot + ".proton_energy}",
                     "--out-json", "{run.dir}/artifacts/schema/" + energy + "_newid.json"],
            "requirements": ["module:numpy", "module:uproot", "file:{inputs.root_" + energy + "}"],
            "outputs": ["{run.dir}/artifacts/schema/" + energy + "_newid.json"],
            "fail_patterns": ["Traceback"],
        })
        prepare = next((task for task in tasks if task["name"] == "prepare_" + energy), None)
        if prepare is None:
            raise runner.ConfigurationError("Expected prepare_" + energy + " in the installed statistics pipeline")
        prepare["deps"] = list(dict.fromkeys(prepare.get("deps", []) + [name]))
    dependencies = {
        "compare_id_9x130_samples": ["validate_9x130_newid"],
        "compare_id_9x275_samples": ["validate_9x275_newid"],
        "compare_newid_energies_samples": ["validate_9x130_newid", "validate_9x275_newid"],
        "compare_9x275_keep_yneg_samples": ["validate_9x275_legacy_full", "validate_9x275_keep_yneg"],
        "compare_9x275_keep_xpos_samples": ["validate_9x275_legacy_full", "validate_9x275_keep_xpos"],
    }
    # The supplied half-plane producer supports legacy electron ID only. Match
    # that ID in the full-detector reference so these compare the mask alone.
    for scenario, input_key in (("full", "baseline_9x275_root"),
                                ("keep_yneg", "keep_yneg_9x275_root"),
                                ("keep_xpos", "keep_xpos_9x275_root")):
        tag = "legacy_full" if scenario == "full" else scenario
        name = "validate_9x275_" + tag
        if name in names:
            raise runner.ConfigurationError("Task already exists: " + name)
        args = ["--input", "{inputs." + input_key + "}", "--tree", "{analysis.tree}",
                "--electron-energy", "{beam.high.electron_energy}",
                "--proton-energy", "{beam.high.proton_energy}",
                "--electron-id", "legacy", "--femc-scenario", scenario,
                "--allow-missing-legacy-id-metadata", "--hash-input",
                "--out-json", "{run.dir}/artifacts/schema/9x275_" + tag + ".json"]
        if scenario == "full":
            args.append("--allow-missing-full-femc-metadata")
        tasks.append({
            "name": name, "stage": "schema",
            "description": "Check 9x275 legacy ID / " + scenario + " metadata and generated denominator",
            "script": "scripts/validation/validate_newid_sample.py", "args": args,
            "requirements": ["module:numpy", "module:uproot", "file:{inputs." + input_key + "}"],
            "outputs": ["{run.dir}/artifacts/schema/9x275_" + tag + ".json"],
            "fail_patterns": ["Traceback"],
        })
    for original in comparison["task"]:
        task = deepcopy(original)
        if task["name"] in names:
            raise runner.ConfigurationError("Comparison task already exists: " + task["name"])
        task["deps"] = dependencies[task["name"]]
        tasks.append(task)
    groups = result.setdefault("groups", {})
    for name, members in comparison["groups"].items():
        if name in groups:
            raise runner.ConfigurationError("Comparison group already exists: " + name)
        groups[name] = members
    groups["newid_inputs"] = ["validate_" + energy + "_newid" for energy in BEAMS]
    groups["halfplane_inputs_9x275"] = ["validate_9x275_legacy_full", "validate_9x275_keep_yneg", "validate_9x275_keep_xpos"]
    groups["newid_all"] = list(dict.fromkeys(groups["note_all"] + groups["compare_ele_id_all"]))
    return result


def add_halfplane_workflows(config):
    """Clone the installed 9x275 note workflow into independent scenario paths.

    Split the physics renderer first, so each cloned section keeps only its
    own upstream dependencies. The source scripts and numerical arguments are
    inherited unchanged. Each scenario prepares its complete generated tuple
    using the same assigned cross section and its own generated denominator.
    """
    result = continuation.transform_configuration(config)
    groups = result.setdefault("groups", {})
    if "note_9x275" not in groups:
        raise runner.ConfigurationError("Expected note_9x275 in the installed statistics pipeline")
    targets = runner.expand_targets(result, ["note_9x275"])
    closure = runner.dependency_closure(result, targets)
    index = runner.task_index(result)
    strict_gate = "validate_9x275_newid"
    if strict_gate not in closure:
        raise runner.ConfigurationError("The 9x275 workflow is missing its newID input gate")
    source_names = [name for name in closure if name != strict_gate]

    def install_group(name, members):
        members = list(dict.fromkeys(members))
        if name in groups and groups[name] != members:
            raise runner.ConfigurationError("Conflicting half-plane group: " + name)
        groups[name] = members

    combined = []
    scenario_metadata = {}
    for scenario, geometry in (("keep_yneg", "y<0"), ("keep_xpos", "x>=0")):
        gate = "validate_9x275_" + scenario
        if gate not in index:
            raise runner.ConfigurationError("Missing half-plane input gate: " + gate)
        rename = {name: name + "_" + scenario for name in source_names}
        rename[strict_gate] = gate
        prefix = "{run.dir}/artifacts/scenarios/" + scenario + "/"

        def replace_value(value):
            if not isinstance(value, str):
                return value
            value = value.replace("{inputs.root_9x275}", "{inputs." + scenario + "_9x275_root}")
            # --run-artifacts receives the bare root; input/output paths use
            # descendants. Both must point at the same scenario tree.
            value = re.sub(r"\{run\.dir\}/artifacts(?=/|$)",
                           lambda match: prefix.rstrip("/"), value)
            if value == "{beam.high.label}":
                value += ", legacy electron ID, fEMC " + geometry
            return value

        for name in source_names:
            original = index[name]
            task = deepcopy(original)
            task["name"] = rename[name]
            task["description"] = original.get("description", name) + " [legacy ID, " + scenario + "]"
            task["deps"] = [rename[dependency] for dependency in original.get("deps", [])]
            for key in ("args", "requirements", "outputs"):
                if key in task:
                    task[key] = [replace_value(value) for value in task[key]]
            if "cwd" in task:
                task["cwd"] = replace_value(task["cwd"])
            if task["name"] in index:
                # Loading and continuation may transform a configuration more
                # than once. Accept only an identical existing clone.
                if index[task["name"]] != task:
                    raise runner.ConfigurationError("Conflicting half-plane task: " + task["name"])
            else:
                result["task"].append(task)
                index[task["name"]] = task
        members = [rename[name] for name in targets]
        install_group("note_9x275_" + scenario, members)
        combined.extend(members)
        scenario_metadata[scenario] = {
            "electron_id": "legacy", "ele_id_mode": 0,
            "femc_scenario": scenario, "femc_mode": 4 if scenario == "keep_yneg" else 3,
            "input_key": scenario + "_9x275_root",
            "artifacts_subdirectory": "artifacts/scenarios/" + scenario,
            "normalization_policy": "Inherit beam.high visible cross section and normalize each complete generated sample separately",
        }
    install_group("note_halfplanes_9x275", combined)
    install_group("halfplanes_all_9x275", combined + groups["compare_halfplanes_9x275"])
    result["halfplane_workflows"] = scenario_metadata
    runner.task_index(result)
    return result


@contextmanager
def install_newid_defaults():
    original_load = runner.load_configuration
    cached_fingerprints = {}
    def load_configuration(args, run_id):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,159}", run_id):
            raise runner.ConfigurationError("Use a run ID containing only letters, digits, '.', '_' or '-'")
        base_args = copy(args)
        base_args.config = None
        base_args.set_values = []
        config = original_load(base_args, run_id)
        comparison = runner.read_toml(ROOT / "config/compare_newid.toml")
        defaults = {key: value for key, value in comparison.items()
                    if key in ("inputs", "comparison")}
        config = runner.deep_merge(config, defaults)
        config = runner.deep_merge(config, runner.read_toml(ROOT / "config/newid.toml"))
        if args.config:
            config = runner.deep_merge(config, runner.read_toml(runner.resolve_path(args.config)))
        for raw in args.set_values or []:
            key, value = runner.parse_override(raw)
            runner.set_nested(config, key, value)
        try:
            runner.validate_flat_normalization(runner.flatten(config))
        except ValueError as error:
            raise runner.ConfigurationError(str(error)) from error
        config = add_halfplane_workflows(add_newid_tasks(config, comparison))
        # Include actual selected source identities in the statistical cache
        # contract, so replacing a merged ROOT at the same path invalidates it.
        if args.command == "run" and not getattr(args, "dry_run", False):
            selected_config = continuation.transform_configuration(config)
            targets = runner.expand_targets(selected_config, args.targets)
            closure = runner.dependency_closure(selected_config, targets)
            relevant = set()
            selected_tasks = runner.task_index(selected_config)
            for name in closure:
                task = selected_tasks[name]
                for value in task.get("args", []) + task.get("requirements", []):
                    if isinstance(value, str):
                        relevant.update(key for key in re.findall(r"\{inputs\.([A-Za-z0-9_-]+)\}", value)
                                        if key.startswith("root_") or key.endswith("_root"))
            # A subset retry must retain and recheck the original input scope.
            # Adding previously unbound inputs still requires a fresh run ID.
            contract_path = Path(config["run"]["dir"]) / "statistical_contract.json"
            if contract_path.is_file():
                prior_settings = json.loads(contract_path.read_text(encoding="utf-8")).get("settings", {})
                prefix = "newid_source_fingerprints."
                relevant.update(key[len(prefix):].split(".", 1)[0]
                                for key in prior_settings if key.startswith(prefix))
            fingerprints = {}
            for key in sorted(relevant):
                path = runner.resolve_path(config["inputs"][key])
                try:
                    stat = path.stat()
                    signature = (str(path), stat.st_size, stat.st_mtime_ns, bool(args.hash_inputs))
                    if signature not in cached_fingerprints:
                        cached_fingerprints[signature] = runner.file_fingerprint(path, include_hash=args.hash_inputs)
                    fingerprints[key] = cached_fingerprints[signature]
                except OSError:
                    fingerprints[key] = {"path": str(path), "available": False}
            config["newid_source_fingerprints"] = fingerprints
        return config
    runner.load_configuration = load_configuration
    try:
        yield
    finally:
        runner.load_configuration = original_load


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Use standard runner syntax, but default to the full two-energy note group
    # and a distinct newID run name. Explicit targets/run IDs remain supported.
    parser = runner.build_parser()
    parser.set_defaults(pipeline=str(ROOT / "config/pipeline_statistics.toml"))
    args = parser.parse_args(argv)
    if args.command in ("run", "plan", "validate"):
        if not args.targets:
            command_at = argv.index(args.command)
            argv.insert(command_at + 1, "note_all")
        if args.run_id is None:
            stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            argv += ["--run-id", "newid-" + stamp]
    with install_newid_defaults():
        return continuation.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())

