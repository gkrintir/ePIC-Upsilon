#!/usr/bin/env python3
"""Run the complete note workflow for one pure Upsilon state at both energies.

The state selects its own generated sample, visible rate, branching fraction,
Q2 phase space and plot labels. A run directory contains exactly one state.
Use run_all_states.py to execute the three state workflows as one batch.
"""
from contextlib import contextmanager
from copy import copy, deepcopy
import argparse
import datetime as dt
import hashlib
import math
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
        raise SystemExit("Python <3.11 needs tomli: python3 -m pip install --user tomli")
    sys.modules["tomllib"] = tomllib

import run_statistics_continue as continuation
from normalization_contract import ALL_STATES_PROFILE, STATE_REFERENCE

runner = continuation.runner
ROOT = Path(__file__).resolve().parents[1]
BEAMS = {"9x130": "low", "9x275": "high"}
RENDERERS = {"render_note_kinematics.py", "render_note_mass_efficiency.py",
             "render_note_physics.py", "render_note_t_studies.py"}


def remove_options(arguments, options):
    result = []
    skip = False
    for argument in arguments:
        if skip:
            skip = False
            continue
        if argument in options:
            skip = True
        elif not any(isinstance(argument, str) and argument.startswith(option + "=") for option in options):
            result.append(argument)
    return result


def replace_option(arguments, option, value):
    return remove_options(arguments, {option}) + [option, str(value)]


def audit_source_fingerprints(project_root):
    folder = Path(project_root) / "tools/normalization_audit"
    return {str(path.relative_to(project_root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(folder.glob("*.py"))}


def state_configuration(base, state, defaults):
    """Select one pure state without altering the legacy pipeline on disk."""
    if state not in STATE_REFERENCE:
        raise runner.ConfigurationError("Unsupported state: " + state)
    config = runner.deep_merge(deepcopy(base), defaults)
    selected = STATE_REFERENCE[state]
    analysis = config.setdefault("analysis", {})
    analysis.update(normalization_profile=ALL_STATES_PROFILE, states=state,
                    generator_branching_fraction=selected["br_ee"],
                    branching_fractions="1S:0.0238,2S:0.0191,3S:0.0218",
                    q2_flux_max=selected["q2_max"], q2_event_max=selected["q2_max"])
    pattern = config.get("all_states", {}).get("input_pattern")
    if not pattern:
        raise runner.ConfigurationError("config/all_states.toml must define all_states.input_pattern")
    for energy, slot in BEAMS.items():
        config["inputs"]["root_" + energy] = pattern.format(state=state, energy=energy, directory=config["all_states"]["input_directories"][state])
        beam = config["beam"][slot]
        beam.pop("visible_xsec_1s_fb", None)
        beam["visible_xsec_fb"] = selected["visible_xsec_fb"][energy]
        beam["normalization_source"] = (
            "Reproduced eSTARlight 939b11a with Jul26 Upsilon_Sim_Cards f62965a "
            "(%s, Upsilon(%s), Q2 0-%s GeV^2): visible dielectron %s fb; generator B_ee=%s"
            % (energy, state, selected["q2_max"], beam["visible_xsec_fb"], selected["br_ee"]))
    return config


def transform_tasks(config, state):
    """Select state-aware CLI arguments and outputs, keeping pure-state caches."""
    result = deepcopy(config)
    token = state.lower()
    def replace(value):
        if not isinstance(value, str):
            return value
        for slot in BEAMS.values():
            value = value.replace("{beam." + slot + ".visible_xsec_1s_fb}",
                                  "{beam." + slot + ".visible_xsec_fb}")
        value = value.replace("1S", state).replace("1s", token)
        return value
    for task in result["task"]:
        for key in ("args", "requirements", "outputs"):
            if key in task:
                task[key] = [replace(value) for value in task[key]]
        task["description"] = task.get("description", task["name"]) + " [Upsilon(" + state + ")]"
        script = Path(task["script"]).name
        args = task.get("args", [])
        if script in RENDERERS:
            args = replace_option(args, "--state", state)
        if "--min-eff" in args:
            args = replace_option(args, "--min-eff", "{analysis.min_eff}")
        if script == "extract_gamma_p_cross_sections_bin_by_bin.py":
            # Record the complete pure-state generated denominator. The truth
            # mass/Q2/W masks remain selections, never a renormalization.
            args = replace_option(args, "--normalization-scope", "all")
        if state != "1S":
            args = remove_options(args, {"--zeus-csv", "--zeus-mode"})
            task["requirements"] = [x for x in task.get("requirements", [])
                                    if "zeus_step_csv" not in x]
            task["outputs"] = [x for x in task.get("outputs", [])
                               if not x.endswith(("_zeus_overlay.png", "_zeus_overlay.pdf"))
                               and Path(x).name not in {
                                   "momentum_transfer_zeus_shape_covariance.csv",
                                   "impact_parameter_zeus_covariance.csv"}]
        task["args"] = args
    # The same section splitter used by the installed continuation workflow is
    # applied after state-specific output pruning. Repeated transforms are safe.
    result = continuation.transform_configuration(result)
    add_input_and_audit_tasks(result, state)
    runner.task_index(result)
    return result


def add_input_and_audit_tasks(config, state):
    tasks = config["task"]
    names = {task["name"] for task in tasks}
    audits = []
    gates = []
    for energy, slot in BEAMS.items():
        gate = "validate_" + energy + "_state_input"
        audit = "audit_" + energy + "_normalization"
        if gate in names or audit in names:
            raise runner.ConfigurationError("State task name already exists: " + gate)
        gates.append(gate)
        tasks.append({
            "name": gate, "stage": "schema", "description": "Validate pure " + state + " generated input and Q2 phase space",
            "script": "scripts/validation/validate_state_sample.py",
            "args": ["--input", "{inputs.root_" + energy + "}", "--tree", "{analysis.tree}",
                     "--state", state, "--q2-max", "{analysis.q2_flux_max}",
                     "--out", "{run.dir}/artifacts/schema/" + energy + "_state_input.json"],
            "requirements": ["module:numpy", "module:uproot", "file:{inputs.root_" + energy + "}"],
            "outputs": ["{run.dir}/artifacts/schema/" + energy + "_state_input.json"],
            "fail_patterns": ["Traceback"],
        })
        prepare = next((task for task in tasks if task["name"] == "prepare_" + energy), None)
        if prepare is None:
            raise runner.ConfigurationError("Expected preparation task for " + energy)
        prepare["deps"] = list(dict.fromkeys(prepare.get("deps", []) + [gate]))
        outdir = "{run.dir}/artifacts/normalization_audit/" + energy
        # The bin-by-bin extractor uses linspace(w_min, w_max, n_w_bins+1).
        # Its audit must follow that definition even if unfolding w_edges differ.
        beam = config["beam"][slot]
        count = int(beam["n_w_bins"])
        lower, upper = float(beam["w_min"]), float(beam["w_max"])
        if count <= 0 or count != float(beam["n_w_bins"]) or not all(map(math.isfinite, (lower, upper))) or upper <= lower:
            raise runner.ConfigurationError("Invalid bin-by-bin W range or bin count for " + energy)
        audit_edges = ",".join(format(lower + (upper-lower)*i/count, ".17g") for i in range(count+1))
        tasks.append({
            "name": audit, "stage": "normalization", "description": "Audit raw/prepared " + state + " weights, flux and W central values",
            "script": "tools/normalization_audit/audit_upsilon_normalization.py",
            "deps": ["extract_" + energy + "_sigma_gamma_p"],
            "args": ["--state", state, "--merged", "{inputs.root_" + energy + "}",
                     "--run-dir", "{run.dir}", "--tree", "{analysis.tree}",
                     "--electron-energy", "{beam." + slot + ".electron_energy}",
                     "--proton-energy", "{beam." + slot + ".proton_energy}",
                     "--assigned-visible-xsec-fb", "{beam." + slot + ".visible_xsec_fb}",
                     "--generator-visible-xsec-fb", "{beam." + slot + ".visible_xsec_fb}",
                     "--generator-source", "{beam." + slot + ".normalization_source}",
                     "--lumi-invfb", "{analysis.tuple_lumi_invfb}",
                     "--br-ee", "{analysis.generator_branching_fraction}",
                     "--q2-cut", "{analysis.q2_event_max}", "--q2-max", "{analysis.q2_flux_max}",
                     "--true-mass-window", "{analysis.truth_mass_window}",
                     "--mass-window-sigma", "{analysis.reco_mass_window_sigma}", "--min-eff", "{analysis.min_eff}",
                     "--w-edges", audit_edges, "--out-dir", outdir, "--hash-inputs"],
            "requirements": ["module:numpy", "module:uproot"],
            "outputs": [outdir + "/audit.json", outdir + "/REPORT.md", outdir + "/saved_csv_comparison.csv"],
            "fail_patterns": ["Traceback"],
        })
        audits.append(audit)
    groups = config.setdefault("groups", {})
    groups["state_inputs"] = gates
    groups["normalization_audits"] = audits
    groups["state_all"] = list(dict.fromkeys(groups["note_all"] + audits))
    for energy in BEAMS:
        groups["state_" + energy] = list(dict.fromkeys(groups["note_" + energy] + ["audit_" + energy + "_normalization"]))


@contextmanager
def install_state_defaults(state):
    original_load = runner.load_configuration
    cached_fingerprints = {}
    def load_configuration(args, run_id):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,159}", run_id):
            raise runner.ConfigurationError("Use a run ID containing only letters, digits, '.', '_' or '-'")
        base_args = copy(args)
        base_args.config = None
        base_args.set_values = []
        base = original_load(base_args, run_id)
        overrides = runner.read_toml(runner.resolve_path(args.config)) if args.config else {}
        for raw in args.set_values or []:
            key, value = runner.parse_override(raw)
            runner.set_nested(overrides, key, value)
        defaults = runner.deep_merge(runner.read_toml(ROOT / "config/all_states.toml"),
                                     {"all_states": overrides.get("all_states", {})})
        config = runner.deep_merge(state_configuration(base, state, defaults), overrides)
        # The general cache contract already binds scripts/ and workflow/;
        # the independent normalization auditor also lives under tools/.
        config["state_audit_source_sha256"] = audit_source_fingerprints(ROOT)
        if config.get("analysis", {}).get("states") != state:
            raise runner.ConfigurationError("Select the state with --state; analysis.states must match it")
        if config.get("analysis", {}).get("normalization_profile") != ALL_STATES_PROFILE:
            raise runner.ConfigurationError("run_state.py requires the validated pure-state July normalization profile")
        runner.validate_flat_normalization(runner.flatten(config))
        config = transform_tasks(config, state)
        runner.check_category_configuration(config)
        # Bind actual raw-input identities before continuation checks caches.
        if args.command == "run" and not getattr(args, "dry_run", False):
            selected = runner.dependency_closure(config, runner.expand_targets(config, args.targets))
            index = runner.task_index(config)
            relevant = set()
            for name in selected:
                for value in index[name].get("args", []) + index[name].get("requirements", []):
                    if isinstance(value, str):
                        relevant.update(re.findall(r"\{inputs\.(root_9x(?:130|275))\}", value))
            contract = Path(config["run"]["dir"]) / "statistical_contract.json"
            if contract.exists():
                prefix = "state_source_fingerprints."
                settings = json.loads(contract.read_text()).get("settings", {})
                relevant.update(key[len(prefix):].split(".", 1)[0]
                                for key in settings if key.startswith(prefix))
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
            config["state_source_fingerprints"] = fingerprints
        return config
    runner.load_configuration = load_configuration
    try:
        yield
    finally:
        runner.load_configuration = original_load


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    selector = argparse.ArgumentParser(add_help=False)
    selector.add_argument("--state", choices=tuple(STATE_REFERENCE), default="1S")
    selected, remaining = selector.parse_known_args(argv)
    parser = runner.build_parser()
    parser.description = __doc__
    parser.epilog = "State selection: add --state 1S, --state 2S or --state 3S before the command (default 1S)."
    parser.set_defaults(pipeline=str(ROOT / "config/pipeline_statistics.toml"))
    args = parser.parse_args(remaining)
    if args.command in ("run", "plan", "validate"):
        if not args.targets:
            remaining.insert(remaining.index(args.command) + 1, "state_all")
        if args.run_id is None:
            stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            remaining += ["--run-id", "state-" + selected.state.lower() + "-" + stamp]
    with install_state_defaults(selected.state):
        return continuation.main(remaining)


if __name__ == "__main__":
    raise SystemExit(main())
