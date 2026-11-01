"""Validate July legacy/pure-state normalization and guard existence-only caches."""
import hashlib
import json
import math
from pathlib import Path

PROFILE = "jul26-estarlight-939b11a"
ALL_STATES_PROFILE = "jul26-purestates-estarlight-939b11a"
STATE_REFERENCE = {
    "1S": {"br_ee": 0.0238, "q2_max": 1.0,
           "visible_xsec_fb": {"9x130": 119.685, "9x275": 210.136}},
    "2S": {"br_ee": 0.0191, "q2_max": 0.01,
           "visible_xsec_fb": {"9x130": 29.297, "9x275": 53.133}},
    "3S": {"br_ee": 0.0218, "q2_max": 0.01,
           "visible_xsec_fb": {"9x130": 23.447, "9x275": 42.969}},
}


def task_definitions_sha256(task_definitions):
    """Bind command arguments and dependencies independently of target selection."""
    encoded = json.dumps(task_definitions, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def configuration_json_sha256(project_root):
    """Bind external sample/rate catalogs used by configured analysis tasks."""
    root = Path(project_root)
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((root / "config").rglob("*.json"))}


def validate_flat_normalization(values):
    profile = values.get("analysis.normalization_profile")
    if profile not in (PROFILE, ALL_STATES_PROFILE):
        return
    def require_number(key, expected):
        try:
            value = float(values[key])
        except (KeyError, TypeError, ValueError):
            raise ValueError("Verified July normalization requires " + key)
        if not math.isfinite(value) or not math.isclose(value, expected, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("%s requires %s=%s, got %s. For a deliberately different calculation, select an explicit custom normalization_profile." % (profile, key, expected, value))
    state = str(values.get("analysis.states", "1S"))
    if profile == PROFILE:
        if state != "1S":
            raise ValueError("The latest July workflow prepares 1S-only tuples. Use run_state.py for pure 2S/3S inputs with their own Q2 settings.")
        for slot, ep, rate in (("low", 130., 119.685), ("high", 275., 210.136)):
            require_number("beam.%s.electron_energy" % slot, 9.)
            require_number("beam.%s.proton_energy" % slot, ep)
            require_number("beam.%s.visible_xsec_1s_fb" % slot, rate)
        require_number("analysis.generator_branching_fraction", 0.0238)
    else:
        if state not in STATE_REFERENCE:
            raise ValueError("The pure-state July workflow requires exactly one of 1S, 2S or 3S")
        reference = STATE_REFERENCE[state]
        for slot, ep, energy in (("low", 130., "9x130"), ("high", 275., "9x275")):
            require_number("beam.%s.electron_energy" % slot, 9.)
            require_number("beam.%s.proton_energy" % slot, ep)
            require_number("beam.%s.visible_xsec_fb" % slot, reference["visible_xsec_fb"][energy])
        require_number("analysis.generator_branching_fraction", reference["br_ee"])
        require_number("analysis.q2_event_max", reference["q2_max"])
        require_number("analysis.q2_flux_max", reference["q2_max"])
    try:
        entries = [item.strip().split(":") for item in str(values["analysis.branching_fractions"]).split(",")]
        branching = {state: float(value) for state, value in entries}
        if len(branching) != len(entries):
            raise ValueError("duplicate state")
    except (KeyError, TypeError, ValueError):
        raise ValueError("Verified July normalization requires a valid analysis.branching_fractions mapping")
    if branching != {"1S": 0.0238, "2S": 0.0191, "3S": 0.0218}:
        raise ValueError("Verified July normalization requires generator branching fractions 1S:0.0238,2S:0.0191,3S:0.0218")


def check_normalization_cache(config, project_root):
    """Guard direct run_analysis runs; statistics wrappers have a broader contract."""
    if not config.get("analysis", {}).get("normalization_profile"):
        return "LEGACY"
    selected = {key: value for key, value in config.items() if key not in ("run", "paths", "groups", "task")}
    sources = {}
    for folder in ("scripts", "workflow"):
        for path in sorted((Path(project_root) / folder).rglob("*.py")):
            sources[str(path.relative_to(project_root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    contract = {"schema_version": 1, "configuration": selected, "analysis_source_sha256": sources,
                "task_definitions_sha256": task_definitions_sha256(config.get("task", [])),
                "configuration_json_sha256": configuration_json_sha256(project_root)}
    run_dir = Path(config["run"]["dir"])
    path = run_dir / "normalization_contract.json"
    if path.exists():
        if json.loads(path.read_text()) != contract:
            raise ValueError("Normalization settings or code changed. Use a fresh --run-id to regenerate normalized inputs, yields and projected uncertainties.")
        return "MATCH"
    if (run_dir / "run_manifest.json").exists() or any((run_dir / "artifacts").glob("**/*")):
        raise ValueError("Existing results have no current normalization contract. Use a fresh --run-id.")
    run_dir.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(contract, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return "CREATED"
