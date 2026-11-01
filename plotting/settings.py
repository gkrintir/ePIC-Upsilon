"""One explicit profile for the SWAN notebooks and command-line renderers.

A plotting profile selects saved numerical products; it does not change their
selection, luminosity, statistical inference, or validation status.
"""
from copy import deepcopy
import json
import os
from pathlib import Path
import re
import tempfile

import numpy as np
from upsilon_esr import Config

STATES = ("1S", "2S", "3S")
DEFAULT_PROFILE = Path(__file__).resolve().parents[1] / "configs/plots.local.json"


def _path(value, parent, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(label + " must be a nonempty path")
    expanded = os.path.expandvars(os.path.expanduser(value))
    if "$" in expanded:
        raise ValueError(label + " contains an unresolved environment variable")
    path = Path(expanded)
    return (path if path.is_absolute() else parent / path).resolve()


def _run(value, label):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,159}", value):
        raise ValueError(label + " must be one explicit run ID, not an absolute path or latest-run search")
    return value


def load_profile(path=None):
    """Read a local profile; never discover an alternative run automatically."""
    chosen = path if path is not None else os.environ.get("UPSILON_PLOT_PROFILE", DEFAULT_PROFILE)
    source = Path(chosen).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(
            f"Plotting profile not found: {source}. Copy configs/plots_9x275.example.json "
            "to configs/plots.local.json, edit its paths/run IDs, or set UPSILON_PLOT_PROFILE.")
    profile = json.loads(source.read_text())
    if not isinstance(profile, dict) or profile.get("schema_version") != 1:
        raise ValueError("Plotting profile requires schema_version=1")
    if profile.get("energy") not in ("9x130", "9x275"):
        raise ValueError("Plotting profile energy must be 9x130 or 9x275")
    lumi = profile.get("projection_lumi_invfb")
    if isinstance(lumi, bool) or not isinstance(lumi, (int, float)) or not np.isfinite(lumi) or lumi <= 0:
        raise ValueError("projection_lumi_invfb must be finite and positive")
    # This migrated package is the supplied proton-tagged notebook category.
    # An inclusive plot cannot be obtained by turning off provenance checks.
    if profile.get("proton_selection") != "required":
        raise ValueError("These preserved notebooks require proton_selection=required. "
                         "Do not relabel inclusive inputs or bypass the tagged metadata checks.")
    for name in ("full_runs", "q2_max", "categories"):
        if not isinstance(profile.get(name), dict) or set(profile[name]) != set(STATES):
            raise ValueError(name + " must explicitly specify 1S, 2S and 3S; one-state runs need the status/audit view")
    for state in STATES:
        _run(profile["full_runs"][state], "full_runs." + state)
        q2 = profile["q2_max"][state]
        if isinstance(q2, bool) or not isinstance(q2, (int, float)) or not np.isfinite(q2) or q2 <= 0:
            raise ValueError("q2_max." + state + " must be finite and positive")
        if not isinstance(profile["categories"][state], str) or not profile["categories"][state].strip():
            raise ValueError("categories." + state + " must explicitly name the saved sample category")
    _run(profile.get("reoptimization_run"), "reoptimization_run")
    if profile.get("efficiency_2d_run") is not None:
        _run(profile["efficiency_2d_run"], "efficiency_2d_run")
    if profile.get("map_sample") is not None:
        _run(profile["map_sample"], "map_sample")
    if type(profile.get("allow_diagnostic_results", False)) is not bool:
        raise ValueError("allow_diagnostic_results must be true or false")
    profile["allow_diagnostic_results"] = profile.get("allow_diagnostic_results", False)
    merges = profile.get("display_rebin", {})
    if not isinstance(merges, dict) or any(s != "1S" for s in merges):
        raise ValueError("display_rebin supports only explicit 1S W/y display merges")
    for state, observables in merges.items():
        if not isinstance(observables, dict) or any(v not in ("w", "y") for v in observables):
            raise ValueError("display_rebin.1S supports only w and y")
        for observable, edges in observables.items():
            if edges is False or edges is None or edges == "keep_outer_bins":
                continue
            if not isinstance(edges, list):
                raise ValueError("Display rebin must be false, null, keep_outer_bins or saved bin edges")
            a = np.asarray(edges, float)
            if a.ndim != 1 or len(a) < 2 or not np.isfinite(a).all() or np.any(np.diff(a) <= 0):
                raise ValueError("Display bin edges must be finite and increasing")
    profile["display_rebin"] = merges
    profile["results_root"] = _path(profile.get("results_root"), source.parent, "results_root")
    profile["output_root"] = _path(profile.get("output_root"), source.parent, "output_root")
    if profile.get("root_input_override") is not None:
        profile["root_input_override"] = _path(profile["root_input_override"], source.parent, "root_input_override")
    if not isinstance(profile.get("root_tree", "UpsilonTree"), str) or not profile.get("root_tree", "UpsilonTree"):
        raise ValueError("root_tree must be a nonempty tree name")
    profile["profile_path"] = source
    return profile


def _profile(profile):
    return load_profile() if profile is None else (load_profile(profile) if isinstance(profile, (str, Path)) else deepcopy(profile))


def new_output_dir(profile, section):
    """Create an isolated rendering directory without replacing old figures."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", section):
        raise ValueError("Invalid figure section")
    parent = Path(profile["output_root"]) / section
    parent.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="render-", dir=parent))


def make_config(framework_root=None, run_ids=None, energy=None, lumi_invfb=None,
                output_dir=None, *, profile=None):
    """Full-framework performance input; optional positional arguments stay compatible."""
    p = _profile(profile)
    root = (Path(framework_root).expanduser() / "results/generated/runs"
            if framework_root is not None else p["results_root"])
    selected = p["full_runs"] if run_ids is None else run_ids
    return Config(runs={state: root / name for state, name in selected.items()},
                  energy=p["energy"] if energy is None else energy,
                  lumi_invfb=p["projection_lumi_invfb"] if lumi_invfb is None else lumi_invfb,
                  output_dir=Path(output_dir) if output_dir is not None else p["output_root"],
                  campaign=p.get("campaign", "26.07.1 Campaign"),
                  require_proton_tag=True, allow_diagnostic_results=p["allow_diagnostic_results"],
                  categories=deepcopy(p["categories"]), q2_max=deepcopy(p["q2_max"]))


def make_reoptimization_config(framework_root=None, run_id=None, energy=None,
                               lumi_invfb=None, output_dir=None, *, profile=None):
    """Saved inference at its actual luminosity; interval rescaling is forbidden."""
    p = _profile(profile)
    root = (Path(framework_root).expanduser() / "results/generated/runs"
            if framework_root is not None else p["results_root"])
    run = root / (p["reoptimization_run"] if run_id is None else run_id)
    cfg = make_config(energy=energy, lumi_invfb=lumi_invfb, output_dir=output_dir, profile=p)
    cfg.runs = {state: run for state in STATES}
    cfg.reoptimization_run = run
    return cfg


def apply_display_rebin(config, profile, section):
    """Only merge adjacent saved edges; false really keeps the original binning."""
    config.rebin_edges = {}
    wanted = ("w",) if section == "w" else (("y",) if section == "differential" else ())
    for state, choices in profile.get("display_rebin", {}).items():
        for observable, requested in choices.items():
            if observable not in wanted or requested is False or requested is None:
                continue
            from upsilon_esr.reoptimization_inputs import load_case
            from upsilon_esr.rebinning import merge_matrix
            table, _, _, _ = load_case(config, state, observable)
            old = table.x_low.tolist() + [float(table.x_high.iloc[-1])]
            if requested == "keep_outer_bins":
                new = [old[0], old[1], old[-2], old[-1]] if len(old) > 4 else old
            else:
                new = list(requested)
            merge_matrix(old, new)
            config.rebin_edges.setdefault(state, {})[observable] = new
    return config


def map_source(profile):
    if not profile.get("efficiency_2d_run") or not profile.get("map_sample"):
        raise ValueError("The combined efficiency notebooks require efficiency_2d_run and map_sample in the shared profile")
    return profile["results_root"] / profile["efficiency_2d_run"] / "artifacts/final_efficiency_2d" / profile["map_sample"]


def file_link_path(path):
    """Notebook-relative link; the printed absolute output path remains authoritative."""
    return os.path.relpath(Path(path).resolve(), Path.cwd())
