#!/usr/bin/env python3
"""Optimize differential bins from prepared ROOT tuples and validate on heldout MC.

This additive runner does not overwrite legacy plots or mutate input ROOT files.
It makes no claim of complete experimental uncertainty or mass-state extraction.
"""

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import traceback

import numpy as np

FRAMEWORK = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FRAMEWORK / "scripts"))
sys.path.insert(0, str(FRAMEWORK / "workflow"))
from proton_selection import proton_selection_metadata
from proton_category_isolation import enforce_category_workspace
from reoptimization.data import load_sample, integrated_flux, cross_section_factor, observable_cross_section_scaling
from reoptimization.response import optimize, evaluate_edges, build_model
from reoptimization.inference import analyze, analyze_spectrum, fit_slope, poisson_unfold
from reoptimization.stress import build_stress_tests, propagate_status
from reoptimization.stress import validation_gates as separated_validation_gates
from reoptimization.reporting import write_report


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (np.ndarray, list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, value):
    enforce_category_workspace(FRAMEWORK, path)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def digest(value):
    return hashlib.sha256(json.dumps(json_safe(value), sort_keys=True).encode()).hexdigest()


def print_input_audit(state, sample):
    metadata = sample["metadata"]
    if "missing_truth_rows" in metadata:
        print("TRUTH", state, "available", metadata["available_truth_rows"],
              "missing -999", metadata["missing_truth_rows"],
              "fraction", "%.6g" % metadata["missing_truth_fraction"],
              "projected missing yield", "%.6g" % metadata["missing_truth_target_yield"],
              "(all-generated normalization unchanged)", flush=True)
    for warning in metadata.get("warnings", []):
        print("WARN", state, warning, flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "validate", "run"))
    parser.add_argument("--config", type=Path,
                        default=FRAMEWORK / "config/reoptimization_9x275_wy.json")
    parser.add_argument("--framework-root", type=Path, default=FRAMEWORK,
                        help="Framework root containing results/generated/runs")
    parser.add_argument("--states", nargs="+", choices=("1S", "2S", "3S"),
                        default=("1S", "2S", "3S"))
    parser.add_argument("--observables", nargs="+", choices=("pt2", "abs_t_via_pt2", "w", "y"),
                        default=None, help="Defaults to config default_observables; old configs retain pt2 and abs_t_via_pt2")
    parser.add_argument("--input", action="append", default=[], metavar="STATE=ROOT_PATH",
                        help="Override a prepared normalized ROOT input, repeat per state")
    parser.add_argument("--lumi", type=float, help="Override target fb^-1 (default 2.5)")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Explicit run output directory, containing artifacts/reoptimization")
    parser.add_argument("--resume", action="store_true",
                        help="Reuse complete state/observable results only with identical fingerprints")
    parser.add_argument("--hash-inputs", action="store_true")
    parser.add_argument("--common-binning", action="store_true",
                        help="Require one passing edge set across all requested states")
    parser.add_argument("--n-toys", type=int, help="Override extraction and coverage toys")
    parser.add_argument("--n-response-bootstrap", type=int)
    return parser.parse_args(argv)


def prepare_config(args):
    category = enforce_category_workspace(args.framework_root, args.out_dir)
    cfg = json.loads(args.config.expanduser().resolve().read_text())
    if category is not None:
        cfg["proton_category"] = category
    if cfg.get("schema_version") != 1:
        raise ValueError("Unsupported reoptimization config schema_version")
    if not np.isclose(cfg["electron_energy"], 9) or not np.isclose(cfg["proton_energy"], 275):
        raise ValueError("This configured update is for 9x275. Prepare a separately audited config for another beam.")
    if args.lumi is not None:
        cfg["target_lumi_invfb"] = args.lumi
    if not np.isfinite(cfg["target_lumi_invfb"]) or cfg["target_lumi_invfb"] <= 0:
        raise ValueError("Target luminosity must be positive and finite")
    cfg["hash_inputs"] = bool(args.hash_inputs or cfg.get("hash_inputs", False))
    cfg["common_binning"] = bool(args.common_binning or cfg.get("common_binning", False))
    if args.n_toys is not None:
        if args.n_toys < 1:
            raise ValueError("--n-toys must be positive")
        cfg["inference"]["n_toys"] = args.n_toys
        cfg["inference"]["n_coverage_toys"] = args.n_toys
    if args.n_response_bootstrap is not None:
        if args.n_response_bootstrap < 0:
            raise ValueError("--n-response-bootstrap cannot be negative")
        cfg["inference"]["n_response_bootstrap"] = args.n_response_bootstrap
    if args.observables is None:
        args.observables = list(cfg.get("default_observables", ["pt2", "abs_t_via_pt2"]))
    if not args.observables or any(x not in ("pt2", "abs_t_via_pt2", "w", "y") for x in args.observables):
        raise ValueError("Unknown or empty observable selection")
    if len(set(args.observables)) != len(args.observables):
        raise ValueError("Observable selection must not contain duplicates")
    cfg["requested_observables"] = list(args.observables)
    active_proton_policy = proton_selection_metadata()
    if "proton_selection" in cfg and cfg["proton_selection"] != active_proton_policy:
        raise ValueError("Config proton_selection differs from YPSILON_PROTON_SELECTION; use a consistent explicit policy")
    cfg["proton_selection"] = active_proton_policy
    overrides = {}
    for item in args.input:
        state, separator, path = item.partition("=")
        if not separator or state not in args.states or not path or state in overrides:
            raise ValueError("--input needs one unique selected STATE=ROOT_PATH")
        overrides[state] = Path(path).expanduser().resolve()
    root = args.framework_root.expanduser().resolve()
    inputs = {}
    for state in args.states:
        scfg = cfg["states"][state]
        default_path = root / "results/generated/runs" / scfg["run"] / "artifacts/inputs/9x275_normalized.root"
        inputs[state] = overrides.get(state, default_path)
    return cfg, inputs


def select_common(optimizations):
    """Intersect passing candidates; never substitute a per-state winner silently."""
    lookup = {}
    for state, result in optimizations.items():
        lookup[state] = {tuple(row["edges"]): row for row in result.get("candidates", [])
                         if row.get("passes", False)}
    common = set.intersection(*(set(rows) for rows in lookup.values())) if lookup else set()
    if not common:
        return None
    def score(edges):
        rows = [lookup[state][edges] for state in lookup]
        return (-(len(edges) - 1),
                max(float(row.get("condition_number", float("inf"))) for row in rows),
                -min(float(row.get("count_balance", 0)) for row in rows), edges)
    return list(min(common, key=score))


def source_fingerprint(cfg, inputs):
    files = {}
    for state, path in inputs.items():
        if not path.is_file():
            files[state] = {"path": str(path), "missing": True}
            continue
        stat = path.stat()
        record = {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
        if cfg.get("hash_inputs"):
            h = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                    h.update(chunk)
            record["sha256"] = h.hexdigest()
        files[state] = record
    code = {}
    for path in sorted((FRAMEWORK / "scripts/reoptimization").glob("*.py")):
        code[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    code["run_reoptimization.py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    code["proton_selection.py"] = hashlib.sha256((FRAMEWORK / "workflow/proton_selection.py").read_bytes()).hexdigest()
    return {"config": cfg, "inputs": files, "code_sha256": code}


def run(args, cfg, inputs):
    run_id = args.run_id or "reopt-9x275-L{}-{}".format(
        str(cfg["target_lumi_invfb"]).replace(".", "p"), datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise ValueError("run-id may contain only letters, digits, dot, underscore or hyphen")
    output_candidate = args.out_dir if args.out_dir else (
        args.framework_root / "results/generated/runs" / run_id)
    enforce_category_workspace(args.framework_root, output_candidate)
    out = output_candidate.expanduser().resolve()
    artifacts = out / "artifacts/reoptimization"
    fingerprint = source_fingerprint(cfg, inputs)
    fingerprint["requested_states"] = list(args.states)
    fingerprint["requested_observables"] = list(args.observables)
    fp = digest(fingerprint)
    manifest_path = out / "run_manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if not args.resume:
            raise FileExistsError("Run exists. Use --resume with unchanged settings, or a fresh --run-id.")
        if old.get("fingerprint") != fp:
            raise ValueError("Resume fingerprint differs (input/config/code). Use a fresh run ID.")
    elif args.resume:
        raise FileNotFoundError("--resume requires an existing run_manifest.json")
    manifest = {"run_id": run_id, "status": "RUNNING", "fingerprint": fp,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "scope": cfg["uncertainty_scope"], "publication_ready": False,
                "fingerprint_details": fingerprint, "tasks": {}}
    write_json(manifest_path, manifest)
    write_json(out / "effective_config.json", cfg)
    samples, results, optimizations = {}, {}, {}
    for state in args.states:
        results[state] = {}
        try:
            print("LOAD", state, inputs[state], flush=True)
            sample = load_sample(inputs[state], cfg["states"][state], cfg)
            samples[state] = sample
            write_json(artifacts / state / "input_audit.json", sample["metadata"])
            print_input_audit(state, sample)
        except Exception as exc:
            traceback.print_exc()
            for observable in args.observables:
                results[state][observable] = {"status": "INPUT_ERROR", "error": str(exc),
                    "optimization": {"status": "BLOCKED_INPUT", "selected_edges": None},
                    "complete": False, "publication_ready": False}
            continue
        for observable in args.observables:
            try:
                unavailable = sample.get("metadata", {}).get("unavailable_observables", {})
                if observable in unavailable:
                    raise ValueError(unavailable[observable])
                print("OPTIMIZE", state, observable, flush=True)
                optimizations[state, observable] = optimize(sample, observable, cfg)
            except Exception as exc:
                traceback.print_exc()
                results[state][observable] = {"status": "OPTIMIZATION_ERROR", "error": str(exc),
                    "optimization": {"status": "OPTIMIZATION_ERROR", "selected_edges": None},
                    "complete": False, "publication_ready": False}
    if cfg.get("common_binning"):
        for observable in args.observables:
            all_present = all((state, observable) in optimizations for state in args.states)
            edges = select_common({state: optimizations[state, observable] for state in args.states}) if all_present else None
            for state in args.states:
                if (state, observable) in optimizations:
                    opt = optimizations[state, observable]
                    opt["individual_selected_edges"] = opt.get("selected_edges")
                    opt["selected_edges"] = edges
                    opt["common_binning"] = True
                    opt["status"] = "SELECTED" if edges else "NO_COMMON_FEASIBLE_BINNING"
    for state in args.states:
        for observable in args.observables:
            task = state + "/" + observable
            directory = artifacts / state / observable
            directory.mkdir(parents=True, exist_ok=True)
            result_path = directory / "result.json"
            if args.resume and result_path.exists():
                cached = json.loads(result_path.read_text())
                if cached.get("complete") and cached.get("fingerprint") == fp:
                    results[state][observable] = cached
                    manifest["tasks"][task] = {"status": cached["status"], "resumed": True}
                    print("REUSE", task, flush=True)
                    continue
            if state not in samples or (state, observable) not in optimizations:
                entry = results[state][observable]
                entry.update(state=state, observable=observable)
                propagate_status(entry)
                manifest["tasks"][task] = entry
                write_json(result_path, entry)
                continue
            sample = samples[state]
            opt = optimizations[state, observable]
            entry = {"state": state, "observable": observable, "optimization": opt,
                     "selection_metadata": sample["metadata"], "state_config": cfg["states"][state],
                     "proton_selection": sample["metadata"].get("proton_selection", {"mode": "unrecorded"}),
                     "cross_section_factor": (cross_section_factor(cfg, cfg["states"][state])
                                              if observable in ("pt2", "abs_t_via_pt2") else None),
                     "flux_integral": integrated_flux(cfg, cfg["states"][state]),
                     "fingerprint": fp, "publication_ready": False,
                     "scope": cfg["uncertainty_scope"]}
            results[state][observable] = entry
            try:
                edges = opt.get("selected_edges")
                if edges is None:
                    entry["status"] = opt.get("status", "NO_FEASIBLE_SLOPE_BINNING")
                    entry["complete"] = True
                    print(entry["status"], task, flush=True)
                else:
                    print("VALIDATE_RESPONSE_AND_INFERENCE", task, "edges", edges, flush=True)
                    model = build_model(sample, observable, edges, cfg, split="train")
                    validation = build_model(sample, observable, edges, cfg, split="test")
                    model["coarse"] = evaluate_edges(sample, observable, edges, cfg, split="train")
                    validation["coarse"] = evaluate_edges(sample, observable, edges, cfg, split="test")
                    entry["model"], entry["validation_model"] = model, validation
                    entry["cross_section_scaling"] = observable_cross_section_scaling(
                        cfg, cfg["states"][state], observable, edges)
                    inference = (analyze(model, validation, cfg) if observable in ("pt2", "abs_t_via_pt2")
                                 else analyze_spectrum(model, validation, cfg))
                    entry["inference"] = inference
                    entry["stress_tests"] = build_stress_tests(
                        sample, observable, edges, cfg, model, validation, inference)
                    if observable in ("pt2", "abs_t_via_pt2"):
                        fine_cfg = deepcopy(cfg)
                        fine_cfg["optimization"]["fine_truth_step"] *= 0.5
                        if "fine_truth_step" in fine_cfg["observables"].get(observable, {}):
                            fine_cfg["observables"][observable]["fine_truth_step"] *= 0.5
                        finer = build_model(sample, observable, edges, fine_cfg, split="train")
                        entry["fine_cell_check"] = {
                            "step_GeV2": finer["metadata"].get("max_fine_truth_width"),
                            "fit": fit_slope(finer, validation["reco_expected"], cfg),
                            "interpretation": "Response cell-size sensitivity, not independent data"}
                    else:
                        entry["fine_cell_check"] = {
                            "status": "NOT_APPLICABLE",
                            "interpretation": "Nonparametric coarse-bin spectrum inference uses no finite-cell exponential slope model"}
                    validation_gates(entry, cfg)
                    entry["complete"] = True
                    np.savez_compressed(directory / "response_matrices.npz", **{
                        split + "_" + key: value
                        for split, current in (("train", model), ("test", validation))
                        for key, value in current.items() if isinstance(value, np.ndarray)})
            except Exception as exc:
                traceback.print_exc()
                entry["status"] = "ANALYSIS_ERROR"
                entry["error"] = str(exc)
                entry["complete"] = False
            propagate_status(entry)
            write_json(result_path, entry)
            manifest["tasks"][task] = {"status": entry["status"], "result": str(result_path)}
            write_json(manifest_path, manifest)
    write_json(artifacts / "all_results.json", results)
    try:
        report_config = deepcopy(cfg)
        if any(sample.get("metadata", {}).get("normalization_metadata", {}).get("synthetic_test_input", False)
               for sample in samples.values()):
            report_config["synthetic"] = True
        write_report(artifacts, results, report_config)
    except Exception as exc:
        traceback.print_exc()
        manifest["report_error"] = str(exc)
    errors = any(not entry.get("complete", False) for per_state in results.values() for entry in per_state.values())
    manifest["status"] = "PARTIAL" if errors or "report_error" in manifest else "COMPLETE"
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["statistical_validation_statuses"] = {task: data["status"] for task, data in manifest["tasks"].items()}
    write_json(manifest_path, manifest)
    print("Run status:", manifest["status"], flush=True)
    print("Outputs:", artifacts, flush=True)
    print("A COMPLETE execution is not a claim of validated full experimental uncertainties.", flush=True)
    return 1 if manifest["status"] == "PARTIAL" else 0


def validation_gates(entry, cfg):
    # Legacy dictionaries remain interpretable for archived diagnostics. New
    # executions always produce the three explicitly separated stress families.
    legacy = any(case.get("kind") in ("truth_shape", "feed_in_normalization_sensitivity")
                 for case in entry.get("stress_tests", []))
    if legacy or "observable" not in entry:
        _legacy_validation_gates(entry, cfg)
        propagate_status(entry)
    else:
        separated_validation_gates(entry, cfg)


def _legacy_validation_gates(entry, cfg):
    """Gate sensitivity checks without mistaking an expected slope shift for bias."""
    inf = entry["inference"]
    nominal = inf.get("slope", {})
    b = nominal.get("b")
    low, high = nominal.get("interval_low"), nominal.get("interval_high")
    width = (float(high) - float(low)) / 2 if low is not None and high is not None else None
    usable = b is not None and width is not None and np.isfinite(width) and width > 0
    finer = entry["fine_cell_check"]["fit"]
    refinement_ok = usable and finer.get("status") == "OK_CONDITIONAL" and finer.get("b") is not None
    fraction = abs(finer["b"] - b) / width if refinement_ok else None
    threshold = float(cfg.get("validation", {}).get("max_fine_cell_shift_fraction", 0.2))
    entry["fine_cell_check"].update(shift_over_nominal_halfwidth=fraction,
        threshold=threshold, passes=bool(refinement_ok and fraction <= threshold))
    shape_checks, shape_fit_checks, feed_checks = [], [], []
    for stress in entry["stress_tests"]:
        if stress["kind"] == "truth_shape":
            closure = stress["coarse_closure"]
            bias = np.asarray(closure.get("bias_over_nominal_sampling_sd", []), float)
            tolerance = float(cfg.get("inference", {}).get("max_abs_standardized_bias", 0.3))
            passed = closure.get("status") == "OK_CONDITIONAL" and len(bias) > 0 and np.isfinite(bias).all() and np.all(np.abs(bias) <= tolerance)
            stress["closure_passes"] = bool(passed)
            stress["closure_tolerance_fraction_of_nominal_sampling_sd"] = tolerance
            shape_checks.append(bool(passed))
            fitted = stress.get("fit", {})
            shape_fit_checks.append(fitted.get("status") == "OK_CONDITIONAL" and fitted.get("interval_high") is not None)
        elif stress["kind"] == "feed_in_normalization_sensitivity":
            fitted = stress["fit"]
            ok = usable and fitted.get("status") == "OK_CONDITIONAL" and fitted.get("b") is not None
            shift = abs(fitted["b"] - b) / width if ok else None
            tolerance = float(cfg.get("validation", {}).get("max_feed_in_shift_fraction", 0.5))
            stress["shift_over_nominal_halfwidth"] = shift
            stress["sensitivity_flag"] = not (ok and shift <= tolerance)
            feed_checks.append(not stress["sensitivity_flag"])
    gates = {"fine_cell_refinement": entry["fine_cell_check"]["passes"],
             "truth_shape_stress_closure": bool(shape_checks) and all(shape_checks),
             "truth_shape_slope_identifiable": bool(shape_fit_checks) and all(shape_fit_checks),
             "feed_in_sensitivity": bool(feed_checks) and all(feed_checks)}
    entry["additional_validation_gates"] = gates
    entry["status"] = ("PASS_CONDITIONAL" if inf.get("status") == "PASS_CONDITIONAL" and all(gates.values())
                       else "NOT_VALIDATED")


def main(argv=None):
    args = parse_args(argv)
    cfg, inputs = prepare_config(args)
    if args.command == "plan":
        print("9x275 target luminosity:", cfg["target_lumi_invfb"], "fb^-1")
        print("Observables:", ", ".join(args.observables))
        print("Input ROOT files are read-only. Existing outputs are not changed.")
        for state, path in inputs.items():
            print(state, "FOUND" if path.is_file() else "MISSING", path)
        print("Steps: normalization validation, training-only bin scan, heldout Poisson inference, bin/slope coverage, response-MC variation, separate shape/tail/combined stresses, appendix reports.")
        print("W reports flux-bin-averaged sigma_gamma_p. Rapidity reports inclusive ep dSigma/dy and a labeled total-flux-normalized companion.")
        return 0
    if args.command == "validate":
        failures = 0
        for state, path in inputs.items():
            try:
                sample = load_sample(path, cfg["states"][state], cfg)
                print("PASS", state, "rows", len(sample["weight"]),
                      "source luminosity", sample["source_lumi"], "flux", integrated_flux(cfg, cfg["states"][state]))
                print_input_audit(state, sample)
                unavailable = sample.get("metadata", {}).get("unavailable_observables", {})
                for observable in args.observables:
                    if observable in unavailable:
                        print("FAIL", state + "/" + observable, unavailable[observable], file=sys.stderr)
                        failures += 1
                    else:
                        print("READY", state + "/" + observable)
            except Exception as exc:
                print("FAIL", state, str(exc), file=sys.stderr)
                failures += 1
        return int(failures > 0)
    return run(args, cfg, inputs)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, FileNotFoundError, FileExistsError, KeyError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        sys.exit(2)
