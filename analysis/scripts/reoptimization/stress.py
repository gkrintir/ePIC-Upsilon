"""Separate in-range shape, outside-range yield, and combined robustness tests.

These variations diagnose a fixed calibration model. They are not assigned
systematic uncertainties and never automatically upgrade a failed case.
"""

import numpy as np

from .response import build_model
from .inference import fit_slope, poisson_unfold


SLOPE_OBSERVABLES = ("pt2", "abs_t_via_pt2")
STRESS_MODES = ("shape_only", "tail_only", "combined")


def closure_result(model, counts, reference, sampling_sd):
    reference = np.asarray(reference, float)
    result = poisson_unfold(model["response"], counts, model["feed_in"], reference)
    result["reference_truth"] = reference.tolist()
    if result.get("truth") is not None:
        difference = np.asarray(result["truth"]) - reference
        result["relative_bias"] = np.divide(difference, reference,
            out=np.full_like(difference, np.nan), where=reference > 0).tolist()
        sd = np.asarray(sampling_sd, float)
        if sd.shape == difference.shape:
            result["bias_over_nominal_sampling_sd"] = np.divide(difference, sd,
                out=np.full_like(difference, np.nan), where=sd > 0).tolist()
    return result


def build_stress_tests(sample, observable, edges, config, model, validation, inference):
    """Freeze training calibration and vary only independent held-out means."""
    modes = tuple(config.get("stress_test_modes", STRESS_MODES))
    if set(modes) != set(STRESS_MODES) or len(modes) != len(STRESS_MODES):
        raise ValueError("stress_test_modes must contain shape_only, tail_only and combined exactly once")
    tilts = config.get("shape_stress_tilts", [-1.0, 1.0])
    scales = config.get("feed_in_stress_scales", [.5, 1.5])
    if not tilts or not scales:
        raise ValueError("Both shape tilts and tail scales are required for robustness validation")
    if any(not np.isfinite(v) or v == 0 for v in tilts):
        raise ValueError("Shape tilts must be finite and nonzero")
    if any(not np.isfinite(v) or v < 0 or v == 1 for v in scales):
        raise ValueError("Tail scales must be finite, nonnegative, and different from one")
    sd = inference.get("coarse", {}).get("sampling_sd", [])
    stresses = []
    nominal_feed = np.asarray(validation["feed_in"], float)
    for mode in modes:
        if mode == "tail_only":
            for scale in scales:
                means = (np.asarray(validation["reco_expected"], float)
                         + (float(scale) - 1.) * nominal_feed)
                reference = validation["coarse"]["truth_expected"]
                case = {"kind": mode, "scale": float(scale),
                        "normalization": {
                            "in_range_truth_unchanged": True,
                            "feed_in_expected_before": float(nominal_feed.sum()),
                            "feed_in_expected_after": float(scale * nominal_feed.sum()),
                            "feed_in_scale": float(scale)},
                        "applicability": ("APPLICABLE" if nominal_feed.sum() > 0 else
                                          "NO_FEED_IN_WITHIN_CONFIGURED_TRUTH_FIDUCIAL"),
                        "truth_expected": validation["truth_expected"],
                        "reco_expected": means,
                        "interpretation": "Outside-variable-range feed-in yield varied; in-range truth and its detector response unchanged"}
                case["coarse_closure"] = closure_result(model["coarse"], means, reference, sd)
                if observable in SLOPE_OBSERVABLES:
                    case["fit"] = fit_slope(model, means, config)
                stresses.append(case)
        else:
            for tilt in tilts:
                shifted = build_model(sample, observable, edges, config, split="test",
                                      shape_tilt=float(tilt), shape_mode=mode)
                means = shifted["reco_expected"]
                case = {"kind": mode, "normalization": shifted["metadata"],
                        "truth_expected": shifted["truth_expected"], "reco_expected": means,
                        "interpretation": ("In-range shape varied with its integral preserved and outside-range weights unchanged"
                                           if mode == "shape_only" else
                                           "Clipped global truth-shape tilt with all-selected-truth normalization preserved; outside-range feed-in also changes")}
                if observable in SLOPE_OBSERVABLES:
                    case.update(tilt_GeVminus2=float(tilt), tilt_units="GeV^-2")
                else:
                    case.update(tilt_dimensionless=float(tilt), tilt_units="dimensionless",
                                tilt_coordinate="u=(x-fit_min)/(fit_max-fit_min)")
                case["coarse_closure"] = closure_result(model["coarse"], means,
                    shifted["coarse"]["truth_expected"], sd)
                if observable in SLOPE_OBSERVABLES:
                    case["fit"] = fit_slope(model, means, config)
                stresses.append(case)
    return stresses


def validation_gates(entry, config):
    """Keep component decisions and every stress family visible and conservative."""
    inference = entry["inference"]
    slope_case = entry["observable"] in SLOPE_OBSERVABLES
    closure_tolerance = float(config.get("inference", {}).get("max_abs_standardized_bias", .3))
    groups = {mode: [] for mode in STRESS_MODES}
    fit_groups = {mode: [] for mode in STRESS_MODES}
    nominal = inference.get("slope", {})
    low, high = nominal.get("interval_low"), nominal.get("interval_high")
    halfwidth = (float(high) - float(low)) / 2 if low is not None and high is not None else None
    usable_slope = (nominal.get("b") is not None and halfwidth is not None
                    and np.isfinite(halfwidth) and halfwidth > 0)
    for case in entry.get("stress_tests", []):
        kind = case["kind"]
        if kind not in groups:
            raise ValueError("Unrecognized stress kind: " + str(kind))
        closure = case.get("coarse_closure", {})
        bias = np.asarray(closure.get("bias_over_nominal_sampling_sd", []), float)
        ok = (closure.get("status") == "OK_CONDITIONAL" and len(bias) > 0
              and np.isfinite(bias).all() and np.all(np.abs(bias) <= closure_tolerance))
        case["closure_passes"] = bool(ok)
        case["closure_tolerance_fraction_of_nominal_sampling_sd"] = closure_tolerance
        case["stage_status"] = "PASS_CONDITIONAL" if ok else "NOT_VALIDATED"
        groups[kind].append(bool(ok))
        if slope_case:
            fit = case.get("fit", {})
            identified = (fit.get("status") == "OK_CONDITIONAL"
                          and fit.get("interval_high") is not None)
            if kind == "tail_only":
                shift = abs(fit["b"] - nominal["b"]) / halfwidth if identified and usable_slope else None
                threshold = float(config.get("validation", {}).get("max_feed_in_shift_fraction", .5))
                case.update(shift_over_nominal_halfwidth=shift,
                            sensitivity_flag=not (shift is not None and shift <= threshold))
                identified = identified and not case["sensitivity_flag"]
            case["slope_check_passes"] = bool(identified)
            fit_groups[kind].append(bool(identified))
    gates = {mode + "_closure": bool(groups[mode]) and all(groups[mode]) for mode in STRESS_MODES}
    slope_gates = {}
    if slope_case:
        finer = entry.get("fine_cell_check", {}).get("fit", {})
        ok = usable_slope and finer.get("status") == "OK_CONDITIONAL" and finer.get("b") is not None
        shift = abs(finer["b"] - nominal["b"]) / halfwidth if ok else None
        threshold = float(config.get("validation", {}).get("max_fine_cell_shift_fraction", .2))
        refinement_pass = bool(ok and shift <= threshold)
        entry["fine_cell_check"].update(shift_over_nominal_halfwidth=shift,
                                       threshold=threshold, passes=refinement_pass)
        slope_gates = {mode + "_slope_check": bool(fit_groups[mode]) and all(fit_groups[mode])
                       for mode in STRESS_MODES}
        slope_gates["fine_cell_refinement"] = refinement_pass
        gates.update(slope_gates)
    entry["additional_validation_gates"] = gates
    inference_gates = inference.get("gates", {})
    failures = ["inference." + key for key, ok in inference_gates.items() if not ok]
    if inference.get("status") != "PASS_CONDITIONAL" and not failures:
        failures.append("inference.status")
    failures += ["robustness." + key for key, ok in gates.items() if not ok]
    entry["failed_checks"] = failures
    entry["status"] = "PASS_CONDITIONAL" if not failures else "NOT_VALIDATED"
    differential_ok = (inference.get("coarse", {}).get("status") == "PASS_CONDITIONAL"
                       and all(bool(groups[k]) and all(groups[k]) for k in STRESS_MODES))
    if not slope_case:
        differential_ok = differential_ok and inference.get("status") == "PASS_CONDITIONAL"
    components = {"differential": "PASS_CONDITIONAL" if differential_ok else "NOT_VALIDATED"}
    if slope_case:
        # The case-wide outcome still requires coarse closure. This separate
        # component does not silently promote a failed differential extraction.
        base_slope_gates = {k: v for k, v in inference_gates.items() if k != "coarse_extraction"}
        slope_ok = bool(base_slope_gates) and all(base_slope_gates.values()) and all(slope_gates.values())
        components["slope"] = "PASS_CONDITIONAL" if slope_ok else "NOT_VALIDATED"
    else:
        components["bin_intervals"] = inference.get("status", "NOT_VALIDATED")
    entry["component_validation"] = components
    propagate_status(entry)


def propagate_status(entry):
    """Keep historical substage statuses explicit beside the final case status."""
    overall = entry.get("status", "NOT_VALIDATED")
    entry["overall_status"] = overall
    for value in entry.get("inference", {}).values():
        if isinstance(value, dict) and "status" in value:
            value["stage_status"] = value["status"]
            value["overall_status"] = overall
        elif isinstance(value, list):
            for row in value:
                if isinstance(row, dict):
                    if "status" in row:
                        row["stage_status"] = row["status"]
                    row["overall_status"] = overall
    if "inference" in entry:
        entry["inference"]["stage_status"] = entry["inference"].get("status", "NOT_RUN")
        entry["inference"]["overall_status"] = overall
    for case in entry.get("stress_tests", []):
        case["overall_status"] = overall
        for key in ("fit", "coarse_closure"):
            if key in case:
                case[key]["stage_status"] = case[key].get("status", "NOT_RUN")
                case[key]["overall_status"] = overall
