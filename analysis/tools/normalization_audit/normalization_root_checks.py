#!/usr/bin/env python3
"""Read-only ROOT inspection, ordered tuple comparison, and beam-formula checks.

No generated-event selection or normalization is performed here. Histogram
integrals are producer diagnostics, not inferred physical rates. Formula
agreement identifies stored beam constants, not original generator provenance.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import uproot

WEIGHTS = {"weight", "weight_proj"}
BEAM_FIELDS = ("true_W", "true_Q2", "true_scat_eMinus_pT", "true_scat_eMinus_eta", "true_scat_eMinus_E")
MAX_EXAMPLES = 10
MAX_NESTED_VALUES_PER_CHUNK = 1000000


def sha256(path):
    """Hash a mounted/local file without changing it."""
    result = hashlib.sha256()
    with Path(path).expanduser().open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def _path(path):
    return str(Path(path).expanduser().resolve())


def _float(value):
    value = float(value)
    return value if np.isfinite(value) else None


def _error(exc):
    return type(exc).__name__ + ": " + str(exc)


def _metadata_text(obj, classname):
    return str(obj) if classname == "TObjString" else str(obj.member("fTitle"))


def inspect_root(path, tree="UpsilonTree", hash_input=False):
    """Inspect file metadata and histogram bookkeeping without reading event arrays."""
    result = {"status": "OK", "path": _path(path), "tree": tree, "tree_class": None,
              "entries": None, "branches": {}, "objects": {}, "histograms": {},
              "framework_normalization": None, "string_metadata": {}, "errors": [],
              "histogram_interpretation": "Stored producer diagnostics only; integrals are not identified as physical rates"}
    try:
        with uproot.open(result["path"]) as source:
            result["objects"] = dict(source.classnames(recursive=True, cycle=False))
            if tree not in source:
                result["errors"].append("Missing tree " + tree)
            else:
                events = source[tree]
                result["tree_class"] = events.classname
                result["entries"] = int(events.num_entries)
                for name in events.keys():
                    branch = events[name]
                    interpretation = getattr(branch, "interpretation", None)
                    result["branches"][name] = {"typename": str(getattr(branch, "typename", "unknown")),
                                                "interpretation": str(interpretation) if interpretation is not None else None}
            for name, classname in result["objects"].items():
                if classname.startswith(("TH1", "TH2", "TH3")):
                    item = {"class": classname, "status": "OK", "diagnostic_only": True}
                    try:
                        obj = source[name]
                        shape = [len(axis.edges()) + 1 for axis in obj.axes]
                        item["bins_with_flow"] = shape
                        if int(np.prod(shape)) > 2000000:
                            item.update(status="SKIP", reason="Histogram exceeds two million bins; bounded inspection")
                        else:
                            values = np.asarray(obj.values(flow=True), dtype=float)
                            core = values[tuple(slice(1, -1) for _ in range(values.ndim))]
                            item.update(title=str(obj.title), sum_without_flow=_float(np.sum(core)),
                                        sum_with_flow=_float(np.sum(values)),
                                        flow_sum=_float(np.sum(values) - np.sum(core)),
                                        nonfinite_bins=int(np.count_nonzero(~np.isfinite(values))))
                            if values.ndim == 1:
                                item["underflow"] = _float(values[0])
                                item["overflow"] = _float(values[-1])
                                if core.size <= 200:
                                    item["bin_contents"] = [_float(value) for value in core]
                                    labels = obj.axis().labels()
                                    item["bin_labels"] = list(labels) if labels is not None else None
                    except Exception as exc:
                        item.update(status="ERROR", error=_error(exc))
                    result["histograms"][name] = item
                elif classname in ("TObjString", "TNamed"):
                    try:
                        text = _metadata_text(source[name], classname)
                        if len(text) > 1048576:
                            value = {"status": "SKIP", "reason": "String metadata exceeds one MiB", "characters": len(text)}
                        else:
                            try:
                                def reject_nonfinite(token):
                                    raise ValueError("Nonfinite JSON constant " + token)
                                value = json.loads(text, parse_constant=reject_nonfinite)
                            except (ValueError, TypeError):
                                value = text
                        result["string_metadata"][name] = value
                        if name == "framework_normalization":
                            result["framework_normalization"] = value
                    except Exception as exc:
                        result["errors"].append("Cannot read string metadata " + name + ": " + _error(exc))
        if hash_input:
            result["sha256"] = sha256(result["path"])
    except Exception as exc:
        result["errors"].append(_error(exc))
    if result["errors"]:
        result["status"] = "ERROR"
    return result


def _step(step_size):
    value = int(step_size)
    if value < 1:
        raise ValueError("step_size must be a positive number of rows")
    return value


def _chunks(rows, step_size):
    for start in range(0, rows, _step(step_size)):
        yield start, min(rows, start + _step(step_size))


def _read(branch, start, stop):
    try:
        return branch.array(entry_start=start, entry_stop=stop, library="np")
    except Exception:
        import awkward as ak
        return ak.to_list(branch.array(entry_start=start, entry_stop=stop, library="ak"))


def _nested_equal(a, b, rtol, atol, budget):
    """Compare one bounded nested entry; return exact, close, max finite delta."""
    budget[0] -= 1
    if budget[0] < 0:
        raise ValueError("Nested entry comparison exceeds bounded values per chunk")
    if isinstance(a, (np.ndarray, np.void)):
        a = a.tolist()
    if isinstance(b, (np.ndarray, np.void)):
        b = b.tolist()
    if isinstance(a, dict) or isinstance(b, dict):
        if not isinstance(a, dict) or not isinstance(b, dict) or set(a) != set(b):
            return False, False, None
        comparisons = [_nested_equal(a[key], b[key], rtol, atol, budget) for key in a]
    elif isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        if not isinstance(a, (list, tuple)) or not isinstance(b, (list, tuple)) or len(a) != len(b):
            return False, False, None
        comparisons = [_nested_equal(x, y, rtol, atol, budget) for x, y in zip(a, b)]
    else:
        if a is None or b is None:
            return a is b, a is b, None
        if isinstance(a, (float, int, complex, np.number, bool)) and isinstance(b, (float, int, complex, np.number, bool)):
            exact = bool(a == b or (np.isnan(a) and np.isnan(b)))
            close = bool(np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True))
            delta = _float(abs(complex(a) - complex(b))) if np.isfinite(a) and np.isfinite(b) else None
            return exact, close, delta
        equal = bool(a == b)
        return equal, equal, None
    deltas = [item[2] for item in comparisons if item[2] is not None]
    return all(item[0] for item in comparisons), all(item[1] for item in comparisons), max(deltas) if deltas else None


def _compare_arrays(a, b, rtol, atol):
    """Return entry-level exact/tolerance masks and the maximum finite difference."""
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray) and a.dtype.kind in "biufc" and b.dtype.kind in "biufc":
        if a.shape != b.shape:
            return np.zeros(len(a), dtype=bool), np.zeros(len(a), dtype=bool), None
        exact = (a == b) | (np.isnan(a) & np.isnan(b))
        close = np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
        valid = np.isfinite(a) & np.isfinite(b)
        with np.errstate(over="ignore", invalid="ignore"):
            delta = np.abs(a.astype(np.clongdouble if a.dtype.kind == "c" or b.dtype.kind == "c" else np.longdouble) - b)
        maximum = _float(np.max(delta[valid])) if valid.any() else None
        if a.ndim > 1:
            exact, close = exact.reshape(len(a), -1).all(axis=1), close.reshape(len(a), -1).all(axis=1)
        return exact, close, maximum
    if len(a) != len(b):
        return np.zeros(len(a), dtype=bool), np.zeros(len(a), dtype=bool), None
    budget = [MAX_NESTED_VALUES_PER_CHUNK]
    comparisons = [_nested_equal(x, y, rtol, atol, budget) for x, y in zip(a, b)]
    deltas = [item[2] for item in comparisons if item[2] is not None]
    return (np.array([item[0] for item in comparisons], dtype=bool),
            np.array([item[1] for item in comparisons], dtype=bool), max(deltas) if deltas else None)


def _branch_pair(first, second, rows, step_size, rtol, atol):
    first_type, second_type = str(getattr(first, "typename", "unknown")), str(getattr(second, "typename", "unknown"))
    result = {"status": "PASS", "compared_entries": 0, "changed_entries": 0,
              "outside_tolerance_entries": 0, "max_absolute_difference": None,
              "first_changed_indices": [], "first_outside_tolerance_indices": [],
              "merged_type": first_type, "normalized_type": second_type, "type_equal": first_type == second_type}
    try:
        for start, stop in _chunks(rows, step_size):
            a, b = _read(first, start, stop), _read(second, start, stop)
            exact, close, maximum = _compare_arrays(a, b, rtol, atol)
            result["compared_entries"] += stop - start
            result["changed_entries"] += int(np.count_nonzero(~exact))
            result["outside_tolerance_entries"] += int(np.count_nonzero(~close))
            if maximum is not None:
                prior = result["max_absolute_difference"]
                result["max_absolute_difference"] = maximum if prior is None else max(prior, maximum)
            for key, mask in (("first_changed_indices", ~exact), ("first_outside_tolerance_indices", ~close)):
                remaining = MAX_EXAMPLES - len(result[key])
                result[key].extend((np.flatnonzero(mask)[:remaining] + start).tolist())
        if result["changed_entries"] or not result["type_equal"]:
            result["status"] = "CHANGED"
    except Exception as exc:
        result.update(status="SKIP", reason=_error(exc))
    return result


def _summary():
    return {"count": 0, "sum": 0., "sum_squares": 0., "min": None, "max": None}


def _accumulate(summary, values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return
    summary["count"] += int(values.size)
    summary["sum"] += float(np.sum(values, dtype=np.longdouble))
    summary["sum_squares"] += float(np.sum(values.astype(np.longdouble) ** 2))
    summary["min"] = float(np.min(values)) if summary["min"] is None else min(summary["min"], float(np.min(values)))
    summary["max"] = float(np.max(values)) if summary["max"] is None else max(summary["max"], float(np.max(values)))


def _finish(summary):
    n = summary["count"]
    summary["mean"] = _float(summary["sum"] / n) if n else None
    summary["rms"] = _float(np.sqrt(summary["sum_squares"] / n)) if n else None
    for key in ("sum", "sum_squares"):
        summary[key] = _float(summary[key])
    return summary


def _ratio(numerator, denominator):
    numerator, denominator = np.asarray(numerator, dtype=float), np.asarray(denominator, dtype=float)
    if numerator.ndim != 1 or denominator.ndim != 1:
        raise ValueError("Weight branches must be scalar per event")
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ratio = numerator / denominator
    categories = np.zeros(len(ratio), dtype=np.int8)
    categories[~np.isfinite(ratio)] = 1
    categories[denominator == 0] = 2
    categories[~np.isfinite(numerator)] = 3
    categories[~np.isfinite(denominator)] = 4
    return ratio, categories


def _weights_pair(first, second, rows, step_size, rtol, atol):
    ratio_report = {"status": "OK", "definition": "normalized weight / merged weight", "values": _summary(),
                    "invalid_ratio_rows": 0, "raw_zero_weight_rows": 0, "normalized_zero_weight_rows": 0,
                    "raw_nonfinite_weight_rows": 0, "normalized_nonfinite_weight_rows": 0,
                    "first_invalid_indices": []}
    projective = {"status": "PASS", "definition": "weight_proj/weight compared event by event",
                  "valid_ratio_rows": 0, "undefined_ratio_rows": 0, "changed_entries": 0,
                  "exact_changed_entries": 0, "changed_invalid_classification_rows": 0,
                  "changed_undefined_input_rows": 0, "first_changed_indices": [],
                  "merged_values": _summary(), "normalized_values": _summary(),
                  "max_absolute_difference": None}
    if "weight" not in first or "weight" not in second:
        reason = "Missing weight branch in one or both trees"
        return {"status": "SKIP", "reason": reason}, {"status": "SKIP", "reason": reason}
    has_projective = "weight_proj" in first and "weight_proj" in second
    if not has_projective:
        projective = {"status": "SKIP", "reason": "Missing weight_proj branch in one or both trees"}
    try:
        for start, stop in _chunks(rows, step_size):
            a = np.asarray(_read(first["weight"], start, stop), dtype=float)
            b = np.asarray(_read(second["weight"], start, stop), dtype=float)
            ratio, category = _ratio(b, a)
            valid = category == 0
            _accumulate(ratio_report["values"], ratio[valid])
            ratio_report["invalid_ratio_rows"] += int(np.count_nonzero(~valid))
            ratio_report["raw_zero_weight_rows"] += int(np.count_nonzero(a == 0))
            ratio_report["normalized_zero_weight_rows"] += int(np.count_nonzero(b == 0))
            ratio_report["raw_nonfinite_weight_rows"] += int(np.count_nonzero(~np.isfinite(a)))
            ratio_report["normalized_nonfinite_weight_rows"] += int(np.count_nonzero(~np.isfinite(b)))
            remaining = MAX_EXAMPLES - len(ratio_report["first_invalid_indices"])
            ratio_report["first_invalid_indices"].extend((np.flatnonzero(~valid)[:remaining] + start).tolist())
            if not has_projective:
                continue
            pa = np.asarray(_read(first["weight_proj"], start, stop), dtype=float)
            pb = np.asarray(_read(second["weight_proj"], start, stop), dtype=float)
            fa, ca = _ratio(pa, a)
            fb, cb = _ratio(pb, b)
            both = (ca == 0) & (cb == 0)
            undefined = ~both
            exact = (fa == fb) | (np.isnan(fa) & np.isnan(fb))
            close = np.isclose(fa, fb, rtol=rtol, atol=atol, equal_nan=True)
            category_changed = ca != cb
            # Undefined quotients cannot establish invariance. Record changes to
            # their underlying numerator/denominator separately, including NaNs.
            raw_inputs_equal = ((pa == pb) | (np.isnan(pa) & np.isnan(pb))) & ((a == b) | (np.isnan(a) & np.isnan(b)))
            undefined_changed = undefined & ~raw_inputs_equal
            changed = (both & ~close) | category_changed | undefined_changed
            projective["valid_ratio_rows"] += int(np.count_nonzero(both))
            projective["undefined_ratio_rows"] += int(np.count_nonzero(undefined))
            projective["changed_entries"] += int(np.count_nonzero(changed))
            projective["exact_changed_entries"] += int(np.count_nonzero(both & ~exact))
            projective["changed_invalid_classification_rows"] += int(np.count_nonzero(category_changed))
            projective["changed_undefined_input_rows"] += int(np.count_nonzero(undefined_changed))
            remaining = MAX_EXAMPLES - len(projective["first_changed_indices"])
            projective["first_changed_indices"].extend((np.flatnonzero(changed)[:remaining] + start).tolist())
            _accumulate(projective["merged_values"], fa[ca == 0])
            _accumulate(projective["normalized_values"], fb[cb == 0])
            if both.any():
                maximum = float(np.max(np.abs(fa[both] - fb[both])))
                previous = projective["max_absolute_difference"]
                projective["max_absolute_difference"] = maximum if previous is None else max(previous, maximum)
        _finish(ratio_report["values"])
        if ratio_report["invalid_ratio_rows"]:
            ratio_report["status"] = "PARTIAL"
        if has_projective:
            _finish(projective["merged_values"])
            _finish(projective["normalized_values"])
            projective["status"] = "CHANGED" if projective["changed_entries"] else "PARTIAL" if projective["undefined_ratio_rows"] else "PASS"
    except Exception as exc:
        ratio_report.update(status="ERROR", error=_error(exc))
        if has_projective:
            projective.update(status="ERROR", error=_error(exc))
    return ratio_report, projective


def compare_root_pairs(merged, normalized, tree="UpsilonTree", step_size=100000, rtol=1e-10, atol=1e-12):
    """Compare ordered common branches, excluding deliberate weight replacements.

    Exact differences and tolerance violations are reported separately. Numeric
    nonweight branches must be exactly unchanged for PASS, including paired NaNs.
    Nested/string branches are compared within a bounded per-chunk work budget.
    """
    result = {"status": "PASS", "merged": _path(merged), "normalized": _path(normalized),
              "tree": tree, "entries": {}, "missing_from_normalized": [], "missing_from_merged": [],
              "branches": {}, "errors": [], "rtol": rtol, "atol": atol,
              "comparison_scope": "Ordered common rows; nonweight exact comparison, weight-ratio diagnostics"}
    try:
        step_size = _step(step_size)
        if rtol < 0 or atol < 0 or not np.isfinite(rtol) or not np.isfinite(atol):
            raise ValueError("Tolerances must be finite and nonnegative")
        with uproot.open(result["merged"]) as raw, uproot.open(result["normalized"]) as prepared:
            first, second = raw[tree], prepared[tree]
            n_first, n_second = int(first.num_entries), int(second.num_entries)
            result["entries"] = {"merged": n_first, "normalized": n_second}
            result["row_count_equal"] = n_first == n_second
            result["compared_entries"] = min(n_first, n_second)
            first_names, second_names = set(first.keys()), set(second.keys())
            result["missing_from_normalized"] = sorted(first_names - second_names)
            result["missing_from_merged"] = sorted(second_names - first_names)
            for name in sorted((first_names & second_names) - WEIGHTS):
                result["branches"][name] = _branch_pair(first[name], second[name], min(n_first, n_second), step_size, rtol, atol)
            result["weight_ratio"], result["projective_factor_invariance"] = _weights_pair(
                first, second, min(n_first, n_second), step_size, rtol, atol)
            statuses = [item["status"] for item in result["branches"].values()]
            projective_status = result["projective_factor_invariance"]["status"]
            observed_branch_change = any(item["changed_entries"] or not item["type_equal"] for item in result["branches"].values())
            if n_first != n_second or first_names != second_names or observed_branch_change or projective_status == "CHANGED":
                result["status"] = "CHANGED"
            elif any(status != "PASS" for status in statuses) or projective_status != "PASS" or result["weight_ratio"]["status"] != "OK":
                result["status"] = "PARTIAL"
    except Exception as exc:
        result["errors"].append(_error(exc))
        result["status"] = "ERROR"
    return result


def _beam_setting(settings, name):
    for key in (name, name + "_GeV"):
        if key in settings:
            value = float(settings[key])
            if not np.isfinite(value) or value <= 0:
                raise ValueError("Beam constants must be finite and positive")
            return value
    raise ValueError("Missing beam setting " + name)


def beam_formula_check(path, settings, tree="UpsilonTree", step_size=100000):
    """Compare stored truth W/Q2 with requested and standard beam formulas.

    settings accepts electron_energy and proton_energy (or *_GeV aliases).
    The latter is the proton momentum convention of the supplied producer.
    Missing or invalid kinematics are reported, without event selection or
    attribution of actual generator beam energies.
    """
    result = {"status": "OK", "path": _path(path), "tree": tree, "checked_rows": 0,
              "unavailable_rows": 0, "missing_branches": [], "candidates": {}, "errors": [],
              "scope": "Stored formula agreement only; actual generator beam/rate/phase-space provenance unverified",
              "absolute_tolerance": 0.001, "relative_tolerance": 0.001,
              "pz_method": "sign(eta)*sqrt(E^2-pT^2-m_e^2), m_e=0.000511 GeV"}
    try:
        step_size = _step(step_size)
        electron = _beam_setting(settings, "electron_energy")
        proton = _beam_setting(settings, "proton_energy")
        requested = "{:g}x{:g}".format(electron, proton)
        result["requested_candidate"] = requested
        candidates = {(9., 130.), (9., 275.), (10., 100.), (18., 275.), (electron, proton)}
        for e_beam, p_beam in sorted(candidates):
            result["candidates"]["{:g}x{:g}".format(e_beam, p_beam)] = {
                "electron_energy_GeV": e_beam, "proton_momentum_GeV": p_beam,
                "true_W": {"matching_rows": 0, "absolute_differences": _summary()},
                "true_Q2": {"matching_rows": 0, "absolute_differences": _summary()}}
        with uproot.open(result["path"]) as source:
            events = source[tree]
            result["entries"] = int(events.num_entries)
            result["missing_branches"] = sorted(set(BEAM_FIELDS) - set(events.keys()))
            if result["missing_branches"]:
                result.update(status="SKIP", reason="Missing stored truth branches for optional beam-formula check")
                return result
            for start, stop in _chunks(int(events.num_entries), step_size):
                arrays = {name: np.asarray(_read(events[name], start, stop), dtype=float) for name in BEAM_FIELDS}
                if any(values.ndim != 1 for values in arrays.values()):
                    raise ValueError("Beam check requires scalar truth branches")
                pt, eta, energy = (arrays[name] for name in BEAM_FIELDS[2:])
                usable = np.logical_and.reduce([np.isfinite(arrays[name]) & (arrays[name] != -999) for name in BEAM_FIELDS])
                with np.errstate(over="ignore", invalid="ignore"):
                    pz2 = energy * energy - pt * pt - 0.000511 ** 2
                usable &= (pt >= 0) & (energy > 0) & np.isfinite(pz2) & (pz2 >= 0)
                result["checked_rows"] += int(np.count_nonzero(usable))
                result["unavailable_rows"] += int(np.count_nonzero(~usable))
                if not usable.any():
                    continue
                pt, energy = pt[usable], energy[usable]
                pz = np.sign(eta[usable]) * np.sqrt(pz2[usable])
                for candidate in result["candidates"].values():
                    e_beam, p_beam = candidate["electron_energy_GeV"], candidate["proton_momentum_GeV"]
                    q0, qz = e_beam - energy, -e_beam - pz
                    q2 = pt ** 2 + qz ** 2 - q0 ** 2
                    w2 = (np.sqrt(p_beam ** 2 + 0.938272 ** 2) + q0) ** 2 - pt ** 2 - (p_beam + qz) ** 2
                    for branch, calculated in (("true_W", np.sign(w2) * np.sqrt(np.abs(w2))), ("true_Q2", q2)):
                        difference = np.abs(calculated - arrays[branch][usable])
                        item = candidate[branch]
                        item["matching_rows"] += int(np.count_nonzero(np.isclose(calculated, arrays[branch][usable], rtol=0.001, atol=0.001)))
                        _accumulate(item["absolute_differences"], difference)
        for candidate in result["candidates"].values():
            for name in ("true_W", "true_Q2"):
                item = candidate[name]
                _finish(item["absolute_differences"])
                item["matching_fraction"] = item["matching_rows"] / result["checked_rows"] if result["checked_rows"] else None
                item["max_absolute_difference"] = item["absolute_differences"]["max"]
                item["rms_absolute_difference"] = item["absolute_differences"]["rms"]
        if not result["checked_rows"]:
            result.update(status="SKIP", reason="No usable stored truth rows for optional beam-formula check")
    except Exception as exc:
        result["errors"].append(_error(exc))
        result["status"] = "ERROR"
    return result
