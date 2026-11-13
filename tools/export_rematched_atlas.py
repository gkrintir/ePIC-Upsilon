#!/usr/bin/env python3
"""Export the verified ECAL-rematched results and figures to the current Atlas."""

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

from PIL import Image


STATES = ("1S", "2S", "3S")
FRAMEWORK_COMMIT = "072c1fdc1c642a5cf21242b0363f24bcd89f84d5"
FIGURES = {
    "mass": "paper_figures/mass_efficiency/upsilon_mass_and_efficiency_w_xv",
    "w": "paper_figures/spectra/w/upsilon_es_w_cross_section",
    "differential": "paper_figures/spectra/differential/upsilon_es_differential",
}


class VerifiedArchive:
    """Read only source files whose bytes match the archive's recorded digest."""

    def __init__(self, root):
        self.root = root.resolve()
        manifest = (self.root / "SHA256SUMS.txt").read_bytes()
        self.manifest_sha256 = hashlib.sha256(manifest).hexdigest()
        self.expected = {}
        for line in manifest.decode("utf-8").splitlines():
            checksum, name = line.split("  ", 1)
            if name in self.expected:
                raise ValueError("Duplicate archive checksum entry for " + name)
            self.expected[name] = checksum
        self.checked = {}

    def read(self, relative):
        path = (self.root / relative).resolve()
        if not path.is_relative_to(self.root) or relative not in self.expected:
            raise ValueError("Unrecorded archive source " + relative)
        payload = path.read_bytes()
        checksum = hashlib.sha256(payload).hexdigest()
        if checksum != self.expected[relative]:
            raise ValueError("Archive SHA-256 mismatch for " + relative)
        self.checked[relative] = {"path": relative, "sha256": checksum, "bytes": len(payload)}
        return payload

    def json(self, relative):
        return json.loads(self.read(relative))

    def csv(self, relative):
        return list(csv.DictReader(io.StringIO(self.read(relative).decode("utf-8"))))


def finite(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Expected a finite numerical result")
    return number


def count(value):
    number = finite(value)
    if number < 0 or not number.is_integer():
        raise ValueError("Expected a nonnegative integer count")
    return int(number)


def boolean(value):
    if str(value).lower() not in ("true", "false"):
        raise ValueError("Expected an explicit boolean")
    return str(value).lower() == "true"


def mass_results(archive):
    fit = archive.json("paper_figures/mass_efficiency/mass_fit_intervals.json")
    if fit["target_luminosity_invfb"] != 2.5 or set(fit["states"]) != set(STATES):
        raise ValueError("Expected the three-state 9x275, 2.5 inverse-femtobarn fit")
    mass = {}
    for state in STATES:
        mass[state] = {}
        for label, key in (("yield", "fitted_yield"), ("width", "fitted_sigma_core_MeV")):
            values = fit["states"][state][key]
            mass[state][label] = [finite(values[field]) for field in ("value", "error_low", "error_high")]
            if values["status"] != "TWO_SIDED" or min(mass[state][label]) < 0:
                raise ValueError("Unsupported mass-fit interval for " + state)
    keys = (
        "method", "nominal_coverage", "coverage_validated", "publication_ready", "interpretation",
        "assumptions", "resolution_correlation", "upstream_shape_mismatch_flags",
    )
    return mass, {key: fit[key] for key in keys}


def efficiency_results(archive):
    crosscheck = archive.json("paper_figures/mass_efficiency/independent_map_crosscheck.json")
    totals = crosscheck["all_eligible_and_selected_totals"]
    fields = ("denominator", "selected_before_proton", "selected")
    cells = []
    for row in archive.csv("paper_figures/mass_efficiency/efficiency_map_cells.csv"):
        cell = {"map": row["map"]}
        cell.update({key: finite(row[key]) for key in ("x_low", "x_high", "y_low", "y_high")})
        cell.update({key: count(row[key]) for key in fields})
        if not 0 <= cell["selected"] <= cell["selected_before_proton"] <= cell["denominator"]:
            raise ValueError("Inconsistent efficiency-map counts")
        cells.append(cell)
    for name, expected in crosscheck["maps"].items():
        selected = [cell for cell in cells if cell["map"] == name]
        sums = {key: sum(cell[key] for cell in selected) for key in fields}
        if len(selected) != expected["cells_verified"] or sums != expected["in_displayed_range"]:
            raise ValueError("Efficiency cells disagree with independent check for " + name)
        if any(sums[key] + expected["outside_displayed_range"][key] != totals[key] for key in fields):
            raise ValueError("Efficiency overflow accounting disagrees for " + name)
    return {
        "eligible": count(totals["denominator"]), "pretag": count(totals["selected_before_proton"]),
        "tagged": count(totals["selected"]), "cells": cells,
    }, crosscheck


def point(quantity, row, edges, values):
    result = {"quantity": quantity, "state": row["state"], "status": row["overall_status"]}
    result.update({key: finite(row[field]) for key, field in zip(("bin_low", "bin_high"), edges)})
    result.update({key: finite(row[field]) for key, field in zip(("central", "lower", "upper"), values)})
    if (result["state"] not in STATES or result["bin_low"] >= result["bin_high"]
            or min(result["central"], result["lower"]) < 0 or result["upper"] < result["lower"]):
        raise ValueError("Invalid displayed spectrum bin")
    return result


def spectral_results(archive):
    prefix = "paper_figures/spectra/"
    w_audit = archive.json(prefix + "w/upsilon_es_w_cross_section_audit.json")
    differential_audit = archive.json(prefix + "differential/upsilon_es_differential_audit.json")
    if w_audit["energy"] != "9x275" or w_audit["target_luminosity_invfb"] != 2.5:
        raise ValueError("Unexpected W-figure beam or luminosity")
    if differential_audit["energy"] != "9x275" or differential_audit["target_lumi_invfb"] != 2.5:
        raise ValueError("Unexpected differential-figure beam or luminosity")

    spectra = []
    for row in archive.csv(prefix + "w/upsilon_es_w_cross_section_points.csv"):
        if boolean(row["displayed"]):
            spectra.append(point("W", row, ("W_low", "W_high"),
                                 ("sigma_gamma_p_pb", "interval_low_pb", "interval_high_pb")))
    for row in archive.csv(prefix + "differential/upsilon_es_differential_points.csv"):
        if boolean(row["displayed"]):
            spectra.append(point(row["observable"], row, ("x_low", "x_high"), ("mean", "lower", "upper")))

    availability = []
    for state in STATES:
        metadata = w_audit["states"][state]["sample_metadata"]
        reason = "Adaptive-binning W diagnostic with failed injected-spectrum coverage checks."
        if state in ("2S", "3S"):
            reason += " The separate fixed-bin W unfolding is withheld because a truth-populated bin has zero reconstruction efficiency."
        availability.append({
            "quantity": "W", "state": state, "available": bool(metadata["available"]),
            "status": metadata["overall_status"], "reason": reason, "failed_checks": metadata["failed_checks"],
        })
    for case in differential_audit["case_availability"]:
        reason = case["reason"]
        if case["state"] == "2S" and case["observable"] == "pt2":
            reason = "Diagnostic spectrum and slope. The injected-slope coverage check failed."
        availability.append({
            "quantity": case["observable"], "state": case["state"], "available": bool(case["available"]),
            "status": case["overall_status"], "reason": reason, "failed_checks": case["failed_checks"],
        })

    slopes = {}
    for row in archive.csv(prefix + "differential/upsilon_es_slopes.csv"):
        state = row["state"]
        case = next(item for item in availability if item["state"] == state and item["quantity"] == "pt2")
        slope = {key: case[key] for key in ("available", "status", "reason", "failed_checks")}
        for key, column in (
            ("value", "b_slope_GeVminus2"), ("error_low", "b_slope_err_low"),
            ("error_high", "b_slope_err_high"), ("lower", "b_interval_low"), ("upper", "b_interval_high"),
        ):
            slope[key] = finite(row[column]) if case["available"] else None
        if row["overall_status"] != case["status"]:
            raise ValueError("Slope and availability statuses disagree for " + state)
        slopes[state] = slope
        if case["available"]:
            spectra.append({
                "quantity": "slope", "state": state, "bin_low": 0.0, "bin_high": 1.0,
                "central": slope["value"], "lower": slope["lower"], "upper": slope["upper"],
                "status": slope["status"],
            })
    if set(slopes) != set(STATES):
        raise ValueError("Missing state-specific slope status")
    for case in availability:
        rows = [row for row in spectra if row["state"] == case["state"] and row["quantity"] == case["quantity"]]
        if bool(rows) != case["available"] or any(row["status"] != case["status"] for row in rows):
            raise ValueError("Displayed bins disagree with their availability or validation status")
    return spectra, slopes, availability


def input_provenance(archive):
    inputs = archive.json("sources/verified_inputs.json")
    audit = archive.json("audit/input_audit.json")
    by_hash = {row["original_file_sha256"]: row for row in audit}
    verified = []
    for source in inputs:
        row = by_hash[source["sha256"]]
        if source["bytes"] != row["size_bytes"] or row["tuple_rows"] != row["h_EventCount"]:
            raise ValueError("Input identity or generated-event accounting disagrees")
        verified.append({
            **source, "state": row["state"], "beam": row["beam"], "generated_rows": row["tuple_rows"],
            "eligible_truth_rows": row["eligible_truth_rows"], "visible_cross_section_fb": row["rate_fb"],
            "electron_id_mode": row["ele_id_mode"], "femc_mode": row["femc_mode"],
        })
    if {(row["state"], row["beam"]) for row in verified} != {(state, beam) for state in STATES for beam in ("9x130", "9x275")}:
        raise ValueError("Expected all six corrected state/beam inputs")
    return verified


def validation_provenance(archive):
    tasks = archive.csv("framework/native_task_status.csv")
    failures = [{key: row[key] for key in ("state", "energy", "task", "status")}
                for row in tasks if row["status"] in ("INVALID_OUTPUT", "BLOCKED_DEPENDENCY")]
    validations = []
    for row in archive.csv("framework/native_reoptimization_status.csv"):
        validations.append({
            "energy": row["energy"], "luminosity_invfb": finite(row["lumi_invfb"]),
            "state": row["state"], "quantity": row["observable"], "status": row["overall_status"],
            "failed_checks": json.loads(row["failed_checks"]), "publication_ready": boolean(row["publication_ready"]),
        })
    comparisons = []
    for row in archive.csv("audit/selected_old_vs_rematched.csv"):
        comparisons.append({key: row[key] if key in ("state", "beam", "category") else finite(row[key]) for key in row})
    return {"fixed_bin_withheld_tasks": failures, "reoptimization": validations, "selected_old_vs_rematched": comparisons}


def export(results, destination):
    archive = VerifiedArchive(results)
    if FRAMEWORK_COMMIT not in archive.read("README.md").decode("utf-8"):
        raise ValueError("Archive analysis-commit provenance differs from this rematched dataset")
    mass, mass_validation = mass_results(archive)
    efficiency, map_check = efficiency_results(archive)
    spectra, slopes, availability = spectral_results(archive)
    inputs = input_provenance(archive)
    validation = validation_provenance(archive)
    runs = archive.json("paper_figures/spectra/input_runs_by_state.json")
    data = {
        "schema_version": 2, "dataset": "rematch3011", "energy": "9x275", "luminosity_invfb": 2.5,
        "mass": mass, "efficiency": efficiency, "spectra": spectra,
        "slope": [slopes["1S"][key] for key in ("value", "error_low", "error_high")],
        "slopes": slopes, "availability": availability, "mass_validation": mass_validation,
        "images": {name: f"plots/{name}.webp" for name in FIGURES},
        "pdfs": {name: f"plots/{name}.pdf" for name in FIGURES},
        "image_dimensions": {}, "provenance": "data/paper-provenance.json",
    }
    provenance = {
        "schema_version": 1, "dataset": data["dataset"], "energy": data["energy"],
        "luminosity_invfb": data["luminosity_invfb"], "framework_commit": FRAMEWORK_COMMIT,
        "matching": "Downstream EICrecon issue 3011 ECAL rematching, full forward ECAL",
        "proton_selection": "reco_topology in {1,2} and finite positive reco_tagged_proton_p",
        "exporter": {"path": "tools/export_rematched_atlas.py",
                     "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "archive_manifest_sha256": archive.manifest_sha256,
        "verified_inputs": inputs, "input_runs_by_state": {state: Path(path).name for state, path in runs.items()},
        "efficiency_map_check": map_check, "mass_validation": mass_validation,
        "availability": availability, "slopes": slopes, **validation,
    }
    destination.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".rematched-atlas-", dir=destination) as folder:
        stage = Path(folder)
        (stage / "data").mkdir()
        (stage / "plots").mkdir()
        for name, stem in FIGURES.items():
            with Image.open(io.BytesIO(archive.read(stem + ".png"))) as source:
                source.thumbnail((1600, 1600), Image.Resampling.LANCZOS)
                source.convert("RGB").save(stage / data["images"][name], "WEBP", quality=92, method=6)
                data["image_dimensions"][name] = list(source.size)
            (stage / data["pdfs"][name]).write_bytes(archive.read(stem + ".pdf"))
        for group in ("images", "pdfs"):
            for name, relative in data[group].items():
                version = hashlib.sha256((stage / relative).read_bytes()).hexdigest()[:12]
                data[group][name] = relative + "?v=" + version
        content = json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n"
        (stage / "data/paper.json").write_text(content, encoding="utf-8")
        (stage / "data/paper.js").write_text("window.PAPER=" + content.rstrip() + ";\n", encoding="utf-8")
        provenance["source_files"] = sorted(archive.checked.values(), key=lambda row: row["path"])
        provenance["exported_files"] = [
            {"path": str(path.relative_to(stage)), "bytes": path.stat().st_size,
             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in sorted(stage.rglob("*")) if path.is_file()
        ]
        (stage / "data/paper-provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        for path in sorted(stage.rglob("*")):
            if path.is_file():
                target = destination / path.relative_to(stage)
                target.parent.mkdir(parents=True, exist_ok=True)
                path.replace(target)
    print("Exported corrected Atlas figures and numerical results to", destination)
    print("Spectrum rows by quantity", dict(Counter(row["quantity"] for row in spectra)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Verified rematched results archive directory")
    parser.add_argument("--out", type=Path, required=True, help="Atlas directory containing data/ and plots/")
    args = parser.parse_args()
    export(args.results.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
