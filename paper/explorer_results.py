"""Build the portable numerical results for the paper figure explorer."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from fingerprints import REFERENCE, cell_rows, mass_rows, spectrum_rows, verify


STATES = ("1S", "2S", "3S")
IMAGES = {"mass": "plots/mass.webp", "w": "plots/w.webp", "differential": "plots/differential.webp"}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def results_json(replay, receipt, reference=REFERENCE):
    """Serialize the checked fit, map cells and spectrum bins for an explorer."""
    if not receipt["spectra"]["all_pass"] or receipt["mass_values_and_errors_checked"] != 18:
        raise ValueError("Use a completed paper replay")
    figure = (replay / receipt["figure"]).resolve()
    if not figure.is_relative_to(replay) or not figure.is_file():
        raise ValueError("Missing reconstruction figure under the replay directory")

    fits = json.loads((figure.parent / "mass_fit_intervals.json").read_text())["states"]
    quantities = {"yield": "fitted_yield", "width": "fitted_sigma_core_MeV"}
    fitted = {state: {key: fits[state][key] for key in quantities.values()} for state in STATES}
    verify("mass_fit", mass_rows({"states": fitted}), reference)
    mass = {
        state: {label: [fitted[state][key][k] for k in ("value", "error_low", "error_high")]
                for label, key in quantities.items()}
        for state in STATES
    }

    cells = []
    edges = ("x_low", "x_high", "y_low", "y_high")
    counts = ("denominator", "selected_before_proton", "selected")
    for row in read_csv(replay / "reconstruction/efficiency_cells.csv"):
        cells.append({"map": row["map"], **{k: float(row[k]) for k in edges},
                      **{k: int(row[k]) for k in counts}})
    verify("efficiency_cells", cell_rows(cells), reference)
    totals = receipt["xv_w_totals"]
    for key in counts:
        if sum(cell[key] for cell in cells if cell["map"] == "xv_w") != totals[key]:
            raise ValueError(f"Map total disagrees with the replay for {key}")

    spectra = []
    fields = ("bin_low", "bin_high", "central", "lower", "upper")
    for row in read_csv(replay / "spectra/displayed_points.csv"):
        point = {k: row[k] for k in ("quantity", "state", "status")}
        point.update({k: float(row[k]) for k in fields})
        spectra.append(point)
    verify("spectra", spectrum_rows(spectra), reference)
    slope, = [row for row in spectra if row["quantity"] == "slope"]
    data = {
        "mass": mass,
        "efficiency": {"eligible": totals["denominator"], "pretag": totals["selected_before_proton"],
                       "tagged": totals["selected"], "cells": cells},
        "spectra": spectra,
        "slope": [slope["central"], slope["central"] - slope["lower"], slope["upper"] - slope["central"]],
        "images": IMAGES,
    }
    return json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n"


def write_atomic(path, text):
    with TemporaryDirectory(dir=path.parent) as staging:
        staged = Path(staging) / path.name
        staged.write_text(text, encoding="utf-8")
        staged.replace(path)
