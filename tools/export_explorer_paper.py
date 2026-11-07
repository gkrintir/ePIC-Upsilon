"""Prepare the local interactive figure data from a successful paper replay."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from PIL import Image

STATES = ("1S", "2S", "3S")
IMAGES = {"mass": "plots/mass.webp", "w": "plots/w.webp", "differential": "plots/differential.webp"}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def check(reference, rows, name):
    payload = json.dumps(sorted(rows), separators=(",", ":"), ensure_ascii=True).encode("ascii")
    if hashlib.sha256(payload).hexdigest() != reference[name]:
        raise ValueError(f"{name} differs from the paper reference")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path(__file__).resolve().parents[1] / "paper/reference/checksums.json")
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text())
    replay = args.replay.resolve()
    receipt = json.loads((replay / "paper_replay.json").read_text())
    if not receipt["spectra"]["all_pass"] or receipt["mass_values_and_errors_checked"] != 18:
        raise ValueError("Use a completed paper replay")
    fit = (replay / receipt["figure"]).resolve()
    if not fit.is_relative_to(replay) or not fit.is_file():
        raise ValueError("Missing reconstruction figure under the replay directory")
    fits = json.loads((fit.parent / "mass_fit_intervals.json").read_text())["states"]
    check(reference, [[s, q, k, format(fits[s][q][k], ".1f")] for s in STATES
                      for q in ("fitted_yield", "fitted_sigma_core_MeV") for k in ("value", "error_low", "error_high")], "mass_fit")
    mass = {s: {label: [fits[s][q][k] for k in ("value", "error_low", "error_high")]
                for label, q in (("yield", "fitted_yield"), ("width", "fitted_sigma_core_MeV"))} for s in STATES}

    cells = []
    edges = ("x_low", "x_high", "y_low", "y_high")
    counts = ("denominator", "selected_before_proton", "selected")
    for row in read_csv(replay / "reconstruction/efficiency_cells.csv"):
        cell = {"map": row["map"], **{k: float(row[k]) for k in edges}, **{k: int(row[k]) for k in counts}}
        cells.append(cell)
    check(reference, [[r["map"], *[format(r[k], ".12g") for k in edges], *[r[k] for k in counts]] for r in cells], "efficiency_cells")
    totals = receipt["xv_w_totals"]
    for key in counts:
        if sum(r[key] for r in cells if r["map"] == "xv_w") != totals[key]:
            raise ValueError(f"Map total disagrees with the replay for {key}")

    spectra = []
    fields = ("bin_low", "bin_high", "central", "lower", "upper")
    for row in read_csv(replay / "spectra/displayed_points.csv"):
        point = {k: row[k] for k in ("quantity", "state", "status")}
        point.update({k: float(row[k]) for k in fields})
        spectra.append(point)
    check(reference, [[p["quantity"], p["state"], *[format(p[k], ".7g") for k in fields], p["status"]] for p in spectra], "spectra")
    slope, = [r for r in spectra if r["quantity"] == "slope"]
    data = {"mass": mass, "efficiency": {"eligible": totals["denominator"], "pretag": totals["selected_before_proton"],
          "tagged": totals["selected"], "cells": cells}, "spectra": spectra,
          "slope": [slope["central"], slope["central"]-slope["lower"], slope["upper"]-slope["central"]], "images": IMAGES}
    out = args.out
    (out / "data").mkdir(parents=True, exist_ok=True)
    (out / "plots").mkdir(exist_ok=True)
    sources = {"mass": fit.with_suffix(".png"), "w": replay / "spectra/w/upsilon_es_w_cross_section.png",
               "differential": replay / "spectra/differential/upsilon_es_differential.png"}
    for key, source in sources.items():
        with Image.open(source) as image:
            image.thumbnail((1600, 1600), Image.Resampling.LANCZOS)
            image.convert("RGB").save(out / IMAGES[key], "WEBP", quality=92, method=6)
    content = json.dumps(data, separators=(",", ":"), allow_nan=False)
    (out / "data/paper.json").write_text(content + "\n")
    (out / "data/paper.js").write_text("window.PAPER=" + content + ";\n")
    print(f"Verified {len(cells)} map cells, {len(spectra)-1} spectrum bins, slope interval and all mass-fit results")


if __name__ == "__main__":
    main()
