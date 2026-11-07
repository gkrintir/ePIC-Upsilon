"""Export the inclusive two-energy data and source PDFs for AN Figs. 10–12."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil


STATES = ("1S", "2S", "3S")
GROUPS = {
    "quarkonium": ("projected_{state}_upsilon_kinematics", 10),
    "decay": ("projected_{state}_decay_leptons", 11),
    "scattered": ("projected_{state}_scattered_electron_tracking", 11),
    "estimators": ("t_estimators_projective_transverse", 12),
}
BEAMS = {"9x130 GeV": "130", "9x275 GeV": "275"}


def endpoint(value, method):
    result = float(value)
    if math.isfinite(result):
        return result
    if method == "unavailable_empty_bin_correction" and math.isnan(result):
        return None
    raise ValueError(f"Unexpected non-finite endpoint for {method}")


def export_state(state, source, out):
    audit = json.loads((source / "kinematic_selection_audit.json").read_text())
    expected = {"state": state, "stat_mode": "projected", "projection_luminosity_invfb": 10,
                "source_luminosity_invfb": 10, "q2_cut": 1 if state == "1S" else 0.01,
                "projected_pages_mask": "selected", "count_interval": "median_garwood"}
    if any(audit.get(key) != value for key, value in expected.items()) or audit["proton_selection"]["mode"] != "inclusive":
        raise ValueError(f"{state}: wrong beam comparison or selection")
    groups = {key: (stem.format(state=state.lower()), fig) for key, (stem, fig) in GROUPS.items()
              if state == "1S" or key != "estimators"}
    page_lookup = {stem: key for key, (stem, _) in groups.items()}
    result = {}
    for key, (stem, fig) in groups.items():
        pdf = source / (stem + ".pdf")
        rel = Path("plots") / "beams" / f"fig{fig}-{state}-{key}.pdf"
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(pdf, out / rel)
        result[key] = {"pdf": str(rel), "sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
                       "observables": {}}
    with (source / "kinematic_bins.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["page"] not in page_lookup:
                continue
            group = result[page_lookup[row["page"]]]
            series = group["observables"].setdefault(row["observable"], {"130": [], "275": []})
            method = row["interval_method"]
            values = [float(row[k]) for k in ("bin_low", "bin_high", "display_value")]
            if not all(math.isfinite(x) for x in values) or values[0] >= values[1]:
                raise ValueError(f"Malformed bin for {row['observable']}")
            values += [endpoint(row[k], method) for k in ("display_interval_lower", "display_interval_upper")]
            series[BEAMS[row["sample"]]].append(values)
    for group in result.values():
        if not group["observables"]:
            raise ValueError(f"{state}: empty figure")
        for series in group["observables"].values():
            low, high = series["130"], series["275"]
            if not low or [x[:2] for x in low] != [x[:2] for x in high]:
                raise ValueError(f"{state}: incompatible bin edges")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--runs", type=Path, help="Workspace results/generated/runs directory")
    sources.add_argument("--input-dirs", type=Path, nargs=3, metavar=("DIR_1S", "DIR_2S", "DIR_3S"),
                         help="Three comparison directories containing kinematic_bins.csv and the PDFs")
    parser.add_argument("--batch-id", default="note-two-energy-L10")
    parser.add_argument("--out", type=Path, required=True, help="Site directory")
    args = parser.parse_args()
    directories = args.input_dirs or [args.runs / f"{args.batch_id}-{s.lower()}" /
                                      "artifacts/note/compare_9x130_9x275" for s in STATES]
    data = {s: export_state(s, directory, args.out) for s, directory in zip(STATES, directories)}
    target = args.out / "data" / "beams.js"
    target.parent.mkdir(parents=True, exist_ok=True)
    json_text = json.dumps(data, separators=(",", ":"), allow_nan=False)
    target.write_text("window.BEAM_FIGURES=" + json_text + ";\n")
    target.with_suffix(".json").write_text(json_text + "\n")


if __name__ == "__main__":
    main()
