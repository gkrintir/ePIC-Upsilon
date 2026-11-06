"""Keep saved category keys separate from verified tuple reconstruction modes."""
import csv
from pathlib import Path


MANIFEST = Path(__file__).resolve().parents[2] / "inputs/paper_inputs.csv"


def electron_id_metadata(saved_category, raw_sha256, state, energy):
    metadata = {"saved_category": saved_category, "category": saved_category,
                "category_source": "saved sample category; tuple electron-ID mode unverified"}
    if not raw_sha256 or not MANIFEST.is_file():
        return metadata
    with MANIFEST.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["sha256"] != raw_sha256:
                continue
            expected = Path(row["relative_path"]).parts[:2]
            if expected != ("Upsilon" + state, energy):
                raise ValueError(f"{state}/{energy}: raw tuple checksum belongs to {expected}")
            mode = int(row["ele_id_mode"])
            metadata.update(category="reconstructed track–cluster matching" if mode == 1 else f"producer mode {mode}",
                            category_source="raw SHA-256 matched inputs/paper_inputs.csv",
                            ele_id_mode=mode, femc_mode=int(row["femc_mode"]), raw_sha256=raw_sha256)
            break
    return metadata
