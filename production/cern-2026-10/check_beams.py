"""Check the beam kinematics and inspect raw weights in the six reference tuples."""
import argparse
import csv
from pathlib import Path

import numpy as np
import uproot


MANIFEST = Path(__file__).resolve().parents[2] / "inputs/paper_inputs.csv"
ELECTRON_MASS = 0.000511
PROTON_MASS = 0.938272


def check(path, record):
    with path.open("rb") as stream, uproot.open(stream) as source:
        tree = source["UpsilonTree"]
        generated = int(record["generated_rows"])
        if tree.num_entries != generated or source["h_EventCount"].values().sum() != generated:
            raise ValueError(f"Generated-event count disagrees with the manifest: {path}")
        fields = ["ele_id_mode", "femc_mode", "weight", "true_mass", "true_W",
                  "true_scat_eMinus_pT", "true_scat_eMinus_eta", "true_scat_eMinus_E"]
        rows = tree.arrays(fields, library="np")
    for field in ("ele_id_mode", "femc_mode"):
        if not np.all(rows[field] == int(record[field])):
            raise ValueError(f"Unexpected {field}: {path}")

    state, beam = Path(record["relative_path"]).parts[:2]
    electron, proton = map(float, beam.split("x"))
    energy = rows["true_scat_eMinus_E"]
    pt = rows["true_scat_eMinus_pT"]
    eta = rows["true_scat_eMinus_eta"]
    w = rows["true_W"]
    # The producer writes sentinel values when no generated decay is identified.
    valid = (rows["true_mass"] > 0) & (w > 0) & (energy > 0) & (pt >= 0)
    valid &= np.isfinite(energy) & np.isfinite(pt) & np.isfinite(eta) & np.isfinite(w)
    if not np.any(valid):
        raise ValueError(f"No usable generated decays: {path}")
    energy, pt, eta = energy[valid], pt[valid], eta[valid]
    # Use eta only for the sign. E and pT retain precision near the beam axis.
    pz = np.sign(eta) * np.sqrt(np.maximum(0, energy**2 - ELECTRON_MASS**2 - pt**2))
    e_total = np.hypot(proton, PROTON_MASS) + electron - energy
    z_total = proton - electron - pz
    w2 = e_total**2 - z_total**2 - pt**2
    computed = np.sign(w2) * np.sqrt(np.abs(w2))
    error = np.max(np.abs(computed - w[valid]))
    if error > 1.e-8:
        raise ValueError(f"Stored true_W disagrees with the {beam} beam setting: {path}")
    weights = np.unique(rows["weight"])
    print(f"{state} {beam}: {valid.sum():,}/{generated:,} rows with generated kinematics, "
          f"max |ΔW| = {error:.3g} GeV, raw weight(s) = {weights.tolist()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, help="Verified input directory containing Upsilon1S, Upsilon2S and Upsilon3S")
    args = parser.parse_args()
    with MANIFEST.open(newline="") as handle:
        for entry in csv.DictReader(handle):
            check(args.inputs / entry["relative_path"], entry)
