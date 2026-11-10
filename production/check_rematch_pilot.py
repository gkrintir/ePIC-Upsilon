"""Compare stored and rematched tuple outputs for either documented 100-event pilot."""
import argparse
import json
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import uproot


# Zero-based selected rows from an independent event replay of the two EDM4eic inputs.
# Each pair lists the rows before and after requiring the recoil proton.
PILOTS = {
    "f2701e1d496d931b117afd5760efe15c909061668890e2d67d6fd28d29856908": {
        "reco": ([4, 9, 18, 27, 38, 39, 45, 63, 71, 76, 85, 86, 93],
                 [4, 9, 39, 45, 63, 71, 76, 86]),
        "reco-rematch3011": ([4, 5, 8, 9, 18, 24, 27, 36, 38, 39, 45, 48, 59, 60, 63, 66, 67, 71, 76, 78, 85, 86, 93],
                             [4, 9, 24, 36, 39, 45, 60, 63, 67, 71, 76, 78, 86]),
    },
    "b5f8450deba55cc4555532f11f87be1d81713ef3cefc51774f84aa7b756a01ff": {
        "reco": ([1, 6, 11, 15, 20, 21, 23, 24, 26, 34, 35, 45, 46, 49, 55, 60, 63, 65, 67, 68, 72, 76, 97],
                 [1, 6, 49, 60, 65]),
        "reco-rematch3011": ([1, 6, 11, 15, 19, 20, 21, 22, 23, 24, 26, 31, 34, 35, 36, 44, 45, 46, 49, 55, 60, 61, 63, 65, 67, 68, 72, 76, 81, 84, 85, 91, 97],
                             [1, 6, 19, 36, 49, 60, 65, 85, 91]),
    },
}


def selected_rows(tree):
    names = ["reco_mass", "reco_W", "sigma_epz", "reco_topology", "reco_tagged_proton_p"]
    data = tree.arrays(names, library="np")
    mass, w = data["reco_mass"], data["reco_W"]
    selected = ((mass >= 8) & (mass <= 12) & (w >= 20) & (w <= 100)
                & (np.abs(data["sigma_epz"] - 18) <= 1.8))
    tagged = np.isin(data["reco_topology"], [1, 2]) & (data["reco_tagged_proton_p"] > 0)
    return np.flatnonzero(selected).tolist(), np.flatnonzero(selected & tagged).tolist()


def verify(stored_job, rematched_job):
    with ExitStack() as files:
        trees = []
        first_digest = None
        for mode, path in [("reco", stored_job), ("reco-rematch3011", rematched_job)]:
            job = json.loads(path.read_text())
            digest = job["inputs"][0].get("sha256")
            if digest not in PILOTS or (first_digest is not None and first_digest != digest):
                raise ValueError("Supply both jobs for the same documented pilot, planned with --hash-inputs.")
            if job["electron_id"] != mode or job["femc"] != "full":
                raise ValueError(f"Expected {mode} and full fEMC in {path}")
            first_digest = digest
            stream = files.enter_context(Path(job["output"]).open("rb"))
            root = files.enter_context(uproot.open(stream))
            tree = root["UpsilonTree"]
            if tree.num_entries != 100 or root["h_EventCount"].values().sum() != 100:
                raise ValueError(f"Expected 100 generated events for {mode}")
            rows = selected_rows(tree)
            if rows != PILOTS[digest][mode]:
                raise ValueError(f"{mode}: selected event indices differ from the independent replay: {rows}")
            print(f"{job['state']} {job['beam']}, {mode}: {len(rows[0])} before proton tag, {len(rows[1])} tagged")
            trees.append(tree)
        stored, rematched = trees
        if set(stored.keys()) != set(rematched.keys()):
            raise ValueError("Tuple schema changed")
        for key in [k for k in stored.keys() if k.startswith("true_")] + ["weight"]:
            original = stored[key].array(library="np")
            corrected = rematched[key].array(library="np")
            if original.dtype != corrected.dtype or original.tobytes() != corrected.tobytes():
                raise ValueError(f"Generated rows or input weights differ in {key}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stored_job", type=Path, help="Completed reco job.json")
    parser.add_argument("rematched_job", type=Path, help="Completed reco-rematch3011 job.json")
    args = parser.parse_args()
    verify(args.stored_job, args.rematched_job)
