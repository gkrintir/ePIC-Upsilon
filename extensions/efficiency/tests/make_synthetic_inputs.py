#!/usr/bin/env python3
"""Create explicitly synthetic ROOT fixtures for diagnostics integration tests.

The injected efficiency changes are intentionally artificial. These files are
not a production model and must never be used in the physics note.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
import uproot

MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    root = args.destination.resolve()
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(223314)
    samples = []
    common = None
    for key, state, mode, cap, n in [
            ("1s_newid", "1S", "new", 1., 6000),
            ("1s_standard", "1S", "standard", 1., 6000),
            ("2s_standard", "2S", "standard", .01, 4000),
            ("3s_standard", "3S", "standard", .01, 3200)]:
        if state == "1S" and common is not None:
            arrays = copy.deepcopy(common)
        else:
            pt = np.clip(rng.exponential(.65, n), 0, 3.4)
            q2 = 10 ** rng.uniform(-6, np.log10(cap), n)
            if state == "1S":
                tail = rng.random(n) < .1
                pt[tail] = rng.uniform(1.6, 3.2, tail.sum())
                q2[tail] = 10 ** rng.uniform(-1.6, 0, tail.sum())
            arrays = dict(
                true_mass=MASSES[state] + rng.normal(0, .002, n),
                true_Q2=q2,
                true_W=np.clip(30 + 18 * pt + rng.normal(0, 13, n), 15, 110),
                true_y=np.clip(1.4 - .35 * pt + rng.normal(0, .45, n), -1.5, 3.5),
                true_pT=pt,
            )
            if state == "1S":
                common = copy.deepcopy(arrays)
        pt, q2 = arrays["true_pT"], arrays["true_Q2"]
        high = q2 > .01
        probability = np.where(high, .32 + .17 * pt, .55 - .12 * pt)
        if mode == "new":
            probability *= .72
        probability = np.clip(probability, .03, .9)
        arrays["reco_mass"] = MASSES[state] + rng.normal(0, .06, n)
        arrays["reco_mass"][rng.random(n) > probability] = -999.
        arrays["reco_W"] = arrays["true_W"] + rng.normal(0, 3, n)
        arrays["reco_W"][rng.random(n) < .03 + .025 * pt] = 115.
        arrays["sigma_epz"] = rng.normal(18, .4, n)
        arrays["sigma_epz"][rng.random(n) < .04] = 22.
        # Exercise common schema failures without changing nominal normalization.
        for branch in ["true_mass", "true_Q2", "true_W", "true_y", "true_pT"]:
            arrays[branch][:3] = -999.
        arrays["true_Q2"][3] = -.1
        arrays["true_pT"][4] = np.nan
        arrays["true_y"][5] = np.nan
        # Synthetic coherent producer topology. Tagging never changes truth rows.
        combo = np.isfinite(arrays["reco_mass"]) & (arrays["reco_mass"] > 0)
        tagged = combo & (rng.random(n) < .65)
        central = rng.random(n) < .5
        arrays["reco_topology"] = np.where(combo, np.where(tagged, np.where(central, 1, 2), np.where(central, 3, 4)), 0).astype(np.int32)
        arrays["reco_tagged_proton_p"] = np.where(tagged, rng.uniform(150., 275., n), -999.)
        arrays["weight"] = np.full(n, .1234)
        file = root / (key + ".root")
        with uproot.recreate(file) as output:
            output.mktree("UpsilonTree", {name: "int32" if name == "reco_topology" else "float64" for name in arrays})
            output["UpsilonTree"].extend(arrays)
        samples.append(dict(key=key, label="SYNTHETIC " + state + " " + mode,
                            state=state, id_mode=mode, q2_max=cap, path=str(file)))
    config = dict(schema_version=1, energy="9x275", synthetic=True, samples=samples,
                  settings=dict(proton_selection="required", electron_energy=9., w_range=[20., 100.],
                                true_mass_window=.2, q2_split=.01,
                                coverage=.6826894921370859, min_denominator=50,
                                pt_edges=np.linspace(0, 3, 16).tolist(),
                                w_edges=np.linspace(20, 100, 11).tolist(),
                                y_edges=np.linspace(-1, 3, 17).tolist(),
                                q2_edges=[0, 1e-6, 1e-5, 1e-4, .001, .01, .03, .1, .3, 1],
                                tail_edges=[0, 1, 1.6, 2.2, 2.6, 3],
                                strata_w_edges=[20, 40, 60, 80, 100],
                                strata_y_edges=[-1, 0, 1, 2, 3]))
    config_path = root / "synthetic_config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    print(config_path)


if __name__ == "__main__":
    main()
