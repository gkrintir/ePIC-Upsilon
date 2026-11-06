"""Render the two paper spectrum figures and check their numerical fingerprints."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

from fingerprints import spectrum_rows, verify


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reopt-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    run, out = args.reopt_run.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    profile = json.loads((repo/"configs/plots_9x275.paper.json").read_text())
    profile.update(results_root=str(run.parent), reoptimization_run=run.name, output_root=str(out))
    profile_path = out/"paper_profile.json"
    profile_path.write_text(json.dumps(profile, indent=2)+"\n")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MPLBACKEND="Agg", MPLCONFIGDIR=str(out/"mpl"))
    for section in ("w","differential"):
        cmd = [sys.executable, str(repo/"plotting/render_all.py"), "--section", section,
               "--profile", str(profile_path), "--allow-diagnostics", "--out", str(out/section)]
        with (out/(section+".log")).open("w") as log:
            subprocess.run(cmd, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)

    quantities = [
        ("W", out/"w/upsilon_es_w_cross_section_points.csv", ("W_low","W_high"),
         ("sigma_gamma_p_pb","interval_low_pb","interval_high_pb")),
        ("y", out/"differential/upsilon_es_differential_points.csv", ("x_low","x_high"),
         ("mean","lower","upper")),
        ("pt2", out/"differential/upsilon_es_differential_points.csv", ("x_low","x_high"),
         ("mean","lower","upper")),
    ]
    checks = []
    for name, path, edges, values in quantities:
        for row in read_csv(path):
            if row["displayed"].lower() == "true" and (name == "W" or row["observable"] == name):
                checks.append(dict(quantity=name,state=row["state"],bin_low=row[edges[0]],bin_high=row[edges[1]],
                                   central=row[values[0]],lower=row[values[1]],upper=row[values[2]],
                                   status=row["overall_status"]))
    slopes = {row["state"]: row for row in read_csv(out/"differential/upsilon_es_slopes.csv")}
    for state in ("2S","3S"):
        if slopes[state]["overall_status"] != "NO_FEASIBLE_SLOPE_BINNING" or slopes[state]["b_slope_GeVminus2"]:
            raise ValueError(f"{state}: unexpected slope availability")
    slope = slopes["1S"]
    checks.append(dict(quantity="slope",state="1S",bin_low=0,bin_high=1,central=slope["b_slope_GeVminus2"],
                       lower=slope["b_interval_low"],upper=slope["b_interval_high"],status=slope["overall_status"]))
    with (out/"displayed_points.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(checks[0]))
        writer.writeheader()
        writer.writerows(checks)
    sha = verify("spectra", spectrum_rows(checks))
    summary = {"all_pass": True, "numbers": 3*len(checks), "significant_digits": 7, "sha256": sha}
    (out/"comparison.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(f"Verified {summary['numbers']} paper values and the corresponding validation statuses")


if __name__ == "__main__":
    main()
