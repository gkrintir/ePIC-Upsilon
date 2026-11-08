"""Reproduce the three proton-tagged Upsilon paper figures from merged tuples."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from explorer_results import results_json, write_atomic
from maps import render_and_check_maps
from fingerprints import mass_rows, verify

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
BACKEND = REPO / "analysis"
PROFILE = "jul26-purestates-estarlight-939b11a"

def run(script, arguments, env):
    subprocess.run([sys.executable, str(script), *map(str, arguments)], check=True, env=env)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw_base, out = args.inputs.resolve(), args.output.resolve()
    config = BACKEND / "config/reoptimization_9x275_wy.json"
    states = json.loads(config.read_text())["states"]
    runs = {state: "paper-" + state.lower() for state in states}
    env = dict(os.environ, YPSILON_PROTON_SELECTION="required", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg", MPLCONFIGDIR=str(out/"mplcache"))
    raw_files = {}
    with (REPO / "inputs/paper_inputs.csv").open(newline="") as manifest:
        for item in csv.DictReader(manifest):
            rel = Path(item["relative_path"])
            if rel.parts[1] != "9x275": continue
            src = raw_base / rel
            if src.stat().st_size != int(item["size_bytes"]):
                raise ValueError(f"Incorrect size: {src}")
            sha = hashlib.sha256()
            with src.open("rb") as stream:
                for chunk in iter(lambda: stream.read(4 * 1024**2), b""):
                    sha.update(chunk)
            if sha.hexdigest() != item["sha256"]:
                raise ValueError(f"Incorrect SHA-256: {src}")
            raw_files[rel.parts[0].removeprefix("Upsilon")] = src
    if set(raw_files) != set(states):
        raise ValueError("Expected one 9x275 input for each of 1S, 2S and 3S")
    out.mkdir(parents=True, exist_ok=False)
    inputs = []
    for state, parameters in states.items():
        folder = out/"results/generated/runs"/runs[state]
        artifacts = folder/"artifacts"
        folder.mkdir(parents=True)
        norm = artifacts/"inputs/9x275_normalized.root"
        rate, br, qmax = (parameters[k] for k in ("visible_xsec_fb","br_ee","q2_max"))
        run(BACKEND/"scripts/validation/prepare_latest_sample.py", [
            "--input",raw_files[state],"--cross-section-fb",rate,"--generator-branching-fraction",br,
            "--lumi-invfb",10,"--normalization-profile",PROFILE,
            "--normalization-source",f"Reproduced eSTARlight 939b11a with Jul26 Upsilon_Sim_Cards f62965a (9x275, Upsilon({state}), Q2 0-{qmax} GeV^2): visible dielectron {rate} fb; generator B_ee={br}",
            "--out-root",norm,"--out-json",artifacts/"normalization/9x275_input_audit.json"], env)
        perf = artifacts/"performance/9x275"
        run(BACKEND/"scripts/validation/plot_latest_sample_performance.py", [
            "--input",norm,"--state",state,"--electron-energy",9,"--energy-label","9x275 GeV",
            "--q2-cut",qmax,"--w-min",20,"--w-max",100,"--out-dir",perf,"--stat-mode","projected",
            "--count-interval","median_garwood","--interval-coverage",.6826894921370859,
            "--projection-lumi-invfb",2.5,"--projection-shape-mode","floating","--lumi-invfb",10], env)
        selection = json.loads((perf/"performance_audit.json").read_text())["proton_selection"]
        (folder/"proton_selection.json").write_text(json.dumps({"proton_selection":selection},indent=2)+"\n")
        inputs.extend(["--input",f"{state}={norm}"])
    result = render_and_check_maps(raw_files["1S"],out,runs,env)
    figure = out/result["figure"]
    fresh = json.loads((figure.parent/"mass_fit_intervals.json").read_text())
    verify("mass_fit", mass_rows(fresh))
    reopt = out/"results/generated/runs/paper-reopt"
    run(BACKEND/"workflow/run_reoptimization.py",["run","--config",config,
        "--framework-root",BACKEND,*inputs,"--states",*states,"--observables","pt2","abs_t_via_pt2","w","y",
        "--lumi",2.5,"--hash-inputs","--run-id",reopt.name,"--out-dir",reopt], env)
    run(HERE/"check_spectra.py",["--reopt-run",reopt,"--output",out/"spectra"],env)
    result["mass_values_and_errors_checked"] = 18
    result["spectra"] = json.loads((out/"spectra/comparison.json").read_text())
    write_atomic(out/"paper.json", results_json(out, result))
    write_atomic(out/"paper_replay.json", json.dumps(result,indent=2)+"\n")
    print(f"Paper reproduction complete: {out/'paper_replay.json'}")
    print(f"Figure explorer results: {out/'paper.json'}")

if __name__ == "__main__":
    main()
