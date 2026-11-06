"""Build the published mass/efficiency figure and compare both maps cell by cell."""
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
from fingerprints import map_rows, verify
EDGES = ["x_low", "x_high", "y_low", "y_high"]
COUNTS = ["denominator", "selected_before_proton", "selected"]


def render_and_check_maps(raw_1s, output_root, full_runs, env):
    """Run the efficiency extension followed by notebook 07 on the same tagged tuples."""
    work = Path(output_root).resolve()
    runs = work / "results/generated/runs"
    profile = json.loads((REPO / "configs/plots_9x275.paper.json").read_text())
    sample_key = profile["map_sample"]
    extension = REPO / "extensions/efficiency"
    config = json.loads((extension / "config_9x275.json").read_text())
    sample, = [s for s in config["samples"] if s["key"] == sample_key]
    sample["path"] = str(Path(raw_1s).resolve())
    sample.pop("reference_efficiency_csv", None)
    config["samples"] = [sample]
    config_path = work / "paper_efficiency_config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    map_run = "paper-1s-tagged-maps"
    subprocess.run([
        sys.executable, str(extension / "run_final_efficiency_2d.py"), "run",
        "--framework-root", str(work), "--config", str(config_path),
        "--sample", sample_key, "--run-id", map_run, "--proton-selection", "required"
    ], check=True, env=env)

    plot_root = work / "reconstruction"
    plot_root.mkdir()
    profile.update(results_root=str(runs), full_runs=full_runs,
                   efficiency_2d_run=map_run, root_input_override=sample["path"],
                   output_root=str(plot_root))
    profile_path = work / "reconstruction_profile.json"
    profile_path.write_text(json.dumps(profile, indent=2) + "\n")
    # IPython executes the same cells as SWAN and propagates a failed cell.
    # The captured paths are those returned by this execution's plotting code.
    notebook = "07_mass_and_efficiency_w_xv.ipynb"
    receipt = plot_root / "outputs.json"
    execute = (
        "import json, sys\n"
        "from pathlib import Path\n"
        "from IPython.core.interactiveshell import InteractiveShell\n"
        "shell = InteractiveShell.instance()\n"
        "shell.safe_execfile_ipy(sys.argv[1], raise_exceptions=True)\n"
        "paths = {k: str(shell.user_ns[k]) for k in ('OUTPUT_DIR', 'wxv_cache')}\n"
        "Path(sys.argv[2]).write_text(json.dumps(paths, indent=2))\n"
    )
    subprocess.run([sys.executable, "-c", execute, notebook, str(receipt)],
                   cwd=REPO / "plotting", check=True,
                   env=dict(env, UPSILON_PLOT_PROFILE=str(profile_path)))
    paths = json.loads(receipt.read_text())
    figure = Path(paths["OUTPUT_DIR"]) / "upsilon_mass_and_efficiency_w_xv.pdf"
    map_dir = runs / map_run / "artifacts/final_efficiency_2d" / sample_key
    pt_y = pd.read_csv(map_dir / "tables/final_efficiency_2d.csv")
    xv_w = pd.read_csv(Path(paths["wxv_cache"]) / "xv_w_efficiency.csv")
    actual = pd.concat([pt_y.loc[pt_y["map"] == "pt_y"], xv_w], ignore_index=True)
    order = ["map", *EDGES]
    actual = actual.sort_values(order).reset_index(drop=True)
    verify("efficiency_cells", map_rows(actual))
    if not figure.is_file():
        raise FileNotFoundError(figure)
    actual[order + COUNTS].to_csv(plot_root / "efficiency_cells.csv", index=False)
    return {
        "figure": str(figure.relative_to(work)),
        "cells_checked": actual.groupby("map").size().to_dict(),
        "xv_w_totals": xv_w[COUNTS].sum().to_dict(),
    }
