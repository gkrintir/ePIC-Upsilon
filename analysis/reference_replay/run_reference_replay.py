#!/usr/bin/env python3
"""Execute frozen note-era scripts on original tuples; audit every archived CSV.

--backend root executes their original PyROOT readers. --backend uproot installs
the explicit limited I/O adapter; no physics code or archived script is changed.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import importlib.metadata
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parent
EDGE_STRING = "20,37.1429,45.7143,54.2857,62.8571,71.4286,80,88.5714,97.1429,105.714,114.286,122.857,140"
EXPECTED = {
    "low": "f714ae48f8903ffc83d34bddcf2cdcb56d409d26e84a430b2e875c35119658b8",
    "high": "d3eab02fc1aa72f3467f36eed31ea2970ce2dc85d545c96cfa81fa680dea7dda",
}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def worker(backend, script, arguments):
    os.environ["MPLBACKEND"] = "Agg"
    from atomic_plot_io import install
    install()
    if backend == "uproot":
        import uproot_io_adapter
        sys.modules["ROOT"] = uproot_io_adapter.module()
    elif backend != "root":
        raise ValueError(backend)
    sys.argv = [script, *arguments]
    runpy.run_path(script, run_name="__main__")


def audit_csv(ref, new, out):
    import numpy as np
    import pandas as pd
    a, b = pd.read_csv(ref), pd.read_csv(new)
    issues, rows = [], []
    if list(a.columns) != list(b.columns):
        issues.append({"column_difference": [list(a.columns), list(b.columns)]})
    if len(a) != len(b):
        issues.append({"rows": [len(a), len(b)]})
    maxrel, maxabs = 0., 0.
    for col in a.columns.intersection(b.columns):
        if col == "path":
            continue  # Input locations change; content hashes are checked separately.
        if len(a) != len(b):
            continue
        if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col]):
            av, bv = a[col].to_numpy(float), b[col].to_numpy(float)
            ok = np.isclose(av, bv, rtol=2e-9, atol=1e-9, equal_nan=True)
            finite = np.isfinite(av) & np.isfinite(bv)
            delta = bv - av
            rel = np.divide(delta, av, out=np.full_like(av, np.nan), where=av != 0)
            if finite.any():
                maxabs = max(maxabs, float(np.max(np.abs(delta[finite]))))
                finite_rel = finite & np.isfinite(rel)
                if finite_rel.any(): maxrel = max(maxrel, float(np.max(np.abs(rel[finite_rel]))))
            for i in range(len(a)):
                rows.append({"row": i, "column": col, "archived": av[i], "rerun": bv[i],
                             "difference": delta[i], "relative_difference": rel[i], "pass": bool(ok[i])})
            if not ok.all(): issues.append({"column": col, "failed_rows": np.flatnonzero(~ok).tolist()})
        else:
            ok = a[col].fillna("<NA>").astype(str).eq(b[col].fillna("<NA>").astype(str))
            if not ok.all(): issues.append({"column": col, "failed_rows": np.flatnonzero(~ok).tolist()})
    pd.DataFrame(rows).to_csv(out, index=False)
    return {"reference": str(ref), "rerun": str(new), "numeric_diff_csv": str(out),
            "pass": not issues, "issues": issues, "max_absolute_difference": maxabs,
            "max_relative_difference": maxrel, "byte_identical": sha(ref) == sha(new),
            "rtol": 2e-9, "atol": 1e-9, "ignored_columns": ["path"] if "path" in a else []}


def independent_checks(inputs):
    import numpy as np
    import uproot
    from uproot_io_adapter import RDataFrame
    rows = []
    for sample, file in inputs.items():
        with uproot.open(file) as f:
            a = f["UpsilonTree"].arrays(["true_mass", "true_Q2", "weight"], library="np")
            entries = f["UpsilonTree"].num_entries
        for state, mass in [("1S", 9.46030), ("2S", 10.02326), ("3S", 10.35520)]:
            mask = (a["true_mass"] >= mass - .2) & (a["true_mass"] <= mass + .2)
            expression = f"true_mass >= {mass-.2:.12g} && true_mass <= {mass+.2:.12g}"
            frame = RDataFrame("UpsilonTree", file).Filter(expression)
            direct = float(np.sum(a["weight"][mask], dtype=np.float64))
            qmask = mask & (a["true_Q2"] >= 0) & (a["true_Q2"] <= .001)
            count_pass = int(mask.sum()) == frame.Count().GetValue()
            sum_pass = np.isclose(direct, frame.Sum("weight").GetValue(), rtol=1e-14)
            rows.append(dict(sample=sample, state=state, entries=entries, state_entries=int(mask.sum()),
                             state_sumw_raw=direct, state_q2lt001_entries=int(qmask.sum()),
                             state_q2lt001_sumw_raw=float(a["weight"][qmask].sum()),
                             adapter_count_matches=bool(count_pass), adapter_sum_matches=bool(sum_pass)))
            if not (count_pass and sum_pass): raise AssertionError(rows[-1])
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--low-root", required=True, type=Path)
    p.add_argument("--high-root", required=True, type=Path)
    p.add_argument("--backend", choices=["auto", "root", "uproot"], default="auto")
    p.add_argument("--out-dir", type=Path, default=BASE / "results")
    p.add_argument("--skip-appendix", action="store_true")
    a = p.parse_args()
    backend = a.backend
    if backend == "auto": backend = "root" if importlib.util.find_spec("ROOT") else "uproot"
    inputs = {"low": a.low_root.resolve(), "high": a.high_root.resolve()}
    for k, f in inputs.items():
        if sha(f) != EXPECTED[k]:
            raise ValueError(f"{k} input is not the verified historical tuple: {f}")
    out = a.out_dir.resolve()
    for folder in [out, out / "logs", out / "differences"]: folder.mkdir(parents=True, exist_ok=True)
    manifest = {"backend": backend, "backend_description": "Original PyROOT readers" if backend == "root" else "Byte-preserved scripts with explicit uproot I/O adapter; NOT a PyROOT execution",
                "python_version": sys.version,
                "package_versions": {name: importlib.metadata.version(name) for name in ["numpy", "pandas", "matplotlib", "uproot", "Pillow"]},
                "mplhep_present": importlib.util.find_spec("mplhep") is not None,
                "inputs": {k: {"path": str(v), "sha256": sha(v)} for k, v in inputs.items()},
                "code_sha256": {str(f.relative_to(BASE)): sha(f) for f in [Path(__file__), BASE / "uproot_io_adapter.py", BASE / "atomic_plot_io.py", *sorted((BASE / "frozen_legacy").glob("*.py"))]},
                "plot_output_io": "Matplotlib PNGs serialize to memory, decode, write/fsync/verify temporary file and atomically replace target. Frozen plotting algorithms unchanged.",
                "runs": [], "checks": [], "independent_selection_checks": independent_checks(inputs),
                "conventions": {"lumi_invfb": 10., "tuple_lumi_invfb": 10., "flux_q2_max": .001,
                    "nominal_extraction_event_q2_cut": None, "truth_comparison_event_q2_cut": .001,
                    "normalization_scope": "truth_state", "weights": "visible dielectron; divide by BR exactly once",
                    "historical_uncertainties": "Preserved exactly; not the newer projected-count intervals",
                    "low_selected_note_recipe": "Recovered from frozen producer and original tuple: W10-25GeV/10bins,104.793fb,NO optional event Q2 cut,fluxQ2max=.001. Produces pixel-identical native overlay. No original command or matching low numerical CSV survives; separate raster verification establishes reproduction.",
                    "archived_low_selected_csv": "Quarantined: target304.793, no reconstructed events. Never used as analysis input.",
                    "flux_validation": "Archived CSV rerender; no new eSTARlight executable run"}}
    if manifest["mplhep_present"]:
        manifest["package_versions"]["mplhep"] = importlib.metadata.version("mplhep")
    def run(name, script, args):
        cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", backend, str(BASE / "frozen_legacy" / script), *map(str, args)]
        start = time.monotonic()
        with open(out / "logs" / f"{name}.log", "w") as log:
            completed = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env={**os.environ, "MPLBACKEND": "Agg"})
        manifest["runs"].append({"name": name, "command": cmd, "returncode": completed.returncode, "seconds": time.monotonic()-start})
        (out / "replay_manifest.json").write_text(json.dumps(manifest, indent=2))
        print(f"{name}: exit {completed.returncode}", flush=True)
        if completed.returncode: raise RuntimeError(f"Failed {name}: see {out/'logs'/f'{name}.log'}")
    common = ["--q2-max", ".001", "--lumi-invfb", "10", "--tuple-lumi-invfb", "10", "--weight-convention", "visible_decay", "--normalization-scope", "truth_state"]
    high = [inputs["high"], "--states", "1S,2S,3S", "--electron-energy", "18", "--proton-energy", "275", "--w-min", "20", "--w-max", "140", "--n-w-bins", "14", "--visible-xsec-1s-fb", "327.584", "--visible-xsec-2s-fb", "84.780", "--visible-xsec-3s-fb", "68.929", *common]
    low = [inputs["low"], "--states", "1S", "--electron-energy", "10", "--proton-energy", "100", "--visible-xsec-1s-fb", "104.793", *common]
    for name, version, params in [
        ("sigma_gp_family_18x275_fixed", "", [*high, "--min-eff", ".05"]),
        ("sigma_gp_family_18x275_fixed_unfolded", "_v2", [*high, "--include-eff-mc-unc"]),
        ("sigma_gp_1S_10x100_fixed", "", [*low, "--w-min", "20", "--w-max", "50", "--n-w-bins", "14"]),
        ("low_selected_candidate_no_event_q2_cut", "_v2", [*low, "--w-min", "10", "--w-max", "26.5", "--n-w-bins", "10", "--include-eff-mc-unc"]),
        ("low_selected_candidate_q2lt001", "_v2", [*low, "--w-min", "10", "--w-max", "26.5", "--n-w-bins", "10", "--include-eff-mc-unc", "--q2-cut", ".001"]),
        ("low_visual_geometry_hypothesis_W10_25", "_v2", [*low, "--w-min", "10", "--w-max", "25", "--n-w-bins", "10", "--include-eff-mc-unc", "--q2-cut", ".001"]),
        ("recovered_note_low_10x100", "_v2", [*low, "--w-min", "10", "--w-max", "25", "--n-w-bins", "10", "--include-eff-mc-unc"]),
    ]:
        run(name, f"extract_sigma_gamma_p_family_fixed_swan{version}.py", [*params, "--out-prefix", out/name])
    truth = "compare_10x100_18x275_truth_W20to26p5_Q2lt0p001_stateNorm"
    run(truth, "compare_sigma_gp_two_energies_swan_v4.py", ["--low-file", inputs["low"], "--high-file", inputs["high"], "--low-ee", "10", "--low-ep", "100", "--high-ee", "18", "--high-ep", "275", "--low-xsec-visible-fb", "104.793", "--high-xsec-visible-fb", "327.584", "--w-min", "20", "--w-max", "26.5", "--n-w-bins", "3", "--q2-max-flux", ".001", "--q2-cut", ".001", "--out-prefix", out/truth])
    if not a.skip_appendix:
        run("unfold_1S_18x275", "unfold_sigma_gp_W_swan.py", [inputs["high"], "--state", "1S", "--electron-energy", "18", "--proton-energy", "275", "--q2-max", ".001", "--q2-cut", ".001", "--w-edges", EDGE_STRING, "--visible-xsec-fb", "327.584", "--train-fraction", ".5", "--seed", "12345", "--n-iterations", "4", "--n-toys", "300", "--out-prefix", out/"unfold_1S_18x275"])
        run("Wres_1S_18x275_full", "study_W_resolution_swan.py", [inputs["high"], "--state", "1S", "--electron-energy", "18", "--q2-cut", ".001", "--w-edges", EDGE_STRING, "--out-prefix", out/"Wres_1S_18x275_full"])
    run("selected_note_candidate", "build_gamma_p_validation_figure_set.py", ["--indir", out, "--low-csv", "low_selected_candidate_q2lt001.csv", "--high-csv", "sigma_gp_family_18x275_fixed_unfolded.csv", "--truth-compare-csv", truth+".csv", "--flux-csv", BASE/"archived_reference/flux_ratio_scan_W20to26p5_Q2lt0p001.csv", "--min-eff", ".05", "--out-prefix", out/"selected_note_candidate"])
    run("selected_note_visual_hypothesis", "build_gamma_p_validation_figure_set.py", ["--indir", out, "--low-csv", "low_visual_geometry_hypothesis_W10_25.csv", "--high-csv", "sigma_gp_family_18x275_fixed_unfolded.csv", "--truth-compare-csv", truth+".csv", "--flux-csv", BASE/"archived_reference/flux_ratio_scan_W20to26p5_Q2lt0p001.csv", "--min-eff", ".05", "--out-prefix", out/"selected_note_visual_hypothesis"])
    run("recovered_note", "build_gamma_p_validation_figure_set.py", ["--indir", out, "--low-csv", "recovered_note_low_10x100.csv", "--high-csv", "sigma_gp_family_18x275_fixed_unfolded.csv", "--truth-compare-csv", truth+".csv", "--flux-csv", BASE/"archived_reference/flux_ratio_scan_W20to26p5_Q2lt0p001.csv", "--min-eff", ".05", "--out-prefix", out/"recovered_note"])
    for new in sorted(out.glob("*.csv")):
        ref = BASE/"archived_reference"/new.name
        if ref.exists(): manifest["checks"].append(audit_csv(ref, new, out/"differences"/(new.stem+"_differences.csv")))
    manifest["all_archival_csv_checks_pass"] = all(x["pass"] for x in manifest["checks"])
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from atomic_plot_io import install
    install()
    from PIL import Image
    manifest["native_image_checks"] = []
    for stem in ["overlay_1S_family", "18x275_family_efficiency", "closure", "truth_low_high_ratio", "flux_validation"]:
        reference_image = BASE/"archived_reference"/f"upsilon_selected_{stem}.png"
        rerun_image = out/f"recovered_note_{stem}.png"
        old_pixels = np.asarray(Image.open(reference_image).convert("RGBA"))
        new_pixels = np.asarray(Image.open(rerun_image).convert("RGBA"))
        same_shape = old_pixels.shape == new_pixels.shape
        identical = same_shape and np.array_equal(old_pixels,new_pixels)
        manifest["native_image_checks"].append({"reference": str(reference_image), "rerun": str(rerun_image),
            "reference_sha256": sha(reference_image), "rerun_sha256": sha(rerun_image), "pixel_identical": bool(identical),
            "reference_rgba_shape": list(old_pixels.shape), "rerun_rgba_shape": list(new_pixels.shape),
            "max_absolute_channel_difference": int(np.max(np.abs(old_pixels.astype(int)-new_pixels.astype(int)))) if same_shape else None,
            "note": "Decoded pixel comparison; PNG metadata/compression may differ. Flux panel is archival CSV rerender, other panels use recomputed tuple results."})
    manifest["all_native_reference_images_pixel_identical"] = all(c["pixel_identical"] for c in manifest["native_image_checks"])
    (out/"native_image_comparison.json").write_text(json.dumps(manifest["native_image_checks"],indent=2))
    frozen_plot = runpy.run_path(str(BASE/"frozen_legacy/build_gamma_p_validation_figure_set.py"))
    recovered_figure = frozen_plot["plot_final_overlay"](pd.read_csv(out/"recovered_note_low_10x100.csv"), pd.read_csv(out/"sigma_gp_family_18x275_fixed_unfolded.csv"),.05)
    recovered_figure.savefig(out/"recovered_note_overlay_1S_family.pdf")
    plt.close(recovered_figure)
    pd.DataFrame([{"file": Path(c["reference"]).name, "pass": c["pass"], "byte_identical": c["byte_identical"],
                   "max_absolute_difference": c["max_absolute_difference"], "max_relative_difference": c["max_relative_difference"]}
                  for c in manifest["checks"]]).to_csv(out/"reference_comparison_summary.csv", index=False)
    (out/"reference_comparison_summary.json").write_text(json.dumps(manifest["checks"], indent=2))
    candidates = []
    for stem in ["low_selected_candidate_no_event_q2_cut", "low_selected_candidate_q2lt001", "low_visual_geometry_hypothesis_W10_25", "recovered_note_low_10x100"]:
        table = pd.read_csv(out/(stem+".csv"));table.insert(0,"recipe",stem)
        table["passes_note_min_eff"] = table.eff_reco_over_true > .05
        table["archival_numeric_validation"] = "Native figure pixel identical; no matching surviving low CSV" if stem == "recovered_note_low_10x100" else "No matching surviving archival CSV"
        candidates.append(table)
    pd.concat(candidates, ignore_index=True).to_csv(out/"low_recipe_candidates.csv", index=False)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    for name, label, color in [("sigma_gp_1S_10x100_fixed", "10×100 1S", "black"), ("sigma_gp_family_18x275_fixed_unfolded", "18×275", None)]:
        old, new = pd.read_csv(BASE/"archived_reference"/(name+".csv")), pd.read_csv(out/(name+".csv"))
        for state, subset in new.groupby("state"):
            ix = subset.index; valid = subset.eff_reco_over_true > .05
            title = label if color else f"{label} {state}"
            ax[0, 0].plot(subset.W_center[valid], subset.sigma_gamma_p_pb[valid], "o", label=title, color=color)
            ax[0, 0].plot(old.loc[ix, "W_center"][valid], old.loc[ix, "sigma_gamma_p_pb"][valid], "x", color=color)
            ax[0, 1].plot(subset.W_center, subset.eff_reco_over_true, "o-", label=title, color=color)
            delta = (subset.sigma_gamma_p_pb-old.loc[ix,"sigma_gamma_p_pb"])/old.loc[ix,"sigma_gamma_p_pb"]
            ax[1, 0].plot(subset.W_center, delta, "o-", label=title, color=color)
            if "closure_unfold_over_truth" in subset:
                ax[1, 1].plot(subset.W_center[valid], subset.closure_unfold_over_truth[valid], "o-", label=title)
    ax[0,0].set_yscale("log");ax[0,0].set_ylabel("Cross section [pb]");ax[0,0].set_title("Rerun circles; archived CSV crosses")
    ax[0,1].set_ylabel("Reco/truth correction");ax[0,1].set_title("Recomputed normalization and efficiency")
    ax[1,0].set_ylabel("(rerun − archived) / archived");ax[1,0].set_title("Numerical reproduction residual")
    ax[1,1].set_ylabel("Corrected / truth");ax[1,1].set_title("Bin-by-bin algebraic closure (same sample)")
    for axis in ax.flat:
        axis.set_xlabel("W [GeV]", fontsize=11);axis.legend(fontsize=8);axis.grid(alpha=.2)
        axis.yaxis.label.set_size(10);axis.title.set_size(11);axis.tick_params(labelsize=9)
    fig.tight_layout();fig.savefig(out/"archived_vs_rerun_validation.png",dpi=180);fig.savefig(out/"archived_vs_rerun_validation.pdf");plt.close(fig)
    # A successful producer exit is insufficient: load the decoded pixels,
    # which also detects empty/truncated PNGs before publishing their hashes.
    png_checks = []
    for f in sorted(out.glob("*.png")):
        with Image.open(f) as decoded:
            decoded.load()
            png_checks.append({"file": f.name, "size_bytes": f.stat().st_size, "dimensions": list(decoded.size), "decoded": True})
    manifest["png_integrity_checks"] = png_checks
    manifest["all_generated_pngs_decodable"] = True
    manifest["outputs"] = [{"file": str(f.relative_to(out)), "sha256": sha(f)} for f in sorted(out.rglob("*")) if f.is_file() and f.name != "replay_manifest.json"]
    (out/"replay_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({"backend": backend, "checks": len(manifest["checks"]), "all_pass": manifest["all_archival_csv_checks_pass"]}))
    if not manifest["all_archival_csv_checks_pass"]: sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1] == "--worker": worker(sys.argv[2], sys.argv[3], sys.argv[4:])
    else: main()
