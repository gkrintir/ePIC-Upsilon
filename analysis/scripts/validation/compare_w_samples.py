#!/usr/bin/env python3
"""Recompute historical and July W comparisons from actual UpsilonTree rows.

These are truth-spectrum normalization diagnostics with conditional weighted-MC
errors. Reconstruction support is shown separately; R/(R/T)=T is not an
independent closure test. No projected Poisson errors or fitted scale factors.
The photon flux and mass/exclusivity selections follow the framework's
extract_gamma_p_cross_sections_bin_by_bin.py. Historical normalization is fixed
from the truth-state sum BEFORE Q2/W cuts; July normalization is target*L/Ngen.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from branching_fractions import positive_branching_fraction

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813
BRANCHES = ("weight", "true_mass", "true_W", "true_Q2", "reco_mass", "reco_W", "sigma_epz")
ERROR_CONVENTION = "sqrt(sum normalized_weight^2), conditional on the fixed normalization; no normalization, flux, luminosity, branching-fraction, response, or model uncertainty"


def photon_flux_dW(w, electron_energy, proton_energy, q2_max):
    s_minus_mp2 = 4.0 * electron_energy * proton_energy
    y = (w * w - M_P * M_P) / s_minus_mp2
    if not 0.0 < y < 1.0:
        return 0.0
    q2_min = M_E * M_E * y * y / (1.0 - y)
    if not 0.0 < q2_min < q2_max:
        return 0.0
    dndy = ALPHA_EM / (2.0 * math.pi) * (
        (1.0 + (1.0 - y)**2) / y * math.log(q2_max / q2_min)
        - 2.0 * (1.0 - y) / y * (1.0 - q2_min / q2_max))
    return dndy * 2.0 * w / s_minus_mp2


def integrated_flux(low, high, electron_energy, proton_energy, q2_max):
    grid = np.linspace(low, high, 500)
    values = np.array([photon_flux_dW(w, electron_energy, proton_energy, q2_max) for w in grid])
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    return float(integrate(values, grid))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def state_mask(data, settings):
    return np.abs(data["true_mass"] - settings["mass_gev"]) <= settings["true_mass_halfwidth_gev"]


def sample_branching_fraction(spec, settings):
    """Allow explicit campaign conventions while retaining archived global fallback."""
    value = spec["branching_fraction_ee"] if "branching_fraction_ee" in spec else settings["branching_fraction_ee"]
    return positive_branching_fraction(value)


def normalize_weights(data, spec, settings):
    """Normalize before any analysis Q2/W/reconstruction selection."""
    raw = np.asarray(data["weight"], float)
    if raw.ndim != 1 or not len(raw) or not np.isfinite(raw).all():
        raise ValueError("Input weights must be a nonempty finite one-dimensional array")
    if np.any(raw < 0):
        raise ValueError("This normalization diagnostic requires nonnegative input weights")
    state = state_mask(data, settings)
    target = float(spec["visible_xsec_fb"]) * float(settings["lumi_invfb"])
    if not np.isfinite(target) or target <= 0:
        raise ValueError("Visible cross section and luminosity must give a positive finite target")
    mode = spec["normalization_mode"]
    if mode == "truth_state":
        denominator = float(raw[state].sum())
        if denominator <= 0:
            raise ValueError("Truth-state weight sum before analysis cuts is nonpositive")
        scale = target / denominator
        weights = raw * scale
        constant = None
    elif mode == "all_generated":
        denominator = len(raw)
        constant = target / denominator
        weights = np.full(len(raw), constant, dtype=float)
        scale = None
    else:
        raise ValueError("Unknown normalization_mode: " + str(mode))
    audit = {
        "mode": mode, "normalization_source": spec["normalization_source"],
        "normalization_warning": spec.get("normalization_warning", ""),
        "normalization_before_analysis_cuts": True,
        "lumi_invfb": float(settings["lumi_invfb"]), "visible_xsec_fb": float(spec["visible_xsec_fb"]),
        "branching_fraction_ee": sample_branching_fraction(spec, settings),
        "branching_fraction_convention": "sample_override" if "branching_fraction_ee" in spec else "global_setting",
        "branching_fraction_uncertainty_included": False,
        "target_visible_event_sum": target, "generated_rows": len(raw),
        "truth_state_rows_before_analysis": int(state.sum()),
        "raw_weight_sum_all": float(raw.sum()), "raw_weight_sumw2_all": float(np.square(raw).sum()),
        "raw_weight_sum_truth_state": float(raw[state].sum()),
        "normalization_denominator": denominator, "raw_weight_scale": scale,
        "constant_weight_per_generated_row": constant,
        "normalized_weight_sum_all": float(weights.sum()),
        "normalized_weight_sumw2_all": float(np.square(weights).sum()),
        "normalized_weight_sum_truth_state": float(weights[state].sum()),
        "normalized_weight_sumw2_truth_state": float(np.square(weights[state]).sum()),
        "uncertainty_convention": ERROR_CONVENTION,
    }
    return weights, audit


def load_sample(spec, settings):
    path = Path(spec["root"]).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(str(path))
    with uproot.open(path) as source:
        if settings["tree"] not in source:
            raise ValueError("Missing tree " + settings["tree"])
        tree = source[settings["tree"]]
        branches = list(BRANCHES) + proton_branches()
        missing = sorted(set(branches) - set(tree.keys()))
        if missing:
            raise ValueError("Missing required branches: " + ", ".join(missing))
        data = {key: np.asarray(tree[key].array(library="np"), float) for key in branches}
    proton_tag_mask(data)
    if any(values.ndim != 1 or len(values) != len(data["weight"]) for values in data.values()):
        raise ValueError("Branches must be aligned scalar arrays")
    weights, audit = normalize_weights(data, spec, settings)
    audit.update({"proton_selection": proton_selection_metadata(), "input_path": str(path), "input_sha256": sha256(path), "input_size_bytes": path.stat().st_size,
                  "nonfinite_branch_rows": {k: int((~np.isfinite(v)).sum()) for k, v in data.items()},
                  "electron_energy": spec["electron_energy"], "proton_energy": spec["proton_energy"]})
    return {"spec": spec, "settings": settings, "data": data, "weights": weights, "normalization": audit}


def histogram(values, weights, mask, edges):
    mask = mask & np.isfinite(values)
    return np.histogram(values[mask], bins=edges, weights=weights[mask])[0].astype(float)


def compute_curve(sample, edges, q2_event_max, q2_flux_max):
    edges = np.asarray(edges, float)
    if len(edges) < 2 or not np.isfinite(edges).all() or not np.all(np.diff(edges) > 0):
        raise ValueError("W edges must be finite and strictly increasing")
    if not np.isfinite([q2_event_max, q2_flux_max]).all() or min(q2_event_max, q2_flux_max) <= 0:
        raise ValueError("Q2 event and flux maxima must be positive finite numbers")
    d, w, cfg, spec = sample["data"], sample["weights"], sample["settings"], sample["spec"]
    state = state_mask(d, cfg)
    q2 = np.isfinite(d["true_Q2"]) & (d["true_Q2"] >= 0) & (d["true_Q2"] <= q2_event_max)
    truth = state & q2
    # Match the framework's reconstructed selection, including migration and
    # possible other-state leakage into the reconstructed 1S mass window.
    reco = (q2 & proton_tag_mask(d) & (np.abs(d["reco_mass"] - cfg["mass_gev"]) <= cfg["reco_mass_halfwidth_gev"])
            & (np.abs(d["sigma_epz"] - 2.0 * spec["electron_energy"]) <= 0.2 * spec["electron_energy"]))
    t, r = histogram(d["true_W"], w, truth, edges), histogram(d["reco_W"], w, reco, edges)
    vt, vr = histogram(d["true_W"], w*w, truth, edges), histogram(d["reco_W"], w*w, reco, edges)
    nt = histogram(d["true_W"], np.ones(len(w)), truth, edges).astype(int)
    nr = histogram(d["reco_W"], np.ones(len(w)), reco, edges).astype(int)
    raw_t = histogram(d["true_W"], d["weight"], truth, edges)
    raw_r = histogram(d["reco_W"], d["weight"], reco, edges)
    flux = np.array([integrated_flux(a, b, spec["electron_energy"], spec["proton_energy"], q2_flux_max)
                     for a, b in zip(edges[:-1], edges[1:])])
    br_ee = sample_branching_fraction(spec, cfg)
    factor = np.divide(1.0, cfg["lumi_invfb"] * 1000.0 * br_ee * flux,
                       out=np.full_like(flux, np.nan), where=flux > 0)
    efficiency = np.divide(r, t, out=np.full_like(t, np.nan), where=t > 0)
    valid = (t > 0) & (flux > 0) & np.isfinite(efficiency) & (r > 0)
    supported = valid & (efficiency > cfg["min_reconstruction_ratio"])
    rows = []
    for i in range(len(t)):
        rows.append({"proton_selection": proton_selection_metadata()["mode"],
                     "sample": spec["key"], "W_low": edges[i], "W_high": edges[i+1],
                     "branching_fraction_ee": br_ee,
                     "W_center": (edges[i]+edges[i+1])/2, "q2_event_max": q2_event_max, "q2_flux_max": q2_flux_max,
                     "true_rows": nt[i], "reco_rows": nr[i], "true_raw_weight_sum": raw_t[i], "reco_raw_weight_sum": raw_r[i],
                     "true_visible_weight_sum": t[i], "true_visible_sumw2": vt[i],
                     "reco_visible_weight_sum": r[i], "reco_visible_sumw2": vr[i],
                     "reco_over_truth_matching_bins": efficiency[i], "flux_int": flux[i],
                     "sigma_gamma_p_truth_pb": t[i]*factor[i], "sigma_gamma_p_truth_mc_err_pb": math.sqrt(vt[i])*factor[i],
                     "truth_supported": bool(t[i] > 0 and flux[i] > 0),
                     "reco_correction_defined": bool(valid[i]), "reco_supported": bool(supported[i]),
                     "min_reconstruction_ratio_strict": cfg["min_reconstruction_ratio"],
                     "normalization_mode": spec["normalization_mode"], "visible_xsec_fb": spec["visible_xsec_fb"],
                     "normalization_warning": spec.get("normalization_warning", "")})
    selection = {
        "proton_selection": proton_selection_metadata(),
        "branching_fraction_ee": br_ee,
        "q2_event_max": q2_event_max, "q2_flux_max": q2_flux_max, "w_edges": edges.tolist(),
        "truth_state_rows_after_q2": int(truth.sum()),
        "truth_state_visible_sum_after_q2": float(w[truth].sum()),
        "truth_state_visible_sumw2_after_q2": float(np.square(w[truth]).sum()),
        "truth_visible_sum_in_w": float(t.sum()), "truth_visible_sumw2_in_w": float(vt.sum()),
        "reco_visible_sum_in_w": float(r.sum()), "reco_visible_sumw2_in_w": float(vr.sum()),
        "reco_epz_center": 2.0*spec["electron_energy"], "reco_epz_halfwidth": 0.2*spec["electron_energy"],
        "mass_gev": cfg["mass_gev"], "true_mass_halfwidth_gev": cfg["true_mass_halfwidth_gev"],
        "reco_mass_halfwidth_gev": cfg["reco_mass_halfwidth_gev"],
        "reco_selection_includes_truth_state_requirement": False,
        "normalization_recomputed_after_selection": False,
        "undefined_ratios": "NaN; never replaced with zero",
        "support_definition": "R(W_reco)/T(W_true) > threshold; includes migration and is not a response-row acceptance",
    }
    return {"sample": sample, "rows": rows, "selection": selection, "edges": edges,
            "q2_event_max": q2_event_max, "q2_flux_max": q2_flux_max}


def compare_curves(left, right):
    """Central-value ratio only; refuse mismatched bins or event/flux cuts."""
    if not np.array_equal(left["edges"], right["edges"]):
        raise ValueError("Cannot form a ratio with different W-bin edges")
    if any(left[key] != right[key] for key in ["q2_event_max", "q2_flux_max"]):
        raise ValueError("Cannot form a ratio with different Q2 selections or flux bounds")
    if left["selection"]["proton_selection"] != right["selection"]["proton_selection"]:
        raise ValueError("Cannot compare reconstruction support with different proton selections")
    rows = []
    for a, b in zip(left["rows"], right["rows"]):
        numerator, denominator = a["sigma_gamma_p_truth_pb"], b["sigma_gamma_p_truth_pb"]
        ok = a["truth_supported"] and b["truth_supported"] and denominator > 0
        rows.append({"W_low": a["W_low"], "W_high": a["W_high"], "W_center": a["W_center"],
                     "numerator_sample": a["sample"], "denominator_sample": b["sample"],
                     "ratio": numerator/denominator if ok else np.nan, "ratio_defined": bool(ok),
                     "both_reco_supported": a["reco_supported"] and b["reco_supported"],
                     "ratio_uncertainty": "not assigned; no independence assumption",
                     "same_input_sha256": left["sample"]["normalization"]["input_sha256"] == right["sample"]["normalization"]["input_sha256"]})
    return rows


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(data), indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_curve(curve, directory, name):
    write_csv(directory/(name+".csv"), curve["rows"])
    write_json(directory/(name+"_audit.json"), {"normalization": curve["sample"]["normalization"], "selection": curve["selection"]})


def plot_curves(curves, directory, stem, title, ratios=None):
    if ratios is not None:
        fig, (ax, ratio_ax) = plt.subplots(2, 1, figsize=(11, 9), sharex=True,
                                         gridspec_kw={"height_ratios": [3, 1]})
    else:
        fig, ax = plt.subplots(figsize=(11, 7.5))
        ratio_ax = None
    warnings = []
    for i, curve in enumerate(curves):
        rows, spec = curve["rows"], curve["sample"]["spec"]
        x = np.array([r["W_center"] for r in rows])
        dx = np.diff(curve["edges"])/2
        y = np.array([r["sigma_gamma_p_truth_pb"] for r in rows])
        dy = np.array([r["sigma_gamma_p_truth_mc_err_pb"] for r in rows])
        supported = np.array([r["reco_supported"] for r in rows])
        valid = np.isfinite(y) & (y > 0)
        color = "C"+str(i)
        label = (spec["label"] + "; Q²≤" + str(curve["q2_event_max"])
                 + "; flux Q²max=" + str(curve["q2_flux_max"]))
        ax.errorbar(x[valid], y[valid], xerr=dx[valid], yerr=dy[valid], fmt="o", color=color,
                    markerfacecolor="none", capsize=2, label=label)
        mask = valid & supported
        ax.scatter(x[mask], y[mask], color=color, marker="s", s=15, zorder=4)
        if spec.get("normalization_warning"):
            warnings.append(spec["normalization_warning"])
    ax.set_yscale("log")
    ax.set_ylabel("Truth-derived σ(γp → Υ(1S)p) [pb]")
    ax.set_title(title + "\nWeighted-MC normalization diagnostic; not an independent closure test", fontsize=12)
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="best")
    ax.text(0.02, 0.02, "Open circles: truth; small squares: bins passing the reconstruction-ratio threshold", transform=ax.transAxes, fontsize=8)
    if ratio_ax is not None:
        valid = np.array([r["ratio_defined"] for r in ratios])
        x = np.array([r["W_center"] for r in ratios])
        y = np.array([r["ratio"] for r in ratios])
        ratio_ax.plot(x[valid], y[valid], "o", color="black")
        ratio_ax.axhline(1, color="gray", linestyle="--")
        ratio_ax.set_ylabel("First / second")
        ratio_ax.grid(alpha=0.2)
        ratio_ax.set_xlabel("Wγp [GeV]")
    else:
        ax.set_xlabel("Wγp [GeV]")
    footer = "MC bars: √Σw² with normalization fixed. Ratios have no assigned independent error."
    if warnings:
        footer += "\n" + "\n".join(dict.fromkeys(warnings))
    fig.text(0.5, 0.01, footer, ha="center", va="bottom", fontsize=8)
    fig.tight_layout(rect=(0, 0.065 if warnings else 0.04, 1, 1))
    directory.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf"]:
        fig.savefig(directory/(stem+"."+suffix), dpi=180)
    plt.close(fig)


def make_edges(spec):
    count = int(spec["n_bins"])
    if count <= 0 or count != spec["n_bins"]:
        raise ValueError("n_bins must be a positive integer")
    return np.linspace(float(spec["w_min"]), float(spec["w_max"]), count+1)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).resolve().parents[2]/"config"/"w_samples.json"
    parser.add_argument("--config", type=Path, default=default)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--set-root", action="append", default=[], metavar="KEY=PATH")
    parser.add_argument("--only-comparisons", action="append", default=[], metavar="ID[,ID...]",
                        help="Run only these comparison IDs and their inputs; repeatable. Explicitly requested optional comparisons become required.")
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text())
    for override in args.set_root:
        key, separator, value = override.partition("=")
        if not separator or key not in config["samples"] or not value:
            parser.error("--set-root requires an existing sample key and nonempty path: KEY=PATH")
        config["samples"][key]["root"] = value
    requested = {name.strip() for value in args.only_comparisons for name in value.split(",") if name.strip()}
    if args.only_comparisons and not requested:
        parser.error("--only-comparisons needs at least one comparison ID")
    if requested:
        known = {item["id"] for item in config["comparisons"]}
        unknown = requested - known
        if unknown:
            parser.error("Unknown comparison IDs: " + ", ".join(sorted(unknown)))
        config["comparisons"] = [item for item in config["comparisons"] if item["id"] in requested]
        needed = {key for item in config["comparisons"] for key in item["samples"]}
        config["samples"] = {key: spec for key, spec in config["samples"].items() if key in needed}
        for item in config["comparisons"]:
            item["optional"] = False
        for spec in config["samples"].values():
            spec["optional"] = False
    # Prevent a partial rerun from leaving an older successful figure looking current.
    if args.out_dir.is_symlink() or (args.out_dir.exists() and not args.out_dir.is_dir()):
        parser.error("--out-dir must be a real directory, not a file or symlink")
    if args.out_dir.exists() and any(path.is_symlink() or not path.is_dir()
                                   for path in args.out_dir.rglob("*")):
        parser.error("--out-dir must contain no files or symlinks; empty directory trees are allowed")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "status": "RUNNING",
                "config_path": str(args.config.resolve()), "config_sha256": sha256(args.config),
                "script_sha256": sha256(__file__), "resolved_configuration": config,
                "requested_comparisons": sorted(requested) if requested else "all",
                "uncertainty_convention": ERROR_CONVENTION, "samples": {}, "comparisons": {},
                "notes": ["No source distributions are digitized or fitted to force agreement.",
                          "Truth-state normalization precedes all Q2, W, and reconstruction cuts.",
                          "A reconstruction support flag is a migration-inclusive matching-bin ratio, not acceptance.",
                          "Native Q2 comparisons are overlays only; ratios require identical event and flux Q2 bounds."]}
    loaded, native = {}, []
    for key, spec in config["samples"].items():
        spec["key"] = key
        try:
            sample = load_sample(spec, config)
            selection = spec["native"]
            curve = compute_curve(sample, make_edges(selection), selection["q2_event_max"], selection["q2_flux_max"])
            save_curve(curve, args.out_dir/"native", key)
            loaded[key] = sample
            native.append(curve)
            manifest["samples"][key] = {"status": "DONE", **sample["normalization"]}
            print("DONE sample", key)
        except (OSError, ValueError, KeyError) as error:
            status = "SKIPPED_OPTIONAL" if spec.get("optional") and isinstance(error, FileNotFoundError) else "BLOCKED"
            manifest["samples"][key] = {"status": status, "reason": str(error), "error_type": type(error).__name__}
            print(status + " sample", key, str(error))
        write_json(args.out_dir/"manifest.json", manifest)
    if native:
        plot_curves(native, args.out_dir/"native", "native_overlay", "Available samples with their native selections")
    for comparison in config["comparisons"]:
        name, keys = comparison["id"], comparison["samples"]
        missing = [key for key in keys if key not in loaded]
        if missing:
            status = ("SKIPPED_OPTIONAL" if comparison.get("optional") and
                      all(manifest["samples"][key]["status"] == "SKIPPED_OPTIONAL" for key in missing) else "BLOCKED")
            manifest["comparisons"][name] = {"status": status, "missing_samples": missing}
            print(status + " comparison", name, ", ".join(missing))
            continue
        try:
            curves = []
            for key in keys:
                q2 = loaded[key]["spec"]["native"] if comparison.get("use_native_q2") else comparison
                curve = compute_curve(loaded[key], make_edges(comparison), q2["q2_event_max"], q2["q2_flux_max"])
                save_curve(curve, args.out_dir/"comparisons"/name, key)
                curves.append(curve)
            ratios = None if comparison.get("use_native_q2") else compare_curves(curves[0], curves[1])
            if ratios is not None:
                write_csv(args.out_dir/"comparisons"/name/"ratio.csv", ratios)
            plot_curves(curves, args.out_dir/"comparisons"/name, "comparison", name.replace("_", " "), ratios)
            manifest["comparisons"][name] = {"status": "DONE", "samples": keys,
                                             "ratio_produced": ratios is not None,
                                             "selections": {key: c["selection"] for key, c in zip(keys, curves)}}
            print("DONE comparison", name)
        except (OSError, ValueError, KeyError) as error:
            manifest["comparisons"][name] = {"status": "FAILED", "reason": str(error), "error_type": type(error).__name__}
            print("FAILED comparison", name, str(error))
        write_json(args.out_dir/"manifest.json", manifest)
    incomplete = any(r["status"] not in ["DONE", "SKIPPED_OPTIONAL"]
                     for group in ["samples", "comparisons"] for r in manifest[group].values())
    manifest["status"] = "PARTIAL" if incomplete else "COMPLETE"
    write_json(args.out_dir/"manifest.json", manifest)
    if not incomplete:
        # The framework must list this marker as an expected output. Existing
        # CSVs/figures from a PARTIAL run are not evidence of successful work.
        write_json(args.out_dir/"comparison_complete.json", {
            "status": "COMPLETE",
            "manifest_sha256": sha256(args.out_dir/"manifest.json"),
            "script_sha256": manifest["script_sha256"],
            "required_samples_and_comparisons_completed": True,
            "optional_skips": {group: [key for key, record in manifest[group].items()
                                        if record["status"] == "SKIPPED_OPTIONAL"]
                               for group in ["samples", "comparisons"]},
        })
    print("Manifest:", args.out_dir/"manifest.json")
    return 1 if incomplete else 0


if __name__ == "__main__":
    sys.exit(main())
