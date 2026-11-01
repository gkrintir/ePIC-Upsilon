#!/usr/bin/env python3
"""Compare two UpsilonTree samples without assuming cross-section normalization."""

from __future__ import annotations

import argparse

# Explicit reconstructed-proton category; truth comparisons stay unchanged.
import sys as _proton_sys
from pathlib import Path as _ProtonPath
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
REQUIRED = ("true_mass", "reco_mass", "reco_W", "sigma_epz")
# The final field excludes nonphysical negative values and negative sentinels.
TRUTH = [
    ("true_mass", "Truth dilepton mass [GeV]", True),
    ("true_W", "Truth W [GeV]", True),
    ("true_Q2", "Truth Q² [GeV²]", True),
    ("true_t", "Truth |t| [GeV²]", True),
    ("true_pT", "Truth Upsilon pT [GeV]", True),
    ("true_y", "Truth Upsilon rapidity", False),
    ("true_ePlus_eta", "Truth positron eta", False),
    ("true_decay_eMinus_eta", "Truth decay electron eta", False),
    ("true_scat_eMinus_E", "Truth scattered electron energy [GeV]", True),
]
RECO = [
    ("reco_mass", "Reconstructed dilepton mass [GeV]", True),
    ("reco_W", "Reconstructed W [GeV]", True),
    ("reco_pT", "Reconstructed Upsilon pT [GeV]", True),
    ("reco_y", "Reconstructed Upsilon rapidity", False),
    ("sigma_epz", "Sum(E − pz) [GeV]", True),
    ("reco_scat_eMinus_E", "Reconstructed scattered electron energy [GeV]", True),
    ("reco_ePlus_E", "Reconstructed positron energy [GeV]", True),
    ("reco_decay_eMinus_E", "Reconstructed decay electron energy [GeV]", True),
    ("reco_t_eX", "Reconstructed |t|, eX [GeV²]", True),
    ("reco_t_BABE", "Reconstructed |t|, BABE [GeV²]", True),
    ("reco_t_eXBABE", "Reconstructed |t|, eXBABE [GeV²]", True),
    ("reco_t_proj", "Reconstructed |t|, valid projective [GeV²]", True),
]


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference-file", required=True, type=Path)
    p.add_argument("--candidate-file", required=True, type=Path)
    p.add_argument("--reference-label", default="10x100 GeV")
    p.add_argument("--candidate-label", default="9x130 GeV")
    p.add_argument("--reference-electron-energy", default=10.0, type=float)
    p.add_argument("--candidate-electron-energy", default=9.0, type=float)
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--state", choices=MASSES, default="1S")
    p.add_argument("--weight-branch", help="Optional explicit shape weights; defaults to unweighted")
    p.add_argument("--w-min", type=float, help="Optional lower truth-W bound for both samples")
    p.add_argument("--w-max", type=float, help="Optional upper truth-W bound for both samples")
    p.add_argument("--bins", default=60, type=int)
    p.add_argument("--out-dir", required=True, type=Path)
    a = p.parse_args()
    if a.bins < 1 or min(a.reference_electron_energy, a.candidate_electron_energy) <= 0:
        p.error("bins and both electron energies must be positive")
    if a.w_min is not None and a.w_max is not None and a.w_min >= a.w_max:
        p.error("--w-min must be smaller than --w-max")
    return a


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path.resolve()), "size_bytes": path.stat().st_size,
            "sha256": digest.hexdigest()}


def load_sample(path, label, energy, role, args):
    path = path.expanduser()
    needed = set(REQUIRED) | set(proton_branches())
    if args.w_min is not None or args.w_max is not None:
        needed.add("true_W")
    if args.weight_branch:
        needed.add(args.weight_branch)
    requested = needed | {v[0] for v in TRUTH + RECO} | {"weight_proj"}
    with uproot.open(path) as root_file:
        if args.tree not in root_file:
            raise ValueError(f"{role} ({label}): missing tree {args.tree!r}")
        tree = root_file[args.tree]
        available = set(tree.keys())
        missing = sorted(needed - available)
        if missing:
            raise ValueError(f"{role} ({label}): missing required branches {missing}")
        data = {}
        for key in sorted(requested & available):
            values = np.asarray(tree[key].array(library="np"), dtype=float)
            if values.ndim != 1:
                raise ValueError(f"{role} ({label}): branch {key} must be scalar per row")
            data[key] = values
        n_rows = int(tree.num_entries)
    mass = MASSES[args.state]
    denominator = np.isfinite(data["true_mass"]) & (np.abs(data["true_mass"] - mass) <= 0.2)
    if args.w_min is not None or args.w_max is not None:
        denominator &= np.isfinite(data["true_W"]) & (data["true_W"] >= 0)
        if args.w_min is not None:
            denominator &= data["true_W"] >= args.w_min
        if args.w_max is not None:
            denominator &= data["true_W"] <= args.w_max
    reco = denominator.copy()
    for branch in ("reco_mass", "reco_W", "sigma_epz"):
        reco &= np.isfinite(data[branch])
    reco &= (data["reco_mass"] >= 8) & (data["reco_mass"] <= 12)
    reco &= (data["reco_W"] >= 0) & (np.abs(data["sigma_epz"] - 2 * energy) <= 0.2 * energy)
    reco &= proton_tag_mask(data)
    if not np.any(denominator):
        raise ValueError(f"{role} ({label}): empty {args.state} tuple denominator after truth selections")
    if not np.any(reco):
        raise ValueError(f"{role} ({label}): empty reconstructed selection for {args.state}")
    weights = data[args.weight_branch] if args.weight_branch else np.ones(n_rows)
    if args.weight_branch and np.any(weights < 0):
        raise ValueError(f"{role} ({label}): negative weights in {args.weight_branch!r} are unsupported")
    n_gen, n_reco = int(denominator.sum()), int(reco.sum())
    summary = {"proton_selection": proton_selection_metadata(), "sample": role, "label": label, "tree_rows": n_rows,
               "truth_selected_tuple_rows": n_gen, "reco_selected_tuple_rows": n_reco,
               "unweighted_reco_fraction_of_tuple_denominator": n_reco / n_gen,
               "electron_energy_GeV": energy, "sigma_epz_center_GeV": 2 * energy,
               "sigma_epz_halfwidth_GeV": 0.2 * energy}
    metadata = {**fingerprint(path), **summary, "missing_optional_branches": sorted(requested - available),
                "nonfinite_weight_rows": int(np.count_nonzero(~np.isfinite(weights)))}
    return {"data": data, "label": label, "role": role, "weights": weights,
            "truth_mask": denominator, "reco_mask": reco, "summary": summary, "metadata": metadata}


def plot_family(family, specifications, samples, args, bin_rows, statuses):
    prepared = []
    for branch, xlabel, nonnegative in specifications:
        status = {"family": family, "branch": branch}
        required_for_plot = {branch}
        if branch == "reco_t_proj":
            required_for_plot |= {"weight_proj", "reco_t_eXBABE"}
        absent = [s["role"] for s in samples if not required_for_plot.issubset(s["data"])]
        if absent:
            statuses.append({**status, "status": "skipped", "reason": "missing branch", "samples": absent})
            continue
        arrays = []
        for s in samples:
            values = s["data"][branch]
            valid = s[family + "_mask"] & np.isfinite(values) & np.isfinite(s["weights"])
            if branch == "reco_t_proj":
                valid &= np.isfinite(s["data"]["weight_proj"]) & (s["data"]["weight_proj"] > 0)
                valid &= np.isfinite(s["data"]["reco_t_eXBABE"]) & (s["data"]["reco_t_eXBABE"] > 0)
            # Negative physical magnitudes and |eta/y|>=900 are tuple sentinels.
            valid &= values >= 0 if nonnegative else np.abs(values) < 900
            arrays.append((values[valid], s["weights"][valid]))
        if any(len(v) == 0 or float(w.sum()) <= 0 for v, w in arrays):
            statuses.append({**status, "status": "skipped", "reason": "empty valid values or zero weight sum"})
            continue
        lo = min(float(v.min()) for v, _ in arrays)
        hi = max(float(v.max()) for v, _ in arrays)
        if branch == "true_Q2":
            lo, hi = 0.0, 1.0
        else:
            if hi == lo:
                span = max(abs(lo) * 0.01, 0.5)
                lo, hi = lo - span, hi + span
        edges = np.linspace(lo, hi, args.bins + 1)
        density_values = []
        for s, (values, weights) in zip(samples, arrays):
            counts, _ = np.histogram(values, bins=edges)
            sums, _ = np.histogram(values, bins=edges, weights=weights)
            density = sums / (float(sums.sum()) * np.diff(edges))
            density_values.append(density)
            for i in range(args.bins):
                bin_rows.append({"family": family, "branch": branch, "sample": s["role"],
                                 "label": s["label"], "bin_low": float(edges[i]),
                                 "bin_high": float(edges[i + 1]), "row_count": int(counts[i]),
                                 "weight_sum": float(sums[i]), "unit_density": float(density[i])})
        statuses.append({**status, "status": "plotted", "range": [lo, hi],
                         "valid_rows": {s["role"]: len(v) for s, (v, _) in zip(samples, arrays)}})
        prepared.append((xlabel, edges, density_values))
    if not prepared:
        raise ValueError(f"No comparable {family} variables survived validation")
    ncols = 3
    nrows = (len(prepared) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.5 * nrows), squeeze=False)
    colors = ("#2473B8", "#C35B24")
    for ax, (xlabel, edges, densities) in zip(axes.flat, prepared):
        for s, density, color in zip(samples, densities, colors):
            ax.stairs(density, edges, label=s["label"], color=color, linewidth=1.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Unit probability density")
        ax.ticklabel_format(axis="both", style="sci", scilimits=(-3, 4), useOffset=False)
        if "|t|" or "Q²" in xlabel:
            ax.set_yscale("log")
        ax.grid(alpha=0.15)
        ax.legend(fontsize=8, frameon=False)
    for ax in list(axes.flat)[len(prepared):]:
        ax.set_visible(False)
    mode = f"weighted by {args.weight_branch}" if args.weight_branch else "unweighted"
    fig.suptitle(f"Upsilon({args.state}) {family} shapes · {mode}", fontsize=16)
    foot = "Each curve normalized separately over its finite, valid values. No cross-section normalization."
    if family == "reco":
        foot += "\nReconstruction: mass 8–12 GeV, W ≥ 0, |Sum(E−pz)−2Ee| ≤ 0.2Ee."
    if args.w_min is not None or args.w_max is not None:
        foot += f"\nCommon truth W selection: {args.w_min if args.w_min is not None else '0'} to {args.w_max if args.w_max is not None else 'unbounded'} GeV."
    fig.text(0.5, 0.012, foot, ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, 0.085, 1, 0.96))
    for extension in ("png", "pdf"):
        fig.savefig(args.out_dir / f"{family}_shapes.{extension}", dpi=180)
    plt.close(fig)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = arguments()
    args.out_dir = args.out_dir.expanduser()
    samples = [load_sample(args.reference_file, args.reference_label, args.reference_electron_energy, "reference", args),
               load_sample(args.candidate_file, args.candidate_label, args.candidate_electron_energy, "candidate", args)]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    bin_rows, statuses = [], []
    for family, specifications in (("truth", TRUTH), ("reco", RECO)):
        plot_family(family, specifications, samples, args, bin_rows, statuses)
    write_csv(args.out_dir / "selection_summary.csv", [s["summary"] for s in samples])
    write_csv(args.out_dir / "bin_values.csv", bin_rows)
    manifest = {
        "schema_version": 1, "tree": args.tree, "state": args.state,
        "truth_mass_center_GeV": MASSES[args.state], "truth_mass_halfwidth_GeV": 0.2,
        "common_truth_w_min_GeV": args.w_min, "common_truth_w_max_GeV": args.w_max,
        "shape_weight_branch": args.weight_branch, "bins": args.bins,
        "density_definition": "bin weight sum / (valid selected weight sum * bin width)",
        "selection_fraction_definition": "unweighted reconstructed selected rows / truth-selected stored tuple rows",
        "denominator_caveat": "Tuple selection fraction, not a detector efficiency relative to all generated events.",
        "reco_selection": "8 <= reco_mass <= 12; reco_W >= 0; abs(sigma_epz - 2*Ee) <= 0.2*Ee; finite required reco fields",
        "truth_selection": "abs(true_mass - state_mass) <= 0.2; optional inclusive common truth-W bounds",
        "sentinel_policy": "Per variable: finite values; nonnegative magnitudes; abs(eta or rapidity) < 900; projective t also requires finite positive weight_proj and reco_t_eXBABE",
        "samples": [s["metadata"] for s in samples], "variables": statuses,
        "software": {"numpy": np.__version__, "matplotlib": matplotlib.__version__, "uproot": uproot.__version__},
        "outputs": ["truth_shapes.png", "truth_shapes.pdf", "reco_shapes.png", "reco_shapes.pdf",
                    "selection_summary.csv", "bin_values.csv", "comparison_manifest.json"],
    }
    (args.out_dir / "comparison_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for s in samples:
        m = s["summary"]
        print(f"{s['label']}: {m['reco_selected_tuple_rows']}/{m['truth_selected_tuple_rows']} tuple rows pass reco selection ({m['unweighted_reco_fraction_of_tuple_denominator']:.4%})")
    for item in statuses:
        if item["status"] == "skipped":
            print(f"SKIP {item['family']} {item['branch']}: {item['reason']}")
    print(f"Outputs: {args.out_dir.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        raise SystemExit(f"Comparison failed: {error}")
