#!/usr/bin/env python3
"""Rebuild the ten report TH2 observables for a configured latest tuple.

The original histogram definitions and bin edges are reused. Beam exclusivity,
truth state and Q2 selections are explicit. Projected mode preserves the supplied
weights without further normalization; shape mode is an explicit optional view.
No comparison to a historical stored histogram is claimed. Every ROOT histogram
retains underflow and overflow bins, and NPZ stores counts, sumw, sumw2 and shape.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

# Reconstructed categories use the shared producer-tag contract only.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH2x

from reconstruct_report_histograms_from_tuple import SPECS, required_branches, npz_name

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reporting"))
from note_style import set_style, finish
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
LABELS = {
    "true_decay_eMinus_eta": r"True decay $e^-$ $\eta$",
    "true_ePlus_eta": r"True decay $e^+$ $\eta$",
    "true_decay_eMinus_pT": r"True decay $e^-$ $p_T$ [GeV]",
    "true_ePlus_pT": r"True decay $e^+$ $p_T$ [GeV]",
    "true_y": r"True $\Upsilon$ rapidity $y$",
    "true_W": r"True $W_{\gamma p}$ [GeV]",
    "true_t": r"True $|t|$ [GeV$^2$]",
    "reco_nTracks_Central": "Central track count",
    "reco_nTracks_Tagger": "Low-$Q^2$ tagger track count",
    "reco_decay_eMinus_eta": r"Reco decay $e^-$ $\eta$",
    "reco_ePlus_eta": r"Reco decay $e^+$ $\eta$",
    "reco_decay_eMinus_pT": r"Reco decay $e^-$ $p_T$ [GeV]",
    "reco_ePlus_pT": r"Reco decay $e^+$ $p_T$ [GeV]",
    "reco_t_eX": r"Reco $|t|_{eX}$ [GeV$^2$]",
    "reco_t_eXBABE": r"Reco $|t|_{eXBABE}$ [GeV$^2$]",
}
TITLES = {
    "GenLevel/h2_Gen_DecayLeps_Eta_vs_Eta": "Generated decay-lepton pseudorapidity",
    "GenLevel/h2_Gen_DecayLeps_pT_vs_pT": "Generated decay-lepton transverse momenta",
    "GenLevel/h2_Gen_DecayLeps_Eta_vs_pT": "Generated decay-lepton kinematics",
    "GenLevel/h2_Gen_W_vs_Rapidity": r"Generated $W_{\gamma p}$ and rapidity",
    "GenLevel/h2_Gen_t_vs_W": r"Generated $|t|$ and $W_{\gamma p}$",
    "h2_nTracks_Central_vs_Tagger": "Central and tagger track multiplicities",
    "h2_eta_ePlus_eMinus": "Reconstructed decay-lepton pseudorapidity",
    "h2_pT_ePlus_eMinus": "Reconstructed decay-lepton transverse momenta",
    "h2_Reco_DecayLeps_Eta_vs_pT": "Reconstructed decay-lepton kinematics",
    "h2_t_Correlation": "Reconstructed momentum-transfer estimators",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def masks_for(arrays, args):
    weight = arrays[args.weight_branch]
    finite_weight = np.isfinite(weight) & (weight >= 0)
    if not finite_weight.all():
        raise ValueError("TH2 reconstruction requires finite, nonnegative input weights")
    truth = np.isfinite(arrays["true_mass"]) & (
        np.abs(arrays["true_mass"] - MASSES[args.state]) <= args.true_mass_window
    )
    truth &= np.isfinite(arrays["true_Q2"]) & (arrays["true_Q2"] >= 0)
    if args.q2_cut is not None:
        truth &= arrays["true_Q2"] <= args.q2_cut
    if args.w_min is not None:
        truth &= np.isfinite(arrays["true_W"]) & (arrays["true_W"] >= args.w_min)
    if args.w_max is not None:
        truth &= np.isfinite(arrays["true_W"]) & (arrays["true_W"] <= args.w_max)
    reco = truth & np.isfinite(arrays["reco_mass"]) & (
        (arrays["reco_mass"] >= 8.0) & (arrays["reco_mass"] <= 12.0)
    )
    reco &= np.isfinite(arrays["sigma_epz"]) & (
        np.abs(arrays["sigma_epz"] - 2 * args.electron_energy) <= 0.2 * args.electron_energy
    )
    reco &= proton_tag_mask(arrays)
    # Retain the existing estimator-validity condition; it is not the proton tag.
    tagged = reco & np.isfinite(arrays["reco_t_BABE"]) & (arrays["reco_t_BABE"] > 0)
    return {"all": finite_weight, "truth": truth, "reco": reco, "tagged": tagged}


def fill_with_flow(spec, arrays, mask, weights):
    xe, ye = np.asarray(spec.x_edges), np.asarray(spec.y_edges)
    counts = np.zeros((len(xe) + 1, len(ye) + 1), dtype=np.int64)
    sumw = np.zeros(counts.shape, dtype=np.float64)
    sumw2 = np.zeros(counts.shape, dtype=np.float64)
    excluded_nonfinite, sentinel_fills = 0, 0
    for xname, yname in spec.pairs:
        x, y = arrays[xname], arrays[yname]
        finite = np.isfinite(x) & np.isfinite(y)
        excluded_nonfinite += int(np.count_nonzero(mask & ~finite))
        keep = mask & finite
        x, y, w = x[keep], y[keep], weights[keep]
        sentinel_fills += int(np.count_nonzero((x <= -900) | (y <= -900)))
        # ROOT convention: a value on the final upper edge is overflow.
        ix, iy = np.searchsorted(xe, x, side="right"), np.searchsorted(ye, y, side="right")
        np.add.at(counts, (ix, iy), 1)
        np.add.at(sumw, (ix, iy), w)
        np.add.at(sumw2, (ix, iy), w * w)
    return counts, sumw, sumw2, excluded_nonfinite, sentinel_fills


def axis_labels(spec):
    if len(spec.pairs) == 2:
        prefix = "True" if spec.selection == "truth" else "Reco"
        return prefix + r" decay-lepton $p_T$ [GeV]", prefix + r" decay-lepton $\eta$"
    x, y = spec.pairs[0]
    return LABELS[x], LABELS[y]


def root_histogram(spec, values, variances, entries, xlabel, ylabel):
    xe, ye = np.asarray(spec.x_edges), np.asarray(spec.y_edges)
    xx, yy = np.meshgrid((xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2, indexing="ij")
    v, var = values[1:-1, 1:-1], variances[1:-1, 1:-1]
    return to_TH2x(
        fName=spec.key.rsplit("/", 1)[-1], fTitle=TITLES[spec.key],
        data=values.T.ravel(), fEntries=float(entries),
        fTsumw=float(v.sum()), fTsumw2=float(var.sum()),
        fTsumwx=float((v * xx).sum()), fTsumwx2=float((v * xx * xx).sum()),
        fTsumwy=float((v * yy).sum()), fTsumwy2=float((v * yy * yy).sum()),
        fTsumwxy=float((v * xx * yy).sum()), fSumw2=variances.T.ravel(),
        fXaxis=to_TAxis("xaxis", xlabel, len(xe) - 1, xe[0], xe[-1], xe),
        fYaxis=to_TAxis("yaxis", ylabel, len(ye) - 1, ye[0], ye[-1], ye),
    )


def render_map(spec, values, record, args):
    regular = values[1:-1, 1:-1]
    positive = regular[regular > 0]
    if positive.size:
        lo, hi = float(positive.min()), float(positive.max())
        norm = LogNorm(vmin=lo, vmax=max(hi, lo * 1.01))
        shown = np.ma.masked_less_equal(regular, 0)
    else:
        norm, shown = Normalize(vmin=0, vmax=1), regular
    fig, ax = plt.subplots(figsize=(8.7, 7.4))
    xlabel, ylabel = axis_labels(spec)
    image = ax.pcolormesh(spec.x_edges, spec.y_edges, shown.T, shading="auto", cmap="viridis", norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    label = {"projected": "Expected entries per bin", "shape": "Fraction of all entries per bin", "unit": "MC entries per bin"}[args.normalization]
    fig.colorbar(image, ax=ax, label=label)
    selection_labels = {"all": "all tuple rows", "truth": args.state + ", truth",
                        "reco": args.state + ", reconstructed selection",
                        "tagged": args.state + ", reconstructed selection; valid BABE estimator"}
    selected = selection_labels[spec.selection]
    if spec.selection in ("reco", "tagged"):
        selected += " | " + proton_selection_metadata()["label"]
    q2 = "" if spec.selection == "all" or args.q2_cut is None else rf", $Q^2\leq {args.q2_cut:g}$ GeV$^2$"
    subtitle = args.energy_label + " | " + selected + q2
    if args.normalization == "projected":
        subtitle += rf" | $\mathcal{{L}}={args.lumi_invfb:g}$ fb$^{{-1}}$"
    subtitle += "\n" + f"MC fills in view {record['visible_fills']:,}; outside view {record['flow_fills']:,}"
    return finish(fig, args.plot_dir, spec.key.replace("/", "_"), TITLES[spec.key], subtitle)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--state", choices=MASSES, default="1S")
    parser.add_argument("--true-mass-window", type=float, default=0.2)
    parser.add_argument("--electron-energy", type=float, required=True)
    parser.add_argument("--energy-label", default="")
    parser.add_argument("--q2-cut", type=float)
    parser.add_argument("--w-min", type=float)
    parser.add_argument("--w-max", type=float)
    parser.add_argument("--weight-branch", default="weight")
    parser.add_argument("--normalization", choices=["projected", "shape", "unit"], default="projected")
    parser.add_argument("--lumi-invfb", type=float, default=10.0, help="Label only; weights are not rescaled")
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--out-npz", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--plot-dir", type=Path)
    args = parser.parse_args()
    if args.electron_energy <= 0 or args.true_mass_window <= 0 or args.lumi_invfb <= 0:
        parser.error("Beam energy, truth mass window, and luminosity must be positive")
    if args.q2_cut is not None and args.q2_cut < 0:
        parser.error("--q2-cut must be nonnegative")
    if args.w_min is not None and args.w_max is not None and args.w_min >= args.w_max:
        parser.error("--w-min must be below --w-max")
    source = args.input.expanduser().resolve()
    destinations = [args.out_root.expanduser().resolve(), args.out_npz.expanduser().resolve(), args.out_json.expanduser().resolve()]
    if source in destinations or len(set(destinations)) != len(destinations):
        parser.error("Outputs must be separate from the input and from each other")
    branches = sorted(set(required_branches()) | set(proton_branches()) | {"true_Q2", "true_W", args.weight_branch})
    with uproot.open(source) as root_file:
        tree = root_file[args.tree]
        missing = sorted(set(branches) - set(tree.keys()))
        if missing:
            raise ValueError("Missing required branches: " + ", ".join(missing))
        arrays = tree.arrays(branches, library="np")
        normalization_metadata = json.loads(str(root_file["framework_normalization"])) if "framework_normalization" in root_file else None
    if normalization_metadata is not None and args.normalization == "projected":
        if not np.isclose(float(normalization_metadata["luminosity_invfb"]), args.lumi_invfb):
            raise ValueError("Plot luminosity differs from input normalization; regenerate normalized input")
    for path in destinations:
        path.parent.mkdir(parents=True, exist_ok=True)
    if args.plot_dir:
        set_style()
    masks = masks_for(arrays, args)
    weights = arrays[args.weight_branch]
    records, payload = [], {}
    report = {
        "input": str(source), "sha256": sha256(source), "tree": args.tree,
        "input_normalization": normalization_metadata,
        "proton_selection": proton_selection_metadata(),
        "proton_selection_scope": "Reconstructed reco/tagged masks only; all/truth masks and supplied weights unchanged",
        "histogram_definition_source": "reconstruct_report_histograms_from_tuple.py:SPECS",
        "historical_reference_compared": False,
        "normalization": args.normalization,
        "projected_mode": "Uses supplied weights without renormalization after any selection or axis window",
        "shape_mode": "Bin sumw divided by all finite-coordinate sumw, including underflow and overflow",
        "weight_branch": args.weight_branch, "state": args.state,
        "electron_energy_GeV": args.electron_energy,
        "q2_cut_GeV2": args.q2_cut, "truth_w_min_GeV": args.w_min, "truth_w_max_GeV": args.w_max,
        "selection_definitions": {
            "all": "All tuple rows; no state, Q2, W, or reconstruction selection",
            "truth": "Configured truth-state mass window, finite nonnegative Q2 and optional Q2/true-W bounds",
            "reco": "truth plus 8 <= reco_mass <= 12, abs(sigma_epz - 2*Ee) <= 0.2*Ee, and the declared proton-selection policy",
            "tagged": "reco plus finite reco_t_BABE > 0 (retained estimator-validity condition, not a substitute proton tag)",
        },
        "selected_rows": {name: int(mask.sum()) for name, mask in masks.items()},
        "coordinate_policy": "Nonfinite coordinate fills excluded and counted; finite sentinel values retained in ROOT flow bins",
        "histograms": records,
    }
    with uproot.recreate(destinations[0]) as output:
        for spec in SPECS:
            counts, sumw, sumw2, nonfinite, sentinel = fill_with_flow(spec, arrays, masks[spec.selection], weights)
            total = float(sumw.sum())
            shape = sumw / total if total > 0 else np.zeros_like(sumw)
            if args.normalization == "shape":
                values, variances = shape, sumw2 / total**2 if total > 0 else np.zeros_like(sumw2)
            elif args.normalization == "unit":
                values, variances = counts.astype(float), counts.astype(float)
            else:
                values, variances = sumw, sumw2
            key = npz_name(spec.key)
            payload.update({key + "__counts_flow": counts, key + "__sumw_flow": sumw,
                            key + "__sumw2_flow": sumw2, key + "__shape_flow": shape,
                            key + "__x_edges": np.asarray(spec.x_edges), key + "__y_edges": np.asarray(spec.y_edges)})
            record = {
                "key": spec.key, "selection": spec.selection, "pairs": [list(p) for p in spec.pairs],
                "proton_selection_applied": spec.selection in ("reco", "tagged") and proton_selection_metadata()["mode"] == "required",
                "fills": int(counts.sum()), "visible_fills": int(counts[1:-1, 1:-1].sum()),
                "flow_fills": int(counts.sum() - counts[1:-1, 1:-1].sum()),
                "nonfinite_coordinate_fills_excluded": nonfinite, "sentinel_coordinate_fills": sentinel,
                "sumw": total, "sumw_visible": float(sumw[1:-1, 1:-1].sum()),
                "sumw_flow": float(total - sumw[1:-1, 1:-1].sum()),
                "report_use": spec.report_use,
            }
            output[spec.key] = root_histogram(spec, values, variances, record["fills"], *axis_labels(spec))
            if args.plot_dir:
                record["plots"] = render_map(spec, values, record, args)
            records.append(record)
        output["latest_report_histogram_provenance"] = json.dumps(report, sort_keys=True)
    np.savez_compressed(destinations[1], **payload)
    destinations[2].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Reconstructed {len(records)} TH2 observables in {args.normalization} mode; no historical reference comparison")
    print(f"ROOT: {destinations[0]}")
    print(f"Audit: {destinations[2]}")


if __name__ == "__main__":
    main()
