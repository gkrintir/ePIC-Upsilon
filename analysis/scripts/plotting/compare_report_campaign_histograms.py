#!/usr/bin/env python3
"""Render the report-era TH1 and TH2 campaign histogram comparisons.

The default TH1 style preserves raw counts and the historical ratio layout.
Unlike the retained notebook cell, this implementation counts only files that
were actually written and also handles all TH2 histograms.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import uproot


def safe_name(key: str) -> str:
    return key.replace("/", "_").replace(" ", "_")


def values_for_mode(values: np.ndarray, mode: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if mode == "shape" and values.sum() > 0:
        return values / values.sum()
    return values


def plot_th1(
    old_hist,
    new_hist,
    key: str,
    output: Path,
    old_label: str,
    new_label: str,
    mode: str,
) -> bool:
    old_values, edges = old_hist.to_numpy(flow=False)
    new_values, new_edges = new_hist.to_numpy(flow=False)
    if not np.array_equal(edges, new_edges):
        raise ValueError(f"TH1 axis mismatch for {key}")
    old_values = values_for_mode(old_values, mode)
    new_values = values_for_mode(new_values, mode)
    if old_values.sum() == 0 and new_values.sum() == 0:
        return False

    centers = (edges[:-1] + edges[1:]) / 2
    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(8, 6),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )
    plt.subplots_adjust(hspace=0.05)
    ax1.hist(
        centers,
        bins=edges,
        weights=old_values,
        histtype="step",
        label=old_label,
        color="blue",
        linewidth=1.5,
    )
    ax1.hist(
        centers,
        bins=edges,
        weights=new_values,
        histtype="step",
        label=new_label,
        color="red",
        linestyle="--",
        linewidth=1.5,
    )
    ax1.set_ylabel("Events" if mode == "raw" else "Normalized events")
    ax1.legend()
    ax1.set_title(key)

    ratio = np.full_like(old_values, np.nan, dtype=float)
    mask = old_values != 0
    ratio[mask] = new_values[mask] / old_values[mask]
    ax2.plot(centers, ratio, "ko", markersize=3)
    ax2.axhline(1.0, color="gray", linestyle="--")
    ax2.set_ylabel("Ratio (F2 / F1)")
    ax2.set_xlabel("Observable Value")
    valid_ratio = ratio[np.isfinite(ratio)]
    if len(valid_ratio) > 0:
        rmin, rmax = np.min(valid_ratio), np.max(valid_ratio)
        if rmax - rmin < 0.1:
            ax2.set_ylim(0.9, 1.1)
        else:
            ax2.set_ylim(max(0.0, rmin - 0.2), min(3.0, rmax + 0.2))

    fig.savefig(output, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return True


def positive_norm(old_values: np.ndarray, new_values: np.ndarray):
    positive = np.concatenate((old_values[old_values > 0], new_values[new_values > 0]))
    if positive.size == 0:
        return None
    low = max(float(positive.min()), 1.0)
    high = float(positive.max())
    return colors.LogNorm(vmin=low, vmax=max(low, high))


def plot_th2(
    old_hist,
    new_hist,
    key: str,
    output: Path,
    old_label: str,
    new_label: str,
) -> bool:
    old_values, x_edges, y_edges = old_hist.to_numpy(flow=False)
    new_values, new_x_edges, new_y_edges = new_hist.to_numpy(flow=False)
    if not np.array_equal(x_edges, new_x_edges) or not np.array_equal(y_edges, new_y_edges):
        raise ValueError(f"TH2 axis mismatch for {key}")
    if old_values.sum() == 0 and new_values.sum() == 0:
        return False

    ratio = np.divide(
        new_values,
        old_values,
        out=np.full_like(new_values, np.nan, dtype=float),
        where=old_values != 0,
    )
    norm = positive_norm(old_values, new_values)
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    image = None
    for ax, values, label in zip(axes[:2], (old_values, new_values), (old_label, new_label)):
        image = ax.pcolormesh(x_edges, y_edges, values.T, shading="auto", norm=norm, cmap="viridis")
        ax.set_title(label)
        ax.set_xlabel(old_hist.axes[0].member("fTitle") or "x")
        ax.set_ylabel(old_hist.axes[1].member("fTitle") or "y")
    if image is not None:
        fig.colorbar(image, ax=axes[:2], label="Events")

    ratio_image = axes[2].pcolormesh(
        x_edges,
        y_edges,
        ratio.T,
        shading="auto",
        cmap="coolwarm",
        vmin=0.5,
        vmax=1.5,
    )
    axes[2].set_title(f"{new_label} / {old_label}")
    axes[2].set_xlabel(old_hist.axes[0].member("fTitle") or "x")
    axes[2].set_ylabel(old_hist.axes[1].member("fTitle") or "y")
    fig.colorbar(ratio_image, ax=axes[2], label="Ratio")
    fig.suptitle(key)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return True


def plot_efficiency_candidate(root_file, label: str, output: Path) -> bool:
    truth_key = "GenLevel/h2_Gen_DecayLeps_Eta_vs_pT"
    reco_key = "h2_Reco_DecayLeps_Eta_vs_pT"
    if truth_key not in root_file or reco_key not in root_file:
        return False
    truth, x_edges, y_edges = root_file[truth_key].to_numpy(flow=False)
    reco, reco_x, reco_y = root_file[reco_key].to_numpy(flow=False)
    if not np.array_equal(x_edges, reco_x) or not np.array_equal(y_edges, reco_y):
        raise ValueError("Truth and reconstructed eta-versus-pT axes differ")
    efficiency = np.divide(
        reco,
        truth,
        out=np.full_like(reco, np.nan, dtype=float),
        where=truth > 0,
    )
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    image = ax.pcolormesh(
        x_edges,
        y_edges,
        efficiency.T,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    fig.colorbar(image, ax=ax, label="Reco / truth")
    ax.set_xlabel(root_file[truth_key].axes[0].member("fTitle") or "pT")
    ax.set_ylabel(root_file[truth_key].axes[1].member("fTitle") or "eta")
    ax.set_title(f"Decay-lepton efficiency candidate: {label}")
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-root", required=True)
    parser.add_argument("--new-root", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--old-label", default="merged_elecs_test (1)")
    parser.add_argument("--new-label", default="merged_elecs_test_")
    parser.add_argument("--th1-mode", choices=("raw", "shape"), default="raw")
    parser.add_argument("--skip-th1", action="store_true")
    parser.add_argument("--skip-th2", action="store_true")
    args = parser.parse_args()

    old_path = Path(args.old_root).expanduser().resolve()
    new_path = Path(args.new_root).expanduser().resolve()
    output_dir = Path(args.out_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, object] = {
        "old_root": str(old_path),
        "new_root": str(new_path),
        "th1_mode": args.th1_mode,
        "written": [],
        "skipped_empty": [],
    }
    written: list[str] = manifest["written"]  # type: ignore[assignment]
    skipped: list[str] = manifest["skipped_empty"]  # type: ignore[assignment]

    with uproot.open(old_path) as old_file, uproot.open(new_path) as new_file:
        old_classes = old_file.classnames(cycle=False, recursive=True)
        new_classes = new_file.classnames(cycle=False, recursive=True)
        common = sorted(set(old_classes) & set(new_classes))
        if not args.skip_th1:
            for key in common:
                if not old_classes[key].startswith("TH1") or not new_classes[key].startswith("TH1"):
                    continue
                output = output_dir / f"{safe_name(key)}.png"
                if plot_th1(
                    old_file[key],
                    new_file[key],
                    key,
                    output,
                    args.old_label,
                    args.new_label,
                    args.th1_mode,
                ):
                    written.append(str(output))
                else:
                    skipped.append(key)
        if not args.skip_th2:
            for key in common:
                if not old_classes[key].startswith("TH2") or not new_classes[key].startswith("TH2"):
                    continue
                output = output_dir / f"{safe_name(key)}_TH2.png"
                if plot_th2(
                    old_file[key],
                    new_file[key],
                    key,
                    output,
                    args.old_label,
                    args.new_label,
                ):
                    written.append(str(output))
                else:
                    skipped.append(key)
            for root_file, suffix, label in (
                (old_file, "old", args.old_label),
                (new_file, "new", args.new_label),
            ):
                output = output_dir / f"DecayLepton_Efficiency_Eta_vs_pT_{suffix}.png"
                if plot_efficiency_candidate(root_file, label, output):
                    written.append(str(output))

    manifest["written_count"] = len(written)
    manifest["skipped_empty_count"] = len(skipped)
    manifest_path = output_dir / "campaign_histogram_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Rendered {len(written)} files to {output_dir}; skipped {len(skipped)} empty histograms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
