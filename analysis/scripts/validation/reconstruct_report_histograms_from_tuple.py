#!/usr/bin/env python3
"""Reconstruct the report-era TH2 observables directly from ``UpsilonTree``.

The rules in this module were recovered by comparing the supplied
``merged_elecs_test_tuples.root`` and ``merged_elecs_test_.root`` files.  They
reproduce every regular and flow bin of all ten stored TH2D histograms with
unit weights.  No reconstruction macro is loaded or called.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Histogram2DSpec:
    key: str
    selection: str
    pairs: tuple[tuple[str, str], ...]
    x_edges: tuple[float, ...]
    y_edges: tuple[float, ...]
    report_use: str


def regular_edges(low: float, high: float, bins: int) -> tuple[float, ...]:
    return tuple(np.linspace(low, high, bins + 1, dtype=float))


SPECS: tuple[Histogram2DSpec, ...] = (
    Histogram2DSpec(
        "GenLevel/h2_Gen_DecayLeps_Eta_vs_Eta",
        "truth",
        (("true_decay_eMinus_eta", "true_ePlus_eta"),),
        regular_edges(-5.0, 5.0, 100),
        regular_edges(-5.0, 5.0, 100),
        "Original report p. 7 truth decay-lepton eta correlation",
    ),
    Histogram2DSpec(
        "GenLevel/h2_Gen_DecayLeps_pT_vs_pT",
        "truth",
        (("true_decay_eMinus_pT", "true_ePlus_pT"),),
        regular_edges(0.0, 12.0, 100),
        regular_edges(0.0, 12.0, 100),
        "Original report p. 7 decay-lepton pT correlation",
    ),
    Histogram2DSpec(
        "GenLevel/h2_Gen_DecayLeps_Eta_vs_pT",
        "truth",
        (
            ("true_decay_eMinus_pT", "true_decay_eMinus_eta"),
            ("true_ePlus_pT", "true_ePlus_eta"),
        ),
        regular_edges(0.0, 10.0, 100),
        regular_edges(-5.0, 5.0, 100),
        "Original report pp. 7 and 9 truth eta-versus-pT map",
    ),
    Histogram2DSpec(
        "GenLevel/h2_Gen_W_vs_Rapidity",
        "truth",
        (("true_y", "true_W"),),
        regular_edges(-5.0, 5.0, 100),
        regular_edges(0.0, 150.0, 100),
        "Original report pp. 6 and 12 W-versus-rapidity panels",
    ),
    Histogram2DSpec(
        "GenLevel/h2_Gen_t_vs_W",
        "truth",
        (("true_W", "true_t"),),
        regular_edges(10.0, 150.0, 100),
        regular_edges(0.0, 4.0, 100),
        "Original report p. 6 generated |t|-versus-W panel",
    ),
    Histogram2DSpec(
        "h2_nTracks_Central_vs_Tagger",
        "all",
        (("reco_nTracks_Central", "reco_nTracks_Tagger"),),
        regular_edges(0.0, 15.0, 15),
        regular_edges(0.0, 15.0, 15),
        "Follow-up report p. 4 raw central-versus-tagger track map",
    ),
    Histogram2DSpec(
        "h2_eta_ePlus_eMinus",
        "reco",
        (("reco_decay_eMinus_eta", "reco_ePlus_eta"),),
        regular_edges(-5.0, 5.0, 100),
        regular_edges(-5.0, 5.0, 100),
        "Original report p. 9 reconstructed decay-lepton eta correlation",
    ),
    Histogram2DSpec(
        "h2_pT_ePlus_eMinus",
        "reco",
        (("reco_decay_eMinus_pT", "reco_ePlus_pT"),),
        regular_edges(0.0, 12.0, 100),
        regular_edges(0.0, 12.0, 100),
        "Reconstructed decay-lepton pT correlation",
    ),
    Histogram2DSpec(
        "h2_Reco_DecayLeps_Eta_vs_pT",
        "reco",
        (
            ("reco_decay_eMinus_pT", "reco_decay_eMinus_eta"),
            ("reco_ePlus_pT", "reco_ePlus_eta"),
        ),
        regular_edges(0.0, 10.0, 100),
        regular_edges(-5.0, 5.0, 100),
        "Original report p. 9 reconstructed eta-versus-pT map",
    ),
    Histogram2DSpec(
        "h2_t_Correlation",
        "tagged",
        (("reco_t_eX", "reco_t_eXBABE"),),
        regular_edges(0.0, 2.0, 50),
        regular_edges(0.0, 2.0, 50),
        "Reconstructed |t| estimator correlation",
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def required_branches(specs: Iterable[Histogram2DSpec] = SPECS) -> list[str]:
    branches = {
        "true_mass",
        "reco_mass",
        "sigma_epz",
        "reco_t_BABE",
    }
    branches.update(proton_branches())
    for spec in specs:
        for x_branch, y_branch in spec.pairs:
            branches.update((x_branch, y_branch))
    return sorted(branches)


def selection_mask(arrays: dict[str, np.ndarray], name: str) -> np.ndarray:
    size = len(next(iter(arrays.values())))
    if name == "all":
        return np.ones(size, dtype=bool)
    if name == "truth":
        return arrays["true_mass"] > 0.0
    reco = (
        (arrays["reco_mass"] >= 8.0)
        & (arrays["reco_mass"] <= 12.0)
        & (np.abs(arrays["sigma_epz"] - 36.0) <= 3.6)
    )
    # Preserve the explicit inclusive switch for historical report replay.
    reco &= proton_tag_mask(arrays)
    if name == "reco":
        return reco
    if name == "tagged":
        if proton_selection_metadata()["mode"] == "inclusive":
            return reco & (arrays["reco_t_BABE"] > 0.0)
        return reco
    raise ValueError(f"Unknown selection {name!r}")


def histogram2d_with_flow(
    x: np.ndarray,
    y: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> np.ndarray:
    """Return ROOT-style bins ordered as underflow, regular bins, overflow."""

    x_bin = np.searchsorted(x_edges, x, side="right")
    y_bin = np.searchsorted(y_edges, y, side="right")
    x_bin = np.where(np.isnan(x), len(x_edges), x_bin)
    y_bin = np.where(np.isnan(y), len(y_edges), y_bin)
    counts = np.zeros((len(x_edges) + 1, len(y_edges) + 1), dtype=np.int64)
    np.add.at(counts, (x_bin, y_bin), 1)
    return counts


def reconstruct(spec: Histogram2DSpec, arrays: dict[str, np.ndarray]) -> np.ndarray:
    mask = selection_mask(arrays, spec.selection)
    x_edges = np.asarray(spec.x_edges)
    y_edges = np.asarray(spec.y_edges)
    result = np.zeros((len(x_edges) + 1, len(y_edges) + 1), dtype=np.int64)
    for x_branch, y_branch in spec.pairs:
        result += histogram2d_with_flow(
            arrays[x_branch][mask],
            arrays[y_branch][mask],
            x_edges,
            y_edges,
        )
    return result


def npz_name(key: str) -> str:
    return key.replace("/", "__")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tuple-root", required=True)
    parser.add_argument("--reference-root")
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--out-npz")
    args = parser.parse_args()

    try:
        import uproot
    except ImportError as exc:
        raise SystemExit("uproot is required") from exc

    tuple_path = Path(args.tuple_root).expanduser().resolve()
    if not tuple_path.is_file():
        raise SystemExit(f"Tuple ROOT file does not exist: {tuple_path}")

    with uproot.open(tuple_path) as root_file:
        if args.tree not in root_file:
            raise SystemExit(f"Tree {args.tree!r} is missing from {tuple_path}")
        tree = root_file[args.tree]
        available = set(tree.keys())
        missing = sorted(set(required_branches()) - available)
        if missing:
            raise SystemExit(f"Tuple is missing required branches: {', '.join(missing)}")
        arrays = tree.arrays(required_branches(), library="np")

    reference_path = Path(args.reference_root).expanduser().resolve() if args.reference_root else None
    reference_file = None
    if reference_path:
        if not reference_path.is_file():
            raise SystemExit(f"Reference ROOT file does not exist: {reference_path}")
        reference_file = uproot.open(reference_path)

    records: list[dict[str, object]] = []
    payload: dict[str, np.ndarray] = {}
    all_exact = True
    try:
        for spec in SPECS:
            observed = reconstruct(spec, arrays)
            key_name = npz_name(spec.key)
            payload[key_name] = observed[1:-1, 1:-1]
            payload[f"{key_name}__x_edges"] = np.asarray(spec.x_edges)
            payload[f"{key_name}__y_edges"] = np.asarray(spec.y_edges)

            record: dict[str, object] = {
                "key": spec.key,
                "selection": spec.selection,
                "pairs": [list(pair) for pair in spec.pairs],
                "x_bins": len(spec.x_edges) - 1,
                "x_range": [spec.x_edges[0], spec.x_edges[-1]],
                "y_bins": len(spec.y_edges) - 1,
                "y_range": [spec.y_edges[0], spec.y_edges[-1]],
                "fills": int(observed.sum()),
                "visible_fills": int(observed[1:-1, 1:-1].sum()),
                "flow_fills": int(observed.sum() - observed[1:-1, 1:-1].sum()),
                "report_use": spec.report_use,
                "unit_weight": True,
            }

            if reference_file is not None:
                if spec.key not in reference_file:
                    record.update({"exact": False, "error": "missing reference histogram"})
                    all_exact = False
                else:
                    reference_hist = reference_file[spec.key]
                    expected = np.asarray(reference_hist.values(flow=True))
                    expected_x = np.asarray(reference_hist.axes[0].edges())
                    expected_y = np.asarray(reference_hist.axes[1].edges())
                    edges_exact = np.array_equal(expected_x, spec.x_edges) and np.array_equal(
                        expected_y, spec.y_edges
                    )
                    if expected.shape != observed.shape:
                        exact = False
                        l1 = None
                        max_abs = None
                        differing = None
                    else:
                        delta = observed.astype(np.float64) - expected
                        exact = bool(edges_exact and np.array_equal(observed, expected))
                        l1 = float(np.abs(delta).sum())
                        max_abs = float(np.abs(delta).max(initial=0.0))
                        differing = int(np.count_nonzero(delta))
                    record.update(
                        {
                            "exact": exact,
                            "edges_exact": bool(edges_exact),
                            "differing_bins": differing,
                            "l1_difference": l1,
                            "max_abs_difference": max_abs,
                        }
                    )
                    all_exact &= exact
            records.append(record)
    finally:
        if reference_file is not None:
            reference_file.close()

    truth_key = npz_name("GenLevel/h2_Gen_DecayLeps_Eta_vs_pT")
    reco_key = npz_name("h2_Reco_DecayLeps_Eta_vs_pT")
    efficiency = np.divide(
        payload[reco_key],
        payload[truth_key],
        out=np.full_like(payload[reco_key], np.nan, dtype=float),
        where=payload[truth_key] > 0,
    )
    payload["derived__decay_lepton_efficiency_eta_pt"] = efficiency

    report = {
        "proton_selection": proton_selection_metadata(),
        "tuple": {
            "path": str(tuple_path),
            "size_bytes": tuple_path.stat().st_size,
            "sha256": sha256(tuple_path),
            "tree": args.tree,
            "entries": len(next(iter(arrays.values()))),
        },
        "reference": (
            {
                "path": str(reference_path),
                "size_bytes": reference_path.stat().st_size,
                "sha256": sha256(reference_path),
            }
            if reference_path
            else None
        ),
        "selection_counts": {
            name: int(selection_mask(arrays, name).sum())
            for name in ("all", "truth", "reco", "tagged")
        },
        "histograms": records,
        "all_exact": all_exact if reference_path else None,
        "derived_efficiency_note": (
            "The eta-versus-pT efficiency is reconstructed as reco/truth for occupied truth bins. "
            "The input maps and binning are exact; the precise report-era ratio convention was not retained."
        ),
    }

    out_json = Path(args.out_json).expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.out_npz:
        out_npz = Path(args.out_npz).expanduser().resolve()
        out_npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(out_npz, **payload)

    print(f"Reconstructed {len(records)} TH2 histograms from {tuple_path}")
    if reference_path:
        exact_count = sum(bool(record.get("exact")) for record in records)
        print(f"Reference comparison: {exact_count}/{len(records)} exact, including flow bins")
        if not all_exact:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
