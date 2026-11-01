#!/usr/bin/env python3
"""Validate that a ROOT file contains the branches required by an analysis profile."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


PROFILES = {
    "core": {"weight", "true_mass", "true_W", "reco_mass", "reco_W", "sigma_epz"},
    "w": {"weight", "true_mass", "true_W", "true_Q2", "reco_mass", "reco_W", "sigma_epz"},
    "t": {
        "weight", "true_mass", "true_W", "true_Q2", "true_t", "reco_mass", "reco_W", "sigma_epz",
        "reco_t_eX", "reco_t_BABE", "reco_t_eXBABE", "reco_t_proj", "reco_pT",
    },
    "differential": {
        "weight", "true_mass", "true_W", "true_Q2", "true_pT", "true_y", "true_t",
        "reco_mass", "reco_W", "reco_pT", "reco_y", "sigma_epz",
    },
    "correlation": {
        "weight", "true_mass", "true_W", "reco_mass", "reco_W", "sigma_epz",
        "reco_nTracks_Central", "reco_nTracks_Tagger",
    },
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("--tree", default="UpsilonTree")
    parser.add_argument("--profile", choices=sorted(PROFILES), required=True)
    parser.add_argument("--branch", action="append", default=[], help="Additional required branch")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    try:
        import uproot
    except ImportError as exc:
        raise SystemExit("uproot is required for ROOT schema validation") from exc

    path = Path(os.path.expandvars(os.path.expanduser(args.input))).resolve()
    if not path.is_file():
        raise SystemExit(f"ROOT input does not exist: {path}")
    with uproot.open(path) as root_file:
        if args.tree not in root_file:
            raise SystemExit(f"Tree {args.tree!r} is missing from {path}")
        tree = root_file[args.tree]
        available = sorted(tree.keys())
        entries = int(tree.num_entries)
    required = sorted(PROFILES[args.profile] | set(args.branch))
    missing = sorted(set(required) - set(available))
    record = {
        "input": str(path),
        "tree": args.tree,
        "entries": entries,
        "profile": args.profile,
        "required_branches": required,
        "available_branches": available,
        "missing_branches": missing,
        "valid": not missing,
    }
    out = Path(args.out).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if missing:
        raise SystemExit(f"Missing branches for profile {args.profile}: {', '.join(missing)}")
    print(f"Validated {path}: {len(available)} branches, {record['entries']} entries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
