"""Explicit reconstructed proton category for the audited July Upsilon producer.

Never apply this mask to generated rows, truth denominators or input weights.
The optional inclusive mode reproduces the previous no-proton-requirement category.
"""
import os

PROTON_BRANCHES = ("reco_topology", "reco_tagged_proton_p")
CONTRACT = "jul26-reco-topology-v1"


def proton_selection_metadata():
    mode = os.environ.get("YPSILON_PROTON_SELECTION", "required")
    if mode not in ("required", "inclusive"):
        raise ValueError("YPSILON_PROTON_SELECTION must be required or inclusive")
    return {
        "mode": mode,
        "contract": CONTRACT,
        "label": "Proton tagged" if mode == "required" else "No proton requirement",
        "predicate": "reco_topology in {1,2} and finite reco_tagged_proton_p > 0" if mode == "required" else "none",
        "applies_to": "reconstructed selection only",
        "truth_denominator": "unchanged; untagged generated events retained",
        "normalization": "all generated rows; unchanged cross section, luminosity and flux",
        "tag_definition": "producer RP/off-momentum candidate; not a truth-match or exclusivity guarantee",
    }


def proton_branches():
    return list(PROTON_BRANCHES) if proton_selection_metadata()["mode"] == "required" else []


def _size(arrays):
    import numpy as np
    if hasattr(arrays, "columns"):
        return len(arrays)
    if isinstance(arrays, np.ndarray):
        return len(arrays)
    return len(next(iter(arrays.values())))


def proton_tag_mask(arrays):
    import numpy as np
    if proton_selection_metadata()["mode"] == "inclusive":
        return np.ones(_size(arrays), dtype=bool)
    missing = [key for key in PROTON_BRANCHES if key not in arrays]
    if missing:
        raise ValueError("Proton-tagged selection requires branches: " + ", ".join(missing)
                         + ". No silent untagged fallback is permitted.")
    topology = np.asarray(arrays[PROTON_BRANCHES[0]], dtype=float)
    momentum = np.asarray(arrays[PROTON_BRANCHES[1]], dtype=float)
    if topology.ndim != 1 or momentum.shape != topology.shape:
        raise ValueError("Proton tag branches must be scalar, equal-length event arrays")
    valid_topology = np.isfinite(topology) & np.isin(topology, (0, 1, 2, 3, 4))
    tagged = np.isin(topology, (1, 2))
    coherent = valid_topology & np.where(tagged, np.isfinite(momentum) & (momentum > 0), momentum == -999.)
    if not np.all(coherent):
        where = np.flatnonzero(~coherent)
        raise ValueError("Proton producer/schema mismatch in %d rows; first indices %s. "
                         "Expected topology 1/2 with positive finite proton momentum, "
                         "or topology 0/3/4 with -999 sentinel. No repair is applied."
                         % (len(where), where[:10].tolist()))
    return tagged


def guard_proton_run(config, root, dry_run=False):
    """Prevent cached inclusive outputs from masquerading as tagged outputs."""
    import hashlib
    import json
    from pathlib import Path
    root = Path(root).resolve()
    run_dir = Path(config["run"]["dir"]).resolve()
    policy = proton_selection_metadata()
    config.setdefault("analysis", {})["proton_selection"] = policy
    files = [p for folder in ("workflow", "scripts", "tools")
             for p in (root / folder).rglob("*.py") if "__pycache__" not in p.parts]
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(files)}
    contract = {"proton_selection": policy, "source_sha256": hashes,
                "run_context": config.get("proton_run_context")}
    receipt = run_dir / "proton_selection.json"
    if run_dir.exists() and any(run_dir.iterdir()):
        if not receipt.exists():
            raise ValueError("Existing run has no proton-selection contract. Use a fresh run ID; old outputs are unchanged.")
        if json.loads(receipt.read_text()) != contract:
            raise ValueError("Run selection or analysis source differs from its recorded contract. Use a fresh run ID; old outputs are unchanged.")
    if not dry_run:
        run_dir.mkdir(parents=True, exist_ok=True)
        receipt.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    return contract
