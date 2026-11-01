"""Reject missing or inclusive provenance in the separate tagged notebooks.

These checks authenticate agreement between saved framework metadata, not the
upstream producer or an experimental efficiency/coverage validation.
"""
from pathlib import Path
import hashlib
import json

EXPECTED_POLICY = {
    "mode": "required",
    "contract": "jul26-reco-topology-v1",
    "predicate": "reco_topology in {1,2} and finite reco_tagged_proton_p > 0",
    "applies_to": "reconstructed selection only",
    "truth_denominator": "unchanged; untagged generated events retained",
    "normalization": "all generated rows; unchanged cross section, luminosity and flux",
}


def require_policy(policy, source):
    if not isinstance(policy, dict):
        raise ValueError(f"{source}: missing proton-selection metadata. Use tagged numerical outputs; no inclusive fallback.")
    disagreements = [key for key, value in EXPECTED_POLICY.items() if policy.get(key) != value]
    if disagreements:
        raise ValueError(f"{source}: incompatible proton selection ({', '.join(disagreements)}). "
                         "This notebook requires the recorded reconstructed-proton tag and unchanged truth normalization.")
    return dict(policy)


def _read(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing tagged numerical input: {path}. Copy the complete tagged run, including JSON metadata.")
    return json.loads(path.read_text())


def full_run_policy(config, state, performance_audit):
    root = Path(config.runs[state])
    root = root.parent if root.name == "artifacts" else root
    receipt = root / "proton_selection.json"
    recorded = require_policy(_read(receipt).get("proton_selection"), receipt)
    base = root / "artifacts/performance" / config.energy
    for data, path in ((performance_audit, base / "performance_audit.json"),
                       (_read(base / "mass_fit.json"), base / "mass_fit.json")):
        policy = require_policy(data.get("proton_selection"), path)
        if policy != recorded:
            raise ValueError(f"{path}: proton-selection metadata disagrees with {receipt}.")
    return {"proton_selection": recorded, "proton_selection_mode": recorded["mode"],
            "proton_contract_source": str(receipt),
            "proton_contract_sha256": hashlib.sha256(receipt.read_bytes()).hexdigest()}


def reoptimization_policy(source_config, entry, config_path, result_path):
    recorded = require_policy(source_config.get("proton_selection"), config_path)
    policies = [(entry.get("proton_selection"), str(result_path)),
                ((entry.get("selection_metadata") or {}).get("proton_selection"),
                 str(result_path) + ":selection_metadata")]
    for raw, source in policies:
        policy = require_policy(raw, source)
        if policy != recorded:
            raise ValueError(f"{source}: proton-selection metadata disagrees with {config_path}.")
    return recorded
