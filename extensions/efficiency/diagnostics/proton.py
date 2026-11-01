"""Explicit reconstructed-proton category for the audited Upsilon tuple schema."""
import numpy as np

PROTON_BRANCHES = ("reco_topology", "reco_tagged_proton_p")


def policy(settings):
    value = settings.get("proton_selection", "required")
    if value not in ("required", "inclusive"):
        raise ValueError("proton_selection must be required or inclusive")
    return value


def proton_mask(arrays, settings, size):
    """Validate every stored tag row; never repair codes or proton momenta.

    Audited producer topology: 0=no valid lepton combination, 1/2=tagged
    central/tagger electron categories, 3/4=the corresponding untagged cases.
    A valid lepton combination is required before proton finding is attempted.
    """
    if policy(settings) == "inclusive":
        return np.ones(size, dtype=bool)
    if isinstance(arrays, np.ndarray) and arrays.dtype.names:
        arrays = {name: arrays[name] for name in arrays.dtype.names}
    missing = [name for name in PROTON_BRANCHES if name not in arrays]
    if missing:
        raise ValueError("Proton selection requires branches: " + ", ".join(missing))
    values = {}
    for name in PROTON_BRANCHES:
        value = np.asarray(arrays[name])
        if value.ndim != 1 or len(value) != size or value.dtype.kind not in "iuf":
            raise ValueError("Invalid scalar proton branch: " + name)
        values[name] = value
    topology, momentum = (values[name] for name in PROTON_BRANCHES)
    if np.any(~np.isfinite(topology) | ~np.isin(topology, (0, 1, 2, 3, 4))):
        raise ValueError("Unknown, noninteger or nonfinite reco_topology; producer/schema audit required")
    tagged = np.isin(topology, (1, 2))
    inconsistent = ((tagged & (~np.isfinite(momentum) | (momentum <= 0)))
                    | (~tagged & (momentum != -999)))
    if np.any(inconsistent):
        raise ValueError("Inconsistent proton schema: topology1/2 requires finite p>0; topology0/3/4 requires p==-999")
    return tagged


def policy_metadata(settings):
    value = policy(settings)
    return {
        "policy": value,
        "category": "proton tagged" if value == "required" else "proton inclusive",
        "predicate": "reco_topology in {1,2} and finite reco_tagged_proton_p > 0" if value == "required" else "No reconstructed proton requirement",
        "schema_validation": "All rows checked for known topology and topology/momentum consistency" if value == "required" else "NOT_REQUESTED; proton branches are not required or read",
        "producer_convention": "Topology0/3/4 have reco_tagged_proton_p==-999; topology1/2 have a positive finite momentum",
        "truth_denominator": "Unchanged; reconstructed proton requirement applies only to the numerator",
        "yield_convention": "Before/after yields are raw unweighted selected MC row counts; no projected luminosity is applied",
    }
