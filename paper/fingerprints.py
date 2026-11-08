"""Canonical digests of the numerical quantities shown in the paper figures."""
import hashlib
import json
from pathlib import Path


REFERENCE = Path(__file__).with_name("reference") / "checksums.json"


def digest(rows):
    # Scientific values are formatted by each caller before hashing. Use JSON
    # arrays to preserve field boundaries and ASCII for platform-independent bytes.
    payload = json.dumps(sorted(rows), separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def verify(name, rows, reference=REFERENCE):
    actual = digest(rows)
    expected = json.loads(reference.read_text())[name]
    if actual != expected:
        raise ValueError(f"{name}: numerical fingerprint differs from the paper reference ({actual})")
    return actual


def map_rows(frame):
    return cell_rows(frame.to_dict("records"))


def cell_rows(records):
    return [[str(row["map"]), *[format(row[col], ".12g") for col in ("x_low","x_high","y_low","y_high")],
             *[int(row[col]) for col in ("denominator","selected_before_proton","selected")]]
            for row in records]


def mass_rows(result):
    return [[state, quantity, key, format(values[key], ".1f")]
            for state, states in result["states"].items()
            for quantity, values in states.items()
            for key in ("value", "error_low", "error_high")]


def spectrum_rows(rows):
    return [[r["quantity"], r["state"], *[format(float(r[k]), ".7g") for k in
             ("bin_low","bin_high","central","lower","upper")], r["status"]] for r in rows]
