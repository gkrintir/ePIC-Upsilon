#!/usr/bin/env python3
"""Write a separate exposure-corrected historical W-unfolding CSV.

The frozen replay remains unchanged. A test sample drawn with probability f
from tuples weighted to L has expected exposure f*L. Historical cross-section
columns divided the unscaled test yields by L; this sidecar multiplies those
four columns by 1/f. No response, yield, closure, or other field is changed.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path


CROSS_SECTION_COLUMNS = (
    "sigma_unfolded_pb", "sigma_unfolded_err_pb",
    "sigma_truth_test_pb", "sigma_truth_test_err_pb",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        fields, rows = reader.fieldnames, list(reader)
    if not fields or len(fields) != len(set(fields)) or not rows:
        raise ValueError("Input must have a unique header and at least one data row")
    missing = set(CROSS_SECTION_COLUMNS) - set(fields)
    if missing:
        raise ValueError(f"Missing historical cross-section columns: {sorted(missing)}")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Malformed CSV row")
    return fields, rows


def correction(input_csv, out_dir, recorded_lumi_invfb, test_fraction):
    input_csv, out_dir = Path(input_csv), Path(out_dir)
    if input_csv.stem.endswith("_split_luminosity_corrected"):
        raise ValueError("Input is already labeled split-luminosity corrected; use the historical replay CSV")
    lumi, fraction = Decimal(str(recorded_lumi_invfb)), Decimal(str(test_fraction))
    if not lumi.is_finite() or lumi <= 0 or not fraction.is_finite() or not 0 < fraction <= 1:
        raise ValueError("Recorded luminosity must be positive and 0 < test fraction <= 1")
    factor, exposure = 1 / fraction, lumi * fraction
    fields, original = read_csv(input_csv)
    original_hash = digest(input_csv)
    output_csv = out_dir / (input_csv.stem + "_split_luminosity_corrected.csv")
    if output_csv.resolve() == input_csv.resolve():
        raise ValueError("Corrected output must differ from the historical input")
    corrected = []
    scaled_cells, nonfinite_cells = 0, []
    for i, row in enumerate(original):
        changed = row.copy()
        for column in CROSS_SECTION_COLUMNS:
            text = row[column].strip()
            if not text:
                nonfinite_cells.append(dict(row=i, column=column, value=row[column]))
                continue
            try:
                value = Decimal(text)
            except InvalidOperation as error:
                raise ValueError(f"Invalid numeric value in row {i}, {column}: {text}") from error
            if not value.is_finite():
                nonfinite_cells.append(dict(row=i, column=column, value=row[column]))
                continue
            changed[column] = str(value * factor)
            scaled_cells += 1
        corrected.append(changed)
    out_dir.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(corrected)
    read_fields, written = read_csv(output_csv)
    unchanged_columns = [column for column in fields if column not in CROSS_SECTION_COLUMNS]
    preserved = read_fields == fields and len(written) == len(original) and all(
        old[column] == new[column] for old, new in zip(original, written) for column in unchanged_columns)
    scale_verified = all(
        (new[column] == old[column] if (not old[column].strip() or not Decimal(old[column]).is_finite())
         else Decimal(new[column]) == Decimal(old[column]) * factor)
        for old, new in zip(original, written) for column in CROSS_SECTION_COLUMNS)
    input_unchanged = digest(input_csv) == original_hash
    unsupported = []
    for i, row in enumerate(original):
        if "efficiency" in row:
            try:
                efficiency = Decimal(row["efficiency"])
                if not efficiency.is_finite() or efficiency <= 0:
                    unsupported.append(i)
            except InvalidOperation:
                unsupported.append(i)
    audit = dict(
        validator_sha256=digest(__file__),
        input_csv=str(input_csv.resolve()), input_sha256=original_hash,
        output_csv=str(output_csv.resolve()), output_sha256=digest(output_csv),
        recorded_luminosity_invfb=float(lumi), test_fraction=float(fraction),
        expected_test_exposure_invfb=float(exposure), cross_section_scale_factor=float(factor),
        corrected_columns=list(CROSS_SECTION_COLUMNS), unchanged_columns=unchanged_columns,
        rows=len(original), scaled_numeric_cells=scaled_cells, preserved_nonfinite_cells=nonfinite_cells,
        checks=dict(all_four_cross_section_columns_scaled=scale_verified,
                    all_other_column_values_preserved_exactly=preserved,
                    historical_input_bytes_unchanged=input_unchanged),
        numerical_checks_passed=bool(scale_verified and preserved and input_unchanged),
        unsupported_response_rows=unsupported,
        interpretation="Separate correction of the historical missing test-subsampling exposure factor; frozen replay is retained as reproduction evidence, not certified absolute cross sections.",
        exposure_convention="Expected exposure uses the configured Bernoulli test probability, not the observed random test event fraction.",
        unchanged_physics="Response, efficiencies, purity, stability, raw weighted train/test and unfolded yields, and dimensionless closure values/errors are unchanged.",
        limitations=[
            "Recorded luminosity and test fraction are supplied run settings; this sidecar checks their conversion arithmetic, not independent provenance of those settings.",
            "Does not validate the assigned generator rate, photon flux, event selection, response model, or historical uncertainty prescription.",
            "Zero-efficiency bins remain unidentifiable. Retained zero placeholders are not measured zero cross sections and should not be interpreted as physical results.",
            "This correction does not apply to the separate full-sample bin-by-bin sigma_gp_family_*_fixed_unfolded.csv products.",
        ],
    )
    audit_path = out_dir / "split_luminosity_correction_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    if not audit["numerical_checks_passed"]:
        raise RuntimeError("Written normalization sidecar failed preservation/scale checks")
    return output_csv, audit_path, audit


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--recorded-lumi-invfb", type=Decimal, default=Decimal("10"))
    parser.add_argument("--test-fraction", type=Decimal, default=Decimal("0.5"))
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        output, audit_path, audit = correction(args.input, args.out_dir, args.recorded_lumi_invfb, args.test_fraction)
    except (ValueError, InvalidOperation) as error:
        parser.error(str(error))
    print(f"PASS: {audit['rows']} rows; four cross-section columns multiplied by {audit['cross_section_scale_factor']:g}; every other column preserved")
    print(f"Expected test exposure: {audit['expected_test_exposure_invfb']:g} fb^-1")
    print(f"Corrected CSV: {output.resolve()}")
    print(f"Audit: {audit_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
