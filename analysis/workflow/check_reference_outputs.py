#!/usr/bin/env python3
"""Audit committed reference CSVs for known scientific-output failures.

The checker intentionally uses only the Python standard library.  It can run in
CI, and its unit tests do not need ROOT, pandas, or any analysis input files.

Issue identifiers are part of the public interface.  They are composed from a
stable rule name, the path relative to the reference directory, and (where
needed) a state discriminator.  Row numbers and absolute paths are deliberately
excluded so that harmless CSV reordering does not invalidate the baseline.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterable, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE_DIR = PROJECT_ROOT / "results" / "reference"
DEFAULT_BASELINE = PROJECT_ROOT / "config" / "reference_quality_baseline.json"

# These are the visible dielectron cross sections used throughout the reference
# production scripts.  A target for the wrong beam energy/state is a scientific
# error; a raw or inferred cross section that is subsequently rescaled is not.
EXPECTED_VISIBLE_XSEC_FB: dict[str, dict[str, float]] = {
    "10x100": {"1S": 104.793, "2S": 29.297, "3S": 23.447},
    "18x275": {"1S": 327.584, "2S": 84.780, "3S": 68.929},
}

W_LABEL_RE = re.compile(
    r"(?:^|_)W(?P<low>\d+(?:p\d+)?)to(?P<high>\d+(?:p\d+)?)(?:_|\.)"
)
ENERGY_RE = re.compile(r"(?:^|_)(?P<energy>10x100|18x275)(?:_|\.)")
UNFOLD_METHOD_RE = re.compile(
    r"^unfold_t_[^_]+_(?:10x100|18x275)_"
    r"(?P<method>eX|pT2|BABE|eXBABE|proj)(?:_|\.)",
    re.IGNORECASE,
)

FLOAT_ABS_TOLERANCE = 1.0e-6
TARGET_ABS_TOLERANCE_FB = 1.0e-3
TRUE_VALUES = {"1", "true", "t", "yes", "y"}


class BaselineError(ValueError):
    """Raised for an absent or malformed issue baseline."""


@dataclass(frozen=True)
class Issue:
    """One deterministic reference-quality finding."""

    issue_id: str
    rule: str
    severity: str
    path: str
    message: str
    evidence: dict[str, Any]


@dataclass(frozen=True)
class AuditResult:
    """Detected issues split by their baseline status."""

    issues: tuple[Issue, ...]
    known_issues: tuple[Issue, ...]
    new_issues: tuple[Issue, ...]
    resolved_baseline_ids: tuple[str, ...]

    @property
    def new_errors(self) -> tuple[Issue, ...]:
        return tuple(issue for issue in self.new_issues if issue.severity == "error")


def _float(value: object) -> float | None:
    """Return a finite float, or ``None`` for blanks/NaN/invalid values."""

    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _bool(value: object) -> bool:
    return str(value).strip().lower() in TRUE_VALUES


def _label_number(value: str) -> float:
    return float(value.replace("p", "."))


def _relative(path: Path, reference_dir: Path) -> str:
    try:
        return path.relative_to(reference_dir).as_posix()
    except ValueError:
        return path.name


def _issue_id(rule: str, relative_path: str, discriminator: str | None = None) -> str:
    parts = [rule, relative_path]
    if discriminator:
        parts.append(discriminator)
    return "::".join(parts)


def _issue(
    rule: str,
    relative_path: str,
    message: str,
    evidence: dict[str, Any],
    *,
    discriminator: str | None = None,
) -> Issue:
    return Issue(
        issue_id=_issue_id(rule, relative_path, discriminator),
        rule=rule,
        severity="error",
        path=relative_path,
        message=message,
        evidence=evidence,
    )


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def _numeric_values(rows: Iterable[dict[str, str]], column: str) -> list[float]:
    return [number for row in rows if (number := _float(row.get(column))) is not None]


def _actual_w_range(
    fieldnames: set[str], rows: list[dict[str, str]]
) -> tuple[float, float, str, str] | None:
    # Prefer explicit analysis selection metadata over per-bin edges when both
    # are present.  The latter still works for binned W output tables.
    for low_column, high_column in (
        ("w_min", "w_max"),
        ("common_W_low", "common_W_high"),
        ("W_low", "W_high"),
    ):
        if low_column not in fieldnames or high_column not in fieldnames:
            continue
        lows = _numeric_values(rows, low_column)
        highs = _numeric_values(rows, high_column)
        if lows and highs:
            return min(lows), max(highs), low_column, high_column
    return None


def _w_label_issue(
    relative_path: str, fieldnames: set[str], rows: list[dict[str, str]]
) -> Issue | None:
    match = W_LABEL_RE.search(Path(relative_path).name)
    if not match:
        return None
    actual = _actual_w_range(fieldnames, rows)
    if actual is None:
        return None
    expected_low = _label_number(match.group("low"))
    expected_high = _label_number(match.group("high"))
    actual_low, actual_high, low_column, high_column = actual
    if math.isclose(actual_low, expected_low, abs_tol=FLOAT_ABS_TOLERANCE) and math.isclose(
        actual_high, expected_high, abs_tol=FLOAT_ABS_TOLERANCE
    ):
        return None
    return _issue(
        "W_RANGE_MISMATCH",
        relative_path,
        (
            f"Filename advertises W={expected_low:g}-{expected_high:g}, but CSV "
            f"content spans W={actual_low:g}-{actual_high:g}."
        ),
        {
            "filename_range": [expected_low, expected_high],
            "content_range": [actual_low, actual_high],
            "range_columns": [low_column, high_column],
        },
    )


def _has_positive_truth(fieldnames: set[str], rows: list[dict[str, str]]) -> bool | None:
    truth_columns = [
        column
        for column in (
            "n_true_weighted_visible",
            "n_truth_visible",
            "n_true",
            "sum_truth_visible",
        )
        if column in fieldnames
    ]
    if not truth_columns:
        return None
    return any(
        (number := _float(row.get(column))) is not None and number > FLOAT_ABS_TOLERANCE
        for column in truth_columns
        for row in rows
    )


def _zero_efficiency_and_plot_issues(
    relative_path: str, fieldnames: set[str], rows: list[dict[str, str]]
) -> list[Issue]:
    if not rows:
        return []

    positive_truth = _has_positive_truth(fieldnames, rows)
    # A state with explicitly zero truth content is diagnosed by the metadata
    # rule below, rather than being counted again as a detector-efficiency fault.
    if positive_truth is False:
        return []

    issues: list[Issue] = []
    efficiency_column = next(
        (name for name in ("eff_reco_over_true", "efficiency") if name in fieldnames),
        None,
    )
    if efficiency_column:
        efficiencies = _numeric_values(rows, efficiency_column)
        if efficiencies and all(abs(value) <= FLOAT_ABS_TOLERANCE for value in efficiencies):
            issues.append(
                _issue(
                    "ZERO_EFFICIENCY",
                    relative_path,
                    f"All {len(efficiencies)} finite efficiency values are zero.",
                    {
                        "efficiency_column": efficiency_column,
                        "finite_values": len(efficiencies),
                        "positive_truth_content": positive_truth,
                    },
                )
            )

    plotted_column = next(
        (name for name in ("plotted", "passes_min_eff") if name in fieldnames),
        None,
    )
    if plotted_column and not any(_bool(row.get(plotted_column)) for row in rows):
        issues.append(
            _issue(
                "NO_PLOTTED_BINS",
                relative_path,
                f"No row is selected by the {plotted_column!r} flag.",
                {
                    "selection_column": plotted_column,
                    "rows": len(rows),
                    "positive_truth_content": positive_truth,
                },
            )
        )
    return issues


def _energy_from_filename(relative_path: str) -> str | None:
    match = ENERGY_RE.search(Path(relative_path).name)
    return match.group("energy") if match else None


def _target_normalization_issues(
    relative_path: str, fieldnames: set[str], rows: list[dict[str, str]]
) -> list[Issue]:
    if "target_visible_xsec_fb" not in fieldnames or "state" not in fieldnames:
        return []
    energy = _energy_from_filename(relative_path)
    if energy is None:
        return []

    mismatches: dict[str, tuple[float, float]] = {}
    expected_by_state = EXPECTED_VISIBLE_XSEC_FB[energy]
    for row in rows:
        state = str(row.get("state", "")).strip()
        expected = expected_by_state.get(state)
        observed = _float(row.get("target_visible_xsec_fb"))
        if expected is None or observed is None:
            continue
        if not math.isclose(
            observed,
            expected,
            rel_tol=1.0e-9,
            abs_tol=TARGET_ABS_TOLERANCE_FB,
        ):
            mismatches.setdefault(state, (observed, expected))

    return [
        _issue(
            "TARGET_NORMALIZATION_MISMATCH",
            relative_path,
            (
                f"{state} target is {observed:g} fb, but {energy} reference "
                f"normalization is {expected:g} fb."
            ),
            {
                "energy": energy,
                "state": state,
                "observed_target_visible_xsec_fb": observed,
                "expected_target_visible_xsec_fb": expected,
            },
            discriminator=state,
        )
        for state, (observed, expected) in sorted(mismatches.items())
    ]


def _empty_state_issues(
    relative_path: str, fieldnames: set[str], rows: list[dict[str, str]]
) -> list[Issue]:
    # Metadata is the authoritative one-row-per-state summary.  Restricting this
    # rule to metadata avoids six redundant findings for each empty pT2/y table.
    if "metadata" not in Path(relative_path).stem.lower():
        return []
    required = {"state", "raw_visible_xsec_fb", "target_visible_xsec_fb"}
    if not required.issubset(fieldnames):
        return []

    state_rows: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        state = str(row.get("state", "")).strip()
        if state:
            state_rows.setdefault(state, []).append(row)

    issues: list[Issue] = []
    for state, grouped_rows in sorted(state_rows.items()):
        raw_values = _numeric_values(grouped_rows, "raw_visible_xsec_fb")
        target_values = _numeric_values(grouped_rows, "target_visible_xsec_fb")
        if not raw_values or not target_values:
            continue
        if all(abs(value) <= FLOAT_ABS_TOLERANCE for value in raw_values) and any(
            value > FLOAT_ABS_TOLERANCE for value in target_values
        ):
            issues.append(
                _issue(
                    "EMPTY_STATE_OUTPUT",
                    relative_path,
                    (
                        f"{state} has zero raw visible cross section despite a "
                        f"nonzero target ({max(target_values):g} fb)."
                    ),
                    {
                        "state": state,
                        "raw_visible_xsec_fb": raw_values,
                        "target_visible_xsec_fb": target_values,
                    },
                    discriminator=state,
                )
            )
    return issues


def _normalized_method(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _method_filename_issue(
    relative_path: str, fieldnames: set[str], rows: list[dict[str, str]]
) -> Issue | None:
    if "method" not in fieldnames:
        return None
    match = UNFOLD_METHOD_RE.search(Path(relative_path).name)
    if not match:
        return None
    filename_method = match.group("method")
    content_methods = sorted(
        {str(row.get("method", "")).strip() for row in rows if str(row.get("method", "")).strip()}
    )
    if not content_methods or all(
        _normalized_method(value) == _normalized_method(filename_method)
        for value in content_methods
    ):
        return None
    return _issue(
        "METHOD_FILENAME_MISMATCH",
        relative_path,
        (
            f"Filename method is {filename_method!r}, but CSV method value(s) are "
            f"{content_methods!r}."
        ),
        {"filename_method": filename_method, "content_methods": content_methods},
    )


def _row_is_low_energy(row: dict[str, str]) -> bool:
    evidence = " ".join(
        str(row.get(column, ""))
        for column in ("energy_label", "energy", "sample", "sample_label", "source_csv")
    ).lower()
    compact = re.sub(r"[\s$\\{}_-]+", "", evidence)
    return "low" in evidence or "10x100" in compact or "10times100" in compact


def _low_energy_overlay_issue(
    relative_path: str, fieldnames: set[str], rows: list[dict[str, str]]
) -> Issue | None:
    filename = Path(relative_path).name.lower()
    if "overlay" not in filename or "10x100" not in filename or "18x275" not in filename:
        return None
    evidence_columns = {
        "energy_label",
        "energy",
        "sample",
        "sample_label",
        "source_csv",
    }
    if not fieldnames.intersection(evidence_columns):
        return None
    low_rows = sum(_row_is_low_energy(row) for row in rows)
    if low_rows:
        return None
    return _issue(
        "LOW_ENERGY_OVERLAY_EMPTY",
        relative_path,
        "Two-energy overlay contains no identifiable 10x100/low-energy rows.",
        {
            "rows": len(rows),
            "low_energy_rows": low_rows,
            "evidence_columns": sorted(fieldnames.intersection(evidence_columns)),
        },
    )


def scan_reference_outputs(reference_dir: Path | str) -> list[Issue]:
    """Scan all reference CSVs and return issues sorted by stable identifier."""

    reference_dir = Path(reference_dir)
    if not reference_dir.is_dir():
        raise FileNotFoundError(f"Reference directory does not exist: {reference_dir}")

    issues_by_id: dict[str, Issue] = {}
    for path in sorted(reference_dir.rglob("*.csv")):
        relative_path = _relative(path, reference_dir)
        try:
            fieldname_list, rows = _read_csv(path)
        except (OSError, csv.Error, UnicodeError) as exc:
            issue = _issue(
                "CSV_READ_ERROR",
                relative_path,
                f"CSV could not be read: {exc}",
                {"exception_type": type(exc).__name__},
            )
            issues_by_id[issue.issue_id] = issue
            continue

        fieldnames = set(fieldname_list)
        if not fieldnames:
            issue = _issue(
                "CSV_EMPTY",
                relative_path,
                "CSV has no header.",
                {"size_bytes": path.stat().st_size},
            )
            issues_by_id[issue.issue_id] = issue
            continue

        candidates: list[Issue | None] = [
            _w_label_issue(relative_path, fieldnames, rows),
            _method_filename_issue(relative_path, fieldnames, rows),
            _low_energy_overlay_issue(relative_path, fieldnames, rows),
        ]
        candidates.extend(_zero_efficiency_and_plot_issues(relative_path, fieldnames, rows))
        candidates.extend(_target_normalization_issues(relative_path, fieldnames, rows))
        candidates.extend(_empty_state_issues(relative_path, fieldnames, rows))
        for issue in candidates:
            if issue is not None:
                issues_by_id[issue.issue_id] = issue

    return [issues_by_id[issue_id] for issue_id in sorted(issues_by_id)]


def load_baseline(path: Path | str) -> set[str]:
    """Load and validate a reference-quality baseline JSON file."""

    path = Path(path)
    if not path.is_file():
        raise BaselineError(f"Baseline does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BaselineError(f"Could not read baseline {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise BaselineError("Baseline root must be a JSON object")
    if payload.get("schema_version") != 1:
        raise BaselineError("Baseline schema_version must be 1")
    issue_ids = payload.get("known_issue_ids")
    if not isinstance(issue_ids, list) or not all(isinstance(item, str) for item in issue_ids):
        raise BaselineError("Baseline known_issue_ids must be an array of strings")
    if len(issue_ids) != len(set(issue_ids)):
        raise BaselineError("Baseline known_issue_ids contains duplicates")
    return set(issue_ids)


def classify_issues(issues: Sequence[Issue], baseline_ids: Iterable[str]) -> AuditResult:
    baseline = set(baseline_ids)
    ordered = tuple(sorted(issues, key=lambda issue: issue.issue_id))
    known = tuple(issue for issue in ordered if issue.issue_id in baseline)
    new = tuple(issue for issue in ordered if issue.issue_id not in baseline)
    detected_ids = {issue.issue_id for issue in ordered}
    resolved = tuple(sorted(baseline - detected_ids))
    return AuditResult(ordered, known, new, resolved)


def audit_reference_outputs(
    reference_dir: Path | str, baseline_path: Path | str
) -> AuditResult:
    """Scan a directory and classify its findings against a JSON baseline."""

    return classify_issues(scan_reference_outputs(reference_dir), load_baseline(baseline_path))


def result_payload(result: AuditResult) -> dict[str, Any]:
    status_by_id = {
        issue.issue_id: "known" for issue in result.known_issues
    } | {issue.issue_id: "new" for issue in result.new_issues}
    return {
        "schema_version": 1,
        "summary": {
            "detected": len(result.issues),
            "known": len(result.known_issues),
            "new": len(result.new_issues),
            "new_errors": len(result.new_errors),
            "resolved_baseline": len(result.resolved_baseline_ids),
        },
        "issues": [
            {**asdict(issue), "baseline_status": status_by_id[issue.issue_id]}
            for issue in result.issues
        ],
        "resolved_baseline_ids": list(result.resolved_baseline_ids),
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _print_text(result: AuditResult) -> None:
    print(
        "Reference quality: "
        f"{len(result.issues)} detected, {len(result.known_issues)} known, "
        f"{len(result.new_issues)} new ({len(result.new_errors)} new errors)."
    )
    known_ids = {issue.issue_id for issue in result.known_issues}
    for issue in result.issues:
        status = "KNOWN" if issue.issue_id in known_ids else "NEW"
        print(f"[{status}] {issue.issue_id}: {issue.message}")
    if result.resolved_baseline_ids:
        print(f"Resolved/stale baseline entries: {len(result.resolved_baseline_ids)}")
        for issue_id in result.resolved_baseline_ids:
            print(f"[RESOLVED] {issue_id}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help=f"Directory of CSV outputs (default: {DEFAULT_REFERENCE_DIR})",
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=DEFAULT_BASELINE,
        help=f"Known-issue baseline JSON (default: {DEFAULT_BASELINE})",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        help="Optionally write the complete machine-readable audit report here.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the text report (errors still use stderr).",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = audit_reference_outputs(args.reference_dir, args.baseline)
    except (BaselineError, FileNotFoundError) as exc:
        print(f"reference-quality configuration error: {exc}", file=sys.stderr)
        return 2
    if args.json_output:
        _write_json(args.json_output, result_payload(result))
    if not args.quiet:
        _print_text(result)
    return 1 if result.new_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
