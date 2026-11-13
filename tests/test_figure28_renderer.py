"""Saved-result integrity and reproducibility for the Figure 28 renderer."""

import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from tools import render_paper_figure28 as figure28


ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "paper" / "figure28_inputs"
MASS_AUDIT = figure28.MASS_PREFIX + "mass_reference/upsilon_es_reconstruction_audit.json"
POINTS = figure28.DIFFERENTIAL_PREFIX + "upsilon_es_differential_points.csv"
DIFFERENTIAL_AUDIT = figure28.DIFFERENTIAL_PREFIX + "upsilon_es_differential_audit.json"


def replace_verified_source(archive, relative, payload):
    """Give a modified test input a matching digest to exercise semantic checks."""
    (archive / relative).write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    manifest = archive / "SHA256SUMS.txt"
    lines = []
    for line in manifest.read_text().splitlines():
        checksum, name = line.split("  ", 1)
        lines.append(f"{digest if name == relative else checksum}  {name}")
    manifest.write_text("\n".join(lines) + "\n")


class Figure28Renderer(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def test_bundled_models_preserve_points_intervals_and_diagnostics(self):
        archive = figure28.VerifiedArchive(BUNDLE)
        self.assertIn(figure28.FRAMEWORK_COMMIT, archive.read("README.md").decode())
        mass = figure28.load_mass(archive)
        differential = figure28.load_differential(archive)
        original_intervals = json.loads((BUNDLE / figure28.MASS_PREFIX / "mass_fit_intervals.json").read_text())
        self.assertEqual(mass["intervals"], original_intervals)
        self.assertEqual(mass["bin_model_max_abs_difference"], {"1S": 0.0, "2S": 0.0, "3S": 0.0})
        self.assertEqual(differential["bin_model_max_abs_difference"], {"1S": 0.0, "2S": 0.0})

        with (BUNDLE / POINTS).open(newline="") as source:
            rows = list(csv.DictReader(source))
        for variable, states in differential["points"].items():
            for state, points in states.items():
                originals = [row for row in rows if row["displayed"] == "True"
                             and row["variable"] == variable and row["state"] == state]
                self.assertEqual(len(points), len(originals))
                for point, original in zip(points, originals):
                    for key in ("x_low", "x_high", "mean", "lower", "upper"):
                        self.assertEqual(point[key], float(original[key]))
                    self.assertEqual(point["status"], original["overall_status"])

        self.assertEqual(differential["slopes"]["1S"]["status"], "PASS_CONDITIONAL")
        self.assertEqual(differential["slopes"]["2S"]["status"], "NOT_VALIDATED")
        self.assertEqual(differential["points"]["pT2"]["3S"], [])
        provenance = json.loads((BUNDLE / "provenance.json").read_text())
        self.assertEqual(list(archive.checked.values()), provenance["copied_files"])

    def test_corrupted_input_is_rejected_before_writing_figures(self):
        archive = self.directory / "inputs"
        shutil.copytree(BUNDLE, archive)
        path = archive / figure28.MASS_PREFIX / "mass_reference/upsilon_es_reconstruction_mass.csv"
        path.write_bytes(path.read_bytes() + b"\n")
        output = self.directory / "output"
        with patch.object(figure28, "make_figures") as drawing:
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                figure28.render(archive, output)
        drawing.assert_not_called()
        self.assertFalse(output.exists())

    def test_inconsistent_metadata_is_rejected_before_writing_figures(self):
        calibration = "framework/full_state_runs/r275-1only-L2p5-1s/artifacts/performance/9x275/mass_fit.json"
        cases = (
            (MASS_AUDIT, ("energy",), "9x130", "9x275"),
            (DIFFERENTIAL_AUDIT, ("target_lumi_invfb",), 100.0, "2.5 fb"),
            (calibration, ("state",), "2S", "Invalid saved mass calibration"),
            (DIFFERENTIAL_AUDIT, ("slope_fits", "2S", "overall_status"),
             "PASS_CONDITIONAL", "Slope validation status disagrees"),
        )
        for index, (relative, keys, value, message) in enumerate(cases):
            with self.subTest(relative=relative, keys=keys):
                archive = self.directory / f"inputs-{index}"
                shutil.copytree(BUNDLE, archive)
                data = json.loads((archive / relative).read_text())
                parent = data
                for key in keys[:-1]:
                    parent = parent[key]
                parent[keys[-1]] = value
                replace_verified_source(archive, relative, json.dumps(data).encode())
                output = self.directory / f"output-{index}"
                with patch.object(figure28, "make_figures") as drawing:
                    with self.assertRaisesRegex(ValueError, message):
                        figure28.render(archive, output)
                drawing.assert_not_called()
                self.assertFalse(output.exists())

    def test_wrong_archive_provenance_is_rejected(self):
        archive = self.directory / "inputs"
        shutil.copytree(BUNDLE, archive)
        payload = (archive / "README.md").read_bytes().replace(figure28.FRAMEWORK_COMMIT.encode(), b"other-commit")
        replace_verified_source(archive, "README.md", payload)
        output = self.directory / "output"
        with self.assertRaisesRegex(ValueError, "analysis-commit provenance"):
            figure28.render(archive, output)
        self.assertFalse(output.exists())

    def test_existing_output_is_preserved(self):
        output = self.directory / "output"
        output.mkdir()
        sentinel = output / "preserve.txt"
        sentinel.write_text("keep")
        with self.assertRaisesRegex(FileExistsError, "new output directory"):
            figure28.render(BUNDLE, output)
        self.assertEqual(sentinel.read_text(), "keep")
        self.assertEqual(list(output.iterdir()), [sentinel])

    def test_cli_bundle_render_is_deterministic_from_another_directory(self):
        audits = []
        for name in ("first", "second"):
            output = self.directory / name
            command = [sys.executable, str(ROOT / "tools" / "render_paper_figure28.py"),
                       "--results-dir", str(BUNDLE), "--output-dir", str(output)]
            result = subprocess.run(command, cwd=self.directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            audit = json.loads((output / "figure28_audit.json").read_text())
            self.assertEqual(len(audit["verified_sources"]), 12)
            self.assertEqual(len(audit["outputs"]), 6)
            self.assertFalse(audit["refit_performed"])
            self.assertFalse(audit["luminosity_rescaling_applied"])
            for artifact in audit["outputs"]:
                actual = hashlib.sha256((output / artifact["file"]).read_bytes()).hexdigest()
                self.assertEqual(actual, artifact["sha256"])
            audits.append(audit)
        self.assertEqual(audits[0], audits[1])
