"""Small identity and transfer regressions for the paper-input workflow."""
import csv
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from inputs import fetch_paper_inputs as inputs

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "plotting"))
from upsilon_esr.electron_id import electron_id_metadata


class PaperInputs(unittest.TestCase):
    def test_verified_id_is_bound_to_each_raw_input(self):
        with inputs.MANIFEST.open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(len(rows), 6)
        for row in rows:
            state_dir, energy = Path(row["relative_path"]).parts[:2]
            state = state_dir.removeprefix("Upsilon")
            saved = "new electron ID" if state == "1S" else "standard electron ID"
            with self.subTest(state=state, energy=energy):
                result = electron_id_metadata(saved, row["sha256"], state, energy)
                self.assertEqual(result["ele_id_mode"], 1)
                self.assertEqual(result["femc_mode"], 0)
                self.assertEqual(result["saved_category"], saved)
                self.assertEqual(result["category"], "reconstructed track–cluster matching")
        self.assertNotIn("ele_id_mode", electron_id_metadata("control", "unknown", "1S", "9x275"))
        with self.assertRaises(ValueError):
            electron_id_metadata("control", rows[0]["sha256"], "3S", "9x275")

    def test_failed_hash_never_installs_partial_download(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(inputs, "read_chunk", return_value=b"abcd"):
            target = Path(folder) / "sample.root"
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                inputs.fetch("sample.root", target, size=4, checksum=hashlib.sha256(b"wxyz").hexdigest())
            self.assertEqual(list(Path(folder).iterdir()), [])

    def test_completed_download_preserves_conflicting_destination(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(inputs, "read_chunk", return_value=b"abcd"):
            target = Path(folder) / "sample.root"
            target.write_bytes(b"keep")
            with self.assertRaises(ValueError):
                inputs.fetch("sample.root", target, size=4, checksum=hashlib.sha256(b"abcd").hexdigest())
            self.assertEqual(target.read_bytes(), b"keep")
            self.assertEqual(list(Path(folder).iterdir()), [target])
