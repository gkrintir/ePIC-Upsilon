import json
from pathlib import Path
import tempfile
import unittest
import warnings
from upsilon_esr.common import Config, sample_metadata


class MetadataTests(unittest.TestCase):
    def build(self, directory, status):
        root = Path(directory)
        config = Config(runs={s: root/s for s in ("1S", "2S", "3S")})
        base = root / "1S/artifacts"
        perf = base / "performance/9x275/performance_audit.json"
        perf.parent.mkdir(parents=True)
        perf.write_text(json.dumps({"state": "1S", "energy_label": "9x275 GeV",
            "truth_selection": {"Q2_GeV2": [0, 1], "W_GeV": [20, 100]},
            "statistical_configuration": {"active_luminosity_invfb": 10}}))
        audit = base / "normalization_audit/9x275/audit.json"
        audit.parent.mkdir(parents=True)
        audit.write_text(json.dumps({"status": status, "settings": {"q2_cut": 1, "q2_max": 1},
            "checks": [{"id": "normalization_test", "status": status, "detail": "fixture"}]}))
        return config

    def test_failed_audit_is_not_published(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.build(directory, "FAIL")
            with self.assertRaisesRegex(ValueError, "normalization_test"):
                sample_metadata(config, "1S")

    def test_warning_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.build(directory, "WARN")
            with warnings.catch_warnings(record=True) as messages:
                result = sample_metadata(config, "1S")
            self.assertTrue(messages)
            self.assertEqual(result["normalization_audit_status"], "WARN")
            self.assertEqual(result["normalization_audit_warnings"][0]["id"], "normalization_test")


if __name__ == "__main__":
    unittest.main()
