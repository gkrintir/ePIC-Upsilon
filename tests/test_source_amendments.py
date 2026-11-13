"""Keep imported source hashes auditable when a reviewed script is amended."""

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('check_repository', REPO / 'tools/check_repository.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class SourceAmendmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'analysis').mkdir()
        self.path = 'analysis/renderer.py'
        original = b'COUNTER = 1\n'
        self.current = b'COUNTER = 2\n'
        for name in ('renderer.py', 'preserved.py'):
            (self.root / 'analysis' / name).write_bytes(original)
        baseline = hashlib.sha256(original).hexdigest()
        self.provenance = {'copied_files': [
            {'path': self.path, 'sha256': baseline},
            {'path': 'analysis/preserved.py', 'sha256': baseline},
        ]}
        self.record = {
            'path': self.path,
            'baseline_sha256': baseline,
            'new_sha256': hashlib.sha256(self.current).hexdigest(),
            'reason': 'Correct the producer-counter chain labels.',
        }

    def amend(self):
        (self.root / self.path).write_bytes(self.current)
        self.provenance['amended_files'] = [self.record]

    def test_original_snapshot_without_amendments(self):
        counts, errors = checker.check_analysis_snapshot(self.root, self.provenance)
        self.assertEqual(counts, {'preserved_analysis_files': 2, 'amended_analysis_files': 0})
        self.assertEqual(errors, [])

    def test_recorded_amendment_and_preserved_file(self):
        self.amend()
        counts, errors = checker.check_analysis_snapshot(self.root, self.provenance)
        self.assertEqual(counts, {'preserved_analysis_files': 1, 'amended_analysis_files': 1})
        self.assertEqual(errors, [])

    def test_undeclared_change_is_rejected(self):
        self.amend()
        (self.root / 'analysis/preserved.py').write_bytes(self.current)
        _, errors = checker.check_analysis_snapshot(self.root, self.provenance)
        self.assertEqual(errors, ['analysis/preserved.py: differs from preserved analysis baseline'])

    def test_changed_or_missing_amended_file_is_rejected(self):
        self.amend()
        source = self.root / self.path
        source.write_bytes(b'COUNTER = 3\n')
        for missing in (False, True):
            with self.subTest(missing=missing):
                if missing:
                    source.unlink()
                _, errors = checker.check_analysis_snapshot(self.root, self.provenance)
                self.assertEqual(errors, [self.path + ': differs from recorded analysis amendment'])

    def test_baseline_hash_must_match_original_record(self):
        self.amend()
        self.record['baseline_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'baseline differs from copied_files'):
            checker.check_analysis_snapshot(self.root, self.provenance)

    def test_unknown_path_is_rejected(self):
        self.amend()
        self.record['path'] = 'analysis/unknown.py'
        with self.assertRaisesRegex(ValueError, 'Unknown amended path'):
            checker.check_analysis_snapshot(self.root, self.provenance)

    def test_duplicate_amendment_is_rejected(self):
        self.amend()
        self.provenance['amended_files'].append(dict(self.record))
        with self.assertRaisesRegex(ValueError, 'duplicate amendment'):
            checker.check_analysis_snapshot(self.root, self.provenance)

    def test_duplicate_baseline_is_rejected(self):
        self.provenance['copied_files'].append(dict(self.provenance['copied_files'][0]))
        with self.assertRaisesRegex(ValueError, 'copied_files contains duplicate paths'):
            checker.check_analysis_snapshot(self.root, self.provenance)

    def test_missing_or_unknown_fields_are_rejected(self):
        self.amend()
        invalid_records = [
            {key: value for key, value in self.record.items() if key != 'reason'},
            {**self.record, 'unknown': True},
        ]
        for record in invalid_records:
            with self.subTest(record=record):
                self.provenance['amended_files'] = [record]
                with self.assertRaisesRegex(ValueError, 'Each amendment requires exactly'):
                    checker.check_analysis_snapshot(self.root, self.provenance)

    def test_reason_must_be_nonempty_text(self):
        self.amend()
        for reason in ('', '  \n', None, 1):
            with self.subTest(reason=reason):
                self.record['reason'] = reason
                with self.assertRaisesRegex(ValueError, 'nonempty reason'):
                    checker.check_analysis_snapshot(self.root, self.provenance)

    def test_hashes_must_be_sha256_hex(self):
        self.amend()
        for field in ('baseline_sha256', 'new_sha256'):
            for digest in ('f' * 63, 'g' * 64, None):
                with self.subTest(field=field, digest=digest):
                    record = {**self.record, field: digest}
                    self.provenance['amended_files'] = [record]
                    with self.assertRaisesRegex(ValueError, 'invalid ' + field):
                        checker.check_analysis_snapshot(self.root, self.provenance)

    def test_main_reports_separate_counts_and_rejects_invalid_provenance(self):
        self.amend()
        (self.root / 'README.md').write_text('# Test repository\n')
        (self.root / 'provenance').mkdir()
        manifest = self.root / 'provenance/analysis_snapshot.json'
        for invalid in (False, True):
            with self.subTest(invalid=invalid):
                if invalid:
                    self.record['reason'] = ''
                manifest.write_text(json.dumps(self.provenance))
                output, errors = io.StringIO(), io.StringIO()
                with patch.object(checker, 'ROOT', self.root), redirect_stdout(output), redirect_stderr(errors):
                    status = checker.main()
                self.assertEqual(status, int(invalid))
                if invalid:
                    self.assertIn('nonempty reason', errors.getvalue())
                else:
                    self.assertIn('preserved_analysis_files: 1', output.getvalue())
                    self.assertIn('amended_analysis_files: 1', output.getvalue())


if __name__ == '__main__':
    unittest.main()
