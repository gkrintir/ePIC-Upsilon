"""Exercise the grouped retry command with real plans and a simulated scheduler."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
RETRY = ROOT / 'production' / 'retry_rematch_1s_9x275.py'
GENERIC_RETRY = ROOT / 'production' / 'retry_rematch.py'
SPEC = importlib.util.spec_from_file_location('production_jobs', ROOT / 'production' / 'jobs.py')
jobs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(jobs)


class GroupedRetryTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix='upsilon_retry_')
        self.addCleanup(tmp.cleanup)
        self.base = Path(tmp.name)
        self.full = self.base / 'workspaces' / 'condor-1S-9x275-full-grouped-v1'
        self.batch = self.full / 'batch'
        src = self.base / 'Upsilon1S' / '9x275'
        src.mkdir(parents=True)
        inputs = []
        for i in range(6):
            p = src / f'input_{i}.root'
            p.write_text(str(100 + i))
            inputs.append(p)
        self.manifest = self.base / 'input.txt'
        self.manifest.write_text('\n'.join(map(str, inputs)) + '\n')
        self.plan(self.batch)
        self.record = self.full / 'submission.txt'
        self.record.write_text('3 job(s) submitted to cluster 7301.\n')
        self.bin = self.base / 'bin'
        self.bin.mkdir()
        self.calls = self.base / 'submits'
        self.qrows = self.base / 'active.json'
        self.qrows.write_text('[]')
        self.env = {**os.environ, 'PATH': str(self.bin) + os.pathsep + os.environ['PATH'],
                    'RETRY_SUBMITS': str(self.calls), 'RETRY_QUEUE': str(self.qrows)}
        self.stub('condor_q', 'cat "$RETRY_QUEUE"\nexit "${RETRY_QUERY_RC:-0}"')
        self.stub('condor_submit', """echo "$@" >> "$RETRY_SUBMITS"
if [ "${RETRY_UNKNOWN:-0}" = 1 ]; then exit 1; fi
sequence=$(wc -l < "$RETRY_SUBMITS")
echo "2 job(s) submitted to cluster $((7301 + sequence))."
""")

    def plan(self, batch, *extra):
        with redirect_stdout(io.StringIO()):
            jobs.main(['plan', '--input-list', str(self.manifest), '--work-dir', str(batch),
                       '--state', '1S', '--beam', '9x275', '--electron-id', jobs.REMATCH,
                       '--environment-label', 'eic_xl:26.07.1-stable', '--femc', 'full',
                       '--job-flavour', 'longlunch', '--files-per-job', '2', '--write', *extra])

    def stub(self, name, body):
        p = self.bin / name
        p.write_text('#!/bin/sh\n' + body + '\n')
        p.chmod(0o755)

    def unit(self, i):
        return self.batch / 'jobs' / f'job_{i:06d}'

    def complete(self, i, batch=None):
        path = (batch or self.batch) / 'jobs' / f'job_{i:06d}' / 'job.json'
        def run(command, setup, **unused):
            Path(json.loads(path.read_text())['output']).write_text(str(100 + i))
            return SimpleNamespace(returncode=0)
        with patch.object(jobs, 'run_with_setup', side_effect=run), \
                patch.object(jobs, 'root_counts', side_effect=lambda p, _: int(Path(p).read_text())), \
                redirect_stdout(io.StringIO()):
            jobs.worker(SimpleNamespace(job=path, execute=True))

    def retry(self, *args, ok=True):
        r = subprocess.run([sys.executable, str(RETRY), '--batch', str(self.batch), *args],
                           capture_output=True, text=True, env=self.env)
        if ok:
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        else:
            self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        return r

    def test_generic_retry_requires_explicit_batch_and_supports_all_six_cells(self):
        no_batch = subprocess.run([sys.executable, str(GENERIC_RETRY), '--submit'],
                                  capture_output=True, text=True, env=self.env)
        self.assertNotEqual(no_batch.returncode, 0)
        self.assertFalse(self.calls.exists())
        for state in ('1S', '2S', '3S'):
            for beam in ('9x130', '9x275'):
                with self.subTest(state=state, beam=beam):
                    cell = self.base / 'matrix' / f'Upsilon{state}' / beam
                    cell.mkdir(parents=True)
                    src = cell / 'input.root'
                    src.write_text('100')
                    manifest = cell / 'input.txt'
                    manifest.write_text(str(src) + '\n')
                    batch = cell / 'batch'
                    self.plan(batch, '--state', state, '--beam', beam, '--input-list', str(manifest))
                    (cell / 'submission.txt').write_text('1 job(s) submitted to cluster 8001.\n')
                    command = [sys.executable, str(GENERIC_RETRY), '--batch', str(batch)]
                    if (state, beam) == ('3S', '9x275'):
                        command.append('--submit')
                    result = subprocess.run(command, capture_output=True, text=True, env=self.env)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn('Retry: 1 files across 1 Condor groups.', result.stdout)
        calls = self.calls.read_text().splitlines()
        self.assertEqual(len(calls), 1)
        self.assertIn('/Upsilon3S/9x275/batch/retries/attempt-0001/retry.sub', calls[0])

    def test_preview_then_submit_only_unfinished_groups(self):
        for i in (0, 1, 2):
            self.complete(i)
        complete = (self.unit(2) / 'result.json').read_bytes()
        # Failure, interrupted attempt, and an input that was never started.
        failed = self.unit(3)
        for name, data in [('result.json', '{"status": "FAILED"}'), ('output.root', 'partial tuple'),
                           ('root.log', 'ROOT error')]:
            (failed / name).write_text(data)
        interrupted = self.unit(4)
        (interrupted / 'execution.lock').touch()
        (interrupted / 'result.json').write_text('{"status":')
        initial_log = self.batch / 'logs' / 'group_000001.err'
        initial_log.write_text('original error')
        self.retry()
        self.assertFalse(self.calls.exists())
        self.assertEqual((failed / 'output.root').read_text(), 'partial tuple')
        self.assertTrue((interrupted / 'execution.lock').exists())
        self.retry('--submit')
        attempt = self.batch / 'retries' / 'attempt-0001'
        submit = (attempt / 'retry.sub').read_text()
        self.assertEqual([s for s in submit.splitlines() if s.isdigit()], ['000001', '000002'])
        self.assertIn('group_u_CMS.CAF.COMM', submit)
        self.assertIn(str(attempt / 'logs'), submit)
        argv = self.calls.read_text()
        self.assertIn('-maxjobs 2', argv)
        self.assertIn('eic_xl:26.07.1-stable', argv)
        self.assertIn('getenv = False', argv)
        self.assertIn('7302', (attempt / 'submission.txt').read_text())
        self.assertEqual((attempt / 'files' / failed.name / 'output.root').read_text(), 'partial tuple')
        self.assertTrue((attempt / 'files' / interrupted.name / 'execution.lock').is_file())
        self.assertFalse((failed / 'result.json').exists())
        self.assertFalse((interrupted / 'execution.lock').exists())
        self.assertEqual((self.unit(2) / 'result.json').read_bytes(), complete)
        self.assertEqual(initial_log.read_text(), 'original error')

        self.qrows.write_text(json.dumps([{'ClusterId': 7302, 'ProcId': 0, 'JobStatus': 2}]))
        self.retry('--submit', ok=False)
        self.assertEqual(len(self.calls.read_text().splitlines()), 1)
        self.qrows.write_text('[]')
        for i in (3, 4, 5):
            self.complete(i)
        self.retry('--submit')
        self.assertEqual(len(self.calls.read_text().splitlines()), 1)

    def test_verified_outputs_adopted_from_removed_cluster_are_preserved(self):
        old = self.base / 'old-batch'
        with redirect_stdout(io.StringIO()):
            jobs.main(['plan', '--input-list', str(self.base / 'input.txt'), '--work-dir', str(old),
                       '--state', '1S', '--beam', '9x275', '--electron-id', jobs.REMATCH,
                       '--environment-label', 'eic_xl:26.07.1-stable', '--femc', 'full',
                       '--job-flavour', 'longlunch', '--write'])
        for i in (0, 1):
            self.complete(i, old)
        subprocess.run([sys.executable, str(ROOT / 'production' / 'reuse_removed_outputs.py'),
                        str(old), str(self.batch)], check=True, capture_output=True)
        imported = self.unit(0) / 'result.json'
        original_bytes = imported.read_bytes()
        # The first grouped batch can itself be a source for an EOS replan.
        grouped = self.batch
        self.complete(2, grouped)
        self.batch = self.full / 'batch-eos'
        tuples = self.base / 'tuples'
        self.plan(self.batch, '--tuple-output-dir', str(tuples))
        subprocess.run([sys.executable, str(ROOT / 'production' / 'reuse_removed_outputs.py'),
                        str(grouped), str(self.batch)], check=True, capture_output=True)
        migrated = json.loads((self.unit(0) / 'result.json').read_text())
        middle = migrated['reused_from']
        self.assertEqual(migrated['command'], json.loads(original_bytes)['command'])
        self.assertEqual(middle['job'], str(imported.with_name('job.json')))
        self.assertEqual(middle['result_sha256'], jobs.sha256(imported))
        self.assertEqual(middle['reused_from'], json.loads(original_bytes)['reused_from'])
        self.assertEqual((tuples / 'job_000002.root').read_bytes(), (grouped / 'jobs/job_000002/output.root').read_bytes())
        self.assertEqual(list(self.batch.rglob('*.root')), [])
        self.retry('--submit')
        submit = (self.batch / 'retries/attempt-0001/retry.sub').read_text()
        self.assertEqual([s for s in submit.splitlines() if s.isdigit()], ['000001', '000002'])
        self.assertEqual(imported.read_bytes(), original_bytes)

    def test_external_retry_archives_payloads_on_tuple_filesystem(self):
        self.batch = self.full / 'batch-external'
        tuples = self.base / 'tuples'
        self.plan(self.batch, '--tuple-output-dir', str(tuples))
        self.complete(0)
        verified = (tuples / 'job_000000.root').read_bytes()
        (tuples / 'job_000001.root').write_text('orphaned partial')
        (tuples / 'job_000002.root').write_text('failed partial')
        (self.unit(2) / 'result.json').write_text('{"status": "FAILED"}')
        (self.unit(2) / 'root.log').write_text('failed to finish')
        preview = self.retry()
        self.assertIn('Verified 1 DONE. Unfinished inputs: 1 failed, 1 partial, 3 untouched.', preview.stdout)
        self.assertFalse((tuples / '.retries').exists())

        self.retry('--submit')
        attempt = self.batch / 'retries/attempt-0001'
        saved = tuples / '.retries/attempt-0001'
        self.assertEqual((attempt / 'tuple-archive.txt').read_text().strip(), str(saved))
        self.assertEqual((saved / 'job_000001.root').read_text(), 'orphaned partial')
        self.assertEqual((saved / 'job_000002.root').read_text(), 'failed partial')
        self.assertFalse((tuples / 'job_000002.root').exists())
        self.assertEqual((tuples / 'job_000000.root').read_bytes(), verified)
        self.assertEqual(list(self.batch.rglob('*.root')), [])
        self.assertEqual((attempt / 'files/job_000002/root.log').read_text(), 'failed to finish')

    def test_serial_and_parallel_retry_scans_select_the_same_files(self):
        for i in (0, 2, 3):
            self.complete(i)
        (self.unit(4) / 'execution.lock').touch()
        (self.unit(1) / 'result.json').write_text('{"status":"FAILED"}')
        serial = self.retry('--check-workers', '1').stdout
        self.assertEqual(self.retry('--check-workers', '3').stdout, serial)
        self.assertFalse(self.calls.exists())
        self.retry('--check-workers', '0', '--submit', ok=False)
        self.assertFalse((self.batch / 'retries').exists())
        self.retry('--check-workers', '3', '--submit')
        pending = (self.batch / 'retries/attempt-0001/file_jobs.txt').read_text().splitlines()
        self.assertEqual(pending, [str(self.unit(i) / 'job.json') for i in (1, 4, 5)])

    def test_a_second_failed_retry_keeps_the_preceding_attempt(self):
        for i in range(5):
            self.complete(i)
        failed = self.unit(5) / 'root.log'
        failed.write_text('initial error')
        self.retry('--submit')
        first = self.batch / 'retries/attempt-0001'
        saved = (first / 'submission.txt').read_bytes()
        failed.write_text('new attempt error')
        self.retry('--submit')
        second = self.batch / 'retries/attempt-0002'
        self.assertEqual((first / 'files/job_000005/root.log').read_text(), 'initial error')
        self.assertEqual((second / 'files/job_000005/root.log').read_text(), 'new attempt error')
        self.assertEqual((first / 'submission.txt').read_bytes(), saved)
        self.assertEqual(len(self.calls.read_text().splitlines()), 2)

    def test_held_or_removing_cluster_prevents_state_changes(self):
        lock = self.unit(0) / 'execution.lock'
        lock.touch()
        for status in (5, 3):
            with self.subTest(status=status):
                self.qrows.write_text(json.dumps([{'ClusterId': 7301, 'ProcId': 0, 'JobStatus': status,
                                                   'HoldReason': 'fixture'}]))
                self.retry('--submit', ok=False)
                self.assertTrue(lock.is_file())
                self.assertFalse(self.calls.exists())

    def test_invalid_done_stops_before_other_evidence_is_archived(self):
        failed = self.unit(0) / 'result.json'
        failed.write_text('{"status":"FAILED"}')
        self.complete(5)
        corrupt = self.unit(5) / 'output.root'
        corrupt.write_text('damaged')
        self.retry('--submit', ok=False)
        self.assertEqual(corrupt.read_text(), 'damaged')
        self.assertTrue(failed.is_file())
        self.assertFalse(self.calls.exists())

    def test_inconclusive_queue_or_submission_is_not_retried_automatically(self):
        self.env['RETRY_QUERY_RC'] = '1'
        self.retry('--submit', ok=False)
        self.assertFalse(self.calls.exists())
        self.env.pop('RETRY_QUERY_RC')
        self.record.write_text('Submitting jobs...\n')
        self.retry('--submit', ok=False)
        self.assertFalse(self.calls.exists())
        self.record.write_text('3 job(s) submitted to cluster 7301.\n')
        self.env['RETRY_UNKNOWN'] = '1'
        self.retry('--submit', ok=False)
        self.env.pop('RETRY_UNKNOWN')
        self.retry('--submit', ok=False)
        self.assertEqual(len(self.calls.read_text().splitlines()), 1)


if __name__ == '__main__':
    unittest.main()
