"""Regression checks for carrying a removed cluster's validated tuples into a grouped plan."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('production_jobs', ROOT / 'production/jobs.py')
jobs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(jobs)


class RemovedBatchReuseTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='upsilon_removed_')
        self.addCleanup(temp.cleanup)
        self.base = Path(temp.name)
        inputs = self.base / 'Upsilon1S/9x275'
        inputs.mkdir(parents=True)
        self.inputs = [inputs / f'file_{i}.root' for i in range(5)]
        for i, path in enumerate(self.inputs):
            path.write_text(str(i + 1))
        self.manifest = self.base / 'inputs.txt'
        self.manifest.write_text('\n'.join(map(str, self.inputs)) + '\n')
        self.old = self.base / 'old'
        self.new = self.base / 'new'
        self.plan(self.old, 1)
        self.plan(self.new, 2)

    def plan(self, work, width, *extra):
        with redirect_stdout(io.StringIO()):
            jobs.main(['plan', '--input-list', str(self.manifest), '--work-dir', str(work),
                       '--beam', '9x275', '--state', '1S', '--electron-id', jobs.REMATCH,
                       '--files-per-job', str(width), '--environment-label', 'pilot-test',
                       '--job-flavour', 'longlunch', '--write', *extra])

    @staticmethod
    def directory(batch, index):
        return batch / 'jobs' / f'job_{index:06d}'

    def complete_old(self, index):
        # Run the actual frozen-job validation and receipt-writing path. Only ROOT
        # input/output is replaced by a small scalar fixture with a per-file count.
        job = self.directory(self.old, index) / 'job.json'
        def run(command, setup, **kwargs):
            (job.parent / 'output.root').write_bytes(self.inputs[index].read_bytes())
            return SimpleNamespace(returncode=0)
        with patch.object(jobs, 'run_with_setup', side_effect=run), \
                patch.object(jobs, 'root_counts', side_effect=lambda p, _: int(Path(p).read_text())), \
                redirect_stdout(io.StringIO()):
            jobs.worker(SimpleNamespace(job=job, execute=True))

    def adopt(self, expect_success=True):
        result = subprocess.run([sys.executable, str(ROOT / 'production/reuse_removed_outputs.py'),
                                 str(self.old), str(self.new)], capture_output=True, text=True)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def queued_groups(self):
        lines = (self.new / 'after-removal.sub').read_text().splitlines()
        return [line for line in lines if line.isdigit()]

    def test_no_completed_jobs_queues_all_groups(self):
        self.adopt()
        self.assertEqual(self.queued_groups(), ['000000', '000001', '000002'])
        self.assertEqual(list(self.new.glob('jobs/*/result.json')), [])

    def test_mixed_results_reuse_and_idempotency(self):
        # Two fully completed inputs fill group 0. Group 1 is partly complete.
        for i in (0, 1, 3):
            self.complete_old(i)
        self.directory(self.old, 2).joinpath('result.json').write_text('{"status":"FAILED"}')
        self.directory(self.old, 4).joinpath('execution.lock').touch()
        self.directory(self.old, 4).joinpath('output.root').write_bytes(b'incomplete')
        original = {p: p.read_bytes() for p in self.old.glob('jobs/*/*')}
        self.adopt()
        self.assertEqual(self.queued_groups(), ['000001', '000002'])
        for i in (0, 1, 3):
            old_dir, new_dir = self.directory(self.old, i), self.directory(self.new, i)
            result = json.loads((new_dir / 'result.json').read_text())
            self.assertEqual(result['job'], str(new_dir / 'job.json'))
            self.assertEqual(result['output'], str(new_dir / 'output.root'))
            self.assertEqual(result['reused_from']['job'], str(old_dir / 'job.json'))
            self.assertEqual(result['reused_from']['result_sha256'], jobs.sha256(old_dir / 'result.json'))
            self.assertEqual((old_dir / 'output.root').read_bytes(), (new_dir / 'output.root').read_bytes())
            with patch.object(jobs, 'root_counts', return_value=i + 1):
                self.assertEqual(jobs.verified_output_rows(jobs.load_job(new_dir / 'job.json')), i + 1)
        for i in (2, 4):
            self.assertFalse((self.directory(self.new, i) / 'result.json').exists())
            self.assertFalse((self.directory(self.new, i) / 'output.root').exists())
        adopted = {p: p.read_bytes() for p in self.new.glob('jobs/*/result.json')}
        self.adopt()
        self.assertEqual({p: p.read_bytes() for p in adopted}, adopted)
        self.assertEqual({p: p.read_bytes() for p in original}, original)

    def test_all_completed_leaves_no_condor_jobs(self):
        for i in range(len(self.inputs)):
            self.complete_old(i)
        self.adopt()
        self.assertEqual(self.queued_groups(), [])
        self.assertEqual(len(list(self.new.glob('jobs/*/result.json'))), len(self.inputs))

    def test_corrupt_done_refuses_reuse(self):
        self.complete_old(0)
        output = self.directory(self.old, 0) / 'output.root'
        output.write_text('different')
        self.adopt(False)
        self.assertEqual(output.read_text(), 'different')
        self.assertFalse((self.directory(self.new, 0) / 'result.json').exists())

    def test_foreign_target_output_is_preserved(self):
        self.complete_old(0)
        output = self.directory(self.new, 0) / 'output.root'
        output.write_text('existing target')
        self.adopt(False)
        self.assertEqual(output.read_text(), 'existing target')

    def test_interrupted_import_finishes_without_rerunning_root(self):
        self.complete_old(0)
        self.adopt()
        dst = self.directory(self.new, 0)
        wanted = (dst / 'result.json').read_bytes()
        for has_output in (True, False):
            with self.subTest(has_output=has_output):
                (dst / 'result.json').rename(dst / 'reuse-pending.json')
                if not has_output:
                    (dst / 'output.root').unlink()
                self.adopt()
                self.assertEqual((dst / 'result.json').read_bytes(), wanted)
                self.assertFalse((dst / 'reuse-pending.json').exists())

    def test_receipt_with_other_input_command_is_rejected(self):
        self.complete_old(0)
        receipt = self.directory(self.old, 0) / 'result.json'
        result = json.loads(receipt.read_text())
        result['command'][-1] = result['command'][-1].replace(str(self.inputs[0]), str(self.inputs[1]))
        receipt.write_text(json.dumps(result))
        self.adopt(False)
        self.assertFalse((self.directory(self.new, 0) / 'result.json').exists())

    def test_different_configuration_and_inventory_are_refused(self):
        new = self.new
        self.new = self.base / 'wrong-configuration'
        self.plan(self.new, 2, '--eta-split', '1.5')
        self.adopt(False)
        self.new = self.base / 'missing-input'
        self.manifest.write_text('\n'.join(map(str, self.inputs[:-1])) + '\n')
        self.plan(self.new, 2)
        self.adopt(False)
        self.new = new


if __name__ == '__main__':
    unittest.main()
