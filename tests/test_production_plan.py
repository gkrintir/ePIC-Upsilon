"""Plan/ownership regressions; no ROOT executable or Condor submission required."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('production_jobs', ROOT / 'production/jobs.py')
jobs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(jobs)
HAS_UPROOT = importlib.util.find_spec('uproot') is not None


class ProductionPlanTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='upsilon_production_')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.input = self.base / 'Upsilon1S/input.root'
        self.input.parent.mkdir()
        self.input.write_bytes(b'planning does not open ROOT')
        self.inputs = self.base / 'inputs.txt'
        self.inputs.write_text(str(self.input) + '\n')
        self.work = self.base / 'batch'

    def args(self, *extra):
        return jobs.parser().parse_args([
            'plan', '--input-list', str(self.inputs), '--work-dir', str(self.work),
            '--beam', '9x275', '--state', '1S', '--electron-id', 'reco',
            '--environment-label', 'test-only-no-ROOT', *extra])

    def create(self, *extra):
        with redirect_stdout(io.StringIO()):
            jobs.plan(self.args('--write', '--hash-inputs', *extra))
        return self.work / 'jobs/job_000000/job.json'

    def test_preserved_sources_hashes(self):
        source = ROOT / 'production'
        for record in json.loads((source / 'sources.json').read_text())['files']:
            self.assertEqual(jobs.sha256(source / record['path']), record['sha256'])
        for beam in ('9x130', '9x275'):
            jobs.verify_macro_sources(beam)

    def test_plan_is_read_only(self):
        with patch.object(jobs.subprocess, 'run') as execute, redirect_stdout(io.StringIO()):
            self.assertEqual(jobs.plan(self.args()), 0)
        execute.assert_not_called()
        self.assertFalse(self.work.exists())

    def make_counted_inputs(self, count):
        sources = []
        for i in range(count):
            source = self.input.with_name(f'events-{i}.root')
            source.write_text(str(i + 1))
            sources.append(source)
        self.inputs.write_text('\n'.join(map(str, sources)) + '\n')
        return sources

    def group(self, index, *options):
        return jobs.main(['group-worker', '--batch', str(self.work), '--index', str(index), *options])

    def test_grouped_plan_validates_size_and_queues_stable_indices(self):
        self.make_counted_inputs(5)
        for invalid in ('0', '-1'):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                jobs.plan(self.args('--files-per-job', invalid, '--write'))
            self.assertFalse(self.work.exists())
        self.create('--files-per-job', '2')
        plan = json.loads((self.work / 'plan.json').read_text())
        self.assertEqual((plan['job_count'], plan['files_per_job'], plan['condor_job_count']), (5, 2, 3))
        self.assertEqual(len(plan['jobs']), 5)
        submit = (self.work / 'submit.sub').read_text()
        self.assertTrue(submit.endswith('queue batch_index from (\n000000\n000001\n000002\n)\n'))
        self.assertIn('arguments = $(batch_index)', submit)
        for index in (-1, 3):
            with self.subTest(index=index), self.assertRaises(ValueError):
                self.group(index)

    def test_group_workers_checkpoint_each_input_and_merge_in_input_order(self):
        sources = self.make_counted_inputs(5)
        self.create('--files-per-job', '2')
        records = [Path(p) for p in json.loads((self.work / 'plan.json').read_text())['jobs']]
        outputs = [p.with_name('output.root') for p in records]
        calls = []

        # Use the actual worker, locks and receipts. Only ROOT I/O is replaced
        # by a small count fixture with a distinct count for each input.
        def fake_run(command, setup, **kwargs):
            if command[0] == 'root':
                spec = json.loads((Path(kwargs['cwd']) / 'job.json').read_text())
                source = Path(spec['inputs'][0]['path'])
                calls.append(source)
                Path(spec['output']).write_bytes(source.read_bytes())
            else:
                self.assertEqual(command[0], 'hadd')
                inputs = [Path(p) for p in Path(command[2][1:]).read_text().splitlines()]
                self.assertEqual(inputs, outputs)
                Path(command[1]).write_text(str(sum(int(p.read_text()) for p in inputs)))
            return SimpleNamespace(returncode=0)

        with patch.object(jobs, 'root_counts', side_effect=lambda p, _: int(Path(p).read_text())), \
             patch.object(jobs, 'run_with_setup', side_effect=fake_run), redirect_stdout(io.StringIO()):
            self.group(0)
            self.assertEqual(calls, [])
            # A previously finished individual file is part of an incomplete group.
            jobs.worker(SimpleNamespace(job=records[0], execute=True))
            self.group(1, '--execute')
            jobs.pending(SimpleNamespace(batch=self.work, write=True))
            self.assertTrue((self.work / 'pending.sub').read_text().endswith(
                'queue batch_index from (\n000000\n000002\n)\n'))
            self.group(2, '--execute')  # final group contains only one file
            self.group(0, '--execute')
            receipts = [p.with_name('result.json').read_bytes() for p in records]
            self.group(0, '--execute')
            self.assertEqual([p.with_name('result.json').read_bytes() for p in records], receipts)
            self.assertEqual(calls, [sources[i] for i in (0, 2, 3, 4, 1)])
            merged = self.base / 'merged.root'
            jobs.merge(SimpleNamespace(batch=self.work, output=merged, execute=True))
            self.assertEqual(merged.read_text(), '15')
            # A corrupt checkpoint must not count as DONE or get overwritten.
            outputs[-1].write_text('99')
            with self.assertRaises((ValueError, RuntimeError)):
                self.group(2, '--execute')
            self.assertEqual(outputs[-1].read_text(), '99')
            self.assertEqual(len(calls), 5)

    def test_group_preserves_locked_and_changed_files_and_finishes_independent_input(self):
        sources = self.make_counted_inputs(3)
        self.create('--files-per-job', '3')
        records = [Path(p) for p in json.loads((self.work / 'plan.json').read_text())['jobs']]
        sources[0].write_text('different input')
        lock = records[1].with_name('execution.lock')
        lock.write_text('held by another worker')

        def produce(command, setup, **kwargs):
            self.assertEqual(Path(kwargs['cwd']), records[2].parent)
            records[2].with_name('output.root').write_text('3')
            return SimpleNamespace(returncode=0)

        with patch.object(jobs, 'root_counts', side_effect=lambda p, _: int(Path(p).read_text())), \
             patch.object(jobs, 'run_with_setup', side_effect=produce) as executed, \
             redirect_stdout(io.StringIO()), self.assertRaises(RuntimeError):
            self.group(0, '--execute')
        self.assertEqual(executed.call_count, 1)
        self.assertFalse(records[0].with_name('output.root').exists())
        self.assertFalse(records[1].with_name('output.root').exists())
        self.assertEqual(lock.read_text(), 'held by another worker')
        self.assertEqual(json.loads(records[2].with_name('result.json').read_text())['status'], 'DONE')

    def test_rematching_freezes_corrected_reader_and_archived_input(self):
        for beam in ('9x130', '9x275'):
            with self.subTest(beam=beam):
                self.work = self.base / ('rematch-' + beam)
                path = self.create('--beam', beam, '--electron-id', jobs.REMATCH)
                job = jobs.load_job(path)
                macros = Path(job['macro_dir'])
                self.assertEqual((macros / jobs.PRODUCER).read_bytes(),
                                 (jobs.PRODUCTION / jobs.PRODUCERS[beam]).read_bytes())
                self.assertEqual((macros / jobs.ARCHIVED_READER).read_bytes(),
                                 (jobs.PRODUCTION / 'macros/ElectronFinder.C').read_bytes())
                self.assertNotEqual((macros / jobs.ARCHIVED_READER).read_bytes(),
                                    (macros / 'ElectronFinder.C').read_bytes())
                self.assertIn(str(macros / jobs.REMATCH_ENTRY) + '(', jobs.root_command(job)[-1])
                self.assertEqual(job['electron_id'], jobs.REMATCH)
                wrapper = (macros / jobs.REMATCH_ENTRY).read_text()
                self.assertLess(wrapper.index('#define UPSILON_REMATCH_3011'), wrapper.index('#include'))
                self.assertIn('femc_eta_split, "reco")', wrapper)
                with (macros / jobs.REMATCH_ENTRY).open('a') as output:
                    output.write('// changed\n')
                with self.assertRaisesRegex(ValueError, 'snapshot hash mismatch'):
                    jobs.load_job(path)

    def test_preserved_electron_ids_do_not_apply_the_patch(self):
        for mode in ('legacy', 'reco'):
            with self.subTest(mode=mode), patch.object(jobs, 'build_rematched_reader') as rematch:
                self.work = self.base / mode
                path = self.create('--electron-id', mode)
                job = jobs.load_job(path)
                self.assertTrue(jobs.root_command(job)[-1].endswith(',"' + mode + '")'))
                self.assertFalse((Path(job['macro_dir']) / jobs.REMATCH_ENTRY).exists())
                rematch.assert_not_called()

    def test_rematch_condor_accounting_in_initial_and_pending_submissions(self):
        self.create('--electron-id', jobs.REMATCH)
        account = 'group_u_CMS.CAF.COMM'
        self.assertEqual(json.loads((self.work / 'plan.json').read_text())['accounting_group'], account)
        with redirect_stdout(io.StringIO()):
            jobs.pending(SimpleNamespace(batch=self.work, write=True))
        for filename in ('submit.sub', 'pending.sub'):
            lines = (self.work / filename).read_text().splitlines()
            setting = f'+AccountingGroup = "{account}"'
            self.assertEqual(lines.count(setting), 1)
            self.assertLess(lines.index(setting), lines.index('queue batch_index from ('))

    def test_unsupported_id_and_batch_combinations_fail_before_writes(self):
        for extra in [('--femc', 'keep_xpos'), ('--femc', 'keep_yneg'),
                      ('--job-flavour', 'bad\nqueue 50')]:
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                jobs.plan(self.args('--write', *extra))
            self.assertFalse(self.work.exists())

    def test_duplicate_and_wrong_state_inputs_rejected(self):
        self.inputs.write_text((str(self.input) + '\n') * 2)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            jobs.plan(self.args('--write'))
        self.inputs.write_text(str(self.input) + '\n')
        with self.assertRaisesRegex(ValueError, 'identify 2S'):
            jobs.plan(self.args('--write', '--state', '2S'))
        self.assertFalse(self.work.exists())

    def test_beam_labels_must_match_before_planning(self):
        for label in ('9x130', '9To130', '10on100'):
            with self.subTest(label=label):
                file = self.input.with_name(f'Upsilon1S_{label}_q2_0to0.01.root')
                file.touch()
                self.inputs.write_text(str(file) + '\n')
                with self.assertRaisesRegex(ValueError, 'conflicts with 9x275'):
                    jobs.plan(self.args('--write'))
                self.assertFalse(self.work.exists())

    def test_130_job_freezes_dedicated_macro_under_original_name(self):
        input130 = self.input.with_name('Upsilon1S_9x130_q2_0to1.root')
        self.input.rename(input130)
        self.inputs.write_text(str(input130) + '\n')
        path = self.create('--beam', '9x130')
        job = jobs.load_job(path)
        source = ROOT / 'production/macros/9x130/DVMP_Upsilon_Analysis.C'
        snapshot = Path(job['macro_dir'])
        self.assertEqual(job['beam'], '9x130')
        self.assertEqual((snapshot / jobs.PRODUCER).read_bytes(), source.read_bytes())
        self.assertEqual(job['source_hashes'][jobs.PRODUCER], jobs.sha256(source))
        self.assertNotEqual(job['source_hashes'][jobs.PRODUCER], jobs.sha256(source.parent.parent / jobs.PRODUCER))
        for helper in ('ElectronFinder.C', 'ePICStyle.h', 'ePICStyle.C'):
            self.assertEqual((snapshot / helper).read_bytes(), (source.parent.parent / helper).read_bytes())
        with patch.object(jobs.subprocess, 'run') as execute, redirect_stdout(io.StringIO()):
            jobs.worker(SimpleNamespace(job=path, execute=False))
        execute.assert_not_called()
        self.assertIn(str(snapshot / jobs.PRODUCER), jobs.root_command(job)[-1])

        # A relabelled plan cannot run a frozen source with the other energy.
        for record in (path, self.work / 'plan.json'):
            spec = json.loads(record.read_text())
            spec['beam'] = '9x275'
            record.write_text(json.dumps(spec))
        with self.assertRaisesRegex(ValueError, 'beam constant E_p_beam conflicts'):
            jobs.load_job(path)

    def test_write_freezes_sources_and_worker_dryrun_does_not_execute(self):
        path = self.create()
        job = jobs.load_job(path)
        for name in jobs.MACROS:
            self.assertEqual((self.work / 'snapshot/macros' / name).read_bytes(),
                             (ROOT / 'production/macros' / name).read_bytes())
        with patch.object(jobs.subprocess, 'run') as execute, redirect_stdout(io.StringIO()):
            jobs.worker(SimpleNamespace(job=path, execute=False))
        execute.assert_not_called()
        self.assertFalse(Path(job['output']).exists())
        command = jobs.root_command(job)
        self.assertEqual(command[:4], ['root', '-l', '-b', '-q'])
        self.assertTrue(command[-1].endswith(',"reco")'))
        self.assertIn('should_transfer_files = NO', (self.work / 'submit.sub').read_text())

    def test_existing_work_is_untouched(self):
        self.work.mkdir()
        protected = self.work / 'keep.txt'
        protected.write_text('old result')
        with self.assertRaises(FileExistsError):
            jobs.plan(self.args('--write'))
        self.assertEqual(protected.read_text(), 'old result')
        self.assertEqual(list(self.work.iterdir()), [protected])

    def test_relocated_condor_launcher_reaches_frozen_worker(self):
        job = self.create()
        slot = self.base / 'execute'
        slot.mkdir()
        executable = slot / 'condor_exec.exe'
        shutil.copy2(self.work / 'snapshot/run_job.sh', executable)
        # Exercise the real worker up to its overwrite guard without invoking ROOT.
        output = job.parent / 'output.root'
        output.write_bytes(b'keep this tuple')
        result = subprocess.run([str(executable), str(job)], cwd=slot, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Existing output/result is preserved', result.stderr)
        self.assertEqual(output.read_bytes(), b'keep this tuple')

    def test_changed_inputs_or_macros_block_execution(self):
        path = self.create()
        self.input.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'Input changed'):
            jobs.load_job(path)
        macro = self.work / 'snapshot/macros/ElectronFinder.C'
        macro.write_text(macro.read_text() + '\n// local edit\n')
        with self.assertRaisesRegex(ValueError, 'snapshot hash mismatch'):
            jobs.load_job(path)

    def test_pending_never_submits_and_preserves_its_file(self):
        self.create()
        with patch.object(jobs.subprocess, 'run') as execute, redirect_stdout(io.StringIO()):
            jobs.pending(SimpleNamespace(batch=self.work, write=True))
        execute.assert_not_called()
        self.assertTrue((self.work / 'pending.sub').read_text().endswith(
            'queue batch_index from (\n000000\n)\n'))
        with redirect_stdout(io.StringIO()), self.assertRaises(FileExistsError):
            jobs.pending(SimpleNamespace(batch=self.work, write=True))

    def test_pending_queue_preserves_original_nonconsecutive_job_indices(self):
        inputs = [self.input, self.input.with_name('second.root'), self.input.with_name('third.root')]
        for source in inputs[1:]:
            source.write_bytes(b'planning input')
        self.inputs.write_text('\n'.join(map(str, inputs)) + '\n')
        self.create()
        completed = self.work / 'jobs/job_000001'
        output = completed / 'output.root'
        output.write_bytes(b'completed output')
        (completed / 'result.json').write_text(json.dumps(
            {'status': 'DONE', 'output_sha256': jobs.sha256(output)}))
        with redirect_stdout(io.StringIO()):
            jobs.pending(SimpleNamespace(batch=self.work, write=True))
        full_file = (self.work / 'submit.sub').read_text()
        pending_file = (self.work / 'pending.sub').read_text()
        for text, ids in [(full_file, ['000000', '000001', '000002']),
                          (pending_file, ['000000', '000002'])]:
            header, rows = text.split('queue batch_index from (\n')
            self.assertEqual(rows.splitlines(), [*ids, ')'])
            self.assertIn(f'initialdir = {self.work}/jobs/job_$(batch_index)\n', header)
            self.assertIn(f'arguments = {self.work}/jobs/job_$(batch_index)/job.json\n', header)
            for field, suffix in [('output', 'out'), ('error', 'err'), ('log', 'log')]:
                self.assertIn(f'{field} = {self.work}/logs/$(batch_index).{suffix}\n', header)

    def test_job_cannot_drift_from_batch_selection(self):
        path = self.create()
        data = json.loads(path.read_text())
        data['electron_id'] = 'legacy'
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'batch configuration: electron_id'):
            jobs.load_job(path)

    @unittest.skipUnless(HAS_UPROOT, 'uproot needed for synthetic ROOT bookkeeping checks')
    def test_worker_and_merge_bookkeeping_with_synthetic_ROOT(self):
        import numpy as np
        import uproot
        with uproot.recreate(self.input) as output:
            output.mktree('events', {'event': 'int64'})
            output['events'].extend({'event': np.arange(4)})
        path = self.create('--beam', '9x130')
        job = jobs.load_job(path)

        def write_tuple(destination, count=4, generated=4):
            with uproot.recreate(destination) as output:
                output.mktree('UpsilonTree', {'event': 'int64', 'weight': 'float64'})
                output['UpsilonTree'].extend({'event': np.arange(count), 'weight': np.full(count, 5.)})
                output['h_EventCount'] = (np.array([float(generated)]), np.array([0., 1.]))

        def fake_root(command, setup, **kwargs):
            self.assertEqual(command[0], 'root')
            write_tuple(job['output'])
            return SimpleNamespace(returncode=0)

        with patch.object(jobs, 'run_with_setup', side_effect=fake_root), redirect_stdout(io.StringIO()):
            jobs.worker(SimpleNamespace(job=path, execute=True))
        result = json.loads((path.parent / 'result.json').read_text())
        self.assertEqual(result['status'], 'DONE')
        self.assertEqual(result['output_rows'], 4)
        with redirect_stdout(io.StringIO()):
            jobs.pending(SimpleNamespace(batch=self.work, write=True))
        self.assertFalse((self.work / 'pending.sub').exists())
        with redirect_stdout(io.StringIO()), self.assertRaises(FileExistsError):
            jobs.worker(SimpleNamespace(job=path, execute=True))

        out = self.base / 'merged.root'
        with patch.object(jobs, 'run_with_setup') as execute, redirect_stdout(io.StringIO()):
            jobs.merge(SimpleNamespace(batch=self.work, output=out, execute=False))
        execute.assert_not_called()
        self.assertFalse(out.exists())
        with self.assertRaisesRegex(ValueError, 'conflicts with 9x130'):
            jobs.merge(SimpleNamespace(batch=self.work, output=self.base / '9x275.root', execute=False))

        def fake_hadd(command, setup, **kwargs):
            self.assertEqual(command[0:2], ['hadd', str(out)])
            self.assertEqual(command[2], '@' + str(out) + '.inputs.txt')
            self.assertEqual(Path(command[2][1:]).read_text().splitlines(), [job['output']])
            self.assertNotIn('-f', command)
            write_tuple(out)
            return SimpleNamespace(returncode=0)

        with patch.object(jobs, 'run_with_setup', side_effect=fake_hadd), redirect_stdout(io.StringIO()):
            jobs.merge(SimpleNamespace(batch=self.work, output=out, execute=True))
        merge_record = json.loads(out.with_name(out.name + '.merge.json').read_text())
        self.assertEqual(merge_record['status'], 'DONE')
        self.assertEqual(merge_record['beam'], '9x130')
        self.assertEqual(merge_record['source_hashes'], job['source_hashes'])
        with uproot.open(out) as source:
            self.assertEqual(source['UpsilonTree']['weight'].array(library='np').sum(), 20.)
        with redirect_stdout(io.StringIO()), self.assertRaises(FileExistsError):
            jobs.merge(SimpleNamespace(batch=self.work, output=out, execute=True))
        write_tuple(self.base / 'bad.root', generated=5)
        with self.assertRaisesRegex(ValueError, 'h_EventCount'):
            jobs.root_counts(self.base / 'bad.root', 'UpsilonTree')

    @unittest.skipUnless(HAS_UPROOT, 'uproot needed for synthetic ROOT bookkeeping checks')
    def test_lost_rows_mark_job_failed_and_prevent_merge(self):
        import numpy as np
        import uproot
        with uproot.recreate(self.input) as output:
            output.mktree('events', {'event': 'int64'})
            output['events'].extend({'event': np.arange(4)})
        path = self.create()
        job = jobs.load_job(path)

        def fake_root(command, setup, **kwargs):
            with uproot.recreate(job['output']) as output:
                output.mktree('UpsilonTree', {'event': 'int64'})
                output['UpsilonTree'].extend({'event': np.arange(3)})
                output['h_EventCount'] = (np.array([3.]), np.array([0., 1.]))
            return SimpleNamespace(returncode=0)

        with patch.object(jobs, 'run_with_setup', side_effect=fake_root), redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, 'Generated rows were lost'):
                jobs.worker(SimpleNamespace(job=path, execute=True))
        self.assertEqual(json.loads((path.parent / 'result.json').read_text())['status'], 'FAILED')
        self.assertTrue(Path(job['output']).exists())
        self.assertFalse((path.parent / 'execution.lock').exists())
        with self.assertRaisesRegex(ValueError, 'validated DONE'):
            jobs.merge(SimpleNamespace(batch=self.work, output=self.base / 'merged.root', execute=True))
        self.assertFalse((self.base / 'merged.root').exists())


if __name__ == '__main__':
    unittest.main()
