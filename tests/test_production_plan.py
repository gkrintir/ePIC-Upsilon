"""Plan/ownership regressions; no ROOT executable or Condor submission required."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
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
            self.assertLess(lines.index(setting), lines.index('queue 1'))

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
        self.assertEqual((self.work / 'pending.sub').read_text().count('queue 1'), 1)
        with redirect_stdout(io.StringIO()), self.assertRaises(FileExistsError):
            jobs.pending(SimpleNamespace(batch=self.work, write=True))

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
