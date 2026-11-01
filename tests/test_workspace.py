"""Fresh temporary workspaces exercise the preserved backend without CERN ROOT files."""
import argparse
from contextlib import redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('workspace',REPO/'tools/workspace.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)


class WorkspaceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.repo=Path(self.temp.name)/'repo';self.repo.mkdir()
        shutil.copytree(REPO/'analysis',self.repo/'analysis',ignore=shutil.ignore_patterns('__pycache__'))
        self.settings={'workspace_root':str(self.repo/'workspaces'),'input_base':str(self.repo/'missing-eos')}

    def init(self,selection='proton-tagged',name='test'):
        with redirect_stdout(io.StringIO()):return w.init_workspace(self.repo,self.settings,name,selection)

    def test_initialization_selection_and_no_overwrite(self):
        tagged=self.init()
        _,meta=w.load_workspace(self.repo,self.settings,'test')
        self.assertEqual(meta['selection'],'proton-tagged')
        self.assertEqual(json.loads((tagged/'.proton_category.json').read_text())['category_root'],str(tagged))
        with self.assertRaises(FileExistsError):self.init()
        inclusive=self.init('inclusive','inclusive')
        self.assertFalse((inclusive/'.proton_category.json').exists())
        w.load_workspace(self.repo,self.settings,'inclusive')
        self.assertFalse((tagged/'results').exists())

    def test_dry_init_writes_nothing(self):
        with redirect_stdout(io.StringIO()):w.init_workspace(self.repo,self.settings,'test','inclusive',True)
        self.assertFalse((self.repo/'workspaces').exists())

    def test_mutated_workspace_rejected(self):
        workspace=self.init()
        (workspace/'workflow/run_all_states.py').write_text('changed')
        with self.assertRaisesRegex(ValueError,'source changed'):w.load_workspace(self.repo,self.settings,'test')
        # Status remains available read-only after code drift.
        w.load_workspace(self.repo,self.settings,'test',verify=False)

    def test_mutated_template_and_path_escape_rejected(self):
        self.init()
        (self.repo/'analysis/workflow/run_all_states.py').write_text('changed')
        with self.assertRaisesRegex(ValueError,'source changed'):w.load_workspace(self.repo,self.settings,'test')
        with self.assertRaises(ValueError):w.workspace_path(self.repo,self.settings,'../escape')
        bad={**self.settings,'workspace_root':str(self.repo/'analysis')}
        with self.assertRaises(ValueError):w.workspace_path(self.repo,bad,'test')

    def test_symlink_workspace_rejected(self):
        self.init();link=self.repo/'workspaces/link';link.symlink_to(self.repo/'workspaces/test',target_is_directory=True)
        with self.assertRaisesRegex(ValueError,'Symlink'):w.load_workspace(self.repo,self.settings,'link')

    def test_real_backend_plans_both_categories(self):
        import subprocess,sys
        for selection in ('inclusive','proton-tagged'):
            workspace=self.init(selection,selection)
            args=argparse.Namespace(states=['1S'],batch_id='smoke-'+selection,action='plan',energy='9x275',stat_mode='projected',lumi=2.5,hash_inputs=False)
            command=w.analysis_command(args,self.settings,workspace)
            env={**os.environ,'YPSILON_PROTON_SELECTION':w.SELECTIONS[selection],'PYTHONDONTWRITEBYTECODE':'1'}
            result=subprocess.run(command,cwd=workspace,env=env,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('prepare_9x275',result.stdout)
            self.assertIn(str(self.repo/'missing-eos'),result.stdout)
            self.assertFalse((workspace/'results').exists())
            w.load_workspace(self.repo,self.settings,selection)

    def test_unavailable_status_is_reported(self):
        workspace=self.init()
        run=workspace/'results/generated/runs/example'
        (run/'artifacts/reoptimization').mkdir(parents=True)
        (run/'run_manifest.json').write_text(json.dumps({'status':'COMPLETE'}))
        (run/'artifacts/reoptimization/validation_status.csv').write_text('state,observable,overall_status,inference_status\n1S,pt2,NO_FEASIBLE_SLOPE_BINNING,NOT_RUN\n')
        buf=io.StringIO()
        with redirect_stdout(buf):w.show_status(workspace,'example')
        self.assertIn('NO_FEASIBLE_SLOPE_BINNING (inference NOT_RUN)',buf.getvalue())
        self.assertIn('does not mean physics validation passed',buf.getvalue())

    def test_reopt_unrelated_output_is_preserved(self):
        workspace=self.init()
        out=workspace/'results/generated/runs/existing'
        out.mkdir(parents=True);(out/'keep.txt').write_text('keep')
        args=argparse.Namespace(states=['1S'],run_id='existing',action='run')
        with self.assertRaisesRegex(FileExistsError,'nonempty'):
            w.reopt_command(args,self.repo,workspace,{'selection':'proton-tagged'})
        self.assertEqual((out/'keep.txt').read_text(),'keep')

    def test_import_check_detects_discoverable_broken_scipy(self):
        # A package spec exists, but its compiled dependency fails when linalg
        # loads. Checking find_spec('scipy') alone would miss this failure.
        fake=Path(self.temp.name)/'fake-dependencies'
        scipy=fake/'scipy';scipy.mkdir(parents=True)
        (scipy/'__init__.py').write_text("__version__='synthetic-test'\n")
        (scipy/'linalg.py').write_text("raise ImportError('libopenblas.so: cannot read file data')\n")
        (scipy/'optimize.py').write_text('')
        with patch.dict(os.environ,{'PYTHONPATH':str(fake)}):
            ok, version=w.import_check('scipy')
            self.assertTrue(ok);self.assertEqual(version,'synthetic-test')
            ok, error=w.import_check('scipy',('scipy.linalg','scipy.optimize'))
        self.assertFalse(ok)
        self.assertEqual(error,'ImportError: libopenblas.so: cannot read file data')

    def test_doctor_fails_required_import_but_root_is_optional(self):
        def fake_import(module, additional=()):
            if module=='scipy':
                self.assertEqual(additional,('scipy.linalg','scipy.optimize'))
                return False,'ImportError: libopenblas.so: cannot read file data'
            return (False,'ModuleNotFoundError: ROOT') if module=='ROOT' else (True,'test-version')
        with patch.object(w,'import_check',side_effect=fake_import),redirect_stdout(io.StringIO()) as text:
            status=w.doctor(self.repo,self.settings)
        self.assertEqual(status,1)
        self.assertIn('IMPORT FAILED scipy',text.getvalue())
        self.assertIn('OPTIONAL IMPORT UNAVAILABLE ROOT',text.getvalue())
        with patch.object(w,'import_check',side_effect=lambda module,additional=():
                          (False,'ModuleNotFoundError: ROOT') if module=='ROOT' else (True,'test-version')):
            with redirect_stdout(io.StringIO()):status=w.doctor(self.repo,self.settings)
        self.assertEqual(status,0)

    def test_reopt_source_contract_and_explicit_config(self):
        workspace=self.init()
        run=workspace/'results/generated/runs/source-1s';rootfile=run/'artifacts/inputs/9x275_normalized.root'
        rootfile.parent.mkdir(parents=True);rootfile.write_bytes(b'synthetic path test, not physics ROOT')
        conf={'analysis':{'states':'1S','generator_branching_fraction':.0238,'q2_event_max':1.,'q2_flux_max':1.},
              'beam':{'high':{'electron_energy':9.,'proton_energy':275.,'visible_xsec_fb':210.136}}}
        (run/'run_manifest.json').write_text(json.dumps({'configuration':conf,'tasks':[{'task':'prepare_9x275','status':'DONE'}]}))
        receipt=run/'proton_selection.json';receipt.write_text(json.dumps({'proton_selection':{'mode':'required'}}))
        args=argparse.Namespace(states=['1S'],run_id='reopt-test',source_batch='source',energy='9x275',action='plan',observables=['pt2'],lumi=2.5,hash_inputs=False,resume=False)
        command=w.reopt_command(args,self.repo,workspace,{'selection':'proton-tagged'})
        self.assertIn('1S='+str(rootfile),command)
        self.assertEqual(args.generated_reopt_config['states']['1S']['run'],str(run))
        self.assertFalse('/afs/' in json.dumps(args.generated_reopt_config))
        receipt.write_text(json.dumps({'proton_selection':{'mode':'inclusive'}}))
        with self.assertRaisesRegex(ValueError,'selection differs'):w.reopt_command(args,self.repo,workspace,{'selection':'proton-tagged'})

if __name__=='__main__':unittest.main()
