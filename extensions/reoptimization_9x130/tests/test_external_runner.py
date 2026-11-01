"""Integration checks use synthetic ROOT data, never EOS production inputs."""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import uproot

PACKAGE = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ['YPSILON_TEST_FRAMEWORK']).resolve()
RUNNER = PACKAGE / 'run_reoptimization_9x130.py'


class ExternalRunner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name) / 'tagged'
        cls.root.mkdir()
        for folder in ('scripts/reoptimization', 'workflow'):
            (cls.root/folder).mkdir(parents=True, exist_ok=True)
        for p in (SOURCE/'scripts/reoptimization').glob('*.py'):
            shutil.copy2(p, cls.root/'scripts/reoptimization'/p.name)
        for name in ('proton_selection.py', 'proton_category_isolation.py'):
            shutil.copy2(SOURCE/'workflow'/name, cls.root/'workflow'/name)
        marker = dict(schema_version=1, mode='required', installation_status='COMPLETE', package_version=3,
                      source_root=str(Path(cls.temp.name)/'original'), category_root=str(cls.root))
        (cls.root/'.proton_category.json').write_text(json.dumps(marker))
        # Extract policy via the real helper in a separate process, no hardcoded metadata shortcuts.
        r = subprocess.run([sys.executable, '-B', '-c',
            'import sys,json;sys.path.insert(0,sys.argv[1]);from proton_selection import proton_selection_metadata;print(json.dumps(proton_selection_metadata()))',
            str(cls.root/'workflow')], capture_output=True,text=True,check=True)
        cls.policy=json.loads(r.stdout)
        cls.source_run=cls.root/'results/generated/runs/synthetic-1s'
        cls.input=cls.source_run/'artifacts/inputs/9x130_normalized.root'
        cls.input.parent.mkdir(parents=True)
        cls.cfg=json.loads((PACKAGE/'reoptimization_9x130.json').read_text())
        cls.cfg['synthetic']=True
        cls.config=Path(cls.temp.name)/'test_config.json'
        cls.config.write_text(json.dumps(cls.cfg))
        n=1000
        weight=1196.85/n
        arr=dict(weight=np.full(n,weight), true_mass=np.full(n,9.4603), true_t=np.linspace(.005,.995,n),
                 true_pT=np.sqrt(np.linspace(.005,.995,n)), true_W=np.full(n,40.), true_Q2=np.full(n,.1),
                 reco_pT=np.sqrt(np.linspace(.005,.995,n)),reco_mass=np.full(n,9.4603),
                 reco_W=np.full(n,40.),sigma_epz=np.full(n,18.),
                 reco_topology=np.zeros(n),reco_tagged_proton_p=np.full(n,-999.))
        meta=dict(normalization_profile='jul26-purestates-estarlight-939b11a',
                  normalization_scope='All generated events, before truth-state, Q2, W, and reconstruction cuts',
                  luminosity_invfb=10., entries=n,generated_events_h_EventCount=n,
                  cross_section_visible_ee_fb=119.685,generator_branching_fraction_ee=.0238,
                  target_sum_weight=1196.85,corrected_event_weight=weight,written_sum_weight=1196.85)
        with uproot.recreate(cls.input) as f:
            f['UpsilonTree']=arr
            f['framework_normalization']=json.dumps(meta)
            f['h_EventCount']=(np.array([float(n)]),np.array([0.,1.]))
        conf=dict(analysis=dict(states='1S',generator_branching_fraction=.0238,q2_event_max=1.,q2_flux_max=1.,projection_lumi_invfb=1.),
                  beam=dict(low=dict(electron_energy=9.,proton_energy=130.,w_min=20.,w_max=70.,visible_xsec_fb=119.685)))
        task=dict(task='prepare_9x130',status='DONE',command=['python3','prepare_latest_sample.py',
              '--cross-section-fb','119.685','--lumi-invfb','10.0','--generator-branching-fraction','0.0238','--out-root',str(cls.input)])
        cls.manifest=dict(configuration=conf,tasks=[task],run_status='PARTIAL')
        (cls.source_run/'run_manifest.json').write_text(json.dumps(cls.manifest))
        (cls.source_run/'proton_selection.json').write_text(json.dumps(dict(proton_selection=cls.policy)))
        cls.source_hash=cls.input.read_bytes()

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def invoke(self,command,*extra):
        return subprocess.run([sys.executable,'-B',str(RUNNER),command,'--framework-root',str(self.root),
                '--source-batch','synthetic','--states','1S','--observables','pt2','--config',str(self.config),*extra],
                capture_output=True,text=True,env={**os.environ,'MPLBACKEND':'Agg'})

    def test_01_plan_validate_partial_source(self):
        p=self.invoke('plan');self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        self.assertIn('9x130 target luminosity: 1.0',p.stdout)
        p=self.invoke('validate');self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        self.assertIn('PASS 1S rows 1000',p.stdout)
        self.assertIn('0.068801067',p.stdout)

    def test_02_source_energy_rejected(self):
        bad=copy.deepcopy(self.manifest);bad['configuration']['beam']['low']['proton_energy']=275.
        path=self.source_run/'run_manifest.json'
        try:
            path.write_text(json.dumps(bad))
            p=self.invoke('plan');self.assertNotEqual(p.returncode,0);self.assertIn('proton energy differs',p.stderr)
        finally:path.write_text(json.dumps(self.manifest))

    def test_03_missing_contract_rejected(self):
        path=self.source_run/'proton_selection.json';text=path.read_text()
        try:
            path.write_text('{}');p=self.invoke('plan')
            self.assertNotEqual(p.returncode,0);self.assertIn('proton selection contract',p.stderr)
        finally:path.write_text(text)

    def test_04_no_feasible_result_stays_unavailable(self):
        p=self.invoke('run','--run-id','synthetic-no-support','--hash-inputs')
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        run=self.root/'results/generated/runs/synthetic-no-support'
        out=run/'artifacts/reoptimization'
        self.assertIn('NO_FEASIBLE_SLOPE_BINNING',(out/'validation_status.csv').read_text())
        self.assertIn(r'9\times130',(out/'appendix_tables.tex').read_text())
        self.assertNotIn(r'9\times275',(out/'appendix_tables.tex').read_text())
        self.assertIn('Measurement scope:',(out/'summary.md').read_text())
        cfg=json.loads((run/'effective_config.json').read_text())
        self.assertEqual(cfg['target_lumi_invfb'],1.)
        self.assertEqual(cfg['proton_selection']['mode'],'required')
        self.assertEqual(self.input.read_bytes(),self.source_hash)
        p=self.invoke('run','--run-id','synthetic-no-support','--hash-inputs')
        self.assertNotEqual(p.returncode,0);self.assertIn('Run exists',p.stderr)
        p=self.invoke('run','--run-id','synthetic-no-support','--hash-inputs','--resume')
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        p=self.invoke('run','--run-id','synthetic-no-support','--hash-inputs','--resume','--lumi','2.5')
        self.assertNotEqual(p.returncode,0);self.assertIn('fingerprint differs',p.stderr)

    def test_05_unrelated_output_preserved(self):
        run=self.root/'results/generated/runs/unrelated';run.mkdir()
        (run/'keep.txt').write_text('keep')
        p=self.invoke('run','--run-id','unrelated')
        self.assertNotEqual(p.returncode,0);self.assertIn('nonempty',p.stderr)
        self.assertEqual((run/'keep.txt').read_text(),'keep')
        self.assertFalse((run/'run_manifest.json').exists())

if __name__=='__main__':unittest.main()
