"""Synthetic tests of tagged provenance; these do not validate CERN physics."""
import json
from pathlib import Path
import tempfile
import unittest
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import test_reconstruction
from settings import make_config as _make_config, make_reoptimization_config as _make_reoptimization_config, load_profile

# Explicit synthetic-test profile; no production default or local user settings.
_TEST_PROFILE = load_profile(Path(__file__).resolve().parents[2] / "configs/plots_9x275.example.json")
def make_config(*args, **kwargs):
    return _make_config(*args, profile=_TEST_PROFILE, **kwargs)
def make_reoptimization_config(*args, **kwargs):
    return _make_reoptimization_config(*args, profile=_TEST_PROFILE, **kwargs)
from upsilon_esr.common import sample_metadata
from upsilon_esr.proton_inputs import EXPECTED_POLICY, require_policy
from upsilon_esr.reoptimization_inputs import load_case, saved_slope
import test_metadata
from test_reoptimization_inputs import make_fixture


def patch_json(path, edit):
    value = json.loads(path.read_text())
    edit(value)
    path.write_text(json.dumps(value))


class ProtonInputsTests(unittest.TestCase):
    def test_settings_always_require_tag(self):
        self.assertTrue(make_config().require_proton_tag)
        self.assertTrue(make_reoptimization_config().require_proton_tag)
        self.assertIn('ProtonTagged', str(make_config().runs['1S']))

    def test_missing_inclusive_or_changed_denominator_rejected(self):
        for bad in (None, {}, dict(EXPECTED_POLICY, mode='inclusive'),
                    dict(EXPECTED_POLICY, truth_denominator='tagged events only'),
                    dict(EXPECTED_POLICY, predicate='reco_topology > 0')):
            with self.subTest(policy=bad), self.assertRaises(ValueError):
                require_policy(bad, 'fixture')

    def full_fixture(self, root):
        cfg = test_metadata.MetadataTests().build(root, 'PASS')
        cfg.require_proton_tag = True
        run = cfg.runs['1S']
        perf = run / 'artifacts/performance/9x275'
        policy = dict(EXPECTED_POLICY, label='Proton tagged')
        (run / 'proton_selection.json').write_text(json.dumps({'proton_selection': policy}))
        patch_json(perf / 'performance_audit.json', lambda v:v.update(proton_selection=policy))
        (perf / 'mass_fit.json').write_text(json.dumps({'proton_selection':policy}))
        return cfg, run, perf

    def test_full_run_requires_all_three_saved_contracts(self):
        with tempfile.TemporaryDirectory() as root:
            cfg, run, perf = self.full_fixture(root)
            result = sample_metadata(cfg, '1S')
            self.assertEqual(result['proton_selection_mode'], 'required')
            self.assertEqual(len(result['proton_contract_sha256']), 64)
            (run / 'proton_selection.json').unlink()
            with self.assertRaisesRegex(FileNotFoundError, 'proton_selection.json'):
                sample_metadata(cfg, '1S')
        for name in ('performance_audit.json','mass_fit.json'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as root:
                cfg, run, perf = self.full_fixture(root)
                patch_json(perf / name, lambda v:v['proton_selection'].update(mode='inclusive'))
                with self.assertRaisesRegex(ValueError, 'incompatible proton selection'):
                    sample_metadata(cfg, '1S')

    def test_guarded_reconstruction_renders_and_retains_policy(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            cfg = make_config(root, {s:s for s in ('1S','2S','3S')},
                              output_dir=root/'synthetic-output')
            cfg.synthetic = True
            example = test_reconstruction.CombinedMassFitTests()
            example.setUp()
            policy = dict(EXPECTED_POLICY, label='Proton tagged')
            for state, run in cfg.runs.items():
                perf = run/'artifacts/performance/9x275'
                perf.mkdir(parents=True)
                (run/'proton_selection.json').write_text(json.dumps({'proton_selection':policy}))
                audit = {'state':state, 'energy_label':'9x275 GeV',
                         'truth_selection':{'Q2_GeV2':[0,cfg.q2_max[state]],'W_GeV':[20,100]},
                         'reco_selection':{'proton_selection':policy},
                         'proton_selection':policy,
                         'statistical_configuration':{'active_luminosity_invfb':2.5}}
                (perf/'performance_audit.json').write_text(json.dumps(audit))
                fit = dict(example.calibrations[state], state=state,
                           proton_selection=policy, constant_selected_weights=True,
                           statistical_configuration={'active_luminosity_invfb':2.5},
                           selected_projected_yield=float(example.means[state].sum()))
                (perf/'mass_fit.json').write_text(json.dumps(fit))
                pd.DataFrame({'low':example.edges[:-1], 'high':example.edges[1:],
                              'projected_yield':example.means[state]}).to_csv(perf/'mass_bins.csv',index=False)
                rows=[]
                for var, edges in [('true_W',[20,40,60,100]),('true_pT',[0,.5,1,2])]:
                    for lo,hi in zip(edges[:-1],edges[1:]):
                        rows.append(dict(variable=var,low=lo,high=hi,denominator=100,
                                         selected=20,efficiency=.2,wilson_low=.16,wilson_high=.24))
                pd.DataFrame(rows).to_csv(perf/'selection_efficiency.csv',index=False)
            from upsilon_esr.reconstruction import plot_reconstruction
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fig, audit = plot_reconstruction(cfg)
            plt.close(fig)
            self.assertTrue((cfg.output_dir/'upsilon_es_reconstruction.pdf').is_file())
            for state in cfg.runs:
                self.assertEqual(audit['input_validation'][state]['proton_selection_mode'],'required')
                self.assertAlmostEqual(audit['yield_scale_factors'][state],1.)

    def reopt_fixture(self, root):
        cfg, _ = make_fixture(root)
        cfg.require_proton_tag = True
        run = cfg.reoptimization_run
        policy = dict(EXPECTED_POLICY, label='Proton tagged')
        patch_json(run / 'effective_config.json', lambda v:v.update(proton_selection=policy))
        for path in (run / 'artifacts/reoptimization').glob('*/*/result.json'):
            patch_json(path, lambda v:v.update(proton_selection=policy,
                selection_metadata={'proton_selection':policy}))
        return cfg

    def test_reoptimization_preserves_status_and_saved_slope(self):
        with tempfile.TemporaryDirectory() as root:
            cfg = self.reopt_fixture(root)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                frame, meta, _, _ = load_case(cfg, '1S', 'w')
            self.assertEqual(meta['overall_status'], 'NOT_VALIDATED')
            self.assertEqual(meta['proton_selection_mode'], 'required')
            self.assertEqual(set(frame.proton_selection_mode), {'required'})
            fit, predictions, profile = saved_slope(cfg, '1S', (0.,1.))
            self.assertEqual(fit['proton_selection']['mode'], 'required')
            self.assertEqual(set(profile.proton_selection_mode), {'required'})

    def test_reoptimization_rejects_untagged_or_mixed_input(self):
        for where in ('config', 'result', 'selection_metadata'):
            with self.subTest(where=where), tempfile.TemporaryDirectory() as root:
                cfg = self.reopt_fixture(root)
                run=cfg.reoptimization_run
                if where=='config':
                    patch_json(run/'effective_config.json',lambda v:v.pop('proton_selection'))
                else:
                    path=run/'artifacts/reoptimization/1S/w/result.json'
                    def edit(v):
                        target=v if where=='result' else v['selection_metadata']
                        target['proton_selection']['mode']='inclusive'
                    patch_json(path,edit)
                with self.assertRaises(ValueError):
                    load_case(cfg,'1S','w')


if __name__ == '__main__':
    unittest.main()
