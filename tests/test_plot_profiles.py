"""Profiles, no-fallback behavior and explicit diagnostic rendering policy."""
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'plotting'))
sys.path.insert(0, str(ROOT / 'plotting/tests'))
from settings import load_profile, make_config, make_reoptimization_config, new_output_dir, apply_display_rebin
from upsilon_esr.reoptimization_inputs import load_case
import test_proton_inputs


class PlotProfiles(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.values = json.loads((ROOT/'configs/plots_9x275.example.json').read_text())
        self.values.update(results_root='external-results', output_root='figures')
        self.path = self.directory/'profile.json'

    def profile(self, values=None):
        self.path.write_text(json.dumps(self.values if values is None else values))
        return load_profile(self.path)

    def test_paths_and_explicit_run_ids(self):
        p = self.profile()
        self.assertEqual(p['results_root'], self.directory/'external-results')
        cfg = make_config(profile=p)
        self.assertEqual(cfg.runs['1S'], p['results_root']/self.values['full_runs']['1S'])
        self.assertEqual(cfg.q2_max, self.values['q2_max'])
        self.assertTrue(cfg.require_proton_tag)
        reopt = make_reoptimization_config(profile=p)
        self.assertEqual(reopt.reoptimization_run, p['results_root']/self.values['reoptimization_run'])
        self.assertFalse(reopt.allow_diagnostic_results)
        self.assertEqual(p['results_root'].exists(), False)

    def test_environment_profile_and_no_search(self):
        p = self.profile()
        with mock.patch.dict(os.environ, {'UPSILON_PLOT_PROFILE': str(self.path)}):
            self.assertEqual(load_profile()['profile_path'], self.path)
        with mock.patch.dict(os.environ, {'UPSILON_PLOT_PROFILE': str(self.directory/'absent.json')}):
            with self.assertRaises(FileNotFoundError):
                load_profile()
        # An existing sibling folder cannot replace an absent explicitly selected run.
        (p['results_root']/'newest-run').mkdir(parents=True)
        cfg = make_reoptimization_config(profile=p)
        with self.assertRaises(FileNotFoundError):
            load_case(cfg, '1S', 'w')

    def test_invalid_profile_selections_rejected(self):
        for edit in [dict(proton_selection='inclusive'), dict(projection_lumi_invfb=0),
                     dict(projection_lumi_invfb=True), dict(allow_diagnostic_results='yes'),
                     dict(reoptimization_run='../other-run'), dict(full_runs={'1S':'one-run'})]:
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                self.profile(dict(self.values, **edit))
        bad = deepcopy(self.values); bad['categories']['1S'] = ''
        with self.assertRaises(ValueError): self.profile(bad)

    def test_output_directories_are_fresh(self):
        p = self.profile()
        first = new_output_dir(p, 'w'); (first/'preserve.txt').write_text('keep')
        second = new_output_dir(p, 'w')
        self.assertNotEqual(first, second)
        self.assertEqual((first/'preserve.txt').read_text(), 'keep')
        self.assertTrue(second.is_dir())

    def test_false_rebin_preserves_saved_bins(self):
        p = self.profile()
        cfg = make_reoptimization_config(profile=p)
        with mock.patch('upsilon_esr.reoptimization_inputs.load_case', side_effect=AssertionError('must not read data')):
            apply_display_rebin(cfg, p, 'w')
        self.assertEqual(cfg.rebin_edges, {})

    def test_merge_uses_only_saved_edges(self):
        p = self.profile(); p['display_rebin'] = {'1S':{'w':'keep_outer_bins'}}
        cfg = make_reoptimization_config(profile=p)
        frame = pd.DataFrame({'x_low':[20.,30.,40.,50.,60.], 'x_high':[30.,40.,50.,60.,70.]})
        with mock.patch('upsilon_esr.reoptimization_inputs.load_case', return_value=(frame, {}, [], {})):
            apply_display_rebin(cfg, p, 'w')
            self.assertEqual(cfg.rebin_edges, {'1S':{'w':[20.,30.,60.,70.]}})
            p['display_rebin'] = {'1S':{'w':[20.,30.,65.,70.]}}
            with self.assertRaisesRegex(ValueError, 'not a saved edge'):
                apply_display_rebin(cfg, p, 'w')

    def test_legacy_category_field_keeps_its_original_meaning(self):
        old = deepcopy(self.values)
        old.pop("saved_categories", None)
        self.assertEqual(self.profile(old)["categories"], old["categories"])
        old["saved_categories"] = dict(old["categories"], **{"2S": "different"})
        with self.assertRaisesRegex(ValueError, "disagree"):
            self.profile(old)

    def test_paper_profile_keeps_outer_bins_for_both_1s_observables(self):
        profile = load_profile(ROOT / 'configs/plots_9x275.paper.json')
        self.assertEqual(profile['display_rebin'], {'1S': {'w': 'keep_outer_bins', 'y': 'keep_outer_bins'}})
        self.assertFalse(profile['allow_diagnostic_results'])
        config = make_reoptimization_config(profile=profile)
        for section, variable, edges, paper in (
            ('w', 'w', [20, 40, 45, 50, 60, 100], [20, 40, 60, 100]),
            ('differential', 'y', [-1, .25, .5, .75, 1, 3], [-1, .25, 1, 3]),
        ):
            frame = pd.DataFrame({'x_low': edges[:-1], 'x_high': edges[1:]})
            with self.subTest(variable=variable), mock.patch(
                'upsilon_esr.reoptimization_inputs.load_case', return_value=(frame, {}, [], {})
            ):
                apply_display_rebin(config, profile, section)
                self.assertEqual(config.rebin_edges, {'1S': {variable: paper}})

    def test_diagnostics_require_explicit_permission(self):
        cfg = test_proton_inputs.ProtonInputsTests().reopt_fixture(self.directory/'fixture')
        cfg.allow_diagnostic_results = False
        with self.assertRaisesRegex(ValueError, 'diagnostic rendering is disabled'):
            load_case(cfg, '1S', 'w')
        cfg.allow_diagnostic_results = True
        with self.assertWarnsRegex(RuntimeWarning, 'NOT_VALIDATED'):
            frame, metadata, _, _ = load_case(cfg, '1S', 'w')
        self.assertEqual(metadata['overall_status'], 'NOT_VALIDATED')
        self.assertTrue(metadata['diagnostic_rendering_explicit'])
        self.assertFalse(metadata['publication_ready'])

    def test_wrong_luminosity_and_category_rejected(self):
        cfg = test_proton_inputs.ProtonInputsTests().reopt_fixture(self.directory/'fixture')
        cfg.lumi_invfb = 1.
        with self.assertRaisesRegex(ValueError, 'Saved inference is at'):
            load_case(cfg, '1S', 'w')
        cfg.lumi_invfb = 2.5
        cfg.categories['1S'] = 'unverified ID'
        with self.assertRaisesRegex(ValueError, 'category disagrees'):
            load_case(cfg, '1S', 'w')

    def test_notebooks_are_clean_and_use_profile(self):
        for p in (ROOT/'plotting').glob('*.ipynb'):
            nb = json.loads(p.read_text())
            code = '\n'.join(''.join(c.get('source',[])) for c in nb['cells'] if c['cell_type']=='code')
            with self.subTest(notebook=p.name):
                self.assertIn('PROFILE = load_profile()', code)
                self.assertNotIn('/eos/', code)
                self.assertNotIn('/afs/', code)
                self.assertNotIn('FRAMEWORK_ROOT = Path(', code)
                for c in nb['cells']:
                    if c['cell_type']=='code':
                        self.assertIsNone(c['execution_count'])
                        self.assertEqual(c['outputs'],[])
                        compile(''.join(c.get('source',[])),str(p),'exec')

if __name__=='__main__': unittest.main()
