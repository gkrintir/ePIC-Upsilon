"""Independent small-fixture checks of normalization audit; no EOS access."""
from copy import deepcopy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'audit_upsilon_normalization.py'
M_E = 0.00051099895
M_P = 0.9382720813


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_helper(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def physics_settings():
    return {'electron_energy': 9., 'proton_energy': 130., 'q2_cut': 1., 'q2_max': 1.,
            'w_edges': [20., 25., 30., 35., 40.], 'lumi_invfb': 10., 'br_ee': 0.0234,
            'mass': 9.46030, 'mass_sigma': .053, 'mass_window_sigma': 3.,
            'true_mass_window': .2, 'sigma_epz_center': 18., 'sigma_epz_halfwidth': 1.8,
            'min_eff': .05}


def fixture_arrays(sigma=119.685, luminosity=10., electron=9., proton=130.):
    """Two unequal-size unweighted jobs, each carrying sigma*L/N_job weights."""
    w = np.array([22., 23., 24., 26., 27., 28., 32., 33., 34., 36., 37., 38.])
    q2 = np.full(len(w), 0.001)
    proton_energy = np.sqrt(proton**2 + M_P**2)
    epz = (q2 + M_E**2) / (2 * electron)
    energy = electron - (w*w - M_P*M_P + q2 - 2*proton*epz) / (2*(proton_energy + proton))
    pz = epz - energy
    pt = np.sqrt(energy*energy - pz*pz - M_E*M_E)
    reco_mass = np.full(len(w), 9.46030)
    reco_mass[[1, 4, 9]] = -999.
    raw = np.r_[np.full(3, sigma*luminosity/3), np.full(9, sigma*luminosity/9)]
    return {
        'weight': raw,
        'weight_proj': raw * np.where(np.arange(len(w)) % 2, 6., 0.),
        'true_mass': np.full(len(w), 9.46030), 'true_W': w, 'true_Q2': q2,
        'true_scat_eMinus_pT': pt, 'true_scat_eMinus_eta': np.arcsinh(pz/pt),
        'true_scat_eMinus_E': energy,
        'reco_mass': reco_mass, 'reco_W': w.copy(), 'sigma_epz': np.full(len(w), 18.),
        'ele_id_mode': np.ones(len(w), dtype=np.int32),
    }


def write_root(path, arrays, generated=None, metadata=None):
    generated = len(arrays['weight']) if generated is None else generated
    with uproot.recreate(path) as root:
        root.mktree('UpsilonTree', {name: value.dtype for name, value in arrays.items()})
        root['UpsilonTree'].extend(arrays)
        root['h_EventCount'] = (np.array([float(generated)]), np.array([0., 1.]))
        if metadata is not None:
            root['framework_normalization'] = json.dumps(metadata)


class NormalizationAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.folder = Path(self.temporary.name)
        self.merged = self.folder / 'merged.root'
        self.normalized = self.folder / 'normalized.root'
        self.raw = fixture_arrays()
        write_root(self.merged, self.raw)
        self.normal = {key: value.copy() for key, value in self.raw.items()}
        ratios = self.raw['weight_proj'] / self.raw['weight']
        self.normal['weight'] = np.full(len(self.raw['weight']), 119.685*10/len(self.raw['weight']))
        self.normal['weight_proj'] = self.normal['weight'] * ratios
        self.metadata = {
            'input': str(self.merged), 'sha256': digest(self.merged), 'tree': 'UpsilonTree',
            'entries': len(self.raw['weight']), 'generated_events_h_EventCount': len(self.raw['weight']),
            'cross_section_visible_ee_fb': 119.685, 'luminosity_invfb': 10.,
            'target_sum_weight': 1196.85, 'written_sum_weight': 1196.85,
            'corrected_event_weight': 1196.85/len(self.raw['weight']),
            'normalization_source': 'Synthetic fixture rate, not production validation',
            'output': str(self.normalized),
        }
        write_root(self.normalized, self.normal, metadata=self.metadata)

    def tearDown(self):
        self.temporary.cleanup()

    def run_audit(self, extra=(), normalized=True, historical_branching=True):
        output = self.folder / 'audit_output'
        args = [sys.executable, str(SCRIPT), '--merged', str(self.merged),
                '--electron-energy', '9', '--proton-energy', '130',
                '--assigned-visible-xsec-fb', '119.685', '--lumi-invfb', '10',
                '--w-edges', '20,25,30,35,40', '--out-dir', str(output), '--hash-inputs']
        if historical_branching:
            args += ['--br-ee', '0.0234']
        if normalized:
            args += ['--normalized', str(self.normalized)]
        args += list(extra)
        before = {path: digest(path) for path in (self.merged, self.normalized) if path.exists()}
        result = subprocess.run(args, capture_output=True, text=True)
        after = {path: digest(path) for path in before}
        self.assertEqual(before, after, 'The read-only audit modified an input ROOT file')
        report_path = output / 'audit.json'
        self.assertTrue(report_path.is_file(), result.stdout + result.stderr)
        return result, json.loads(report_path.read_text())

    def check(self, report, identifier):
        items = [check for check in report['checks'] if check['id'] == identifier]
        self.assertEqual(len(items), 1, identifier + ': ' + json.dumps(report['checks']))
        return items[0]

    def write_csv(self, sigma_factor=1., truth_factor=1., yield_lumi=10.):
        physics = load_helper('normalization_physics').physics_summary(self.normal, physics_settings())
        rows = []
        for source in physics['bins']:
            rows.append({
                'state': '1S', 'W_low': source['W_low'], 'W_high': source['W_high'], 'W_center': source['W_center'],
                'flux_int': source['flux_int'],
                'n_true_weighted_visible': source['truth_sumw'] * truth_factor * yield_lumi/10,
                'n_reco_weighted_visible': source['reco_sumw'] * yield_lumi/10,
                'eff_reco_over_true': source['eff_reco_over_true'],
                'n_unfolded_visible': source['corrected_visible'] * yield_lumi/10,
                'sigma_gamma_p_pb': source['sigma_gamma_p_pb'] * sigma_factor,
                'sigma_gamma_p_truth_closure_pb': source['sigma_gamma_p_truth_closure_pb'] * sigma_factor,
                'closure_unfold_over_truth': source['closure_unfold_over_truth'], 'plotted': source['plotted'],
                'tuple_lumi_invfb': 10., 'yield_lumi_invfb': yield_lumi, 'projection_lumi_invfb': yield_lumi,
                'yield_luminosity_scale': yield_lumi/10., 'stat_mode': 'projected',
            })
        path = self.folder / '9x130_sigma_gamma_p_1s.csv'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_full_cli_preserves_inputs_and_accepts_correct_uniform_preparation(self):
        result, report = self.run_audit()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(report['status'], ('PASS', 'WARN'))
        for identifier in ('normalized_target_sum', 'normalized_uniform_weights',
                           'generated_count_merged', 'metadata_source_sha256',
                           'nonweight_branch_preservation'):
            self.assertEqual(self.check(report, identifier)['status'], 'PASS', identifier)
        self.assertEqual(report['schema_version'], 1)
        self.assertIsNotNone(report['files']['merged']['physics'])
        self.assertIsNotNone(report['files']['normalized']['physics'])
        self.assertFalse(report['physical_rate_provenance_verified'])
        self.assertAlmostEqual(report['files']['merged']['weights']['weight']['sumw'], 2*1196.85)
        self.assertAlmostEqual(report['files']['normalized']['weights']['weight']['sumw'], 1196.85)
        self.assertTrue(np.allclose([row['uniform_vs_global_shape_ratio'] for row in report['bin_comparison']],
                                   [.5, 1.5, 1.5, 1.5]))

    def test_jul26_default_and_historical_branching_override_change_only_inclusive_values(self):
        current_result, current = self.run_audit(historical_branching=False)
        legacy_result, legacy = self.run_audit()
        alias_result, alias = self.run_audit(('--branching-fraction', '0.0234'), historical_branching=False)
        for result in (current_result, legacy_result, alias_result):
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(current['settings']['br_ee'], .0238)
        self.assertEqual(legacy['settings']['br_ee'], .0234)
        self.assertEqual(alias['settings']['br_ee'], .0234)
        for current_bin, legacy_bin, alias_bin in zip(
                current['files']['normalized']['physics']['bins'],
                legacy['files']['normalized']['physics']['bins'],
                alias['files']['normalized']['physics']['bins']):
            for key in ('truth_sumw', 'reco_sumw', 'corrected_visible', 'flux_int'):
                self.assertEqual(current_bin[key], legacy_bin[key])
            self.assertAlmostEqual(current_bin['sigma_gamma_p_pb']/legacy_bin['sigma_gamma_p_pb'], .0234/.0238)
            self.assertEqual(alias_bin['sigma_gamma_p_pb'], legacy_bin['sigma_gamma_p_pb'])

    def test_recorded_generator_branching_fraction_detects_mixed_conventions(self):
        self.metadata['generator_branching_fraction_ee'] = .0238
        write_root(self.normalized, self.normal, metadata=self.metadata)
        result, report = self.run_audit(historical_branching=False)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.check(report, 'metadata_generator_branching_fraction_ee')['status'], 'PASS')
        result, report = self.run_audit()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'metadata_generator_branching_fraction_ee')['status'], 'FAIL')

    def test_corrected_target_rejects_prepared_root_with_old_high_rate(self):
        self.raw = fixture_arrays(sigma=327.584, proton=275.)
        write_root(self.merged, self.raw)
        self.normal = {key: value.copy() for key, value in self.raw.items()}
        self.normal['weight'] = np.full(12, 327.584*10/12)
        self.normal['weight_proj'] = self.normal['weight'] * (self.raw['weight_proj']/self.raw['weight'])
        self.metadata.update(sha256=digest(self.merged), cross_section_visible_ee_fb=327.584,
                             target_sum_weight=3275.84, written_sum_weight=3275.84,
                             corrected_event_weight=3275.84/12)
        write_root(self.normalized, self.normal, metadata=self.metadata)
        result, report = self.run_audit(('--proton-energy', '275', '--assigned-visible-xsec-fb', '210.136'),
                                        historical_branching=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'normalized_target_sum')['status'], 'FAIL')
        self.assertEqual(self.check(report, 'metadata_cross_section_visible_ee_fb')['status'], 'FAIL')

    def test_generated_denominator_mismatch_is_not_fixed_or_renormalized(self):
        write_root(self.merged, self.raw, generated=100)
        result, report = self.run_audit(normalized=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'generated_count_merged')['status'], 'FAIL')

    def test_same_target_total_but_nonuniform_prepared_weights_rejected(self):
        self.normal['weight'] = self.raw['weight'] / 2
        self.normal['weight_proj'] = self.normal['weight'] * (self.raw['weight_proj'] / self.raw['weight'])
        write_root(self.normalized, self.normal, metadata=self.metadata)
        result, report = self.run_audit()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'normalized_target_sum')['status'], 'PASS')
        self.assertEqual(self.check(report, 'normalized_uniform_weights')['status'], 'FAIL')

    def test_branching_luminosity_or_unit_scaled_preparation_fails_target_check(self):
        for factor in (1/0.0234, 10., 0.001):
            with self.subTest(factor=factor):
                changed = {key: value.copy() for key, value in self.normal.items()}
                changed['weight'] *= factor
                changed['weight_proj'] *= factor
                write_root(self.normalized, changed, metadata=self.metadata)
                result, report = self.run_audit()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.check(report, 'normalized_target_sum')['status'], 'FAIL')

    def test_nonweight_event_data_mutation_is_rejected(self):
        self.normal['true_W'][0] += 0.25
        write_root(self.normalized, self.normal, metadata=self.metadata)
        result, report = self.run_audit()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'nonweight_branch_preservation')['status'], 'FAIL')

    def test_nonweight_dtype_mutation_and_projective_factor_mutation_are_rejected(self):
        for mutation in ('dtype', 'projective'):
            with self.subTest(mutation=mutation):
                changed = {key: value.copy() for key, value in self.normal.items()}
                if mutation == 'dtype':
                    changed['ele_id_mode'] = changed['ele_id_mode'].astype(np.int64)
                    identifier = 'nonweight_branch_preservation'
                else:
                    changed['weight_proj'][1] *= 2
                    identifier = 'projective_factor_preservation'
                write_root(self.normalized, changed, metadata=self.metadata)
                result, report = self.run_audit()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.check(report, identifier)['status'], 'FAIL')

    def test_missing_truth_q2_blocks_physics_instead_of_ignoring_cut(self):
        del self.raw['true_Q2']
        write_root(self.merged, self.raw)
        result, report = self.run_audit(normalized=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'required_branches_merged')['status'], 'FAIL')
        self.assertIsNone(report['files']['merged']['physics'])

    def test_wrong_source_metadata_hash_rejected(self):
        self.metadata['sha256'] = '0' * 64
        write_root(self.normalized, self.normal, metadata=self.metadata)
        result, report = self.run_audit()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'metadata_source_sha256')['status'], 'FAIL')

    def test_output_collision_cannot_overwrite_supplied_reference_report(self):
        output = self.folder / 'collision_output'
        output.mkdir()
        reference = output / 'audit.json'
        contents = b'{"sentinel":"This supplied reference must remain unchanged."}\n'
        reference.write_bytes(contents)
        result = subprocess.run([sys.executable, str(SCRIPT), '--merged', str(self.merged),
            '--electron-energy', '9', '--proton-energy', '130',
            '--assigned-visible-xsec-fb', '119.685', '--w-edges', '20,25,30,35,40',
            '--reference-report', str(reference), '--out-dir', str(output)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(reference.read_bytes(), contents)

    def test_external_generator_value_is_compared_but_not_authenticated(self):
        result, report = self.run_audit(('--generator-visible-xsec-fb', '119.685',
                                       '--generator-source', 'Synthetic independent-rate fixture'))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.check(report, 'external_generator_rate')['status'], 'PASS')
        self.assertFalse(report['external_generator_rate']['provenance_authenticated'])
        self.assertFalse(report['physical_rate_provenance_verified'])
        result, report = self.run_audit(('--generator-visible-xsec-fb', '200',
                                       '--generator-source', 'Synthetic rate mismatch fixture'))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.check(report, 'external_generator_rate')['status'], 'FAIL')

    def test_cross_energy_report_decomposes_assigned_rate_population_and_flux(self):
        result, low = self.run_audit()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        reference = self.folder / 'low_reference.json'
        reference.write_text(json.dumps(low))
        self.merged = self.folder / 'high_merged.root'
        self.normalized = self.folder / 'high_normalized.root'
        self.raw = fixture_arrays(sigma=327.584, proton=275.)
        write_root(self.merged, self.raw)
        self.normal = {key: value.copy() for key, value in self.raw.items()}
        self.normal['weight'] = np.full(12, 327.584*10/12)
        self.normal['weight_proj'] = self.normal['weight'] * (self.raw['weight_proj']/self.raw['weight'])
        self.metadata.update(input=str(self.merged), sha256=digest(self.merged), output=str(self.normalized),
                             cross_section_visible_ee_fb=327.584, target_sum_weight=3275.84,
                             written_sum_weight=3275.84, corrected_event_weight=3275.84/12)
        write_root(self.normalized, self.normal, metadata=self.metadata)
        result, high = self.run_audit(('--proton-energy', '275', '--assigned-visible-xsec-fb', '327.584',
                                      '--reference-report', str(reference)))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.check(high, 'energy_comparison_settings')['status'], 'PASS')
        self.assertEqual(len(high['cross_energy']['rows']), 4)
        for row, low_bin, high_bin in zip(high['cross_energy']['rows'], low['bin_comparison'], high['bin_comparison']):
            expected = (327.584/119.685) * low_bin['flux_int']/high_bin['flux_int']
            self.assertAlmostEqual(row['assigned_visible_xsec_ratio'], 327.584/119.685)
            self.assertAlmostEqual(row['generated_truth_fraction_ratio'], 1.)
            self.assertAlmostEqual(row['actual_normalized_truth_sigma_ratio'], expected)
            self.assertAlmostEqual(row['predicted_uniform_truth_sigma_ratio'], expected)
            self.assertAlmostEqual(row['closure_actual_over_uniform_prediction'], 1.)

    def test_saved_cross_section_csv_matches_actual_root_estimator(self):
        path = self.write_csv()
        result, report = self.run_audit(('--csv', str(path)))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(report['saved_csv_check']['status'], 'PASS')

    def test_projected_40fb_csv_keeps_cross_section_and_invalid_luminosity_is_rejected(self):
        path = self.write_csv(yield_lumi=40.)
        result, report = self.run_audit(('--csv', str(path)))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(report['saved_csv_check']['status'], 'PASS')
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        rows[0]['yield_lumi_invfb'] = 'nan'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        result, report = self.run_audit(('--csv', str(path)))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(report['saved_csv_check']['status'], 'FAIL')

    def test_saved_csv_double_br_luminosity_units_and_stale_yields_rejected(self):
        for sigma_factor, truth_factor in ((1/.0234, 1.), (10., 1.), (.001, 1.), (1., 1.2)):
            with self.subTest(sigma_factor=sigma_factor, truth_factor=truth_factor):
                path = self.write_csv(sigma_factor, truth_factor)
                result, report = self.run_audit(('--csv', str(path)))
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(report['saved_csv_check']['status'], 'FAIL')

    def test_physics_exposes_uniform_vs_global_rescale_shape_difference(self):
        physics = load_helper('normalization_physics')
        settings = physics_settings()
        raw = physics.physics_summary(self.raw, settings)
        uniform = physics.physics_summary(self.normal, settings)
        global_arrays = {key: value.copy() for key, value in self.raw.items()}
        global_arrays['weight'] *= 0.5
        global_result = physics.physics_summary(global_arrays, settings)
        self.assertAlmostEqual(raw['cutflow']['all']['sumw'], 2*1196.85)
        self.assertAlmostEqual(uniform['cutflow']['all']['sumw'], 1196.85)
        self.assertAlmostEqual(global_result['cutflow']['all']['sumw'], 1196.85)
        self.assertEqual([r['truth_count'] for r in uniform['bins']], [3, 3, 3, 3])
        self.assertEqual([r['reco_count'] for r in uniform['bins']], [2, 2, 3, 2])
        ratios = [a['sigma_gamma_p_pb']/b['sigma_gamma_p_pb']
                  for a, b in zip(global_result['bins'], uniform['bins'])]
        self.assertTrue(np.allclose(ratios, [2., 2/3., 2/3., 2/3.]))
        for row in uniform['bins']:
            expected = (1196.85/4)/(1000*10*.0234*row['flux_int'])
            self.assertAlmostEqual(row['sigma_gamma_p_pb'], expected)
            self.assertAlmostEqual(row['closure_unfold_over_truth'], 1.)
            self.assertLess(abs(row['flux_relative_difference']), 1e-6)
        accumulated = physics.PhysicsAccumulator(settings)
        for first, last in ((0, 5), (5, 12)):
            accumulated.update({key: value[first:last] for key, value in self.normal.items()})
        streamed = accumulated.result()
        self.assertTrue(np.allclose([r['sigma_gamma_p_pb'] for r in streamed['bins']],
                                   [r['sigma_gamma_p_pb'] for r in uniform['bins']]))

    def test_pair_reader_checks_all_nonweight_fields_and_projective_factors(self):
        root_checks = load_helper('normalization_root_checks')
        first = root_checks.compare_root_pairs(self.merged, self.normalized, step_size=5)
        self.assertEqual(first['status'], 'PASS', first)
        self.assertEqual(first['projective_factor_invariance']['status'], 'PASS')
        self.assertAlmostEqual(first['weight_ratio']['values']['min'], .25)
        self.assertAlmostEqual(first['weight_ratio']['values']['max'], .75)
        self.assertTrue(all(row['changed_entries'] == 0 for row in first['branches'].values()))
        self.normal['weight_proj'][1] *= 2
        write_root(self.normalized, self.normal, metadata=self.metadata)
        changed = root_checks.compare_root_pairs(self.merged, self.normalized, step_size=5)
        self.assertEqual(changed['status'], 'CHANGED')
        self.assertEqual(changed['projective_factor_invariance']['changed_entries'], 1)
        self.assertEqual(changed['projective_factor_invariance']['first_changed_indices'], [1])


if __name__ == '__main__':
    unittest.main()
