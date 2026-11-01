"""State selection, card phase space and saved-artifact regression checks."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from test_normalization_audit import (SCRIPT, M_E, M_P, digest, fixture_arrays,
                                      load_helper, physics_settings, write_root)


STATES = {
    '1S': (9.46030, .053, .0238, 1., 119.685, 210.136),
    '2S': (10.02326, .056, .0191, .01, 29.297, 53.133),
    '3S': (10.35520, .059, .0218, .01, 23.447, 42.969),
}


class AllStateAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.folder = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def sample(self, state, proton=130.):
        mass, width, br, q2max, low_rate, high_rate = STATES[state]
        rate = low_rate if proton == 130 else high_rate
        folder = self.folder / (state + '_' + str(proton))
        folder.mkdir(exist_ok=True)
        merged, normalized = folder / 'merged.root', folder / 'normalized.root'
        arrays = fixture_arrays(sigma=rate, proton=proton)
        arrays['true_mass'][:] = mass
        arrays['reco_mass'][arrays['reco_mass'] > 0] = mass
        # Deliberate event above the excited-state Q2 boundary; keep kinematics consistent.
        arrays['true_Q2'][0] = .02
        q2, w = arrays['true_Q2'], arrays['true_W']
        epz = (q2 + M_E**2) / 18.
        energy = 9. - (w*w - M_P*M_P + q2 - 2*proton*epz) / (2*(np.sqrt(proton**2 + M_P**2) + proton))
        pz = epz - energy
        pt = np.sqrt(energy*energy - pz*pz - M_E*M_E)
        arrays['true_scat_eMinus_E'] = energy
        arrays['true_scat_eMinus_pT'] = pt
        arrays['true_scat_eMinus_eta'] = np.arcsinh(pz/pt)
        write_root(merged, arrays)
        normal = {key: value.copy() for key, value in arrays.items()}
        normal['weight'] = np.full(12, rate*10/12)
        normal['weight_proj'] = normal['weight'] * arrays['weight_proj']/arrays['weight']
        metadata = dict(input=str(merged), output=str(normalized), sha256=digest(merged),
                        tree='UpsilonTree', entries=12, generated_events_h_EventCount=12,
                        cross_section_visible_ee_fb=rate, luminosity_invfb=10.,
                        target_sum_weight=rate*10, written_sum_weight=rate*10,
                        corrected_event_weight=rate*10/12, generator_branching_fraction_ee=br,
                        normalization_source='Synthetic state fixture, not EOS validation')
        write_root(normalized, normal, metadata=metadata)
        return dict(state=state, proton=proton, rate=rate, merged=merged,
                    normalized=normalized, normal=normal, folder=folder)

    def run_audit(self, sample, extra=(), explicit_normalized=True):
        output = sample['folder'] / 'audit'
        args = [sys.executable, str(SCRIPT), '--state', sample['state'],
                '--merged', str(sample['merged']), '--electron-energy', '9',
                '--proton-energy', str(sample['proton']), '--assigned-visible-xsec-fb', str(sample['rate']),
                '--w-edges', '20,25,30,35,40', '--out-dir', str(output), '--hash-inputs']
        if explicit_normalized:
            args += ['--normalized', str(sample['normalized'])]
        result = subprocess.run(args + list(extra), capture_output=True, text=True)
        path = output / 'audit.json'
        self.assertTrue(path.exists(), result.stdout + result.stderr)
        return result, json.loads(path.read_text())

    @staticmethod
    def check(report, identifier):
        return next(row for row in report['checks'] if row['id'] == identifier)

    def write_csv(self, sample, path, wrong_state_rows=False):
        state = sample['state']
        mass, width, br, q2max, _, _ = STATES[state]
        settings = physics_settings()
        settings.update(state=state, proton_energy=sample['proton'], mass=mass, mass_sigma=width,
                        br_ee=br, q2_cut=q2max, q2_max=q2max)
        physics = load_helper('normalization_physics').physics_summary(sample['normal'], settings)
        mapping = {'n_true_weighted_visible': 'truth_sumw', 'n_reco_weighted_visible': 'reco_sumw',
                   'n_unfolded_visible': 'corrected_visible'}
        rows = []
        for source in physics['bins']:
            row = {key: source[key] for key in ('W_low', 'W_high', 'flux_int', 'eff_reco_over_true', 'sigma_gamma_p_pb')}
            row.update({target: source[origin] for target, origin in mapping.items()})
            row['state'] = state
            rows.append(row)
        if wrong_state_rows:
            rows = [dict(row, state='1S', sigma_gamma_p_pb=-1000) for row in rows] + rows
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def test_state_defaults_select_mass_branching_and_q2_phase_space(self):
        for state, values in STATES.items():
            with self.subTest(state=state):
                sample = self.sample(state)
                result, report = self.run_audit(sample)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                mass, width, br, q2max, rate, _ = values
                for name, expected in dict(state=state, mass=mass, mass_sigma=width, br_ee=br,
                                           q2_cut=q2max, q2_max=q2max).items():
                    self.assertEqual(report['settings'][name], expected)
                physics = report['files']['normalized']['physics']
                self.assertEqual(physics['cutflow']['truthmass']['count'], 12)
                self.assertEqual(physics['cutflow']['truthmassQ2W']['count'], 12 if state == '1S' else 11)
                self.assertAlmostEqual(physics['cutflow']['all']['sumw'], rate*10)
                self.assertEqual(self.check(report, 'independent_flux_integration')['status'], 'PASS')
                if state == '1S':
                    self.assertEqual(physics['conventions']['model_status'], 'AVAILABLE')
                    self.assertTrue(all(row['model_average_pb'] > 0 for row in physics['bins']))
                else:
                    self.assertEqual(self.check(report, 'real_photon_benchmark')['status'], 'NOT_APPLICABLE')
                    for row in report['bin_comparison']:
                        self.assertIsNone(row['model_average_pb'])
                        self.assertIsNone(row['normalized_sigma_over_model'])

    def test_run_directory_discovers_correct_state_csv_and_filters_mixed_rows(self):
        for state in ('2S', '3S'):
            with self.subTest(state=state):
                sample = self.sample(state)
                run_dir = sample['folder'] / 'run'
                inputs = run_dir / 'artifacts' / 'inputs'
                inputs.mkdir(parents=True)
                shutil.copyfile(sample['normalized'], inputs / '9x130_normalized.root')
                csv_path = run_dir / 'artifacts' / 'cross_sections' / ('9x130_sigma_gamma_p_' + state.lower() + '.csv')
                self.write_csv(sample, csv_path, wrong_state_rows=True)
                (csv_path.parent / '9x130_sigma_gamma_p_1s.csv').write_text('invalid stale 1S CSV\n')
                result, report = self.run_audit(sample, ('--run-dir', str(run_dir)), explicit_normalized=False)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(report['saved_csv_check']['path'], str(csv_path))
                self.assertEqual(report['saved_csv_check']['status'], 'PASS')

    def test_missing_state_csv_rows_fail_instead_of_using_legacy_1s(self):
        sample = self.sample('2S')
        path = sample['folder'] / 'no_state.csv'
        self.write_csv(sample, path)
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            del row['state']
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        result, report = self.run_audit(sample, ('--csv', str(path)))
        self.assertEqual(result.returncode, 1)
        self.assertEqual(report['saved_csv_check']['status'], 'FAIL')

    def test_wrong_state_input_fails_population_check(self):
        sample = self.sample('1S')
        result, report = self.run_audit(sample, ('--state', '2S'), explicit_normalized=False)
        self.assertEqual(result.returncode, 1)
        self.assertEqual(self.check(report, 'selected_state_population')['status'], 'FAIL')

    def test_explicit_q2_override_is_preserved(self):
        sample = self.sample('2S')
        result, report = self.run_audit(sample, ('--q2-cut', '.05', '--q2-max', '.05'))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(report['settings']['q2_cut'], .05)
        self.assertEqual(report['settings']['q2_max'], .05)
        self.assertEqual(report['files']['normalized']['physics']['cutflow']['truthmassQ2W']['count'], 12)

    def test_excited_state_cross_energy_comparison_leaves_model_ratio_null(self):
        low = self.sample('2S')
        result, low_report = self.run_audit(low)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        high = self.sample('2S', proton=275.)
        reference = low['folder'] / 'audit' / 'audit.json'
        result, report = self.run_audit(high, ('--reference-report', str(reference)))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.check(report, 'energy_comparison_settings')['status'], 'PASS')
        for row in report['cross_energy']['rows']:
            self.assertIsNone(row['real_photon_model_ratio'])
            self.assertAlmostEqual(row['closure_actual_over_uniform_prediction'], 1.)
        low_report['settings']['state'] = '1S'
        reference.write_text(json.dumps(low_report))
        result, report = self.run_audit(high, ('--reference-report', str(reference)))
        self.assertEqual(result.returncode, 1)
        self.assertEqual(self.check(report, 'energy_comparison_settings')['status'], 'FAIL')


if __name__ == '__main__':
    unittest.main()
