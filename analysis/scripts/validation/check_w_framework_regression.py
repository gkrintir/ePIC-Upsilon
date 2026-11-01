#!/usr/bin/env python3
"""Run the installed extractor on historical ROOTs and compare note CSVs.

Central-value regression only: historical and current MC error conventions differ.
The bundled current snapshot is used only when no installed extractor is found.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PACKAGE = Path(__file__).resolve().parents[2]
CENTRAL = ['flux_int', 'n_true_weighted_visible', 'n_reco_weighted_visible',
           'eff_reco_over_true', 'n_unfolded_visible', 'sigma_gamma_p_pb',
           'sigma_gamma_p_truth_closure_pb', 'closure_unfold_over_truth']
KEYS = ['state', 'W_low', 'W_high']


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def compare(reference, rerun, out_prefix):
    old = pd.read_csv(reference)
    new = pd.read_csv(rerun)
    required = set(KEYS + ['flux_int', 'n_true_weighted_visible',
                          'n_reco_weighted_visible', 'eff_reco_over_true', 'sigma_gamma_p_pb'])
    for label, frame in [('reference', old), ('rerun', new)]:
        if required - set(frame.columns):
            raise ValueError(f'{label} lacks required columns: {sorted(required-set(frame.columns))}')
    for frame in (old, new):
        if frame.duplicated(KEYS).any():
            raise ValueError('Duplicate state/bin keys')
    if len(old) != len(new):
        raise ValueError('Different row counts')
    old = old.sort_values(KEYS).reset_index(drop=True)
    new = new.sort_values(KEYS).reset_index(drop=True)
    if not old.state.equals(new.state) or not np.allclose(old[KEYS[1:]], new[KEYS[1:]], atol=1e-10, rtol=0):
        raise ValueError('State or bin edges differ; no interpolation is allowed')
    rows, checks = [], []
    for column in CENTRAL:
        if column not in old or column not in new:
            checks.append({'column': column, 'status': 'NOT_IN_REFERENCE' if column not in old else 'FAIL_MISSING_COLUMN'})
            continue
        a, b = old[column].to_numpy(float), new[column].to_numpy(float)
        finite = np.isfinite(a) & np.isfinite(b)
        representation_only = np.zeros(len(a), dtype=bool)
        if column == 'eff_reco_over_true':
            # Archived extractor encoded a zero truth denominator as 0.
            # Current code correctly records an undefined ratio. The reco
            # numerator can be nonzero due to migration; it is checked above.
            empty = ((old.n_true_weighted_visible == 0) & (new.n_true_weighted_visible == 0)).to_numpy()
            representation_only = empty & (a == 0) & np.isnan(b)
        patterns = np.array_equal(np.isfinite(a)[~representation_only], np.isfinite(b)[~representation_only])
        ok = patterns and np.allclose(a[finite], b[finite], rtol=1e-8, atol=1e-10)
        diff = np.abs(b[finite] - a[finite])
        rel = diff / np.maximum(np.abs(a[finite]), 1e-20)
        checks.append({'column': column, 'status': 'PASS' if ok else 'FAIL',
                       'finite_pattern_equal_on_defined_bins': patterns,
                       'empty_efficiency_zero_to_nan_bins': int(representation_only.sum()),
                       'max_abs_difference': float(diff.max()) if len(diff) else 0.,
                       'max_relative_difference': float(rel.max()) if len(rel) else 0.})
        for i in range(len(a)):
            rows.append({'state': old.state.iloc[i], 'W_low': old.W_low.iloc[i], 'W_high': old.W_high.iloc[i],
                         'column': column, 'reference': a[i], 'rerun': b[i], 'difference': b[i]-a[i]})
    pd.DataFrame(rows).to_csv(str(out_prefix)+'_differences.csv', index=False)
    fig, (ax, ratio) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={'height_ratios':[3,1]})
    for state, color in zip(old.state.unique(), ['royalblue','crimson','forestgreen']):
        sel = (old.state == state).to_numpy() & np.isfinite(old.sigma_gamma_p_pb) & (old.sigma_gamma_p_pb > 0)
        if 'eff_reco_over_true' in old:
            sel &= old.eff_reco_over_true.to_numpy() > .05
        x = .5*(old.W_low.to_numpy()+old.W_high.to_numpy())
        ax.plot(x[sel], old.sigma_gamma_p_pb.to_numpy()[sel], 'o', markerfacecolor='none', color=color, label=f'{state} archived note data')
        ax.plot(x[sel], new.sigma_gamma_p_pb.to_numpy()[sel], 'x', color=color, label=f'{state} fresh framework')
        ratio.plot(x[sel], new.sigma_gamma_p_pb.to_numpy()[sel]/old.sigma_gamma_p_pb.to_numpy()[sel], 'o', color=color)
    ax.set_yscale('log'); ax.set_ylabel(r'$\sigma_{\gamma p}$ [pb]')
    ax.set_title('Same ROOT input and historical settings')
    ax.legend(fontsize=8); ratio.axhline(1, color='gray', lw=1)
    ratio.set_ylim(.999999, 1.000001); ratio.ticklabel_format(axis='y', style='plain', useOffset=False)
    ratio.set_ylabel('Fresh / archived'); ratio.set_xlabel(r'$W_{\gamma p}$ [GeV]')
    fig.tight_layout(); fig.savefig(str(out_prefix)+'.png', dpi=170); fig.savefig(str(out_prefix)+'.pdf'); plt.close(fig)
    return {'status': 'FAIL' if any(c['status'].startswith('FAIL') for c in checks) else 'PASS', 'checks': checks,
            'reference_sha256': sha256(reference), 'rerun_sha256': sha256(rerun), 'rows': len(old)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old-low', required=True)
    p.add_argument('--old-high', required=True)
    p.add_argument('--extractor')
    p.add_argument('--out-dir', required=True)
    args = p.parse_args()
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError('Choose a new output directory to avoid mixing regression runs')
    installed = PACKAGE/'scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py'
    snapshot = PACKAGE/'validation_snapshot/scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py'
    extractor = Path(args.extractor).resolve() if args.extractor else installed if installed.exists() else snapshot
    result = {'extractor': str(extractor), 'extractor_sha256': sha256(extractor),
              'using_bundled_snapshot': extractor == snapshot,
              'comparison_scope': 'Central values only; current same-MC uncertainties intentionally differ from archived fixed-efficiency uncertainties.',
              'cases': {}}
    cases = [
        ('old_high_family', args.old_high, 18, 275, 20, 140, '1S,2S,3S', 327.584, 'sigma_gp_family_18x275_fixed_unfolded'),
        ('old_low', args.old_low, 10, 100, 20, 50, '1S', 104.793, 'sigma_gp_1S_10x100_fixed'),
    ]
    for name, source, ee, ep, wlo, whi, states, xsec, refstem in cases:
        prefix = out/name
        cmd = [sys.executable, str(extractor), str(Path(source).resolve()), '--states', states,
               '--electron-energy', str(ee), '--proton-energy', str(ep), '--w-min', str(wlo), '--w-max', str(whi),
               '--n-w-bins', '14', '--q2-max', '0.001', '--lumi-invfb', '10', '--tuple-lumi-invfb', '10',
               '--visible-xsec-1s-fb', str(xsec), '--visible-xsec-2s-fb', '84.780', '--visible-xsec-3s-fb', '68.929',
               '--normalization-scope', 'truth_state', '--state-norm-mode', 'target_visible_xsec',
               '--weight-convention', 'visible_decay', '--stat-mode', 'mc', '--min-eff', '0.05', '--out-prefix', str(prefix)]
        record = {'command': cmd, 'input': str(Path(source).resolve()), 'input_sha256': sha256(source),
                  'event_q2_cut': None, 'flux_q2_max': .001}
        try:
            completed = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            (out/(name+'.log')).write_text(completed.stdout)
            if completed.returncode:
                raise RuntimeError(f'Extractor returned {completed.returncode}; see {name}.log')
            record.update(compare(PACKAGE/'reference_data'/(refstem+'.csv'), str(prefix)+'.csv', out/(name+'_regression')))
        except Exception as exc:
            record.update(status='FAIL', error=str(exc))
        result['cases'][name] = record
    result['status'] = 'PASS' if all(x['status']=='PASS' for x in result['cases'].values()) else 'FAIL'
    (out/'framework_regression.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status':result['status'], 'cases':{k:v['status'] for k,v in result['cases'].items()}}, indent=2))
    return 0 if result['status']=='PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
