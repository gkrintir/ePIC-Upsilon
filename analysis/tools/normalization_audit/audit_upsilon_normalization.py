#!/usr/bin/env python3
"""Read-only state-specific Upsilon normalization audit. Python 3.9+, numpy, uproot.

This program never writes a ROOT file. Assigned cross sections are hypotheses,
not measurements recovered from the tuple. See README.md for assumptions.
"""
import argparse
import csv
import datetime
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
import uproot

from normalization_physics import PhysicsAccumulator, required_branches, validate_settings
from normalization_root_checks import inspect_root, compare_root_pairs, beam_formula_check


STATE_DEFAULTS = {
    '1S': dict(mass=9.46030, mass_sigma=.053, br_ee=.0238, q2_cut=1., q2_max=1.),
    '2S': dict(mass=10.02326, mass_sigma=.056, br_ee=.0191, q2_cut=.01, q2_max=.01),
    '3S': dict(mass=10.35520, mass_sigma=.059, br_ee=.0218, q2_cut=.01, q2_max=.01),
}


def resolve_state_defaults(args):
    """Use each July card's phase space unless the user explicitly overrides it."""
    for name, value in STATE_DEFAULTS[args.state].items():
        if getattr(args, name) is None:
            setattr(args, name, value)
    return args


def finite(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def ratio(a, b):
    return a / b if a is not None and b is not None and b != 0 else None


def close(a, b, rtol=1e-7, atol=1e-10):
    a, b = finite(a), finite(b)
    return (a is None and b is None) or (a is not None and b is not None and abs(a-b) <= atol+rtol*abs(b))


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return finite(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def check(report, key, status, detail, evidence=None):
    row = dict(id=key, status=status, detail=detail)
    if evidence is not None:
        row['evidence'] = clean(evidence)
    report['checks'].append(row)


def numeric_check(report, key, actual, expected, detail, rtol=1e-7, atol=1e-10):
    passed = finite(actual) is not None and finite(expected) is not None and close(actual, expected, rtol, atol)
    check(report, key, 'PASS' if passed else 'FAIL', detail, dict(actual=actual, expected=expected))


def csv_write(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=keys)
        writer.writeheader()
        writer.writerows(clean(rows))


class WeightSummary:
    def __init__(self, group_limit=5000):
        self.rows = self.finite_count = self.negative = self.zero = 0
        self.sumw = self.sumw2 = 0.0
        self.minimum, self.maximum = math.inf, -math.inf
        self.groups, self.group_limit, self.truncated = {}, group_limit, False
        self.outside_target = 0

    def update(self, values, target=None):
        values = np.asarray(values, dtype=float)
        self.rows += len(values)
        valid = values[np.isfinite(values)]
        self.finite_count += len(valid)
        self.negative += int(np.sum(valid < 0))
        self.zero += int(np.sum(valid == 0))
        self.sumw += float(np.sum(valid, dtype=np.float64))
        self.sumw2 += float(np.sum(valid*valid, dtype=np.float64))
        if len(valid):
            self.minimum = min(self.minimum, float(np.min(valid)))
            self.maximum = max(self.maximum, float(np.max(valid)))
        if target is not None:
            self.outside_target += int(np.sum(~np.isclose(values, target, rtol=1e-7, atol=1e-12)))
        unique, counts = np.unique(valid, return_counts=True)
        for v, n in zip(unique, counts):
            v = float(v)
            if v in self.groups or len(self.groups) < self.group_limit:
                self.groups[v] = self.groups.get(v, 0)+int(n)
            else:
                self.truncated = True

    def result(self):
        return clean(dict(rows=self.rows, finite_count=self.finite_count,
                          nonfinite_count=self.rows-self.finite_count, negative_count=self.negative,
                          zero_count=self.zero, sumw=self.sumw, sumw2=self.sumw2,
                          minimum=self.minimum, maximum=self.maximum,
                          mean=ratio(self.sumw, self.finite_count),
                          effective_entries=ratio(self.sumw*self.sumw, self.sumw2) if not self.negative else None,
                          unique_groups_retained=len(self.groups), groups_truncated=self.truncated,
                          outside_uniform_target_count=self.outside_target,
                          groups=[dict(weight=w, count=n, sumw=w*n) for w, n in sorted(self.groups.items())]))


def scan(path, metadata, settings, args):
    names = metadata.get('branches', {})
    required = required_branches(settings)
    missing = sorted(set(required)-set(names))
    physics = None if missing else PhysicsAccumulator(settings)
    weights = {name: WeightSummary() for name in ('weight', 'weight_proj') if name in names}
    target = args.assigned_visible_xsec_fb*args.lumi_invfb/max(metadata.get('entries', 0), 1)
    wanted = sorted(set(weights) | (set(required) if not missing else set()))
    with uproot.open(path) as root:
        if wanted:
            for arrays in root[args.tree].iterate(wanted, step_size=args.chunk_size, library='np', how=dict):
                for name, stats in weights.items():
                    stats.update(arrays[name], target if name == 'weight' else None)
                if physics is not None:
                    physics.update(arrays)
    return dict(metadata=metadata, weights={k:v.result() for k,v in weights.items()},
                missing_required_branches=missing, physics=physics.result() if physics else None)


def file_checks(report, label, data, args):
    meta, weights = data['metadata'], data['weights']
    check(report, 'nonempty_tree_'+label, 'PASS' if meta.get('entries',0)>0 else 'FAIL',
          'A positive generated-event population is required for normalization.',meta.get('entries'))
    check(report, 'required_branches_'+label, 'FAIL' if data['missing_required_branches'] else 'PASS',
          'All W-extraction branches, including true_Q2 and the declared proton category branches, must exist.', data['missing_required_branches'])
    hist = meta.get('histograms', {}).get('h_EventCount')
    if hist:
        numeric_check(report, 'generated_count_'+label, hist.get('sum_without_flow'), meta['entries'],
                      'h_EventCount regular-bin sum must equal all generated tuple rows.', rtol=0, atol=1e-6)
        check(report, 'generated_count_flow_'+label, 'WARN' if abs(hist.get('flow_sum') or 0)>1e-8 else 'PASS',
              'Nonzero under/overflow in h_EventCount makes its denominator convention ambiguous.', hist)
    else:
        check(report, 'generated_count_'+label, 'NOT_CHECKED', 'No h_EventCount; tuple completeness cannot be checked.')
    for name, stats in weights.items():
        bad = stats['nonfinite_count']+stats['negative_count']
        check(report, 'valid_'+name+'_'+label, 'FAIL' if bad else 'PASS',
              'Nonfinite/negative weights are incompatible with this positive-weight preparation/extraction.',
              {k:stats[k] for k in ('nonfinite_count','negative_count','zero_count')})
    if 'weight' in weights:
        s = weights['weight']
        s['literal_sumw_over_lumi_fb'] = s['sumw']/args.lumi_invfb
        s['sumw_over_assigned_target'] = s['sumw']/(args.assigned_visible_xsec_fb*args.lumi_invfb)
        if label == 'merged':
            check(report, 'merged_total_vs_assigned', 'PASS' if close(s['sumw_over_assigned_target'],1) else 'WARN',
                  'Literal raw sum(weight)/(assigned sigma*L). A mismatch is diagnostic, not proof that either rate is correct.',
                  s['sumw_over_assigned_target'])
            if s['zero_count']:
                check(report, 'merged_zero_weights', 'WARN', 'Current prepare_latest_sample requires strictly positive input weight.')
            if meta.get('framework_normalization'):
                check(report, 'merged_is_original', 'WARN', 'The supplied merged file already contains framework_normalization metadata.')
        else:
            numeric_check(report, 'normalized_target_sum', s['sumw'], args.assigned_visible_xsec_fb*args.lumi_invfb,
                          'Prepared weights must sum to the explicitly assigned visible sigma times tuple luminosity.')
            check(report, 'normalized_uniform_weights', 'PASS' if s['outside_uniform_target_count']==0 else 'FAIL',
                  'Every prepared weight must match assigned sigma*L/N_rows under the current unweighted-generator assumption.',
                  dict(expected=args.assigned_visible_xsec_fb*args.lumi_invfb/max(meta['entries'],1),
                       different_events=s['outside_uniform_target_count']))


def metadata_checks(report, args):
    normalized = report['files'].get('normalized')
    if not normalized:
        check(report, 'normalized_comparison', 'NOT_CHECKED', 'No normalized ROOT file supplied or discovered.')
        return
    meta = normalized['metadata'].get('framework_normalization')
    if not isinstance(meta, dict):
        check(report, 'normalization_metadata', 'NOT_CHECKED', 'No readable framework_normalization JSON in normalized ROOT.')
        return
    raw = report['files']['merged']
    pairs = [('cross_section_visible_ee_fb', args.assigned_visible_xsec_fb), ('luminosity_invfb',args.lumi_invfb),
             ('entries',raw['metadata']['entries']),('generated_events_h_EventCount',raw['metadata']['entries']),
             ('target_sum_weight',args.assigned_visible_xsec_fb*args.lumi_invfb),
             ('corrected_event_weight',args.assigned_visible_xsec_fb*args.lumi_invfb/max(raw['metadata']['entries'],1))]
    for field, expected in pairs:
        if field in meta:
            numeric_check(report, 'metadata_'+field, meta[field], expected, 'Embedded normalization metadata versus supplied inputs/settings.',
                          rtol=0 if field in ('entries','generated_events_h_EventCount') else 1e-7,
                          atol=1e-6 if field in ('entries','generated_events_h_EventCount') else 1e-10)
        else:
            check(report, 'metadata_'+field, 'NOT_CHECKED', 'Embedded metadata field absent.')
    if meta.get('generator_branching_fraction_ee') is not None:
        numeric_check(report, 'metadata_generator_branching_fraction_ee',
                      meta['generator_branching_fraction_ee'], args.br_ee,
                      'Prepared input generator branching convention versus inclusive audit conversion.',
                      rtol=1e-10, atol=1e-12)
    observed = raw['metadata'].get('sha256')
    if args.hash_inputs and observed and meta.get('sha256'):
        check(report, 'metadata_source_sha256', 'PASS' if observed==meta['sha256'] else 'FAIL',
              'Prepared-file recorded source hash must match the original merged file bytes.',
              dict(recorded=meta['sha256'], actual=observed))
    else:
        check(report, 'metadata_source_sha256', 'NOT_CHECKED', 'Use --hash-inputs; a recorded source sha256 is also required.')
    check(report, 'metadata_source_path', 'PASS' if os.path.realpath(str(meta.get('input','')))==os.path.realpath(args.merged) else 'WARN',
          'Source paths can differ after moving a file; hashes are the stronger identity check.',
          dict(recorded=meta.get('input'), supplied=args.merged))
    if args.normalization_audit:
        external = json.loads(Path(args.normalization_audit).read_text())
        differences = {k:dict(embedded=meta.get(k), external=external.get(k)) for k in set(meta)|set(external)
                       if meta.get(k)!=external.get(k)}
        check(report, 'external_audit_vs_embedded', 'FAIL' if differences else 'PASS',
              'External preparation JSON must agree with the normalized ROOT metadata.', differences)


def build_bins(report, settings, args):
    merged = report['files']['merged']
    if not merged.get('physics'):
        return []
    normalized = report['files'].get('normalized', {}).get('physics')
    rawsum = merged['weights']['weight']['sumw']
    entries = merged['metadata']['entries']
    per_event = args.assigned_visible_xsec_fb*args.lumi_invfb/max(entries,1)
    global_scale = ratio(args.assigned_visible_xsec_fb*args.lumi_invfb, rawsum)
    result = []
    for index, raw in enumerate(merged['physics']['bins']):
        norm = normalized['bins'][index] if normalized else {}
        denom = 1000*args.lumi_invfb*args.br_ee*raw['flux_int']
        uniform_truth = raw['truth_count']*per_event
        uniform_reco = raw['reco_count']*per_event
        row = {key:raw[key] for key in ('W_low','W_high','W_center','flux_int','independent_flux','model_average_pb')}
        row.update(truth_count=raw['truth_count'], reco_count=raw['reco_count'],
                   upper_edge_inclusive=raw.get('upper_edge_inclusive'),
                   truth_count_exact_upper_edge=raw.get('truth_count_exact_upper_edge'),
                   reco_count_exact_upper_edge=raw.get('reco_count_exact_upper_edge'),
                   truth_fraction_all_generated=ratio(raw['truth_count'],entries),
                   raw_truth_sumw=raw['truth_sumw'],raw_reco_sumw=raw['reco_sumw'],
                   raw_sigma_pb=raw['sigma_gamma_p_pb'], raw_truth_sigma_pb=ratio(raw['truth_sumw'],denom),
                   normalized_truth_sumw=norm.get('truth_sumw'), normalized_reco_sumw=norm.get('reco_sumw'),
                   normalized_sigma_pb=norm.get('sigma_gamma_p_pb'),
                   normalized_truth_sigma_pb=ratio(norm.get('truth_sumw'),denom),
                   normalized_eff_reco_over_true=norm.get('eff_reco_over_true'),
                   uniform_target_truth_sumw=uniform_truth, uniform_target_reco_sumw=uniform_reco,
                   uniform_target_truth_sigma_pb=ratio(uniform_truth,denom),
                   global_rescale_truth_sigma_pb=ratio(raw['truth_sumw']*global_scale,denom) if global_scale is not None else None,
                   uniform_vs_global_shape_ratio=ratio(uniform_truth,raw['truth_sumw']*global_scale) if global_scale else None,
                   assigned_visible_xsec_fb=args.assigned_visible_xsec_fb,tuple_lumi_invfb=args.lumi_invfb,br_ee=args.br_ee,
                   correction_reco_over_truth=norm.get('eff_reco_over_true',raw['eff_reco_over_true']),
                   normalized_sigma_over_model=ratio(norm.get('sigma_gamma_p_pb'),raw['model_average_pb']),
                   uniform_sigma_over_model=ratio(ratio(uniform_truth,denom),raw['model_average_pb']),
                   normalized_plotted=norm.get('plotted'),
                   flux_relative_difference=ratio(raw['flux_int']-raw['independent_flux'],raw['independent_flux']))
        result.append(row)
    worst = max((abs(r['flux_relative_difference']) for r in result if r['flux_relative_difference'] is not None),default=0.)
    check(report,'independent_flux_integration','PASS' if worst<=1e-3 else 'WARN' if worst<.01 else 'FAIL',
          '500-point framework W integration versus independent 2D Gauss integration in log(y), log(Q2); this tests the EPA implementation, not production-generator provenance.',
          dict(max_relative_difference=worst))
    shape = max((abs(r['uniform_vs_global_shape_ratio']-1) for r in result if r['uniform_vs_global_shape_ratio'] is not None),default=0.)
    check(report,'weight_replacement_shape','WARN' if shape>1e-5 else 'PASS',
          'Uniform replacement versus a global rescaling of original weights: equal totals need not give equal W shapes.',
          dict(max_bin_relative_shape_change=shape))
    if normalized:
        unsupported = [i for i,b in enumerate(normalized['bins']) if b['truth_sumw']>0 and b['reco_sumw']<=0]
        check(report,'unsupported_reco_bins','WARN' if unsupported else 'PASS',
              'Positive truth with no reconstructed support: central correction is undefined; truth diagnostic remains available.',unsupported)
    return result


def csv_check(report, args):
    if not args.csv:
        check(report,'saved_csv','NOT_CHECKED','No W cross-section CSV supplied or discovered.')
        return None
    data = report['files'].get('normalized')
    if not data or not data.get('physics'):
        check(report,'saved_csv','NOT_CHECKED','Supply a valid normalized ROOT to compare its saved full-sample bin-by-bin CSV.')
        return None
    with open(args.csv, newline='') as handle:
        source_rows = [r for r in csv.DictReader(handle)
                       if r.get('state', '1S') == args.state]
    indexed = {}
    for r in source_rows:
        key = (float(r['W_low']),float(r['W_high']))
        if key in indexed:
            raise ValueError('Duplicate W bin in saved CSV: '+str(key))
        indexed[key] = r
    output, failures, missing = [], [], []
    mapping = {'n_true_weighted_visible':'truth_sumw','n_reco_weighted_visible':'reco_sumw',
               'eff_reco_over_true':'eff_reco_over_true','n_unfolded_visible':'corrected_visible',
               'flux_int':'flux_int','sigma_gamma_p_pb':'sigma_gamma_p_pb'}
    for b in data['physics']['bins']:
        key=(b['W_low'],b['W_high'])
        r=indexed.pop(key,None)
        if r is None:
            missing.append(key)
            continue
        # Count columns are at yield_lumi; cross-section columns divide by that same L.
        yield_lumi=finite(r.get('yield_lumi_invfb'))
        if 'yield_lumi_invfb' in r and (yield_lumi is None or yield_lumi<=0):
            failures.append(dict(W_low=b['W_low'],column='yield_lumi_invfb',saved=r['yield_lumi_invfb'],reason='invalid luminosity metadata'))
        if yield_lumi is None:
            yield_lumi=args.csv_yield_lumi_invfb or args.lumi_invfb
        scale=yield_lumi/args.lumi_invfb
        if yield_lumi<=0:
            failures.append(dict(W_low=b['W_low'],column='yield_lumi_invfb',saved=yield_lumi))
        optional_mapping={k:k for k in ('sigma_gamma_p_truth_closure_pb','closure_unfold_over_truth') if k in r}
        for column, field in dict(mapping,**optional_mapping).items():
            if column not in r:
                missing.append(column)
                continue
            expected=b[field]
            if column.startswith('n_') and expected is not None:
                expected*=scale
            actual=finite(r[column])
            passed=close(actual,expected)
            item=dict(W_low=b['W_low'],W_high=b['W_high'],column=column,saved=actual,recomputed=expected,
                      saved_over_recomputed=ratio(actual,expected),status='PASS' if passed else 'FAIL')
            output.append(item)
            if not passed:
                failures.append(item)
        if 'plotted' in r:
            wanted=str(r['plotted']).lower() in ('1','true')
            if wanted != b['plotted']:
                failures.append(dict(W_low=b['W_low'],column='plotted',saved=wanted,recomputed=b['plotted']))
        for key, expected in [('tuple_lumi_invfb',args.lumi_invfb)]:
            if key in r and (finite(r[key]) is None or not close(r[key],expected)):
                failures.append(dict(W_low=b['W_low'],column=key,saved=r[key],recomputed=expected))
        if 'yield_luminosity_scale' in r and (finite(r['yield_luminosity_scale']) is None or not close(r['yield_luminosity_scale'],scale)):
            failures.append(dict(W_low=b['W_low'],column='yield_luminosity_scale',saved=r['yield_luminosity_scale'],recomputed=scale))
        if r.get('stat_mode')=='projected' and 'projection_lumi_invfb' in r:
            if finite(r['projection_lumi_invfb']) is None or not close(r['projection_lumi_invfb'],yield_lumi):
                failures.append(dict(W_low=b['W_low'],column='projection_lumi_invfb',saved=r['projection_lumi_invfb'],recomputed=yield_lumi))
    if indexed:
        missing.extend(['extra bin '+str(k) for k in indexed])
    check(report,'saved_csv','FAIL' if failures or missing else 'PASS',
          'Recompute full-sample '+args.state+' bin-by-bin central values from actual normalized weights, with recorded yield luminosity.',
          dict(failed_fields=len(failures),missing=missing))
    check(report,'saved_csv_provenance','NOT_CHECKED',
          'Numerical agreement is checked. A CSV without an authenticated ROOT hash does not prove source identity.')
    sigma_ratios=[r['saved_over_recomputed'] for r in output if r['column']=='sigma_gamma_p_pb' and r['saved_over_recomputed'] is not None]
    hints=[]
    if sigma_ratios:
        med=float(np.median(sigma_ratios))
        if np.allclose(sigma_ratios,med,rtol=.01,atol=1e-10) and not close(med,1,rtol=.01):
            candidates={'factor 2 / half-sample luminosity':2., 'factor 1/2 / half-sample luminosity':.5,
                        'fb-to-pb factor':1000., 'inverse fb-to-pb factor':.001,
                        'missing division by B_ee':args.br_ee,'extra division by B_ee':1/args.br_ee}
            hints=[name for name,val in candidates.items() if close(med,val,rtol=.03)]
            hints.append('Approximately constant saved/recomputed factor: %.9g; inspect assigned rate, L, B and unit conventions.'%med)
    return dict(status='FAIL' if failures or missing else 'PASS',path=args.csv,rows=output,failures=failures,missing=missing,factor_hints=hints)


def energy_comparison(report,args):
    if not args.reference_report:
        return None
    reference=json.loads(Path(args.reference_report).read_text())
    reference_settings=reference['settings']
    fields=('q2_cut','q2_max','mass','mass_sigma','mass_window_sigma','true_mass_window','proton_selection')
    mismatched={k:[reference_settings.get(k),report['settings'].get(k)] for k in fields
                if reference_settings.get(k)!=report['settings'].get(k)}
    if reference_settings.get('state', '1S') != report['settings'].get('state', '1S'):
        mismatched['state'] = [reference_settings.get('state', '1S'), report['settings'].get('state', '1S')]
    check(report,'energy_comparison_settings','FAIL' if mismatched else 'PASS',
          'Common-bin energy comparison requires the same state, mass, Q2 and proton-selection definitions.',mismatched)
    reco_differences={k:[reference_settings.get(k),report['settings'].get(k)]
                      for k in ('sigma_epz_center','sigma_epz_halfwidth','min_eff')
                      if reference_settings.get(k)!=report['settings'].get(k)}
    check(report,'energy_reco_selection_settings','WARN' if reco_differences else 'PASS',
          'Reco exclusivity and plotting cuts are recorded separately from the truth-rate decomposition.',reco_differences)
    refs={(b['W_low'],b['W_high']):b for b in reference.get('bin_comparison',[])}
    output=[]
    edge_differences=[]
    for b in report['bin_comparison']:
        r=refs.get((b['W_low'],b['W_high']))
        if r is None:
            continue
        if b.get('upper_edge_inclusive')!=r.get('upper_edge_inclusive'):
            edge_differences.append(dict(W_low=b['W_low'],W_high=b['W_high'],current_inclusive=b.get('upper_edge_inclusive'),
                                         reference_inclusive=r.get('upper_edge_inclusive'),
                                         current_truth_at_edge=b.get('truth_count_exact_upper_edge'),
                                         reference_truth_at_edge=r.get('truth_count_exact_upper_edge')))
        cross=ratio(b['assigned_visible_xsec_fb'],r['assigned_visible_xsec_fb'])
        fraction=ratio(b['truth_fraction_all_generated'],r['truth_fraction_all_generated'])
        flux=ratio(r['flux_int'],b['flux_int'])
        br=ratio(r['br_ee'],b['br_ee'])
        expected=cross*fraction*flux*br if all(v is not None for v in (cross,fraction,flux,br)) else None
        actual=ratio(b['normalized_truth_sigma_pb'],r['normalized_truth_sigma_pb'])
        output.append(dict(W_low=b['W_low'],W_high=b['W_high'],ratio_direction='current/reference',
                           assigned_visible_xsec_ratio=cross,generated_truth_fraction_ratio=fraction,
                           inverse_flux_ratio=flux,inverse_branching_ratio=br,
                           predicted_uniform_truth_sigma_ratio=expected,actual_normalized_truth_sigma_ratio=actual,
                           actual_normalized_reco_corrected_sigma_ratio=ratio(b['normalized_sigma_pb'],r['normalized_sigma_pb']),
                           closure_actual_over_uniform_prediction=ratio(actual,expected),
                           real_photon_model_ratio=ratio(b['model_average_pb'],r['model_average_pb'])))
    check(report,'energy_common_bins','PASS' if output else 'FAIL','Compare only identical W-bin boundaries.',dict(common_bins=len(output)))
    check(report,'energy_bin_edge_policy','WARN' if edge_differences else 'PASS',
          'A final bin includes its upper edge; an interior bin excludes it. Exact-edge populations can differ even with matching boundaries.',edge_differences)
    if reference.get('status') in ('FAIL','ERROR'):
        check(report,'reference_report_status','WARN','Reference audit contains failures; inspect before interpreting energy ratios.',reference.get('status'))
    return dict(reference_report=args.reference_report,reference_status=reference.get('status'),settings_mismatches=mismatched,rows=output)


def write_report(report,out):
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    (out/'audit.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n')
    csv_write(out/'w_bin_audit.csv',report.get('bin_comparison',[]))
    weights,groups,cutflow,branches=[],[],[],[]
    for label,data in report['files'].items():
        for branch,s in data.get('weights',{}).items():
            weights.append(dict(file=label,branch=branch,**{k:v for k,v in s.items() if k!='groups'}))
            groups.extend(dict(file=label,branch=branch,**row) for row in s['groups'])
        if data.get('physics'):
            cutflow.extend(dict(file=label,selection=k,**v) for k,v in data['physics']['cutflow'].items())
    for name,b in (report.get('pair') or {}).get('branches',{}).items():
        branches.append(dict(branch=name,**b))
    csv_write(out/'weight_summary.csv',weights)
    csv_write(out/'weight_groups.csv',groups)
    csv_write(out/'cutflow.csv',cutflow)
    csv_write(out/'branch_comparison.csv',branches)
    csv_write(out/'saved_csv_comparison.csv',(report.get('saved_csv_check') or {}).get('rows',[]))
    csv_write(out/'energy_comparison.csv',(report.get('cross_energy') or {}).get('rows',[]))
    lines=['# Upsilon normalization audit','', '**Internal consistency status: '+report['status']+'**','',
           'This is a read-only audit. No ROOT file or framework setting was changed.', '',
           'The assigned cross section is an external hypothesis. PASS does not certify its production provenance.', '',
           '## What enters the calculation','',
           '`weight` is used for sigma(W). `weight_proj` is inspected but does not enter this W extraction.',
           'The prepared file replaces original `weight` with assigned visible sigma × tuple luminosity / all generated rows.',
           'For supported bins the same-MC correction R/(R/T) equals T. This closure cannot independently validate normalization.',
           'The cross section in pb is T / (1000 × L[fb^-1] × B_ee × integrated photon flux).', '',
           'Uniform replacement and a global rescaling of original weights are separate diagnostic hypotheses.',
           'Uniform replacement assumes unweighted generator events and complete generated rows. The tuple alone cannot prove that.', '',
           '## Inputs and settings','', '```json',json.dumps(report['settings'],indent=2),'```','']
    for label,data in report['files'].items():
        s=data.get('weights',{}).get('weight',{})
        lines.extend(['- **'+label+'**: `'+str(data['metadata'].get('path'))+'`; '+str(data['metadata'].get('entries'))+
                      ' rows; sum(weight)='+str(s.get('sumw'))+'; literal sum(weight)/L='+str(s.get('literal_sumw_over_lumi_fb'))+' fb.'])
    lines+=['','## Checks','','| Status | Check | Result |','|---|---|---|']
    for c in report['checks']:
        lines.append('| '+c['status']+' | '+c['id']+' | '+c['detail'].replace('|','/')+' |')
        if c['status'] in ('WARN','FAIL','ERROR') and c.get('evidence') is not None:
            lines.append('| | evidence | `'+json.dumps(c['evidence'],sort_keys=True).replace('|','/')[:2400]+'` |')
    lines+=['','## Interpretation','',
            '- `weight_summary.csv`: totals, invalid weights and implied rates; these are arithmetic checks, not independently generated cross sections.',
            '- `weight_groups.csv`: repeated weight values and populations. Groups do not uniquely identify generator jobs.',
            '- `w_bin_audit.csv`: raw, actual normalized, uniform-target and globally rescaled predictions; flux and branching factors are explicit.',
            '- `branch_comparison.csv`: row-by-row preservation checks for non-weight branches.',
            '- `cutflow.csv`: generated denominator and every truth/reco selection.',
            '- `saved_csv_comparison.csv`: saved central values compared with recomputation, when provided.',
            '- `energy_comparison.csv`: cross-energy factor decomposition, when a reference report is supplied.',
            '- `audit.json`: complete metadata, beam hypotheses, field-level differences and all findings.', '',
            'The current-source real-photon eSTARlight benchmark is available only for 1S. Its fields are null for 2S/3S. Agreement does not verify production revision, card ranges, or finite-Q2 generation.',
            'The independent flux integration tests the framework EPA formula. It does not run eSTARlight or read its production flux table.',
            'To establish the physical rate, retain the generator revision, exact card and generator rate/log for each energy, including whether the rate includes B_ee.',
            'This audit targets normalization and W central values. It does not revalidate unfolding responses or statistical confidence intervals.','']
    if report.get('saved_csv_check',{}):
        lines += ['CSV factor hints: '+json.dumps(report['saved_csv_check'].get('factor_hints',[])), '']
    (out/'REPORT.md').write_text('\n'.join(lines))


def parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--merged',required=True)
    p.add_argument('--state',choices=tuple(STATE_DEFAULTS),default='1S',help='State in this single-state generated sample.')
    p.add_argument('--normalized')
    p.add_argument('--run-dir',help='Discover exact standard artifacts/inputs, cross_sections and normalization files in one run.')
    p.add_argument('--csv',help='Full-sample bin-by-bin sigma(W) CSV for --state; not a train/test unfolded CSV.')
    p.add_argument('--normalization-audit',help='External preparation audit JSON.')
    p.add_argument('--reference-report',help='audit.json from the other energy; ratios are current/reference.')
    p.add_argument('--out-dir',required=True)
    p.add_argument('--tree',default='UpsilonTree')
    p.add_argument('--electron-energy',type=float,required=True)
    p.add_argument('--proton-energy',type=float,required=True)
    p.add_argument('--assigned-visible-xsec-fb',type=float,required=True,help='Explicit current framework rate hypothesis, including ee branching.')
    p.add_argument('--generator-visible-xsec-fb',type=float,help='Optional independently supplied generator visible rate, same generated phase space.')
    p.add_argument('--generator-source',help='Provenance for the optional generator rate; recorded, not authenticated by this script.')
    p.add_argument('--lumi-invfb',type=float,default=10.)
    p.add_argument('--csv-yield-lumi-invfb',type=float,help='Only for legacy CSV lacking yield_lumi_invfb; otherwise metadata wins.')
    p.add_argument('--w-edges',help='Comma-separated bin edges; default current 9x130:20..70; 9x275:20..100, step5.')
    p.add_argument('--q2-cut',type=float,help='Event Q2 upper limit in GeV^2; July card defaults: 1 for 1S, 0.01 for 2S/3S.')
    p.add_argument('--q2-max',type=float,help='Flux Q2 upper limit in GeV^2; July card defaults: 1 for 1S, 0.01 for 2S/3S.')
    p.add_argument('--br-ee', '--branching-fraction', dest='br_ee', type=float,
                   help='Visible-to-inclusive branching fraction. Defaults: 1S=0.0238, 2S=0.0191, 3S=0.0218; --br-ee 0.0234 reproduces the archived 1S convention.')
    p.add_argument('--mass',type=float,help='State mass in GeV; defaults: 1S=9.46030, 2S=10.02326, 3S=10.35520.')
    p.add_argument('--mass-sigma',type=float,help='Reconstruction selection width in GeV; defaults: 1S=0.053, 2S=0.056, 3S=0.059.')
    p.add_argument('--mass-window-sigma',type=float,default=3.)
    p.add_argument('--true-mass-window',type=float,default=.2)
    p.add_argument('--sigma-epz-center',type=float)
    p.add_argument('--sigma-epz-halfwidth',type=float)
    p.add_argument('--min-eff',type=float,default=.05)
    p.add_argument('--chunk-size',type=int,default=100000)
    p.add_argument('--hash-inputs',action='store_true',help='SHA256 entire ROOT files; may take time over EOS.')
    return p


def main(argv=None):
    args=resolve_state_defaults(parser().parse_args(argv))
    if args.assigned_visible_xsec_fb<=0 or not math.isfinite(args.assigned_visible_xsec_fb) or args.chunk_size<=0:
        raise SystemExit('Cross section and chunk size must be positive.')
    if args.csv_yield_lumi_invfb is not None and (not math.isfinite(args.csv_yield_lumi_invfb) or args.csv_yield_lumi_invfb<=0):
        raise SystemExit('--csv-yield-lumi-invfb must be finite and positive.')
    if args.generator_visible_xsec_fb is not None and (not math.isfinite(args.generator_visible_xsec_fb) or args.generator_visible_xsec_fb<=0 or not args.generator_source):
        raise SystemExit('--generator-visible-xsec-fb requires a positive value and --generator-source.')
    energy='%gx%g'%(args.electron_energy,args.proton_energy)
    discovery=[]
    if args.run_dir:
        base=Path(args.run_dir)/'artifacts'
        for name,path in [('normalized',base/'inputs'/(energy+'_normalized.root')),
                          ('csv',base/'cross_sections'/(energy+'_sigma_gamma_p_'+args.state.lower()+'.csv')),
                          ('normalization_audit',base/'normalization'/(energy+'_input_audit.json'))]:
            if not getattr(args,name):
                if path.is_file():
                    setattr(args,name,str(path))
                else:
                    discovery.append(str(path))
    edges=[float(v) for v in args.w_edges.split(',')] if args.w_edges else list(np.arange(20, (70 if args.proton_energy==130 else 100)+.1,5))
    settings={k:getattr(args,k) for k in ('state','electron_energy','proton_energy','q2_cut','q2_max','lumi_invfb','br_ee','mass','mass_sigma','mass_window_sigma','true_mass_window','min_eff')}
    settings.update(w_edges=edges,sigma_epz_center=args.sigma_epz_center if args.sigma_epz_center is not None else 2*args.electron_energy,
                    sigma_epz_halfwidth=args.sigma_epz_halfwidth if args.sigma_epz_halfwidth is not None else .2*args.electron_energy)
    settings=validate_settings(settings)
    settings['assigned_visible_xsec_fb']=args.assigned_visible_xsec_fb
    report=dict(schema_version=1,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                status='ERROR',settings=settings,checks=[],files={},pair=None,bin_comparison=[],saved_csv_check=None,
                cross_energy=None,physical_rate_provenance_verified=False,
                scope='read-only original-versus-prepared normalization and full-sample W central-value audit')
    inputs=[args.merged]+([args.normalized] if args.normalized else [])
    before={}
    # Validate before the report-writing exception handler: never overwrite an input.
    out=Path(args.out_dir).resolve()
    reserved={'audit.json','REPORT.md','w_bin_audit.csv','weight_summary.csv','weight_groups.csv','cutflow.csv',
              'branch_comparison.csv','saved_csv_comparison.csv','energy_comparison.csv'}
    for path in inputs+[p for p in (args.csv,args.normalization_audit,args.reference_report) if p]:
        resolved=Path(path).resolve()
        if resolved.parent==out and resolved.name in reserved:
            raise SystemExit('Output would overwrite a supplied input: '+str(path))
    try:
        if args.normalized and os.path.realpath(args.merged)==os.path.realpath(args.normalized):
            raise ValueError('Merged and normalized paths refer to the same file.')
        for label,path in [('merged',args.merged),('normalized',args.normalized)]:
            if not path:
                continue
            print('Reading %s: %s'%(label,path),flush=True)
            st=os.stat(path)
            before[path]=(st.st_size,st.st_mtime_ns)
            meta=inspect_root(path,tree=args.tree,hash_input=args.hash_inputs)
            if meta.get('status')=='ERROR':
                report['files'][label]=dict(metadata=meta)
                raise ValueError('Cannot inspect '+label+': '+str(meta.get('errors')))
            report['files'][label]=scan(path,meta,settings,args)
            file_checks(report,label,report['files'][label],args)
        if discovery:
            check(report,'run_discovery','FAIL' if not args.normalized else 'WARN','Some expected run artifacts were absent; specify scenario paths explicitly when needed.',discovery)
        if args.normalized:
            print('Comparing all ordered non-weight branches and weight ratios...',flush=True)
            report['pair']=compare_root_pairs(args.merged,args.normalized,tree=args.tree,step_size=args.chunk_size)
            pair=report['pair']
            changed=[k for k,v in pair.get('branches',{}).items() if k not in ('weight','weight_proj') and (v.get('changed_entries',0) or v.get('type_equal') is False)]
            bad=changed or pair.get('missing_from_normalized') or pair.get('missing_from_merged') or pair.get('entries',{}).get('merged')!=pair.get('entries',{}).get('normalized')
            skipped=[k for k,v in pair.get('branches',{}).items() if v.get('status') in ('SKIP','ERROR')]
            check(report,'nonweight_branch_preservation','FAIL' if bad or pair.get('status')=='ERROR' else 'WARN' if skipped else 'PASS',
                  'Prepared file must retain all non-weight event values, branch names and row order.',
                  dict(changed=changed,missing_from_normalized=pair.get('missing_from_normalized'),missing_from_merged=pair.get('missing_from_merged'),skipped=skipped))
            invariance=pair.get('projective_factor_invariance',{})
            mapped={'PASS':'PASS','CHANGED':'FAIL','PARTIAL':'WARN','SKIP':'NOT_CHECKED','ERROR':'ERROR'}
            check(report,'projective_factor_preservation',mapped.get(invariance.get('status'),'NOT_CHECKED'),
                  'weight_proj/weight should be preserved; weight_proj is not used in sigma(W).',invariance)
            ratios=pair.get('weight_ratio',{})
            check(report,'weight_ratio_inspection',{'OK':'PASS','PARTIAL':'WARN','SKIP':'NOT_CHECKED','ERROR':'ERROR'}.get(ratios.get('status'),'WARN'),
                  'Original-to-normalized event weight ratios; replacement factors need not be constant.',ratios)
        metadata_checks(report,args)
        print('Checking W yields, independent flux integration and beam hypotheses...',flush=True)
        report['bin_comparison']=build_bins(report,settings,args)
        merged_physics = report['files']['merged'].get('physics') or {}
        truth_mass = merged_physics.get('cutflow', {}).get('truthmass')
        if truth_mass is not None:
            check(report, 'selected_state_population', 'PASS' if truth_mass['count'] > 0 else 'FAIL',
                  'The supplied single-state sample must contain events in the requested '+args.state+' truth-mass window.',
                  dict(state=args.state, mass=args.mass, true_mass_window=args.true_mass_window,
                       count=truth_mass['count'], fraction_of_all_rows=truth_mass.get('fraction_of_all_rows')))
        benchmark = merged_physics.get('conventions', {}).get('model_status')
        if benchmark == 'NOT_APPLICABLE':
            check(report, 'real_photon_benchmark', 'NOT_APPLICABLE',
                  'The available real-photon benchmark describes Upsilon(1S). It is not evaluated for the selected state or a different mass convention; model ratios are null.')
        report['beam_formula_check']=beam_formula_check(args.merged,settings,tree=args.tree,step_size=args.chunk_size)
        beam=report['beam_formula_check']
        requested=beam.get('candidates',{}).get(beam.get('requested_candidate'),{})
        deficient=any((requested.get(k,{}).get('matching_fraction') or 0)<.99 for k in ('true_W','true_Q2'))
        check(report,'beam_reconstruction', 'NOT_CHECKED' if beam.get('status') in ('SKIP','ERROR') else 'WARN' if deficient else 'PASS',
              'Compare supplied beams with alternate beam hypotheses using stored truth scattered-electron kinematics; see audit.json. This is a kinematic diagnostic.',report['beam_formula_check'])
        check(report,'q2_event_flux_consistency','PASS' if close(args.q2_cut,args.q2_max) else 'WARN',
              'Event Q2 cut and flux Q2 upper limit should describe the same selected phase space.',
              dict(event_q2_max=args.q2_cut,flux_q2_max=args.q2_max))
        report['saved_csv_check']=csv_check(report,args)
        report['cross_energy']=energy_comparison(report,args)
        if args.generator_visible_xsec_fb is not None:
            numeric_check(report,'external_generator_rate',args.assigned_visible_xsec_fb,args.generator_visible_xsec_fb,
                          'Assigned versus user-supplied independent visible generator rate; same generated phase space is required.',rtol=1e-4)
            report['external_generator_rate']=dict(visible_xsec_fb=args.generator_visible_xsec_fb,source=args.generator_source,provenance_authenticated=False)
        else:
            check(report,'external_generator_rate','NOT_CHECKED','No independent production generator rate provided. Original/normalized weights alone cannot certify the assigned cross section.')
        check(report,'rate_provenance','WARN','Assigned visible rate and generator-source description are explicit inputs. This audit checks their numerical consistency but does not authenticate the executable used for the supplied EOS production. The rate and branching fraction must belong to the selected '+args.state+' state and its generated Q2 phase space.')
        after={p:(os.stat(p).st_size,os.stat(p).st_mtime_ns) for p in before}
        check(report,'inputs_unchanged','PASS' if before==after else 'FAIL','Input ROOT size and modification time are unchanged during this read-only audit.')
        statuses={c['status'] for c in report['checks']}
        report['status']='ERROR' if 'ERROR' in statuses else 'FAIL' if 'FAIL' in statuses else 'WARN' if 'WARN' in statuses else 'PASS'
    except Exception as exc:
        check(report,'audit_exception','ERROR',str(exc))
        report['status']='ERROR'
    write_report(report,args.out_dir)
    print('Audit status: '+report['status'])
    print('Report: '+str(Path(args.out_dir)/'REPORT.md'))
    print('Machine-readable evidence: '+str(Path(args.out_dir)/'audit.json'))
    return 2 if report['status']=='ERROR' else 1 if report['status']=='FAIL' else 0


if __name__=='__main__':
    sys.exit(main())
