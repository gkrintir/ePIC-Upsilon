#!/usr/bin/env python3
"""Run optional fEMC x/Q2 loss maps without rerunning the full note workflow."""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / 'scripts/validation/plot_femc_xq2.py'
SCENARIOS = ('nominal', 'remove_inner', 'remove_outer')
LABELS = {'nominal': 'Nominal fEMC', 'remove_inner': 'Remove inner', 'remove_outer': 'Remove outer'}


def positive(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError('must be positive and finite')
    return result


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def resolve_configuration(args, catalog, energy):
    state = catalog['states'][args.state]
    beam = catalog['beams'][energy]
    config = {key: catalog[key] for key in ('tree', 'electron_energy', 'truth_mass_window',
              'q2_edges', 'xb_edges', 'xv_edges', 'normalization_source', 'production_family')}
    config.update(beam)
    config.update(state=args.state, energy_label=energy + ' GeV',
                  visible_xsec_fb=state['visible_xsec_fb'][energy], q2_max=state['q2_max'],
                  lumi_invfb=args.lumi if args.lumi is not None else catalog['lumi_invfb'],
                  selection_scope=args.scope or catalog['selection_scope'],
                  min_nominal_mc=args.min_nominal_mc if args.min_nominal_mc is not None else catalog['min_nominal_mc'],
                  hash_inputs=args.hash_inputs)
    samples = []
    for key in SCENARIOS:
        explicit = getattr(args, key)
        directory = getattr(args, key + '_dir') or catalog['directories'][key]
        path = explicit or catalog['input_pattern'].format(state=args.state, energy=energy, directory=directory)
        samples.append({'key': key, 'label': LABELS[key], 'path': str(Path(path).expanduser().resolve())})
    config['samples'] = samples
    if config['selection_scope'] == 'full' and (args.w_min is not None or args.w_max is not None):
        raise ValueError('--w-min/--w-max apply only to the analysis scope')
    if args.w_min is not None:
        config['w_min'] = args.w_min
    if args.w_max is not None:
        config['w_max'] = args.w_max
    if not config['w_min'] < config['w_max']:
        raise ValueError('w_min must be smaller than w_max')
    if config['min_nominal_mc'] < 1:
        raise ValueError('min_nominal_mc must be at least 1')
    if len({sample['path'] for sample in samples}) != len(samples):
        raise ValueError('Nominal and removal scenarios must refer to distinct ROOT files')
    return config


def outputs_complete(folder):
    try:
        audit = json.loads((folder / 'audit.json').read_text())
    except (OSError, ValueError):
        return False
    if audit.get('status') not in ('PASS', 'WARN'):
        return False
    required = ['audit.json', 'bins.csv']
    for coord in ('xb', 'xv'):
        required.append(coord + '_q2_bins.csv')
        for kind in ('comparison', 'loss_remove_inner', 'loss_remove_outer'):
            required.extend(coord + '_q2_' + kind + '.' + extension for extension in ('png', 'pdf'))
    return all((folder / name).is_file() and (folder / name).stat().st_size > 0 for name in required)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/femc_xq2.json')
    parser.add_argument('--energy', choices=('all', '9x130', '9x275'), default='all')
    parser.add_argument('--state', choices=('1S', '2S', '3S'), default='1S')
    parser.add_argument('--lumi', type=positive, help='Projected luminosity in fb^-1; catalog default 10')
    parser.add_argument('--scope', choices=('analysis', 'full'), help='analysis applies W cuts; full retains generated threshold region')
    parser.add_argument('--w-min', type=positive)
    parser.add_argument('--w-max', type=positive)
    parser.add_argument('--min-nominal-mc', type=int)
    parser.add_argument('--run-id', help='New output directory ID; default contains state and timestamp')
    parser.add_argument('--plan', action='store_true', help='Print resolved inputs/settings without reading ROOT or writing outputs')
    parser.add_argument('--hash-inputs', action='store_true')
    parser.add_argument('--force', action='store_true', help='Rerun completed maps only if their configuration/input contract still matches')
    for key in SCENARIOS:
        option = key.replace('_', '-')
        parser.add_argument('--' + option, help='Explicit ROOT path, requires a single --energy')
        parser.add_argument('--' + option + '-dir', help='Override this scenario directory name at both energies')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.energy == 'all' and any(getattr(args, key) for key in SCENARIOS):
        parser.error('Explicit ROOT paths require --energy 9x130 or --energy 9x275')
    if args.scope == 'full' and (args.w_min is not None or args.w_max is not None):
        parser.error('--w-min/--w-max apply to --scope analysis; full intentionally omits those W cuts')
    try:
        catalog = json.loads(args.config.read_text())
        energies = list(catalog['beams']) if args.energy == 'all' else [args.energy]
        configurations = {energy: resolve_configuration(args, catalog, energy) for energy in energies}
    except (KeyError, ValueError, OSError) as error:
        parser.error(str(error))
    run_id = args.run_id or ('femc-xq2-' + args.state.lower() + '-' + dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,149}', run_id):
        parser.error('--run-id must contain letters, digits, underscores, dots or hyphens')
    run_dir = ROOT / 'results/generated/runs' / run_id
    if args.plan:
        for energy, config in configurations.items():
            print(energy + '\n' + json.dumps(config, indent=2))
            print('Output: ' + str(run_dir / 'artifacts/femc_xq2' / energy))
        return 0
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / '.femc_xq2.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error('This run ID is already active; choose another ID or wait for it to finish')
        manifest_path = run_dir / 'femc_xq2_manifest.json'
        manifest = {'run_id': run_id, 'state': args.state, 'status': 'RUNNING', 'tasks': {}}
        if manifest_path.exists():
            previous = json.loads(manifest_path.read_text())
            if previous.get('state') != args.state:
                parser.error('This run ID belongs to another Upsilon state; use a fresh --run-id')
            manifest['tasks'].update(previous.get('tasks', {}))
        source_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in (ENGINE, Path(__file__))}
        for energy, config in configurations.items():
            output = run_dir / 'artifacts/femc_xq2' / energy
            print('fEMC x/Q2 ' + args.state + ', ' + energy, flush=True)
            missing = [sample['path'] for sample in config['samples'] if not Path(sample['path']).is_file()]
            if missing:
                manifest['tasks'][energy] = {'status': 'BLOCKED_INPUT', 'missing': missing}
                print('Missing ROOT files:\n' + '\n'.join(missing), flush=True)
                write_json(manifest_path, manifest)
                continue
            fingerprints = []
            for sample in config['samples']:
                path = Path(sample['path']); stat = path.stat()
                fingerprint = {'key': sample['key'], 'path': str(path), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
                if args.hash_inputs:
                    fingerprint['sha256'] = sha256(path)
                fingerprints.append(fingerprint)
            contract = {'configuration': config, 'inputs': fingerprints, 'source_sha256': source_hashes}
            contract_path = output / 'contract.json'
            if contract_path.exists() and json.loads(contract_path.read_text()) != contract:
                manifest['tasks'][energy] = {'status': 'BLOCKED_CONTRACT', 'error': 'Settings, inputs or code changed; use a fresh --run-id'}
                print('BLOCKED_CONTRACT: use a fresh --run-id', flush=True)
                write_json(manifest_path, manifest)
                continue
            if not contract_path.exists() and output.exists() and any(output.iterdir()):
                manifest['tasks'][energy] = {'status': 'BLOCKED_CONTRACT', 'error': 'Existing outputs lack a contract; use a fresh --run-id'}
                write_json(manifest_path, manifest)
                continue
            if outputs_complete(output) and not args.force:
                manifest['tasks'][energy] = {'status': 'CACHED', 'output': str(output)}
                print('CACHED ' + str(output), flush=True)
                write_json(manifest_path, manifest)
                continue
            write_json(contract_path, contract)
            resolved_path = output / 'resolved_config.json'
            write_json(resolved_path, config)
            command = [sys.executable, str(ENGINE), '--config', str(resolved_path), '--out-dir', str(output)]
            print(shlex.join(command), flush=True)
            log_path = output / 'run.log'
            with log_path.open('w') as log:
                result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            done = result.returncode == 0 and outputs_complete(output)
            manifest['tasks'][energy] = {'status': 'DONE' if done else 'FAILED', 'returncode': result.returncode,
                                         'output': str(output), 'log': str(log_path)}
            print(manifest['tasks'][energy]['status'] + ' ' + str(output), flush=True)
            if not done:
                print(log_path.read_text()[-5000:], flush=True)
            write_json(manifest_path, manifest)
        manifest['status'] = 'COMPLETE' if all(task['status'] in ('DONE', 'CACHED') for task in manifest['tasks'].values()) else 'PARTIAL'
        write_json(manifest_path, manifest)
        print('Manifest: ' + str(manifest_path), flush=True)
        return 0 if manifest['status'] == 'COMPLETE' else 1


if __name__ == '__main__':
    raise SystemExit(main())
