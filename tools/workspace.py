#!/usr/bin/env python3
"""Create isolated analysis workspaces and invoke the preserved Upsilon backend."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_BASE = '/eos/cms/store/group/phys_top/gkrintir/upsilon_data'
SELECTIONS = {'inclusive': 'inclusive', 'proton-tagged': 'required'}
STATES = ('1S', '2S', '3S')
SKIP_DIRS = {'results', '__pycache__', '.git', '.ipynb_checkpoints'}
SKIP_FILES = {'.proton_category.json', '.workspace.json', 'images.tar'}


def safe_name(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,139}', str(value)):
        raise ValueError('Use one name containing letters, digits, dot, underscore or hyphen (maximum 140 characters)')
    return str(value)


def inside(path, parent):
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def no_links(path):
    path = Path(path).absolute()
    for part in [path, *path.parents]:
        if part.is_symlink():
            raise ValueError('Symlink workspace/code path is not supported: ' + str(part))
    return path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_files(source):
    files = {}
    source = no_links(source).resolve()
    for folder, dirs, names in os.walk(source, followlinks=False):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
        for dirname in dirs:
            if (Path(folder)/dirname).is_symlink():
                raise ValueError('Source contains a symlink: ' + str(Path(folder)/dirname))
        for name in sorted(names):
            if name in SKIP_FILES or name.endswith(('.pyc', '.bak', '~')) or '.bak.' in name:
                continue
            path = Path(folder)/name
            if path.is_symlink():
                raise ValueError('Source contains a symlink: ' + str(path))
            if path.suffix == '.root':
                raise ValueError('Analysis source contains a ROOT data file: ' + str(path))
            files[str(path.relative_to(source))] = path
    if not files:
        raise ValueError('Analysis source directory is empty: ' + str(source))
    return files


def read_settings(repo, config=None):
    repo = Path(repo).resolve()
    path = Path(config).expanduser() if config else repo/'configs/local.json'
    if config and not path.is_absolute():
        path = Path.cwd()/path
    if not path.exists() and not config:
        path = repo/'configs/local.example.json'
    settings = json.loads(path.read_text()) if path.exists() else {}
    if config and not path.exists():
        raise FileNotFoundError(path)
    settings.setdefault('workspace_root', 'workspaces')
    settings.setdefault('input_base', DEFAULT_INPUT_BASE)
    root = Path(os.path.expandvars(os.path.expanduser(str(settings['workspace_root']))))
    if not root.is_absolute():
        root = repo/root
    settings['workspace_root'] = str(no_links(root).resolve())
    base = Path(os.path.expandvars(os.path.expanduser(str(settings['input_base']))))
    if not base.is_absolute():
        raise ValueError('input_base must be an absolute path (e.g. the EOS upsilon_data directory)')
    settings['input_base'] = str(base)
    return settings


def workspace_path(repo, settings, name):
    source = Path(repo).resolve()/'analysis'
    path = no_links(Path(settings['workspace_root'])/safe_name(name)).resolve()
    if inside(path, source) or inside(source, path):
        raise ValueError('Workspace and analysis source must be distinct, disjoint directories')
    return path


def init_workspace(repo, settings, name, selection, dry_run=False):
    repo = Path(repo).resolve()
    source = repo/'analysis'
    destination = workspace_path(repo, settings, name)
    if selection not in SELECTIONS:
        raise ValueError('Unknown selection')
    if destination.exists():
        raise FileExistsError('Workspace already exists; it will not be overwritten: ' + str(destination))
    files = source_files(source)
    for required in ('workflow/run_all_states.py', 'workflow/proton_selection.py', 'workflow/proton_category_isolation.py'):
        if required not in files:
            raise FileNotFoundError('Incomplete analysis backend: ' + required)
    hashes = {name: sha(path) for name, path in files.items()}
    print('Source:', source)
    print('Workspace:', destination)
    print('Selection:', selection, '(' + SELECTIONS[selection] + ')')
    print('Files:', len(files), '; existing code and results are unchanged')
    if dry_run:
        print('DRY RUN: no directories created')
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    no_links(destination.parent)
    destination.mkdir(exist_ok=False)
    receipt = {'schema_version': 1, 'installation_status': 'INSTALLING',
               'workspace': str(destination), 'source': str(source), 'selection': selection,
               'created_utc': datetime.now(timezone.utc).isoformat(),
               'source_sha256': hashes, 'orchestrator_sha256': sha(Path(__file__))}
    marker = destination/'.workspace.json'
    marker.write_text(json.dumps(receipt, indent=2)+'\n')
    # Leave an explicit INSTALLING receipt on failure. Never delete a workspace.
    for relative, original in files.items():
        target = destination/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(original.read_bytes())
        target.chmod(original.stat().st_mode & 0o777)
        if sha(target) != hashes[relative]:
            raise ValueError('Source changed while copying: ' + relative)
    if selection == 'proton-tagged':
        category = {'schema_version': 1, 'package_version': 3,
                    'installation_status': 'COMPLETE', 'mode': 'required',
                    'source_root': str(source), 'category_root': str(destination),
                    'created_by': 'ePIC-Upsilon tools/workspace.py'}
        (destination/'.proton_category.json').write_text(json.dumps(category, indent=2)+'\n')
    receipt['installation_status'] = 'COMPLETE'
    marker.write_text(json.dumps(receipt, indent=2)+'\n')
    print('Initialized:', destination)
    return destination


def load_workspace(repo, settings, name, verify=True):
    path = workspace_path(repo, settings, name)
    marker = path/'.workspace.json'
    no_links(marker)
    if not marker.is_file():
        raise FileNotFoundError('Workspace is not initialized: ' + str(path) + '. Use init first.')
    meta = json.loads(marker.read_text())
    if meta.get('schema_version') != 1 or meta.get('installation_status') != 'COMPLETE':
        raise ValueError('Workspace initialization is incomplete; inspect it and initialize a fresh name')
    if meta.get('workspace') != str(path) or meta.get('source') != str(Path(repo).resolve()/'analysis'):
        raise ValueError('Workspace was moved or belongs to another source checkout; initialize a fresh workspace')
    if meta.get('selection') not in SELECTIONS:
        raise ValueError('Unknown recorded workspace selection')
    if verify:
        current = {k: sha(v) for k, v in source_files(Path(repo)/'analysis').items()}
        copied = {k: sha(v) for k, v in source_files(path).items()}
        if current != meta['source_sha256'] or copied != meta['source_sha256']:
            raise ValueError('Workspace or repository analysis source changed. Preserve this workspace and initialize a fresh name.')
        if meta.get('orchestrator_sha256') != sha(Path(__file__)):
            raise ValueError('Workspace launcher version changed; initialize a fresh workspace before running')
    category = path/'.proton_category.json'
    if meta['selection'] == 'proton-tagged':
        no_links(category)
        data = json.loads(category.read_text())
        if (data.get('category_root') != str(path) or data.get('source_root') != meta['source']
                or data.get('mode') != 'required' or data.get('installation_status') != 'COMPLETE'):
            raise ValueError('Proton category marker is inconsistent')
    elif category.exists():
        raise ValueError('Inclusive workspace must not contain a proton-tagged category marker')
    no_links(path/'results/generated')
    return path, meta


def selected_states(values):
    result = list(STATES) if values == ['all'] else [value.upper() for value in values]
    if not result or len(set(result)) != len(result) or any(value not in STATES for value in result):
        raise ValueError('--states must be all or unique states separated by spaces, e.g. 1S 2S 3S')
    return result


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return number


def execute(command, cwd, selection, dry_run=False):
    print('Selection:', selection)
    print('Command:', shlex.join([str(x) for x in command]), flush=True)
    if dry_run:
        print('DRY RUN: command not executed')
        return 0
    env = dict(os.environ, YPSILON_PROTON_SELECTION=SELECTIONS[selection], PYTHONDONTWRITEBYTECODE='1')
    return subprocess.run([str(x) for x in command], cwd=str(cwd), env=env).returncode


def analysis_command(args, settings, workspace):
    states = selected_states(args.states)
    safe_name(args.batch_id)
    command = [sys.executable, '-B', str(workspace/'workflow/run_all_states.py'), args.action,
               '--states', ','.join(states), '--energy', args.energy,
               '--stat-mode', args.stat_mode, '--projection-lumi', str(args.lumi),
               '--tuple-lumi', '10.0', '--batch-id', args.batch_id,
               '--set', 'all_states.input_pattern=' + settings['input_base'] + '/Upsilon{state}/{energy}/hiAcc/{directory}/merged_Upsilon_final.root']
    if args.action == 'run':
        command += ['--keep-going']
        if args.hash_inputs:
            command += ['--hash-inputs']
    return command


def check_source_batch(workspace, batch, states, energy, mode):
    safe_name(batch)
    slot = 'low' if energy == '9x130' else 'high'
    rates = {'1S': (119.685,210.136,.0238,1.), '2S': (29.297,53.133,.0191,.01), '3S': (23.447,42.969,.0218,.01)}
    inputs = {}
    for state in states:
        run = workspace/'results/generated/runs'/(batch+'-'+state.lower())
        no_links(run)
        manifest = json.loads((run/'run_manifest.json').read_text())
        policy = json.loads((run/'proton_selection.json').read_text()).get('proton_selection', {})
        if policy.get('mode') != mode:
            raise ValueError('Source batch selection differs for ' + state)
        conf = manifest['configuration']; beam = conf['beam'][slot]; analysis = conf['analysis']
        rate = rates[state][0 if energy == '9x130' else 1]
        checks = [(beam['electron_energy'], 9.), (beam['proton_energy'], 130. if energy == '9x130' else 275.),
                  (beam['visible_xsec_fb'],rate), (analysis['generator_branching_fraction'],rates[state][2]),
                  (analysis['q2_event_max'],rates[state][3]), (analysis['q2_flux_max'],rates[state][3])]
        if analysis.get('states') != state or not all(math.isclose(float(a),b,rel_tol=1e-9,abs_tol=1e-12) for a,b in checks):
            raise ValueError('Source batch state, beam, rate or Q2 provenance mismatch for ' + state)
        stage_path = run/'stages'/('prepare_'+energy)/'manifest.json'
        tasks = [t for t in manifest.get('tasks',[]) if t.get('task')=='prepare_'+energy]
        task = json.loads(stage_path.read_text()) if stage_path.exists() else (tasks[0] if len(tasks)==1 else {})
        if task.get('status') not in ('DONE','CACHED'):
            raise ValueError('Prepared input has no successful task record for '+state+'/'+energy)
        path = run/'artifacts/inputs'/(energy+'_normalized.root')
        no_links(path)
        if not path.is_file():
            raise FileNotFoundError('Missing normalized input: '+str(path))
        inputs[state] = path
    return inputs


def reopt_command(args, repo, workspace, meta):
    states = selected_states(args.states)
    safe_name(args.run_id)
    mode = SELECTIONS[meta['selection']]
    output = workspace/'results/generated/runs'/args.run_id
    no_links(output)
    if args.action == 'run' and output.exists():
        if not output.is_dir():
            raise FileExistsError('Output path is not a directory: '+str(output))
        if not (output/'run_manifest.json').is_file() and any(output.iterdir()):
            raise FileExistsError('Output directory is nonempty without a reoptimization manifest; use a fresh run ID')
    inputs = check_source_batch(workspace,args.source_batch,states,args.energy,mode)
    if args.energy == '9x130':
        if mode != 'required':
            raise ValueError('The current 9x130 reoptimization extension requires proton tagging')
        runner = Path(repo)/'extensions/reoptimization_9x130/run_reoptimization_9x130.py'
        command = [sys.executable,'-B',str(runner),args.action,'--framework-root',str(workspace),
                   '--source-batch',args.source_batch]
    else:
        command = [sys.executable,'-B',str(workspace/'workflow/run_reoptimization.py'),args.action,
                   '--framework-root',str(workspace)]
        for state,path in inputs.items():
            command += ['--input',state+'='+str(path)]
        # Remove old AFS defaults from the effective configuration as well as
        # overriding the actual input paths. This temporary config is persisted
        # by the backend as effective_config.json in the new run.
        cfg = json.loads((workspace/'config/reoptimization_9x275_wy.json').read_text())
        for state in STATES:
            cfg['states'][state]['run'] = str(workspace/'results/generated/runs'/(args.source_batch+'-'+state.lower()))
        cfg['source_batch'] = args.source_batch
        cfg['source_records_sha256'] = {state: sha(path.parents[2]/'run_manifest.json') for state,path in inputs.items()}
        args.generated_reopt_config = cfg
    command += ['--states',*states,'--observables',*args.observables,'--lumi',str(args.lumi),'--run-id',args.run_id]
    if args.hash_inputs:
        command += ['--hash-inputs']
    if args.resume:
        if args.action != 'run':
            raise ValueError('--resume applies only to run')
        command += ['--resume']
    return command


def show_status(workspace, run_id=None):
    base = workspace/'results/generated/runs'
    runs = [base/safe_name(run_id)] if run_id else sorted(base.glob('*'))
    found = False
    for run in runs:
        manifest_path = run/'run_manifest.json'
        if not manifest_path.is_file():
            continue
        found = True
        data = json.loads(manifest_path.read_text())
        print(run.name + ': ' + str(data.get('status',data.get('run_status','UNKNOWN'))))
        status = run/'artifacts/reoptimization/validation_status.csv'
        if status.is_file():
            with status.open(newline='') as stream:
                for row in csv.DictReader(stream):
                    print('  {}/{}: {} (inference {})'.format(row['state'],row['observable'],row.get('overall_status',row.get('status')),row.get('inference_status',row.get('stage_status','UNKNOWN'))))
        summary = data.get('continuation_summary',{})
        if summary.get('incomplete_tasks'):
            print('  Incomplete: '+', '.join(summary['incomplete_tasks']))
    if not found:
        print('No run manifests found in '+str(base))
    print('Execution COMPLETE does not mean physics validation passed.')
    return 0


def import_check(module, additional=()):
    """Actually load extension modules in a fresh process, including BLAS users."""
    code = ("import importlib,json,sys; "
            "main=importlib.import_module(sys.argv[1]); "
            "[importlib.import_module(name) for name in sys.argv[2:]]; "
            "print(json.dumps({'version':str(getattr(main,'__version__','unknown'))}))")
    try:
        result = subprocess.run([sys.executable, '-B', '-c', code, module, *additional],
                                capture_output=True, text=True, timeout=20,
                                env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
    except subprocess.TimeoutExpired:
        return False, 'import timed out after 20 seconds'
    if result.returncode:
        lines = (result.stderr or result.stdout or 'import process failed').strip().splitlines()
        return False, lines[-1][:500] if lines else 'import process failed'
    try:
        version = json.loads(result.stdout.strip().splitlines()[-1])['version']
    except (ValueError, KeyError, IndexError):
        version = 'unknown'
    return True, str(version)[:100]


def doctor(repo, settings, name=None):
    failed = False
    print('Python:',sys.version.split()[0],sys.executable)
    if sys.version_info < (3,9):
        print('MISSING: Python >=3.9');failed=True
    for module in ('numpy','scipy','pandas','matplotlib','mplhep','uproot'):
        extra = ('scipy.linalg', 'scipy.optimize') if module == 'scipy' else ()
        ok, detail = import_check(module, extra)
        label = module + (' (including scipy.linalg and scipy.optimize)' if extra else '')
        print(('IMPORT OK ' if ok else 'IMPORT FAILED ')+label+': '+detail)
        failed |= not ok
    if sys.version_info < (3,11):
        ok, detail = import_check('tomli')
        print(('IMPORT OK ' if ok else 'IMPORT FAILED ')+'tomli: '+detail);failed |= not ok
    root_ok, root_detail = import_check('ROOT')
    print(('OPTIONAL IMPORT OK ' if root_ok else 'OPTIONAL IMPORT UNAVAILABLE ')+'ROOT: '+root_detail)
    print('ROOT executable (production and legacy tasks):',shutil.which('root') or 'not found')
    print('Input base:',settings['input_base'],'(available)' if Path(settings['input_base']).is_dir() else '(not mounted here)')
    if name:
        path,meta=load_workspace(repo,settings,name)
        print('Workspace verified:',path,meta['selection'])
    print('Dependency imports were executed; ROOT event data and CERN production jobs were not run.')
    return int(failed)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',help='Local JSON settings (default configs/local.json or local.example.json)')
    sub=p.add_subparsers(dest='command',required=True)
    for name in ('init','doctor','analysis','reopt','status'):
        s=sub.add_parser(name)
        s.add_argument('--config',default=argparse.SUPPRESS,help='Local JSON settings')
        if name=='init':
            s.add_argument('--name',required=True);s.add_argument('--selection',choices=tuple(SELECTIONS),required=True)
            s.add_argument('--dry-run',action='store_true')
        elif name=='doctor':s.add_argument('--workspace')
        else:s.add_argument('--workspace',required=True)
        if name in ('analysis','reopt'):
            s.add_argument('action',choices=('plan','validate','run'))
            s.add_argument('--energy',choices=('9x130','9x275','all') if name=='analysis' else ('9x130','9x275'),required=True)
            s.add_argument('--states',nargs='+',default=['all'])
            s.add_argument('--lumi',type=positive,required=True,help='Target luminosity in fb^-1; always explicit')
            s.add_argument('--hash-inputs',action='store_true')
            s.add_argument('--dry-run',action='store_true',help='Print the command without executing it')
        if name=='analysis':
            s.add_argument('--batch-id',required=True)
            s.add_argument('--stat-mode',choices=('mc','projected'),default='projected')
        elif name=='reopt':
            s.add_argument('--source-batch',required=True)
            s.add_argument('--run-id',required=True)
            s.add_argument('--observables',nargs='+',choices=('pt2','abs_t_via_pt2','w','y'),default=['pt2','abs_t_via_pt2'])
            s.add_argument('--resume',action='store_true')
        elif name=='status':s.add_argument('--run-id')
    return p


def main(argv=None, repo=None):
    args=parser().parse_args(argv)
    repo=Path(repo or REPO).resolve()
    try:
        settings=read_settings(repo,args.config)
        if args.command=='init':
            init_workspace(repo,settings,args.name,args.selection,args.dry_run);return 0
        if args.command=='doctor':return doctor(repo,settings,args.workspace)
        workspace,meta=load_workspace(repo,settings,args.workspace,verify=args.command!='status')
        if args.command=='status':return show_status(workspace,args.run_id)
        command=(analysis_command(args,settings,workspace) if args.command=='analysis'
                 else reopt_command(args,repo,workspace,meta))
        if getattr(args, 'generated_reopt_config', None) is not None:
            if args.dry_run:
                print('A temporary 9x275 config will bind the listed source batch; saved in the run as effective_config.json.')
                return execute(command + ['--config', '<generated-9x275-config>'],workspace,meta['selection'],True)
            with tempfile.TemporaryDirectory(prefix='upsilon-reopt-config-') as folder:
                path=Path(folder)/'reoptimization_9x275.json'
                path.write_text(json.dumps(args.generated_reopt_config,indent=2)+'\n')
                return execute(command + ['--config',str(path)],workspace,meta['selection'])
        return execute(command,workspace,meta['selection'],args.dry_run)
    except (ValueError,FileNotFoundError,FileExistsError,KeyError,OSError) as error:
        print('ERROR: '+str(error),file=sys.stderr);return 2


if __name__=='__main__':
    raise SystemExit(main())
