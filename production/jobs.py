#!/usr/bin/env python3
"""New, explicit shared-filesystem production planner; never submits Condor jobs.

Each beam uses its own recovered, byte-exact tuple macro. Each batch freezes
the selected macro, its shared helpers and the worker script.
"""
import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
from itertools import islice
import io
import json
import math
import multiprocessing
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import NamedTuple
import zlib

PRODUCTION = Path(__file__).resolve().parent
REPO = PRODUCTION.parent
PRODUCER = 'DVMP_Upsilon_Analysis.C'
PRODUCERS = {'9x275': 'macros/' + PRODUCER, '9x130': 'macros/9x130/' + PRODUCER}
MACROS = (PRODUCER, 'ElectronFinder.C', 'ePICStyle.C', 'ePICStyle.h')
REMATCH = 'reco-rematch3011'
REMATCH_ENTRY = 'RunRematch3011.C'
REMATCH_PATCH = 'ElectronFinder.patch'
ARCHIVED_READER = 'ElectronFinder.stored.C'
SCENARIOS = ('full', 'remove_inner', 'remove_outer', 'keep_xpos', 'keep_yneg')
DEFAULT_CHECK_WORKERS = min(4, len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1))
IN_MEMORY_TUPLE_BYTES = 4 * 1024 * 1024


def ordered_checks(check, items, workers, *, executor=ThreadPoolExecutor, **executor_options):
    """Overlap filesystem reads, preserving manifest order and bounding outstanding work."""
    if workers < 1:
        raise ValueError('The number of check workers must be positive')
    if workers == 1:
        yield from map(check, items)
        return
    items = iter(items)
    with executor(max_workers=workers, **executor_options) as pool:
        pending = deque(pool.submit(check, item) for item in islice(items, 2 * workers))
        try:
            for item in items:
                yield pending.popleft().result()
                pending.append(pool.submit(check, item))
            while pending:
                yield pending.popleft().result()
        finally:
            for future in pending:
                future.cancel()


def add_check_workers(parser):
    parser.add_argument('--check-workers', type=int, default=DEFAULT_CHECK_WORKERS,
                        help='Concurrent per-file checks (default: %(default)s; 1 for serial checks)')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, value):
    path = Path(path)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def safe_path(value):
    path = Path(value).expanduser().absolute()
    # A simple shared-filesystem contract avoids ambiguous Condor/ROOT syntax.
    if not re.fullmatch(r'[A-Za-z0-9_./+\-]+', str(path)):
        raise ValueError('Paths must use letters, digits, underscore, dot, slash, plus or hyphen: ' + str(path))
    return path


def validate_work_location(path, repository=None):
    resolved = Path(path).resolve()
    repository = Path(repository or REPO).resolve()
    if resolved.is_relative_to(repository) and not resolved.is_relative_to((repository / 'workspaces').resolve()):
        raise ValueError('Inside the repository, generated work must live under ignored workspaces/')


def expected_output_path(batch, index):
    name = f'job_{index:06d}'
    if 'tuple_output_dir' in batch:
        return Path(batch['tuple_output_dir']) / (name + '.root')
    return Path(batch['work_dir']) / 'jobs' / name / 'output.root'


def submit_text(info, jobs, *, log_dir=None):
    work = Path(info['work_dir'])
    lines = ['# Newly generated wrapper, not the historical CERN submission.',
             'universe = vanilla', 'should_transfer_files = NO',
             'getenv = PATH,LD_LIBRARY_PATH,PYTHONPATH,ROOTSYS,ROOT_INCLUDE_PATH,CPATH,CPLUS_INCLUDE_PATH',
             'executable = ' + str(work / 'snapshot/run_job.sh'),
             'request_cpus = 1', 'request_memory = %dMB' % info['memory_mb'],
             'request_disk = %dMB' % info['disk_mb']]
    if info.get('accounting_group'):
        lines.append('+AccountingGroup = "' + info['accounting_group'] + '"')
    if info.get('job_flavour'):
        lines.append('+JobFlavour = "' + info['job_flavour'] + '"')
    # Indexes identify the original file or group, including in a sparse pending.sub.
    group_size = info.get('files_per_job', 1)
    job_dir = work / 'jobs' / 'job_$(batch_index)'
    initialdir, argument = (job_dir, f'{job_dir}/job.json') if group_size == 1 else (work, '$(batch_index)')
    log_dir = Path(log_dir) if log_dir is not None else work / 'logs'
    log_prefix = log_dir / ('$(batch_index)' if group_size == 1 else 'group_$(batch_index)')
    lines += [f'initialdir = {initialdir}', f'arguments = {argument}',
              f'output = {log_prefix}.out', f'error = {log_prefix}.err',
              f'log = {log_prefix}.log', 'queue batch_index from (']
    file_ids = (int(Path(path).parent.name.removeprefix('job_')) for path in jobs)
    group_ids = sorted({index // group_size for index in file_ids})
    lines.extend(f'{index:06d}' for index in group_ids)
    lines.append(')')
    return '\n'.join(lines) + '\n'


def check_beam_label(path, beam):
    # Campaign names include both 9x130 and 9To130, beside cuts such as q2_0to1.
    for electron, proton in re.findall(r'(?<!\d)(5|9|10|18)(?:x|on|to)(41|100|130|250|275)(?!\d)', str(path), re.I):
        if f'{int(electron)}x{int(proton)}' != beam:
            raise ValueError(f'Beam label {electron}x{proton} conflicts with {beam}: {path}')


def input_records(path, state, hash_inputs=False, *, beam=None):
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    if not lines:
        raise ValueError('Input list is empty')
    paths = [safe_path(line).resolve() for line in lines]
    if len(set(paths)) != len(paths):
        raise ValueError('Duplicate input ROOT files would duplicate events')
    records = []
    for line, p in zip(lines, paths):
        if beam is not None:
            check_beam_label(line, beam)
            check_beam_label(p, beam)
        if not p.is_file() or p.suffix != '.root':
            raise ValueError('Expected a mounted, existing .root file: ' + str(p))
        # Match the preserved producer path-based state inference exactly.
        inferred = next((s for s in ('1S', '2S', '3S') if s.lower() in str(p).lower()), None)
        if inferred != state:
            raise ValueError('Producer infers state from input path; ' + str(p) + ' does not identify ' + state)
        stat = p.stat()
        record = dict(path=str(p), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
        if hash_inputs:
            record['sha256'] = sha256(p)
        records.append(record)
    return records


def verify_beam_constants(path, beam):
    text = path.read_text()
    electron, proton = map(float, beam.split('x'))
    for name, value in [('E_e_beam', electron), ('E_p_beam', proton)]:
        match = re.search(r'^\s*const double\s+' + name + r'\s*=\s*([0-9.]+)\s*;', text, re.M)
        if not match or float(match.group(1)) != value:
            raise ValueError(f'Producer beam constant {name} conflicts with {beam}: {path}')


def verify_macro_sources(beam='9x275', electron_id='reco'):
    sources = json.loads((PRODUCTION / 'sources.json').read_text())
    expected = {row['path']: row['sha256'] for row in sources['files']}
    selected = {name: 'macros/' + name for name in MACROS}
    selected[PRODUCER] = PRODUCERS[beam]
    if electron_id == REMATCH:
        selected[ARCHIVED_READER] = selected.pop('ElectronFinder.C')
        for name in (REMATCH_ENTRY, REMATCH_PATCH):
            selected[name] = 'rematch3011/' + name
    for relative in selected.values():
        if sha256(PRODUCTION / relative) != expected[relative]:
            raise ValueError('Recovered macro changed; audit a new producer revision before planning: ' + relative)
    verify_beam_constants(PRODUCTION / selected[PRODUCER], beam)
    hashes = {name: expected[relative] for name, relative in selected.items()}
    if electron_id == REMATCH:
        hashes['ElectronFinder.C'] = sources['rematch3011_reader_sha256']
    return selected, hashes


def build_rematched_reader(selected, expected_hash):
    # Emit the patched bytes to stdout. Neither source in the checkout is modified.
    result = subprocess.run(
        ['patch', '--batch', '--fuzz=0', '--silent', '--output=-', '--reject-file=-',
         str(PRODUCTION / selected[ARCHIVED_READER]), str(PRODUCTION / selected[REMATCH_PATCH])],
        capture_output=True, check=False)
    if result.returncode or hashlib.sha256(result.stdout).hexdigest() != expected_hash:
        raise ValueError('Cannot reproduce the reviewed EICrecon #3011 reader from the archived source and patch')
    return result.stdout


def plan(args):
    if args.files_per_job < 1:
        raise ValueError('--files-per-job must be at least 1')
    if args.femc in ('keep_xpos', 'keep_yneg') and args.electron_id != 'legacy':
        raise ValueError('The preserved producer supports half-plane modes only with legacy electron ID')
    if not math.isfinite(args.eta_split) or args.eta_split <= 0:
        raise ValueError('--eta-split must be positive and finite')
    if args.memory_mb < 1 or args.disk_mb < 1:
        raise ValueError('Resource requests must be positive')
    if not args.environment_label.strip():
        raise ValueError('Provide the exact tested ROOT/podio environment label')
    if args.job_flavour and not re.fullmatch(r'[A-Za-z0-9_-]+', args.job_flavour):
        raise ValueError('Invalid job flavour')
    work = safe_path(args.work_dir)
    if work.exists() or work.is_symlink():
        raise FileExistsError('Use a fresh work directory; existing work is unchanged: ' + str(work))
    validate_work_location(work)
    tuple_dir = args.tuple_output_dir
    if tuple_dir is not None:
        if not tuple_dir.is_absolute():
            raise ValueError('--tuple-output-dir requires an absolute mounted filesystem path')
        tuple_dir = safe_path(tuple_dir)
        if tuple_dir.is_symlink() or tuple_dir.exists():
            raise FileExistsError('Use a fresh tuple output directory: ' + str(tuple_dir))
        tuple_dir = safe_path(tuple_dir.resolve())
        check_beam_label(tuple_dir, args.beam)
        validate_work_location(tuple_dir)
        control_dir = work.resolve()
        if tuple_dir.is_relative_to(control_dir) or control_dir.is_relative_to(tuple_dir):
            raise ValueError('Tuple output and work directories must be separate, non-overlapping paths')
    records = input_records(args.input_list, args.state, args.hash_inputs, beam=args.beam)
    selected_sources, source_hashes = verify_macro_sources(args.beam, args.electron_id)
    rematched_reader = None
    if args.electron_id == REMATCH:
        rematched_reader = build_rematched_reader(selected_sources, source_hashes['ElectronFinder.C'])
    setup = None
    if args.setup_script:
        setup_path = safe_path(args.setup_script).resolve()
        if not setup_path.is_file():
            raise FileNotFoundError(setup_path)
        setup = dict(path=str(setup_path), sha256=sha256(setup_path))
    info = dict(schema_version=1, kind='new_shared_filesystem_condor_plan',
                producer_identity='archived source: ' + selected_sources[PRODUCER],
                beam=args.beam, state=args.state, data_version=args.data_version,
                electron_id=args.electron_id, femc=args.femc, eta_split=args.eta_split,
                memory_mb=args.memory_mb, disk_mb=args.disk_mb, job_flavour=args.job_flavour,
                environment_label=args.environment_label, setup_script=setup,
                input_list=str(Path(args.input_list).resolve()), inputs=records,
                work_dir=str(work), repository_root=str(REPO), source_hashes=source_hashes,
                worker_sha256=sha256(__file__), created_utc=datetime.now(timezone.utc).isoformat(),
                normalization='Historical per-job weights; merge does not normalize. Run framework preparation/audit before yields.',
                job_count=len(records), files_per_job=args.files_per_job,
                condor_job_count=(len(records) + args.files_per_job - 1) // args.files_per_job)
    if tuple_dir is not None:
        info['tuple_output_dir'] = str(tuple_dir)
    if args.electron_id == REMATCH:
        info['accounting_group'] = 'group_u_CMS.CAF.COMM'
        info['producer_identity'] = selected_sources[PRODUCER] + ' with downstream EICrecon #3011 rematching'
        info['source_paths'] = selected_sources
    print(json.dumps(info, indent=2))
    if not args.write:
        print('PLAN ONLY. Add --write to create the fresh work directory. No ROOT or Condor command executed.')
        return 0
    work.mkdir(parents=True, exist_ok=False)
    if tuple_dir is not None:
        try:
            tuple_dir.mkdir(parents=True, exist_ok=False)
        except OSError:
            work.rmdir()  # Only our newly created, still-empty control directory.
            raise
    (work / 'snapshot/macros').mkdir(parents=True)
    for name, relative in selected_sources.items():
        shutil.copyfile(PRODUCTION / relative, work / 'snapshot/macros' / name)
    if rematched_reader is not None:
        (work / 'snapshot/macros/ElectronFinder.C').write_bytes(rematched_reader)
    worker_path = work / 'snapshot/jobs.py'
    shutil.copyfile(__file__, worker_path)
    launcher = work / 'snapshot/run_job.sh'
    # Condor may stage the launcher in the execution slot. Its worker stays on shared storage.
    command = 'worker "$1"' if args.files_per_job == 1 else f'group-worker --batch {shlex.quote(str(work))} --index "$1"'
    launcher.write_text('#!/usr/bin/env bash\nset -euo pipefail\n'
                        f'exec python3 -B {shlex.quote(str(worker_path))} {command} --execute\n')
    launcher.chmod(0o755)
    (work / 'logs').mkdir()
    jobs = []
    for index, record in enumerate(records):
        job_dir = work / 'jobs' / ('job_%06d' % index)
        job_dir.mkdir(parents=True)
        job = {**info, 'inputs': [record], 'job_id': index, 'job_dir': str(job_dir),
               'macro_dir': str(work / 'snapshot/macros'), 'output': str(expected_output_path(info, index))}
        job_path = job_dir / 'job.json'
        save_json(job_path, job)
        jobs.append(str(job_path))
    info['jobs'] = jobs
    save_json(work / 'plan.json', info)
    (work / 'submit.sub').write_text(submit_text(info, jobs))
    print('CREATED:', work)
    print(f'Condor processes: {info["condor_job_count"]} for {len(records)} inputs (up to {args.files_per_job} files per process).')
    print('Inspect a job: python3 -B %s worker %s' % (work / 'snapshot/jobs.py', jobs[0]))
    print('After a successful local smoke job, submit remaining pending jobs explicitly; see docs/production.md.')
    print('Condor description:', work / 'submit.sub')
    return 0


def pending(args):
    """Create a submit file for untouched jobs after a local smoke run."""
    work = Path(args.batch).expanduser().resolve()
    info = json.loads((work / 'plan.json').read_text())
    jobs = []
    for path in info['jobs']:
        job = load_job(path, batch=info)
        directory = Path(job['job_dir'])
        result_path = directory / 'result.json'
        if result_path.exists():
            result = json.loads(result_path.read_text())
            if (result.get('status') != 'DONE' or result.get('output') != job['output']
                    or sha256(job['output']) != result.get('output_sha256')):
                raise ValueError('Failed or changed job requires a fresh attempt; refusing to silently skip it: ' + str(path))
            continue
        if (directory / 'execution.lock').exists() or Path(job['output']).exists() or (directory / 'root.log').exists():
            raise ValueError('A job has started without a validated result; inspect before submitting: ' + str(path))
        jobs.append(path)
    print('Untouched jobs:', len(jobs))
    if not jobs:
        print('No pending jobs; all outputs are DONE.')
        return 0
    if args.write:
        destination = work / 'pending.sub'
        with destination.open('x') as stream:
            stream.write(submit_text(info, jobs))
        print('CREATED:', destination)
        print('Submit explicitly only once, after checking the batch queue: condor_submit ' + shlex.quote(str(destination)))
    else:
        print('PLAN ONLY. Add --write to create pending.sub. This checks files, not the live Condor queue.')
    return 0


class BatchSnapshot(NamedTuple):
    plan: dict
    work: Path


def validate_batch_snapshot(batch, *, worker_hash=None):
    """Check the artifacts shared by all jobs before processing the batch."""
    work = Path(batch['work_dir']).resolve()
    jobs_dir = work / 'jobs'
    if jobs_dir.resolve() != jobs_dir:
        raise ValueError('Job directory differs from its frozen location')
    if 'tuple_output_dir' in batch:
        tuple_dir = Path(batch['tuple_output_dir'])
        if tuple_dir.resolve() != tuple_dir:
            raise ValueError('Tuple output directory differs from its frozen destination')
        if not tuple_dir.is_dir():
            raise FileNotFoundError('Tuple output directory is unavailable: ' + str(tuple_dir))
    macro_dir = work / 'snapshot/macros'
    if macro_dir.resolve() != macro_dir:
        raise ValueError('Macro snapshot path mismatch')
    expected_sources = set(MACROS)
    if batch['electron_id'] == REMATCH:
        expected_sources.update((ARCHIVED_READER, REMATCH_ENTRY, REMATCH_PATCH))
    if set(batch['source_hashes']) != expected_sources:
        raise ValueError('Unexpected source list for electron ID ' + batch['electron_id'])
    for name in expected_sources:
        if sha256(macro_dir / name) != batch['source_hashes'][name]:
            raise ValueError('Macro snapshot hash mismatch: ' + name)
    verify_beam_constants(macro_dir / PRODUCER, batch['beam'])
    if worker_hash is None:
        worker_hash = sha256(work / 'snapshot/jobs.py')
    if worker_hash != batch['worker_sha256']:
        raise ValueError('Worker differs from the recorded plan')
    setup = batch.get('setup_script')
    if setup and sha256(setup['path']) != setup['sha256']:
        raise ValueError('Environment setup script changed after planning')
    return BatchSnapshot(batch, work)


def load_job(path, *, batch=None, worker_hash=None, snapshot=None):
    path = Path(path)
    # The snapshot already resolved the common parents. Only the two job-specific
    # components can add a symlink for canonical manifests written by the planner.
    if (snapshot is None or path.parent.parent != snapshot.work / 'jobs'
            or path.parent.is_symlink() or path.is_symlink()):
        path = path.resolve()
    job = json.loads(path.read_text())
    if job.get('schema_version') != 1 or job.get('kind') != 'new_shared_filesystem_condor_plan':
        raise ValueError('Unrecognized job format')
    if job['beam'] not in PRODUCERS:
        raise ValueError('Unknown production beam: ' + str(job['beam']))
    if snapshot is None:
        if batch is None:
            batch = json.loads((Path(job['work_dir']) / 'plan.json').read_text())
        if worker_hash is None:
            worker_hash = sha256(__file__)
        snapshot = validate_batch_snapshot(batch, worker_hash=worker_hash)
    else:
        if batch is not None and batch is not snapshot.plan:
            raise ValueError('Snapshot belongs to another batch')
        batch = snapshot.plan
    work = snapshot.work
    index = job['job_id']
    if not isinstance(index, int) or index < 0 or index >= len(batch['jobs']):
        raise ValueError('Invalid job index in batch')
    declared_path = Path(batch['jobs'][index])
    if ((declared_path != path and declared_path.resolve() != path)
            or batch['job_count'] != len(batch['jobs'])):
        raise ValueError('Job does not belong to the declared complete batch')
    for key, value in batch.items():
        if key not in ('inputs', 'jobs') and job.get(key) != value:
            raise ValueError('Job differs from its batch configuration: ' + key)
    if job['inputs'] != [batch['inputs'][index]]:
        raise ValueError('Job input differs from its batch manifest')
    expected = work / 'jobs' / ('job_%06d' % job['job_id'])
    directory = Path(job['job_dir'])
    if path != expected / 'job.json' or (directory != expected and directory.resolve() != expected):
        raise ValueError('Job paths do not match the declared work directory')
    output = Path(job['output'])
    if 'tuple_output_dir' in batch:
        planned_output = expected_output_path(batch, index)
        if output != planned_output or output.is_symlink():
            raise ValueError('Tuple output path differs from its frozen destination')
    else:
        planned_output = expected / 'output.root'
        redirected = output.is_symlink() if output == planned_output else output.resolve() != planned_output
        if redirected or 'tuple_output_dir' in job:
            raise ValueError('Output must be this job output.root')
    macro_dir = work / 'snapshot/macros'
    if Path(job['macro_dir']) != macro_dir and Path(job['macro_dir']).resolve() != macro_dir:
        raise ValueError('Macro snapshot path mismatch')
    record = job['inputs'][0]
    stat = Path(record['path']).stat()
    if (stat.st_size, stat.st_mtime_ns) != (record['size_bytes'], record['mtime_ns']):
        raise ValueError('Input changed after planning')
    if record.get('sha256') and sha256(record['path']) != record['sha256']:
        raise ValueError('Input SHA256 mismatch')
    return job


def root_command(job):
    quote = lambda text: json.dumps(str(text), ensure_ascii=True)
    args = [quote(job['inputs'][0]['path']), quote(job['output']), quote(job['data_version']),
            quote(job['femc']), format(job['eta_split'], '.17g')]
    entry = PRODUCER
    if job['electron_id'] == REMATCH:
        entry = REMATCH_ENTRY
    else:
        args.append(quote(job['electron_id']))
    call = str(Path(job['macro_dir']) / entry) + '(' + ','.join(args) + ')'
    return ['root', '-l', '-b', '-q', call]


def run_with_setup(command, setup, **kwargs):
    if setup:
        return subprocess.run(['bash', '-c', 'set -euo pipefail; source "$1"; shift; exec "$@"',
                               'upsilon-env', setup['path'], *command], **kwargs)
    return subprocess.run(command, **kwargs)


def root_counts(path, tree, **options):
    import uproot
    if hasattr(path, 'read'):
        options.setdefault('handler', uproot.source.object.ObjectSource)
        options.setdefault('use_threads', False)
    with uproot.open(path, **options) as source:
        n = int(source[tree].num_entries)
        if n <= 0:
            raise ValueError('Empty event tree in ' + str(path))
        if tree == 'UpsilonTree':
            values = source['h_EventCount'].values(flow=True)
            if any(not math.isfinite(float(x)) or x < 0 for x in values) or values[0] != 0 or values[-1] != 0:
                raise ValueError('Malformed generated-event count histogram')
            if not math.isclose(float(sum(values)), n, rel_tol=0., abs_tol=1e-6):
                raise ValueError('h_EventCount differs from all tuple rows')
        return n


def worker(args, *, batch=None, resume=False):
    job = load_job(args.job, batch=batch)
    command = root_command(job)
    print('ROOT command:', shlex.join(command))
    print('Environment:', job['environment_label'])
    if not args.execute:
        print('PLAN ONLY. Add --execute to process this one job. Existing files remain unchanged.')
        return 0
    job_dir = Path(job['job_dir'])
    output = Path(job['output'])
    if resume:
        if (job_dir / 'execution.lock').exists():
            raise FileExistsError(f'File job is locked: {job_dir}')
        if (job_dir / 'result.json').exists():
            verified_output_rows(job)
            print('ALREADY DONE:', output, flush=True)
            return 0
        if (job_dir / 'root.log').exists():
            raise FileExistsError(f'Partial file job has a log without a receipt: {job_dir}')
    if output.exists() or output.is_symlink() or (job_dir / 'result.json').exists():
        raise FileExistsError('Existing output/result is preserved; plan a fresh attempt directory')
    lock = job_dir / 'execution.lock'
    fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.close(fd)
    result = dict(schema_version=1, status='FAILED', job=str(Path(args.job).resolve()),
                  beam=job['beam'], state=job['state'], electron_id=job['electron_id'], source_hashes=job['source_hashes'],
                  command=command, environment_label=job['environment_label'],
                  normalization=job['normalization'], started_utc=datetime.now(timezone.utc).isoformat())
    stage = 'input_read'
    try:
        n_input = root_counts(job['inputs'][0]['path'], 'events')
        stage = 'tuple_production'
        # Reserve our output exclusively. ROOT RECREATE writes this owned placeholder.
        out_fd = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        os.close(out_fd)
        with (job_dir / 'root.log').open('x') as log:
            process = run_with_setup(command, job.get('setup_script'), cwd=job_dir, stdout=log, stderr=subprocess.STDOUT, check=False)
        result['root_returncode'] = process.returncode
        if process.returncode != 0:
            raise RuntimeError('ROOT failed; inspect root.log; partial output retained')
        n_output = root_counts(output, 'UpsilonTree')
        if n_input != n_output:
            raise ValueError('Generated rows were lost: input events count differs from tuple count')
        result.update(status='DONE', input_rows=n_input, output_rows=n_output,
                      output=str(output), output_sha256=sha256(output))
    except Exception as exc:
        result.update(error=str(exc), failure_stage=stage, failure_type=type(exc).__name__)
        if getattr(exc, 'errno', None) is not None:
            result['failure_errno'] = exc.errno
        raise
    finally:
        result['finished_utc'] = datetime.now(timezone.utc).isoformat()
        try:
            save_json(job_dir / 'result.json', result)
        finally:
            lock.unlink()
    print('DONE:', output)
    return 0


def group_worker(args):
    work = Path(args.batch).resolve()
    batch = json.loads((work / 'plan.json').read_text())
    width = batch.get('files_per_job', 1)
    count = len(batch['jobs'])
    if (Path(batch['work_dir']).resolve() != work or type(width) is not int or width < 1
            or count != batch['job_count'] or count != len(batch['inputs'])):
        raise ValueError('Invalid grouped batch manifest')
    groups = (count + width - 1) // width
    if groups != batch.get('condor_job_count', count) or not 0 <= args.index < groups:
        raise ValueError('Invalid Condor group index or count')
    start = args.index * width
    members = batch['jobs'][start:start + width]
    print(f'Condor group {args.index}: {len(members)} file jobs', flush=True)
    failures = 0
    for path in members:
        try:
            worker(argparse.Namespace(job=path, execute=args.execute), batch=batch, resume=True)
        except (OSError, ValueError, RuntimeError, ImportError, KeyError) as exc:
            failures += 1
            print(f'FAILED {path}: {exc}', file=sys.stderr, flush=True)
    if failures:
        raise RuntimeError(f'Condor group {args.index}: {failures} of {len(members)} file jobs failed')
    return 0


def verified_output_rows(job, result=None):
    if result is None:
        result = json.loads((Path(job['job_dir']) / 'result.json').read_text())
    path = Path(job['output'])
    if result.get('status') != 'DONE' or result.get('output') != str(path):
        detail = result.get('error') or f"recorded output {result.get('output')!r}, expected {path}"
        raise ValueError(f"Input {job['job_id']:06d} has no validated DONE output "
                         f"(status={result.get('status')!r}): {detail}")
    for key in ('beam', 'state', 'source_hashes'):
        if result.get(key) != job[key]:
            raise ValueError('Job result differs from the planned ' + key)
    if job['electron_id'] == REMATCH and result.get('electron_id') != REMATCH:
        raise ValueError('Corrected batch requires confirmed rematching in every job result')
    rows = result['output_rows']
    with path.open('rb') as stream:
        data = stream.read(IN_MEMORY_TUPLE_BYTES + 1)
        digest = hashlib.sha256(data)
        if len(data) <= IN_MEMORY_TUPLE_BYTES:
            # Typical intermediate tuples fit here: inspect the exact bytes already read from EOS.
            source = io.BytesIO(data)
        else:
            for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                digest.update(chunk)
            stream.seek(0)
            source = stream
        if digest.hexdigest() != result['output_sha256']:
            raise ValueError('Output changed after job validation: ' + str(path))
        if result.get('input_rows') != rows or root_counts(source, 'UpsilonTree') != rows:
            raise ValueError('Output count changed or differs from input rows: ' + str(path))
    return rows


def is_short_read(message):
    match = re.search(r'expected Chunk of length (\d+),\s*received (\d+) bytes', message)
    return bool(match and int(match[2]) < int(match[1]))


def confirm_input_damage(record):
    """Read a separate local copy. Source access and interrupted copies must fail the merge."""
    import uproot

    src = Path(record['path'])
    with tempfile.TemporaryDirectory(prefix='upsilon-merge-check-') as tmp:
        copy = Path(tmp) / 'input.root'
        shutil.copyfile(src, copy)
        size = copy.stat().st_size
        stat = src.stat()
        if size != record['size_bytes'] or (stat.st_size, stat.st_mtime_ns) != (size, record['mtime_ns']):
            raise ValueError(f'Cannot confirm a complete copy of the planned input: {src}')
        digest = sha256(copy)
        if record.get('sha256') and digest != record['sha256']:
            raise ValueError(f'Diagnostic copy differs from the planned input checksum: {src}')
        with copy.open('rb') as stream:
            header = stream.read(20)
            if len(header) < 20 or header[:4] != b'root':
                reason = 'Missing or incomplete ROOT header'
            else:
                width = 8 if int.from_bytes(header[4:8], 'big') >= 1_000_000 else 4
                end = int.from_bytes(header[12:12 + width], 'big')
                stream.seek(0)
                try:
                    root_counts(stream, 'events', handler=uproot.source.object.ObjectSource)
                except (OSError, EOFError, ValueError, zlib.error) as exc:
                    if getattr(exc, 'errno', None) is not None:
                        raise
                    if end > size:
                        reason = f'ROOT header declares {end} bytes; file has {size}'
                    elif is_short_read(str(exc)) or isinstance(exc, EOFError):
                        reason = 'Unexpected end of file reading the events tree from a verified local copy'
                    elif isinstance(exc, zlib.error) and str(exc).startswith('Error -3 '):
                        reason = f'Corrupt compressed ROOT data in the local copy: {exc}'
                    else:
                        raise ValueError(f'Cannot confirm input corruption in {src}: {exc}') from exc
                else:
                    raise ValueError(f'Original input now has a readable events tree; retry this job: {src}')
    return dict(confirmation=reason, size_bytes=size, sha256=digest)


def excluded_input(job, result):
    folder = Path(job['job_dir'])
    record = job['inputs'][0]
    error = result.get('error')
    identity = ('beam', 'state', 'electron_id', 'source_hashes', 'environment_label', 'normalization')
    if (not isinstance(error, str) or not error or result.get('job') != str(folder / 'job.json')
            or any(result.get(key) != job[key] for key in identity) or result.get('command') != root_command(job)):
        raise ValueError(f'Cannot verify failed worker receipt: {folder}')
    stage = result.get('failure_stage')
    if stage is None:
        # Old workers did not label the read stage. Accept only uproot errors naming this exact input.
        src = record['path']
        named_input = any(line.strip() in (f'for file path {src}', f'in file {src}') for line in error.splitlines())
        if named_input and (is_short_read(error) or 'not a ROOT file' in error):
            stage = 'input_read'
    evidence = (folder / 'execution.lock', folder / 'reuse-pending.json', folder / 'root.log', Path(job['output']))
    if (stage != 'input_read' or result.get('failure_errno') is not None
            or result.get('failure_type') in ('ImportError', 'ModuleNotFoundError')
            or any(os.path.lexists(p) for p in evidence)
            or any(key in result for key in ('output', 'input_rows', 'root_returncode'))):
        raise ValueError(f'Only a failed original-input read can be excluded: {folder}: {error}')
    confirmation = confirm_input_damage(record)
    return dict(job_id=job['job_id'], input=record['path'], error=error, **confirmation)


def check_merge_input(index, snapshot, skip_unreadable):
    job = load_job(snapshot.plan['jobs'][index], snapshot=snapshot)
    if job['job_id'] != index:
        raise ValueError('Unexpected job order in merge plan')
    result = json.loads((Path(job['job_dir']) / 'result.json').read_text())
    if result.get('status') == 'FAILED' and skip_unreadable:
        return None, 0, excluded_input(job, result)
    return Path(job['output']), verified_output_rows(job, result), None


def _init_merge_process(snapshot, skip_unreadable):
    # Transfer the full plan once per process. Individual tasks carry just a file index.
    global _merge_context
    _merge_context = (snapshot, skip_unreadable)


def _check_in_merge_process(index):
    return check_merge_input(index, *_merge_context)


def merge(args):
    started = time.monotonic()
    check_workers = getattr(args, 'check_workers', DEFAULT_CHECK_WORKERS)
    hadd_jobs = getattr(args, 'hadd_jobs', 1)
    if min(check_workers, hadd_jobs) < 1:
        raise ValueError('--check-workers and --hadd-jobs must both be positive')
    plan_path = Path(args.batch).expanduser().resolve() / 'plan.json'
    plan_bytes = plan_path.read_bytes()
    plan_hash = hashlib.sha256(plan_bytes).hexdigest()
    plan_data = json.loads(plan_bytes)
    if Path(plan_data['work_dir']).resolve() != plan_path.parent:
        raise ValueError('Merge plan belongs to a different batch')
    out = safe_path(args.output)
    check_beam_label(out, plan_data['beam'])
    validate_work_location(out, plan_data['repository_root'])
    if out.exists() or out.is_symlink():
        raise FileExistsError('Merged output already exists; it will not be replaced')
    if not out.parent.is_dir():
        raise FileNotFoundError('Create the intended output directory first: ' + str(out.parent))
    input_list = out.with_name(out.name + '.inputs.txt')
    result_path = out.with_name(out.name + '.merge.json')
    lock = out.with_name(out.name + '.merge.lock')
    if args.execute and any(os.path.lexists(p) for p in (input_list, result_path, lock)):
        raise FileExistsError('Existing merge records or lock: inspect the previous attempt or use a fresh output name')
    hadd_temp = Path(getattr(args, 'hadd_temp_dir', None) or out.parent).expanduser().resolve()
    if not hadd_temp.is_dir():
        raise FileNotFoundError(f'Missing hadd temporary directory: {hadd_temp}')
    snapshot = validate_batch_snapshot(plan_data)
    paths, excluded, expected = [], [], 0
    total = plan_data['job_count']
    if len(plan_data['inputs']) != total or len(plan_data['jobs']) != total:
        raise ValueError('Merge plan inputs and jobs differ in length')
    # Small pilots finish faster in the current process than starting additional Python interpreters.
    check_workers = min(check_workers, max(1, total // 32))
    skip_unreadable = getattr(args, 'skip_unreadable_inputs', False)
    if check_workers == 1:
        checks = (check_merge_input(i, snapshot, skip_unreadable) for i in range(total))
    else:
        checks = ordered_checks(_check_in_merge_process, range(total), check_workers,
                                executor=ProcessPoolExecutor, mp_context=multiprocessing.get_context('spawn'),
                                initializer=_init_merge_process, initargs=(snapshot, skip_unreadable))
    print(f"Checking {total:,} planned inputs with {check_workers} validation process(es).", flush=True)
    last_progress = started
    for count, (path, rows, skipped) in enumerate(checks, 1):
        if skipped is not None:
            excluded.append(skipped)
            print(f"EXCLUDED {skipped['job_id']:06d}: {skipped['input']} ({skipped['confirmation']})", flush=True)
        else:
            paths.append(path)
            expected += rows
        now = time.monotonic()
        if count % 1000 == 0 or now - last_progress >= 15:
            print(f"Checked {count:,}/{total:,} inputs in {now - started:.1f} s; {len(excluded)} excluded.", flush=True)
            last_progress = now
    if not paths or len(set(paths)) != len(paths):
        raise ValueError('Empty or duplicate merge inputs')
    if sha256(plan_path) != plan_hash:
        raise ValueError('Batch plan changed during merge validation')
    validate_batch_snapshot(plan_data)
    # For a small list, launching a subprocess per few inputs costs more than the merge.
    processes = min(hadd_jobs, max(1, len(paths) // 3))
    parallel_flags = ['-j', str(processes), '-d', str(hadd_temp)] if processes > 1 else []
    command = ['hadd', *parallel_flags, str(out), '@' + str(input_list)]
    validation_seconds = time.monotonic() - started
    print(f"Validated {len(paths):,} tuples in {validation_seconds:.1f} s.", flush=True)
    print('Merge command:', shlex.join(command), flush=True)
    print('Expected all-generated rows:', expected)
    if excluded:
        print(f'Merging {len(paths)} of {plan_data["job_count"]} inputs; {len(excluded)} confirmed unreadable.')
    if not args.execute:
        print('PLAN ONLY. Add --execute to merge. Weights and histograms are summed, not renormalized.')
        return 0
    # hadd without -f refuses an existing output; never overwrite an older file.
    if any(os.path.lexists(p) for p in (out, result_path, input_list)):
        raise FileExistsError('Existing merge provenance is preserved; use a fresh output name')
    lock_fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    os.close(lock_fd)
    result = dict(status='FAILED', command=command, expected_rows=expected, inputs=list(map(str, paths)),
                  beam=plan_data['beam'], state=plan_data['state'], source_hashes=plan_data['source_hashes'],
                  electron_id=plan_data['electron_id'],
                  normalization='Unchanged raw per-job weights. Normalize the complete merged sample in the framework.',
                  plan_sha256=plan_hash, merge_driver_sha256=sha256(__file__),
                  planned_input_count=plan_data['job_count'], included_input_count=len(paths), excluded_inputs=excluded,
                  check_workers=check_workers, hadd_jobs=processes, validation_seconds=validation_seconds)
    try:
        with input_list.open('x') as stream:
            stream.write('\n'.join(map(str, paths)) + '\n')
        result['input_list_sha256'] = sha256(input_list)
        hadd_started = time.monotonic()
        process = run_with_setup(command, plan_data.get('setup_script'), check=False)
        result['hadd_seconds'] = time.monotonic() - hadd_started
        if process.returncode != 0:
            raise RuntimeError('hadd failed; any partial output remains for inspection')
        print(f"hadd finished in {result['hadd_seconds']:.1f} s. Checking merged output.", flush=True)
        final_check_started = time.monotonic()
        count = root_counts(out, 'UpsilonTree')
        if count != expected:
            raise ValueError('Merged tuple/generated count differs from the sum of all declared jobs')
        result.update(status='DONE_WITH_EXCLUSIONS' if excluded else 'DONE',
                      output=str(out), output_rows=count, sha256=sha256(out))
        result['final_check_seconds'] = time.monotonic() - final_check_started
    except Exception as exc:
        result['error'] = str(exc)
        raise
    finally:
        try:
            save_json(result_path, result)
        finally:
            lock.unlink()
    print(result['status'] + ':', out)
    print('Next: framework normalization and input audit; this ROOT file still has raw per-job weights.')
    return 0


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('plan', help='Inspect a production plan; --write creates job/Condor files')
    a.add_argument('--input-list', type=Path, required=True)
    a.add_argument('--work-dir', type=Path, required=True)
    a.add_argument('--tuple-output-dir', type=Path,
                   help='Fresh absolute directory on a mounted filesystem for per-input ROOT files; '
                        'Condor logs and job records stay under --work-dir')
    a.add_argument('--beam', choices=sorted(PRODUCERS), required=True)
    a.add_argument('--state', choices=('1S', '2S', '3S'), required=True)
    a.add_argument('--electron-id', choices=('legacy', 'reco', REMATCH), required=True)
    a.add_argument('--femc', choices=SCENARIOS, default='full')
    a.add_argument('--eta-split', type=float, default=1.68)
    a.add_argument('--data-version', choices=('old', 'new'), default='new')
    a.add_argument('--environment-label', required=True)
    a.add_argument('--setup-script', type=Path)
    a.add_argument('--memory-mb', type=int, default=4096)
    a.add_argument('--disk-mb', type=int, default=8192)
    a.add_argument('--job-flavour', default=None, help='Optional CERN JobFlavour, explicitly supplied')
    a.add_argument('--files-per-job', type=int, default=1, help='Maximum number of sequential input files per Condor process')
    a.add_argument('--hash-inputs', action='store_true')
    a.add_argument('--write', action='store_true')
    a = sub.add_parser('worker', help='Inspect one frozen job; --execute invokes ROOT')
    a.add_argument('job', type=Path)
    a.add_argument('--execute', action='store_true')
    a = sub.add_parser('group-worker', help='Run a group of frozen per-file workers in original input order')
    a.add_argument('--batch', type=Path, required=True)
    a.add_argument('--index', type=int, required=True, help='Zero-based group index')
    a.add_argument('--execute', action='store_true')
    a = sub.add_parser('pending', help='List untouched jobs; --write creates pending.sub without submission')
    a.add_argument('--batch', type=Path, required=True)
    a.add_argument('--write', action='store_true')
    a = sub.add_parser('merge', help='Verify a complete batch; --execute invokes hadd without -f')
    a.add_argument('--batch', type=Path, required=True)
    a.add_argument('--output', type=Path, required=True)
    add_check_workers(a)
    a.add_argument('--hadd-jobs', type=int, default=DEFAULT_CHECK_WORKERS,
                   help='Maximum ROOT merge processes (default: %(default)s; 1 for serial hadd)')
    a.add_argument('--hadd-temp-dir', type=Path,
                   help='Existing directory for temporary parallel hadd files (default: final output directory)')
    a.add_argument('--skip-unreadable-inputs', action='store_true',
                   help='Exclude only failed original-input reads with confirmed ROOT file damage; record each exclusion')
    a.add_argument('--execute', action='store_true')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    return {'plan': plan, 'worker': worker, 'group-worker': group_worker,
            'pending': pending, 'merge': merge}[args.command](args)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, ImportError, RuntimeError, KeyError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        sys.exit(2)
