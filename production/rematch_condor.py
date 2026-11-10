"""Locate the grouped campaign's recorded submissions and query their live scheduler state."""
import json
import os
from pathlib import Path
import pwd
import re
import subprocess


def default_batch():
    workspaces = Path(__file__).resolve().parents[1] / 'workspaces'
    afs = workspaces / 'condor-1S-9x275-full-grouped-v1' / 'batch'
    eos = workspaces / 'condor-1S-9x275-full-grouped-eos-v1' / 'batch'
    if (eos.parent / 'submission.txt').exists() or (eos / 'after-removal.sub').exists():
        return eos
    if (afs / 'plan.json').is_file() and (
            (afs.parent / 'submission.txt').exists() or not (eos / 'plan.json').is_file()):
        return afs
    return eos


def recorded_clusters(batch):
    """Read grouped submissions only; an ambiguous scheduler acknowledgement needs investigation."""
    batch = Path(batch).resolve()
    initial = batch.parent / 'submission.txt'
    records = ([initial] if initial.exists() else []) + sorted(
        (batch / 'retries').glob('attempt-*/submission.txt'))
    clusters = set()
    for record in records:
        ids = re.findall(r'\b[1-9][0-9]*\s+job\(s\) submitted to cluster ([1-9][0-9]*)\.',
                         record.read_text(), re.IGNORECASE)
        if not ids:
            raise ValueError(f'No confirmed Condor cluster ID in {record}. Check the scheduler and submission output before retrying.')
        clusters.update(map(int, ids))
    return sorted(clusters)


def queue_rows(batch, clusters):
    """Get a narrow view of the current user's live jobs on the original submit schedd."""
    batch = Path(batch).resolve()
    owner = pwd.getpwuid(os.getuid()).pw_name
    jobs_prefix = str(batch / 'jobs') + '/'
    matches = [f'ClusterId == {int(cid)}' for cid in clusters]
    matches.extend([f'Iwd =?= {json.dumps(str(batch))}',
                    f'substr(Iwd, 0, {len(jobs_prefix)}) =?= {json.dumps(jobs_prefix)}',
                    f'Cmd =?= {json.dumps(str(batch / "snapshot" / "run_job.sh"))}'])
    constraint = f'Owner =?= {json.dumps(owner)} && (' + ' || '.join(matches) + ')'
    command = ['condor_q', '-nobatch', '-constraint', constraint, '-json',
               '-attributes', 'ClusterId,ProcId,JobStatus,HoldReason,Iwd,Cmd,Arguments,Args']
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f'Cannot query Condor on this host: {exc}') from exc
    if result.returncode:
        raise RuntimeError(f'Could not query the current Condor scheduler: {result.stderr.strip() or result.returncode}')
    # Some condor_q releases print no JSON bytes for an empty selection.
    if not result.stdout.strip() and result.stderr.strip():
        raise RuntimeError(f'Could not confirm an empty queue: {result.stderr.strip()}')
    rows = json.loads(result.stdout.strip() or '[]')
    if not isinstance(rows, list) or any(not isinstance(r, dict) or
            any(type(r.get(k)) is not int for k in ('ClusterId', 'ProcId', 'JobStatus')) for r in rows):
        raise ValueError('Unexpected condor_q output; the queue state is unknown.')
    return rows
