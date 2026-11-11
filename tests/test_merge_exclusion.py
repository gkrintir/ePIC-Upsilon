"""Only verified, unreadable original ROOT files may be omitted from a merge."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'production'))
import jobs


def short_read(path):
    return f"expected Chunk of length 561,\nreceived 0 bytes from FSSpecSource\nfor file path {path}"


@pytest.fixture
def campaign(tmp_path, monkeypatch, request):
    uproot = pytest.importorskip('uproot')
    import numpy as np

    def write_tree(path, tree, rows):
        with uproot.recreate(path) as output:
            output.mktree(tree, {'event': 'int32'})
            output[tree].extend({'event': np.arange(rows, dtype=np.int32)})
            if tree == 'UpsilonTree':
                output['h_EventCount'] = np.array([float(rows)]), np.array([0., 1.])

    # Include another state/beam to exercise the generic merge entry point.
    dataset = tmp_path / 'Upsilon2S/9x130'
    dataset.mkdir(parents=True)
    sources = [dataset / f'file{i}.root' for i in range(4)]
    broken = getattr(request, 'param', (1, 3))
    for i, path in enumerate(sources):
        write_tree(path, 'events', i + 1)
        if i in broken:
            with path.open('r+b') as stream:
                stream.truncate(1000)  # Retain a genuine ROOT header declaring its original EOF.
    inventory = tmp_path / 'inputs.txt'
    inventory.write_text('\n'.join(map(str, sources)) + '\n')
    batch = tmp_path / 'controls'
    tuples = tmp_path / 'tuples'
    revision = tmp_path / 'earlier-jobs.py'
    revision.write_bytes(Path(jobs.__file__).read_bytes() + b'\n# earlier frozen worker\n')
    with monkeypatch.context() as worker:
        worker.setattr(jobs, '__file__', str(revision))
        jobs.main(['plan', '--input-list', str(inventory), '--work-dir', str(batch),
                   '--tuple-output-dir', str(tuples), '--state', '2S', '--beam', '9x130',
                   '--electron-id', 'reco', '--environment-label', 'test-only', '--write'])
        real_count = jobs.root_counts

        def count(path, tree):
            if tree == 'events' and Path(path) in (sources[i] for i in broken):
                raise OSError(short_read(path))
            return real_count(path, tree)

        def produce(command, setup, *, cwd, **kwargs):
            job = json.loads((Path(cwd) / 'job.json').read_text())
            write_tree(job['output'], 'UpsilonTree', job['job_id'] + 1)
            return SimpleNamespace(returncode=0)

        worker.setattr(jobs, 'root_counts', count)
        worker.setattr(jobs, 'run_with_setup', produce)
        for i in range(4):
            job = batch / 'jobs' / f'job_{i:06d}' / 'job.json'
            if i in broken:
                with pytest.raises(OSError, match='expected Chunk'):
                    jobs.main(['worker', str(job), '--execute'])
                # Existing campaigns have no failure-stage or input-diagnostic fields.
                result = job.with_name('result.json')
                record = json.loads(result.read_text())
                legacy_fields = ('schema_version', 'status', 'job', 'beam', 'state', 'electron_id',
                                 'source_hashes', 'command', 'environment_label', 'normalization',
                                 'started_utc', 'error', 'finished_utc')
                result.write_text(json.dumps({key: record[key] for key in legacy_fields}))
            else:
                jobs.main(['worker', str(job), '--execute'])

    def hadd(command, setup, **kwargs):
        assert command[0] == 'hadd'
        inputs = [Path(path) for path in Path(command[2][1:]).read_text().splitlines()]
        assert inputs == [tuples / f'job_{i:06d}.root' for i in range(4) if i not in broken]
        write_tree(command[1], 'UpsilonTree', sum(jobs.root_counts(p, 'UpsilonTree') for p in inputs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(jobs, 'run_with_setup', hadd)
    return batch, sources, tmp_path / 'merged.root'


def merge(batch, output, *args):
    return jobs.main(['merge', '--batch', str(batch), '--output', str(output), *args])


def test_opt_in_records_each_exclusion_and_retains_frozen_receipts(campaign):
    batch, sources, output = campaign
    evidence = [batch / 'plan.json', batch / 'snapshot/jobs.py', *batch.glob('jobs/*/result.json')]
    before = {p: p.read_bytes() for p in evidence}
    with pytest.raises(ValueError):
        merge(batch, output, '--execute')
    merge(batch, output, '--skip-unreadable-inputs')
    assert not output.exists()
    merge(batch, output, '--skip-unreadable-inputs', '--execute')
    receipt = json.loads(output.with_suffix('.root.merge.json').read_text())
    assert receipt['status'] == 'DONE_WITH_EXCLUSIONS'
    assert (receipt['included_input_count'], receipt['planned_input_count'], receipt['output_rows']) == (2, 4, 4)
    assert [(r['job_id'], r['input'], r['error']) for r in receipt['excluded_inputs']] == [
        (i, str(sources[i]), short_read(sources[i])) for i in (1, 3)]
    assert all(p.read_bytes() == content for p, content in before.items())
    with pytest.raises(FileExistsError):
        merge(batch, output, '--skip-unreadable-inputs', '--execute')


@pytest.mark.parametrize('campaign', [()], indirect=True)
def test_opt_in_keeps_complete_merge_status(campaign):
    batch, _, output = campaign
    merge(batch, output, '--skip-unreadable-inputs', '--execute')
    receipt = json.loads(output.with_suffix('.root.merge.json').read_text())
    assert (receipt['status'], receipt['output_rows']) == ('DONE', 10)


@pytest.mark.parametrize('campaign', [(0, 1, 2, 3)], indirect=True)
def test_all_unreadable_refuses_empty_merge(campaign):
    batch, _, output = campaign
    with pytest.raises(ValueError):
        merge(batch, output, '--skip-unreadable-inputs', '--execute')
    assert not output.exists()


@pytest.mark.parametrize('problem', ['unfinished', 'partial', 'permission', 'wrong_path',
                                     'producer_error', 'readable', 'changed_done_tuple'])
def test_opt_in_still_blocks_other_problems(campaign, problem):
    batch, sources, output = campaign
    folder = batch / 'jobs/job_000001'
    receipt = folder / 'result.json'
    if problem == 'unfinished':
        receipt.unlink()
    elif problem == 'partial':
        (folder / 'execution.lock').touch()
    elif problem == 'changed_done_tuple':
        (output.parent / 'tuples/job_000000.root').write_bytes(b'changed tuple')
    else:
        if problem == 'readable':
            # This original is readable. A previous short read alone does not justify omission.
            folder = batch / 'jobs/job_000000'
            receipt = folder / 'result.json'
            record = json.loads((batch / 'jobs/job_000001/result.json').read_text())
            spec = json.loads((folder / 'job.json').read_text())
            record.update(job=str(folder / 'job.json'), command=jobs.root_command(spec), error=short_read(sources[0]))
            (folder / 'root.log').unlink()
            Path(spec['output']).unlink()
        else:
            record = json.loads(receipt.read_text())
            if problem == 'permission':
                record['error'] = f"[Errno 13] Permission denied: '{sources[1]}'"
            elif problem == 'wrong_path':
                record['error'] = short_read(str(sources[1]) + '.other')
            elif problem == 'producer_error':
                record['root_returncode'] = 1
                (folder / 'root.log').touch()
        receipt.write_text(json.dumps(record))
    with pytest.raises((ValueError, OSError, RuntimeError)):
        merge(batch, output, '--skip-unreadable-inputs', '--execute')
    assert not output.exists()
