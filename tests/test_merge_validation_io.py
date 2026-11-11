"""Checks retained when the merge shares snapshot reads and buffers small tuples."""
import io
import json
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'production'))
import jobs


def test_snapshot_checked_once_and_each_manifest_and_input_still_checked(tmp_path):
    original = tmp_path / 'Upsilon1S.root'
    original.write_bytes(b'unit test source')
    inventory = tmp_path / 'inventory.txt'
    inventory.write_text(f'{original}\n')
    setup = tmp_path / 'setup.sh'
    setup.write_text('true\n')
    batch = tmp_path / 'batch'
    jobs.main(['plan', '--input-list', str(inventory), '--work-dir', str(batch), '--beam', '9x275',
               '--state', '1S', '--electron-id', 'reco', '--environment-label', 'test',
               '--setup-script', str(setup), '--hash-inputs', '--write'])
    plan = json.loads((batch / 'plan.json').read_text())
    manifest = Path(plan['jobs'][0])
    with patch.object(jobs, 'sha256', wraps=jobs.sha256) as sha:
        snapshot = jobs.validate_batch_snapshot(plan)
        one_time = [Path(c.args[0]) for c in sha.call_args_list]
        assert one_time.count(setup) == 1
        assert one_time.count(batch / 'snapshot/jobs.py') == 1
        assert set(one_time) == {setup, batch / 'snapshot/jobs.py',
                                 *(batch / 'snapshot/macros' / name for name in plan['source_hashes'])}
        sha.reset_mock()
        with patch.object(Path, 'resolve', side_effect=AssertionError('repeated ancestor traversal')):
            job = jobs.load_job(manifest, snapshot=snapshot)
        assert [Path(c.args[0]) for c in sha.call_args_list] == [original]

    output = Path(job['output'])
    output.symlink_to(original)
    with pytest.raises(ValueError, match='Output must be'):
        jobs.load_job(manifest, snapshot=snapshot)
    output.unlink()
    job['source_hashes']['ElectronFinder.C'] = '0' * 64
    manifest.write_text(json.dumps(job))
    with pytest.raises(ValueError, match='batch configuration: source_hashes'):
        jobs.load_job(manifest, snapshot=snapshot)
    setup.write_text('false\n')
    with pytest.raises(ValueError, match='Environment setup script changed'):
        jobs.validate_batch_snapshot(plan)


@pytest.mark.parametrize('buffer_bytes', [4 << 20, 16])
def test_buffer_and_stream_verify_real_root_contents_in_one_open(tmp_path, monkeypatch, buffer_bytes):
    uproot = pytest.importorskip('uproot')
    import numpy as np

    path = tmp_path / 'tuple.root'
    with uproot.recreate(path) as f:
        f.mktree('UpsilonTree', {'i': 'int64'}).extend({'i': np.arange(3)})
        f['h_EventCount'] = (np.array([3.]), np.array([0., 1.]))
    job = dict(output=str(path), beam='9x130', state='2S', electron_id='reco', source_hashes={})
    receipt = dict(job, status='DONE', output_rows=3, input_rows=3, output_sha256=jobs.sha256(path))
    monkeypatch.setattr(jobs, 'IN_MEMORY_TUPLE_BYTES', buffer_bytes)
    open_file = Path.open
    opens = []

    def count_open(candidate, *args, **kwargs):
        if candidate == path:
            opens.append(candidate)
        return open_file(candidate, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', count_open)
    with patch.object(jobs, 'root_counts', wraps=jobs.root_counts) as counts:
        assert jobs.verified_output_rows(job, receipt) == 3
        assert isinstance(counts.call_args.args[0], io.BytesIO) == (path.stat().st_size <= buffer_bytes)
    assert opens == [path]
    with pytest.raises(ValueError, match='Output changed'):
        jobs.verified_output_rows(job, dict(receipt, output_sha256='0' * 64))
    with pytest.raises(ValueError, match='Output count changed'):
        jobs.verified_output_rows(job, dict(receipt, input_rows=2, output_rows=2))
