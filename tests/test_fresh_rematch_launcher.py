"""Exercise the full-sample shell launcher without a Condor scheduler or EOS mount."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest


PRODUCTION = Path(__file__).resolve().parents[1] / 'production'
PREFIX = 'condor-1S-9x275-full-grouped-eos-'


@pytest.fixture
def campaign(tmp_path, monkeypatch):
    repo = tmp_path / 'repo'
    shutil.copytree(PRODUCTION, repo / 'production')
    work = repo / 'workspaces'
    pilot = work / 'condor-1S-9x275-check-v2/batch/jobs/job_000000/result.json'
    pilot.parent.mkdir(parents=True)
    pilot.write_text(json.dumps(dict(status='DONE', state='1S', beam='9x275',
                                     output_rows=100, electron_id='reco-rematch3011')))
    data = tmp_path / 'data'
    src = data / 'Upsilon1S/9x275/hiAcc'
    src.mkdir(parents=True)
    for i in range(3):
        name = f'eSTARlight1.3.0-1.0_Upsilon1S_9x275_q2_0to0.01_hiAcc_ab_run000.{i:04d}.eicrecon.edm4eic.root'
        (src / name).write_bytes(b'root')
    commands = tmp_path / 'bin'
    commands.mkdir()
    for name, body in {
            'condor_q': 'cat "$MOCK_QUEUE"',
            'condor_submit': 'echo "$*" >> "$MOCK_SUBMITS"\necho "2 job(s) submitted to cluster 77."',
    }.items():
        tool = commands / name
        tool.write_text('#!/bin/sh\n' + body + '\n')
        tool.chmod(0o755)
    for key, value in dict(PATH=str(commands) + os.pathsep + os.environ['PATH'],
                           UPSILON_DATA=data, CORRECTED_BASE=tmp_path / 'corrected',
                           FILES_PER_JOB='2', MOCK_QUEUE=tmp_path / 'queue',
                           MOCK_SUBMITS=tmp_path / 'submits').items():
        monkeypatch.setenv(key, str(value))
    monkeypatch.delenv('TUPLE_OUTPUT_DIR', raising=False)
    (tmp_path / 'queue').touch()
    return repo


def launch(repo, *args):
    return subprocess.run(['bash', 'production/submit_1s_9x275.sh', *args], cwd=repo,
                          text=True, capture_output=True, timeout=30)


def test_named_campaign_submits_all_inputs_with_earlier_campaign_active(campaign, tmp_path):
    work = campaign / 'workspaces'
    old = work / 'condor-1S-9x275-full-grouped-v1'
    old.mkdir()
    previous = '11 job(s) submitted to cluster 44.\n'
    (old / 'submission.txt').write_text(previous)
    (old / 'input.txt').write_text('earlier inventory\n')
    # The second active name shares the new run's prefix.
    other = work / (PREFIX + 'test-extra')
    (tmp_path / 'queue').write_text(f'44 {old}/batch worker\n55 {other}/batch worker\n')
    preview = launch(campaign, '--new-run', 'test')
    assert preview.returncode == 0, preview.stderr
    assert '--new-run test --submit' in preview.stdout
    assert not (tmp_path / 'submits').exists()
    result = launch(campaign, '--submit', '--new-run', 'test')
    assert result.returncode == 0, result.stderr
    full = work / (PREFIX + 'test')
    batch = full / 'batch'
    plan = json.loads((batch / 'plan.json').read_text())
    output = tmp_path / 'corrected/Upsilon1S/9x275/hiAcc/output_newID/runs/test/intermediate'
    assert plan['tuple_output_dir'] == str(output)
    assert output.is_dir() and not list(output.iterdir())
    assert (plan['job_count'], plan['condor_job_count']) == (3, 2)
    for i, jobfile in enumerate(plan['jobs']):
        assert json.loads(Path(jobfile).read_text())['output'] == str(output / f'job_{i:06d}.root')
    assert (batch / 'submit.sub').read_text().endswith('000000\n000001\n)\n')
    assert (full / 'submission.txt').read_text() == '2 job(s) submitted to cluster 77.\n'
    assert (old / 'submission.txt').read_text() == previous
    again = launch(campaign, '--new-run', 'test', '--submit')
    assert again.returncode != 0
    assert len((tmp_path / 'submits').read_text().splitlines()) == 1


def test_same_named_campaign_in_queue_blocks_even_without_receipt(campaign, tmp_path):
    full = campaign / 'workspaces' / (PREFIX + 'test')
    (tmp_path / 'queue').write_text(f'44 {full}/batch {full}/batch/snapshot/run_job.sh\n')
    result = launch(campaign, '--submit', '--new-run', 'test')
    assert result.returncode != 0
    assert not (tmp_path / 'submits').exists()
    assert not (full / 'submission.txt').exists()


@pytest.mark.parametrize('args', [('--new-run', 'v1'), ('--new-run', '../other'),
                                  ('--new-run', 'test', '--after-removal')])
def test_invalid_named_or_mixed_campaign_rejected(campaign, tmp_path, args):
    result = launch(campaign, *args, '--submit')
    assert result.returncode != 0
    assert not (tmp_path / 'submits').exists()
    assert not (tmp_path / 'corrected').exists()
