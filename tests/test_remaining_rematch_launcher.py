"""Run the July matrix launcher with real frozen plans and a simulated scheduler."""
import json
from pathlib import Path
import subprocess

import pytest

from test_fresh_rematch_launcher import campaign


REMAINING = [('1S', '9x130'), ('2S', '9x130'), ('2S', '9x275'),
             ('3S', '9x130'), ('3S', '9x275')]
RUN = 'rest-01'


def control(repo, state, beam):
    return repo / 'workspaces' / f'condor-{state}-{beam}-full-grouped-eos-{RUN}'


def launch(repo, *options):
    return subprocess.run(
        ['bash', 'production/submit_rematch.sh', '--new-run', RUN, *options], cwd=repo,
        text=True, capture_output=True, timeout=60)


@pytest.fixture
def matrix(campaign, tmp_path):
    for state, beam in REMAINING:
        src = tmp_path / 'data' / f'Upsilon{state}' / beam / 'hiAcc'
        src.mkdir(parents=True)
        for run in ('003.0019', '003.0017', '003.0021'):
            name = f'eSTARlight1.3.0-1.0_Upsilon{state}_{beam}_q2_0to0.01_hiAcc_ab_run{run}.eicrecon.edm4eic.root'
            (src / name).write_bytes(b'root')
    (tmp_path / 'bin/condor_submit').write_text(
        '#!/bin/sh\n'
        'echo "$*" >> "$MOCK_SUBMITS"\n'
        'cluster=$(wc -l < "$MOCK_SUBMITS")\n'
        'echo "2 job(s) submitted to cluster $((100 + cluster))."\n')
    return campaign


def submissions(tmp_path):
    log = tmp_path / 'submits'
    return log.read_text().splitlines() if log.exists() else []


def test_full_matrix_plans_then_submits_each_remaining_pair_once(matrix, tmp_path):
    existing = matrix / 'workspaces/condor-1S-9x275-full-grouped-eos-fresh-01'
    existing.mkdir()
    receipt = existing / 'submission.txt'
    receipt.write_text('114 job(s) submitted to cluster 12.\n')
    (tmp_path / 'queue').write_text(f'12 {existing}/batch {existing}/batch/snapshot/run_job.sh\n')
    preview = launch(matrix, '--remaining')
    assert preview.returncode == 0, preview.stderr
    assert not submissions(tmp_path)
    # A partial invocation may have submitted this cell already.
    one = launch(matrix, '--state', '3S', '--beam', '9x130', '--submit')
    assert one.returncode == 0, one.stderr
    assert len(submissions(tmp_path)) == 1
    live = control(matrix, '3S', '9x130') / 'batch'
    with (tmp_path / 'queue').open('a') as queue:
        queue.write(f'101 {live} {live}/snapshot/run_job.sh\n')
    finish = launch(matrix, '--remaining', '--submit')
    assert finish.returncode == 0, finish.stderr
    calls = submissions(tmp_path)
    assert len(calls) == len(REMAINING)
    for state, beam in REMAINING:
        directory = control(matrix, state, beam)
        batch = directory / 'batch'
        plan = json.loads((batch / 'plan.json').read_text())
        assert (plan['state'], plan['beam'], plan['electron_id'], plan['femc']) == (
            state, beam, 'reco-rematch3011', 'full')
        assert (plan['job_count'], plan['files_per_job'], plan['condor_job_count']) == (3, 2, 2)
        src = tmp_path / 'data' / f'Upsilon{state}' / beam / 'hiAcc'
        original_files = sorted(str(path) for path in src.glob('*.root'))
        assert [item['path'] for item in plan['inputs']] == original_files
        assert (directory / 'input.txt').read_text().splitlines() == original_files
        macro = ('macros/9x130/' if beam == '9x130' else 'macros/') + 'DVMP_Upsilon_Analysis.C'
        assert plan['source_paths']['DVMP_Upsilon_Analysis.C'] == macro
        layout = 'output_newID' if state == '1S' else 'output'
        out = tmp_path / 'corrected' / f'Upsilon{state}' / beam / 'hiAcc' / layout / 'runs' / RUN / 'intermediate'
        assert Path(plan['tuple_output_dir']) == out and out.is_dir()
        assert (directory / 'submission.txt').is_file()
        submit = batch / 'submit.sub'
        assert '+AccountingGroup = "group_u_CMS.CAF.COMM"' in submit.read_text()
        assert sum(str(submit) in line and 'eic_xl:26.07.1-stable' in line and '-maxjobs 2' in line
                   for line in calls) == 1
    assert not control(matrix, '1S', '9x275').exists()
    assert receipt.read_text() == '114 job(s) submitted to cluster 12.\n'
    rerun = launch(matrix, '--remaining', '--submit')
    assert rerun.returncode == 0, rerun.stderr
    assert submissions(tmp_path) == calls


def test_unconfirmed_submission_never_submits_that_cell_again(matrix, tmp_path):
    preview = launch(matrix, '--remaining')
    assert preview.returncode == 0, preview.stderr
    ambiguous = control(matrix, *REMAINING[0])
    receipt = ambiguous / 'submission.txt'
    receipt.write_text('Submitting job(s).\n')
    result = launch(matrix, '--remaining', '--submit')
    assert result.returncode != 0
    assert all(str(ambiguous) not in call for call in submissions(tmp_path))
    assert receipt.read_text() == 'Submitting job(s).\n'


def test_matrix_rejects_single_tuple_directory_override(matrix, tmp_path, monkeypatch):
    override = tmp_path / 'shared-tuples'
    monkeypatch.setenv('TUPLE_OUTPUT_DIR', str(override))
    result = launch(matrix, '--remaining', '--submit')
    assert result.returncode != 0
    assert not override.exists()
    assert not submissions(tmp_path)


def test_new_matrix_requires_proven_pilot_before_submitting(matrix, tmp_path):
    pilot = matrix / 'workspaces/condor-1S-9x275-check-v2/batch/jobs/job_000000/result.json'
    pilot.unlink()
    result = launch(matrix, '--remaining', '--submit')
    assert result.returncode != 0
    assert not submissions(tmp_path)
