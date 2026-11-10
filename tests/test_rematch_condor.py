"""Campaign-scoped scheduler queries; never contact a schedd in these tests."""
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'production'))
from rematch_condor import queue_rows, recorded_clusters


class ReadCondorTest(TestCase):
    def test_collect_only_grouped_attempts_and_fail_closed_on_ambiguous_submission(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            batch = root / 'grouped/batch'
            attempt = batch / 'retries/attempt-0001'
            attempt.mkdir(parents=True)
            (batch.parent / 'submission.txt').write_text('114 job(s) submitted to cluster 1234.\n')
            unrelated = root / 'earlier-one-file-campaign'
            unrelated.mkdir()
            (unrelated / 'submission.txt').write_text('999 job(s) submitted to cluster 55.\n')
            self.assertEqual(recorded_clusters(batch), [1234])
            (attempt / 'submission.txt').write_text('2 job(s) submitted to cluster 1237.\n')
            self.assertEqual(recorded_clusters(batch), [1234, 1237])
            (attempt / 'submission.txt').write_text('Submitting job(s).')
            with self.assertRaisesRegex(ValueError, 'No confirmed'):
                recorded_clusters(batch)

    def test_confirm_empty_selection_and_keep_query_scoped(self):
        for response in ('[]', ''):
            with patch('rematch_condor.subprocess.run', return_value=subprocess.CompletedProcess(
                    [], 0, response, '')) as run:
                self.assertEqual(queue_rows(Path('/test/batch'), [72, 75]), [])
                cmd = run.call_args.args[0]
                constraint = cmd[cmd.index('-constraint') + 1]
                for token in ('Owner', 'ClusterId == 72', 'ClusterId == 75', '/test/batch/jobs/', 'snapshot/run_job.sh'):
                    self.assertIn(token, constraint)

    def test_scheduler_failure_or_bad_response_is_not_an_empty_queue(self):
        for code, out, err, error in ((1, '', 'schedd error', RuntimeError), (0, '', 'connection lost', RuntimeError),
                                      (0, '{}', '', ValueError), (0, '[{}]', '', ValueError)):
            with patch('rematch_condor.subprocess.run',
                       return_value=subprocess.CompletedProcess([], code, out, err)), self.assertRaises(error):
                queue_rows(Path('/test/batch'), [])
