"""Exercise campaign receipt and scheduler reporting without a live schedd or ROOT."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


PRODUCTION = Path(__file__).resolve().parents[1] / "production"
SPEC = importlib.util.spec_from_file_location("rematch_status", PRODUCTION / "rematch_status.py")
status = importlib.util.module_from_spec(SPEC)
with patch.object(sys, "path", [str(PRODUCTION), *sys.path]):
    SPEC.loader.exec_module(status)


class RematchStatusTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.batch = Path(temp.name) / "batch"
        self.jobs = [self.batch / "jobs" / f"job_{i:06d}" for i in range(6)]
        for directory in self.jobs:
            directory.mkdir(parents=True)
            (directory / "job.json").write_text("{}")
        plan = {"work_dir": str(self.batch), "state": "1S", "beam": "9x275", "files_per_job": 2,
                "job_count": len(self.jobs), "inputs": [{"path": f"/input/{i}.root"} for i in range(6)],
                "jobs": [str(d / "job.json") for d in self.jobs]}
        (self.batch / "plan.json").write_text(json.dumps(plan))

    def done(self, i, rows):
        job = self.jobs[i]
        (job / "output.root").write_text("placeholder ROOT")
        (job / "result.json").write_text(json.dumps(
            {"status": "DONE", "job": str(job / "job.json"), "output": str(job / "output.root"),
             "input_rows": rows, "output_rows": rows}))

    def printed(self, limit=2, query_queue=True):
        out = io.StringIO()
        with redirect_stdout(out):
            status.report(self.batch, limit, query_queue)
        return out.getvalue()

    def test_summary_and_bounded_failure_logs_across_retry(self):
        self.done(0, 100)
        self.done(1, 101)
        (self.jobs[2] / "execution.lock").touch()
        (self.jobs[3] / "result.json").write_text('{"status":"FAILED","error":"ROOT failed"}')
        (self.jobs[3] / "root.log").touch()
        (self.jobs[4] / "result.json").write_text("{")
        attempt = self.batch / "retries" / "attempt-0001"
        (attempt / "logs").mkdir(parents=True)
        (attempt / "logs/group_000001.err").touch()
        (self.batch.parent / "submission.txt").write_text("3 job(s) submitted to cluster 52.\n")
        (attempt / "submission.txt").write_text("1 job(s) submitted to cluster 53.\n")
        with patch.object(status, "queue_rows", return_value=[
                {"ClusterId": 53, "ProcId": 0, "JobStatus": 5, "HoldReason": "out of memory"}]) as q:
            report = self.printed(1)
        q.assert_called_once_with(self.batch, [52, 53])
        self.assertIn("held 53.0: out of memory", report)
        self.assertIn("2/6 DONE (33.3%), FAILED=1, started/partial=1, untouched=1, unreadable=1", report)
        self.assertIn("201 (full-sample reference 1,137,600)", report)
        self.assertIn("/input/3.root", report)
        self.assertIn(str(self.jobs[3] / "root.log"), report)
        self.assertIn(str(attempt / "logs/group_000001.err"), report)
        self.assertNotIn(str(attempt / "logs/group_000000.err"), report)
        self.assertNotIn("untouched input", report)

    def test_queue_unavailable_preserves_progress(self):
        self.done(5, 100)
        with patch.object(status, "queue_rows", side_effect=FileNotFoundError("condor_q missing")):
            report = self.printed()
        self.assertIn("Condor queue unavailable", report)
        self.assertNotIn("Current Condor processes: 0", report)
        self.assertIn("1/6 DONE", report)
        (self.batch.parent / "submission.txt").touch()
        with patch.object(status, "queue_rows") as q:
            report = self.printed()
        q.assert_not_called()
        self.assertIn("No confirmed Condor cluster ID", report)
        with patch.object(status, "recorded_clusters") as clusters:
            report = self.printed(query_queue=False)
        clusters.assert_not_called()
        self.assertIn("1/6 DONE", report)

    def test_empty_queue_is_not_completion_and_limits_examples(self):
        with patch.object(status, "queue_rows", return_value=[]):
            report = self.printed(1)
        self.assertIn("Current Condor processes: 0", report)
        self.assertIn("0/6 DONE", report)
        self.assertEqual(report.count("untouched input "), 1)
        self.assertIn(str(self.batch / "logs/group_000000.err"), report)

    def test_missing_outputs_do_not_count_done_rows(self):
        self.done(0, 100)
        (self.jobs[0] / "output.root").unlink()
        status_name, rows, _ = status.file_state(self.jobs[0])
        self.assertEqual((status_name, rows), ("unreadable", 0))

    def test_queue_only_does_not_read_worker_directories(self):
        with patch.object(status, "queue_rows", return_value=[]), patch.object(status, "scan_files") as scan:
            with redirect_stdout(io.StringIO()):
                status.report(self.batch, 2, queue_only=True)
        scan.assert_not_called()


if __name__ == "__main__":
    unittest.main()
