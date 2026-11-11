"""Exercise campaign receipt and scheduler reporting without a live schedd or ROOT."""
from contextlib import redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
from threading import Event
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

    def test_queue_overview_accepts_multiple_dirs_and_repeated_batch_flags(self):
        paths = [self.batch.parent / f'camp-{i}/batch' for i in range(3)]
        with patch.object(status, 'report') as report:
            status.main(['--queue-only', '--batch', str(paths[0]), str(paths[1]), '--batch', str(paths[2])])
        self.assertEqual([call.args[0] for call in report.call_args_list], paths)
        self.assertTrue(all(call.args[3] for call in report.call_args_list))

    def test_overview_continues_past_an_unplanned_cell_but_exits_nonzero(self):
        with patch.object(status, 'report', side_effect=[FileNotFoundError('missing plan'), None]) as report, \
             patch.object(sys, 'stderr', io.StringIO()) as errors, self.assertRaises(SystemExit) as exc:
            status.main(['--batch', 'not-yet-planned', str(self.batch), '--queue-only'])
        self.assertEqual(exc.exception.code, 2)
        self.assertEqual(report.call_count, 2)
        self.assertIn('not-yet-planned', errors.getvalue())

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

    def test_external_tuple_health_and_partial_state(self):
        tuples = self.batch.parent / "tuples"
        tuples.mkdir()
        path = self.batch / "plan.json"
        plan = json.loads(path.read_text())
        plan["tuple_output_dir"] = str(tuples)
        path.write_text(json.dumps(plan))
        for i in (0, 1):
            output = tuples / f"job_{i:06d}.root"
            record = {"status": "DONE", "job": str(self.jobs[i] / "job.json"), "output": str(output),
                      "input_rows": 100, "output_rows": 100}
            (self.jobs[i] / "result.json").write_text(json.dumps(record))
        (tuples / "job_000000.root").write_text("verified tuple")
        (tuples / "job_000002.root").write_text("partial tuple without receipt")
        report = self.printed(query_queue=False)
        self.assertIn(str(tuples), report)
        self.assertIn("1/6 DONE (16.7%), FAILED=0, started/partial=1, untouched=3, unreadable=1", report)
        self.assertIn("Reported generated rows: 100 ", report)

    def test_parallel_status_preserves_manifest_order_for_bounded_examples(self):
        plan = json.loads((self.batch / 'plan.json').read_text())
        later = Event()

        def probe(directory, output, external_files):
            index = int(directory.name.split('_')[-1])
            if index == 0:
                self.assertTrue(later.wait(5))
            elif index == 1:
                later.set()
            return "FAILED", 0, f"input {index}"

        with patch.object(status, "file_state", probe):
            counts, rows, examples = status.scan_files(plan, 2, check_workers=2)
        self.assertEqual((counts["FAILED"], rows), (6, 0))
        self.assertEqual(examples["FAILED"], [(0, "input 0"), (1, "input 1")])

    def test_queue_only_does_not_read_worker_directories(self):
        with patch.object(status, "queue_rows", return_value=[]), patch.object(status, "scan_files") as scan:
            with redirect_stdout(io.StringIO()):
                status.report(self.batch, 2, queue_only=True)
        scan.assert_not_called()


if __name__ == "__main__":
    unittest.main()
