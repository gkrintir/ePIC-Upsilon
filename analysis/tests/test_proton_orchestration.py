"""Exercise actual entry-point function bodies with synthetic IO/tasks only."""
import ast
import __future__
from contextlib import contextmanager, redirect_stdout
import datetime as dt
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflow"
sys.path.insert(0, str(WORKFLOW))
from proton_selection import guard_proton_run


def function_from_script(script, name, namespace):
    path = WORKFLOW / script
    parsed = ast.parse(path.read_text(), filename=str(path))
    function = next(node for node in parsed.body if isinstance(node, ast.FunctionDef) and node.name == name)
    namespace = dict(namespace, __file__=str(path), __package__="workflow")
    module = ast.Module(body=[function], type_ignores=[])
    exec(compile(module, str(path), "exec", flags=__future__.annotations.compiler_flag), namespace)
    return namespace[name]


class FakeParser:
    def __init__(self, args): self.args = args
    def set_defaults(self, **kwargs): pass
    def parse_args(self, *args): return self.args
    def error(self, text): raise AssertionError(text)


class ProtonOrchestrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.env = patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"})
        self.env.start()

    def tearDown(self):
        self.env.stop()
        self.temp.cleanup()

    def test_statistical_wrappers_create_proton_contract_before_statistics_file(self):
        for script in ("run_statistics.py", "run_statistics_continue.py"):
            with self.subTest(script=script), redirect_stdout(io.StringIO()):
                run = self.root / script
                cfg = {"run": {"dir": str(run), "id": "synthetic"}}
                events = []
                def check_contract(path, contract):
                    # Real statistical-contract creation makes the directory nonempty.
                    self.assertTrue((path / "proton_selection.json").is_file())
                    path.mkdir(parents=True, exist_ok=True)
                    (path / "statistical_contract.json").write_text(json.dumps(contract))
                    events.append("statistical")
                    return "CREATED"
                def run_command(args):
                    guard_proton_run(cfg, self.root)
                    events.append("runner")
                    return 0
                args = SimpleNamespace(command="run", dry_run=False, run_id="synthetic",
                                       keep_going=False, func=run_command)
                runner = SimpleNamespace(PROJECT_ROOT=self.root, ConfigurationError=RuntimeError,
                    build_parser=lambda: FakeParser(args), load_configuration=lambda *args: cfg,
                    flatten=lambda value: {})
                @contextmanager
                def continuation(**kwargs):
                    yield {"summary": None}
                main = function_from_script(script, "main", {
                    "runner": runner, "Path": Path, "dt": dt,
                    "build_contract": lambda *args: {"synthetic": True},
                    "check_contract": check_contract,
                    "install_continuation": continuation,
                })
                self.assertEqual(main(), 0)
                self.assertEqual(events, ["statistical", "runner"])

    def allstates_module(self):
        spec = importlib.util.spec_from_file_location("synthetic_proton_allstates", WORKFLOW / "run_all_states.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ROOT = self.root
        return module

    def test_old_batch_manifest_is_not_overwritten_and_no_children_launch(self):
        module = self.allstates_module()
        old = self.root / "results/generated/batches/legacy"
        old.mkdir(parents=True)
        manifest = old / "batch_manifest.json"
        content = b'{"original": "keep for bookkeeping"}\n'
        manifest.write_bytes(content)
        with patch.object(module.subprocess, "run") as child, redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                module.main(["run", "--states", "1S", "--batch-id", "legacy"])
        self.assertEqual(manifest.read_bytes(), content)
        self.assertFalse((old / "proton_selection.json").exists())
        child.assert_not_called()

    def test_fresh_batch_and_identical_continuation_work_but_policy_change_fails(self):
        module = self.allstates_module()
        argv = ["run", "--states", "1S", "--batch-id", "fresh"]
        with patch.object(module.subprocess, "run", return_value=SimpleNamespace(returncode=0)) as child, redirect_stdout(io.StringIO()):
            self.assertEqual(module.main(argv), 0)
            self.assertEqual(module.main(argv), 0)
            self.assertEqual(child.call_count, 2)
            path = self.root / "results/generated/batches/fresh/batch_manifest.json"
            before = path.read_bytes()
            with patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "inclusive"}):
                with self.assertRaises(ValueError): module.main(argv)
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(child.call_count, 2)

    def test_package_style_cmd_run_resolves_guard_and_blocks_old_run(self):
        run = self.root / "run"
        run.mkdir()
        artifact = run / "old.csv"
        artifact.write_bytes(b"old numerical outputs")
        cfg = {"run": {"dir": str(run), "id": "synthetic"}}
        args = SimpleNamespace(run_id="synthetic", dry_run=False)
        command = function_from_script("run_analysis.py", "cmd_run", {
            "dt": dt, "sys": sys, "Path": Path, "PROJECT_ROOT": self.root,
            "load_configuration": lambda *args: cfg, "ConfigurationError": RuntimeError,
        })
        before = list(sys.path)
        saved_helper = sys.modules.pop("proton_selection", None)
        try:
            sys.path[:] = [entry for entry in sys.path if Path(entry or '.').resolve() != WORKFLOW.resolve()]
            with self.assertRaisesRegex(RuntimeError, "fresh run ID"):
                command(args)
        finally:
            sys.path[:] = before
            if saved_helper is not None: sys.modules["proton_selection"] = saved_helper
        self.assertEqual(artifact.read_bytes(), b"old numerical outputs")
        self.assertFalse((run / "proton_selection.json").exists())


    def test_dry_run_does_not_poison_a_fresh_run_id(self):
        run = self.root / "dry"
        cfg = {"run": {"dir": str(run), "id": "synthetic"}}
        args = SimpleNamespace(run_id="synthetic", dry_run=True, targets=[])
        def forbidden_manifest(*args):
            self.fail("A dry-run must not write a run manifest")
        command = function_from_script("run_analysis.py", "cmd_run", {
            "dt": dt, "sys": sys, "Path": Path, "PROJECT_ROOT": self.root,
            "load_configuration": lambda *args: cfg, "ConfigurationError": RuntimeError,
            "expand_targets": lambda *args: [], "dependency_closure": lambda *args: [],
            "task_index": lambda *args: {}, "verify_configured_checksums": lambda *args: [],
            "write_run_manifest": forbidden_manifest,
        })
        with redirect_stdout(io.StringIO()):
            self.assertEqual(command(args), 0)
        self.assertFalse(run.exists())
        guard_proton_run(cfg, self.root)
        self.assertTrue((run / "proton_selection.json").is_file())


if __name__ == "__main__":
    unittest.main()
