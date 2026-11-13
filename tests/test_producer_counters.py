"""Producer-counter rendering and notebook losses without a PyROOT dependency."""

import ast
import csv
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import StepPatch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis/scripts/reporting"))
import render_note_kinematics as renderer


COUNTS = np.array([1000, 1000, 900, 1200, 1000, 457, 400, 300, 250, 200], dtype=float)
ZERO_COUNTS = np.array([0, 0, 0, 0, 10, 0, 0, 0, 1, 0], dtype=float)


def notebook_loss_functions():
    path = ROOT / "analysis/notebooks/03_event_selection_truth_and_bt_studies.ipynb"
    notebook = json.loads(path.read_text())
    functions = {}
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell.get("source", []))
        tree = ast.parse(source, filename=f"{path}:cell{index}")
        definitions = [
            node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "calc_pct_lost"
        ]
        if definitions:
            namespace = {"np": np}
            module = ast.Module(body=definitions, type_ignores=[])
            exec(compile(module, str(path), "exec"), namespace)
            functions[index] = namespace["calc_pct_lost"]
    return functions


class ProducerCounterPages(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(plt.close, "all")
        self.output = Path(temporary.name)

    def render_counts(self, counts=COUNTS):
        sample = {
            "label": "9x275 GeV",
            "path": "synthetic.root",
            "stored": {
                "h_EleID_EpRatio": {
                    "values": np.array([10, 20, 30]),
                    "edges": np.array([0, 0.5, 1, 1.5]),
                },
                "h_EleID_IsoFraction": {
                    "values": np.array([4, 5, 6]),
                    "edges": np.array([0, 0.5, 1, 1.5]),
                },
                "h_CutFlow": {
                    "values": counts.copy(),
                    "edges": np.arange(0.5, 11),
                    "labels": [f"Legacy stage {stage}" for stage in range(1, 11)],
                },
            },
        }
        args = SimpleNamespace(
            out_dir=self.output, state="1S", stat_mode="projected",
            lumi_invfb=10, projection_lumi_invfb=2.5,
        )
        rows, audit = [], []
        with mock.patch.object(renderer, "finish") as finish:
            renderer.stored_pages([sample], args, rows, audit)
        figures = {call.args[2]: call.args[0] for call in finish.call_args_list}
        return rows, audit, figures["event_selection_cutflow_numbered"]

    def test_raw_counters_and_pid_entries_survive_projected_arguments(self):
        rows, audit, _ = self.render_counts()
        for histogram, expected in (
            ("h_CutFlow", COUNTS),
            ("h_EleID_EpRatio", [10, 20, 30]),
            ("h_EleID_IsoFraction", [4, 5, 6]),
        ):
            with self.subTest(histogram=histogram):
                actual = [row["count"] for row in rows if row["histogram"] == histogram]
                np.testing.assert_array_equal(actual, expected)
        self.assertFalse(any("missing_stored_histogram" in entry for entry in audit))
        counter_audit = next(entry for entry in audit if entry.get("histogram") == "h_CutFlow")
        self.assertEqual(counter_audit["stored_labels"], [f"Legacy stage {stage}" for stage in range(1, 11)])
        self.assertEqual(counter_audit["nested_stage_sequences"], [[1, 2, 3], [1, 4, 5, 6, 7, 8, 9, 10]])
        self.assertEqual(counter_audit["undefined_adjacent_loss_stages"], [1, 4])

    def test_undefined_losses_are_nan_on_plot_and_blank_in_csv(self):
        rows, _, figure = self.render_counts()
        loss_steps = [patch for patch in figure.axes[1].patches if isinstance(patch, StepPatch)]
        self.assertEqual(len(loss_steps), 1)
        loss = loss_steps[0].get_data().values
        self.assertIsNone(loss_steps[0].get_data().baseline)
        self.assertTrue(np.isnan(loss[0]))
        self.assertTrue(np.isnan(loss[3]))
        expected = np.array([0, 10, 100 / 6, 54.3, 100 * 57 / 457, 25, 100 / 6, 20])
        np.testing.assert_allclose(loss[[1, 2, 4, 5, 6, 7, 8, 9]], expected)

        path = self.output / "producer_counters.csv"
        renderer.write_csv(path, rows)
        with path.open(newline="") as source:
            counters = [row for row in csv.DictReader(source) if row["histogram"] == "h_CutFlow"]
        self.assertEqual(len(counters), 10)
        self.assertEqual(counters[0]["loss_from_previous_percent"], "")
        self.assertEqual(counters[3]["loss_from_previous_percent"], "")
        self.assertAlmostEqual(float(counters[5]["loss_from_previous_percent"]), 54.3)

    def test_zero_denominators_are_undefined_and_complete_losses_remain_100(self):
        rows, _, figure = self.render_counts(ZERO_COUNTS)
        steps = [patch for patch in figure.axes[1].patches if isinstance(patch, StepPatch)]
        loss = steps[0].get_data().values
        np.testing.assert_array_equal(np.flatnonzero(np.isfinite(loss)), [5, 9])
        np.testing.assert_array_equal(loss[[5, 9]], [100, 100])
        counters = [row for row in rows if row["histogram"] == "h_CutFlow"]
        for index in [0, 1, 2, 3, 4, 6, 7, 8]:
            with self.subTest(stage=index + 1):
                self.assertEqual(counters[index]["loss_from_previous_percent"], "")

    def test_count_curve_breaks_between_generator_and_reconstruction_counters(self):
        _, _, figure = self.render_counts()
        steps = [patch for patch in figure.axes[0].patches if isinstance(patch, StepPatch)]
        self.assertEqual(len(steps), 2)
        generator, reconstruction = [patch.get_data() for patch in steps]
        np.testing.assert_array_equal(generator.values, COUNTS[:3])
        np.testing.assert_array_equal(reconstruction.values, COUNTS[3:])
        self.assertEqual(generator.edges[-1], 3.5)
        self.assertEqual(reconstruction.edges[0], 3.5)
        self.assertIsNone(generator.baseline)
        self.assertIsNone(reconstruction.baseline)

    def test_stage_definitions_explain_availability_and_nonnegative_w(self):
        rows, _, figure = self.render_counts()
        counters = [row for row in rows if row["histogram"] == "h_CutFlow"]
        self.assertRegex(counters[5]["definition"].lower(), "availab")
        self.assertIn("central electron", counters[5]["definition"].lower())
        self.assertIn("tagger", counters[5]["definition"].lower())
        self.assertRegex(counters[8]["definition"], r"W.*(?:>=|≥|\\geq).*0")
        key = "\n".join(text.get_text() for text in figure.axes[2].texts)
        self.assertRegex(key.lower(), "availab")
        self.assertRegex(key, r"W.*(?:>=|≥|\\geq).*0")


class NotebookProducerLosses(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.functions = notebook_loss_functions()

    def test_both_cells_keep_non_nested_and_initial_losses_undefined(self):
        self.assertEqual(set(self.functions), {8, 10})
        for cell, calculate in self.functions.items():
            with self.subTest(cell=cell):
                counts = COUNTS.copy()
                loss = calculate(counts)
                self.assertTrue(np.isnan(loss[0]))
                self.assertTrue(np.isnan(loss[3]))
                self.assertAlmostEqual(loss[5], 54.3)
                self.assertAlmostEqual(loss[9], 20)
                np.testing.assert_array_equal(counts, COUNTS)

    def test_both_cells_keep_zero_denominators_undefined(self):
        for cell, calculate in self.functions.items():
            with self.subTest(cell=cell):
                loss = calculate(ZERO_COUNTS.copy())
                np.testing.assert_array_equal(np.flatnonzero(np.isfinite(loss)), [5, 9])
                np.testing.assert_array_equal(loss[[5, 9]], [100, 100])


if __name__ == "__main__":
    unittest.main()
