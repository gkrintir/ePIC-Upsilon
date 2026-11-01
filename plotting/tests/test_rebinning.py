"""Independent tests of display aggregation and requested plotting changes."""
from contextlib import contextmanager
from copy import deepcopy
import sys
import tempfile
import types
import unittest
from unittest import mock
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from test_reoptimization_inputs import make_fixture
from upsilon_esr.reoptimization_inputs import load_case
from upsilon_esr.rebinning import merge_matrix, maybe_rebin


class RebinningTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cfg, _ = make_fixture(self.temp.name)
        self.warning_context = warnings.catch_warnings()
        self.warning_context.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)
        self.addCleanup(self.warning_context.__exit__, None, None, None)
        rng = np.random.default_rng(932)
        switch = 20. * rng.integers(0, 2, 400)
        # A strong anticorrelation makes summing marginal intervals or diagonal
        # variances conspicuously wrong. Toy sums retain a nonnegative tail.
        self.toys = np.column_stack([rng.poisson(3., 400), switch,
                                    20.-switch, rng.poisson(2., 400),
                                    rng.poisson(.5, 400)]).astype(float)
        self.covariance = np.cov(self.toys, rowvar=False, ddof=1)

    def case(self, observable):
        frame, meta, _, entry = load_case(self.cfg, "1S", observable)
        entry = deepcopy(entry)
        entry["inference"]["settings"] = {"coverage": self.cfg.coverage}
        entry["inference"]["coarse"]["covariance"] = self.covariance.tolist()
        frame["toy_mean"] = self.toys.mean(axis=0)
        alpha = 1. - self.cfg.coverage
        frame["toy_lower"], frame["toy_upper"] = np.quantile(
            self.toys, [alpha/2., 1.-alpha/2.], axis=0)
        factors = np.asarray(entry["cross_section_scaling"][
            "alternate_factors" if observable == "y" else "factors"])
        meta["covariance"] = (self.covariance*factors[:,None]*factors[None,:]).tolist()
        return frame, meta, entry

    @contextmanager
    def replay(self):
        # The replay solver has its own integration tests. Isolate aggregation
        # here with known non-Gaussian correlated joint toy estimates.
        module = types.ModuleType("upsilon_esr.rebin_toys")
        module.replay_joint_toys = mock.Mock(return_value=(self.toys, {"synthetic": True}))
        with mock.patch.dict(sys.modules, {"upsilon_esr.rebin_toys": module}):
            yield module.replay_joint_toys

    def test_w_flux_weighting_yield_and_full_covariance(self):
        self.cfg.rebin_edges = {"1S": {"w": [20., 50., 65., 100.]}}
        frame, meta, entry = self.case("w")
        original = frame.copy(deep=True)
        with self.replay() as replay:
            merged, result = maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)
        replay.assert_called_once()
        self.assertEqual(len(merged), 3)
        old_flux = np.asarray(entry["cross_section_scaling"]["bin_flux"])
        expected = np.average(frame["mean"].iloc[1:4], weights=old_flux[1:4])
        self.assertAlmostEqual(merged["mean"].iloc[1], expected)
        self.assertNotAlmostEqual(expected, frame["mean"].iloc[1:4].mean())
        np.testing.assert_allclose(merged.fitted_truth, [frame.fitted_truth.iloc[0],
                                    frame.fitted_truth.iloc[1:4].sum(), frame.fitted_truth.iloc[-1]])
        self.assertAlmostEqual(merged.fitted_truth.sum(), frame.fitted_truth.sum())
        new_flux = np.asarray(result["cross_section_scaling"]["bin_flux"])
        self.assertAlmostEqual(np.dot(merged["mean"], new_flux), np.dot(frame["mean"], old_flux))
        transform = np.asarray(result["rebinning"]["cross_section_merge_matrix"])
        expected_covariance = transform @ np.asarray(meta["covariance"]) @ transform.T
        np.testing.assert_allclose(result["covariance"], expected_covariance)
        diagonal_only = transform @ np.diag(np.diag(meta["covariance"])) @ transform.T
        self.assertGreater(abs(expected_covariance[1,1] - diagonal_only[1,1]), 1.)
        self.assertEqual(set(merged.overall_status), {"NOT_VALIDATED"})
        self.assertFalse(merged.publication_ready.any())
        self.assertEqual(result["overall_status"], "NOT_VALIDATED")
        self.assertFalse(result["rebinning"]["validation_status_changed"])
        pd.testing.assert_frame_equal(frame, original)

    def test_y_width_weighting_primary_and_flux_normalized_companion(self):
        # Merging unequal widths tests a genuinely weighted density average.
        self.cfg.rebin_edges = {"1S": {"y": [-1., .5, 1., 3.]}}
        frame, meta, entry = self.case("y")
        with self.replay():
            merged, result = maybe_rebin(self.cfg, "1S", "y", frame, meta, entry)
        widths = frame.bin_width.to_numpy(float)
        expected = np.average(frame["mean"].iloc[:2], weights=widths[:2])
        self.assertAlmostEqual(merged["mean"].iloc[0], expected)
        self.assertNotAlmostEqual(expected, frame["mean"].iloc[:2].mean())
        np.testing.assert_allclose(merged.fitted_truth_alternate_cross_section,
                                   merged.fitted_truth_cross_section / entry["flux_integral"])
        self.assertAlmostEqual(np.dot(merged["mean"], merged.bin_width), np.dot(frame["mean"], widths))
        self.assertAlmostEqual(np.dot(merged.fitted_truth_cross_section, merged.bin_width),
                               np.dot(frame.fitted_truth_cross_section, widths))
        self.assertEqual(result["displayed_quantity"], "flux_normalized_ep_rapidity")

    def test_joint_toy_quantiles_not_added_endpoints(self):
        self.cfg.rebin_edges = {"1S": {"w": [20., 50., 65., 100.]}}
        frame, meta, entry = self.case("w")
        with self.replay():
            merged, _ = maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)
        alpha = 1. - self.cfg.coverage
        expected = np.quantile(self.toys[:,1:4].sum(axis=1), [alpha/2., 1.-alpha/2.])
        np.testing.assert_allclose([merged.toy_lower.iloc[1], merged.toy_upper.iloc[1]], expected)
        self.assertNotAlmostEqual(merged.toy_lower.iloc[1], frame.toy_lower.iloc[1:4].sum())
        self.assertNotAlmostEqual(merged.toy_upper.iloc[1], frame.toy_upper.iloc[1:4].sum())

    def test_identity_and_disabled_rebin_leave_inputs_exactly_unchanged(self):
        frame, meta, entry = self.case("w")
        with self.replay() as replay:
            result, info = maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)
            self.assertIs(result, frame)
            self.assertIs(info, meta)
            self.cfg.rebin_edges = {"1S": {"w": [20.,50.,55.,60.,65.,100.]}}
            result, info = maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)
            self.assertIs(result, frame)
            self.assertIs(info, meta)
        replay.assert_not_called()

    def test_default_five_to_three_groups_and_invalid_edges(self):
        matrix, edges = merge_matrix([20,50,55,60,65,100], [20,50,65,100])
        np.testing.assert_array_equal(matrix, [[1,0,0,0,0], [0,1,1,1,0], [0,0,0,0,1]])
        np.testing.assert_array_equal(edges, [20,50,65,100])
        for requested in ([20,45,65,100], [20,50,65], [20,65,50,100]):
            with self.assertRaises(ValueError):
                merge_matrix([20,50,55,60,65,100], requested)

    def test_automatic_groups_follow_changed_saved_edges(self):
        for old in ([20,40,55,75,100], [20,35,45,60,80,100],
                    [20,30,40,50,60,80,100]):
            with self.subTest(old=old):
                group, edges = merge_matrix(old, "keep_outer_bins")
                self.assertEqual(group.shape, (3, len(old)-1))
                np.testing.assert_array_equal(edges, [old[0],old[1],old[-2],old[-1]])
                np.testing.assert_array_equal(group.sum(axis=0), np.ones(len(old)-1))
                self.assertEqual(group[0].sum(), 1)
                self.assertEqual(group[-1].sum(), 1)
                self.assertTrue(set(edges).issubset(old))
        with self.assertRaisesRegex(ValueError, "Saved edges"):
            merge_matrix([20,35,45,60,80,100], [20,50,65,100])
        with self.assertRaisesRegex(ValueError, "Unknown"):
            merge_matrix([20,40,60,100], "snap_edges")
        for old in ([20,np.nan,100], [20,20,100], [100,20], [20]):
            with self.assertRaises(ValueError):
                merge_matrix(old, "keep_outer_bins")

    def test_three_or_fewer_saved_bins_are_never_split(self):
        group, edges = merge_matrix([20,40,60,100], "keep_outer_bins")
        np.testing.assert_array_equal(group, np.eye(3))
        for old in ([20,60,100], [20,100]):
            with self.assertWarnsRegex(RuntimeWarning, "keeping them unchanged"):
                group, edges = merge_matrix(old, "keep_outer_bins")
            np.testing.assert_array_equal(group, np.eye(len(old)-1))
            np.testing.assert_array_equal(edges, old)

    def test_automatic_merge_matches_explicit_covariance_quantiles_and_tag_metadata(self):
        frame, meta, entry = self.case("w")
        frame["proton_selection_mode"] = "required"
        frame["proton_selection_contract"] = "jul26-reco-topology-v1"
        frame["x_low"], frame["x_high"] = [20,35,45,60,80], [35,45,60,80,100]
        self.cfg.rebin_edges = {"1S": {"w": "keep_outer_bins"}}
        with self.replay():
            automatic, auto_meta = maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)
        self.cfg.rebin_edges = {"1S": {"w": [20,35,80,100]}}
        with self.replay():
            explicit, explicit_meta = maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)
        pd.testing.assert_frame_equal(automatic, explicit)
        np.testing.assert_array_equal(auto_meta["covariance"], explicit_meta["covariance"])
        self.assertEqual(auto_meta["rebinning"]["edge_selection"], "keep_outer_bins")
        self.assertEqual(set(automatic.proton_selection_mode), {"required"})
        self.assertEqual(set(automatic.proton_selection_contract), {"jul26-reco-topology-v1"})
        frame.loc[2,"proton_selection_mode"] = "inclusive"
        with self.replay(), self.assertRaisesRegex(ValueError, "inconsistent proton_selection_mode"):
            maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)

    def test_preflight_accepts_automatic_rule_for_w_and_y(self):
        from upsilon_esr.reoptimization_inputs import preflight_reoptimization
        self.cfg.rebin_edges = {"1S": {"w": "keep_outer_bins", "y": "keep_outer_bins"}}
        for section, observable in (("w", "w"), ("differential", "y")):
            rows = preflight_reoptimization(self.cfg, section)
            row = rows.loc[(rows.state == "1S") & (rows.observable == observable)].iloc[0]
            self.assertEqual(row["bins"], 5)
            self.assertEqual(row["display_bins"], 3)

    def test_pt2_rebin_and_inconsistent_metadata_are_rejected(self):
        frame, meta, entry = self.case("w")
        self.cfg.rebin_edges = {"1S": {"pt2": [0,.15,.4,1]}}
        with self.assertRaisesRegex(ValueError, "W and rapidity only"):
            maybe_rebin(self.cfg, "1S", "pt2", frame, meta, entry)
        self.cfg.rebin_edges = {"1S": {"w": [20,50,65,100]}}
        frame.loc[1, "overall_status"] = "PASS_CONDITIONAL"
        with self.replay(), self.assertRaisesRegex(ValueError, "inconsistent overall_status"):
            maybe_rebin(self.cfg, "1S", "w", frame, meta, entry)

    def test_w_hides_mc_and_draws_reco_above_zeus(self):
        from upsilon_esr.w_cross_section import plot_w_cross_section
        with mock.patch("upsilon_esr.w_cross_section._steps", side_effect=AssertionError("MC curve drawn")):
            fig, audit = plot_w_cross_section(self.cfg)
        self.addCleanup(plt.close, fig)
        self.assertEqual([t.get_text() for t in fig.axes[0].get_legend().get_texts()], ["RECO", "ZEUS"])
        for axis in fig.axes:
            self.assertNotIn("MC", axis.get_legend_handles_labels()[1])
        containers = {c.get_label(): c for c in fig.axes[0].containers}
        self.assertGreater(containers["RECO"].lines[0].get_zorder(), containers["ZEUS"].lines[0].get_zorder())
        self.assertFalse(audit["display_annotations"]["MC_curves_displayed"])

    def test_pt2_legend_is_upper_left(self):
        from upsilon_esr.differential import plot_differential
        fig, _ = plot_differential(self.cfg)
        self.addCleanup(plt.close, fig)
        legend = fig.axes[0].get_legend()
        self.assertEqual(legend._loc, 2)  # matplotlib upper-left location
        self.assertEqual([t.get_text() for t in legend.get_texts()][-1],
                         r"Fit: $0<p_T^2<1$ GeV$^2$")


if __name__ == "__main__":
    unittest.main()
