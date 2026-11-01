"""Synthetic contract tests for saved reoptimization inputs (never physics data).

``make_fixture`` can also receive the framework's actual ``write_report``
function for an independent producer/consumer interoperability check.
"""
from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from upsilon_esr.common import Config, DEFAULT_COVERAGE, write_json
from upsilon_esr.reoptimization_inputs import load_case, saved_slope, w_inputs


def make_fixture(root, reporter=None):
    """A portable, explicitly synthetic v4 run with variable-width bins."""
    root = Path(root)
    run = root / "synthetic-v4-run"
    base = run / "artifacts/reoptimization"
    base.mkdir(parents=True)
    cfg = Config(allow_diagnostic_results=True, runs={s: root / ("unused-legacy-" + s) for s in ("1S", "2S", "3S")},
                 reoptimization_run=run, output_dir=root / "synthetic-figures", synthetic=True)
    source = dict(electron_energy=9., proton_energy=275., target_lumi_invfb=2.5,
                  w_range=[20., 100.], inference={"coverage": DEFAULT_COVERAGE},
                  uncertainty_scope="SYNTHETIC fixed-response test", synthetic=True)
    states = {}
    edges_by_state = {
        "1S": {"pt2": [0., .05, .15, .3, .6, 1.], "w": [20., 50., 55., 60., 65., 100.], "y": [-1., .25, .5, .75, 1., 3.]},
        "2S": {"pt2": [0., .1, .25, .45, 1.], "w": [20., 45., 60., 100.], "y": [-1., .5, 1., 3.]},
        "3S": {"pt2": [0., .15, .4, 1.], "w": [20., 45., 60., 100.], "y": [-1., .5, 1., 3.]},
    }
    for number, state in enumerate(cfg.runs):
        br = (.0238, .0191, .0218)[number]
        state_cfg = dict(state=state, br_ee=br, q2_max=cfg.q2_max[state],
                         category=cfg.categories[state], visible_xsec_fb=(210.136, 53.133, 42.969)[number],
                         sigma=(.053, .056, .059)[number])
        states[state] = {}
        for observable, listed_edges in edges_by_state[state].items():
            edges = np.asarray(listed_edges, float)
            n = len(edges) - 1
            truth = np.linspace(120., 20., n) / (number + 1)
            fitted = truth * 1.02
            covariance = np.diag(truth * 2.) + .08 * np.sqrt(np.outer(truth, truth))
            # Deliberately include one band that does not bracket its Asimov
            # estimate: a renderer must preserve, not clamp, sampling quantiles.
            lower, upper = fitted * .65, fitted * 1.4
            lower[-1], upper[-1] = fitted[-1] * 1.05, fitted[-1] * 1.5
            flux = .1012186984099357 if state == "1S" else .07446628165521689
            factor = 1. / (1000. * 2.5 * br)
            scaling = dict(luminosity_invfb=2.5, branching_fraction_ee=br, bin_edges=edges.tolist(),
                           units="pb/GeV^2" if observable == "pt2" else "pb")
            if observable == "w":
                bin_flux = np.linspace(.01, .025, n)
                scaling.update(quantity="sigma_gamma_p_w", bin_flux=bin_flux.tolist(),
                               factors=(factor / bin_flux).tolist())
            elif observable == "y":
                scaling.update(quantity="dsigma_ep_dy", factors=(factor / np.diff(edges)).tolist(),
                               alternate_factors=(factor / np.diff(edges) / flux).tolist(),
                               alternate_quantity="flux_normalized_ep_rapidity", alternate_units="pb")
            else:
                scaling.update(quantity="dsigma_gamma_p_dpt2", factors=(factor / flux / np.diff(edges)).tolist())
            reco = truth * .2
            reco_edges, indices = edges.copy(), list(range(n))
            if observable == "pt2":
                reco_edges, reco = np.r_[edges, 5.], np.r_[reco, 7.]
            elif observable == "y":
                reco_edges, reco, indices = np.r_[-2., edges, 4.], np.r_[11., reco, 13.], list(range(1, n+1))
            coarse = dict(truth_edges=edges.tolist(), reco_edges=reco_edges.tolist(),
                          truth_expected=truth.tolist(), reco_expected=reco.tolist(),
                          metadata={"fit_reco_indices": indices})
            extraction = dict(status="PASS_CONDITIONAL", truth_edges=edges.tolist(),
                              fitted_truth=fitted.tolist(), reference_truth=truth.tolist(),
                              toy_mean=fitted.tolist(), toy_lower=lower.tolist(), toy_upper=upper.tolist(),
                              sampling_sd=np.sqrt(np.diag(covariance)).tolist(), covariance=covariance.tolist())
            status = "PASS_CONDITIONAL" if observable == "pt2" else "NOT_VALIDATED"
            inference = dict(status="PASS_CONDITIONAL", coarse=extraction)
            if observable == "pt2":
                inference["slope"] = dict(status="OK_CONDITIONAL", b=2.7+number*.1,
                    interval_low=1.8, interval_high=4.2, amplitude=143./(number+1),
                    profile=[{"b": 0., "delta_deviance": 8.}, {"b": 2.7+number*.1, "delta_deviance": 0.},
                             {"b": 5., "delta_deviance": 3.}],
                    interval_interpretation="Saved conditional synthetic profile")
            entry = dict(state=state, observable=observable, complete=True, status=status, overall_status=status,
                         publication_ready=False, failed_checks=[] if observable == "pt2" else ["synthetic.stress_gate"],
                         state_config=state_cfg, selection_metadata={"source_luminosity_invfb": 10., "target_luminosity_invfb": 2.5},
                         optimization={"status": "SELECTED", "selected_edges": edges.tolist(), "candidates": []},
                         cross_section_scaling=scaling, flux_integral=flux,
                         cross_section_factor=factor/flux if observable == "pt2" else None,
                         model={"truth_edges": edges.tolist(), "coarse": deepcopy(coarse)},
                         validation_model={"coarse": coarse}, inference=inference, scope="SYNTHETIC TEST")
            states[state][observable] = entry
            write_json(base / state / observable / "result.json", entry)
    write_json(run / "effective_config.json", source)
    write_json(base / "all_results.json", states)
    if reporter is not None:
        reporter(base, states, source)
    else:
        # Same named columns/serialization as reporting.write_report. Runtime
        # producer interoperability is checked separately with reporter=... .
        rows, statuses = [], []
        for state, cases in states.items():
            for observable, entry in cases.items():
                status = entry["overall_status"]
                common = dict(state=state, observable=observable, status=status, overall_status=status,
                              stage_status="PASS_CONDITIONAL", publication_ready=False,
                              failed_checks=json.dumps(entry["failed_checks"]))
                statuses.append(common.copy())
                scaling, coarse = entry["cross_section_scaling"], entry["inference"]["coarse"]
                edges = np.asarray(coarse["truth_edges"])
                factors = np.asarray(scaling["factors"])
                alternate = scaling.get("alternate_factors")
                for i in range(len(factors)):
                    row = dict(common, bin=i, x_low=edges[i], x_high=edges[i+1],
                               quantity=scaling["quantity"], quantity_units=scaling["units"],
                               conversion_factor=factors[i],
                               interval_interpretation="Central toy sampling quantiles, not confidence limits")
                    for key in ("fitted_truth", "reference_truth", "toy_mean", "toy_lower", "toy_upper", "sampling_sd"):
                        row[key] = coarse[key][i]
                        row[key+"_cross_section"] = coarse[key][i] * factors[i]
                        if alternate is not None:
                            row[key+"_alternate_cross_section"] = coarse[key][i] * alternate[i]
                    if alternate is not None:
                        row.update(alternate_conversion_factor=alternate[i], alternate_quantity=scaling["alternate_quantity"], alternate_units="pb")
                    rows.append(row)
                covariance = np.asarray(coarse["covariance"])
                np.savetxt(base / f"{state}_{observable}_cross_section_covariance.csv",
                           covariance * factors[:, None] * factors[None, :], delimiter=",")
                if alternate is not None:
                    alternate = np.asarray(alternate)
                    np.savetxt(base / f"{state}_{observable}_alternate_cross_section_covariance.csv",
                               covariance * alternate[:, None] * alternate[None, :], delimiter=",")
                write_json(base / f"{state}_{observable}_covariance_metadata.json", common)
        pd.DataFrame(rows).to_csv(base / "differential_results.csv", index=False)
        pd.DataFrame(statuses).to_csv(base / "validation_status.csv", index=False)
    return cfg, states


class ReoptimizationInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cfg, self.entries = make_fixture(self.temp.name)
        self.warning_context = warnings.catch_warnings()
        self.warning_context.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)
        self.addCleanup(self.warning_context.__exit__, None, None, None)
        self.base = self.cfg.reoptimization_run / "artifacts/reoptimization"

    def test_w_no_second_luminosity_or_bin_width_factor(self):
        frame, projected, metadata, _, _ = w_inputs(self.cfg, "1S")
        expected = 120. * 1.02 / (1000. * 2.5 * .0238 * .01)
        self.assertAlmostEqual(frame.sigma_gamma_p_pb.iloc[0], expected)
        self.assertAlmostEqual(frame.sigma_gamma_p_truth_closure_pb.iloc[0], expected / 1.02)
        self.assertAlmostEqual(projected["mean_count"].iloc[0], 24.)
        self.assertFalse(metadata["yield_luminosity_rescaled"])
        self.assertEqual(metadata["overall_status"], "NOT_VALIDATED")
        self.assertEqual(metadata["stage_status"], "PASS_CONDITIONAL")
        self.assertEqual(frame.W_center.iloc[0], 35.)
        self.assertEqual(len(w_inputs(self.cfg, "3S")[0]), 3)

    def test_different_luminosity_requires_new_inference(self):
        with self.assertRaisesRegex(ValueError, "Rerun reoptimization"):
            load_case(replace(self.cfg, lumi_invfb=10.), "1S", "pt2")

    def test_rapidity_uses_alternate_and_full_covariance(self):
        frame, meta, _, entry = load_case(self.cfg, "2S", "y")
        np.testing.assert_allclose(frame["mean"], frame.fitted_truth_cross_section / entry["flux_integral"])
        self.assertEqual(meta["displayed_quantity"], "flux_normalized_ep_rapidity")
        factors = np.asarray(entry["cross_section_scaling"]["alternate_factors"])
        count_cov = np.asarray(entry["inference"]["coarse"]["covariance"])
        np.testing.assert_allclose(meta["covariance"], count_cov*factors[:, None]*factors[None, :])
        self.assertGreater(meta["covariance"][0][1], 0.)
        # Leading and trailing reconstructed control bins must not shift values.
        np.testing.assert_allclose(frame.mean_count, np.linspace(120., 20., 3) / 2. * .2)
        self.assertGreater(frame.lower.iloc[-1], frame["mean"].iloc[-1])

    def test_saved_slope_amplitude_and_profile(self):
        fit, averages, profile = saved_slope(self.cfg, "1S", (0., 1.))
        factor = 1. / (1000. * 2.5 * .0238 * .1012186984099357)
        self.assertAlmostEqual(fit["amplitude_pb_per_GeV2"], 143. * factor)
        self.assertEqual(fit["b_slope_GeVminus2"], 2.7)
        self.assertEqual(fit["b_interval_low"], 1.8)
        self.assertEqual(fit["b_interval_high"], 4.2)
        np.testing.assert_allclose(profile.b, [0., 2.7, 5.])
        np.testing.assert_allclose(profile.delta_deviance, [8., 0., 3.])
        last = averages.iloc[-1]
        expected = 143.*factor*(np.exp(-2.7*.6)-np.exp(-2.7))/(2.7*.4)
        self.assertAlmostEqual(last.fit_bin_average_pb_per_GeV2, expected)

    def test_overall_status_mismatch_is_rejected(self):
        path = self.base / "validation_status.csv"
        table = pd.read_csv(path)
        table.loc[(table.state == "1S") & (table.observable == "w"), "status"] = "PASS_CONDITIONAL"
        table.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "statuses disagree"):
            load_case(self.cfg, "1S", "w")

    def test_missing_new_input_never_falls_back(self):
        (self.base / "3S/y/result.json").unlink()
        with self.assertRaises(FileNotFoundError):
            load_case(self.cfg, "3S", "y")

    def test_covariance_mismatch_is_rejected(self):
        path = self.base / "1S_y_alternate_cross_section_covariance.csv"
        covariance = np.loadtxt(path, delimiter=",")
        covariance[0, 1] *= 2.
        np.savetxt(path, covariance, delimiter=",")
        with self.assertRaisesRegex(ValueError, "full covariance"):
            load_case(self.cfg, "1S", "y")

    def test_control_bin_mapping_mismatch_is_rejected(self):
        path = self.base / "1S/y/result.json"
        entry = json.loads(path.read_text())
        # The number of indices and their uniqueness remain valid, but this
        # shifts the leading reconstructed control bin into the first point.
        n = len(entry["optimization"]["selected_edges"]) - 1
        entry["validation_model"]["coarse"]["metadata"]["fit_reco_indices"] = list(range(n))
        write_json(path, entry)
        with self.assertRaises(ValueError):
            load_case(self.cfg, "1S", "y")

    def test_mixed_slope_normalization_is_not_hidden_as_missing_fit(self):
        from upsilon_esr.differential import plot_differential
        path = self.base / "1S/pt2/result.json"
        entry = json.loads(path.read_text())
        entry["cross_section_factor"] *= 2.
        write_json(path, entry)
        with self.assertRaisesRegex(ValueError, "slope amplitude normalization"):
            plot_differential(self.cfg)
        plt.close("all")

    def test_unavailable_slope_retains_final_status_in_audit_and_csv(self):
        from upsilon_esr.differential import plot_differential
        path = self.base / "3S/pt2/result.json"
        entry = json.loads(path.read_text())
        entry.update(status="NOT_VALIDATED", overall_status="NOT_VALIDATED")
        entry["inference"]["slope"].update(b=None, amplitude=None, status="UNIDENTIFIABLE_NO_SIGNAL")
        write_json(path, entry)
        for name in ("differential_results.csv", "validation_status.csv"):
            path = self.base / name
            table = pd.read_csv(path)
            mask = (table.state == "3S") & (table.observable == "pt2")
            table.loc[mask, ["status", "overall_status"]] = "NOT_VALIDATED"
            table.to_csv(path, index=False)
        path = self.base / "3S_pt2_covariance_metadata.json"
        metadata = json.loads(path.read_text())
        metadata.update(status="NOT_VALIDATED", overall_status="NOT_VALIDATED")
        write_json(path, metadata)
        fig, audit = plot_differential(self.cfg)
        plt.close(fig)
        row = audit["slope_fits"]["3S"]
        self.assertEqual(row["status"], "NOT_VALIDATED")
        self.assertEqual(row["overall_status"], "NOT_VALIDATED")
        self.assertEqual(row["stage_status"], "UNAVAILABLE")
        self.assertFalse(row["publication_ready"])
        slopes = pd.read_csv(self.cfg.output_dir / "upsilon_es_slopes.csv")
        saved = slopes.loc[slopes.state == "3S"].iloc[0]
        self.assertEqual(saved.overall_status, "NOT_VALIDATED")
        self.assertEqual(saved.stage_status, "UNAVAILABLE")

    def test_renderers_do_not_reproject_or_refit(self):
        from upsilon_esr.w_cross_section import plot_w_cross_section
        from upsilon_esr.differential import plot_differential
        with mock.patch("upsilon_esr.w_cross_section.project_cross_section", side_effect=AssertionError("old reproject")), \
             mock.patch("upsilon_esr.differential.project_cross_section", side_effect=AssertionError("old reproject")), \
             mock.patch("upsilon_esr.differential.fit_projected_exponential", side_effect=AssertionError("old refit")):
            fig, audit_w = plot_w_cross_section(self.cfg)
            plt.close(fig)
            fig, audit_d = plot_differential(self.cfg)
            plt.close(fig)
        self.assertFalse(audit_w["publication_ready"])
        self.assertFalse(audit_d["publication_ready"])
        self.assertEqual(audit_d["slope_fits"]["1S"]["b_slope_GeVminus2"], 2.7)
        points = pd.read_csv(self.cfg.output_dir / "upsilon_es_differential_points.csv")
        self.assertEqual(set(points.loc[points.variable == "y", "overall_status"]), {"NOT_VALIDATED"})
        for name in ("upsilon_es_differential", "upsilon_es_w_cross_section"):
            self.assertTrue((self.cfg.output_dir / (name + ".pdf")).is_file())


if __name__ == "__main__":
    unittest.main()
