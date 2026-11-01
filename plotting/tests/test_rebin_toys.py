"""Independent low-count and correlated-bin checks for joint toy replay."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
from scipy.special import xlogy

from upsilon_esr.rebin_toys import replay_joint_toys, poisson_unfold


def _deviance(counts, means):
    if np.any((counts > 0) & (means == 0)):
        return float("inf")
    return max(0., float(2 * np.sum(xlogy(counts, counts) - xlogy(counts, means) - counts + means)))


def analytic_fixture(migrating=False, n_toys=200):
    """Two-bin MLE is analytic: interior inverse or either physical boundary.

    These saved summaries are built without calling the replay solver. The
    migration case produces off-diagonal covariance and boundary solutions.
    """
    response = (np.array([[.65, .18], [.13, .55]]) if migrating
                else np.diag([.25, .5]))
    truth = np.array([3., 8.])
    means = response @ truth
    settings = dict(seed=20261023, n_toys=n_toys, coverage=.6826894921370859)
    draws = np.random.default_rng(np.random.SeedSequence(settings["seed"]).spawn(3)[0]).poisson(
        means, size=(n_toys, len(means)))

    def solve(counts):
        interior = np.linalg.solve(response, counts)
        if np.all(interior >= 0):
            return interior
        total = counts.sum()
        candidates = [np.array([total / response[:, 0].sum(), 0.]),
                      np.array([0., total / response[:, 1].sum()])]
        return min(candidates, key=lambda fitted: _deviance(counts, response @ fitted))

    values = np.asarray([solve(n) for n in draws])
    covariance = np.cov(values, rowvar=False, ddof=1)
    alpha = (1 - settings["coverage"]) / 2
    lower, upper = np.quantile(values, [alpha, 1-alpha], axis=0)
    model = dict(truth_edges=[20., 50., 100.], reco_edges=[20., 50., 100.],
                 response=response.tolist(), feed_in=[0., 0.],
                 truth_expected=truth.tolist(), reco_expected=means.tolist())
    summary = dict(status="PASS_CONDITIONAL", truth_edges=model["truth_edges"],
                   fitted_truth=truth.tolist(), reference_truth=truth.tolist(),
                   reco_expected=means.tolist(), fitted_reco=means.tolist(), asimov_deviance=0.,
                   n_toys=n_toys, toy_status="COMPLETED", successful_toys=n_toys, failed_fraction=0.,
                   toy_mean=values.mean(axis=0).tolist(), covariance=covariance.tolist(),
                   sampling_sd=np.sqrt(np.diag(covariance)).tolist(),
                   toy_lower=lower.tolist(), toy_upper=upper.tolist(),
                   boundary_fractions=(values <= 1e-8).mean(axis=0).tolist(),
                   zero_reco_fractions=(draws == 0).mean(axis=0).tolist())
    entry = dict(state="1S", observable="w", complete=True, overall_status="NOT_VALIDATED",
                 model={"coarse": deepcopy(model)}, validation_model={"coarse": deepcopy(model)},
                 inference=dict(kind="spectrum", settings=settings, coarse=summary))
    return entry, values


class JointToyReplayTests(unittest.TestCase):
    def test_sparse_zero_observations_match_analytic_diagonal_model(self):
        entry, expected = analytic_fixture()
        actual, meta = replay_joint_toys(entry)
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)
        self.assertTrue(np.any(actual == 0))
        self.assertTrue(meta["saved_summary_reproduced"])
        self.assertEqual(meta["overall_status"], "NOT_VALIDATED")
        self.assertFalse(meta["publication_ready"])
        self.assertEqual(meta["failed_toys"], 0)

    def test_migration_and_boundary_solutions_match_independent_analytic_mle(self):
        entry, expected = analytic_fixture(migrating=True)
        actual, _ = replay_joint_toys(entry)
        np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-6)
        covariance = np.cov(actual, rowvar=False)
        self.assertLess(covariance[0, 1], 0.)
        self.assertTrue(np.any((actual[:, 0] == 0) & (actual[:, 1] > 0)))
        # Sum quantiles differ from the sum of marginal quantiles. The joint
        # fitted samples, rather than interval endpoints, are required.
        q = [.15865525393145705, .8413447460685429]
        self.assertFalse(np.allclose(np.quantile(actual.sum(axis=1), q),
                                     np.quantile(actual, q, axis=0).sum(axis=1)))

    def test_cache_reuses_checked_vectors_without_drawing_new_toys(self):
        entry, expected = analytic_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            actual, first = replay_joint_toys(entry, tmp)
            with mock.patch("upsilon_esr.rebin_toys.np.random.default_rng",
                            side_effect=AssertionError("Cached replay must not draw new toys")):
                cached, second = replay_joint_toys(entry, tmp)
            self.assertFalse(first["cache_hit"])
            self.assertTrue(second["cache_hit"])
            np.testing.assert_array_equal(cached, actual)
            self.assertEqual(len(list(Path(tmp).iterdir())), 1)

    def test_corrupt_cache_rejected(self):
        entry, _ = analytic_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            _, meta = replay_joint_toys(entry, tmp)
            cache = Path(meta["cache_path"])
            with np.load(cache, allow_pickle=False) as source:
                values, metadata = source["values"].copy(), source["metadata"].copy()
            values[0, 0] += 1.
            np.savez_compressed(cache, values=values, metadata=metadata)
            with self.assertRaisesRegex(ValueError, "Invalid joint-toy cache"):
                replay_joint_toys(entry, tmp)

    def test_saved_covariance_or_toy_summary_mismatch_fails(self):
        for key in ("covariance", "toy_lower", "toy_mean"):
            entry, _ = analytic_fixture()
            value = np.asarray(entry["inference"]["coarse"][key])
            value.flat[0] += 1.
            entry["inference"]["coarse"][key] = value.tolist()
            with self.assertRaisesRegex(ValueError, "does not reproduce"):
                replay_joint_toys(entry)

    def test_settings_and_unsuccessful_fit_count_cannot_be_ignored(self):
        entry, _ = analytic_fixture()
        entry["inference"]["coarse"]["successful_toys"] -= 1
        with self.assertRaisesRegex(ValueError, "success count"):
            replay_joint_toys(entry)
        entry, _ = analytic_fixture()
        entry["inference"]["settings"]["n_toys"] += 1
        with self.assertRaisesRegex(ValueError, "settings disagree"):
            replay_joint_toys(entry)

    def test_wrong_observable_and_invalid_seed_are_rejected(self):
        entry, _ = analytic_fixture()
        entry["observable"] = "pt2"
        with self.assertRaisesRegex(ValueError, "only for saved W/y"):
            replay_joint_toys(entry)
        entry, _ = analytic_fixture()
        entry["inference"]["settings"]["seed"] = -1
        with self.assertRaisesRegex(ValueError, "seed"):
            replay_joint_toys(entry)


if __name__ == "__main__":
    unittest.main()
