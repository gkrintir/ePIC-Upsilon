"""Recover the saved W/y joint sampling ensemble before merging bins.

Only the original-bin estimator is replayed. No Gaussian approximation, new
response binning, extra coverage validation or new physical uncertainty is
introduced. The two numerical functions below are copied verbatim from the v4
framework scripts/reoptimization/inference.py supplied with this analysis.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import scipy
from scipy.optimize import brentq, nnls
from scipy.special import xlogy

SOLVER_SOURCE_SHA256 = "bdbd2e989594311fdb7311ea342ef6fb0bbb95a1ab71af39a18b347835699c39"
# The original solver stops at a normalized KKT residual of 1e-7.
# These comparisons allow small platform/linear-algebra roundoff, not a new
# sampling ensemble. Counts are in expected-event units; covariance is events².
REPLAY_RTOL = 2e-6
REPLAY_ATOL = 2e-7


def poisson_deviance(n, mu):
    n, mu = np.broadcast_arrays(np.asarray(n, float), np.asarray(mu, float))
    if np.any(n < 0) or np.any(mu < 0) or not np.isfinite(n).all() or not np.isfinite(mu).all():
        raise ValueError("Poisson inputs must be finite and nonnegative")
    if np.any((n > 0) & (mu == 0)):
        return float("inf")
    return float(max(0.0, 2 * np.sum(xlogy(n, n) - xlogy(n, mu) - n + mu)))


def poisson_unfold(response, counts, feed, initial=None):
    """Nonnegative coarse truth-bin MLE, with a fixed response and feed-in."""
    response, counts, feed = map(lambda x: np.asarray(x, float), (response, counts, feed))
    if response.ndim != 2 or counts.shape != (response.shape[0],) or feed.shape != counts.shape:
        raise ValueError("Poisson extraction dimensions mismatch")
    if any(not np.isfinite(x).all() or np.any(x < 0) for x in (response, counts, feed)):
        raise ValueError("Poisson extraction requires finite nonnegative inputs")
    if np.linalg.matrix_rank(response) < response.shape[1]:
        return dict(status="UNIDENTIFIABLE_RESPONSE_RANK", truth=None)
    if np.any((counts > 0) & (response.sum(axis=1) == 0) & (feed == 0)):
        return dict(status="UNSUPPORTED_RECO_COUNTS", truth=None)
    initial = np.asarray(initial, float) if initial is not None else nnls(response, counts - feed)[0]
    if initial.shape != (response.shape[1],) or not np.isfinite(initial).all() or np.any(initial < 0):
        raise ValueError("Invalid initial truth yields")
    truth = initial.copy()
    efficiency = response.sum(axis=0)
    n_total = float(counts.sum())
    # The convex objective is minimized by exact coordinate updates. Each
    # coordinate has a monotone derivative and a finite bracket. This preserves
    # physical zeros without a discontinuous numerical barrier at mu=0.
    max_sweeps = 3000
    residual = float("inf")
    for sweep in range(max_sweeps):
        for j in range(response.shape[1]):
            column = response[:, j]
            active = column > 0
            rest = feed + response @ truth - column * truth[j]
            # Only roundoff cancellation can make a nonnegative sum negative.
            if np.any(rest < -1e-10 * (1 + n_total)):
                raise RuntimeError("Negative Poisson coordinate remainder")
            rest = np.maximum(rest, 0.)
            def derivative(value):
                mu = rest[active] + column[active] * value
                numerator = counts[active] * column[active]
                if np.any((mu == 0) & (numerator > 0)):
                    return -float("inf")
                ratio = np.divide(numerator, mu, out=np.zeros_like(mu), where=mu > 0)
                return float(efficiency[j] - ratio.sum())
            if n_total == 0 or derivative(0.) >= 0:
                truth[j] = 0.
            else:
                upper = n_total / efficiency[j]
                d_upper = derivative(upper)
                if abs(d_upper) <= 1e-12 * efficiency[j]:
                    truth[j] = upper
                else:
                    truth[j] = brentq(derivative, 0., upper,
                                      xtol=1e-11 * (1 + upper), rtol=1e-11)
        mu = response @ truth + feed
        if np.any((counts > 0) & (mu == 0)):
            continue
        ratio = np.divide(counts, mu, out=np.zeros_like(mu), where=mu > 0)
        gradient = response.T @ (1 - ratio)
        projected = np.where(truth > 1e-9, gradient, np.minimum(gradient, 0.))
        residual = float(np.max(np.abs(projected / efficiency)))
        if residual < 1e-7:
            break
    if residual >= 1e-7:
        return dict(status="FIT_FAILED_KKT", truth=None, kkt_residual=residual,
                    iterations=max_sweeps)
    return dict(status="OK_CONDITIONAL", truth=truth.tolist(),
                fitted_reco=mu.tolist(), deviance=poisson_deviance(counts, mu),
                at_boundary=(truth <= 1e-8).tolist(), kkt_residual=residual,
                iterations=sweep + 1)


def _same(actual, expected, name):
    actual, expected = np.asarray(actual, float), np.asarray(expected, float)
    if (actual.shape != expected.shape or not np.isfinite(actual).all()
            or not np.isfinite(expected).all()
            or not np.allclose(actual, expected, rtol=REPLAY_RTOL, atol=REPLAY_ATOL)):
        raise ValueError(f"Joint-toy replay does not reproduce saved {name}. "
                         "The rebin option cannot proceed with incompatible numerical inputs; "
                         "no Gaussian or marginal-interval fallback is applied.")


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"Saved {name} must be an integer.")
    if not np.isfinite(value) or int(value) != value or value < minimum:
        raise ValueError(f"Saved {name} must be an integer >= {minimum}.")
    return int(value)


def _values_hash(values):
    return hashlib.sha256(np.ascontiguousarray(values, dtype="<f8").tobytes()).hexdigest()


def _verify_ensemble(values, failed, n_toys, coverage, saved, draws=None):
    values = np.asarray(values, float)
    n_bins = len(saved["fitted_truth"])
    if (values.ndim != 2 or values.shape[1] != n_bins or len(values) < 2
            or not np.isfinite(values).all() or np.any(values < 0)):
        raise ValueError("Invalid replayed joint-toy fit array.")
    if len(values) + failed != n_toys:
        raise ValueError("Joint-toy replay lost or duplicated toy outcomes.")
    if len(values) != _integer(saved["successful_toys"], "successful_toys", 2):
        raise ValueError("Joint-toy replay success count differs from the saved run.")
    _same(failed / n_toys, saved["failed_fraction"], "toy failure fraction")
    _same(values.mean(axis=0), saved["toy_mean"], "toy mean")
    covariance = np.atleast_2d(np.cov(values, rowvar=False, ddof=1))
    _same(covariance, saved["covariance"], "full joint covariance")
    _same(np.sqrt(np.diag(covariance)), saved["sampling_sd"], "sampling standard deviations")
    alpha = 1. - coverage
    quantiles = np.quantile(values, [alpha / 2., 1. - alpha / 2.], axis=0)
    _same(quantiles[0], saved["toy_lower"], "lower sampling quantiles")
    _same(quantiles[1], saved["toy_upper"], "upper sampling quantiles")
    if "boundary_fractions" in saved:
        _same((values <= 1e-8).mean(axis=0), saved["boundary_fractions"], "boundary fractions")
    if draws is not None and "zero_reco_fractions" in saved:
        _same((draws == 0).mean(axis=0), saved["zero_reco_fractions"], "zero-count fractions")
    return values


def replay_joint_toys(entry, cache_dir=None):
    """Return original-bin fitted truth vectors and replay provenance.

    The RNG uses the original spectrum inference's first spawned stream. All
    successful fit vectors are kept in their original order; every failed fit
    is counted and checked against the stored result. Cached arrays are checked
    against both their content digest and the saved numerical summaries.
    """
    live_solver = "\n\n\n".join(inspect.getsource(fn).rstrip("\n")
                                for fn in (poisson_deviance, poisson_unfold)) + "\n"
    if hashlib.sha256(live_solver.encode()).hexdigest() != SOLVER_SOURCE_SHA256:
        raise ValueError("Poisson replay solver source differs from the recorded v4 implementation.")
    if entry.get("observable") not in ("w", "y") or entry.get("inference", {}).get("kind") != "spectrum":
        raise ValueError("Joint-toy replay is supported only for saved W/y spectrum inference.")
    if not entry.get("complete"):
        raise ValueError("Cannot replay an incomplete inference case.")
    inference = entry["inference"]
    settings, saved = inference["settings"], inference["coarse"]
    seed = _integer(settings["seed"], "inference seed")
    n_toys = _integer(settings["n_toys"], "n_toys", 2)
    if _integer(saved["n_toys"], "coarse n_toys", 2) != n_toys or saved.get("toy_status") != "COMPLETED":
        raise ValueError("Saved spectrum toys are incomplete or their settings disagree.")
    coverage = float(settings["coverage"])
    if not np.isfinite(coverage) or not 0 < coverage < 1:
        raise ValueError("Invalid saved sampling-quantile coverage.")
    train, validation = entry["model"]["coarse"], entry["validation_model"]["coarse"]
    response = np.asarray(train["response"], float)
    feed = np.asarray(train["feed_in"], float)
    counts = np.asarray(validation["reco_expected"], float)
    initial = np.asarray(validation["truth_expected"], float)
    edges = np.asarray(train["truth_edges"], float)
    if (edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all()
            or np.any(np.diff(edges) <= 0) or response.ndim != 2
            or response.shape != (len(counts), len(edges)-1)
            or feed.shape != counts.shape or initial.shape != (len(edges)-1,)):
        raise ValueError("Saved spectrum model has incompatible dimensions.")
    if (any(not np.isfinite(v).all() or np.any(v < 0) for v in (response, feed, counts, initial))
            or np.any(initial <= 0) or np.any(response.sum(axis=0) <= 0)
            or np.any(response.sum(axis=0) > 1 + 1e-9)):
        raise ValueError("Saved spectrum model has invalid probabilities or expectations.")
    _same(edges, validation["truth_edges"], "training/validation truth edges")
    _same(edges, saved["truth_edges"], "extraction truth edges")
    _same(train["reco_edges"], validation["reco_edges"], "training/validation reconstructed edges")
    _same(counts, saved["reco_expected"], "reconstructed expectations")
    _same(initial, saved["reference_truth"], "held-out truth expectation")
    asimov = poisson_unfold(response, counts, feed, initial)
    if asimov.get("truth") is None:
        raise ValueError(f"Cannot reproduce nominal Poisson fit: {asimov.get('status')}")
    _same(asimov["truth"], saved["fitted_truth"], "nominal fitted truth")
    _same(asimov["fitted_reco"], saved["fitted_reco"], "nominal fitted reconstruction")
    _same(asimov["deviance"], saved["asimov_deviance"], "nominal deviance")
    # Summaries are included in the key so stale/mixed result files never reuse
    # a cache solely because their response happens to have the same shape.
    payload = dict(solver=SOLVER_SOURCE_SHA256, response=response.tolist(),
                   feed=feed.tolist(), counts=counts.tolist(), initial=initial.tolist(),
                   truth_edges=edges.tolist(), seed=seed, n_toys=n_toys,
                   coverage=coverage, summary=saved)
    key = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()
    cache_path = None if cache_dir is None else Path(cache_dir) / ("joint-toys-" + key + ".npz")
    cache_hit = cache_path is not None and cache_path.is_file()
    if cache_hit:
        try:
            with np.load(cache_path, allow_pickle=False) as cached:
                values = np.asarray(cached["values"], float)
                cache_metadata = json.loads(str(cached["metadata"].item()))
            if (cache_metadata.get("key") != key or cache_metadata.get("values_sha256") != _values_hash(values)
                    or cache_metadata.get("solver_source_sha256") != SOLVER_SOURCE_SHA256):
                raise ValueError("Cache content identity mismatch")
            failed = _integer(cache_metadata["failed_toys"], "cached failed_toys")
            values = _verify_ensemble(values, failed, n_toys, coverage, saved)
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise ValueError(f"Invalid joint-toy cache {cache_path}. Remove that cache file and retry. {error}") from error
    else:
        label = f"{entry.get('state', '')}/{entry['observable']}"
        print(f"Replaying {n_toys} saved-model Poisson toys for {label} before merging bins...", flush=True)
        rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(3)[0])
        draws = rng.poisson(counts, size=(n_toys, len(counts)))
        fitted, failed = [], 0
        for i, draw in enumerate(draws):
            result = poisson_unfold(response, draw, feed, initial)
            if result.get("truth") is None:
                failed += 1
            else:
                fitted.append(result["truth"])
            if n_toys >= 500 and (i + 1) % 250 == 0:
                print(f"  {label}: {i + 1}/{n_toys} toys, {failed} failed fits", flush=True)
        values = _verify_ensemble(fitted, failed, n_toys, coverage, saved, draws)
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_metadata = dict(key=key, values_sha256=_values_hash(values),
                                  solver_source_sha256=SOLVER_SOURCE_SHA256, failed_toys=failed)
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(mode="wb", dir=cache_path.parent,
                                                 prefix="joint-toys-", suffix=".npz.tmp", delete=False) as stream:
                    temporary = Path(stream.name)
                    np.savez_compressed(stream, values=values,
                                        metadata=np.asarray(json.dumps(cache_metadata, sort_keys=True)))
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary, cache_path)
            finally:
                if temporary is not None and temporary.exists():
                    temporary.unlink()
    metadata = dict(solver_source_sha256=SOLVER_SOURCE_SHA256, seed=seed,
                    rng="default_rng(SeedSequence(seed).spawn(3)[0])", n_toys=n_toys,
                    successful_toys=len(values), failed_toys=failed, failed_fraction=failed/n_toys,
                    saved_summary_reproduced=True, comparison_rtol=REPLAY_RTOL, comparison_atol=REPLAY_ATOL,
                    original_bin_count=len(initial), coverage=coverage,
                    numpy_version=np.__version__, scipy_version=scipy.__version__,
                    cache_hit=cache_hit, cache_key=key,
                    cache_path=str(cache_path) if cache_path is not None else None,
                    values_sha256=_values_hash(values),
                    overall_status=entry.get("overall_status", entry.get("status", "NOT_VALIDATED")),
                    publication_ready=False,
                    interpretation="Replay of original-bin fixed-response/fixed-feed Poisson sampling ensemble; no new coverage validation")
    return values, metadata
