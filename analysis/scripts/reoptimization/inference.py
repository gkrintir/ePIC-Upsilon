"""Conditional low-count inference for the reoptimization appendix.

The likelihood is at reconstructed-count level. The response and feed-in
template are fixed calibration inputs. Reported profile limits are conditional
Wilks forecasts; independent injected-slope toys explicitly test their coverage.
No function treats fractional Asimov expectations as integer observations.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar, nnls
from scipy.special import xlogy
from scipy.stats import chi2


SCOPE = (
    "Signal-only reconstructed Poisson likelihood with fixed response and fixed "
    "out-of-fit-range feed-in. No mass-component extraction, background, "
    "luminosity, branching-ratio or flux uncertainty. Response-MC variation is "
    "reported separately and is not silently added to the count interval."
)


def _settings(config):
    defaults = dict(n_toys=1000, n_coverage_toys=1000, n_response_bootstrap=200,
                    n_model_toys=500, model_tail_probability_min=0.01,
                    min_toys_for_validation=400,
                    seed=20261023, coverage=0.6826894921370859,
                    b_search_max=40.0, b_limit_max=1000.0,
                    coverage_min=0.65, max_bias_over_sd=0.30,
                    max_failed_fraction=0.01, max_boundary_fraction=0.20,
                    feed_in_scales=[0.5, 1.5], profile_points=121)
    defaults.update({k: config[k] for k in defaults if k in config})
    defaults.update(config.get("inference", {}))
    aliases = {"coverage_slopes": "test_slopes", "max_slope": "b_search_max",
               "min_coverage": "coverage_min", "max_failure_fraction": "max_failed_fraction",
               "max_abs_standardized_bias": "max_bias_over_sd"}
    for supplied, canonical in aliases.items():
        if supplied in defaults:
            defaults[canonical] = defaults.pop(supplied)
    if "feed_in_stress_scales" in config:
        defaults["feed_in_scales"] = config["feed_in_stress_scales"]
    for key in ("n_toys", "n_coverage_toys", "n_response_bootstrap", "n_model_toys"):
        if int(defaults[key]) != defaults[key] or defaults[key] < 0:
            raise ValueError(f"{key} must be a nonnegative integer")
        defaults[key] = int(defaults[key])
    if not 0 < defaults["coverage"] < 1:
        raise ValueError("coverage must be between zero and one")
    if defaults["b_search_max"] <= 0 or defaults["b_limit_max"] < defaults["b_search_max"]:
        raise ValueError("Invalid slope search limits")
    return defaults


def exponential_integrals(edges, b):
    """Exact integrals of exp(-b*x), including the b=0 limit."""
    edges = np.asarray(edges, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all():
        raise ValueError("Truth edges must be finite and one-dimensional")
    if np.any(np.diff(edges) <= 0) or edges[0] < 0:
        raise ValueError("Truth edges must be increasing and nonnegative")
    if not np.isfinite(b) or b < 0:
        raise ValueError("Slope must be finite and nonnegative")
    widths = np.diff(edges)
    if b < 1e-10:
        return widths.copy()
    return np.exp(-b * edges[:-1]) * (-np.expm1(-b * widths)) / b


def poisson_deviance(n, mu):
    n, mu = np.broadcast_arrays(np.asarray(n, float), np.asarray(mu, float))
    if np.any(n < 0) or np.any(mu < 0) or not np.isfinite(n).all() or not np.isfinite(mu).all():
        raise ValueError("Poisson inputs must be finite and nonnegative")
    if np.any((n > 0) & (mu == 0)):
        return float("inf")
    return float(max(0.0, 2 * np.sum(xlogy(n, n) - xlogy(n, mu) - n + mu)))


def _validate_model(model):
    edges = np.asarray(model["truth_edges"], float)
    exponential_integrals(edges, 0.)
    response = np.asarray(model["response"], float)
    expected = np.asarray(model["truth_expected"], float)
    reco = np.asarray(model["reco_expected"], float)
    feed = np.asarray(model["feed_in"], float)
    if response.ndim != 2 or response.shape != (len(reco), len(edges) - 1):
        raise ValueError("Response shape does not match truth/reconstructed bins")
    if expected.shape != (response.shape[1],) or feed.shape != reco.shape:
        raise ValueError("Model count dimensions do not match response")
    for name, values in (("response", response), ("truth_expected", expected),
                         ("reco_expected", reco), ("feed_in", feed)):
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError(f"{name} contains invalid or negative values")
    if np.any(response.sum(axis=0) > 1 + 1e-9):
        raise ValueError("Response probability columns exceed one")
    if np.any(expected <= 0):
        raise ValueError("Every fine truth cell must have positive generated MC support")
    unsupported = (response.sum(axis=1) == 0) & (feed == 0) & (reco > 0)
    if np.any(unsupported):
        raise ValueError("Positive reconstructed expectation without response/feed-in support")
    return edges, response, expected, reco, feed


class ForwardModel:
    """A finite-cell response model with a profiled nonnegative signal yield."""

    def __init__(self, model, settings=None):
        self.edges, self.response, self.truth, self.reco, self.feed = _validate_model(model)
        self.settings = _settings({}) if settings is None else settings
        self._cache = {}
        # A response with constant normalized shape cannot identify a slope.
        p0, _ = self.shape(0.)
        p1, _ = self.shape(float(self.settings["b_search_max"]))
        if np.linalg.norm(p0 - p1) < 1e-8:
            raise ValueError("Slope is not identifiable from this response")

    def shape(self, b):
        key = float(b)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        h = self.response @ exponential_integrals(self.edges, key)
        total = float(h.sum())
        if not np.isfinite(total) or total <= 0:
            raise ValueError("Exponential model has no selected response support")
        value = (h / total, total)
        # Do not let many independent toy optimizations make an unbounded cache.
        if len(self._cache) < 10000:
            self._cache[key] = value
        return value

    def profile(self, b, counts, return_parameters=False):
        counts = np.asarray(counts, float)
        if counts.shape != self.reco.shape or np.any(counts < 0) or not np.isfinite(counts).all():
            raise ValueError("Invalid reconstructed counts")
        p, h_total = self.shape(b)
        if np.any((counts > 0) & (p == 0) & (self.feed == 0)):
            return (float("inf"), 0., self.feed.copy(), 0.) if return_parameters else float("inf")
        n_total = float(counts.sum())
        active = p > 0
        if n_total == 0:
            selected_signal = 0.
        elif np.all(self.feed == 0):
            selected_signal = n_total
        else:
            # dNLL/dS = 1 - sum(n*p/(S*p+g)); monotonic in S.
            # Bins with p=0 contain no information about S. The optimum is at
            # S <= N_active=sum(n[p>0]), since p is a normalized signal shape.
            # An occupied-bin feed term can vanish in sparse Poisson toys even
            # when the complete feed template is nonzero. The optimum then
            # lies exactly at the upper endpoint. The subtractive score above
            # can round to a tiny negative value there and break Brent's bracket.
            n_active = math.fsum(counts[active].tolist())
            def derivative(s):
                mu = s * p[active] + self.feed[active]
                numerator = counts[active] * p[active]
                if np.any((mu == 0) & (numerator > 0)):
                    return -float("inf")
                if s > 0 and s >= n_active / 2.:
                    # Exactly the same score, rearranged near the upper end:
                    # (S-N_active)/S + sum((n/S)*g/(S*p+g)).
                    # At S=N_active every term is nonnegative, with exact zero
                    # when feed contributes only in empty bins. Near S=0 use
                    # the ordinary score to avoid cancellation in this form.
                    fraction = np.divide(self.feed[active], mu,
                                         out=np.zeros_like(mu), where=mu > 0)
                    return ((s - n_active) / s
                            + math.fsum((counts[active] / s * fraction).tolist()))
                ratio = np.divide(numerator, mu, out=np.zeros_like(mu), where=mu > 0)
                return 1. - math.fsum(ratio.tolist())
            if n_active == 0 or derivative(0.) >= 0:
                selected_signal = 0.
            else:
                selected_signal = float(brentq(derivative, 0., n_active,
                                               xtol=1e-10, rtol=1e-10))
        mu = selected_signal * p + self.feed
        deviance = poisson_deviance(counts, mu)
        if return_parameters:
            return deviance, selected_signal / h_total, mu, selected_signal
        return deviance

    def fit(self, counts, limits=True):
        counts = np.asarray(counts, float)
        if counts.sum() == 0:
            return dict(status="UNIDENTIFIABLE_ZERO_COUNT", b=None,
                        interval_low=0., interval_high=None,
                        lower_at_boundary=True, upper_unbounded=True,
                        fitted_reco=self.feed.tolist(), expected_signal_reco=0.)
        upper = float(self.settings["b_search_max"])
        # Locate candidate minima before bounded refinement. Arbitrary migrations
        # need not give a globally unimodal one-dimensional profile.
        grid = np.unique(np.r_[np.linspace(0., min(10., upper), 21),
                               np.geomspace(min(10., upper), upper, 9)])
        values = np.array([self.profile(b, counts) for b in grid])
        candidates = list(zip(grid.tolist(), values.tolist()))
        for j in range(1, len(grid) - 1):
            if values[j] <= values[j - 1] and values[j] <= values[j + 1]:
                fit = minimize_scalar(lambda b: self.profile(b, counts),
                                      bounds=(grid[j - 1], grid[j + 1]),
                                      method="bounded", options={"xatol": 2e-6})
                if fit.success and np.isfinite(fit.fun):
                    candidates.append((float(fit.x), float(fit.fun)))
        # Refine edge brackets, whose minima can lie between a boundary and the
        # nearest coarse-grid point.
        edge_brackets = []
        if values[0] <= values[1]:
            edge_brackets.append((grid[0], grid[1]))
        if values[-1] <= values[-2]:
            edge_brackets.append((grid[-2], grid[-1]))
        for a, z in edge_brackets:
            fit = minimize_scalar(lambda b: self.profile(b, counts), bounds=(a, z),
                                  method="bounded", options={"xatol": 2e-6})
            if fit.success and np.isfinite(fit.fun):
                candidates.append((float(fit.x), float(fit.fun)))
        if not any(np.isfinite(v) for _, v in candidates):
            return dict(status="FIT_FAILED", b=None, message="No finite profile value")
        slope, deviance = min(candidates, key=lambda pair: pair[1])
        _, amplitude, mu, selected_signal = self.profile(slope, counts, True)
        if selected_signal <= 1e-8:
            return dict(status="UNIDENTIFIABLE_NO_SIGNAL", b=None,
                        interval_low=0., interval_high=None,
                        lower_at_boundary=True, upper_unbounded=True,
                        fitted_reco=mu.tolist(), expected_signal_reco=selected_signal,
                        deviance=deviance)
        high_boundary = slope >= upper * (1 - 1e-5)
        result = dict(status="SEARCH_BOUNDARY" if high_boundary else "OK_CONDITIONAL",
                      b=slope, amplitude=amplitude,
                      expected_signal_reco=selected_signal, deviance=deviance,
                      ndof=max(int(np.sum(mu > 0)) - 2, 0),
                      fitted_reco=mu.tolist(), observed_asimov_reco=counts.tolist(),
                      at_lower_boundary=slope <= 1e-7, at_search_boundary=high_boundary,
                      interval_interpretation=("Conditional Asimov profile-deviance forecast; "
                          "chi-square threshold, not a toy-calibrated confidence interval."))
        if not limits:
            return result
        delta = float(chi2.ppf(self.settings["coverage"], 1))
        threshold = deviance + delta
        lower_at_boundary = self.profile(0., counts) <= threshold
        lower = 0. if lower_at_boundary else float(brentq(
            lambda b: self.profile(b, counts) - threshold, 0., slope))
        hi = max(slope + 1., 2 * slope, 1.)
        limit = float(self.settings["b_limit_max"])
        while hi < limit and self.profile(hi, counts) < threshold:
            hi = min(hi * 2, limit)
        unbounded = self.profile(hi, counts) < threshold
        upper_limit = None if unbounded or high_boundary else float(brentq(
            lambda b: self.profile(b, counts) - threshold, slope, hi))
        result.update(interval_low=lower, interval_high=upper_limit,
                      lower_at_boundary=lower_at_boundary,
                      upper_unbounded=False,
                      upper_limit_not_found_within_search=unbounded,
                      upper_limit_search_max=limit,
                      upper_search_limited=high_boundary,
                      delta_deviance_threshold=delta,
                      nominal_coverage=float(self.settings["coverage"]))
        return result


def fit_slope(model, counts, config):
    """Fast public fit/profile API, with no toy ensembles or bootstrap."""
    settings = _settings(config)
    forward = ForwardModel(model, settings)
    result = forward.fit(np.asarray(counts, float), limits=True)
    result["profile"] = _profile_table(forward, np.asarray(counts, float), result, settings)
    result["scope"] = SCOPE
    return result


def _profile_table(forward, counts, fit, settings):
    if fit.get("b") is None:
        return []
    high = fit.get("interval_high")
    xmax = max(5., fit["b"] * 2., (high * 1.3 if high is not None else 0.))
    xmax = min(xmax, float(settings["b_search_max"]))
    grid = np.unique(np.r_[np.linspace(0., xmax, int(settings["profile_points"])), fit["b"]])
    return [dict(b=float(b), delta_deviance=max(0., forward.profile(b, counts) - fit["deviance"]))
            for b in grid]


def _model_adequacy(forward, slope, settings, rng):
    n = int(settings["n_model_toys"])
    if n == 0 or slope.get("b") is None:
        return dict(status="NOT_RUN", n_toys=n)
    fitted = np.asarray(slope["fitted_reco"], float)
    draws = rng.poisson(fitted, size=(n, len(fitted)))
    deviances, failures = [], 0
    for counts in draws:
        fit = forward.fit(counts, limits=False)
        if fit.get("deviance") is not None:
            deviances.append(fit["deviance"])
        elif fit["status"] == "UNIDENTIFIABLE_ZERO_COUNT":
            deviances.append(poisson_deviance(counts, forward.feed))
        else:
            failures += 1
    if not deviances:
        return dict(status="FAIL_NO_SUCCESSFUL_TOYS", n_toys=n, failed_fraction=1.)
    tail = float(np.mean(np.asarray(deviances) >= slope["deviance"]))
    passed = (tail >= settings["model_tail_probability_min"]
              and failures / n <= settings["max_failed_fraction"])
    return dict(status="PASS_CONDITIONAL" if passed else "FAIL_VALIDATION",
                n_toys=n, successful_toys=len(deviances), failed_fraction=failures / n,
                heldout_asimov_minimum_deviance=slope["deviance"],
                fitted_model_toy_deviance_quantiles=np.quantile(deviances, [.5, .9, .95, .99]).tolist(),
                expected_template_discrepancy_tail_probability=tail,
                tail_probability_mc_se=math.sqrt(tail * (1 - tail) / len(deviances)),
                interpretation=("Projected-count model-adequacy diagnostic: compare the minimized "
                    "deviance of the held-out expected template with refitted Poisson experiments "
                    "from its fitted exponential. This is not an observed-data goodness-of-fit p value "
                    "and does not establish the mass model or response model as correct."))


def _slope_coverage(forward, validation, asimov, settings, rng):
    n_toys = int(settings["n_coverage_toys"])
    if not n_toys or asimov.get("b") is None:
        return []
    reference_b = float(asimov["b"])
    injected = settings.get("test_slopes")
    if injected is None:
        injected = [max(.5, .5 * reference_b), reference_b,
                    max(reference_b + 1., 1.5 * reference_b)]
    injected = np.unique(np.asarray(injected, float))
    if not np.isfinite(injected).all() or np.any(injected < 0) or np.any(injected >= settings["b_search_max"]):
        raise ValueError("Injected slopes must be finite, nonnegative and below b_search_max")
    response = np.asarray(validation["response"], float)
    total_truth = float(np.asarray(validation["truth_expected"], float).sum())
    feed = np.asarray(validation["feed_in"], float)
    threshold = float(chi2.ppf(settings["coverage"], 1))
    result = []
    for slope in injected:
        integrals = exponential_integrals(forward.edges, float(slope))
        # Keep the held-out fiducial truth normalization fixed while changing its shape.
        true_yield = total_truth * integrals / integrals.sum()
        means = response @ true_yield + feed
        draws = rng.poisson(means, size=(n_toys, len(means)))
        estimates, qvalues, boundary, failures, uninformative = [], [], 0, 0, 0
        for counts in draws:
            fit = forward.fit(counts, limits=False)
            if fit.get("b") is None:
                if fit["status"].startswith("UNIDENTIFIABLE"):
                    # These experiments constrain no slope. Their full-domain
                    # confidence set covers the injected value; no point estimate
                    # enters bias/SD summaries.
                    uninformative += 1
                    qvalues.append(0.)
                else:
                    failures += 1
                continue
            estimates.append(fit["b"])
            boundary += bool(fit["at_lower_boundary"] or fit["at_search_boundary"])
            qvalues.append(max(0., forward.profile(float(slope), counts) - fit["deviance"]))
        estimates = np.asarray(estimates, float)
        successful = len(estimates)
        if not successful:
            result.append(dict(input_b=float(slope), n_toys=n_toys, successful_toys=0,
                               failed_fraction=failures / n_toys,
                               uninformative_fraction=uninformative / n_toys,
                               coverage=uninformative / n_toys,
                               status="FAIL_NO_IDENTIFIABLE_SLOPES"))
            continue
        # Numerical failures count as not covered, making this a conservative
        # lower-bound diagnostic if any occur. Uninformative full-domain sets
        # cover, and are reported separately from numerical failures.
        covered = int(np.sum(np.asarray(qvalues) <= threshold))
        coverage = covered / n_toys
        conditional_coverage = (covered - uninformative) / successful
        se = math.sqrt(coverage * (1 - coverage) / n_toys)
        std = float(np.std(estimates, ddof=1)) if successful > 1 else 0.
        mean = float(np.mean(estimates))
        bias = mean - float(slope)
        bias_over_sd = bias / std if std > 0 else None
        issues = []
        if coverage < settings["coverage_min"]:
            issues.append("coverage_below_configured_minimum")
        if failures / n_toys > settings["max_failed_fraction"]:
            issues.append("fit_failures")
        if boundary / n_toys > settings["max_boundary_fraction"]:
            issues.append("frequent_boundary_solutions")
        if bias_over_sd is None or abs(bias_over_sd) > settings["max_bias_over_sd"]:
            issues.append("mean_bias_relative_to_sampling_sd")
        result.append(dict(
            input_b=float(slope), n_toys=n_toys, successful_toys=successful,
            expected_reco=means.tolist(), expected_truth_total=total_truth,
            coverage=coverage, coverage_mc_se=se,
            coverage_among_successful=conditional_coverage,
            nominal_coverage=float(settings["coverage"]),
            mean_b=mean, median_b=float(np.median(estimates)),
            sampling_sd=std, bias=bias, bias_over_sd=bias_over_sd,
            b_quantiles=np.quantile(estimates, [.1586552539, .5, .8413447461]).tolist(),
            boundary_fraction=boundary / n_toys, failed_fraction=failures / n_toys,
            uninformative_fraction=uninformative / n_toys,
            zero_count_fraction=float(np.mean(draws.sum(axis=1) == 0)),
            likelihood_ratio_quantiles=np.quantile(qvalues, [.5, settings["coverage"], .95]).tolist(),
            status="PASS_CONDITIONAL" if not issues else "FAIL_VALIDATION",
            issues=issues,
            interpretation=("Independent reconstructed-Poisson experiments generated with held-out "
                "response and fixed truth total, fitted using training response/feed-in. Coverage "
                "tests the chi-square profile set at this injected slope only; no global coverage guarantee.")))
    return result


def _response_bootstrap(model, counts, settings, rng):
    n = int(settings["n_response_bootstrap"])
    if not n:
        return dict(status="NOT_RUN", n_bootstrap=0)
    raw = model.get("raw_categories")
    if raw is None:
        return dict(status="NOT_VERIFIED", reason="No raw response categories supplied")
    raw = np.asarray(raw, float)
    response = np.asarray(model["response"], float)
    if raw.shape != (response.shape[0] + 1, response.shape[1]):
        raise ValueError("raw_categories must have reconstructed rows followed by one missed row")
    if not np.isfinite(raw).all() or np.any(raw < 0) or not np.allclose(raw, np.round(raw)):
        raise ValueError("Response bootstrap requires raw nonnegative integer category counts")
    totals = raw.sum(axis=0).astype(int)
    if np.any(totals == 0):
        return dict(status="FAIL_EMPTY_MC_CELL", n_bootstrap=n)
    probabilities = raw / totals
    if not np.allclose(probabilities[:-1], response, atol=1e-10, rtol=1e-8):
        raise ValueError("Raw multinomial counts disagree with response; weighted/stressed bootstrap is not supported")
    values, failed = [], 0
    for _ in range(n):
        boot = np.column_stack([rng.multinomial(int(total), probabilities[:, j])
                                for j, total in enumerate(totals)])
        candidate = dict(model, response=boot[:-1] / totals)
        try:
            fit = ForwardModel(candidate, settings).fit(counts, limits=False)
        except ValueError:
            failed += 1
            continue
        if fit.get("b") is None or fit.get("at_search_boundary"):
            failed += 1
        else:
            values.append(fit["b"])
    if len(values) < 2:
        return dict(status="FAIL_INSUFFICIENT_BOOTSTRAPS", n_bootstrap=n,
                    successful=len(values), failed_fraction=failed / n)
    return dict(status="PASS_CONDITIONAL" if failed / n <= settings["max_failed_fraction"] else "FAIL_VALIDATION",
                n_bootstrap=n, successful=len(values), failed_fraction=failed / n,
                b_mean=float(np.mean(values)), b_std=float(np.std(values, ddof=1)),
                b_quantiles=np.quantile(values, [.1586552539, .5, .8413447461]).tolist(),
                interpretation=("Multinomial response-MC resampling at fixed generated truth-cell counts, "
                    "including the missed-event category. Data expectation and feed-in template held fixed. "
                    "Separate calibration sensitivity, not combined with the projected count interval. "
                    "Does not include feed-in MC uncertainty or within-cell response variation. "
                    "An empirical multinomial bootstrap cannot populate categories with zero "
                    "observed MC events, so it does not bound unseen migration probabilities."))


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


def _coarse_analysis(model, validation, settings, rng):
    response = np.asarray(model["response"], float)
    ref = np.asarray(validation["truth_expected"], float)
    counts = np.asarray(validation["reco_expected"], float)
    feed = np.asarray(model["feed_in"], float)
    fit = poisson_unfold(response, counts, feed, ref)
    if fit.get("truth") is None:
        return fit
    result = dict(status="NOT_VALIDATED", truth_edges=np.asarray(model["truth_edges"]).tolist(),
                  fitted_truth=fit["truth"], reference_truth=ref.tolist(),
                  reco_expected=counts.tolist(), fitted_reco=fit["fitted_reco"],
                  asimov_deviance=fit["deviance"], asimov_boundary=fit["at_boundary"],
                  asimov_relative_bias=(np.asarray(fit["truth"]) / ref - 1).tolist(),
                  response_condition_number=float(np.linalg.cond(response)),
                  interval_interpretation=("Pointwise central sampling quantiles of the nonnegative "
                     "Poisson truth-yield MLE under fixed response/feed-in. Not confidence intervals."))
    n = int(settings["n_toys"])
    if not n:
        result["toy_status"] = "NOT_RUN"
        return result
    draws = rng.poisson(counts, size=(n, len(counts)))
    values, failed = [], 0
    for draw in draws:
        toy = poisson_unfold(response, draw, feed, ref)
        if toy.get("truth") is None:
            failed += 1
        else:
            values.append(toy["truth"])
    if len(values) < 2:
        result.update(status="FAIL_INSUFFICIENT_TOYS", n_toys=n, failed_fraction=failed / n)
        return result
    values = np.asarray(values)
    mean = values.mean(axis=0)
    covariance = np.atleast_2d(np.cov(values, rowvar=False, ddof=1))
    sd = np.sqrt(np.diag(covariance))
    bias = mean - ref
    bias_over_sd = np.divide(bias, sd, out=np.full_like(bias, np.inf), where=sd > 0)
    alpha = 1 - settings["coverage"]
    quantiles = np.quantile(values, [alpha / 2, 1 - alpha / 2], axis=0)
    problems = []
    if failed / n > settings["max_failed_fraction"]:
        problems.append("fit_failures")
    if np.any(np.abs(bias_over_sd) > settings["max_bias_over_sd"]):
        problems.append("mean_bias_relative_to_sampling_sd")
    result.update(status="PASS_CONDITIONAL" if not problems else "FAIL_VALIDATION",
                  toy_status="COMPLETED", n_toys=n, successful_toys=len(values),
                  failed_fraction=failed / n, toy_mean=mean.tolist(),
                  covariance=covariance.tolist(), toy_lower=quantiles[0].tolist(),
                  toy_upper=quantiles[1].tolist(), sampling_sd=sd.tolist(),
                  bias=bias.tolist(), bias_over_sd=bias_over_sd.tolist(),
                  boundary_fractions=(values <= 1e-8).mean(axis=0).tolist(),
                  zero_reco_fractions=(draws == 0).mean(axis=0).tolist(), issues=problems)
    return result


def analyze(model, validation_model, config):
    """Run independent held-out closure and explicitly conditional statistical tests.

    Optional ``model['coarse']`` and validation counterpart enable binwise
    nonparametric Poisson extraction. Output is JSON serializable.
    """
    settings = _settings(config)
    edges, response, truth, _, feed = _validate_model(model)
    vedges, vresponse, vtruth, vcounts, vfeed = _validate_model(validation_model)
    if response.shape != vresponse.shape or not np.array_equal(edges, vedges):
        raise ValueError("Training and held-out model geometry differ")
    seeds = np.random.SeedSequence(int(settings["seed"])).spawn(4)
    coverage_rng, mc_rng, coarse_rng, model_rng = [np.random.default_rng(s) for s in seeds]
    try:
        forward = ForwardModel(model, settings)
    except ValueError as error:
        return dict(status="FAIL_UNIDENTIFIABLE", scope=SCOPE, reason=str(error))
    slope = forward.fit(vcounts, limits=True)
    slope["profile"] = _profile_table(forward, vcounts, slope, settings)
    slope["fine_truth_edges"] = edges.tolist()
    slope["fit_input"] = "Independent held-out reconstructed expected counts"
    slope["fixed_feed_in"] = feed.tolist()
    slope["heldout_feed_in"] = vfeed.tolist()
    slope["feed_in_fraction"] = float(vfeed.sum() / vcounts.sum()) if vcounts.sum() else None
    slope["zero_efficiency_fine_truth_cells"] = np.flatnonzero(response.sum(axis=0) == 0).tolist()
    if slope.get("b") is not None:
        threshold = slope["delta_deviance_threshold"]
        membership = np.array([p["delta_deviance"] <= threshold for p in slope["profile"]])
        components = int(np.sum(membership & ~np.r_[False, membership[:-1]]))
        slope["profile_set_components_on_display_grid"] = components
        if components > 1:
            slope["status"] = "DISCONNECTED_PROFILE_SET"
    adequacy = _model_adequacy(forward, slope, settings, model_rng)
    coverage = _slope_coverage(forward, validation_model, slope, settings, coverage_rng)
    bootstrap = _response_bootstrap(model, vcounts, settings, mc_rng)
    stress = []
    if slope.get("b") is not None:
        signal = np.maximum(vcounts - vfeed, 0)
        for factor in settings["feed_in_scales"]:
            if factor < 0 or not np.isfinite(factor):
                raise ValueError("Feed-in stress scales must be finite and nonnegative")
            fit = forward.fit(signal + float(factor) * vfeed, limits=False)
            stress.append(dict(scale=float(factor), status=fit["status"], b=fit.get("b"),
                               b_shift=None if fit.get("b") is None else fit["b"] - slope["b"],
                               interpretation="Held-out feed-in varied, training feed-in fixed in fit"))
    coarse_train = model.get("coarse", model.get("coarse_model"))
    coarse_validation = validation_model.get("coarse", validation_model.get("coarse_model"))
    coarse = (_coarse_analysis(coarse_train, coarse_validation, settings, coarse_rng)
              if coarse_train is not None and coarse_validation is not None
              else dict(status="NOT_RUN", reason="Coarse models not supplied"))
    gates = dict(
        slope_identifiable=slope.get("status") == "OK_CONDITIONAL",
        slope_upper_limit_finite=slope.get("interval_high") is not None,
        injected_slope_coverage=bool(coverage) and all(r["status"] == "PASS_CONDITIONAL" for r in coverage),
        response_mc_completed=bootstrap.get("status") == "PASS_CONDITIONAL",
        coarse_extraction=coarse.get("status") == "PASS_CONDITIONAL",
        expected_template_model_adequacy=adequacy.get("status") == "PASS_CONDITIONAL",
        sufficient_validation_toys=(settings["n_coverage_toys"] >= settings["min_toys_for_validation"]
                                    and settings["n_toys"] >= settings["min_toys_for_validation"]),
    )
    return dict(status="PASS_CONDITIONAL" if all(gates.values()) else "NOT_VALIDATED",
                scope=SCOPE, gates=gates, settings=settings, slope=slope,
                slope_toys=coverage, feed_in_stress=stress, response_mc=bootstrap,
                coarse=coarse, model_adequacy=adequacy,
                unresolved=["Per-kinematic-bin mass-component separation and background",
                            "Feed-in normalization/shape uncertainty is stressed, not constrained",
                            "Finite-cell response and truth-shape model dependence",
                            "Statistical coverage established only at tested injected slopes",
                            "Flux, luminosity, branching-ratio and other systematics"])


SPECTRUM_SCOPE = (
    "Nonparametric reconstructed-Poisson spectrum inference with nonnegative truth-bin "
    "yields, fixed training response and fixed outside-fit-range feed-in. Pointwise "
    "profile limits use a chi-square likelihood-ratio threshold and are conditional "
    "forecasts, not a calibrated confidence belt or simultaneous band. Independent "
    "held-out nominal and shape-only injections test coverage at selected spectra. "
    "No exponential slope is fitted. Response-MC variation is reported separately; "
    "background, mass-component separation and systematic uncertainties are excluded."
)


def _spectrum_model(model):
    """Validate a coarse spectrum without restricting the coordinate to positive x."""
    current = model.get("coarse", model)
    edges = np.asarray(current["truth_edges"], float)
    response = np.asarray(current["response"], float)
    truth = np.asarray(current["truth_expected"], float)
    counts = np.asarray(current["reco_expected"], float)
    feed = np.asarray(current["feed_in"], float)
    if (edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all()
            or np.any(np.diff(edges) <= 0)):
        raise ValueError("Spectrum truth edges must be finite and strictly increasing")
    if (response.ndim != 2 or response.shape != (len(counts), len(edges) - 1)
            or truth.shape != (len(edges) - 1,) or feed.shape != counts.shape):
        raise ValueError("Spectrum model dimensions disagree")
    for name, value in (("response", response), ("truth", truth), ("reco", counts), ("feed", feed)):
        if not np.isfinite(value).all() or np.any(value < 0):
            raise ValueError("Spectrum " + name + " must be finite and nonnegative")
    if np.any(truth <= 0):
        raise ValueError("Every fitted spectrum bin needs positive generated truth support")
    if np.any(response.sum(axis=0) > 1 + 1e-9):
        raise ValueError("Spectrum response column sums exceed probability one")
    return current


def profile_truth_bin(response, counts, feed, bin_index, value, initial=None):
    """Profile other nonnegative truth yields at one fixed truth-bin yield.

    An impossible fixed value has infinite deviance, not a numerical-fit failure.
    Zero reconstructed counts and single-truth-bin models are supported.
    """
    response, counts, feed = (np.asarray(x, float) for x in (response, counts, feed))
    if (response.ndim != 2 or counts.shape != (response.shape[0],)
            or feed.shape != counts.shape):
        raise ValueError("Profile-bin dimensions mismatch")
    if not isinstance(bin_index, (int, np.integer)) or not 0 <= bin_index < response.shape[1]:
        raise ValueError("Profile bin index is outside the response")
    if (not np.isfinite(value) or value < 0
            or any(not np.isfinite(x).all() or np.any(x < 0) for x in (response, counts, feed))):
        raise ValueError("Profile inputs must be finite and nonnegative")
    free = np.arange(response.shape[1]) != bin_index
    fixed_feed = feed + float(value) * response[:, bin_index]
    if not free.any():
        return dict(status="OK_CONDITIONAL", deviance=poisson_deviance(counts, fixed_feed),
                    truth=[float(value)], fitted_reco=fixed_feed.tolist())
    free_response = response[:, free]
    if np.any((counts > 0) & (free_response.sum(axis=1) == 0) & (fixed_feed == 0)):
        return dict(status="OK_CONDITIONAL", deviance=float("inf"),
                    truth=None, excluded_fixed_value=True)
    free_initial = None if initial is None else np.asarray(initial, float)[free]
    fitted = poisson_unfold(free_response, counts, fixed_feed, free_initial)
    if fitted.get("truth") is None:
        return fitted
    truth = np.empty(response.shape[1], float)
    truth[bin_index], truth[free] = float(value), np.asarray(fitted["truth"], float)
    fitted["truth"] = truth.tolist()
    return fitted


def spectrum_bin_intervals(model, counts, config, fitted=None):
    """Pointwise conditional profile-deviance intervals for arbitrary truth bins."""
    settings = _settings(config)
    current = _spectrum_model(model)
    response = np.asarray(current["response"], float)
    feed = np.asarray(current["feed_in"], float)
    counts = np.asarray(counts, float)
    fitted = fitted or poisson_unfold(response, counts, feed, current["truth_expected"])
    if fitted.get("truth") is None:
        return []
    point = np.asarray(fitted["truth"], float)
    minimum = float(fitted["deviance"])
    threshold = float(chi2.ppf(settings["coverage"], 1))
    result = []
    for j, estimate in enumerate(point):
        def difference(value):
            item = profile_truth_bin(response, counts, feed, j, float(value), point)
            if item.get("deviance") is None:
                raise RuntimeError("Profile nuisance fit failed: " + item["status"])
            return float(item["deviance"]) - minimum - threshold
        row = dict(bin=j, truth_yield=float(estimate), nominal_coverage=float(settings["coverage"]),
                   delta_deviance_threshold=threshold,
                   interpretation=("Pointwise conditional profile-deviance forecast with other "
                      "truth-bin yields profiled nonnegative; chi-square threshold, not a "
                      "toy-calibrated confidence interval or simultaneous band."))
        try:
            at_zero = difference(0.)
            low = 0. if estimate == 0 or at_zero <= 0 else float(
                brentq(difference, 0., float(estimate), xtol=1e-8, rtol=1e-10))
            efficiency = float(response[:, j].sum())
            if efficiency <= 0:
                raise RuntimeError("No selected support for the profiled truth bin")
            high = float(estimate + max(1., math.sqrt(float(counts.sum()) + 1.) / efficiency))
            for _ in range(50):
                if difference(high) >= 0:
                    break
                high = 2 * high + 1.
            else:
                raise RuntimeError("No finite profile upper endpoint found")
            high = float(brentq(difference, float(estimate), high, xtol=1e-8, rtol=1e-10))
            row.update(status="OK_CONDITIONAL", interval_low=low, interval_high=high,
                       lower_at_boundary=bool(low == 0.), upper_unbounded=False)
        except (ValueError, RuntimeError) as error:
            row.update(status="PROFILE_FAILED", interval_low=None, interval_high=None,
                       reason=str(error))
        result.append(row)
    return result


def _spectrum_injections(validation, coarse, settings):
    """Use held-out fine responses for shape variations, keeping the tail fixed.

    Tilts are dimensionless in z=(x-lo)/(hi-lo), so W and signed rapidity use
    the same bounded-strength deformation. No exponential physics model is
    imposed on the nominal spectrum or its fitted bin yields.
    """
    tilts = np.unique(np.asarray(settings.get("spectrum_coverage_tilts", [-1., 0., 1.]), float))
    if not len(tilts) or not np.isfinite(tilts).all():
        raise ValueError("spectrum_coverage_tilts must be a nonempty finite vector")
    if not np.any(tilts == 0):
        tilts = np.unique(np.r_[tilts, 0.])
    fine_edges = np.asarray(validation["truth_edges"], float)
    coarse_edges = np.asarray(coarse["truth_edges"], float)
    if (fine_edges.ndim != 1 or len(fine_edges) < 2 or np.any(np.diff(fine_edges) <= 0)
            or not np.isfinite(fine_edges).all()
            or not np.allclose(fine_edges[[0, -1]], coarse_edges[[0, -1]], rtol=0, atol=1e-10)
            or any(not np.any(np.isclose(fine_edges, edge, rtol=0, atol=1e-10))
                   for edge in coarse_edges)):
        raise ValueError("Fine held-out truth cells must subdivide the coarse fitted range")
    fine_truth = np.asarray(validation["truth_expected"], float)
    fine_response = np.asarray(validation["response"], float)
    feed = np.asarray(validation["feed_in"], float)
    if (fine_truth.shape != (len(fine_edges) - 1,)
            or fine_response.shape != (len(feed), len(fine_truth))
            or np.any(fine_truth < 0) or np.any(fine_response < 0)
            or not np.isfinite(fine_truth).all() or not np.isfinite(fine_response).all()
            or not np.isfinite(feed).all() or np.any(feed < 0)
            or fine_truth.sum() <= 0):
        raise ValueError("Invalid fine held-out response for spectrum coverage")
    centers = (fine_edges[1:] + fine_edges[:-1]) / 2
    indices = np.searchsorted(coarse_edges, centers, side="right") - 1
    z = (centers - coarse_edges[0]) / (coarse_edges[-1] - coarse_edges[0])
    result = []
    for tilt in tilts:
        if tilt == 0:
            truth = np.asarray(coarse["truth_expected"], float)
            means = np.asarray(coarse["reco_expected"], float)
        else:
            exponent = -float(tilt) * z
            shifted = fine_truth * np.exp(exponent - exponent.max())
            shifted *= fine_truth.sum() / shifted.sum()
            truth = np.bincount(indices, weights=shifted, minlength=len(coarse_edges)-1)
            means = fine_response @ shifted + feed
        result.append(dict(injection="nominal" if tilt == 0 else "shape_only_tilt_%g" % tilt,
                           tilt_dimensionless=float(tilt), truth=truth, means=means))
    return result


def _spectrum_coverage(model, validation, fine_validation, settings, rng):
    n = int(settings["n_coverage_toys"])
    if not n:
        return []
    response, feed = np.asarray(model["response"], float), np.asarray(model["feed_in"], float)
    nbins = response.shape[1]
    threshold = float(chi2.ppf(settings["coverage"], 1))
    rows = []
    for injection in _spectrum_injections(fine_validation, validation, settings):
        reference, means = injection["truth"], injection["means"]
        draws = rng.poisson(means, size=(n, len(means)))
        estimates, qvalues = [[] for _ in range(nbins)], [[] for _ in range(nbins)]
        failures = np.zeros(nbins, dtype=int)
        for draw in draws:
            fit = poisson_unfold(response, draw, feed, reference)
            if fit.get("truth") is None:
                failures += 1
                continue
            for j in range(nbins):
                profile = profile_truth_bin(response, draw, feed, j, float(reference[j]), fit["truth"])
                if profile.get("deviance") is None:
                    failures[j] += 1
                    continue
                estimates[j].append(fit["truth"][j])
                qvalues[j].append(max(0., profile["deviance"] - fit["deviance"]))
        for j in range(nbins):
            values, ratios = np.asarray(estimates[j], float), np.asarray(qvalues[j], float)
            success = len(values)
            coverage = float(np.sum(ratios <= threshold) / n)
            std = float(values.std(ddof=1)) if success > 1 else 0.
            bias = float(values.mean() - reference[j]) if success else None
            standardized = bias / std if std > 0 and bias is not None else None
            boundary = float(np.sum(values <= 1e-8) / n)
            issues = []
            if success < 2:
                issues.append("insufficient_successful_toys")
            if n < settings["min_toys_for_validation"]:
                issues.append("insufficient_validation_toys")
            if coverage < settings["coverage_min"]:
                issues.append("coverage_below_configured_minimum")
            if failures[j] / n > settings["max_failed_fraction"]:
                issues.append("fit_failures")
            if reference[j] > 0 and boundary > settings["max_boundary_fraction"]:
                issues.append("frequent_boundary_solutions")
            if standardized is None or abs(standardized) > settings["max_bias_over_sd"]:
                issues.append("mean_bias_relative_to_sampling_sd")
            rows.append(dict(
                injection=injection["injection"], tilt_dimensionless=injection["tilt_dimensionless"],
                bin=j, injected_truth=float(reference[j]), expected_reco=means.tolist(),
                n_toys=n, successful_toys=success, failed_fraction=float(failures[j] / n),
                coverage=coverage, coverage_mc_se=math.sqrt(coverage * (1-coverage) / n),
                nominal_coverage=float(settings["coverage"]), coverage_min=float(settings["coverage_min"]),
                delta_deviance_threshold=threshold,
                mean=float(values.mean()) if success else None,
                median=float(np.median(values)) if success else None,
                sampling_sd=std, bias=bias, bias_over_sd=standardized,
                boundary_fraction=boundary,
                zero_total_count_fraction=float(np.mean(draws.sum(axis=1) == 0)),
                status="PASS_CONDITIONAL" if not issues else "FAIL_VALIDATION", issues=issues,
                interpretation=("Pointwise profile-set coverage at this held-out injected truth "
                    "spectrum. Other truth-bin yields are profiled, response/feed are fixed to "
                    "training. Shape tilts preserve in-range truth total and hold the held-out "
                    "tail fixed. Fine cells with no held-out MC remain empty; this test does "
                    "not establish support in unpopulated phase space. Numerical failures "
                    "count as not covered. This does not provide "
                    "simultaneous coverage or a calibrated confidence belt.")))
    return rows


def _spectrum_response_bootstrap(model, counts, settings, rng):
    n = int(settings["n_response_bootstrap"])
    if not n:
        return dict(status="NOT_RUN", n_bootstrap=0)
    raw = model.get("raw_categories")
    if raw is None:
        return dict(status="NOT_VERIFIED", reason="No raw response categories supplied")
    raw, response = np.asarray(raw, float), np.asarray(model["response"], float)
    if (raw.shape != (response.shape[0] + 1, response.shape[1])
            or not np.isfinite(raw).all() or np.any(raw < 0)
            or not np.allclose(raw, np.round(raw))):
        raise ValueError("Spectrum bootstrap needs integer reco-plus-missed raw categories")
    totals = raw.sum(axis=0).astype(int)
    if np.any(totals == 0):
        return dict(status="FAIL_EMPTY_MC_CELL", n_bootstrap=n)
    probabilities = raw / totals
    if not np.allclose(probabilities[:-1], response, atol=1e-10, rtol=1e-8):
        raise ValueError("Spectrum multinomial bootstrap requires unweighted raw response probabilities")
    values, failed = [], 0
    for _ in range(n):
        boot = np.column_stack([rng.multinomial(int(total), probabilities[:, j])
                                for j, total in enumerate(totals)])
        fit = poisson_unfold(boot[:-1] / totals, counts, model["feed_in"], model["truth_expected"])
        if fit.get("truth") is None:
            failed += 1
        else:
            values.append(fit["truth"])
    if len(values) < 2:
        return dict(status="FAIL_INSUFFICIENT_BOOTSTRAPS", n_bootstrap=n,
                    successful=len(values), failed_fraction=failed/n)
    values = np.asarray(values, float)
    covariance = np.atleast_2d(np.cov(values, rowvar=False, ddof=1))
    return dict(status="PASS_CONDITIONAL" if failed/n <= settings["max_failed_fraction"] else "FAIL_VALIDATION",
                n_bootstrap=n, successful=len(values), failed_fraction=failed/n,
                truth_mean=values.mean(axis=0).tolist(), covariance=covariance.tolist(),
                sampling_sd=np.sqrt(np.diag(covariance)).tolist(),
                truth_quantiles=np.quantile(values, [.1586552539, .5, .8413447461], axis=0).tolist(),
                interpretation=("Separate multinomial response-MC variation of truth-bin yields "
                    "at fixed generated counts, including missed events. Data means/feed are held "
                    "fixed. Not combined with projected-count intervals; excludes feed-in MC "
                    "uncertainty and unobserved migration categories."))


def analyze_spectrum(model, validation_model, config):
    """Validate W or rapidity bins without fitting an exponential slope.

    The pointwise confidence-set approximation and estimator sampling bands
    are stored separately. No missing toy or failed stage becomes a PASS.
    """
    settings = _settings(config)
    train, validation = _spectrum_model(model), _spectrum_model(validation_model)
    if (np.shape(train["response"]) != np.shape(validation["response"])
            or not np.array_equal(train["truth_edges"], validation["truth_edges"])
            or ("reco_edges" in train and "reco_edges" in validation
                and not np.array_equal(train["reco_edges"], validation["reco_edges"]))):
        raise ValueError("Training and held-out spectrum geometry differs")
    seeds = np.random.SeedSequence(int(settings["seed"])).spawn(3)
    coarse_rng, coverage_rng, mc_rng = [np.random.default_rng(seed) for seed in seeds]
    coarse = _coarse_analysis(train, validation, settings, coarse_rng)
    if coarse.get("fitted_truth") is None:
        return dict(kind="spectrum", status="NOT_VALIDATED", scope=SPECTRUM_SCOPE,
                    settings=settings, coarse=coarse, gates={"coarse_extraction": False},
                    bin_intervals=[], bin_coverage=[], response_mc={"status": "NOT_RUN"})
    fit = dict(truth=coarse["fitted_truth"], deviance=coarse["asimov_deviance"])
    intervals = spectrum_bin_intervals(train, validation["reco_expected"], config, fit)
    coverage = _spectrum_coverage(train, validation, validation_model, settings, coverage_rng)
    bootstrap = _spectrum_response_bootstrap(train, validation["reco_expected"], settings, mc_rng)
    gates = dict(
        coarse_extraction=coarse.get("status") == "PASS_CONDITIONAL",
        profile_intervals=bool(intervals) and all(x["status"] == "OK_CONDITIONAL" for x in intervals),
        injected_spectrum_coverage=bool(coverage) and all(x["status"] == "PASS_CONDITIONAL" for x in coverage),
        response_mc_completed=bootstrap.get("status") == "PASS_CONDITIONAL",
        sufficient_validation_toys=(settings["n_toys"] >= settings["min_toys_for_validation"]
                                    and settings["n_coverage_toys"] >= settings["min_toys_for_validation"]),
    )
    return dict(kind="spectrum", status="PASS_CONDITIONAL" if all(gates.values()) else "NOT_VALIDATED",
                scope=SPECTRUM_SCOPE, settings=settings, gates=gates, coarse=coarse,
                bin_intervals=intervals, bin_coverage=coverage, response_mc=bootstrap,
                unresolved=["Per-bin signal extraction, backgrounds, and truth-fiducial conditioning",
                            "Tail normalization/shape uncertainty is stressed, not constrained",
                            "Pointwise coverage checked only at selected nominal/shape-only injections",
                            "Profile limits are not toy-calibrated confidence belts or simultaneous bands",
                            "Luminosity, branching fraction, flux and other systematics"])
