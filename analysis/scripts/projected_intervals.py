"""Count confidence intervals and projected interval summaries.

Observed counts use the central (equal-tailed) Garwood construction. Projected
means need a forecast, not an invented fractional observation: the default
uses ensemble-median Garwood endpoints, evaluated at the median integer count
of Poisson(mean). The expected-mean point is not rounded. These displayed
interval widths never replace variances in a covariance matrix.
"""
import math

import numpy as np
from scipy.stats import gamma, norm, poisson

DEFAULT_COVERAGE = math.erf(1.0 / math.sqrt(2.0))


def validate_coverage(coverage):
    coverage = float(coverage)
    if not np.isfinite(coverage) or not 0.0 < coverage < 1.0:
        raise ValueError("Interval coverage must be finite and strictly between zero and one")
    return coverage


def add_interval_arguments(parser):
    parser.add_argument('--count-interval', choices=['median_garwood', 'standard_deviation'],
                        default='median_garwood',
                        help='Projected counts: ensemble-median expected Garwood limits, or legacy symmetric SD bands. Does not change covariance.')
    parser.add_argument('--interval-coverage', type=float, default=DEFAULT_COVERAGE,
                        help='Central probability for count intervals and projected sampling bands (default 0.6826894921)')


def garwood_intervals(observed, coverage=DEFAULT_COVERAGE):
    """Exact central Poisson confidence limits for nonnegative integer counts.

    This matches the central TH1 kPoisson convention, including n=0 upper
    bound -log(alpha/2), rather than the one-sided RooHistError convention.
    Discreteness makes coverage conservative. Weighted counts are rejected.
    """
    coverage = validate_coverage(coverage)
    n = np.asarray(observed, dtype=float)
    if np.any(~np.isfinite(n)) or np.any(n < 0) or np.any(n != np.floor(n)):
        raise ValueError('Garwood intervals require nonnegative integer observed counts')
    alpha = 1.0 - coverage
    lower = np.zeros_like(n)
    positive = n > 0
    lower[positive] = gamma.ppf(alpha / 2.0, n[positive])
    upper = gamma.isf(alpha / 2.0, n + 1.0)
    return lower, upper


def expected_count_intervals(means, coverage=DEFAULT_COVERAGE, method='median_garwood'):
    """Forecast bounds without interpreting a fractional mean as observation.

    Monotone Garwood endpoints at the median Poisson count are the respective
    ensemble-median confidence limits. They vary stepwise with the expected
    count. They are not an exact confidence interval computed on fractional
    Asimov data, and their widths are not standard deviations.
    """
    coverage = validate_coverage(coverage)
    mu = np.asarray(means, dtype=float)
    if np.any(~np.isfinite(mu)) or np.any(mu < 0):
        raise ValueError('Projected count means must be finite and nonnegative')
    if method == 'median_garwood':
        reference = poisson.ppf(0.5, mu)
        lower, upper = garwood_intervals(reference, coverage)
    elif method == 'standard_deviation':
        reference = np.full_like(mu, np.nan)
        halfwidth = norm.ppf((1.0 + coverage) / 2.0) * np.sqrt(mu)
        lower, upper = mu - halfwidth, mu + halfwidth
    else:
        raise ValueError('Unknown count interval method: ' + str(method))
    return {'lower': lower, 'upper': upper, 'reference_count': reference,
            'method': method, 'coverage': coverage}


def sampling_interval(toys, coverage=DEFAULT_COVERAGE, axis=0):
    """Absolute central sampling bounds; these are not confidence limits."""
    coverage = validate_coverage(coverage)
    values = np.asarray(toys, dtype=float)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError('Finite nonempty toys are required for a sampling interval')
    tail = (1.0 - coverage) / 2.0
    bounds = np.quantile(values, [tail, 1.0 - tail], axis=axis)
    return bounds[0], bounds[1]


def draw_absolute_intervals(ax, x, y, lower, upper, **kwargs):
    """Draw absolute interval endpoints without clipping to include the mean.

    Sampling intervals can exclude the nominal mean, especially at tiny event
    rates. In that case draw the interval and mean independently, rather than
    modifying endpoints to manufacture nonnegative matplotlib yerr values.
    Returns a normal ErrorbarContainer for legends.
    """
    x, y, lower, upper = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float),
                                            np.asarray(lower, float), np.asarray(upper, float))
    if np.any(np.isfinite(lower) & np.isfinite(upper) & (upper < lower)):
        raise ValueError('Interval upper bound is below lower bound')
    if 'yerr' in kwargs:
        raise ValueError('Pass absolute interval endpoints, not yerr')
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(lower) & np.isfinite(upper)
    inside = (lower <= y) & (y <= upper)
    if np.all(inside[finite]):
        # NaN intervals are intentionally absent, not silently replaced by zero.
        return ax.errorbar(x, y, yerr=np.vstack([y - lower, upper - y]), **kwargs)
    artists = ax.errorbar(x, y, **kwargs)
    color = kwargs.get('ecolor', kwargs.get('color', None))
    if color is None and artists.lines[0] is not None:
        color = artists.lines[0].get_color()
    line_opts = {'color': color, 'alpha': kwargs.get('alpha', None),
                 'linewidth': kwargs.get('elinewidth', kwargs.get('linewidth', kwargs.get('lw', 1.0))),
                 'zorder': kwargs.get('zorder', 2)}
    ax.vlines(x[finite], lower[finite], upper[finite], **line_opts)
    capsize = kwargs.get('capsize', 0)
    if capsize:
        for endpoint in [lower, upper]:
            ax.plot(x[finite], endpoint[finite], linestyle='none', marker='_',
                    markersize=2*capsize, color=color, alpha=kwargs.get('alpha', None),
                    markeredgewidth=kwargs.get('capthick', line_opts['linewidth']),
                    zorder=line_opts['zorder'])
    return artists
