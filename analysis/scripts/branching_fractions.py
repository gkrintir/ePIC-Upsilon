"""Explicit dielectron branching conventions, with archived defaults preserved.

Only inclusive/visible conversions use these numbers. They do not rescale a
supplied visible generator cross section or the prepared event weights.
"""

import argparse
import math


LEGACY_BRANCHING_FRACTIONS = {"1S": 0.0234, "2S": 0.0191, "3S": 0.0218}


def positive_branching_fraction(value):
    """argparse-compatible finite probability, strictly between zero and one."""
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("branching fractions must be numbers") from exc
    if not math.isfinite(result) or not 0.0 < result < 1.0:
        raise argparse.ArgumentTypeError("branching fractions must be finite and strictly between 0 and 1")
    return result


def parse_branching_fractions(value):
    """Validate a state:value list; omitted states retain archived defaults."""
    result = dict(LEGACY_BRANCHING_FRACTIONS)
    if value is None:
        return result
    if isinstance(value, dict):
        pairs = list(value.items())
    else:
        if not isinstance(value, str) or not value.strip():
            raise argparse.ArgumentTypeError("expected a nonempty list such as 1S:0.0238,2S:0.0191,3S:0.0218")
        pairs = []
        for item in value.split(","):
            parts = item.split(":")
            if len(parts) != 2:
                raise argparse.ArgumentTypeError("expected state:value pairs separated by commas")
            pairs.append((parts[0].strip(), parts[1].strip()))
    seen = set()
    for state, fraction in pairs:
        if state not in LEGACY_BRANCHING_FRACTIONS:
            raise argparse.ArgumentTypeError("unknown Upsilon state {!r}; choose 1S, 2S or 3S".format(state))
        if state in seen:
            raise argparse.ArgumentTypeError("duplicate branching fraction for {}".format(state))
        seen.add(state)
        result[state] = positive_branching_fraction(fraction)
    return result


def add_branching_fraction_argument(parser):
    return parser.add_argument(
        "--branching-fractions", type=parse_branching_fractions, default=None,
        metavar="STATE:FRACTION,...",
        help="Dielectron fractions for visible/inclusive conversions. Omitted states retain archived defaults "
             "(1S:0.0234,2S:0.0191,3S:0.0218). Latest generator convention: "
             "1S:0.0238,2S:0.0191,3S:0.0218. No branching uncertainty is added.",
    )


def branching_fractions_from_args(args):
    return parse_branching_fractions(getattr(args, "branching_fractions", None))


def branching_fraction_metadata(args):
    return {
        "branching_fractions_ee": branching_fractions_from_args(args),
        "branching_fraction_convention": "legacy_default" if getattr(args, "branching_fractions", None) is None else "explicit_cli",
        "branching_fraction_uncertainty_included": False,
    }
