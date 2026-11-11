#!/usr/bin/env python3
"""Run the grouped retry command with the historical 1S 9x275 default batch."""
from rematch_condor import default_batch
from retry_rematch import cli


if __name__ == "__main__":
    cli(default_batch())
