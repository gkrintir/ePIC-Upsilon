#!/usr/bin/env python3
"""Run the statistics pipeline after checking its cache contract.

Uses the existing workflow runner without replacing the user's runner file.
The preflight occurs before the runner can reuse outputs or rewrite manifests.
"""
import datetime as dt
from pathlib import Path
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
    sys.modules["tomllib"] = tomllib

import run_analysis as runner
from check_statistical_run import build_contract, check_contract


def main():
    parser = runner.build_parser()
    parser.description = __doc__
    parser.set_defaults(pipeline=str(Path(__file__).resolve().parents[1] / "config/pipeline_statistics.toml"))
    args = parser.parse_args()
    try:
        if args.command == "run":
            args.run_id = args.run_id or dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            config = runner.load_configuration(args, args.run_id)
            from proton_selection import guard_proton_run
            guard_proton_run(config, runner.PROJECT_ROOT, dry_run=args.dry_run)
            settings = runner.flatten({key: value for key, value in config.items()
                                       if key not in ("groups", "task", "run", "paths")})
            settings = [key + "=" + str(value) for key, value in sorted(settings.items())]
            contract = build_contract(settings, runner.PROJECT_ROOT, config.get("task", []))
            if not args.dry_run:
                status = check_contract(Path(config["run"]["dir"]), contract)
                print("Statistical contract " + status, flush=True)
        return int(args.func(args))
    except (runner.ConfigurationError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
