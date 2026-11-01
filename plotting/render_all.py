#!/usr/bin/env python3
"""Render the numerical framework outputs selected by one explicit plot profile."""
import argparse
import json
from pathlib import Path
import sys

from settings import load_profile, make_config, make_reoptimization_config, new_output_dir, apply_display_rebin
from upsilon_esr import preflight
from upsilon_esr.common import write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, help='Defaults to UPSILON_PLOT_PROFILE or configs/plots.local.json')
    parser.add_argument('--section', choices=('all', 'reconstruction', 'w', 'differential'), default='all')
    parser.add_argument('--out', type=Path, help='Fresh output directory; an existing directory is refused')
    parser.add_argument('--allow-diagnostics', action='store_true', help='Explicitly permit NOT_VALIDATED saved diagnostics; preserve their status')
    parser.add_argument('--check', action='store_true', help='Validate numerical inputs without rendering figures')
    args = parser.parse_args(argv)
    profile = load_profile(args.profile)
    if args.allow_diagnostics:
        profile['allow_diagnostic_results'] = True
    if args.out is not None:
        out = args.out.expanduser().resolve()
        if out.exists():
            raise FileExistsError('Output directory already exists; choose a fresh --out: ' + str(out))
    elif args.check:
        out = profile['output_root'] / 'read-only-preflight'
    else:
        out = new_output_dir(profile, args.section)
    sections = ('reconstruction', 'w', 'differential') if args.section == 'all' else (args.section,)
    jobs = []
    for section in sections:
        maker = make_config if section == 'reconstruction' else make_reoptimization_config
        cfg = maker(profile=profile, output_dir=out)
        if section != 'reconstruction':
            apply_display_rebin(cfg, profile, section)
        print(section)
        print(preflight(cfg, section).to_string(index=False))
        jobs.append((section, cfg))
    if args.check:
        print('Input validation completed. No figures or results were written.')
        return 0
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'plot_profile.json', profile)
    for section, cfg in jobs:
        if section == 'reconstruction':
            from upsilon_esr.reconstruction import plot_reconstruction
            plot_reconstruction(cfg, rebin=2, float_common_shape=True,
                                efficiency_variables=('true_W', 'true_pT'), show_mass_intervals=False)
        elif section == 'w':
            from upsilon_esr.w_cross_section import plot_w_cross_section
            plot_w_cross_section(cfg)
        else:
            from upsilon_esr.differential import plot_differential
            plot_differential(cfg, fit_range=(0., 1.))
    print('Figures and numerical audits:', out)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, FileNotFoundError, FileExistsError, KeyError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        sys.exit(2)
