#!/usr/bin/env python3
"""
Compare eSTARlight internal photon-luminosity behaviour with an analytic WW flux
by scanning narrow W_gamma_p bins.

This does NOT need the wylumi file format. It reruns eSTARlight in the same W bins
for two beam-energy settings, parses the printed total cross section, and compares
the low/high cross-section ratio with the analytic photon-flux ratio.

In a narrow W bin, if the same gamma-p model is used in both runs,

    sigma_ep_low(bin) / sigma_ep_high(bin)  ~=  Phi_low(bin) / Phi_high(bin)

because the same sigma(gamma p -> V p) approximately cancels.

Usage:
  python3 scan_estarlight_flux_ratio.py \
    --executable ./e_starlight \
    --template slight.in \
    --low-electron-gamma 19569 --low-target-gamma 106.6 \
    --high-electron-gamma 35225 --high-target-gamma 293.2 \
    --w-edges "10,14.4167,18.8333,23.25,27.6667,32.0833,36.5,40.9167,45.3333,49.75,54.1667,58.5833,63" \
    --q2-max 1.0 \
    --out flux_ratio_scan_10x100_vs_18x275.csv

Notes:
  - The script writes one temporary working directory per bin/energy.
  - It sets W_GP_MIN/W_GP_MAX explicitly and N_EVENTS=1 for speed.
  - It parses "Total cross section: ... picob." from eSTARlight stdout.
  - The absolute cross section includes the gamma-p model, so use the ratio,
    not the absolute number, as the flux comparison.
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813


def parse_edges(s):
    return [float(x) for x in s.split(",") if x.strip()]


def photon_flux_dW(W, e_e, e_p, q2_max):
    s_minus_mp2 = 4.0 * e_e * e_p
    y = (W * W - M_P * M_P) / s_minus_mp2
    if y <= 0.0 or y >= 1.0:
        return 0.0

    q2_min = (M_E * M_E * y * y) / (1.0 - y)
    if q2_min <= 0.0 or q2_min >= q2_max:
        return 0.0

    log_term = math.log(q2_max / q2_min)
    dndy = (ALPHA_EM / (2.0 * math.pi)) * (
        ((1.0 + (1.0 - y) ** 2) / y) * log_term
        - (2.0 * (1.0 - y) / y) * (1.0 - q2_min / q2_max)
    )
    return dndy * (2.0 * W / s_minus_mp2)


def integrated_flux(w_low, w_high, e_e, e_p, q2_max, n_steps=1000):
    grid = np.linspace(w_low, w_high, n_steps)
    vals = np.array([photon_flux_dW(w, e_e, e_p, q2_max) for w in grid])
    return float(np.trapezoid(vals, grid))


def gamma_to_energy(gamma, mass):
    return gamma * mass


def update_param(text, key, value):
    # Replace "KEY = ..." while preserving comments.
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^#\n]*)(.*)$", re.MULTILINE)
    repl = rf"\g<1>{value}\g<3>"
    if pattern.search(text):
        return pattern.sub(repl, text)
    return text + f"\n{key} = {value}\n"


def make_card(template_text, electron_gamma, target_gamma, w_low, w_high, q2_min, q2_max, n_events):
    text = template_text
    replacements = {
        "ELECTRON_BEAM_GAMMA": f"{electron_gamma:.12g}",
        "TARGET_BEAM_GAMMA": f"{target_gamma:.12g}",
        "W_GP_MIN": f"{w_low:.12g}",
        "W_GP_MAX": f"{w_high:.12g}",
        "MIN_GAMMA_Q2": f"{q2_min:.12g}",
        "MAX_GAMMA_Q2": f"{q2_max:.12g}",
        "N_EVENTS": str(int(n_events)),
        "baseFileName": "slight",
    }
    for key, val in replacements.items():
        text = update_param(text, key, val)
    return text


def parse_xsec(stdout):
    """Parse eSTARlight total and generated cross sections.

    Returns both values in pb. Handles picob, nanob, femtob, microb.
    """
    unit_to_pb = {
        "femtob": 1e-3,
        "picob": 1.0,
        "nanob": 1e3,
        "microb": 1e6,
        "millib": 1e9,
        "barn": 1e12,
    }

    total_pb = None
    generated_pb = None

    # Examples:
    # " Total cross section: 13.764 picob."
    # " Total cross section: 0.000 femtob."
    # " The cross section of the generated sample is 327.584 femtob."
    m = re.search(r"Total\s+cross\s+section:\s*([0-9eE+\-.]+)\s*([A-Za-z]+)", stdout)
    if m:
        val = float(m.group(1))
        unit = m.group(2).lower()
        if unit in unit_to_pb:
            total_pb = val * unit_to_pb[unit]

    m = re.search(r"cross\s+section\s+of\s+the\s+generated\s+sample\s+is\s*([0-9eE+\-.]+)\s*([A-Za-z]+)", stdout)
    if m:
        val = float(m.group(1))
        unit = m.group(2).lower()
        if unit in unit_to_pb:
            generated_pb = val * unit_to_pb[unit]

    return total_pb, generated_pb

def run_estarlight(executable, template_text, electron_gamma, target_gamma, w_low, w_high, q2_min, q2_max, n_events, keep_dir=None):
    workdir = Path(tempfile.mkdtemp(prefix="esl_flux_scan_")) if keep_dir is None else Path(keep_dir)
    workdir.mkdir(parents=True, exist_ok=True)

    card = make_card(template_text, electron_gamma, target_gamma, w_low, w_high, q2_min, q2_max, n_events)
    (workdir / "slight.in").write_text(card)

    exe = Path(executable).resolve()
    proc = subprocess.run(
        [str(exe)],
        cwd=str(workdir),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    stdout = proc.stdout
    (workdir / "estarlight_stdout.txt").write_text(stdout)

    total_pb, generated_fb = parse_xsec(stdout)
    if total_pb is None:
        raise RuntimeError(f"Could not parse total cross section for W={w_low}-{w_high}. See {workdir}/estarlight_stdout.txt")
    if total_pb == 0.0:
        print(f"  NOTE: eSTARlight returned zero total cross section for W={w_low:g}-{w_high:g}; keeping bin as zero.")

    return total_pb, generated_fb, str(workdir), proc.returncode


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--executable", required=True, help="Path to e_starlight executable")
    p.add_argument("--template", required=True, help="Template slight.in")
    p.add_argument("--low-electron-gamma", type=float, required=True)
    p.add_argument("--low-target-gamma", type=float, required=True)
    p.add_argument("--high-electron-gamma", type=float, required=True)
    p.add_argument("--high-target-gamma", type=float, required=True)
    p.add_argument("--w-edges", required=True)
    p.add_argument("--q2-min", type=float, default=0.0)
    p.add_argument("--q2-max", type=float, default=1.0)
    p.add_argument("--n-events", type=int, default=1)
    p.add_argument("--out", default="estarlight_flux_ratio_scan.csv")
    p.add_argument("--keep-workdirs", action="store_true")
    args = p.parse_args()

    template_text = Path(args.template).read_text()
    edges = parse_edges(args.w_edges)

    low_ee = gamma_to_energy(args.low_electron_gamma, M_E)
    low_ep = gamma_to_energy(args.low_target_gamma, M_P)
    high_ee = gamma_to_energy(args.high_electron_gamma, M_E)
    high_ep = gamma_to_energy(args.high_target_gamma, M_P)

    print("=== Beam energies inferred from gammas ===")
    print(f"low:  Ee={low_ee:.6g} GeV, Ep={low_ep:.6g} GeV")
    print(f"high: Ee={high_ee:.6g} GeV, Ep={high_ep:.6g} GeV")
    print("")

    rows = []
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        print(f"[{i+1}/{len(edges)-1}] W bin {lo:g} - {hi:g} GeV")

        low_total_pb, low_gen_pb, low_dir, low_code = run_estarlight(
            args.executable, template_text,
            args.low_electron_gamma, args.low_target_gamma,
            lo, hi, args.q2_min, args.q2_max, args.n_events,
            keep_dir=None,
        )

        high_total_pb, high_gen_pb, high_dir, high_code = run_estarlight(
            args.executable, template_text,
            args.high_electron_gamma, args.high_target_gamma,
            lo, hi, args.q2_min, args.q2_max, args.n_events,
            keep_dir=None,
        )

        phi_low = integrated_flux(lo, hi, low_ee, low_ep, args.q2_max)
        phi_high = integrated_flux(lo, hi, high_ee, high_ep, args.q2_max)

        internal_ratio = low_total_pb / high_total_pb if high_total_pb > 0 else np.nan
        analytic_ratio = phi_low / phi_high if phi_high > 0 else np.nan

        rows.append({
            "W_low": lo,
            "W_high": hi,
            "W_center": 0.5 * (lo + hi),
            "low_total_pb": low_total_pb,
            "high_total_pb": high_total_pb,
            "low_generated_pb": low_gen_pb,
            "high_generated_pb": high_gen_pb,
            "internal_ratio_low_over_high": internal_ratio,
            "analytic_flux_low": phi_low,
            "analytic_flux_high": phi_high,
            "analytic_flux_ratio_low_over_high": analytic_ratio,
            "internal_over_analytic_ratio": internal_ratio / analytic_ratio if analytic_ratio and np.isfinite(analytic_ratio) else np.nan,
            "low_workdir": low_dir,
            "high_workdir": high_dir,
            "low_returncode": low_code,
            "high_returncode": high_code,
        })

        if not args.keep_workdirs:
            shutil.rmtree(low_dir, ignore_errors=True)
            shutil.rmtree(high_dir, ignore_errors=True)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)

    print("\n=== Summary ===")
    cols = [
        "W_low", "W_high",
        "internal_ratio_low_over_high",
        "analytic_flux_ratio_low_over_high",
        "internal_over_analytic_ratio",
        "low_total_pb", "high_total_pb",
    ]
    print(df[cols].to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
