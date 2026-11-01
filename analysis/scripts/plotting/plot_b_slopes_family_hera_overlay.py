#!/usr/bin/env python3
"""
Overlay ePIC Υ(1S), Υ(2S), and Υ(3S) b-slope points on a ZEUS Fig. 4 style
comparison plot.

Input
-----
A b-slope summary CSV produced by extract_differential_gamma_p_swan_v4.py, e.g.

  diff_gamma_p_family_18x275_W37to114_bSlope_b_slope_summary.csv

This script drops the lower-energy overlay and instead shows the three Upsilon
states from the same high-energy sample. The x coordinate is

  Q^2 + M_V^2,

so the three states appear at different hard scales:
  Υ(1S): ~ 89.5 GeV^2
  Υ(2S): ~100.5 GeV^2
  Υ(3S): ~107.2 GeV^2

The y uncertainties are the statistical b-fit uncertainties stored in the summary
CSV, obtained from the finite-bin likelihood/chi2 scan in the extraction script.
"""

import argparse
import csv
import math
import os

import numpy as np
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass


UPSILON = {
    "1S": {"mass": 9.46030, "color": "royalblue", "marker": "D"},
    "2S": {"mass": 10.02326, "color": "crimson", "marker": "s"},
    "3S": {"mass": 10.35520, "color": "forestgreen", "marker": "^"},
}


def read_csv(path):
    with open(os.path.expanduser(path), newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def ffloat(row, key, default=np.nan):
    try:
        return float(row.get(key, default))
    except Exception:
        return default


def load_epic_points(summary_csv, q2):
    rows = read_csv(summary_csv)
    out = []
    for row in rows:
        state = row.get("state", "")
        if state not in UPSILON:
            continue
        mass = UPSILON[state]["mass"]
        out.append({
            "state": state,
            "x": mass * mass + q2,
            "b": ffloat(row, "b_slope_GeVminus2"),
            "err_low": ffloat(row, "b_slope_err_low", 0.0),
            "err_high": ffloat(row, "b_slope_err_high", 0.0),
            "chi2": ffloat(row, "chi2"),
            "ndof": ffloat(row, "ndof"),
            "valid_bins": ffloat(row, "valid_bins"),
            "mass_window_min": ffloat(row, "mass_window_min"),
            "mass_window_max": ffloat(row, "mass_window_max"),
        })
    if not out:
        raise RuntimeError(f"No Upsilon rows found in {summary_csv}")
    return out


def hera_reference_points(include_zeus_syst=False):
    """Compact approximate HERA Fig. 4 context points.

    The ZEUS Υ(1S) point is the actual reported value. Other points are visual
    context anchors for a Fig. 4 style comparison, not a precision redigitization.
    """
    pts = []
    raw = [
        ("rho ZEUS/H1", 0.6, 9.0, 3.0, 3.0, "o", "tab:red"),
        ("rho ZEUS/H1", 1.5, 7.2, 1.3, 1.3, "o", "tab:red"),
        ("rho ZEUS/H1", 3.0, 7.5, 0.8, 0.8, "o", "tab:red"),
        ("rho ZEUS/H1", 6.0, 6.2, 0.7, 0.7, "o", "tab:red"),
        ("rho ZEUS/H1", 10.0, 5.5, 0.6, 0.6, "o", "tab:red"),
        ("phi ZEUS/H1", 2.0, 7.4, 1.1, 1.1, "^", "magenta"),
        ("phi ZEUS/H1", 5.0, 6.6, 0.8, 0.8, "^", "magenta"),
        ("phi ZEUS/H1", 13.0, 5.5, 0.8, 0.8, "^", "magenta"),
        ("J/psi ZEUS/H1", 9.6, 4.6, 0.5, 0.5, "s", "tab:green"),
        ("J/psi ZEUS/H1", 11.0, 4.1, 0.4, 0.4, "s", "tab:green"),
        ("J/psi ZEUS/H1", 16.0, 4.7, 0.4, 0.4, "s", "tab:green"),
        ("J/psi ZEUS/H1", 25.0, 4.2, 0.5, 0.5, "s", "tab:green"),
        ("DVCS ZEUS/H1", 4.0, 5.8, 0.9, 0.9, "*", "tab:blue"),
        ("DVCS ZEUS/H1", 8.0, 5.4, 0.8, 0.8, "*", "tab:blue"),
        ("DVCS ZEUS/H1", 15.0, 5.2, 0.7, 0.7, "*", "tab:blue"),
        ("DVCS ZEUS/H1", 25.0, 5.0, 0.9, 0.9, "*", "tab:blue"),
        ("DVCS ZEUS/H1", 40.0, 4.9, 0.9, 0.9, "*", "tab:blue"),
    ]
    for name, x, b, elo, ehi, marker, color in raw:
        pts.append({
            "name": name, "x": x, "b": b, "err_low": elo, "err_high": ehi,
            "marker": marker, "color": color,
        })

    zeus_elo = math.sqrt(1.3**2 + 0.6**2) if include_zeus_syst else 1.3
    zeus_ehi = math.sqrt(2.0**2 + 0.5**2) if include_zeus_syst else 2.0
    pts.append({
        "name": "ZEUS Υ(1S)",
        "x": UPSILON["1S"]["mass"] ** 2,
        "b": 4.3,
        "err_low": zeus_elo,
        "err_high": zeus_ehi,
        "marker": "P",
        "color": "black",
    })
    return pts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, help="b_slope_summary.csv from extract_differential_gamma_p_swan_v4.py")
    p.add_argument("--q2", type=float, default=0.001)
    p.add_argument("--sample-label", default=r"ePIC $18\times275$ GeV")
    p.add_argument("--include-zeus-syst", action="store_true")
    p.add_argument("--out-prefix", default="b_slope_fig4_overlay_three_upsilon_states")
    args = p.parse_args()

    epic = load_epic_points(args.summary, args.q2)
    refs = hera_reference_points(include_zeus_syst=args.include_zeus_syst)

    fig, ax = plt.subplots(figsize=(8.6, 6.5))

    grouped = {}
    for r in refs:
        grouped.setdefault(r["name"], []).append(r)

    for name, rows in grouped.items():
        x = np.asarray([r["x"] for r in rows])
        y = np.asarray([r["b"] for r in rows])
        yerr = np.vstack([
            np.asarray([r["err_low"] for r in rows]),
            np.asarray([r["err_high"] for r in rows]),
        ])
        marker = rows[0]["marker"]
        color = rows[0]["color"]
        ax.errorbar(
            x, y, yerr=yerr,
            fmt=marker, color=color, markerfacecolor="none",
            linestyle="none", capsize=2, markersize=7,
            label=name, alpha=0.75,
        )

    for point in epic:
        st = UPSILON[point["state"]]
        if not np.isfinite(point["b"]):
            continue
        ax.errorbar(
            point["x"], point["b"],
            yerr=np.array([[point["err_low"]], [point["err_high"]]]),
            fmt=st["marker"], color=st["color"], markerfacecolor=st["color"],
            markersize=8, capsize=3,
            label=rf"{args.sample_label} $\Upsilon({point['state']})$",
            zorder=10,
        )
        ax.text(
            point["x"] + 1.0, point["b"],
            rf"$\Upsilon({point['state']})$",
            color=st["color"], fontsize=9, va="center"
        )

    ax.set_xlabel(r"$Q^2 + M_{V}^{2}$ [GeV$^2$]", fontsize=14)
    ax.set_ylabel(r"$b$ [GeV$^{-2}$]", fontsize=14)
    ax.set_xlim(0, 115)
    ax.set_ylim(0, 14)
    ax.set_title(r"Exponential $t$-slope comparison with HERA measurements", fontsize=15)
    ax.text(
        0.04, 0.95,
        r"ePIC points: finite-bin fit to $b\,e^{-b|t|}$",
        transform=ax.transAxes, ha="left", va="top", fontsize=10.5
    )
    ax.text(
        0.04, 0.90,
        r"$Q^2\simeq 10^{-3}$ GeV$^2$ for the ePIC photoproduction sample",
        transform=ax.transAxes, ha="left", va="top", fontsize=10.0
    )

    ax.legend(frameon=False, fontsize=8.0, ncol=2, loc="upper right")
    fig.tight_layout()

    out_png = args.out_prefix + ".png"
    out_pdf = args.out_prefix + ".pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)

    out_csv = args.out_prefix + "_epic_points.csv"
    with open(out_csv, "w", newline="") as f:
        fields = ["state", "x", "b", "err_low", "err_high", "chi2", "ndof", "valid_bins", "mass_window_min", "mass_window_max"]
        wr = csv.DictWriter(f, fieldnames=fields)
        wr.writeheader()
        for row in epic:
            wr.writerow({k: row.get(k, "") for k in fields})

    print("=== ePIC b-slope points ===")
    for point in epic:
        print(
            f"Υ({point['state']}): x={point['x']:.3f}, "
            f"b={point['b']:.3f} +{point['err_high']:.3f}/-{point['err_low']:.3f}, "
            f"chi2/ndof={point['chi2']:.3f}/{point['ndof']:.0f}"
        )
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
