#!/usr/bin/env python3
"""
Overlay extracted ePIC Upsilon(1S) b slopes on a ZEUS Fig. 4 style comparison plot.

Inputs
------
CSV files produced by extract_differential_gamma_p_swan_v4.py, for example:
  diff_gamma_p_1S_10x100_W20to26p5_bSlope_b_slope_summary.csv
  diff_gamma_p_1S_18x275_W37to114_bSlope_b_slope_summary.csv

The x coordinate is the hard scale Q^2 + M_V^2. For photoproduction Q^2 is tiny, so
the default for Upsilon(1S) is approximately M_Upsilon^2 = 9.4603^2 GeV^2.

The script includes a compact, approximate reproduction of the HERA Fig. 4 reference
points to provide context. The ZEUS Upsilon(1S) reference value is set explicitly to
b = 4.3^{+2.0}_{-1.3} GeV^{-2}, with an optional systematic contribution ignored
by default unless --include-zeus-syst is passed.

This is intended as a presentation/analysis-note comparison, not a re-digitization
of the full HERA compilation.
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

M_UPSILON_1S = 9.46030


def read_csv(path):
    path = os.path.expanduser(path)
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def first_float(row, keys, default=np.nan):
    for k in keys:
        if k in row and str(row[k]).strip() != "":
            try:
                return float(row[k])
            except Exception:
                pass
    return default


def load_b_summary(path, label, scale_x, x_shift=0.0):
    rows = read_csv(path)
    # Prefer 1S row if present.
    row = None
    for r in rows:
        if r.get("state", "") == "1S":
            row = r
            break
    if row is None:
        row = rows[0]

    b = first_float(row, ["b_slope_GeVminus2", "b", "B", "B_GeVminus2"])
    elo = first_float(row, ["b_slope_err_low", "b_err_low", "B_err_low"], default=np.nan)
    ehi = first_float(row, ["b_slope_err_high", "b_err_high", "B_err_high"], default=np.nan)

    if not np.isfinite(elo):
        elo = first_float(row, ["b_slope_err", "b_err", "B_err"], default=0.0)
    if not np.isfinite(ehi):
        ehi = first_float(row, ["b_slope_err", "b_err", "B_err"], default=0.0)

    return {
        "label": label,
        "x": scale_x + x_shift,
        "b": b,
        "err_low": elo,
        "err_high": ehi,
        "path": path,
    }


def hera_reference_points(include_zeus_syst=False):
    """Approximate reference points for a Fig. 4 style HERA comparison.

    The values are intentionally compact and approximate except for the ZEUS Upsilon
    point, which is taken from the ZEUS Upsilon t-dependence paper.
    """
    pts = []

    # Light vector mesons and DVCS, approximate visual anchors from ZEUS Fig. 4.
    # name, x=Q2+M^2, b, err_low, err_high, marker
    raw = [
        ("rho ZEUS/H1", 0.6, 9.0, 3.0, 3.0, "o"),
        ("rho ZEUS/H1", 1.5, 7.2, 1.3, 1.3, "o"),
        ("rho ZEUS/H1", 3.0, 7.5, 0.8, 0.8, "o"),
        ("rho ZEUS/H1", 6.0, 6.2, 0.7, 0.7, "o"),
        ("rho ZEUS/H1", 10.0, 5.5, 0.6, 0.6, "o"),
        ("phi ZEUS/H1", 2.0, 7.4, 1.1, 1.1, "^"),
        ("phi ZEUS/H1", 5.0, 6.6, 0.8, 0.8, "^"),
        ("phi ZEUS/H1", 13.0, 5.5, 0.8, 0.8, "^"),
        ("J/psi ZEUS/H1", 9.6, 4.6, 0.5, 0.5, "s"),
        ("J/psi ZEUS/H1", 11.0, 4.1, 0.4, 0.4, "s"),
        ("J/psi ZEUS/H1", 16.0, 4.7, 0.4, 0.4, "s"),
        ("J/psi ZEUS/H1", 25.0, 4.2, 0.5, 0.5, "s"),
        ("DVCS ZEUS/H1", 4.0, 5.8, 0.9, 0.9, "*"),
        ("DVCS ZEUS/H1", 8.0, 5.4, 0.8, 0.8, "*"),
        ("DVCS ZEUS/H1", 15.0, 5.2, 0.7, 0.7, "*"),
        ("DVCS ZEUS/H1", 25.0, 5.0, 0.9, 0.9, "*"),
        ("DVCS ZEUS/H1", 40.0, 4.9, 0.9, 0.9, "*"),
    ]

    for name, x, b, elo, ehi, marker in raw:
        pts.append({"name": name, "x": x, "b": b, "err_low": elo, "err_high": ehi, "marker": marker})

    # ZEUS Upsilon t-slope result.
    x_ups = M_UPSILON_1S ** 2
    stat_low = 1.3
    stat_high = 2.0
    if include_zeus_syst:
        # Combine stat and syst in quadrature: syst -0.6/+0.5
        stat_low = math.sqrt(1.3**2 + 0.6**2)
        stat_high = math.sqrt(2.0**2 + 0.5**2)
    pts.append({"name": "ZEUS Upsilon(1S)", "x": x_ups, "b": 4.3, "err_low": stat_low, "err_high": stat_high, "marker": "P"})

    return pts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--low-summary", required=True, help="10x100 b_slope_summary.csv")
    p.add_argument("--high-summary", required=True, help="18x275 b_slope_summary.csv")
    p.add_argument("--low-label", default=r"ePIC $10\times100$ GeV")
    p.add_argument("--high-label", default=r"ePIC $18\times275$ GeV")
    p.add_argument("--state", default="1S")
    p.add_argument("--q2-low", type=float, default=0.001)
    p.add_argument("--q2-high", type=float, default=0.001)
    p.add_argument("--mass", type=float, default=M_UPSILON_1S)
    p.add_argument("--include-zeus-syst", action="store_true")
    p.add_argument("--out-prefix", default="b_slope_fig4_overlay")
    args = p.parse_args()

    x_low = args.mass**2 + args.q2_low
    x_high = args.mass**2 + args.q2_high

    # Slight visual x-offset so the two ePIC points do not fully overlap.
    low = load_b_summary(args.low_summary, args.low_label, x_low, x_shift=-1.2)
    high = load_b_summary(args.high_summary, args.high_label, x_high, x_shift=+1.2)

    refs = hera_reference_points(include_zeus_syst=args.include_zeus_syst)

    fig, ax = plt.subplots(figsize=(8.4, 6.4))

    # Group approximate HERA references.
    groups = {}
    for r in refs:
        groups.setdefault(r["name"], []).append(r)

    styles = {
        "rho ZEUS/H1": dict(marker="o", color="tab:red", mfc="none"),
        "phi ZEUS/H1": dict(marker="^", color="magenta", mfc="none"),
        "J/psi ZEUS/H1": dict(marker="s", color="tab:green", mfc="none"),
        "DVCS ZEUS/H1": dict(marker="*", color="tab:blue", mfc="none"),
        "ZEUS Upsilon(1S)": dict(marker="P", color="black", mfc="white"),
    }

    for name, pts in groups.items():
        x = np.asarray([r["x"] for r in pts])
        y = np.asarray([r["b"] for r in pts])
        yerr = np.vstack([
            np.asarray([r["err_low"] for r in pts]),
            np.asarray([r["err_high"] for r in pts]),
        ])
        st = styles.get(name, dict(marker="o", color="gray", mfc="none"))
        ax.errorbar(
            x, y, yerr=yerr,
            fmt=st["marker"], color=st["color"], markerfacecolor=st["mfc"],
            linestyle="none", capsize=2, markersize=7, label=name, alpha=0.85
        )

    # ePIC points.
    for point, color, marker in [(low, "tab:blue", "D"), (high, "tab:orange", "D")]:
        ax.errorbar(
            point["x"], point["b"],
            yerr=np.array([[point["err_low"]], [point["err_high"]]]),
            fmt=marker, color=color, markerfacecolor=color,
            markersize=8, capsize=3, label=point["label"], zorder=10
        )

    ax.set_xlabel(r"$Q^2 + M_{V}^{2}$ [GeV$^2$]", fontsize=14)
    ax.set_ylabel(r"$b$ [GeV$^{-2}$]", fontsize=14)
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 14)
    ax.set_title(r"Exponential $t$-slope comparison with HERA measurements", fontsize=15)
    ax.text(
        0.04, 0.95,
        r"$\Upsilon(1S)$ ePIC points: finite-bin fit to $b\,e^{-b|t|}$",
        transform=ax.transAxes, ha="left", va="top", fontsize=10.5
    )
    ax.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper right")
    fig.tight_layout()

    out_png = args.out_prefix + ".png"
    out_pdf = args.out_prefix + ".pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)

    out_csv = args.out_prefix + "_epic_points.csv"
    with open(out_csv, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["label", "x", "b", "err_low", "err_high", "path"])
        wr.writeheader()
        wr.writerow(low)
        wr.writerow(high)

    print("=== ePIC b-slope points ===")
    for pnt in [low, high]:
        print(f"{pnt['label']}: x={pnt['x']:.3f}, b={pnt['b']:.3f} +{pnt['err_high']:.3f}/-{pnt['err_low']:.3f}")
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
