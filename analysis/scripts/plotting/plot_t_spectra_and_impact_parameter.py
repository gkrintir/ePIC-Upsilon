#!/usr/bin/env python3
"""
Overlay the unfolded |t| spectrum obtained with the pT^2 estimator at two beam
energies together with the ZEUS exclusive Upsilon(1S) reference, and show the
corresponding extracted b_T profile in a neighboring panel.

Intended inputs
  - low and high CSVs from unfold_t_distribution_swan_v5.py
    (for method pT2), i.e. files like
      <prefix>_t_unfold_summary.csv
  - ZEUS step CSV from
      zeus_desy_11_186_fig4_exclusive_upsilon_step.csv

The left panel is a shape comparison in |t|.
The right panel shows the corresponding b_T profiles from a discrete
Fourier-Bessel transform of sqrt(dN/d|t|). Since the absolute normalization of
that profile is convention-dependent for this purpose, the b_T curves are
normalized to unit peak by default.
"""

import argparse
import csv
import math
import os
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise RuntimeError(f"Matplotlib import failed: {exc}")

try:
    import mplhep as hep
    hep.style.use(hep.style.ROOT)
except Exception:
    pass

try:
    from scipy.special import j0 as scipy_j0
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
    except Exception:
        ROOT = None

HBARC_GEV_FM = 0.1973269804


def read_csv_numeric(path):
    path = os.path.expanduser(path)
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    out = {}
    for key in rows[0].keys():
        vals = []
        numeric = True
        for r in rows:
            v = r[key]
            try:
                vals.append(float(v))
            except Exception:
                numeric = False
                break
        if numeric:
            out[key] = np.asarray(vals, dtype=float)
        else:
            out[key] = np.asarray([r[key] for r in rows], dtype=object)
    return out


def get_col(d, names, required=True, default=None):
    for n in names:
        if n in d:
            return d[n]
    if required:
        raise KeyError(f"Could not find any of the columns {names}. Available columns: {list(d.keys())}")
    return default


def load_epic_summary(path):
    d = read_csv_numeric(path)
    t_low = get_col(d, ["t_low"])
    t_high = get_col(d, ["t_high"])
    t_center = get_col(d, ["t_center"], required=False, default=0.5 * (t_low + t_high))
    y = get_col(d, ["dNdt_unfolded_visible", "dNdt_truth_test_visible", "dsigma_dt"])
    yerr = get_col(d, ["dNdt_unfolded_visible_err", "dNdt_truth_test_visible_err", "dsigma_dt_err"], required=False, default=np.zeros_like(y))
    return {
        "t_low": np.asarray(t_low, dtype=float),
        "t_high": np.asarray(t_high, dtype=float),
        "t_center": np.asarray(t_center, dtype=float),
        "xerr": 0.5 * (np.asarray(t_high, dtype=float) - np.asarray(t_low, dtype=float)),
        "y": np.asarray(y, dtype=float),
        "yerr": np.asarray(yerr, dtype=float),
    }


def load_zeus_step(path):
    d = read_csv_numeric(path)
    t_low = get_col(d, ["t_low"])
    t_high = get_col(d, ["t_high"])
    t_center = get_col(d, ["t_center"], required=False, default=0.5 * (t_low + t_high))
    y = get_col(d, ["dNdt_GeVminus2", "dNdt"])
    yerr_low = get_col(d, ["err_low"], required=False, default=np.zeros_like(y))
    yerr_high = get_col(d, ["err_high"], required=False, default=np.zeros_like(y))
    return {
        "t_low": np.asarray(t_low, dtype=float),
        "t_high": np.asarray(t_high, dtype=float),
        "t_center": np.asarray(t_center, dtype=float),
        "xerr": 0.5 * (np.asarray(t_high, dtype=float) - np.asarray(t_low, dtype=float)),
        "y": np.asarray(y, dtype=float),
        "yerr_low": np.asarray(yerr_low, dtype=float),
        "yerr_high": np.asarray(yerr_high, dtype=float),
    }


def bessel_j0(x):
    if HAS_SCIPY:
        return scipy_j0(x)
    if ROOT is None:
        # crude fallback
        x = np.asarray(x, dtype=float)
        out = np.ones_like(x)
        small = np.abs(x) < 1e-3
        xs = x[~small]
        out[small] = 1.0 - 0.25 * x[small] * x[small]
        # asymptotic fallback only for rough plotting
        out[~small] = np.sqrt(2.0 / (np.pi * np.maximum(xs, 1e-9))) * np.cos(xs - np.pi / 4.0)
        return out
    vec = np.vectorize(lambda z: ROOT.TMath.BesselJ0(float(z)))
    return vec(x)


def build_step_edges(t_low, t_high, y):
    edges = np.r_[t_low, t_high[-1]]
    vals = np.r_[y, y[-1]]
    return edges, vals


def rebin_value_at_centers(src, target_centers):
    # piecewise constant in source binning
    out = np.full_like(target_centers, np.nan, dtype=float)
    for i, xc in enumerate(target_centers):
        idx = np.where((xc >= src["t_low"]) & (xc < src["t_high"]))[0]
        if len(idx) == 0 and np.isclose(xc, src["t_high"][-1]):
            idx = [len(src["t_low"]) - 1]
        if len(idx) > 0:
            out[i] = src["y"][idx[0]]
    return out


def fit_scale_to_reference(sample, ref):
    ref_at_sample = rebin_value_at_centers(ref, sample["t_center"])
    m = np.isfinite(sample["y"]) & (sample["y"] > 0) & np.isfinite(ref_at_sample) & (ref_at_sample > 0)
    if not np.any(m):
        return 1.0
    # unweighted shape scale is enough here
    num = np.sum(ref_at_sample[m] * sample["y"][m])
    den = np.sum(sample["y"][m] ** 2)
    return num / den if den > 0 else 1.0


def normalize_area(sample):
    area = np.sum(sample["y"] * (sample["t_high"] - sample["t_low"]))
    if area <= 0:
        return 1.0
    return 1.0 / area


def transform_to_bT(sample, b_grid_fm, nsub=80):
    t_low = sample["t_low"]
    t_high = sample["t_high"]
    y = np.maximum(sample["y"], 0.0)

    # F(b) proportional to ∫ dΔ Δ J0(bΔ) sqrt(dN/dt), with Δ=sqrt(t)
    prof = np.zeros_like(b_grid_fm, dtype=float)
    for i in range(len(t_low)):
        if y[i] <= 0 or t_high[i] <= t_low[i]:
            continue
        delta_lo = math.sqrt(max(t_low[i], 0.0))
        delta_hi = math.sqrt(max(t_high[i], 0.0))
        delta = np.linspace(delta_lo, delta_hi, nsub)
        amp = math.sqrt(y[i])
        for ib, bfm in enumerate(b_grid_fm):
            bg = bfm / HBARC_GEV_FM
            integrand = delta * bessel_j0(bg * delta) * amp
            prof[ib] += np.trapezoid(integrand, delta)
    # Use modulus as a shape-like profile
    prof = np.abs(prof)
    return prof


def toy_bT_unc(sample, b_grid_fm, ntoys=200, seed=12345):
    rng = np.random.default_rng(seed)
    nom = transform_to_bT(sample, b_grid_fm)
    if "yerr" not in sample or sample["yerr"] is None:
        return nom, np.zeros_like(nom)
    if np.all(sample["yerr"] <= 0):
        return nom, np.zeros_like(nom)
    toys = []
    for _ in range(ntoys):
        ytoy = rng.normal(sample["y"], sample["yerr"])
        stoy = dict(sample)
        stoy["y"] = np.maximum(ytoy, 0.0)
        toys.append(transform_to_bT(stoy, b_grid_fm))
    toys = np.asarray(toys)
    return nom, np.std(toys, axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--low-csv", required=True, help="10x100 pT2 t-unfold summary CSV")
    ap.add_argument("--high-csv", required=True, help="18x275 pT2 t-unfold summary CSV")
    ap.add_argument("--zeus-csv", required=True, help="ZEUS step CSV")
    ap.add_argument("--low-label", default=r"ePIC 10$\times$100, $p_{T,\Upsilon}^{2}$")
    ap.add_argument("--high-label", default=r"ePIC 18$\times$275, $p_{T,\Upsilon}^{2}$")
    ap.add_argument("--zeus-label", default=r"ZEUS exclusive $\Upsilon(1S)$ component")
    ap.add_argument("--low-note", default=r"$W_{\gamma p}\in[20,26.5]$ GeV, $Q^2<10^{-3}$ GeV$^2$")
    ap.add_argument("--high-note", default=r"$W_{\gamma p}\in[37,114]$ GeV, $Q^2<10^{-3}$ GeV$^2$")
    ap.add_argument("--shape-mode", choices=["scale_to_zeus", "unit_area", "none"], default="scale_to_zeus")
    ap.add_argument("--bt-norm", choices=["peak", "area", "none"], default="peak")
    ap.add_argument("--bt-max-fm", type=float, default=2.5)
    ap.add_argument("--n-bt", type=int, default=160)
    ap.add_argument("--n-toys", type=int, default=200)
    ap.add_argument("--out-prefix", default="overlay_pt2_twoE_ZEUS_bT")
    args = ap.parse_args()

    low = load_epic_summary(args.low_csv)
    high = load_epic_summary(args.high_csv)
    zeus = load_zeus_step(args.zeus_csv)

    # Scale only for the left-panel shape comparison.
    low_scale = 1.0
    high_scale = 1.0
    zeus_scale = 1.0
    if args.shape_mode == "scale_to_zeus":
        low_scale = fit_scale_to_reference(low, zeus)
        high_scale = fit_scale_to_reference(high, zeus)
    elif args.shape_mode == "unit_area":
        low_scale = normalize_area(low)
        high_scale = normalize_area(high)
        zeus_area = np.sum(zeus["y"] * (zeus["t_high"] - zeus["t_low"]))
        zeus_scale = 1.0 / zeus_area if zeus_area > 0 else 1.0

    low_plot = dict(low)
    high_plot = dict(high)
    zeus_plot = dict(zeus)
    low_plot["y"] = low["y"] * low_scale
    low_plot["yerr"] = low["yerr"] * low_scale
    high_plot["y"] = high["y"] * high_scale
    high_plot["yerr"] = high["yerr"] * high_scale
    zeus_plot["y"] = zeus["y"] * zeus_scale
    zeus_plot["yerr_low"] = zeus["yerr_low"] * zeus_scale
    zeus_plot["yerr_high"] = zeus["yerr_high"] * zeus_scale

    bgrid = np.linspace(0.0, args.bt_max_fm, args.n_bt)
    low_bt, low_bt_err = toy_bT_unc(low_plot, bgrid, ntoys=args.n_toys, seed=11)
    high_bt, high_bt_err = toy_bT_unc(high_plot, bgrid, ntoys=args.n_toys, seed=17)
    zeus_for_bt = {
        "t_low": zeus_plot["t_low"], "t_high": zeus_plot["t_high"],
        "t_center": zeus_plot["t_center"], "y": zeus_plot["y"],
        "yerr": 0.5 * (zeus_plot["yerr_low"] + zeus_plot["yerr_high"]),
    }
    zeus_bt, zeus_bt_err = toy_bT_unc(zeus_for_bt, bgrid, ntoys=max(50, args.n_toys // 2), seed=23)

    def norm_bt(y, yerr):
        y = y.copy()
        yerr = yerr.copy()
        if args.bt_norm == "peak":
            s = np.max(y) if np.max(y) > 0 else 1.0
        elif args.bt_norm == "area":
            s = np.trapezoid(y, bgrid) if np.trapezoid(y, bgrid) > 0 else 1.0
        else:
            s = 1.0
        return y / s, yerr / s, s

    low_btn, low_btn_err, low_bts = norm_bt(low_bt, low_bt_err)
    high_btn, high_btn_err, high_bts = norm_bt(high_bt, high_bt_err)
    zeus_btn, zeus_btn_err, zeus_bts = norm_bt(zeus_bt, zeus_bt_err)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.2))

    # Left panel, |t| comparison.
    zedges, zvals = build_step_edges(zeus_plot["t_low"], zeus_plot["t_high"], zeus_plot["y"])
    ax1.step(zedges, zvals, where="post", color="black", linestyle=":", linewidth=2.8, label=args.zeus_label, zorder=3)
    if np.any(zeus_plot["yerr_low"] > 0) or np.any(zeus_plot["yerr_high"] > 0):
        ax1.errorbar(zeus_plot["t_center"], zeus_plot["y"], xerr=zeus_plot["xerr"],
                     yerr=np.vstack([zeus_plot["yerr_low"], zeus_plot["yerr_high"]]),
                     fmt="none", ecolor="black", capsize=2, zorder=3)

    ax1.errorbar(low_plot["t_center"], low_plot["y"], xerr=low_plot["xerr"], yerr=low_plot["yerr"],
                 fmt="o-", color="tab:blue", capsize=3, label=args.low_label, zorder=5)
    ax1.errorbar(high_plot["t_center"], high_plot["y"], xerr=high_plot["xerr"], yerr=high_plot["yerr"],
                 fmt="s-", color="tab:orange", capsize=3, label=args.high_label, zorder=4)

    ax1.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax1.set_ylabel(r"scaled $dN/d|t|$ [arb.]" if args.shape_mode != "none" else r"$dN/d|t|$ [GeV$^{-2}$]")
    ax1.set_yscale("log")
    ax1.set_xlim(0.0, max(low_plot["t_high"].max(), high_plot["t_high"].max(), zeus_plot["t_high"].max()))
    positive = np.concatenate([low_plot["y"][low_plot["y"] > 0], high_plot["y"][high_plot["y"] > 0], zeus_plot["y"][zeus_plot["y"] > 0]])
    ymin = max(np.min(positive) * 0.5, 1e-2) if len(positive) else 1e-2
    ymax = np.max(positive) * 3.0 if len(positive) else 1.0
    ax1.set_ylim(ymin, ymax)
    ax1.set_title(r"$p_{T,\Upsilon}^{2}$-based unfolded $|t|$ comparison")
    ax1.text(0.03, 0.97, "ePIC and ZEUS", transform=ax1.transAxes, va="top", ha="left", fontsize=12, fontweight="bold")
    ax1.text(0.03, 0.90, args.low_note, transform=ax1.transAxes, va="top", ha="left", fontsize=9.5, color="tab:blue")
    ax1.text(0.03, 0.84, args.high_note, transform=ax1.transAxes, va="top", ha="left", fontsize=9.5, color="tab:orange")
    ax1.legend(frameon=False, fontsize=9)

    # Right panel, bT profiles.
    ax2.plot(bgrid, zeus_btn, color="black", linestyle=":", linewidth=2.8, label=args.zeus_label)
    if np.any(zeus_btn_err > 0):
        ax2.fill_between(bgrid, np.maximum(zeus_btn - zeus_btn_err, 0.0), zeus_btn + zeus_btn_err,
                         color="black", alpha=0.10)
    ax2.plot(bgrid, low_btn, color="tab:blue", linewidth=2.0, label=args.low_label)
    if np.any(low_btn_err > 0):
        ax2.fill_between(bgrid, np.maximum(low_btn - low_btn_err, 0.0), low_btn + low_btn_err,
                         color="tab:blue", alpha=0.18)
    ax2.plot(bgrid, high_btn, color="tab:orange", linewidth=2.0, label=args.high_label)
    if np.any(high_btn_err > 0):
        ax2.fill_between(bgrid, np.maximum(high_btn - high_btn_err, 0.0), high_btn + high_btn_err,
                         color="tab:orange", alpha=0.18)

    ax2.set_xlabel(r"$b_T$ [fm]")
    if args.bt_norm == "peak":
        ax2.set_ylabel(r"normalized $b_T$ profile")
    elif args.bt_norm == "area":
        ax2.set_ylabel(r"area-normalized $b_T$ profile")
    else:
        ax2.set_ylabel(r"$b_T$ profile [arb.]")
    ax2.set_xlim(0.0, args.bt_max_fm)
    ax2.set_ylim(bottom=0.0)
    ax2.set_title(r"Extracted $b_T$ profile from the same $|t|$ shapes")
    ax2.text(0.03, 0.97, r"Fourier-Bessel transform of $\,\sqrt{dN/d|t|}$", transform=ax2.transAxes,
             va="top", ha="left", fontsize=10.5)
    ax2.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    png = args.out_prefix + ".png"
    pdf = args.out_prefix + ".pdf"
    fig.savefig(png, dpi=250)
    fig.savefig(pdf)

    # Save bT curves.
    out_csv = args.out_prefix + "_bT_profiles.csv"
    with open(out_csv, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["bT_fm", "low_profile", "low_profile_err", "high_profile", "high_profile_err", "zeus_profile", "zeus_profile_err"])
        for i in range(len(bgrid)):
            wr.writerow([bgrid[i], low_btn[i], low_btn_err[i], high_btn[i], high_btn_err[i], zeus_btn[i], zeus_btn_err[i]])

    print("=== Scaling used for the left-panel |t| overlay ===")
    print(f"shape_mode = {args.shape_mode}")
    print(f"low_scale   = {low_scale:.6g}")
    print(f"high_scale  = {high_scale:.6g}")
    print(f"zeus_scale  = {zeus_scale:.6g}")
    print(f"bt_norm     = {args.bt_norm}")
    print(f"low_bt_norm = {low_bts:.6g}")
    print(f"high_bt_norm= {high_bts:.6g}")
    print(f"zeus_bt_norm= {zeus_bts:.6g}")
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
