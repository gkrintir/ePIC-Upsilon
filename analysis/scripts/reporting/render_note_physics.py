#!/usr/bin/env python3
"""Render supported single-state physics figures in the supplied note layout.

Consumes existing analysis CSVs. No event selections, extraction, or unfolding
are repeated. The transverse profile follows the existing framework's modulus
of the Fourier--Bessel transform of sqrt(dN/dt), with piecewise-constant bins,
finite t support, and peak normalization. It is a shape diagnostic only.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import j0

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projected_intervals import (DEFAULT_COVERAGE, add_interval_arguments,
                                 draw_absolute_intervals, sampling_interval)
from note_style import BLUE, ORANGE, finish, set_style
import matplotlib.pyplot as plt

HBARC_GEV_FM = 0.1973269804
UNCERTAINTY_NOTES = {
    "mc": "Signal-only MC; error bars show MC statistical precision",
    "projected": "Signal-only projected data statistics; fixed MC response/corrections; systematics excluded",
}


def trapz(values, x, axis=-1):
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    return integrate(values, x, axis=axis)


def numeric(frame, name):
    return frame[name].to_numpy(dtype=float)


def good_bool(frame, name):
    return frame[name].astype(str).str.lower().isin(["true", "1"]).to_numpy()


def ratio_limits(values, errors):
    """Show complete conditional error bars, including broad projections."""
    finite = np.isfinite(values) & np.isfinite(errors)
    if not np.any(finite):
        return 0.7, 1.3
    low = min(0.7, float(np.min(values[finite] - errors[finite])))
    high = max(1.3, float(np.max(values[finite] + errors[finite])))
    padding = 0.05 * (high - low)
    return low - padding, high + padding


def absolute_ratio_limits(values, lower, upper):
    finite = np.isfinite(values) & np.isfinite(lower) & np.isfinite(upper)
    if not np.any(finite):
        return 0.7, 1.3
    low = min(0.7, float(np.min(lower[finite])), float(np.min(values[finite])))
    high = max(1.3, float(np.max(upper[finite])), float(np.max(values[finite])))
    padding = 0.05 * (high - low)
    return low - padding, high + padding


def phase_label(args, mc_reference=False):
    label = (rf"ep, {args.energy_label}  |  $W_{{\gamma p}}\in[{args.w_min:g},{args.w_max:g}]$ GeV"
             rf"  |  $Q^2\leq{args.q2_cut:g}$ GeV$^2$")
    if mc_reference:
        return label + "  |  MC reference"
    luminosity = args.projection_lumi_invfb if args.stat_mode == "projected" else args.lumi_invfb
    return label + rf"  |  $L_{{ep}}={luminosity:g}$ fb$^{{-1}}$"


def checked_covariance(covariance, size):
    covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (size, size) or not np.all(np.isfinite(covariance)):
        raise ValueError("Covariance must be a finite square matrix matching the spectrum bins")
    scale = max(float(np.max(np.abs(covariance))), np.finfo(float).tiny)
    if not np.allclose(covariance, covariance.T, rtol=1e-8, atol=scale * 1e-10):
        raise ValueError("Input covariance is not symmetric")
    covariance = 0.5 * (covariance + covariance.T)
    if np.linalg.eigvalsh(covariance).min() < -scale * 1e-8:
        raise ValueError("Input covariance has a materially negative eigenvalue")
    return covariance


def unit_area(sample):
    widths = sample["high"] - sample["low"]
    area = float(np.sum(sample["y"] * widths))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Shape normalization requires a positive finite area")
    covariance = checked_covariance(sample["covariance"], len(widths))
    jacobian = np.eye(len(widths)) / area - np.outer(sample["y"], widths) / area**2
    normalized = {key: value.copy() for key, value in sample.items()}
    normalized["y"] = sample["y"] / area
    normalized["covariance"] = jacobian @ covariance @ jacobian.T
    normalized["err"] = np.sqrt(np.maximum(np.diag(normalized["covariance"]), 0))
    # Keep original upstream toys for the nonlinear transform. Its peak
    # normalization cancels their individual areas exactly, without a second
    # Gaussian approximation to the unit-area covariance.
    return normalized, area


def check_sample(sample):
    for value in sample.values():
        if not np.all(np.isfinite(value)):
            raise ValueError("Nonfinite values in the momentum-transfer shape")
    if np.any(sample["low"] < 0) or np.any(sample["high"] <= sample["low"]):
        raise ValueError("The Fourier--Bessel shape requires nonnegative ordered t bins")
    if np.any(sample["y"] < 0) or np.any(sample["err"] < 0):
        raise ValueError("Negative spectrum density or uncertainty")
    if not np.allclose(sample["high"][:-1], sample["low"][1:], rtol=0, atol=1e-9):
        raise ValueError("Momentum-transfer bins must be contiguous")
    checked_covariance(sample["covariance"], len(sample["y"]))


def bessel_kernel(sample, b_grid_fm, nsub=80):
    """Same per-bin quadrature as plot_t_spectra_and_impact_parameter.py."""
    kernels = []
    for low, high in zip(sample["low"], sample["high"]):
        delta = np.linspace(np.sqrt(low), np.sqrt(high), nsub)
        integrand = delta[None, :] * j0(b_grid_fm[:, None] * delta[None, :] / HBARC_GEV_FM)
        kernels.append(trapz(integrand, delta, axis=1))
    return np.asarray(kernels)


def normalized_toy_statistics(sample, coverage=DEFAULT_COVERAGE):
    """Central sampling band after normalizing each complete upstream toy."""
    toys = np.asarray(sample["toys"], dtype=float)
    widths = sample["high"] - sample["low"]
    if toys.ndim != 2 or toys.shape[1] != len(widths) or not np.all(np.isfinite(toys)):
        raise ValueError("Shape toys must be finite complete spectra matching the bins")
    if np.any(toys < 0):
        raise ValueError("Unit-area projected shape toys must be nonnegative")
    areas = toys @ widths
    valid = np.isfinite(areas) & (areas > 0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("Fewer than two positive-area shape toys; sampling band is not estimable")
    normalized = toys[valid] / areas[valid, None]
    lower, upper = sampling_interval(normalized, coverage=coverage, axis=0)
    return {"lower": lower, "upper": upper, "median": np.median(normalized, axis=0),
            "covariance": np.atleast_2d(np.cov(normalized, rowvar=False, ddof=1)),
            "audit": {"input_toys": len(toys), "retained_toys": int(valid.sum()),
                      "zero_area_toys_excluded": int((~valid).sum()), "interval_coverage": coverage,
                      "interpretation": "Pointwise central sampling interval of individually area-normalized unfolded toys; conditional on positive total area; not a Garwood confidence interval"}}


def profile_with_uncertainty(sample, b_grid_fm, ntoys, seed, coverage=DEFAULT_COVERAGE):
    kernel = bessel_kernel(sample, b_grid_fm)
    nominal = np.abs(np.sqrt(sample["y"]) @ kernel)
    peak = float(np.max(nominal))
    if peak <= 0:
        raise ValueError("Fourier--Bessel profile has no positive normalization")
    if "toys" in sample:
        toys = np.asarray(sample["toys"], dtype=float)
        source = "Original correlated unfolded spectra from upstream toys"
    elif np.any(sample["covariance"]):
        if ntoys < 2:
            raise ValueError("At least two profile toys are required when input covariance is nonzero")
        rng = np.random.default_rng(seed)
        covariance = checked_covariance(sample["covariance"], len(sample["y"]))
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        toys = sample["y"] + rng.standard_normal((ntoys, len(sample["y"]))) @ (
            eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))).T
        source = "Correlated Gaussian approximation to full input covariance"
    else:
        toys = np.tile(sample["y"], (max(2, ntoys), 1))
        source = "Deterministic reference (input covariance is zero)"
    if toys.ndim != 2 or toys.shape[1] != len(sample["y"]) or not np.all(np.isfinite(toys)):
        raise ValueError("Profile input toys must be finite and match the spectrum bins")
    negative_bins = int(np.count_nonzero(toys < 0))
    negative_toys = int(np.count_nonzero(np.any(toys < 0, axis=1)))
    # sqrt(dN/dt) requires nonnegative bins. Record clipping explicitly; a
    # Gaussian fallback with clipping is only an approximate shape diagnostic.
    toy_profiles = np.abs(np.sqrt(np.maximum(toys, 0)) @ kernel)
    peaks = np.max(toy_profiles, axis=1)
    valid = np.isfinite(peaks) & (peaks > 0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("Fewer than two positive finite toy profiles; cannot estimate profile covariance")
    toy_profiles = toy_profiles[valid] / peaks[valid, None]
    covariance = np.atleast_2d(np.cov(toy_profiles, rowvar=False, ddof=1))
    lower, upper = sampling_interval(toy_profiles, coverage=coverage, axis=0)
    q16, q50, q84 = np.quantile(toy_profiles, [0.16, 0.5, 0.84], axis=0)
    audit = {"toy_source": source, "input_toys": len(toys), "retained_toys": int(valid.sum()),
             "zero_peak_toys_excluded": int((~valid).sum()), "negative_bins_clipped": negative_bins,
             "toys_with_negative_bins_clipped": negative_toys,
             "normalization": "Each complete toy is divided by its own maximum absolute profile",
             "nonnegative_condition": "Negative spectrum bins are clipped to zero before sqrt; zero-peak toys excluded",
             "interval_coverage": coverage,
             "band": "Pointwise central sampling interval at the requested coverage; covariance is sample covariance (ddof=1); not a Garwood confidence interval"}
    return {"profile": nominal / peak, "error": np.sqrt(np.maximum(np.diag(covariance), 0)),
            "covariance": covariance, "lower": lower, "upper": upper,
            "q16": q16, "q50": q50, "q84": q84, "audit": audit}


class Renderer:
    def __init__(self, args):
        self.args = args
        self.state = getattr(args, "state", "1S")
        if self.state not in ("1S", "2S", "3S"):
            raise ValueError(f"Unknown Upsilon state: {self.state}")
        self.state_lower = self.state.lower()
        self.include_1s_reference = self.state == "1S"
        self.args.count_interval = getattr(args, "count_interval", "median_garwood")
        self.args.interval_coverage = getattr(args, "interval_coverage", DEFAULT_COVERAGE)
        self.artifacts = args.run_artifacts.resolve()
        self.out = args.out_dir.resolve()
        self.out.mkdir(parents=True, exist_ok=True)
        self.uncertainty_note = UNCERTAINTY_NOTES[args.stat_mode]
        self.stat_label = "MC stat." if args.stat_mode == "mc" else "projected data stat."
        self.coverage_label = f"{100 * args.interval_coverage:.1f}%"
        self.count_label = (self.coverage_label + " median expected Garwood limits"
                            if args.count_interval == "median_garwood" else
                            self.coverage_label + " Gaussian-scaled SD intervals")
        self.sampling_label = (self.coverage_label + " central toy sampling band"
                              if args.count_interval == "median_garwood" else
                              self.coverage_label + " Gaussian-scaled covariance-SD band")
        self.audit = {"energy_key": args.energy_key, "energy_label": args.energy_label,
                      "state": self.state, "input_artifacts": str(self.artifacts),
                      "w_range_GeV": [args.w_min, args.w_max], "q2_cut_GeV2": args.q2_cut,
                      "source_lumi_invfb": args.lumi_invfb, "stat_mode": args.stat_mode,
                      "projection_lumi_invfb": args.projection_lumi_invfb,
                      "count_interval": args.count_interval, "interval_coverage": args.interval_coverage,
                      "active_lumi_invfb": args.projection_lumi_invfb if args.stat_mode == "projected" else args.lumi_invfb,
                      "uncertainty_model": self.uncertainty_note,
                      "renderer_luminosity_rescaling": False,
                      "source_files": [], "figures": [], "omitted": [], "notes": [],
                      "statistical_metadata_checks": [], "secondary_covariances": [],
                      "published_1S_references_enabled": self.include_1s_reference}
        if not self.include_1s_reference:
            self.audit["notes"].append("Published HERA/ZEUS 1S reference overlays are omitted for the selected excited state; no state rescaling is inferred.")

    def read(self, relative, optional=False):
        path = self.artifacts / relative
        if not path.is_file():
            if optional:
                self.audit["omitted"].append(str(relative))
                return None
            raise FileNotFoundError(path)
        self.audit["source_files"].append({"path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        frame = pd.read_csv(path)
        if "state" in frame:
            frame = frame[frame["state"] == self.state].copy()
        return frame

    def validate_statistics(self, frame, source, unfolding=False):
        if "stat_mode" not in frame or frame.empty:
            raise ValueError(f"{source}: statistical metadata missing; regenerate upstream outputs with --stat-mode")
        modes = set(frame["stat_mode"].astype(str))
        if modes != {self.args.stat_mode}:
            raise ValueError(f"{source}: upstream stat_mode {modes} does not match {self.args.stat_mode}")
        check = {"source": source, "stat_mode": self.args.stat_mode}
        if self.args.stat_mode == "projected":
            column = "yield_lumi_invfb" if unfolding else "projection_lumi_invfb"
            if column not in frame:
                raise ValueError(f"{source}: missing {column}; cannot verify projected luminosity")
            values = numeric(frame, column)
            if not np.allclose(values, self.args.projection_lumi_invfb, rtol=1e-9, atol=0):
                raise ValueError(f"{source}: upstream projected luminosity differs from requested target")
            check[column] = float(values[0])
            if "count_interval" not in frame or set(frame["count_interval"].astype(str)) != {self.args.count_interval}:
                raise ValueError(f"{source}: projected count-interval metadata does not match requested method")
            if "interval_coverage" not in frame or not np.allclose(numeric(frame, "interval_coverage"), self.args.interval_coverage, rtol=1e-10, atol=0):
                raise ValueError(f"{source}: projected interval coverage does not match requested coverage")
            check.update(count_interval=self.args.count_interval, interval_coverage=self.args.interval_coverage)
        elif unfolding and "yield_lumi_invfb" in frame:
            check["yield_lumi_invfb"] = numeric(frame, "yield_lumi_invfb").tolist()
            if "train_fraction" in frame:
                expected = self.args.lumi_invfb * (1 - numeric(frame, "train_fraction"))
                if not np.allclose(numeric(frame, "yield_lumi_invfb"), expected, rtol=1e-9, atol=0):
                    raise ValueError(f"{source}: MC luminosity does not match the requested source luminosity and split")
        elif "yield_lumi_invfb" in frame:
            values = numeric(frame, "yield_lumi_invfb")
            if not np.allclose(values, self.args.lumi_invfb, rtol=1e-9, atol=0):
                raise ValueError(f"{source}: MC luminosity does not match --lumi-invfb")
            check["yield_lumi_invfb"] = float(values[0])
        for column in ("uncertainty_model", "include_efficiency_mc_uncertainty"):
            if column in frame:
                check[column] = frame[column].astype(str).unique().tolist()
        self.audit["statistical_metadata_checks"].append(check)

    def extraction_note(self, frame):
        if self.args.stat_mode == "mc":
            return self.uncertainty_note
        note = "Projected " + self.count_label
        if ("include_efficiency_mc_uncertainty" in frame
                and np.any(good_bool(frame, "include_efficiency_mc_uncertainty"))):
            return note + "; correction-MC added in quadrature (approximation)"
        return note + "; fixed corrections"

    def display_interval(self, ax, frame, x, y_column, error_column,
                         low_column, high_column, mask=None, xerr=None, **kwargs):
        y = numeric(frame, y_column)
        selected = np.ones(len(frame), dtype=bool) if mask is None else np.asarray(mask, dtype=bool).copy()
        selected &= np.isfinite(x) & np.isfinite(y)
        if "unfolding_identifiable" in frame:
            identifiable = good_bool(frame, "unfolding_identifiable")
            selected &= identifiable
            if np.any(~identifiable):
                self.audit["notes"].append(f"{y_column}: omitted unidentifiable bins {np.flatnonzero(~identifiable).tolist()}; no zero-error nominal point displayed")
        if self.args.stat_mode == "projected":
            for column in (low_column, high_column):
                if column not in frame:
                    raise ValueError(f"Missing projected absolute interval column {column}; regenerate upstream results")
            lower, upper = numeric(frame, low_column), numeric(frame, high_column)
            selected &= np.isfinite(lower) & np.isfinite(upper)
            artist = draw_absolute_intervals(ax, x[selected], y[selected], lower[selected], upper[selected],
                                             xerr=None if xerr is None else xerr[selected], **kwargs)
        else:
            error = numeric(frame, error_column)
            selected &= np.isfinite(error)
            lower, upper = y - error, y + error
            artist = ax.errorbar(x[selected], y[selected], yerr=error[selected],
                                xerr=None if xerr is None else xerr[selected], **kwargs)
        self.audit.setdefault("displayed_intervals", []).append({
            "observable": y_column, "stat_mode": self.args.stat_mode,
            "lower_column": low_column if self.args.stat_mode == "projected" else error_column,
            "upper_column": high_column if self.args.stat_mode == "projected" else error_column,
            "displayed_bins": np.flatnonzero(selected).tolist(),
            "variance_changed_by_display_interval": False})
        return artist, selected, lower, upper

    def unfolding_metadata(self, prefix):
        frame = self.read(prefix + "_normalization.csv")
        self.validate_statistics(frame, prefix + "_normalization.csv", unfolding=True)
        if "covariance_valid" in frame and not np.all(good_bool(frame, "covariance_valid")):
            raise ValueError(f"{prefix}: upstream covariance is invalid; rerun with at least two toys")
        return frame

    def unfolding_covariance(self, prefix, size):
        frame = self.read(prefix + "_unfolded_covariance.csv")
        if frame.shape == (size, size + 1):
            frame = frame.iloc[:, 1:]
        return checked_covariance(frame.to_numpy(dtype=float), size)

    def unfolding_toys(self, prefix, edges, normalization):
        path = self.artifacts / (prefix + "_unfolded_toys.npz")
        if not path.is_file():
            if self.args.stat_mode == "projected":
                raise FileNotFoundError(f"{path}: projected nonlinear profiles require original unfolded toys")
            self.audit["notes"].append(f"{path.name} missing: MC profile uses a correlated Gaussian covariance approximation")
            return None
        self.audit["source_files"].append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        with np.load(path, allow_pickle=False) as archive:
            mode = str(archive["stat_mode"].item())
            lumi = float(archive["yield_lumi_invfb"].item())
            expected_lumi = float(normalization["yield_lumi_invfb"].iloc[0])
            if mode != self.args.stat_mode or not np.isclose(lumi, expected_lumi, rtol=1e-9, atol=0):
                raise ValueError(f"{path}: toy mode/luminosity does not match the upstream summary")
            if self.args.stat_mode == "projected":
                if ("count_interval" not in archive or str(archive["count_interval"].item()) != self.args.count_interval
                        or "interval_coverage" not in archive or not np.isclose(float(archive["interval_coverage"].item()), self.args.interval_coverage, rtol=1e-10, atol=0)):
                    raise ValueError(f"{path}: toy interval metadata does not match requested method/coverage")
            if not np.array_equal(archive["bin_edges"], edges):
                if not np.allclose(archive["bin_edges"], edges, rtol=0, atol=1e-10):
                    raise ValueError(f"{path}: toy bins do not match the spectrum")
            toys = archive["unfolded_visible"].copy()
        if toys.ndim != 2 or toys.shape[1] != len(edges) - 1 or len(toys) < 2:
            raise ValueError(f"{path}: at least two complete unfolded spectra are required")
        return toys

    def export_covariance(self, stem, covariance, coordinates):
        path = self.out / (stem + "_covariance.csv")
        pd.DataFrame(covariance, index=coordinates, columns=coordinates).to_csv(path, index_label="bin_or_coordinate")
        self.audit["secondary_covariances"].append(str(path))

    def save(self, fig, stem, title, extra=None, handles=None, labels=None, mc_reference=False):
        if extra is None:
            extra = self.uncertainty_note
        subtitle = phase_label(self.args, mc_reference) + ("\n" + extra if extra else "")
        paths = finish(fig, self.out, stem, title, subtitle, handles, labels, legend_ncol=3)
        self.audit["figures"].append({"stem": stem, "files": paths})

    def w_extraction(self):
        key = self.args.energy_key
        frame = self.read(f"cross_sections/{key}_sigma_gamma_p_{self.state_lower}.csv")
        self.validate_statistics(frame, f"cross_sections/{key}_sigma_gamma_p_{self.state_lower}.csv")
        x, dx = numeric(frame, "W_center"), 0.5 * (numeric(frame, "W_high") - numeric(frame, "W_low"))
        valid = good_bool(frame, "plotted") & np.isfinite(numeric(frame, "sigma_gamma_p_pb"))
        fig, ax = plt.subplots(figsize=(11.8, 7.0))
        h1, _, _, _ = self.display_interval(ax, frame, x, "sigma_gamma_p_pb", "sigma_gamma_p_err_pb",
                                            "sigma_gamma_p_interval_low_pb", "sigma_gamma_p_interval_high_pb",
                                            mask=valid, xerr=dx, color=BLUE, fmt="o", capsize=3)
        handles, labels = [h1], [f"ePIC {self.args.energy_label}"]
        reference_note = ""
        if self.include_1s_reference:
            # These published values describe 1S only; no excited-state extrapolation.
            h2 = ax.errorbar([95., 175.], [160., 321.], xerr=[35., 45.], yerr=[70., 144.],
                            color="#777777", fmt="^", capsize=3, markersize=6)
            handles.append(h2)
            labels.append(r"HERA exclusive $\Upsilon(1S)$")
            reference_note = "\nHERA bars retain the published experimental uncertainty"
        ax.set(xlabel=r"$W_{\gamma p}$ [GeV]", ylabel=rf"$\sigma(\gamma p\to\Upsilon({self.state})p)$ [pb]",
               xlim=(max(0, self.args.w_min - 5), max(225, self.args.w_max + 5) if self.include_1s_reference else self.args.w_max + 5))
        ax.set_yscale("log")
        self.save(fig, "cross_section_w", rf"Exclusive $\Upsilon({self.state})$ photoproduction",
                  self.extraction_note(frame) + reference_note, handles, labels)

        fig, ax = plt.subplots(figsize=(11.8, 6.6))
        h = ax.errorbar(x, numeric(frame, "eff_reco_over_true"), xerr=dx,
                       yerr=numeric(frame, "eff_err_mc_approx"), fmt="o-", color=BLUE, capsize=3)
        ax.set(xlabel=r"$W_{\gamma p}$ [GeV]", ylabel=r"$N_{\rm reco}(W_{\rm reco})/N_{\rm true}(W_{\rm true})$",
               xlim=(self.args.w_min, self.args.w_max), ylim=(0, 1))
        self.save(fig, "w_reconstruction_efficiency", rf"$\Upsilon({self.state})$ reconstruction ratio versus $W$",
                  "Ratio in matching W bins includes migration; correlated delta-method MC uncertainty",
                  [h], [rf"$\Upsilon({self.state})$"], mc_reference=True)

        ratio = numeric(frame, "closure_unfold_over_truth")
        fig, ax = plt.subplots(figsize=(11.8, 6.6))
        h = ax.errorbar(x[valid], ratio[valid], xerr=dx[valid], fmt="o", color=BLUE, capsize=3)
        ax.axhline(1, color="#666666", linestyle="--")
        ax.set(xlabel=r"$W_{\gamma p}$ [GeV]", ylabel="Corrected / truth", ylim=(0.8, 1.2))
        self.save(fig, "bin_by_bin_closure", rf"$\Upsilon({self.state})$ bin-by-bin closure",
                  "Same-sample algebraic identity; no independent closure uncertainty is assigned",
                  [h], ["Corrected / truth"], mc_reference=True)
        self.audit["notes"].append("Bin-by-bin closure equals unity algebraically when the same sample supplies reco, truth, and correction. No independent numerator/denominator error is invented; projected data errors are not an MC-closure uncertainty.")

    def differential(self):
        key = self.args.energy_key
        summary = self.read(f"differential/{key}_{self.state_lower}_b_slope_summary.csv")
        self.validate_statistics(summary, f"differential/{key}_{self.state_lower}_b_slope_summary.csv")
        fit = summary.iloc[0]
        b = float(fit["b_slope_GeVminus2"])
        self.audit["b_fit"] = {name: float(fit[name]) for name in ["b_slope_GeVminus2", "b_slope_err_low",
                             "b_slope_err_high", "chi2", "ndof", "b_fit_t_min", "b_fit_t_max"]}
        statistic = str(fit.get("fit_statistic", "chi2"))
        poisson_fit = statistic.startswith("poisson_deviance")
        constrained_fit = "constraint" in statistic
        statistic_value = float(fit.get("fit_statistic_value", fit["chi2"]))
        statistic_label = (r"D+P" if constrained_fit else "D") if poisson_fit else r"\chi^2"
        fit_label = "Reco-Poisson profile fit" if poisson_fit else "Finite-bin GLS fit"
        fit_coverage = float(fit.get("fit_interval_coverage", DEFAULT_COVERAGE))
        self.audit["b_fit"].update(fit_statistic=statistic, fit_statistic_value=statistic_value,
                                  interval_coverage=fit_coverage if poisson_fit else None,
                                  fit_calibration_model=str(fit.get("fit_calibration_model", "unspecified")),
                                  deviance=float(fit.get("deviance", np.nan)),
                                  constraint_penalty=float(fit.get("constraint_penalty", 0)))
        for variable, stem, title, xlabel, ylabel in [
            ("pT2", "pt2_spectra", rf"$\Upsilon({self.state})$ differential $p_T^2$ spectrum", r"$p_{T,\Upsilon}^{2}$ [GeV$^2$]", r"$d\sigma_{\gamma p}/dp_T^2$ [pb/GeV$^2$]"),
            ("y", "rapidity_coverage", rf"$\Upsilon({self.state})$ rapidity coverage", r"$y_\Upsilon$", r"$d\sigma_{\gamma p}/dy$ [pb]")]:
            frame = self.read(f"differential/{key}_{self.state_lower}_{self.state}_{variable}.csv")
            self.validate_statistics(frame, f"differential/{key}_{self.state_lower}_{self.state}_{variable}.csv")
            valid = good_bool(frame, "passes_min_eff") & np.isfinite(numeric(frame, "dsigma_dx_pb"))
            x = numeric(frame, "x_center")
            y = numeric(frame, "dsigma_dx_pb")
            dy = numeric(frame, "dsigma_dx_pb_err")
            if str(frame["output_convention"].iloc[0]) == "visible_decay":
                ylabel = r"$B_{ee}$ " + ylabel
            fig, ax = plt.subplots(figsize=(11.8, 6.8))
            h, _, _, _ = self.display_interval(ax, frame, x, "dsigma_dx_pb", "dsigma_dx_pb_err",
                                               "dsigma_dx_pb_interval_low", "dsigma_dx_pb_interval_high",
                                               mask=valid, xerr=0.5 * numeric(frame, "bin_width"),
                                               fmt="o", color=BLUE, capsize=3, markersize=6)
            handles, labels = [h], [rf"$\Upsilon({self.state})$ (" + self.stat_label + ")"]
            extra = self.extraction_note(frame) + "; bin-by-bin correction"
            if variable == "pT2":
                low, high = float(fit["b_fit_t_min"]), float(fit["b_fit_t_max"])
                if np.isfinite(b) and b > 0:
                    xx = np.linspace(low, high, 500)
                    norm = np.exp(-b * low) - np.exp(-b * high)
                    yy = float(fit["fit_integrated_sigma_pb"]) * b * np.exp(-b * xx) / norm
                    line, = ax.plot(xx, yy, color=BLUE, alpha=0.7)
                    handles.append(line)
                    labels.append(rf"{fit_label}, $b={b:.3f}^{{+{fit['b_slope_err_high']:.3f}}}_{{-{fit['b_slope_err_low']:.3f}}}$ GeV$^{{-2}}$")
                    extra += rf"; ${statistic_label}/\mathrm{{ndf}}={statistic_value:.1f}/{int(fit['ndof'])}$"
                    predictions = frame[["x_low", "x_high", "x_center", "bin_width"]].copy()
                    predictions["fit_bin_average_pb_per_GeV2"] = float(fit["fit_integrated_sigma_pb"]) * (
                        np.exp(-b * numeric(frame, "x_low")) - np.exp(-b * numeric(frame, "x_high"))) / norm / numeric(frame, "bin_width")
                    predictions.to_csv(self.out / "pt2_fit_bin_averages.csv", index=False)
                ax.set_yscale("log")
                ax.set_xlim(0, float(frame["x_high"].max()))
            else:
                ax.set_ylim(bottom=0)
                ax.set_xlim(float(frame["x_low"].min()), float(frame["x_high"].max()))
            ax.set(xlabel=xlabel, ylabel=ylabel)
            self.save(fig, stem, title, extra, handles, labels)

        fig, ax = plt.subplots(figsize=(11.8, 6.6))
        h1 = ax.errorbar([0], [b], yerr=np.array([[fit["b_slope_err_low"]], [fit["b_slope_err_high"]]]),
                        color=BLUE, fmt="o", capsize=5, markersize=7)
        handles, labels = [h1], [f"ePIC {fit_label.lower()}"]
        tick_positions, tick_labels = [0], [f"ePIC {self.args.energy_label}\n" + r"$p_T^2$ fit"]
        if self.include_1s_reference:
            h2 = ax.errorbar([1], [4.3], yerr=[[1.3], [2.0]], color="#555555", fmt="s", capsize=5, markersize=7)
            handles.append(h2)
            labels.append("ZEUS published statistical uncertainty")
            tick_positions.append(1)
            tick_labels.append(r"ZEUS $\Upsilon(1S)$" + "\n" + r"$|t|$ fit")
        ax.set_xticks(tick_positions, tick_labels)
        ax.set(xlim=(-0.5, 1.5 if self.include_1s_reference else 0.5), ylim=(0, 7), ylabel=r"Exponential slope $b$ [GeV$^{-2}$]")
        self.save(fig, "b_slope_comparison", rf"$\Upsilon({self.state})$ exponential-slope comparison",
                  (f"Projected reco-Poisson profile interval ({100 * fit_coverage:.1f}%)" if poisson_fit else self.uncertainty_note)
                  + rf"; ePIC fit ${statistic_label}/\mathrm{{ndf}}={statistic_value:.1f}/{int(fit['ndof'])}$",
                  handles, labels)
        if self.include_1s_reference:
            self.audit["notes"].append("ZEUS slope reference4.3(+2.0,-1.3) GeV^-2 is taken from the existing framework reference. No approximate light-meson reference points are reproduced.")
        if constrained_fit:
            self.audit["notes"].append("The projected slope uses correlated Gaussian log-correction nuisance parameters (D+P is Poisson deviance plus their constraint penalty). This differs from the approximate quadrature combination used for displayed count intervals.")

    def momentum_shapes(self):
        key = self.args.energy_key
        prefix = f"t/{key}_pT2"
        frame = self.read(prefix + "_t_unfold_summary.csv")
        self.validate_statistics(frame, prefix + "_t_unfold_summary.csv", unfolding=True)
        normalization = self.unfolding_metadata(prefix)
        if set(frame["method"]) != {"pT2"}:
            raise ValueError("Expected the pT2-estimator unfolding summary")
        if "unfolding_identifiable" in frame and not np.all(good_bool(frame, "unfolding_identifiable")):
            unsupported = np.flatnonzero(~good_bool(frame, "unfolding_identifiable")).tolist()
            raise ValueError(f"Cannot form a full-support normalized shape/profile with unidentifiable t bins {unsupported}; use explicitly supported narrower or merged binning")
        samples = {
            "epic": {"low": numeric(frame, "t_low"), "high": numeric(frame, "t_high"),
                     "y": numeric(frame, "dNdt_unfolded_visible"), "err": numeric(frame, "dNdt_unfolded_visible_err")}}
        if self.include_1s_reference:
            reference = self.args.zeus_csv.resolve()
            z = pd.read_csv(reference)
            self.audit["source_files"].append({"path": str(reference), "sha256": hashlib.sha256(reference.read_bytes()).hexdigest()})
            samples["zeus"] = {"low": numeric(z, "t_low"), "high": numeric(z, "t_high"),
                               "y": numeric(z, "dNdt_GeVminus2"),
                               "err": 0.5 * (numeric(z, "err_low") + numeric(z, "err_high"))}
            samples["zeus"]["covariance"] = np.diag(samples["zeus"]["err"]**2)
        widths = samples["epic"]["high"] - samples["epic"]["low"]
        integrated_covariance = self.unfolding_covariance(prefix, len(widths))
        samples["epic"]["covariance"] = integrated_covariance / np.outer(widths, widths)
        edges = np.r_[samples["epic"]["low"], samples["epic"]["high"][-1]]
        toys = self.unfolding_toys(prefix, edges, normalization)
        if toys is not None:
            toy_covariance = np.atleast_2d(np.cov(toys, rowvar=False, ddof=1))
            if not np.allclose(toy_covariance, integrated_covariance, rtol=1e-6, atol=max(float(np.max(np.abs(integrated_covariance))), 1e-30) * 1e-8):
                raise ValueError("Upstream toy covariance does not match the saved unfolded covariance")
            samples["epic"]["toys"] = toys / widths
        for sample in samples.values():
            check_sample(sample)
        supports = [(s["low"][0], s["high"][-1]) for s in samples.values()]
        if self.include_1s_reference and not np.allclose(supports[0], supports[1], rtol=0, atol=1e-9):
            raise ValueError("ZEUS and ePIC spectra must have common support for the shape overlay")
        areas = {}
        for name in samples:
            samples[name], areas[name] = unit_area(samples[name])
            self.export_covariance("momentum_transfer_" + name + "_shape", samples[name]["covariance"],
                                   [f"{lo:g}-{hi:g}" for lo, hi in zip(samples[name]["low"], samples[name]["high"])])
        epic, zeus = samples["epic"], samples.get("zeus")
        shape_toys = None
        if "toys" in epic:
            shape_toys = normalized_toy_statistics(epic, self.args.interval_coverage)
            self.export_covariance("momentum_transfer_epic_shape_toy", shape_toys["covariance"],
                                   [f"{lo:g}-{hi:g}" for lo, hi in zip(epic["low"], epic["high"])])
            self.audit["normalized_shape_toys"] = shape_toys["audit"]
        fig, ax = plt.subplots(figsize=(11.8, 6.8))
        if zeus is not None:
            zh = ax.stairs(zeus["y"], np.r_[zeus["low"], zeus["high"][-1]], color="#444444", linestyle=":", linewidth=2.8)
        if self.args.stat_mode == "projected":
            if shape_toys is None:
                raise ValueError("Projected normalized shape requires original complete unfolding toys")
            eh = draw_absolute_intervals(ax, 0.5 * (epic["low"] + epic["high"]), epic["y"],
                                         shape_toys["lower"], shape_toys["upper"],
                                         xerr=0.5 * (epic["high"] - epic["low"]), fmt="o", color=BLUE, capsize=3)
            shape_note = f"Projected {self.coverage_label} central toy sampling intervals; each toy normalized to unit area; fixed response"
        else:
            eh = ax.errorbar(0.5 * (epic["low"] + epic["high"]), epic["y"],
                             xerr=0.5 * (epic["high"] - epic["low"]), yerr=epic["err"], fmt="o", color=BLUE, capsize=3)
            shape_note = "Unit-area shape; finite support; MC precision from full normalization covariance"
            if zeus is not None:
                shape_note = shape_note.replace("finite support", "common finite support")
        ax.set_yscale("log")
        ax.set(xlabel=r"$|t|$ [GeV$^2$]", ylabel=r"Unit-area $dN/d|t|$ [GeV$^{-2}$]", xlim=supports[0])
        handles, labels = [eh], [f"ePIC {self.args.energy_label}, $p_T^2$ estimator"]
        if zeus is not None:
            handles.append(zh)
            labels.append("ZEUS exclusive component (digitized step)")
        shape_title = r"$p_T^2$-based unfolded $|t|$ shape"
        if not self.include_1s_reference:
            shape_title = rf"$\Upsilon({self.state})$ " + shape_title
        self.save(fig, "momentum_transfer_comparison", shape_title,
                  shape_note, handles, labels)
        records = []
        for name, sample in samples.items():
            for i, (lo, hi, value, err) in enumerate(zip(sample["low"], sample["high"], sample["y"], sample["err"])):
                records.append({"sample": name, "state": self.state, "t_low": lo, "t_high": hi, "density": value, "density_error": err,
                                "stat_mode": self.args.stat_mode if name == "epic" else "published_reference",
                                "error_convention": "sqrt(diag(J C J.T)); first-order normalization covariance; not a confidence interval",
                                "density_sampling_interval_low": shape_toys["lower"][i] if name == "epic" and shape_toys is not None else np.nan,
                                "density_sampling_interval_high": shape_toys["upper"][i] if name == "epic" and shape_toys is not None else np.nan,
                                "density_toy_standard_deviation": np.sqrt(max(shape_toys["covariance"][i, i], 0)) if name == "epic" and shape_toys is not None else np.nan,
                                "interval_coverage": self.args.interval_coverage if name == "epic" and shape_toys is not None else np.nan,
                                "displayed_error_convention": "central toy sampling interval" if name == "epic" and self.args.stat_mode == "projected" else "covariance standard deviation"})
        pd.DataFrame(records).to_csv(self.out / "momentum_transfer_shapes.csv", index=False)

        bgrid = np.linspace(0, 2.5, 160)
        profile_coverage = self.args.interval_coverage if self.args.stat_mode == "projected" else .68
        epic_result = profile_with_uncertainty(epic, bgrid, self.args.n_profile_toys, 17, profile_coverage)
        zeus_result = (profile_with_uncertainty(zeus, bgrid, self.args.n_profile_toys, 23, .68)
                       if zeus is not None else None)
        profile, error = epic_result["profile"], epic_result["error"]
        results = {"epic": epic_result}
        if zeus_result is not None:
            results["zeus"] = zeus_result
        for name, result in results.items():
            self.export_covariance("impact_parameter_" + name, result["covariance"], bgrid)
        fig, ax = plt.subplots(figsize=(11.8, 6.8))
        eh, = ax.plot(bgrid, profile, color=BLUE, linewidth=2.5)
        ax.fill_between(bgrid, epic_result["lower"], epic_result["upper"], color=BLUE, alpha=0.15)
        handles, labels = [eh], [f"ePIC ({self.stat_label}; {100 * profile_coverage:.1f}% central sampling band)"]
        if zeus_result is not None:
            zh, = ax.plot(bgrid, zeus_result["profile"], color="#444444", linestyle=":", linewidth=2.8)
            if np.any(zeus_result["error"]):
                ax.fill_between(bgrid, zeus_result["lower"], zeus_result["upper"], color="#444444", alpha=0.1)
            handles.append(zh)
            labels.append("ZEUS exclusive step transform")
        ax.set(xlabel=r"$b_T$ [fm]", ylabel=r"Peak-normalized $|F(b_T)|$", xlim=(0, 2.5), ylim=(0, 1.08))
        profile_title = "Fourier–Bessel transverse profiles"
        if not self.include_1s_reference:
            profile_title = rf"$\Upsilon({self.state})$ " + profile_title
        self.save(fig, "impact_parameter_profiles", profile_title,
                  rf"Shape diagnostic; finite $|t|\in[{supports[0][0]:g},{supports[0][1]:g}]$ GeV$^2$ support; correlated toys, each peak-normalized",
                  handles, labels)
        profile_columns = {"bT_fm": bgrid, "state": self.state, "epic_profile": profile, "epic_profile_error": error,
                      "epic_profile_q16": epic_result["q16"], "epic_profile_median": epic_result["q50"],
                      "epic_profile_q84": epic_result["q84"],
                      "epic_profile_interval_low": epic_result["lower"], "epic_profile_interval_high": epic_result["upper"],
                           "epic_profile_interval_coverage": profile_coverage, "stat_mode": self.args.stat_mode}
        if zeus_result is not None:
            profile_columns.update(zeus_profile=zeus_result["profile"], zeus_profile_error=zeus_result["error"],
                                   zeus_profile_q16=zeus_result["q16"], zeus_profile_q84=zeus_result["q84"])
        pd.DataFrame(profile_columns).to_csv(self.out / "impact_parameter_profiles.csv", index=False)
        self.audit["transverse_profile"] = {
            "formula": "abs(integral dDelta Delta J0(bT*Delta/hbarc)*sqrt(dN/dt)), Delta=sqrt(t)",
            "normalization": "Nominal profile and every complete correlated toy divided by their respective own peaks",
            "sample_areas_before_normalization": areas, "t_support_GeV2": list(supports[0]),
            "integration_subpoints_per_bin": 80, "hbarc_GeV_fm": HBARC_GEV_FM,
            "gaussian_fallback_toys_requested": self.args.n_profile_toys, "fallback_seed_epic": 17,
            "fallback_seed_zeus": 23 if zeus_result is not None else None,
            "epic_toy_audit": epic_result["audit"], "zeus_toy_audit": zeus_result["audit"] if zeus_result is not None else None,
            "shape_covariance": "J C J.T with J=I/area-outer(density,bin_width)/area^2; preserves normalization correlations",
            "uncertainty": "Original complete unfolded toys where available; per-toy nonlinear transform and peak normalization. Empirical full profile covariance, requested-coverage central sampling bounds, and legacy q16/q50/q84 exported; no Garwood construction is applied to transformed spectra.",
            "conditional_model": "Data-only projected statistics with fixed response/prior; finite response MC, backgrounds, flux, luminosity and other systematics excluded. MC mode uses upstream weighted-template toys.",
            "interpretation": "Shape diagnostic only. Coarse bins, missing amplitude phase, finite t support and taking modulus can produce artifacts; no absolute density, proton radius, or gluon-density constraint is inferred.",
            "reference_uncertainty": ("Exclusive ZEUS step CSV carries zero errors; no reference band is invented."
                                      if zeus is not None else "No applicable excited-state ZEUS reference; omitted.")}

    def unfolding(self, coordinate=None):
        key = self.args.energy_key
        if coordinate not in (None, "w", "t"):
            raise ValueError("Unfolding coordinate must be w or t")
        coordinates = ("w", "t") if coordinate is None else (coordinate,)
        closure_note = (self.sampling_label + "; ratio truth fixed; not confidence limits"
                        if self.args.stat_mode == "projected" else
                        "Ratio bars: numerator only, truth fixed; no joint same-event closure uncertainty")
        if "w" in coordinates:
            w = self.read(f"unfolding/{key}_sigma_gamma_p_w_unfolded_summary.csv")
            self.validate_statistics(w, f"unfolding/{key}_sigma_gamma_p_w_unfolded_summary.csv", unfolding=True)
            self.unfolding_metadata(f"unfolding/{key}_sigma_gamma_p_w")
            x, dx = numeric(w, "W_center"), 0.5 * (numeric(w, "W_high") - numeric(w, "W_low"))
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(11.8, 8.2), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
            uh, _, _, _ = self.display_interval(ax, w, x, "sigma_unfolded_pb", "sigma_unfolded_err_pb",
                                                "sigma_unfolded_interval_low_pb", "sigma_unfolded_interval_high_pb",
                                                xerr=dx, fmt="o", color=BLUE, capsize=3)
            th = ax.errorbar(x, numeric(w, "sigma_truth_test_pb"), xerr=dx, yerr=numeric(w, "sigma_truth_test_err_pb"), fmt="s", color="#555555", mfc="white", capsize=3)
            ax.set_yscale("log"); ax.set_ylabel(rf"$\sigma_{{\gamma p\to\Upsilon({self.state})p}}$ [pb]")
            _, selected, lower, upper = self.display_interval(rax, w, x, "closure_unfolded_over_truth", "closure_ratio_err",
                                                              "closure_ratio_interval_low", "closure_ratio_interval_high",
                                                              xerr=dx, fmt="o", color=BLUE)
            rax.axhline(1, color="#777777", linestyle="--")
            rax.set(xlabel=r"$W_{\gamma p}$ [GeV]", ylabel="Unfolded / truth",
                    ylim=absolute_ratio_limits(numeric(w, "closure_unfolded_over_truth")[selected], lower[selected], upper[selected]))
            self.save(fig, "w_unfolding_closure", r"$W$ unfolding with a held-out MC template",
                      self.uncertainty_note + "\n" + closure_note,
                      [uh, th], [f"Unfolded ({self.stat_label})", "Truth template (MC precision)"])
        for coordinate in coordinates:
            if coordinate == "w":
                prefix = f"unfolding/{key}_sigma_gamma_p_w"
                low, high = numeric(w, "W_low"), numeric(w, "W_high")
                title = "$W$ response"
            else:
                prefix = f"t/{key}_pT2"
                title = r"$p_T^2$-estimator response"
                t = self.read(f"t/{key}_pT2_t_unfold_summary.csv")
                low, high = numeric(t, "t_low"), numeric(t, "t_high")
            response = self.read(prefix + "_response_prob.csv")
            if response.columns[0].startswith("Unnamed"):
                response = response.iloc[:, 1:]
            matrix = response.to_numpy(dtype=float)
            edges = np.r_[low, high[-1]]
            if matrix.shape != (len(low), len(low)):
                raise ValueError("Response matrix dimensions disagree with the summary bins")
            fig, ax = plt.subplots(figsize=(9.7, 7.8))
            image = ax.pcolormesh(edges, edges, matrix, cmap="Blues", vmin=0)
            fig.colorbar(image, ax=ax, label=r"$P(\mathrm{reco\ bin}\mid\mathrm{truth\ bin})$")
            label = r"$W$ [GeV]" if coordinate == "w" else r"$|t|$ [GeV$^2$]"
            ax.set(xlabel="Reconstructed " + label, ylabel="Truth " + label)
            self.save(fig, coordinate + "_response_probability", title,
                      "Fixed MC response; includes reconstruction inefficiency; no projected data fluctuations",
                      mc_reference=True)
        if "t" in coordinates:
            t = self.read(f"t/{key}_pT2_t_unfold_summary.csv")
            self.validate_statistics(t, f"t/{key}_pT2_t_unfold_summary.csv", unfolding=True)
            fig, ax = plt.subplots(figsize=(11.8, 6.6))
            h, selected, lower, upper = self.display_interval(ax, t, numeric(t, "t_center"), "closure_unfold_over_truth", "closure_err",
                                                             "closure_interval_low", "closure_interval_high",
                                                             xerr=0.5 * numeric(t, "bin_width"), fmt="o", color=BLUE, capsize=3)
            ax.axhline(1, color="#777777", linestyle="--")
            ax.set(xlabel=r"Truth $|t|$ [GeV$^2$]", ylabel="Unfolded / truth",
                   ylim=absolute_ratio_limits(numeric(t, "closure_unfold_over_truth")[selected], lower[selected], upper[selected]))
            self.save(fig, "t_unfolding_closure", r"$p_T^2$-estimator held-out closure diagnostic",
                      self.uncertainty_note + "\n" + closure_note, [h], [f"Unfolded / fixed truth ({self.stat_label})"])
        self.audit["notes"].append("W/t responses are trained on independent training subsets, but unfolded and truth test spectra share held-out events. Ratio intervals propagate only the numerator with the truth template fixed; no independent truth-denominator error or full joint closure uncertainty is claimed. Projected central toy bounds are sampling intervals, not Garwood confidence limits. Truth bars remain MC template precision in both modes.")

    def resolutions(self, coordinate=None):
        key = self.args.energy_key
        if coordinate not in (None, "w", "t"):
            raise ValueError("Resolution coordinate must be w or t")
        coordinates = ("w", "t") if coordinate is None else (coordinate,)
        if "w" in coordinates:
            w = self.read(f"resolution/{key}_w_W_resolution_summary.csv")
            for stem, column, error, title, ylabel in [
                ("w_resolution", "median_delta_W", "core_sigma_W", r"$W$ reconstruction resolution", r"$W_{\rm reco}-W_{\rm true}$ [GeV]"),
                ("w_relative_resolution", "median_rel_delta_W", "core_sigma_rel", r"Relative $W$ reconstruction resolution", r"$(W_{\rm reco}-W_{\rm true})/W_{\rm true}$")]:
                fig, ax = plt.subplots(figsize=(11.8, 6.6))
                valid = np.isfinite(numeric(w, column)) & np.isfinite(numeric(w, error))
                h = ax.errorbar(numeric(w, "W_center")[valid], numeric(w, column)[valid],
                                xerr=0.5 * (numeric(w, "W_high") - numeric(w, "W_low"))[valid],
                                yerr=numeric(w, error)[valid], fmt="o", color=BLUE, capsize=3)
                ax.axhline(0, color="#777777", linestyle="--")
                ax.set(xlabel=r"Truth $W_{\gamma p}$ [GeV]", ylabel=ylabel, xlim=(self.args.w_min, self.args.w_max))
                self.save(fig, stem, title,
                          "Median residual; vertical bars are half the central 68% interval, not errors on the median",
                          [h], [rf"$\Upsilon({self.state})$ selected MC sample"], mc_reference=True)
        if "t" in coordinates:
            t = self.read(f"t/{key}_pT2_t_resolution_summary.csv")
            fig, ax = plt.subplots(figsize=(11.8, 6.6))
            valid = np.isfinite(numeric(t, "median_delta_t")) & np.isfinite(numeric(t, "core_sigma_t"))
            h = ax.errorbar(numeric(t, "t_center")[valid], numeric(t, "median_delta_t")[valid],
                            xerr=0.5 * (numeric(t, "t_high") - numeric(t, "t_low"))[valid],
                            yerr=numeric(t, "core_sigma_t")[valid], fmt="o", color=BLUE, capsize=3)
            ax.axhline(0, color="#777777", linestyle="--")
            ax.set(xlabel=r"Truth $|t|$ [GeV$^2$]", ylabel=r"$p_{T,\Upsilon,\mathrm{reco}}^2-|t|_{\rm true}$ [GeV$^2$]")
            self.save(fig, "t_resolution", r"$p_T^2$-estimator resolution",
                      "Median residual; vertical bars are half the central 68% interval, not statistical uncertainty",
                      [h], [rf"$\Upsilon({self.state})$ selected MC sample"], mc_reference=True)

    def run(self):
        section = getattr(self.args, "section", "all")
        self.audit["section"] = section
        if section == "all":
            self.w_extraction()
            self.differential()
            self.momentum_shapes()
            self.unfolding()
            self.resolutions()
        else:
            sections = {
                "w_extraction": self.w_extraction,
                "differential": self.differential,
                "momentum_shapes": self.momentum_shapes,
                "w_unfolding": lambda: self.unfolding("w"),
                "t_unfolding": lambda: self.unfolding("t"),
                "w_resolution": lambda: self.resolutions("w"),
                "t_resolution": lambda: self.resolutions("t"),
            }
            sections[section]()
        audit_name = "note_physics_audit.json" if section == "all" else f"note_physics_{section}_audit.json"
        path = self.out / audit_name
        path.write_text(json.dumps(self.audit, indent=2) + "\n", encoding="utf-8")
        print(f"Rendered {len(self.audit['figures'])} physics figure families into {self.out}")
        print(f"Audit: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", choices=["1S", "2S", "3S"], default="1S")
    parser.add_argument("--run-artifacts", required=True, type=Path)
    parser.add_argument("--section", default="all",
                        choices=["all", "w_extraction", "differential", "momentum_shapes",
                                 "w_unfolding", "t_unfolding", "w_resolution", "t_resolution"],
                        help="Render only the selected figure family; default renders all families")
    parser.add_argument("--energy-key", default="9x130")
    parser.add_argument("--energy-label", default="9 × 130 GeV")
    parser.add_argument("--w-min", type=float, default=20)
    parser.add_argument("--w-max", type=float, default=70)
    parser.add_argument("--q2-cut", type=float, default=1)
    parser.add_argument("--lumi-invfb", type=float, default=10)
    parser.add_argument("--stat-mode", choices=["mc", "projected"], default="mc",
                        help="Must match upstream analysis outputs; the renderer never rescales them")
    parser.add_argument("--projection-lumi-invfb", type=float, default=None,
                        help="Projected luminosity to verify against upstream metadata (defaults to --lumi-invfb)")
    add_interval_arguments(parser)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--zeus-csv", type=Path,
                        default=Path(__file__).resolve().parents[2] / "data/reference/zeus_fig4_exclusive_upsilon_step.csv")
    parser.add_argument("--n-profile-toys", type=int, default=200,
                        help="Number of correlated Gaussian fallback toys; all original upstream toys are used when available")
    args = parser.parse_args()
    if args.projection_lumi_invfb is None:
        args.projection_lumi_invfb = args.lumi_invfb
    if not np.isfinite(args.lumi_invfb) or args.lumi_invfb <= 0:
        parser.error("--lumi-invfb must be finite and positive")
    if not np.isfinite(args.projection_lumi_invfb) or args.projection_lumi_invfb <= 0:
        parser.error("--projection-lumi-invfb must be finite and positive")
    if not np.isfinite(args.interval_coverage) or not 0 < args.interval_coverage < 1:
        parser.error("--interval-coverage must be finite and strictly between zero and one")
    if args.n_profile_toys < 0:
        parser.error("--n-profile-toys cannot be negative")
    set_style()
    Renderer(args).run()


if __name__ == "__main__":
    main()
