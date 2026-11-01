"""Mass reconstruction and tuple-conditional efficiencies from framework tables.

The simultaneous fit is performed here; the three independent upstream fits are
used only to calibrate fixed tails and relative widths.  No ROOT file, simulated
background, or previously saved projected fit uncertainty is used.
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.special import erf

from .common import (
    COLORS, MARKERS, FONT_SIZES, annotate, draw_intervals, expected_intervals, read_csv,
    read_json, resolve_artifacts, sample_metadata, save_figure, set_style, state_label, write_json,
)

MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
SHAPE_NAMES = ("mean_GeV", "sigma_core_GeV", "alpha_left", "n_left",
               "alpha_right", "n_right")


def _crystal_ball_integral(z, al, nl, ar, nr):
    """DSCB antiderivative, matching the framework's calibrated function."""
    z = np.asarray(z, dtype=float)
    left_total = np.exp(-.5 * al**2) * nl / al / (nl - 1)
    core_total = np.sqrt(np.pi / 2) * (erf(ar / np.sqrt(2)) + erf(al / np.sqrt(2)))
    right_total = np.exp(-.5 * ar**2) * nr / ar / (nr - 1)
    result = np.empty_like(z)
    left, right = z < -al, z > ar
    core = ~(left | right)
    result[left] = np.exp(nl * np.log(nl / al) - al**2 / 2 - np.log(nl - 1)
                          + (1 - nl) * np.log(nl / al - al - z[left]))
    result[core] = left_total + np.sqrt(np.pi / 2) * (
        erf(z[core] / np.sqrt(2)) + erf(al / np.sqrt(2)))
    result[right] = left_total + core_total + right_total - np.exp(
        nr * np.log(nr / ar) - ar**2 / 2 - np.log(nr - 1)
        + (1 - nr) * np.log(nr / ar - ar + z[right]))
    return result


def mass_probabilities(edges, shape):
    mean, sigma, al, nl, ar, nr = np.asarray(shape, float)
    if sigma <= 0 or min(al, ar) <= 0 or min(nl, nr) <= 1:
        raise ValueError("DSCB requires sigma, alpha > 0 and tail powers > 1.")
    cdf = _crystal_ball_integral((np.asarray(edges) - mean) / sigma, al, nl, ar, nr)
    prob = np.maximum(np.diff(cdf), 1e-16)
    return prob / prob.sum()


def poisson_residual(mean, model):
    """Expected (Asimov) Poisson-deviance residual; accepts fractional means."""
    mean, model = np.asarray(mean, float), np.maximum(np.asarray(model, float), 1e-14)
    term = model - mean
    pos = mean > 0
    term[pos] += mean[pos] * np.log(mean[pos] / model[pos])
    return np.sign(mean - model) * np.sqrt(np.maximum(2 * term, 0))


def fit_combined_mass(edges, means, calibrations, float_common_shape=True):
    """Fit 3 signal yields, with fixed mass differences and calibrated tails.

    Fits expected histogram means, not a generated pseudo-experiment.  A common
    mass offset and width scale float by default.  State-specific tails and
    relative widths are fixed, with no MC-calibration covariance constraint.
    Returns a conditional descriptive fit, not a background-aware yield study.
    """
    states = tuple(MASSES)
    totals = np.array([np.sum(means[state]) for state in states], float)
    if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
        raise ValueError("Each state needs a positive finite expected yield in the mass fit.")
    shapes = {s: np.array([calibrations[s]["parameters"][n]["value"]
                           for n in SHAPE_NAMES], float) for s in states}
    for state in states:
        if not calibrations[state].get("converged", False):
            raise ValueError(f"{state}: the upstream MC shape fit did not converge.")
        mass_probabilities(edges, shapes[state])
    # Fit five dimensionless or similarly scaled parameters: log(N_1,N_2,N_3),
    # common offset in GeV, log(common width scale).
    offsets = [shapes[s][0] - MASSES[s] for s in states]
    offset0 = float(np.average(offsets, weights=totals))
    initial = np.r_[np.log(totals), np.clip(offset0, -.09, .09), 0.]
    lower = np.r_[np.full(3, -20.), -.1, np.log(.5)]
    upper = np.r_[np.log(np.maximum(10 * totals, 100.)), .1, np.log(2.)]
    active = np.arange(5 if float_common_shape else 3)
    expected = np.sum([means[s] for s in states], axis=0)

    def expand(p):
        full = initial.copy()
        full[active] = p
        return full

    def components(p):
        full = expand(p)
        result = {}
        for j, state in enumerate(states):
            shape = shapes[state].copy()
            shape[0] = MASSES[state] + full[3]
            shape[1] *= np.exp(full[4])
            result[state] = np.exp(full[j]) * mass_probabilities(edges, shape)
        return result

    def model(p):
        return np.sum(list(components(p).values()), axis=0)

    fit = least_squares(lambda p: poisson_residual(expected, model(p)), initial[active],
                        bounds=(lower[active], upper[active]), x_scale="jac",
                        ftol=1e-10, xtol=1e-10, gtol=1e-9, max_nfev=5000)
    full = expand(fit.x)
    fitted = components(fit.x)
    summed = np.sum(list(fitted.values()), axis=0)
    deviance = float(np.sum(poisson_residual(expected, summed)**2))
    parameter_order = ["log_yield_1S", "log_yield_2S", "log_yield_3S",
                       "common_mass_offset_GeV", "log_common_width_scale"]
    on_bound = np.isclose(fit.x, lower[active], atol=1e-5) | np.isclose(fit.x, upper[active], atol=1e-5)
    summary = {
        "converged": bool(fit.success), "optimizer_message": str(fit.message),
        "fit_kind": "Simultaneous three-state fit to expected histogram means",
        "objective": "Expected extended-Poisson deviance (Asimov); not observed counts",
        "background": "none; no continuum or background sample supplied",
        "mass_differences": "Fixed to the nominal Upsilon masses recorded below",
        "nominal_masses_GeV": MASSES,
        "fixed_shape_parameters": "State-specific DSCB tails and relative core widths from independent MC fits",
        "common_shape_floating": bool(float_common_shape),
        "floating_parameter_order": [parameter_order[i] for i in active],
        "common_mass_offset_GeV": float(full[3]),
        "common_width_scale": float(np.exp(full[4])),
        "fit_interval_GeV": [float(edges[0]), float(edges[-1])],
        "poisson_deviance": deviance, "nominal_ndof": int(len(expected) - len(active)),
        "at_parameter_bound": bool(np.any(on_bound)),
        "states": {s: {"expected_selected_yield": float(totals[j]),
                        "fitted_yield": float(np.exp(full[j])),
                        "fitted_mean_GeV": float(MASSES[s] + full[3]),
                        "fitted_sigma_core_GeV": float(shapes[s][1] * np.exp(full[4])),
                        "upstream_shape_mismatch_flag": bool(calibrations[s].get("shape_mismatch_flag", False))}
                   for j, s in enumerate(states)},
        "uncertainties": "No fit-parameter errors quoted; fixed-tail conditional model and known MC shape mismatch require validation before yield inference",
        "deviance_interpretation": "Descriptive expected discrepancy; nominal ndof is not a calibrated goodness-of-fit test at low expected counts",
    }
    if not fit.success:
        raise RuntimeError(f"Combined mass fit failed: {fit.message}")
    return fitted, summary


def _load_inputs(config):
    loaded = {}
    for state in MASSES:
        validated_metadata = sample_metadata(config, state)
        path = resolve_artifacts(config.runs[state]) / "performance" / config.energy
        bins = read_csv(path / "mass_bins.csv", required=("low", "high", "projected_yield"))
        calibration = read_json(path / "mass_fit.json")
        audit = read_json(path / "performance_audit.json")
        eff = read_csv(path / "selection_efficiency.csv", required=("variable", "low", "high", "denominator", "selected", "efficiency", "wilson_low", "wilson_high"))
        source_lumi = calibration.get("statistical_configuration", {}).get("active_luminosity_invfb")
        if source_lumi is None or not np.isfinite(source_lumi) or source_lumi <= 0:
            raise ValueError(f"{path / 'mass_fit.json'}: missing positive active_luminosity_invfb.")
        if calibration.get("state") != state or audit.get("state") != state:
            raise ValueError(f"{path}: state metadata does not match requested {state}.")
        if not calibration.get("constant_selected_weights", False):
            raise ValueError(f"{path}: calibrated unweighted DSCB shapes require constant selected event weights for this projection.")
        truth = audit.get("truth_selection", {})
        if "Q2_GeV2" not in truth or "W_GeV" not in truth:
            raise ValueError(f"{path / 'performance_audit.json'}: truth Q2/W metadata missing.")
        configured_q2 = getattr(config, "q2_max", {}).get(state)
        if configured_q2 is not None and not np.isclose(float(truth["Q2_GeV2"][1]), float(configured_q2)):
            raise ValueError(f"{state}: configured Q2 maximum {configured_q2} differs from the performance audit {truth['Q2_GeV2'][1]}.")
        values = bins[["low", "high", "projected_yield"]].to_numpy(float)
        if not np.all(np.isfinite(values)) or np.any(values[:, 2] < 0):
            raise ValueError(f"{path}: mass bins must contain finite, nonnegative expected yields.")
        if np.any(values[:, 1] <= values[:, 0]) or not np.allclose(values[:-1, 1], values[1:, 0]):
            raise ValueError(f"{path}: mass bins must be ordered and contiguous.")
        if not np.isclose(float(bins.projected_yield.sum()), float(calibration["selected_projected_yield"]), rtol=1e-7, atol=1e-10):
            raise ValueError(f"{path}: mass histogram sum does not match the calibrated selected yield.")
        loaded[state] = dict(path=path, bins=bins, calibration=calibration, audit=audit, efficiency=eff, validated_metadata=validated_metadata,
                             source_lumi=float(source_lumi), scale=float(config.lumi_invfb / source_lumi))
    reference = loaded["1S"]["bins"][["low", "high"]].to_numpy(float)
    for state, data in loaded.items():
        test = data["bins"][["low", "high"]].to_numpy(float)
        if test.shape != reference.shape or not np.allclose(test, reference, atol=1e-10, rtol=0):
            raise ValueError(f"Mass binning differs for {state}; regenerate with common bins before combining.")
    return loaded


def plot_reconstruction(config, *, rebin=2, float_common_shape=True,
                        efficiency_variables=("true_W", "true_pT"), show_mass_intervals=False):
    """Write upsilon_es_reconstruction.pdf/png and its numeric audit products.

    ``rebin`` merges adjacent mass bins exactly before fitting; defaults to
    50 MeV for the framework's 25 MeV bins. Efficiencies retain MC Wilson
    intervals and are never rescaled with luminosity. The mass histogram is
    drawn as a line by default; optional mass intervals remain in the CSV.
    """
    if len(efficiency_variables) != 2:
        raise ValueError("Choose exactly two efficiency variables for the two lower panels.")
    if not isinstance(rebin, (int, np.integer)) or rebin < 1:
        raise ValueError("rebin must be a positive integer.")
    set_style()
    inputs = _load_inputs(config)
    frame = inputs["1S"]["bins"]
    if len(frame) % rebin:
        raise ValueError("rebin must divide the number of input mass bins exactly.")
    fine_edges = np.r_[frame.low.to_numpy(float), float(frame.high.iloc[-1])]
    edges = fine_edges[::rebin]
    means = {s: (v["bins"].projected_yield.to_numpy(float) * v["scale"]).reshape(-1, rebin).sum(axis=1)
             for s, v in inputs.items()}
    fitted, fit_summary = fit_combined_mass(edges, means,
        {s: v["calibration"] for s, v in inputs.items()}, float_common_shape=float_common_shape)
    summed, fit_total = np.sum(list(means.values()), axis=0), np.sum(list(fitted.values()), axis=0)
    coverage = getattr(config, "coverage", .6826894921370859)
    interval = expected_intervals(summed, coverage=coverage)
    # Shared interval helper returns absolute limits, preserving the mean even
    # when it lies outside the interval obtained at the median integer count.
    if isinstance(interval, dict):
        lower, upper = np.asarray(interval["lower"]), np.asarray(interval["upper"])
    else:
        lower, upper = np.asarray(interval[0]), np.asarray(interval[1])
    fig = plt.figure(figsize=(7.2, 8.0))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.65, 1), hspace=.43, wspace=.40,
                          left=.14, right=.98, top=.97, bottom=.105)
    ax = fig.add_subplot(grid[0, :])
    ax.stairs(summed, edges, color="black", lw=1.1, zorder=2, label="_nolegend_")
    centers, halfwidth = (edges[:-1] + edges[1:]) / 2, np.diff(edges) / 2
    if show_mass_intervals:
        positive = summed > 0
        draw_intervals(ax, centers[positive], summed[positive], lower[positive], upper[positive],
                       color="black", fmt="none", elinewidth=.7, capsize=1.3, zorder=2)
    # The thick total remains visible around narrow component dashes and in
    # their gaps, including where a single state almost equals the total.
    # Keep the histogram underneath both curves and preserve the fit arrays.
    total_style = dict(color="#ce242b", lw=2.6, zorder=4)
    component_style = dict(lw=.95, linestyle=(0, (4, 3)), zorder=5)
    width = np.diff(edges)
    if np.allclose(width, width[0]):
        # Display smooth curves as density times the plotted bin width.  The
        # optimization above still used exact integrals over the data bins.
        dense_edges = np.linspace(edges[0], edges[-1], 2001)
        dense_centers = (dense_edges[:-1] + dense_edges[1:]) / 2
        curves = {}
        for state in MASSES:
            shape = np.array([inputs[state]["calibration"]["parameters"][n]["value"] for n in SHAPE_NAMES])
            shape[0] = fit_summary["states"][state]["fitted_mean_GeV"]
            shape[1] = fit_summary["states"][state]["fitted_sigma_core_GeV"]
            curves[state] = (fit_summary["states"][state]["fitted_yield"] * mass_probabilities(dense_edges, shape)
                             * width[0] / np.diff(dense_edges))
        ax.plot(dense_centers, np.sum(list(curves.values()), axis=0), label="Fit",
                **total_style)
        for state in MASSES:
            ax.plot(dense_centers, curves[state], color=COLORS[state], label=state_label(state),
                    **component_style)
    else:
        ax.stairs(fit_total, edges, label="Fit", **total_style)
        for state in MASSES:
            ax.stairs(fitted[state], edges, color=COLORS[state], label=state_label(state),
                      **component_style)
    ax.set_ylabel(f"Expected events / {width[0]:g} GeV" if np.allclose(width, width[0]) else "Expected events / bin")
    ax.set_xlabel(r"$m_{e^{+}e^{-}}$ [GeV]", loc="right")
    ax.set_xlim(edges[0], edges[-1])
    peak = max(float(np.max(upper if show_mass_intervals else summed)), float(fit_total.max()), 1.)
    ax.set_ylim(0, peak * 1.6)
    annotate(ax, config, panel="a")
    ax.legend(loc="upper left", bbox_to_anchor=(.025, .73), fontsize=FONT_SIZES["legend"], frameon=False, handlelength=1.6,
              borderaxespad=.4, labelspacing=.25)
    fit_table_rows = [dict(state=state,
                          fitted_yield=float(fit_summary["states"][state]["fitted_yield"]),
                          fitted_sigma_core_MeV=1000.0 * float(fit_summary["states"][state]["fitted_sigma_core_GeV"]))
                      for state in MASSES]
    fit_table = ax.table(
        cellText=[[state_label(r["state"]), f'{r["fitted_yield"]:.1f}',
                   f'{r["fitted_sigma_core_MeV"]:.1f}'] for r in fit_table_rows],
        colLabels=["State", r"$N_{\mathrm{fit}}$", r"$\sigma_{\mathrm{core}}$ [MeV]"],
        colWidths=[.25, .28, .47], cellLoc="center", colLoc="center",
        bbox=[.52, .61, .45, .29], edges="open")
    fit_table.auto_set_font_size(False)
    fit_table.set_fontsize(FONT_SIZES["annotation"])
    for (row, col), cell in fit_table.get_celld().items():
        cell.PAD = .025
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.visible_edges = "B"
            cell.set_edgecolor("#777777")
            cell.set_linewidth(.6)
        elif col == 0:
            cell.set_text_props(color=COLORS[fit_table_rows[row - 1]["state"]])

    labels = {"true_W": r"True $W_{\gamma p}$ [GeV]", "true_y": r"True $y_{\Upsilon}$",
              "true_pT": r"True $p_{T,\Upsilon}$ [GeV]"}
    q2 = {s: float(v["audit"]["truth_selection"]["Q2_GeV2"][1]) for s, v in inputs.items()}
    categories = getattr(config, "categories", {})
    for index, variable in enumerate(efficiency_variables):
        axis = fig.add_subplot(grid[1, index])
        maximum_efficiency = 0.
        for state, data in inputs.items():
            e = data["efficiency"].loc[data["efficiency"].variable == variable]
            if e.empty:
                raise ValueError(f"{data['path']}: efficiency variable {variable!r} unavailable.")
            good = (e.denominator.to_numpy(float) > 0) & np.isfinite(e.efficiency.to_numpy(float))
            e = e.loc[good]
            values = e.efficiency.to_numpy(float)
            lower_eff, upper_eff = e.wilson_low.to_numpy(float), e.wilson_high.to_numpy(float)
            if np.any((values < 0) | (values > 1) | (lower_eff > values) | (upper_eff < values)):
                raise ValueError(f"{data['path']}: invalid efficiency or Wilson interval.")
            if len(upper_eff):
                maximum_efficiency = max(maximum_efficiency, float(upper_eff.max()))
            axis.errorbar((e.low.to_numpy(float) + e.high.to_numpy(float)) / 2, 100 * values,
                          xerr=(e.high.to_numpy(float) - e.low.to_numpy(float)) / 2,
                          yerr=100 * np.array([values - lower_eff, upper_eff - values]),
                          fmt=MARKERS[state], ms=4.0, color=COLORS[state], lw=.9,
                          capsize=1.8, label=state_label(state))
        axis.set_xlabel(labels.get(variable, variable))
        axis.set_ylim(0, min(100, max(50, 10 * np.ceil(16 * maximum_efficiency))))
        axis.text(.04, .96, "(b)" if index == 0 else "(c)", transform=axis.transAxes,
                  va="top", fontsize=FONT_SIZES["panel"])
        axis.set_ylabel("Efficiency [%]")
        if index == 1:
            axis.legend(frameon=False, loc="upper right", fontsize=FONT_SIZES["legend"], labelspacing=.2)
        from matplotlib.ticker import MaxNLocator
        axis.xaxis.set_major_locator(MaxNLocator(nbins=4))
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4))

    # Selection and model details remain in the audit and notebook notes.
    out = Path(config.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    warnings_list = []
    if len(set(q2.values())) > 1 or len(set(categories.get(s, "unspecified") for s in MASSES)) > 1:
        warnings_list.append("The states have different Q2 and/or electron-ID settings. Efficiency differences do not isolate the state-mass dependence.")
    if any(v["calibration"].get("shape_mismatch_flag", False) for v in inputs.values()):
        warnings_list.append("At least one upstream DSCB calibration has significant shape mismatch. The combined fit is a conditional descriptive baseline, not a validated yield-extraction model.")
    warnings_list.append("Efficiencies are conditional on eligible stored tuple rows; a complete generated-event denominator is not independently established here.")
    for item in warnings_list:
        warnings.warn(item, UserWarning, stacklevel=2)
    fit_summary.update({
        "target_luminosity_invfb": float(config.lumi_invfb), "energy": config.energy,
        "input_files": {s: [str(v["path"] / name) for name in
                            ("mass_bins.csv", "mass_fit.json", "selection_efficiency.csv", "performance_audit.json")]
                        for s, v in inputs.items()},
        "source_luminosities_invfb": {s: v["source_lumi"] for s, v in inputs.items()},
        "yield_scale_factors": {s: v["scale"] for s, v in inputs.items()},
        "truth_selections": {s: v["audit"]["truth_selection"] for s, v in inputs.items()},
        "reconstruction_selections": {s: v["audit"]["reco_selection"] for s, v in inputs.items()},
        "electron_id_labels": categories, "mass_rebin_factor": int(rebin),
        "displayed_fit_table": fit_table_rows,
        "displayed_resolution_definition": "Fitted DSCB Gaussian-core sigma, converted from GeV to MeV; not the full RMS. Relative state widths retain the fit constraints.",
        "displayed_fit_table_uncertainties": "No parameter uncertainties supplied by this conditional fit; none invented for display.",
        "publication_footer_displayed": False,
        "input_validation": {s: v["validated_metadata"] for s, v in inputs.items()},
        "displayed_mass_intervals": "Ensemble-median Garwood limits around median integer observations; expected means unchanged" if show_mass_intervals else "not drawn",
        "mass_interval_coverage": float(coverage),
        "zero_mean_mass_bins": "Shown in histogram; interval bars omitted because empty finite-MC bins do not establish a zero physical expectation",
        "efficiency_intervals": "Original finite-MC Wilson intervals; no luminosity rescaling",
        "efficiency_variables": list(efficiency_variables), "warnings": warnings_list,
    })
    table = pd.DataFrame({"low": edges[:-1], "high": edges[1:], "expected_signal": summed,
                          "interval_low": lower, "interval_high": upper, "combined_fit": fit_total,
                          **{f"expected_{s}": means[s] for s in MASSES},
                          **{f"fit_{s}": fitted[s] for s in MASSES}})
    table.to_csv(out / "upsilon_es_reconstruction_mass.csv", index=False)
    write_json(out / "upsilon_es_reconstruction_audit.json", fit_summary)
    save_figure(fig, config, "upsilon_es_reconstruction")
    return fig, fit_summary
