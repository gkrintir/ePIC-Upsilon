#!/usr/bin/env python3
"""Render paper Figure 28 and the corresponding note panels from saved results.

Run from a checkout with the plotting dependencies installed::

    python tools/render_paper_figure28.py --results-dir EXTRACTED_RESULTS \
        --output-dir NEW_OUTPUT_DIRECTORY

The input is the checksum-manifested full-sample ECAL-rematched results archive.
All fits, profile endpoints and sampling bands are read from that archive. The
renderer validates the saved model against its binned predictions before drawing
the smooth curves. It performs no fit, luminosity rescaling or ROOT-file read.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "plotting"))
sys.path.insert(0, str(REPOSITORY / "tools"))

from export_rematched_atlas import FRAMEWORK_COMMIT, VerifiedArchive
from upsilon_esr.common import COLORS, MARKERS, draw_intervals, set_style, state_label
from upsilon_esr.differential import exponential_bin_average
from upsilon_esr.reconstruction import MASSES, SHAPE_NAMES, mass_probabilities


STATES = tuple(MASSES)
CAMPAIGN = "26.07.1 Campaign"
MASS_PREFIX = "paper_figures/mass_efficiency/"
DIFFERENTIAL_PREFIX = "paper_figures/spectra/differential/"
MASS_TABLE_FONT_SIZE = 18.5
REFERENCE_COLOR = "#62676c"


def numeric_column(rows, key):
    values = np.array([float(row[key]) for row in rows])
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite values in {key}")
    return values


def assert_close(actual, expected, name):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=1e-8, atol=1e-10):
        raise ValueError(f"Saved values disagree for {name}")
    return float(np.max(np.abs(actual - expected)))


def load_mass(archive):
    rows = archive.csv(MASS_PREFIX + "mass_reference/upsilon_es_reconstruction_mass.csv")
    fit = archive.json(MASS_PREFIX + "mass_reference/upsilon_es_reconstruction_audit.json")
    intervals = archive.json(MASS_PREFIX + "mass_fit_intervals.json")
    if (fit["energy"] != "9x275" or fit["target_luminosity_invfb"] != 2.5
            or intervals["target_luminosity_invfb"] != 2.5 or not fit["converged"]):
        raise ValueError("Figure 28 requires the converged 9x275, 2.5 fb^-1 mass fit")
    low, high = numeric_column(rows, "low"), numeric_column(rows, "high")
    if len(rows) == 0 or np.any(high <= low):
        raise ValueError("Invalid mass-bin edges")
    assert_close(low[1:], high[:-1], "contiguous mass bins")
    assert_close(high - low, np.full(len(rows), 0.05), "50 MeV mass bins")
    edges = np.r_[low, high[-1]]
    expected = numeric_column(rows, "expected_signal")
    assert_close(expected, sum(numeric_column(rows, "expected_" + s) for s in STATES),
                 "summed expected mass spectrum")
    dense_edges = np.linspace(edges[0], edges[-1], 2001)
    curves, bin_predictions, differences = {}, {}, {}
    for state in STATES:
        calibration_path, = [p for p in fit["input_files"][state] if p.endswith("/mass_fit.json")]
        run_name = calibration_path.split("/runs/", 1)[1].split("/", 1)[0]
        calibration = archive.json(
            f"framework/full_state_runs/{run_name}/artifacts/performance/9x275/mass_fit.json")
        if calibration["state"] != state or not calibration["converged"]:
            raise ValueError(f"Invalid saved mass calibration for {state}")
        shape = np.array([calibration["parameters"][name]["value"] for name in SHAPE_NAMES])
        saved = fit["states"][state]
        shape[0], shape[1] = saved["fitted_mean_GeV"], saved["fitted_sigma_core_GeV"]
        count = saved["fitted_yield"]
        predictions = count * mass_probabilities(edges, shape)
        differences[state] = assert_close(predictions, numeric_column(rows, "fit_" + state),
                                          state + " binned DSCB")
        bin_predictions[state] = predictions
        curves[state] = count * mass_probabilities(dense_edges, shape) * 0.05 / np.diff(dense_edges)
        for quantity, value in (("fitted_yield", count), ("fitted_sigma_core_MeV", 1000 * shape[1])):
            entry = intervals["states"][state][quantity]
            if entry["status"] != "TWO_SIDED":
                raise ValueError(f"Missing two-sided mass interval for {state}/{quantity}")
            assert_close(entry["value"], value, state + " " + quantity)
            assert_close([entry["value"] - entry["lower"], entry["upper"] - entry["value"]],
                         [entry["error_low"], entry["error_high"]], state + " profile endpoints")
    assert_close(sum(bin_predictions.values()), numeric_column(rows, "combined_fit"), "combined mass fit")
    return dict(edges=edges, expected=expected, dense_x=(dense_edges[:-1] + dense_edges[1:]) / 2,
                curves=curves, intervals=intervals, fit=fit, bin_model_max_abs_difference=differences)


def load_differential(archive):
    rows = archive.csv(DIFFERENTIAL_PREFIX + "upsilon_es_differential_points.csv")
    audit = archive.json(DIFFERENTIAL_PREFIX + "upsilon_es_differential_audit.json")
    slope_rows = archive.csv(DIFFERENTIAL_PREFIX + "upsilon_es_slopes.csv")
    if audit["energy"] != "9x275" or audit["target_lumi_invfb"] != 2.5:
        raise ValueError("Figure 28 requires the 9x275, 2.5 fb^-1 differential projection")
    points = {variable: {state: [] for state in STATES} for variable in ("pT2", "y")}
    for row in rows:
        if row["displayed"] not in ("True", "False") or row["valid"] not in ("True", "False"):
            raise ValueError("Saved spectrum validity flags must be explicit booleans")
        if row["displayed"] == "False":
            continue
        if row["valid"] != "True" or row["proton_selection_mode"] != "required":
            raise ValueError("A displayed spectrum bin must have valid proton-tagged input")
        values = {name: float(row[name]) for name in ("x_low", "x_high", "mean", "lower", "upper")}
        if (not np.isfinite(list(values.values())).all() or values["x_high"] <= values["x_low"]
                or values["mean"] < 0 or values["lower"] < 0 or values["upper"] < values["lower"]):
            raise ValueError("Invalid saved differential point")
        values.update(state=row["state"], status=row["overall_status"])
        points[row["variable"]][row["state"]].append(values)
    for case in audit["case_availability"]:
        actual = points[case["variable"]][case["state"]]
        if len(actual) != case["displayed_bins"]:
            raise ValueError("Displayed spectrum count disagrees with its availability audit")
        if any(row["status"] != case["overall_status"] for row in actual):
            raise ValueError("Spectrum point and availability statuses disagree")
    slopes, differences = {}, {}
    for row in slope_rows:
        state = row["state"]
        if not points["pT2"][state]:
            continue
        values = {key: float(row[key]) for key in (
            "b_slope_GeVminus2", "b_interval_low", "b_interval_high", "b_slope_err_low",
            "b_slope_err_high", "amplitude_pb_per_GeV2")}
        if not np.isfinite(list(values.values())).all() or min(values.values()) < 0:
            raise ValueError("Invalid saved slope parameters")
        values["status"] = row["overall_status"]
        saved = audit["slope_fits"][state]
        if values["status"] != saved["overall_status"]:
            raise ValueError(f"Slope validation status disagrees for {state}")
        for key in values:
            if key != "status":
                assert_close(values[key], saved[key], state + " " + key)
        assert_close([values["b_slope_GeVminus2"] - values["b_interval_low"],
                      values["b_interval_high"] - values["b_slope_GeVminus2"]],
                     [values["b_slope_err_low"], values["b_slope_err_high"]], state + " slope interval")
        bins = archive.csv(DIFFERENTIAL_PREFIX + f"upsilon_es_{state.lower()}_pt2_fit_bins.csv")
        fitted = exponential_bin_average(numeric_column(bins, "x_low"), numeric_column(bins, "x_high"),
                                        values["b_slope_GeVminus2"],
                                        amplitude=values["amplitude_pb_per_GeV2"])
        differences[state] = assert_close(fitted, numeric_column(bins, "fit_bin_average_pb_per_GeV2"),
                                          state + " binned exponential")
        values["fit_bins"] = bins
        slopes[state] = values
    if set(slopes) != {state for state in STATES if points["pT2"][state]}:
        raise ValueError("Every displayed pT2 state requires its saved slope")
    references = audit["slope_references"]
    if [entry["experiment"] for entry in references] != ["ZEUS", "CMS"]:
        raise ValueError("Expected the saved ZEUS and CMS slope references")
    for entry in references:
        numeric_column([entry], "slope")
        for key in ("stat_low", "stat_high", "syst_low", "syst_high"):
            if numeric_column([entry], key)[0] < 0:
                raise ValueError("Invalid reference uncertainty")
    return dict(points=points, slopes=slopes, references=references, audit=audit,
                bin_model_max_abs_difference=differences)


def panel_header(ax, panel):
    texts = [
        ax.text(0, 1.285, "ePIC Performance", transform=ax.transAxes, fontsize=17.5,
                weight="bold", ha="left", va="bottom"),
        ax.text(0, 1.19, CAMPAIGN, transform=ax.transAxes, fontsize=14,
                ha="left", va="bottom"),
        ax.text(0, 1.085, r"$e+p$, $9\times275$ GeV, $\mathcal{L}=2.5$ fb$^{-1}$",
                transform=ax.transAxes, fontsize=14, ha="left", va="bottom"),
        ax.text(1, 1.285, f"({panel})", transform=ax.transAxes, fontsize=16,
                ha="right", va="bottom"),
    ]
    return texts


def interval_text(entry):
    return (rf"${entry['value']:.1f}^{{+{entry['error_high']:.1f}}}"
            rf"_{{-{entry['error_low']:.1f}}}$")


def draw_mass(ax, mass, panel="a"):
    ax.stairs(mass["expected"], mass["edges"], color="black", lw=1.05, zorder=2)
    total = sum(mass["curves"].values())
    ax.plot(mass["dense_x"], total, color="#ce242b", lw=2.3, label="Fit", zorder=4)
    for state in STATES:
        ax.plot(mass["dense_x"], mass["curves"][state], color=COLORS[state],
                ls=(0, (4, 3)), lw=1.0, label=state_label(state), zorder=5)
    ax.set_xlim(mass["edges"][0], mass["edges"][-1])
    ax.set_ylim(0, max(total.max(), mass["expected"].max()) * 1.75)
    ax.set_xlabel(r"$m_{e^{+}e^{-}}$ [GeV]", loc="right")
    ax.set_ylabel("Expected events\n/ 0.05 GeV")
    header = panel_header(ax, panel)
    legend = ax.legend(loc="upper left", bbox_to_anchor=(.01, .975), frameon=False,
                       fontsize=14.5, handlelength=1.55, labelspacing=.24)
    body = [[state_label(s), interval_text(mass["intervals"]["states"][s]["fitted_yield"]),
             interval_text(mass["intervals"]["states"][s]["fitted_sigma_core_MeV"])]
            for s in STATES]
    table = ax.table(cellText=body, colLabels=["State", r"$N_{\mathrm{fit}}$",
                                             r"$\sigma_{\mathrm{core}}$ [MeV]"],
                     cellLoc="center", colLoc="center", colWidths=[.23, .33, .44],
                     bbox=[.41, .43, .58, .53], edges="open")
    table.auto_set_font_size(False)
    table.set_fontsize(MASS_TABLE_FONT_SIZE)
    for (row, column), cell in table.get_celld().items():
        cell.PAD = .02
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.visible_edges = "B"
            cell.set_edgecolor("#777777")
            cell.set_linewidth(.6)
        elif column == 0:
            cell.set_text_props(color=COLORS[STATES[row - 1]])
    return dict(header=header, legend=legend, table=table)


def draw_slope_inset(ax, differential):
    states = [state for state in STATES if state in differential["slopes"]]
    inset = ax.inset_axes([.56, .47, .42, .46])
    inset.set_facecolor("white")
    for position, state in enumerate(states):
        fit = differential["slopes"][state]
        inset.errorbar(position, fit["b_slope_GeVminus2"],
                       yerr=[[fit["b_slope_err_low"]], [fit["b_slope_err_high"]]],
                       fmt=MARKERS[state], color=COLORS[state], mfc=COLORS[state],
                       ms=4.6, capsize=2, elinewidth=1.0)
    labels = [state_label(state) for state in states]
    for position, reference in enumerate(differential["references"], len(states)):
        low = np.hypot(reference["stat_low"], reference["syst_low"])
        high = np.hypot(reference["stat_high"], reference["syst_high"])
        inset.errorbar(position, reference["slope"], yerr=[[low], [high]],
                       fmt=MARKERS["1S"], color=REFERENCE_COLOR, mfc=REFERENCE_COLOR,
                       ms=4.6, capsize=2, elinewidth=1.0)
        labels.append(reference["experiment"])
    inset.set_xticks(range(len(labels)), labels, fontsize=10.5)
    for label, state in zip(inset.get_xticklabels(), states):
        label.set_color(COLORS[state])
    inset.set_ylabel(r"$b$ [GeV$^{-2}$]", fontsize=12, labelpad=3)
    inset.set_xlim(-.55, len(labels) - .45)
    inset.set_ylim(0, max([8.6] + [fit["b_interval_high"] * 1.25
                                 for fit in differential["slopes"].values()]))
    inset.tick_params(axis="both", which="both", direction="in", top=True, right=True,
                      labelsize=10.5)
    return inset


def draw_points(ax, points):
    for state, rows in points.items():
        if not rows:
            continue
        low, high = numeric_column(rows, "x_low"), numeric_column(rows, "x_high")
        draw_intervals(ax, (low + high) / 2, numeric_column(rows, "mean"),
                       numeric_column(rows, "lower"), numeric_column(rows, "upper"),
                       xerr=(high - low) / 2, fmt=MARKERS[state], color=COLORS[state],
                       label=state_label(state), ms=4.5, capsize=2, elinewidth=1.0, lw=1.0)


def draw_pt2(ax, differential, panel):
    points = differential["points"]["pT2"]
    draw_points(ax, points)
    for state, fit in differential["slopes"].items():
        x = np.linspace(0, 1, 240)
        ax.plot(x, fit["amplitude_pb_per_GeV2"] * np.exp(-fit["b_slope_GeVminus2"] * x),
                color=COLORS[state], ls="--", lw=1.2, alpha=.9, zorder=1)
        for row in fit["fit_bins"]:
            if row["included_in_fit"] == "True":
                ax.plot([float(row["x_low"]), float(row["x_high"])],
                        [float(row["fit_bin_average_pb_per_GeV2"])] * 2,
                        color=COLORS[state], lw=1.6, alpha=.9, zorder=1)
    ax.set_xlabel(r"$p_{T,\Upsilon}^{2}$ [GeV$^2$]")
    ax.set_ylabel(r"$d\sigma_{\gamma p}/dp_{T,\Upsilon}^{2}$ [pb/GeV$^2$]")
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    all_rows = [row for rows in points.values() for row in rows]
    positive_lower = [row["lower"] for row in all_rows if row["lower"] > 0]
    floor = min(min(row["mean"] for row in all_rows) / 5, min(positive_lower) / 1.5)
    left_upper = max(row["upper"] for row in all_rows if row["x_low"] <= .50)
    right_upper = max(row["upper"] for row in all_rows if row["x_high"] >= .60)
    ceiling = max(floor * (1.15 * left_upper / floor) ** (1 / .60),
                  floor * (1.25 * right_upper / floor) ** (1 / .39))
    ax.set_ylim(floor, ceiling)
    for state, rows in points.items():
        zeros = [row for row in rows if row["lower"] <= 0]
        ax.plot([(row["x_low"] + row["x_high"]) / 2 for row in zeros],
                [floor * 1.08] * len(zeros), "v", color=COLORS[state], ms=3, clip_on=False)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([], [], color="#555555", lw=1.4, ls="--"))
    labels.append(r"Fit: $0<p_T^2<1$ GeV$^2$")
    legend = ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(.01, .985),
                       frameon=False, fontsize=13.5, labelspacing=.25, handlelength=1.4)
    header = panel_header(ax, panel)
    inset = draw_slope_inset(ax, differential)
    return dict(header=header, legend=legend, inset=inset)


def draw_rapidity(ax, differential):
    points = differential["points"]["y"]
    draw_points(ax, points)
    all_rows = [row for rows in points.values() for row in rows]
    ax.set_xlabel(r"$y_{\Upsilon}$ (laboratory)")
    ax.set_ylabel(r"$d\sigma_{\gamma p}/dy_{\Upsilon}$ [pb]")
    ax.set_xlim(min(row["x_low"] for row in all_rows), max(row["x_high"] for row in all_rows))
    ax.set_ylim(0, max(max(row["upper"] for row in all_rows) * 1.35, 1.0))
    ax.text(.025, .95, "(b)", transform=ax.transAxes, va="top", fontsize=15)
    ax.legend(loc="upper right", frameon=False, fontsize=14, labelspacing=.25)


def check_layout(figure, mass_artists=None, pt2_artists=None):
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    checks = {}
    if mass_artists:
        table = mass_artists["table"]
        legend_box = mass_artists["legend"].get_window_extent(renderer)
        table_box = table.get_window_extent(renderer)
        if table_box.overlaps(legend_box):
            raise ValueError("Mass statistics table overlaps the component legend")
        for text in mass_artists["header"]:
            if text.get_window_extent(renderer).overlaps(table_box):
                raise ValueError("Mass statistics table overlaps the panel header")
        for position, cell in table.get_celld().items():
            text_box = cell.get_text().get_window_extent(renderer)
            box = cell.get_window_extent(renderer)
            if text_box.width > box.width + 1 or text_box.height > box.height + 1:
                raise ValueError(f"Mass statistics text exceeds table cell {position}")
        checks["mass_table_clear_of_header_and_legend"] = True
        checks["mass_table_text_fits_cells"] = True
    if pt2_artists:
        inset = pt2_artists["inset"]
        inset_box = inset.get_tightbbox(renderer)
        if inset_box.overlaps(pt2_artists["legend"].get_window_extent(renderer)):
            raise ValueError("Slope inset overlaps the pT2 legend")
        checks["slope_inset_clear_of_legend"] = True
    return checks


def make_figures(mass, differential):
    set_style()
    paper, axes = plt.subplots(2, 1, figsize=(7.2, 10.2))
    paper.subplots_adjust(left=.145, right=.98, bottom=.07, top=.86, hspace=.73)
    mass_artists = draw_mass(axes[0], mass)
    pt2_artists = draw_pt2(axes[1], differential, "b")
    paper_check = check_layout(paper, mass_artists, pt2_artists)

    mass_figure, mass_axis = plt.subplots(figsize=(7.2, 4.9))
    mass_figure.subplots_adjust(left=.145, right=.98, bottom=.16, top=.75)
    mass_check = check_layout(mass_figure, draw_mass(mass_axis, mass))

    note, note_axes = plt.subplots(2, 1, figsize=(7.0, 8.0))
    note.subplots_adjust(left=.16, right=.97, top=.81, bottom=.10, hspace=.40)
    note_artists = draw_pt2(note_axes[0], differential, "a")
    draw_rapidity(note_axes[1], differential)
    note_check = check_layout(note, pt2_artists=note_artists)
    return {"upsilon_es_figure28": (paper, paper_check),
            "note_mass_panel": (mass_figure, mass_check),
            "note_differential": (note, note_check)}


def render(results_dir, output_dir):
    results_dir, output_dir = Path(results_dir).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"Choose a new output directory: {output_dir}")
    archive = VerifiedArchive(results_dir)
    if FRAMEWORK_COMMIT not in archive.read("README.md").decode("utf-8"):
        raise ValueError("Archive analysis-commit provenance differs from this rematched dataset")
    mass, differential = load_mass(archive), load_differential(archive)
    figures = make_figures(mass, differential)
    audit = {
        "figure": "paper_figure28", "energy": "9x275", "luminosity_invfb": 2.5,
        "framework_commit": FRAMEWORK_COMMIT,
        "campaign": CAMPAIGN, "state_colors": COLORS, "state_markers": MARKERS,
        "reference_marker": MARKERS["1S"], "reference_color": REFERENCE_COLOR,
        "reference_marker_fill": "same as marker color", "inset_epic_label": None,
        "mass_table_font_size_pt": MASS_TABLE_FONT_SIZE,
        "mass_parameter_intervals": mass["intervals"],
        "spectrum_points": differential["points"],
        "slope_fits": differential["slopes"], "slope_references": differential["references"],
        "case_availability": differential["audit"]["case_availability"],
        "interval_interpretation": differential["audit"]["uncertainty"],
        "caption_notes": differential["audit"]["caption_notes"],
        "mass_binned_model_max_abs_difference": mass["bin_model_max_abs_difference"],
        "pt2_binned_model_max_abs_difference": differential["bin_model_max_abs_difference"],
        "archive_manifest_sha256": archive.manifest_sha256,
        "verified_sources": list(archive.checked.values()),
        "layout_checks": {name: checks for name, (_, checks) in figures.items()},
        "refit_performed": False, "luminosity_rescaling_applied": False,
        "outputs": [],
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".figure28-", dir=output_dir.parent) as temporary:
        staging = Path(temporary)
        for name, (figure, _) in figures.items():
            for extension in ("pdf", "png"):
                path = staging / f"{name}.{extension}"
                metadata = {"CreationDate": None, "ModDate": None} if extension == "pdf" else None
                figure.savefig(path, bbox_inches="tight", pad_inches=.08, dpi=220, metadata=metadata)
                audit["outputs"].append({"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            plt.close(figure)
        (staging / "figure28_audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
        if output_dir.exists():
            raise FileExistsError(f"Output directory appeared during rendering: {output_dir}")
        staging.rename(output_dir)
    return output_dir / "figure28_audit.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="Extracted corrected-results archive containing SHA256SUMS.txt")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory for figures and audit")
    args = parser.parse_args()
    audit_path = render(args.results_dir, args.output_dir)
    print(f"Rendered Figure 28 and matching note panels. Audit: {audit_path}")


if __name__ == "__main__":
    main()
