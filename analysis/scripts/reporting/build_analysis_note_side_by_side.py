#!/usr/bin/env python3
"""Build the corrected analysis-note figure comparison PDF and index CSV.

The left column is the exact raster included by the LaTeX analysis note.  The
right column is a fresh rerun (or, where explicitly labelled, a rerender from
an encoded CSV).  Figures 1, 2, 19, and 28 are intentionally outside the
comparison scope requested by the analysis team.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import textwrap
from pathlib import Path

from PIL import Image
from reportlab.lib.colors import Color, HexColor, white
from reportlab.lib.pagesizes import landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen import canvas


PAGE = landscape((540, 960))
NAVY = HexColor("#102A43")
BLUE = HexColor("#1E66F5")
TEAL = HexColor("#0F8B8D")
GREEN = HexColor("#218739")
AMBER = HexColor("#C47F00")
RED = HexColor("#B33A3A")
INK = HexColor("#172B4D")
MUTED = HexColor("#52667A")
LIGHT = HexColor("#F4F7FA")
BORDER = HexColor("#CBD5E1")


FIGURES = [
    dict(n=3, title="Generator validation using W and |t|", report=["Upsilon_TrueLevel_Overlay_Dashboard.png"], rerun=["fig03_truth/Upsilon_TrueLevel_Overlay_Dashboard.png"], old="MassFit-Copy4.ipynb C7", new="notebooks/03_event_selection_truth_and_bt_studies.ipynb C7", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix2.root", status="Reproduced", tone="green", note="Fresh tuple rerun."),
    dict(n=4, title="Electron-identification variables", report=["Upsilon_EleID_Variables.png"], rerun=["fig04_eleid_components/h_EleID_EpRatio.png", "fig04_eleid_components/h_EleID_IsoFraction.png"], old="Exact composite producer not retained; MassFit-Copy8.ipynb C0 has components", new="scripts/plotting/compare_report_campaign_histograms.py", inputs="merged_elecs_test.root + merged_elecs_test_.root (stored campaign TH1s)", status="Component-level", tone="amber", note="The two stored histogram components rerun; the exact two-energy composite layout is not recoverable."),
    dict(n=5, title="Central-track and tagger-track correlation", report=["central_tagger_corr_10x100_1S.png", "central_tagger_corr_18x275_1S.png"], rerun=["fig05_corr/central_tagger_corr_10x100_1S.png", "fig05_corr/central_tagger_corr_18x275_1S.png"], old="plot_central_tagger_correlation.py", new="scripts/validation/plot_track_multiplicity_correlation.py", inputs="merged_elecs_test_10t100.root + merged_elecs_test_tuples_fix2.root", status="Data exact", tone="green", note="Map and summary CSV files reproduce byte-for-byte."),
    dict(n=6, title="Preliminary event-selection cut flow", report=["Upsilon_EventSelection_Cutflow_with_Loss.png"], rerun=["fig06_cutflow/Upsilon_EventSelection_Cutflow_with_Loss.png"], old="MassFit-Copy4.ipynb C8", new="notebooks/03_event_selection_truth_and_bt_studies.ipynb C8", inputs="merged_elecs_test_10t100.root + merged_elecs_test_tuples.root", status="Reproduced", tone="green", note="Fresh notebook rerun."),
    dict(n=7, title="Inclusive Upsilon mass fit at 18x275 GeV", report=["mass_fit_summary.png"], rerun=["fig07_fit_high/Fit_Inclusive_Upsilon_DSCB_Nice.png"], old="MassFit.ipynb inclusive-fit cell", new="notebooks/01_mass_fits_inclusive_and_rapidity.ipynb C1 (parameterized high-energy rerun)", inputs="merged_elecs_test_tuples_fix2.root", status="Numerically close", tone="amber", note="Yields agree within 0-4 events and widths within about 0.2 MeV; minimizer reached its call limit."),
    dict(n=8, title="Inclusive mass fits in rapidity bins", report=["mass_fit_summary-Y.png"], rerun=["fig08_fit_rapidity/Fit_Inclusive_Upsilon_RapidityBins_raw_visible_decay.png"], old="MassFit.ipynb C2", new="notebooks/01_mass_fits_inclusive_and_rapidity.ipynb C2", inputs="merged_elecs_test_tuples_fix2.root", status="Source reproduced", tone="amber", note="Fresh rerun produced the panels, but the fit covariance is not positive definite."),
    dict(n=9, title="Inclusive Upsilon mass fit at 10x100 GeV", report=["mass_fit_summary_10100.png"], rerun=["fig09_fit_low/Fit_Inclusive_Upsilon_DSCB_Nice.png"], old="MassFit.ipynb C1", new="notebooks/01_mass_fits_inclusive_and_rapidity.ipynb C1", inputs="merged_elecs_test_10t100_10x100.root", status="Source reproduced", tone="amber", note="Fresh rerun produced the plot, but the minimizer reached its call limit."),
    dict(n=10, title="Projected scaled-yield distributions at 18x275 GeV", report=["scaled_yields_summary.png"], rerun=["fig10_scaled_yields/Upsilon_ThreeState_Overlay_Dashboard_target_visible_xsec_visibleEE.png"], old="MassFit-Copy4.ipynb dashboard cells", new="notebooks/03_event_selection_truth_and_bt_studies.ipynb C4 (explicit input override)", inputs="merged_elecs_test_tuples_fix2.root; target-visible normalization encoded in notebook", status="Reproduced", tone="green", note="Fresh tuple rerun with the report normalization."),
    dict(n=11, title="Projected 1S distributions at two beam energies", report=["Upsilon_1S_TwoEnergy_Overlay_Dashboard_target_visible_xsec_visibleEE.png"], rerun=["fig11_two_energy/Upsilon_1S_TwoEnergy_Overlay_Dashboard_target_visible_xsec_visibleEE.png"], old="MassFit-Copy4.ipynb two-energy dashboard", new="notebooks/03_event_selection_truth_and_bt_studies.ipynb C2", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix2.root", status="Reproduced", tone="green", note="Fresh tuple rerun."),
    dict(n=12, title="Reconstructed |t| estimator comparison", report=["Upsilon_1S_TMethods_TwoEnergy_SideBySide_target_visible_xsec_visibleEE.png"], rerun=["fig12_t_methods/Upsilon_1S_TMethods_TwoEnergy_SideBySide_target_visible_xsec_visibleEE.png"], old="Upsilon_1S_TMethods_TwoEnergy_Overlay_visibleEE.py", new="notebooks/03_event_selection_truth_and_bt_studies.ipynb C3", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix2.root", status="Reproduced", tone="green", note="Fresh tuple rerun; summary CSV also regenerated."),
    dict(n=13, title="Acceptance and reconstruction efficiency summary", report=["Upsilon_1S_Kinematic_Efficiencies_TwoEnergyOverlay.png"], rerun=["fig13_efficiency/Upsilon_1S_Kinematic_Efficiencies_TwoEnergyOverlay.png"], old="MassFit-Copy2.ipynb C3", new="notebooks/02_efficiency_and_two_energy_studies.ipynb C3", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix2.root", status="Reproduced", tone="green", note="Fresh tuple rerun."),
    dict(n=14, title="Upsilon(nS) rapidity coverage", report=["diff_gamma_p_family_18x275_W37to114_bSlope_pt2_y.png"], rerun=["fig14_17_18/diff_gamma_p_family_18x275_W37to114_bSlope_pt2_y.png"], crop="right", old="extract_differential_gamma_p_swan_v4.py", new="scripts/cross_sections/extract_differential_gamma_p_bin_by_bin.py", inputs="merged_elecs_test_tuples_fix2.root", status="Data exact", tone="green", note="All plotted numeric columns reproduce to floating-point precision."),
    dict(n=15, title="Projected exclusive Upsilon(nS) photoproduction cross sections", report=["upsilon_selected_overlay_1S_family.png"], rerun=["fig15_20_23/upsilon_selected_overlay_1S_family.png"], old="extract_sigma_*_swan_v2.py + plot_selected_upsilon_comparisons_v2.py", new="scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py + scripts/plotting/build_gamma_p_validation_figure_set.py", inputs="merged_elecs_test_10t100.root + merged_elecs_test_tuples_fix.root", status="Central values exact", tone="green", note="Normalization CSVs are byte-exact; promoted code adds explicit uncertainty components."),
    dict(n=16, title="|t| spectrum and impact-parameter profile", report=["pt2_twoE_ZEUS_bT.png"], rerun=["fig16/pt2_twoE_ZEUS_bT.png"], crop="top", old="overlay_pt2_twoE_with_ZEUS_and_bT.py", new="scripts/plotting/plot_t_spectra_and_impact_parameter.py", inputs="Fresh low/high t-unfold summary CSVs + digitized ZEUS CSV", status="Fresh unfolding", tone="green", note="The ROOT-derived summary CSVs are explicit inputs; no DVMP_JPsi_Analysis.C dependency."),
    dict(n=17, title="Differential Upsilon(nS) pT-squared spectra", report=["diff_gamma_p_family_18x275_W37to114_bSlope_pt2_y.png"], rerun=["fig14_17_18/diff_gamma_p_family_18x275_W37to114_bSlope_pt2_y.png"], crop="left", old="extract_differential_gamma_p_swan_v4.py", new="scripts/cross_sections/extract_differential_gamma_p_bin_by_bin.py", inputs="merged_elecs_test_tuples_fix2.root", status="Data exact", tone="green", note="All plotted numeric columns reproduce to floating-point precision."),
    dict(n=18, title="Extracted Upsilon(nS) b-slopes", report=["b_slope_fig4_overlay_3states_18x275.png"], rerun=["fig14_17_18/b_slope_fig4_overlay_3states_18x275.png"], old="overlay_b_slopes_three_states_on_hera_fig4.py", new="scripts/plotting/plot_b_slopes_family_hera_overlay.py", inputs="Fresh b-slope summary CSV from merged_elecs_test_tuples_fix2.root", status="Data exact", tone="green", note="Central values and fit diagnostics reproduce to floating-point precision."),
    dict(n=20, title="Photon-flux validation", report=["upsilon_selected_flux_validation.png"], rerun=["fig15_20_23/upsilon_selected_flux_validation.png"], old="plot_selected_upsilon_comparisons_v2.py", new="scripts/plotting/build_gamma_p_validation_figure_set.py", inputs="flux_ratio_scan_W20to26p5_Q2lt0p001.csv", status="Encoded rerender", tone="blue", note="No ROOT dependency; validation values are explicitly stored in CSV."),
    dict(n=21, title="Low/high truth comparison", report=["upsilon_selected_truth_low_high_ratio.png"], rerun=["fig15_20_23/upsilon_selected_truth_low_high_ratio.png"], old="compare_sigma_gp_two_energies_swan_v4.py + selected-plot script", new="scripts/validation/validate_truth_gamma_p_two_energies.py + scripts/plotting/build_gamma_p_validation_figure_set.py", inputs="compare_10x100_18x275_truth_W20to26p5_Q2lt0p001_stateNorm.csv", status="Encoded rerender", tone="blue", note="No hidden macro dependency; comparison values are explicitly stored in CSV."),
    dict(n=22, title="W-dependent reconstruction efficiency", report=["upsilon_selected_18x275_family_efficiency.png"], rerun=["fig15_20_23/upsilon_selected_18x275_family_efficiency.png"], old="extract_sigma_gamma_p_family_fixed_swan_v2.py + selected-plot script", new="scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py + scripts/plotting/build_gamma_p_validation_figure_set.py", inputs="merged_elecs_test_tuples_fix.root", status="Central data exact", tone="green", note="Normalization is byte-exact; efficiencies and central cross sections reproduce."),
    dict(n=23, title="Bin-by-bin closure", report=["upsilon_selected_closure.png"], rerun=["fig15_20_23/upsilon_selected_closure.png"], old="extract_sigma_gamma_p_family_fixed_swan_v2.py + selected-plot script", new="scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py + scripts/plotting/build_gamma_p_validation_figure_set.py", inputs="merged_elecs_test_10t100.root + merged_elecs_test_tuples_fix.root", status="Central data exact", tone="green", note="Promoted output records closure and uncertainty components explicitly."),
    dict(n=24, title="W response-matrix closure", report=["unfold_1S_18x275_closure.png"], rerun=["fig24_25/unfold_1S_18x275_closure.png"], old="unfold_sigma_gp_W_swan.py", new="scripts/unfolding/unfold_gamma_p_cross_section_vs_w.py", inputs="merged_elecs_test_tuples_fix.root", status="Response exact", tone="green", note="Response CSV is byte-exact; toy uncertainty differs slightly with the current numerical stack."),
    dict(n=25, title="W response matrix", report=["unfold_1S_18x275_response_counts.png", "unfold_1S_18x275_response_probability.png"], rerun=["fig24_25/unfold_1S_18x275_response_counts.png", "fig24_25/unfold_1S_18x275_response_probability.png"], old="unfold_sigma_gp_W_swan.py", new="scripts/unfolding/unfold_gamma_p_cross_section_vs_w.py", inputs="merged_elecs_test_tuples_fix.root", status="Data exact", tone="green", note="Count and probability response CSVs reproduce byte-for-byte."),
    dict(n=26, title="W resolution summary", report=["Wres_1S_18x275_full_W_resolution_vs_W.png", "Wres_1S_18x275_full_W_relative_resolution_vs_W.png"], rerun=["fig26_27_w_resolution/Wres_1S_18x275_full_W_resolution_vs_W.png", "fig26_27_w_resolution/Wres_1S_18x275_full_W_relative_resolution_vs_W.png"], old="study_W_resolution_swan.py", new="scripts/validation/validate_w_resolution.py", inputs="merged_elecs_test_tuples_fix.root", status="Data exact", tone="green", note="Resolution summary and response CSVs reproduce byte-for-byte; raster differences are Matplotlib-only."),
    dict(n=27, title="W residual distributions", report=["Wres_1S_18x275_full_deltaW_distributions.png"], rerun=["fig26_27_w_resolution/Wres_1S_18x275_full_deltaW_distributions.png"], old="study_W_resolution_swan.py", new="scripts/validation/validate_w_resolution.py", inputs="merged_elecs_test_tuples_fix.root", status="Data exact", tone="green", note="The added fix ROOT file resolves this figure's exact event input."),
    dict(n=29, title="|t| estimator response metrics", report=["compare_t_methods_1S_twoE_separateW_metrics.png"], rerun=["fig29_36/compare_t_methods_1S_twoE_separateW_metrics.png"], old="compare_t_methods_two_energies_swan_v2.py", new="scripts/validation/benchmark_t_estimators.py", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix.root", status="Corrected near-match", tone="blue", note="0.50% of pixels change because the corrected high tuple adds events; normalized RMSE 0.0377."),
    dict(n=30, title="|t| residual distributions", report=["compare_t_methods_1S_twoE_separateW_residuals.png"], rerun=["fig29_36/compare_t_methods_1S_twoE_separateW_residuals.png"], old="compare_t_methods_two_energies_swan_v2.py", new="scripts/validation/benchmark_t_estimators.py", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix.root", status="Corrected near-match", tone="blue", note="1.07% of pixels change because the corrected high tuple adds events; normalized RMSE 0.0759."),
    dict(n=31, title="|t| unfolding closure", report=["compare_t_methods_1S_twoE_separateW_closure.png"], rerun=["fig29_36/compare_t_methods_1S_twoE_separateW_closure.png"], old="compare_t_methods_two_energies_swan_v2.py", new="scripts/validation/benchmark_t_estimators.py", inputs="merged_elecs_test_10t100_10x100.root + merged_elecs_test_tuples_fix.root", status="Corrected near-match", tone="blue", note="0.94% of pixels change because the corrected high tuple adds events; normalized RMSE 0.0768."),
    dict(n=32, title="t_eX checks at 18x275 GeV", report=["unfold_t_1S_18x275_eX_zeusBins_zeus_overlay.png", "unfold_t_1S_18x275_eX_zeusBins_closure.png", "unfold_t_1S_18x275_eX_zeusBins_t_resolution.png", "unfold_t_1S_18x275_eX_zeusBins_response_counts.png"], rerun=["fig29_36/unfold_t_1S_18x275_eX_zeusBins_zeus_overlay.png", "fig29_36/unfold_t_1S_18x275_eX_zeusBins_closure.png", "fig29_36/unfold_t_1S_18x275_eX_zeusBins_t_resolution.png", "fig29_36/unfold_t_1S_18x275_eX_zeusBins_response_counts.png"], old="unfold_t_distribution_swan_v5.py", new="scripts/unfolding/unfold_t_spectrum.py", inputs="merged_elecs_test_tuples_fix.root", status="Pixel exact", tone="green", note="All four rerun panels are byte-for-byte identical to the report assets."),
    dict(n=33, title="t_eX checks at 10x100 GeV", report=["unfold_t_1S_10x100_eX_zeusBins_zeus_overlay.png", "unfold_t_1S_10x100_eX_zeusBins_closure.png", "unfold_t_1S_10x100_eX_zeusBins_t_resolution.png", "unfold_t_1S_10x100_eX_zeusBins_response_counts.png"], rerun=["fig29_36/unfold_t_1S_10x100_eX_zeusBins_zeus_overlay.png", "fig29_36/unfold_t_1S_10x100_eX_zeusBins_closure.png", "fig29_36/unfold_t_1S_10x100_eX_zeusBins_t_resolution.png", "fig29_36/unfold_t_1S_10x100_eX_zeusBins_response_counts.png"], old="unfold_t_distribution_swan_v5.py", new="scripts/unfolding/unfold_t_spectrum.py", inputs="merged_elecs_test_10t100.root", status="Pixel exact", tone="green", note="All four panels are byte-exact; archived metadata proves method=pT2 despite the eX filename/caption."),
    dict(n=34, title="Merged-bin t_eX checks at 10x100 GeV", report=["unfold_t_1S_10x100_eX_mergedLowT_zeus_overlay.png", "unfold_t_1S_10x100_eX_mergedLowT_closure.png", "unfold_t_1S_10x100_eX_mergedLowT_t_resolution.png", "unfold_t_1S_10x100_eX_mergedLowT_response_counts.png"], rerun=["fig29_36/unfold_t_1S_10x100_eX_mergedLowT_zeus_overlay.png", "fig29_36/unfold_t_1S_10x100_eX_mergedLowT_closure.png", "fig29_36/unfold_t_1S_10x100_eX_mergedLowT_t_resolution.png", "fig29_36/unfold_t_1S_10x100_eX_mergedLowT_response_counts.png"], old="unfold_t_distribution_swan_v5.py", new="scripts/unfolding/unfold_t_spectrum.py", inputs="merged_elecs_test_10t100.root", status="Pixel exact", tone="green", note="All four rerun panels are byte-for-byte identical to the report assets."),
    dict(n=35, title="pT-squared checks at 18x275 GeV", report=["unfold_t_1S_18x275_pT2_mergedLowT_zeus_overlay.png", "unfold_t_1S_18x275_pT2_mergedLowT_closure.png", "unfold_t_1S_18x275_pT2_mergedLowT_t_resolution.png", "unfold_t_1S_18x275_pT2_mergedLowT_response_counts.png"], rerun=["fig29_36/unfold_t_1S_18x275_pT2_mergedLowT_zeus_overlay.png", "fig29_36/unfold_t_1S_18x275_pT2_mergedLowT_closure.png", "fig29_36/unfold_t_1S_18x275_pT2_mergedLowT_t_resolution.png", "fig29_36/unfold_t_1S_18x275_pT2_mergedLowT_response_counts.png"], old="unfold_t_distribution_swan_v5.py", new="scripts/unfolding/unfold_t_spectrum.py", inputs="merged_elecs_test_tuples_fix.root", status="Pixel exact", tone="green", note="All four rerun panels are byte-for-byte identical to the report assets."),
    dict(n=36, title="pT-squared checks at 10x100 GeV", report=["unfold_t_1S_10x100_pT2_mergedLowT_zeus_overlay.png", "unfold_t_1S_10x100_pT2_mergedLowT_closure.png", "unfold_t_1S_10x100_pT2_mergedLowT_t_resolution.png", "unfold_t_1S_10x100_pT2_mergedLowT_response_counts.png"], rerun=["fig29_36/unfold_t_1S_10x100_pT2_mergedLowT_zeus_overlay.png", "fig29_36/unfold_t_1S_10x100_pT2_mergedLowT_closure.png", "fig29_36/unfold_t_1S_10x100_pT2_mergedLowT_t_resolution.png", "fig29_36/unfold_t_1S_10x100_pT2_mergedLowT_response_counts.png"], old="unfold_t_distribution_swan_v5.py", new="scripts/unfolding/unfold_t_spectrum.py", inputs="merged_elecs_test_10t100.root", status="Pixel exact", tone="green", note="All four rerun panels are byte-for-byte identical to the report assets."),
]


def register_fonts() -> tuple[str, str]:
    regular = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
    bold = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")
    if regular.exists() and bold.exists():
        pdfmetrics.registerFont(TTFont("DV", str(regular)))
        pdfmetrics.registerFont(TTFont("DV-Bold", str(bold)))
        return "DV", "DV-Bold"
    return "Helvetica", "Helvetica-Bold"


def wrap(s: str, width: int) -> list[str]:
    return textwrap.wrap(s, width=width, break_long_words=False, break_on_hyphens=False) or [""]


def crop_image(im: Image.Image, mode: str | None) -> Image.Image:
    if not mode:
        return im
    w, h = im.size
    if mode == "left":
        return im.crop((0, 0, int(w * 0.52), h))
    if mode == "right":
        return im.crop((int(w * 0.48), 0, w, h))
    if mode == "top":
        return im.crop((0, int(h * 0.07), w, h))
    return im


def cached_image(path: Path, cache: Path, crop: str | None, max_px: tuple[int, int]) -> Path:
    digest = hashlib.sha256(f"{path}|{crop}|{max_px}".encode()).hexdigest()[:16]
    out = cache / f"{path.stem}_{digest}.jpg"
    if out.exists():
        return out
    with Image.open(path) as src:
        im = src.convert("RGB")
        im = crop_image(im, crop)
        im.thumbnail(max_px, Image.Resampling.LANCZOS)
        im.save(out, "JPEG", quality=92, subsampling=0, optimize=True)
    return out


def draw_fitted(c: canvas.Canvas, path: Path, box: tuple[float, float, float, float], cache: Path, crop: str | None = None) -> None:
    x, y, w, h = box
    cached = cached_image(path, cache, crop, (1700, 1250))
    with Image.open(cached) as im:
        iw, ih = im.size
    scale = min(w / iw, h / ih)
    dw, dh = iw * scale, ih * scale
    c.drawImage(ImageReader(str(cached)), x + (w - dw) / 2, y + (h - dh) / 2, dw, dh, preserveAspectRatio=True, mask="auto")


def draw_grid(c: canvas.Canvas, paths: list[Path], box: tuple[float, float, float, float], cache: Path, crop: str | None = None) -> None:
    x, y, w, h = box
    n = len(paths)
    if n == 1:
        draw_fitted(c, paths[0], box, cache, crop)
        return
    cols = 2
    rows = math.ceil(n / cols)
    gap = 5
    cw = (w - gap * (cols - 1)) / cols
    ch = (h - gap * (rows - 1)) / rows
    for i, path in enumerate(paths):
        col, row = i % cols, i // cols
        bx = x + col * (cw + gap)
        by = y + h - (row + 1) * ch - row * gap
        draw_fitted(c, path, (bx, by, cw, ch), cache, crop)


def tone_color(name: str) -> Color:
    return {"green": GREEN, "amber": AMBER, "red": RED, "blue": BLUE}.get(name, MUTED)


def validate_paths(report_dir: Path, rerun_dir: Path) -> None:
    missing = []
    for f in FIGURES:
        for rel in f["report"]:
            if not (report_dir / rel).exists():
                missing.append(str(report_dir / rel))
        for rel in f["rerun"]:
            if not (rerun_dir / rel).exists():
                missing.append(str(rerun_dir / rel))
    if missing:
        raise FileNotFoundError("Missing comparison inputs:\n" + "\n".join(missing))


def write_index(path: Path) -> None:
    fields = ["figure", "title", "old_script_or_notebook", "new_script_or_notebook", "runtime_inputs", "verification_status", "verification_note"]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for item in FIGURES:
            w.writerow(dict(figure=item["n"], title=item["title"], old_script_or_notebook=item["old"], new_script_or_notebook=item["new"], runtime_inputs=item["inputs"], verification_status=item["status"], verification_note=item["note"]))


def draw_cover(c: canvas.Canvas, regular: str, bold: str) -> None:
    w, h = PAGE
    c.setFillColor(NAVY)
    c.rect(0, 0, w, h, stroke=0, fill=1)
    c.setFillColor(white)
    c.setFont(bold, 30)
    c.drawString(54, h - 86, "Ypsilon ePIC analysis note")
    c.setFont(bold, 24)
    c.drawString(54, h - 123, "report figure vs fresh rerun")
    c.setFont(regular, 13)
    c.setFillColor(HexColor("#CFE3F5"))
    c.drawString(54, h - 153, "Corrected scope: 32 active figures from the newly supplied note")

    cards = [
        ("32", "side-by-side figure pages"),
        ("7", "audited ROOT inputs, including tuples_fix.root"),
        ("0", "promoted scripts loading DVMP_JPsi_Analysis.C"),
    ]
    x0, y0 = 54, h - 270
    for i, (num, text) in enumerate(cards):
        x = x0 + i * 280
        c.setFillColor(HexColor("#173F5F"))
        c.roundRect(x, y0, 245, 78, 8, stroke=0, fill=1)
        c.setFillColor(white)
        c.setFont(bold, 25)
        c.drawString(x + 18, y0 + 39, num)
        c.setFont(regular, 10.5)
        for j, line in enumerate(wrap(text, 28)):
            c.drawString(x + 63, y0 + 46 - 14 * j, line)

    c.setFont(bold, 13)
    c.setFillColor(white)
    c.drawString(54, h - 320, "Scope and interpretation")
    bullets = [
        "Compared: Figures 3-18, 20-27, and 29-36.",
        "Excluded by request: Figures 1, 2, 19, and 28. Figure 37 is inactive LaTeX content.",
        "Left: exact raster included in the note. Right: fresh code rerun or an explicitly labelled CSV-backed rerender.",
        "The added merged_elecs_test_tuples_fix.root is the exact source for the W-resolution tables and several W/t pipelines.",
        "This verifies reproducibility and provenance. It does not promote invalid mass-fit covariances or call-limit fits to final physics results.",
    ]
    c.setFont(regular, 11.5)
    c.setFillColor(HexColor("#DCEAF5"))
    yy = h - 346
    for bullet in bullets:
        c.drawString(62, yy, u"•")
        lines = wrap(bullet, 115)
        for j, line in enumerate(lines):
            c.drawString(80, yy - 15 * j, line)
        yy -= 20 + 15 * (len(lines) - 1)

    c.setFillColor(HexColor("#9FBAD0"))
    c.setFont(regular, 9.5)
    c.drawRightString(w - 42, 26, "Generated by scripts/reporting/build_analysis_note_side_by_side.py")
    c.showPage()


def draw_figure_page(c: canvas.Canvas, item: dict, report_dir: Path, rerun_dir: Path, cache: Path, regular: str, bold: str, page_no: int) -> None:
    w, h = PAGE
    c.setFillColor(white)
    c.rect(0, 0, w, h, stroke=0, fill=1)
    c.setFillColor(NAVY)
    c.rect(0, h - 63, w, 63, stroke=0, fill=1)
    c.setFillColor(white)
    c.setFont(bold, 19)
    c.drawString(28, h - 36, f"Figure {item['n']}")
    c.setFont(regular, 13)
    c.drawString(140, h - 36, item["title"])

    left_x, right_x = 28, 494
    box_y, box_w, box_h = 154, 438, 286
    for x, label, color in [(left_x, "REPORT FIGURE", NAVY), (right_x, "FRESH RERUN / ENCODED RERENDER", TEAL)]:
        c.setFillColor(color)
        c.roundRect(x, h - 91, box_w, 22, 5, stroke=0, fill=1)
        c.setFillColor(white)
        c.setFont(bold, 9.5)
        c.drawCentredString(x + box_w / 2, h - 84, label)
        c.setStrokeColor(BORDER)
        c.setLineWidth(0.8)
        c.roundRect(x, box_y, box_w, box_h, 6, stroke=1, fill=0)

    report_paths = [report_dir / p for p in item["report"]]
    rerun_paths = [rerun_dir / p for p in item["rerun"]]
    draw_grid(c, report_paths, (left_x + 5, box_y + 5, box_w - 10, box_h - 10), cache, item.get("crop"))
    draw_grid(c, rerun_paths, (right_x + 5, box_y + 5, box_w - 10, box_h - 10), cache, item.get("crop"))

    status_y = 124
    color = tone_color(item["tone"])
    c.setFillColor(color)
    c.roundRect(28, status_y, 156, 20, 5, stroke=0, fill=1)
    c.setFillColor(white)
    c.setFont(bold, 9)
    c.drawCentredString(106, status_y + 6, item["status"].upper())
    c.setFillColor(INK)
    c.setFont(regular, 9.4)
    lines = wrap(item["note"], 118)
    for j, line in enumerate(lines[:2]):
        c.drawString(196, status_y + 11 - 12 * j, line)

    c.setStrokeColor(BORDER)
    c.line(28, 113, w - 28, 113)
    c.setFillColor(MUTED)
    c.setFont(bold, 8.5)
    c.drawString(28, 96, "OLD SOURCE")
    c.drawString(28, 70, "PROMOTED SOURCE")
    c.drawString(28, 44, "RUNTIME INPUTS")
    c.setFillColor(INK)
    c.setFont(regular, 8.7)
    for j, line in enumerate(wrap(item["old"], 132)[:2]):
        c.drawString(130, 96 - 10 * j, line)
    for j, line in enumerate(wrap(item["new"], 132)[:2]):
        c.drawString(130, 70 - 10 * j, line)
    for j, line in enumerate(wrap(item["inputs"], 132)[:2]):
        c.drawString(130, 44 - 10 * j, line)
    c.setFillColor(MUTED)
    c.setFont(regular, 8)
    c.drawRightString(w - 28, 18, f"comparison page {page_no} of {len(FIGURES)}")
    c.showPage()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-dir", type=Path, required=True)
    ap.add_argument("--rerun-dir", type=Path, required=True)
    ap.add_argument("--output-pdf", type=Path, required=True)
    ap.add_argument("--output-index", type=Path, required=True)
    ap.add_argument("--cache-dir", type=Path, required=True)
    args = ap.parse_args()
    validate_paths(args.report_dir, args.rerun_dir)
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    args.output_index.parent.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    regular, bold = register_fonts()
    write_index(args.output_index)
    c = canvas.Canvas(str(args.output_pdf), pagesize=PAGE, pageCompression=1)
    c.setTitle("Ypsilon ePIC analysis note: report figure vs fresh rerun")
    c.setAuthor("Ypsilon Fit automated framework")
    draw_cover(c, regular, bold)
    for page_no, item in enumerate(FIGURES, 1):
        draw_figure_page(c, item, args.report_dir, args.rerun_dir, args.cache_dir, regular, bold, page_no)
    c.save()
    print(f"Wrote {args.output_pdf}")
    print(f"Wrote {args.output_index}")


if __name__ == "__main__":
    main()
