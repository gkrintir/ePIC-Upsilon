# Code audit and rename map

## Scope and attribution limit

The supplied source archive contains 40 Python scripts, 14 notebooks, and 5
editor backup files. The original Python scripts collapse into 21 logical
families. The framework adds two recovered-input utilities, for 23 canonical
scripts overall. There is no report PDF/source, README, input ROOT file,
environment lock, or test suite in the original archive. Three separately
supplied PDFs and six later-supplied ROOT files were therefore inspected
against every retained artifact. Producer attribution also used ROOT schemas,
bin-for-bin reconstruction, CSV schemas, cached notebook output, and exact
version comparisons.

The page-level comparison supports one essentially exact retained-image match
(Performance page 7), plus several source or lineage matches. Many report
observables are reconstructed directly from `UpsilonTree` by retained Python
cells. A smaller set exists only as pre-binned ROOT histograms, and some exact
report layouts are absent. None of the retained code calls or loads
`DVMP_JPsi_Analysis.C`. All ten stored TH2 fills are now encoded directly from
`UpsilonTree`. See `REPORT_CODE_MAP.md`, `FIGURE_SOURCE_RENAME_TABLE.md`, and
`UPLOADED_ROOT_CROSSCHECK.md`.

## Promoted scripts

| Analysis role | Original implementation | Semantic path |
|---|---|---|
| eSTARlight flux validation | `scan_estarlight_flux_ratio_v3.py` | `scripts/flux/validate_estarlight_flux_ratio.py` |
| 1S bin-by-bin extraction, legacy | `extract_sigma_gp_1S_renorm_swan_v2.py` | `scripts/cross_sections/extract_gamma_p_cross_section_1s_legacy.py` |
| Family bin-by-bin extraction | `extract_sigma_gamma_p_family_fixed_swan_v2.py` | `scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py` |
| Differential bin-by-bin extraction | `extract_differential_gamma_p_swan_v4.py` | `scripts/cross_sections/extract_differential_gamma_p_bin_by_bin.py` |
| Bayesian W unfolding | `unfold_sigma_gp_W_swan.py` | `scripts/unfolding/unfold_gamma_p_cross_section_vs_w.py` |
| Bayesian t unfolding | `unfold_t_distribution_swan_v5.py` | `scripts/unfolding/unfold_t_spectrum.py` |
| Common-W hybrid t unfolding | `unfold_t_hybrid_commonW_swan_nopandas.py` | `scripts/unfolding/unfold_hybrid_t_common_w.py` |
| Truth gamma-p comparison | `compare_sigma_gp_two_energies_swan_v4.py` | `scripts/validation/validate_truth_gamma_p_two_energies.py` |
| t-estimator benchmark | `compare_t_methods_two_energies_swan_v2.py` | `scripts/validation/benchmark_t_estimators.py` |
| W-resolution study | `study_W_resolution_swan.py` | `scripts/validation/validate_w_resolution.py` |
| Gamma-p normalization diagnostic | `debug_sigma_gp_1S_swan.py` | `scripts/validation/validate_gamma_p_normalization.py` |
| Fit/yield accounting check | `crosscheck_massfit_vs_sigma_yields.py` | `scripts/validation/validate_mass_fit_yield_accounting.py` |
| Track-multiplicity correlation | `plot_central_tagger_correlation.py` | `scripts/validation/plot_track_multiplicity_correlation.py` |
| Gamma-p energy overlay | `overlay_sigma_gp_two_energies_swan.py` | `scripts/plotting/plot_gamma_p_cross_sections_by_energy.py` |
| Selected validation figure set | `plot_selected_upsilon_comparisons_v2.py` | `scripts/plotting/build_gamma_p_validation_figure_set.py` |
| Different-W hybrid t overlay | `make_hybrid_t_low_high_overlay_v3.py` | `scripts/plotting/plot_hybrid_t_different_w.py` |
| t and impact-parameter plots | `overlay_pt2_twoE_with_ZEUS_and_bT.py` | `scripts/plotting/plot_t_spectra_and_impact_parameter.py` |
| 1S b-slope HERA overlay | `overlay_b_slopes_on_hera_fig4.py` | `scripts/plotting/plot_b_slope_1s_hera_overlay.py` |
| Family b-slope HERA overlay | `overlay_b_slopes_three_states_on_hera_fig4.py` | `scripts/plotting/plot_b_slopes_family_hera_overlay.py` |
| Reconstructed t comparison | `compare_reco_t_distributions_two_energies_swan.py` | `scripts/plotting/plot_reco_t_estimators_by_energy.py` |
| Two-energy t-method overlay | `Upsilon_1S_TMethods_TwoEnergy_Overlay_visibleEE.py` | `scripts/plotting/plot_t_estimators_two_energy_overlay.py` |
| Campaign TH1/TH2 comparison | `MassFit-Copy8.ipynb` C0 | `scripts/plotting/compare_report_campaign_histograms.py` |
| Report TH2 tuple reconstruction | Recovered from uploaded tuple/histogram equality | `scripts/validation/reconstruct_report_histograms_from_tuple.py` |

## Promotion rules

- Version suffixes were removed only after selecting the most complete logical
  implementation in each family.
- Exact original files remain in `archive/original_code/` for line-by-line
  comparison.
- Five `.py~` editor backups are archival only. One contains an `np.trapzd`
  typo and none is considered authoritative.
- `unfold_t_distribution_swan_v2.py` and `v3.py` are byte-identical and both
  fail Python compilation because a normal docstring contains an unescaped
  `\Upsilon` sequence.
- `MassFit-CopyN.ipynb` names are notebook forks, not a chronological version
  guarantee. The curated names describe subject matter instead.
- `MassFit-Copy13.ipynb` is a tau/HepMC study and is outside the Upsilon
  pipeline.

## Safety changes to promoted code

The physics calculations were preserved. Narrow validation changes were made:

- The common-W unfolding now rejects absent state weight, nonpositive
  luminosity/target normalization, empty common phase space, and zero selected
  reconstructed events.
- The family and differential bin-by-bin extractors reject absent requested
  states and invalid normalization samples.
- The selected-figure builder now requires both input CSV paths instead of
  silently loading a known-invalid low-energy artifact.
- The two-energy t overlay now receives files and normalizations from the
  command line rather than embedded paths and constants.

These changes prevent successful-looking zero/NaN outputs. They do not certify
the physics content of the retained reference products.

## Producer evidence

The strongest producer matches are fields or files written by only one version:

- State-scoped metadata identifies the v4 two-energy truth comparison.
- Efficiency-uncertainty and closure columns identify the v2 family extractor.
- Finite-bin b-slope tables identify the v4 differential extractor.
- Response-probability CSVs and unique normalization metadata identify the
  no-Pandas common-W implementation.
- Cached notebook cell output identifies the mass-fit, efficiency, cutflow,
  dashboard, and t-binning plots.

The complete file-level disposition is in `code_inventory.csv`; artifact-level
confidence and unresolved outputs are in `report_output_map.csv`.
