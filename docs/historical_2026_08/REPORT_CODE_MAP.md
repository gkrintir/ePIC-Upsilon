# ePIC Y-meson report-to-code map

## Scope and confidence language

The three supplied PDFs were rendered page by page. Visible plot panels were
compared with `results/reference/` by image similarity and visual inspection,
then traced to retained notebook cells or scripts. Six ROOT files supplied
later were cross-checked against the retained outputs and used to recover all
ten TH2 fill rules. Page numbers below are PDF page numbers.

The `match` label describes the **visible report image**, not whether the
underlying physics result is valid:

| Machine label | Meaning |
|---|---|
| `exact` | Same archived panel pixels/image object, allowing PDF encoding, crop, or scaling. |
| `high-near/source` | Near visual match or the same retained plotting source/data, but not pixel-exact. |
| `medium-lineage` | Later descendant or the same histogram/data series in a changed comparison/style. |
| `low-analogue` | Related observable, but the selection, sample, energy, or comparison is materially different. |
| `missing-producer` | No retained code/artifact reproduces the visible panel. |
| `external` | Figure is credited to an external source and is not produced here. |
| `manual` | Hand-authored table/composite; no generator was identified. |

Only **Performance page 7** has an image-level `exact` report match. In particular,
Performance pages **2, 3, 8, and 10 are not exact**. A source cell listed for a
historical artifact is the closest retained source family; it does not prove
that the report-era run state or inputs survive.

## Reports inspected

| Report ID | PDF | Pages | SHA-256 |
|---|---|---:|---|
| `original` | `Y meson photo-production in ePIC.pdf` | 17 | `0d37e7140c76d89876f5712684fea414f7c51d617dbc1f373a70ed098e834b0e` |
| `followup` | `Followup Y meson photo-production in ePIC.pdf` | 11 | `da26460441c80cfc15688486c5e83061ac494b2fc12238c8be82c75186180941` |
| `performance` | `Performance for Y meson photo-production in ePIC.pdf` | 12 | `b43fa1960a5fa9849579a7fec47d1c33937f73ea78eeb3e7f9f2ed88a70ab4d6` |

## Original report

Page 9 names `DVMP_JPsi_Analysis.C` as the historical reconstruction step.
None of the retained notebooks or Python scripts loads or calls that macro.
The retained figure code operates on either a flattened `UpsilonTree` or
prefilled ROOT histograms. The unavailable provenance layer is the conversion
from raw ePIC output to `UpsilonTree` or stored histograms, plus the exact
report-era inputs. It is not a runtime dependency of the retained code.

| Page | Visible panel(s) | Best retained artifact(s) | Closest code | Match | Evidence/caveat |
|---:|---|---|---|---|---|
| 2 | S. Yoo/Berkeley DSCB mass-fit example | — | — | `external` | Prior external study; its quoted ratios/resolutions are not local fit outputs. |
| 4 | Generated Upsilon mass | `results/reference/Upsilon_TrueLevel_Overlay_Dashboard.png`; `GenLevel_h_Gen_Upsilon_Mass.png` | `notebooks/03_event_selection_truth_and_bt_studies.ipynb`, cell 7; notebook 07, cell 0 | `medium-lineage` | Cell 7 fills the distribution directly from `UpsilonTree`; notebook 07 only renders a later pre-binned campaign view. |
| 5 | Generated Upsilon rapidity; generated Upsilon pT | `Upsilon_TrueLevel_Overlay_Dashboard.png`; `GenLevel_h_Gen_Upsilon_y.png`; `GenLevel_h_Gen_Upsilon_pT.png` | notebook 03, cell 7; notebook 07, cell 0 | `medium-lineage` | Tuple-direct descendants survive; the exact single-sample report layout does not. |
| 6 | W vs Upsilon rapidity | `GenLevel_h2_Gen_W_vs_Rapidity_TH2.png` | `scripts/validation/reconstruct_report_histograms_from_tuple.py`; `scripts/plotting/compare_report_campaign_histograms.py` | `medium-lineage` | Exact stored TH2 and exact tuple fill recovered: x=`true_y`, y=`true_W`. Report-layout pixels are not claimed. |
| 6 | Generated \|t\| vs W | `GenLevel_h2_Gen_t_vs_W_TH2.png`; `Upsilon_TrueLevel_Overlay_Dashboard.png` | reconstruction and campaign scripts; notebook 03, cell 7 | `medium-lineage` | Exact stored TH2 and tuple fill recovered; the dashboard is a later two-energy contour view. |
| 7 | Decay-lepton pT correlation; eta vs pT | generator-level TH2 comparison files; `Upsilon_TrueLevel_Overlay_Dashboard.png` | reconstruction and campaign scripts; notebook 03, cell 7 | `medium-lineage` | All three generator-level decay-lepton TH2 objects reproduce bin-for-bin from the tuple. |
| 8 | Generated scattered-electron energy; eta | `Upsilon_TrueLevel_Overlay_Dashboard.png`; `GenLevel_h_Gen_ScatElec_E.png`; `GenLevel_h_Gen_ScatElec_Eta.png` | notebook 03, cell 7; notebook 07, cell 0 | `medium-lineage` | Cell 7 is tuple-direct; notebook 07 only renders later pre-binned campaign histograms. |
| 9 | Truth/reco decay-lepton maps and efficiency map | truth/reco TH2 comparison files; `DecayLepton_Efficiency_Eta_vs_pT_new.png` | reconstruction and campaign scripts | `medium-lineage` | Truth and reco eta-vs-pT maps are exact tuple reconstructions. The renderer's reco/truth ratio is a source-backed candidate; the report's empty-bin convention is unresolved. |
| 10 | Generated vs reconstructed scattered-electron eta | `Upsilon_TrueLevel_Overlay_Dashboard.png`; `Upsilon_Inclusive_Dashboard.png`; `h_eta_ScatElec.png` | notebook 03, cell 7 and cells 0/1; notebook 07, cell 0 | `low-analogue` | Truth and reco distributions are tuple-direct, but the exact overlay is absent; the stored `h_eta` image compares campaigns. |
| 11 | Electron E/p ratio | `h_EleID_EpRatio.png` | notebook 07, cell 0 | `medium-lineage` | Pre-binned renderer only; the event-to-histogram filling step is absent. |
| 12 | Reconstructed invariant mass | `Upsilon_Inclusive_Dashboard.png`; `h_invMass_Reco_All.png` | notebook 03, cells 0/1; notebook 07, cell 0 | `medium-lineage` | Notebook 03 fills the mass distribution from `UpsilonTree`; notebook 07 renders a later pre-binned comparison. |
| 12 | Generated/reconstructed W and Upsilon rapidity | `Upsilon_TrueLevel_Overlay_Dashboard.png`; `Upsilon_Inclusive_Dashboard.png`; `h_Reco_W.png`; `h_UpsilonRapidity.png` | notebook 03, cell 7 and cells 0/1; notebook 07, cell 0 | `low-analogue` | Tuple-direct truth and reco ingredients survive, but not the exact overlays. |
| 12 | W-vs-rapidity inset | `GenLevel_h2_Gen_W_vs_Rapidity_TH2.png` | reconstruction and campaign scripts | `medium-lineage` | The exact TH2 payload and tuple fill survive; the report inset styling is not pixel-matched. |
| 13 | Mass resolution in 1<y<2, -1<y<0, and 0<y<1 | `RapidityBins_h_MassRes_Bin4.png`; `...Bin2.png`; `...Bin3.png` | notebook 07, cell 0 | `medium-lineage` | Pre-binned renderer only; the tuple-level mass-resolution filler is absent. |
| 14 | Inclusive fit without electron ID | `Fit_Inclusive_Upsilon_DSCB_Nice.png` | `notebooks/01_mass_fits_inclusive_and_rapidity.ipynb`, cell 1 | `low-analogue` | Tuple-direct fit analogue; its exact relation to the report's no-ID selection is unresolved. |
| 14 | Inclusive fit with electron ID | `Fit_Inclusive_Upsilon_DSCB.png` | `notebooks/01_mass_fits_inclusive_and_rapidity.ipynb`, cell 0 | `high-near/source` | Strong near visual match (same ROOT-style fit/data); the generic output's exact run state is unresolved and current cell text writes a scaled variant. |
| 17 | Fits in -1<y<0 and 0<y<1 | `Fit_Inclusive_Upsilon_RapidityBins.png`; `..._raw_visible_decay.png` | notebook 01, cell 2 | `medium-lineage` | Same DSCB model family; retained composite differs in style and normalization from the report panels. |

## Follow-up report

These pages are report-era ROOT histogram screenshots. The closest retained
one-dimensional files are mostly later campaign comparisons, not exact report
exports.

| Page | Visible panel | Best retained artifact | Closest code | Match | Evidence/caveat |
|---:|---|---|---|---|---|
| 3 | Number of low-Q2 tagger tracks | `Upsilon_Inclusive_Dashboard.png`; `h_nTracks_Tagger.png` | `notebooks/03_event_selection_truth_and_bt_studies.ipynb`, cells 0/1; notebook 07, cell 0 | `medium-lineage` | Notebook 03 fills the distribution from the tuple; notebook 07 renders a later pre-binned campaign view. |
| 4 | Central-vs-tagger track correlation | `h2_nTracks_Central_vs_Tagger_TH2.png`; `central_tagger_corr_18x275_1S.png` | reconstruction and campaign scripts; `plot_track_multiplicity_correlation.py` | `medium-lineage` | The all-event raw map is now exactly tuple-reconstructed. The selected/weighted 1S view remains a separate descendant. |
| 5 | Number of central tracks | `Upsilon_Inclusive_Dashboard.png`; `h_nTracks_Central.png` | notebook 03, cells 0/1; notebook 07, cell 0 | `medium-lineage` | Notebook 03 is tuple-direct; notebook 07 is a pre-binned campaign renderer. |
| 6 | Tagged-proton momentum | `h_Reco_Proton_p.png` | notebook 07, cell 0 | `medium-lineage` | Pre-binned renderer only; the event-to-histogram filling step is absent. |
| 7 | Reconstructed \|t\| for 80<W<100 GeV | `h_Reco_t_Wbin2.png`; `reco_t_distributions_1S_twoE_separateW.png` | notebook 07, cell 0; `scripts/plotting/plot_reco_t_estimators_by_energy.py` | `low-analogue` | Notebook 07 renders the stored key; the script is a configurable tuple-direct alternative, but not the exact report layout. |
| 8 | Electron-ID isolation fraction | `h_EleID_IsoFraction.png` | notebook 07, cell 0 | `medium-lineage` | Pre-binned renderer only; the event-to-histogram filling step is absent. |
| 9 | Generated W in acceptance at lower beam energy | `GenLevel_h_Gen_W_InAcc.png`; `Upsilon_TrueLevel_Overlay_Dashboard.png` | notebook 07, cell 0; notebook 03, cell 7 | `low-analogue` | The accepted-W view is pre-binned; cell 7 fills raw true W directly but does not encode the exact report acceptance. |

## Performance report

| Page | Visible panel(s) | Best retained artifact(s) | Closest code | Match | Evidence/caveat |
|---:|---|---|---|---|---|
| 2 | Two rapidity-bin mass-resolution campaign comparisons | `RapidityBins_h_MassRes_Bin2.png`; `...Bin3.png` | `scripts/plotting/compare_report_campaign_histograms.py` | `high-near/source` | Uploaded inputs reproduce all 43 archived nonempty TH1 PNGs byte-for-byte in historical raw-count mode; slide panels still are not pixel-exact copies. |
| 3 | Central and low-Q2 tagger track campaign comparisons | `h_nTracks_Central.png`; `h_nTracks_Tagger.png` | campaign comparison script | `high-near/source` | Uploaded inputs reproduce the archived PNGs byte-for-byte; slide panels still differ by crop or encoding. |
| 4 | Raw eSTARlight vs afterburner pT, rapidity, W, and Q2 | — | — | `missing-producer` | No matching artifact or source was found. |
| 5 | Cross-section and branching-ratio table | — | Constants recur in notebook 06 and extraction scripts | `manual` | No table generator was identified. |
| 6 | Raw ROOT event-selection cutflow | `Upsilon_EventSelection_Cutflow.png`; newer `..._with_Loss.png` | `notebooks/03_event_selection_truth_and_bt_studies.ipynb`, cell 8 (preferred over duplicate cell 10) | `medium-lineage` | Cell 8 reads prefilled `h_CutFlow`; the event-to-histogram filling step and exact raw export are absent. |
| 7 | Four 1S efficiency panels versus pT, y, W, and \|t\| | `Upsilon_1S_Kinematic_Efficiencies-Copy1.png` | `notebooks/02_efficiency_and_two_energy_studies.ipynb`, cell 0 | `exact` | `merged_elecs_test_tuples.root` reproduces the 1S, 2S, and 3S `-Copy1` files byte-for-byte. Input and run state are now identified. |
| 8 | 4x4 “Scaled yields” dashboard | `Upsilon_Inclusive_Dashboard.png` | `notebooks/03_event_selection_truth_and_bt_studies.ipynb`, cell 0; cell 1 is the current audited replacement | `medium-lineage` | Grid, variables, and labels correspond, but mass normalization and several distributions differ; not exact. |
| 9 | ZEUS mass spectrum | — | — | `external` | External ZEUS plot. |
| 10 | Projected gamma-p cross sections for 1S/2S/3S | `Upsilon_Family_GammaP.png` | `notebooks/06_gamma_p_cross_section_family.ipynb`, cell 0 | `high-near/source` | `merged_elecs_test_tuples.root` reproduces the retained image byte-for-byte. The report panel still differs in labels, range, or annotation. |

## Reproducibility and validity limits

- `exact` is an image-provenance statement only. It is not validation of the
  inputs, normalization, fit convergence, or physics result.
- Every mass-fit number remains **unvalidated**. Cached notebook output includes
  MIGRAD call-limit/abnormal-termination messages and non-positive-definite
  covariance reports; some plotted components can become negative. See
  `docs/FIT_MODEL.md` and `docs/NOTEBOOK_AUDIT.md`.
- The inclusive-fit code applies STARlight 2S/3S yield factors after fitting in
  at least one retained path. A visually similar report curve therefore does
  not establish the fitted parameter interpretation.
- `notebooks/03_event_selection_truth_and_bt_studies.ipynb`, cells 0/1 and 7
  reconstruct many reco and truth observables directly from `UpsilonTree`.
- `notebooks/07_compare_campaign_histograms.ipynb`, cell 0 does not fill any
  histogram and its saved source has drifted from the archived PNGs. The
  promoted `compare_report_campaign_histograms.py` restores the exact raw-count
  TH1 state and adds TH2 rendering.
- `reconstruct_report_histograms_from_tuple.py` encodes the selections, branch
  order, axes, and unit-weight fills for all ten TH2 objects. All ten match
  `merged_elecs_test_.root` exactly, including flow bins.
- Full ROOT identities and rerun evidence are in
  `docs/UPLOADED_ROOT_CROSSCHECK.md`.
- The concise old-name to new-name index is
  `docs/FIGURE_SOURCE_RENAME_TABLE.md`, with a CSV counterpart in
  `docs/figure_source_rename_map.csv`.
- `config/report_figures.toml` carries the same panel-level labels for automated
  checks. Missing and external panels must not be silently replaced by a later
  analogue.
