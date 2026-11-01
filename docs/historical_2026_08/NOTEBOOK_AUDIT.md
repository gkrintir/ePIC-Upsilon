# Notebook lineage and authority audit

> **Interpretation:** this audit identifies retained code and its output provenance. It does not validate the physics results. In particular, every cached invariant-mass fit number is unvalidated; see `FIT_MODEL.md`.

## Scope and lineage

The source archive contained 14 notebooks named `MassFit.ipynb` and `MassFit-Copy1.ipynb` through `MassFit-Copy13.ipynb`. The copy number is not a chronological version number. All notebooks have the same original first-cell identifier (`810352e0`), and nearly all retain the same trailing identifiers (`870177de` and `ae40a7b7`), while the cell bodies were overwritten for unrelated analyses. File modification times also contradict copy-number order.

The framework keeps descriptive, source-equivalent copies of the useful notebooks. The original notebooks remain available under the archive.

| Framework notebook | Original notebook | Main purpose | Audit disposition |
| --- | --- | --- | --- |
| `01_mass_fits_inclusive_and_rapidity.ipynb` | `MassFit.ipynb` | Inclusive and rapidity-binned DSCB mass fits | Fit source of record; numerical outputs unvalidated |
| `02_efficiency_and_two_energy_studies.ipynb` | `MassFit-Copy2.ipynb` | Kinematic, decay-lepton, and two-energy efficiencies | Mixed authority by cell |
| `03_event_selection_truth_and_bt_studies.ipynb` | `MassFit-Copy4.ipynb` | Dashboards, cutflow, truth overlays, and bT transform | Mixed authority by cell |
| `04_low_q2_tagger_reach.ipynb` | `MassFit-Copy5.ipynb` | Low-Q2 scattered-electron reach | Unexecuted draft; no archived product |
| `05_differential_rapidity_cross_section.ipynb` | `MassFit-Copy6.ipynb` | Truth/visible differential rapidity cross section | Exploratory predecessor |
| `06_gamma_p_cross_section_family.ipynb` | `MassFit-Copy7.ipynb` | Early Upsilon-family gamma-p extraction | Exploratory predecessor |
| `07_compare_campaign_histograms.ipynb` | `MassFit-Copy8.ipynb` | Campaign histogram comparison | Saved source drifted from archived outputs; promoted script is canonical |
| `08_gamma_p_cross_section_10x100.ipynb` | `MassFit-Copy9.ipynb` | Corrected 10x100 gamma-p extraction | Conflicts with later variants; not canonical |
| `09_gamma_p_cross_section_two_energies.ipynb` | `MassFit-Copy10.ipynb` | Gamma-p variants, theory overlay, and Q2 check | Mixed exploratory cells; not canonical extraction |
| `10_t_migration_diagnostics.ipynb` | `MassFit-Copy11.ipynb` | First t-resolution/migration study | Superseded by notebook 11 |
| `11_t_binning_optimization.ipynb` | `MassFit-Copy12.ipynb` | Resolution-aware t-bin optimization | Current notebook source for bin design |
| `12_basic_kinematic_distributions.ipynb` | `MassFit-Copy3.ipynb` | Basic mass, projective-t, and W plots | Prefer source over Copy1; discard cached images |

`MassFit-Copy1.ipynb` is an older version of the basic-distribution study and was not promoted. `MassFit-Copy13.ipynb` is an unrelated tau-decay/HepMC study and belongs outside the Upsilon-fit workflow.

## Authoritative-cell shortlist

Cells are zero-indexed, matching the notebook JSON.

| Framework notebook | Authoritative cell(s) | Use | Conditions and caveats |
| --- | ---: | --- | --- |
| `01_mass_fits_inclusive_and_rapidity.ipynb` | 0--2 | Exact fit implementation and figure provenance | Source is authoritative; cached fit numbers are not validated |
| `02_efficiency_and_two_energy_studies.ipynb` | 0 | Recovered historical per-state kinematic efficiencies | `merged_elecs_test_tuples.root` reproduces all three `-Copy1` images byte-for-byte |
| `02_efficiency_and_two_energy_studies.ipynb` | 3 | Two-energy, per-state kinematic-efficiency overlays | Uses `...10t100_10x100.root` and `...tuples_fix2.root` with energy-specific sigma-EPz cuts |
| `02_efficiency_and_two_energy_studies.ipynb` | 5 | Event-level decay-lepton efficiencies | Strict 18x275 selection; outputs state plots and CSV points under `decay_lepton_eff_notebook/` |
| `03_event_selection_truth_and_bt_studies.ipynb` | 0 and 1 | Legacy and updated inclusive dashboards | Legacy and `fix2` inputs reproduce their respective retained outputs byte-for-byte; cell 1 is the current source mode |
| `03_event_selection_truth_and_bt_studies.ipynb` | 2 | 1S two-energy reconstructed dashboard | Output name records target-visible normalization |
| `03_event_selection_truth_and_bt_studies.ipynb` | 3 | Side-by-side two-energy t-method dashboard | Produces plot and summary CSV |
| `03_event_selection_truth_and_bt_studies.ipynb` | 7 | Truth-level two-energy dashboard | Current truth-overlay source |
| `03_event_selection_truth_and_bt_studies.ipynb` | 8 | Event-selection cutflow | Prefer over near-duplicate cell 10 |
| `03_event_selection_truth_and_bt_studies.ipynb` | 9 and 11 | bT/Hankel-transform implementation and invocation | Keep only when this optional analysis is required; invocation writes externally |
| `11_t_binning_optimization.ipynb` | 0 | Resolution-aware t-bin scan | Stored recommendation is diagnostic, not independently validated |
| `12_basic_kinematic_distributions.ipynb` | 0 | Current source for mass, projective-t, and W plots | Rerun from inputs; embedded plots are stale |

For the two-energy t-method overlay formerly in notebook 02 cell 4, prefer the standalone `Upsilon_1S_TMethods_TwoEnergy_Overlay_visibleEE.py`: it is almost identical but corrects the low-energy W minimum from 10 to 20 GeV.

## Cell-to-output evidence

The strongest report-like provenance comes from exact save statements, matching cached stream text, output modification times shortly before notebook saves, and descriptive performance-plot filenames. No report, TeX source, Word file, or figure manifest was included, so association with a report is inferred rather than directly proven.

| Source | Input/selection summary | Output evidence |
| --- | --- | --- |
| Notebook 01, cell 0 | `~/merged_elecs_test.root`, `h_invMass_Reco_All` | `Fit_Inclusive_Upsilon_DSCB_Scaled.pdf`; cached ROOT stream says it was created |
| Notebook 01, cell 1 | 10x100 tuple, mass 8--12 GeV, `reco_W >= 0`, visible-decay weight | `Fit_Inclusive_Upsilon_DSCB_Nice.png`; cached stream says it was saved |
| Notebook 01, cell 2 | 18x275 fixed tuple, mass/sigma-EPz/W cuts, two rapidity bins | `Fit_Inclusive_Upsilon_RapidityBins_raw_visible_decay.png` under current modes; target-normalized predecessor also exists |
| Notebook 02, cell 0 | Legacy 18x275 tuple and per-state selections | Three `Upsilon_{state}_Kinematic_Efficiencies-Copy1.png` files reproduce byte-for-byte |
| Notebook 02, cell 3 | 10x100 and 18x275 fixed tuples, state windows and energy-specific cuts | `Upsilon_{1S,2S,3S}_Kinematic_Efficiencies_TwoEnergyOverlay.png` |
| Notebook 02, cell 5 | 18x275 tuple, event-level two-lepton acceptance | Three event-efficiency plots, three denominator plots, and three point CSVs |
| Notebook 03, cell 0 | Legacy 18x275 tuple | `Upsilon_Inclusive_Dashboard.png` reproduces byte-for-byte |
| Notebook 03, cell 1 | Fixed 18x275 tuple, raw visible-ee mode | `Upsilon_Inclusive_Dashboard_raw_visibleEE.png` reproduces byte-for-byte |
| Notebook 03, cell 2 | Two-energy 1S reconstruction | `Upsilon_1S_TwoEnergy_Overlay_Dashboard_target_visible_xsec_visibleEE.png` |
| Notebook 03, cell 3 | Two-energy eX/BABE/eXBABE/projective/pT2 comparison | `Upsilon_1S_TMethods_TwoEnergy_SideBySide_target_visible_xsec_visibleEE.png` and summary CSV |
| Notebook 03, cell 7 | Two-energy truth-level comparison | `Upsilon_TrueLevel_Overlay_Dashboard.png` |
| Notebook 03, cell 8 | ROOT `h_CutFlow` histograms | `Upsilon_EventSelection_Cutflow_with_Loss.png` |
| Notebook 06, cell 0 | Legacy 18x275 tuple | `Upsilon_Family_GammaP.png` reproduces byte-for-byte |
| Notebook 07, cell 0 | Two campaign histogram files | Current source does not reproduce archived PNGs; promoted raw-count renderer does |
| Notebook 11, cell 0 | `~/merged_elecs_test_tuples_fix.root`, strict reconstructed selection | Calibration plots, migrations, scan tables, and recommendation CSV under `t_binning_resolution_aware_v4/` |
| Notebook 12, cell 0 | `~/merged_elecs_test_tuples.root`, strict RDataFrame selection | `Upsilon_Mass_CombinedTopologies.png`, `Upsilon_Projective_t.png`, `Upsilon_W_Spectrum.png` after a required fresh run |

Useful timestamp corroboration includes the scaled-fit PDF on May 22, rapidity-fit variants on May 24, current dashboards between 12:35 and 14:36 on May 25, two-energy efficiency overlays at 15:49 on May 25, and the “Nice” mass-fit PNG at 18:12 immediately before the mass-fit notebook was saved. Timestamps support provenance but do not establish numerical validity.

## Excluded, unsafe, or superseded cells

| Source | Reason not to automate as authoritative code |
| --- | --- |
| Notebook 02, cell 1 | Writes the same generic state-efficiency filenames as cell 0 with different cuts. Keep it separate from the now-verified cell-0 `-Copy1` outputs. |
| Notebook 02, cell 2 | Current source reads a 10x100 file while carrying 18x275 labels. Its cached run found only 1S and skipped 2S/3S, although top-level 2S/3S plots exist; those files came from overwritten code. |
| Notebook 02, cell 4 | Superseded by the corrected standalone t-method overlay script. |
| Notebook 03, cell 0 | Historical dashboard source. Its exact input is now known, but cell 1 remains the current audited replacement. |
| Notebook 03, cell 4 | Reads a 10x100 file but hard-codes 18x275 target cross sections and labels. The resulting raw dashboard visibly behaves like 10x100 while claiming 18x275; exclude it. |
| Notebook 03, cell 5 | Source is empty but a cached display remains. It cannot reproduce its output. |
| Notebook 03, cell 6 | Obsolete 4x5 dashboard with duplicated plotting statements; collides with cell 0 on `Upsilon_Inclusive_Dashboard.png`. |
| Notebook 03, cell 10 | Approximately 96% source-similar to cutflow cell 8 and writes the same filename; use cell 8. |
| Notebook 04, cell 0 | Unexecuted low-Q2 reach draft that writes an external ROOT file; no matching archived product. |
| Notebook 05, cell 0 | Exploratory single-energy rapidity cross section superseded by maintained differential scripts. |
| Notebook 06, cell 0 | Primitive gamma-p extraction superseded by corrected extraction scripts; cached plot code combines log-y with `ylim(0, 600)`, generating a warning. |
| Notebooks 08 and 09 | Conflicting input defaults, Q2 maxima, branching-ratio conventions, and theory sources; use the canonical standalone gamma-p pipeline instead. |
| Notebook 10, cell 0 | Older migration study lacks the sigma-EPz cut, fails to expand `~` before creating the output directory, and is superseded by notebook 11. |
| Notebook 12 cached outputs | The archived Copy1 and Copy3 notebooks contain byte-identical embedded images despite materially different source selections. Treat the cache as stale and rerun cell 0. |
| Archived `MassFit-Copy13.ipynb` | Tau-decay/HepMC study unrelated to the Upsilon report workflow. |

## Output collisions and reproducibility hazards

- The old Copy1 and Copy3 notebooks write the same three basic-distribution filenames.
- Notebook 02 cells 0 and 1 overwrite the same state-efficiency filenames.
- Notebook 03 cells 0 and 6, plus an empty cached predecessor, compete for `Upsilon_Inclusive_Dashboard.png`.
- Notebook 03 cells 8 and 10 overwrite the same cutflow plot.
- Gamma-p notebook variants write the same external “corrected” filename while using different conventions.
- Inputs use ambiguous aliases. The six files supplied for this cross-check now
  have explicit roles and checksums in `config/verified_root_inputs.toml`; any
  remaining `old` or `fix` aliases are still unresolved.

Automation should therefore resolve inputs through a dataset manifest, derive plot labels from that manifest, include selection and normalization identifiers in output paths, clear notebook outputs before execution, and emit machine-readable run metadata. A cached notebook display or an existing PNG must never count as proof that the current source reproduced it.

## Fit-result warning

The cached invariant-mass fits include `STATUS=CALL LIMIT`, `STATUS=NOT POSDEF`, `ERR MATRIX NOT POS-DEF`, an abnormal-minimization warning, negative fitted amplitudes, and a negative integrated 2S yield. Consequently, **the archived masses, widths, yields, uncertainties, and covariance-derived quantities are unvalidated and must not be reported as final results**. Only a fresh run that passes explicit minimizer, covariance, physicality, and stability gates may promote those values.
