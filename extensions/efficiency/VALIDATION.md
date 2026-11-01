# Verification record

This package was verified with controlled synthetic inputs. No CERN AFS/EOS ROOT input was available for the production-sample diagnosis in this environment. Synthetic outputs are not included as physics results.

## Numerical checks

Twelve unit tests pass with `python3 -m unittest discover -s tests -v`:

- Original producer selection and W/rapidity/pT histograms reproduced against the frozen original predicate.
- The Q² split includes the boundary only in the low band and partitions all eligible/selected counts exactly.
- Chunking does not change numerical tables.
- Sequential cutflow and N-minus-one counts correctly handle overlapping failures.
- Nonfinite, −999 sentinel, finite lower/upper rejections, and passing entries form disjoint per-cut partitions of eligible truth.
- The original Σ(E−pz) floating-point boundary predicate is retained exactly.
- Empty denominators are undefined; zero selected with nonzero denominator retains a nonzero Wilson upper limit.
- Under/overflow and last-edge handling are consistent with NumPy histograms and original source.
- Q² mixture contributions exactly sum to the total efficiency.
- Adding arbitrary weight branches leaves counts unchanged.
- Missing, jagged/non-numeric, or inconsistent branch arrays are rejected.

## Full workflow and failure handling

A fresh complete run read four actual synthetic ROOT TTrees with 6,000, 6,000, 4,000, and 3,200 rows, using 128 kB chunks and full input hashes. All invariant checks passed, the two deliberately identical 1S truth populations were recognized, and the run generated eight CSV tables plus twelve PDF/PNG figure pairs and JSON/Markdown audits. The final Q²-mixture and mass-failure breakdown additions were included in this full run. Selected plots were visually inspected for labels, interval display, legends, maps, and numerical conventions. All synthetic figures carry a prominent test-input label.

Additional end-to-end checks showed:

- Reusing a run ID is refused; all existing output hashes remain unchanged.
- A missing input produces PARTIAL while retaining successful samples and plots.
- An intentionally changed reference count produces MISMATCH/PARTIAL while retaining the diagnostic evidence.

The complete workflow ran with Python 3.12 and uproot 5.7.7. All source files and notebook code cells also parse under Python 3.9 syntax. This does not claim a separate Python 3.9 runtime test.

## Scientific scope

These checks establish code behavior on controlled inputs, not the cause of the production-sample efficiency curves. COMPLETE is an execution status. Each report states that the ratios are conditional on eligible stored tuple rows and use finite-MC Wilson intervals. No luminosity, normalization, background model, state-yield fit, or weight branch enters these efficiencies. ID pairing and independence are not assumed. Sentinel values do not authenticate a particular upstream tracking/association failure. The original framework scripts and previously delivered publication notebooks were not modified.


## Standalone final-selection 2D update

- All 22 unit tests passed with Python 3.12. Ten new tests cover exact parity with the frozen original performance selection, nominal pT–rapidity map parity, fixed-mass xV versus fluctuating stored mass, proton-mass denominator, xB cell placement, log underflows and zero Q², Q²-cap boundaries for higher states, independence from pT/y validity for x maps, sparse/empty binomial intervals, and streaming chunk invariance.
- Synthetic ROOT-to-output runs completed for 1S (6,000 input rows) and 3S (3,200 rows) through uproot 5.7.7 with 64 kB chunks. Each produced exactly four PDFs, four PNGs, the count/interval CSV and selection audit. Full input hashing was exercised on 1S.
- All map count/flow invariants passed. Existing run IDs were refused and every original output byte remained unchanged.
- All four rendered 1S synthetic maps were visually inspected for readable labels, log axes, fixed efficiency scale, sparse/undefined-cell styling and clipping. Synthetic previews are marked and are not production results.
- All Python sources and the new SWAN notebook's cells parse using Python 3.9 syntax. This is not a Python 3.9 runtime test.
- Original executable scripts, diagnostics engine, JSON configuration, requirements, notebook04 and original tests are byte-identical to the prior delivered package. Changes are additive scripts/modules/tests/notebook plus this validation record and README instructions.
- CERN production ROOT files are not available in this environment. Their four new maps must be produced on lxplus/SWAN; no x–Q² distributions have been inferred from marginal histograms.


## Explicit proton-category update

The earlier selections and figures were proton-inclusive. New runs now require topology1/2 and positive finite reconstructed-proton momentum unless `--proton-selection inclusive` is explicitly supplied. This changes the reconstructed numerator only. It is not a relabeling of previous results.

- All 29 tests passed. Seven new tests cover the required default, missing schema, all five known topology codes, unknown/noninteger/nonfinite codes, inconsistent momentum including topology0, rejection even for rows outside truth selection, explicit inclusive bypass, unchanged denominators, smaller tagged numerators, before/after cutflow reconciliation, streaming and independence from weight branches.
- Synthetic ROOT integration completed for the full diagnostic and four-map entry point in required mode, plus four-map inclusive mode. Identical truth population: 5,507 eligible rows. Inclusive selected count: 1,741. Tagged selected count: 1,115. Both entry points agree exactly. These figures are test values, not CERN production results.
- Required runs read ten scalar branches; inclusive runs read the original eight. Validation checks tag/momentum consistency over all rows. Required/inclusive categories and raw before/after MC yield counts are retained in audits and tables. Legacy CSV comparisons are explicitly pre-proton inclusive checks.
- A repeated run ID was refused with every existing output byte unchanged. Python sources and both SWAN notebook cells parse as Python3.9 syntax. The tagged plot subtitle and Total efficiency [%] colorbar were visually inspected.
- No CERN production inputs were available locally. Production comparisons and normalization/extraction updates require separate fresh runs. No normalization, branching fraction, flux or event-weight factors were changed by this standalone update.


## Separate proton-tagged notebook copy

The `Ypsilon_Efficiency_Diagnostics_ProtonTagged` copy changes only notebook default output roots/run IDs and documentation. Both notebooks retain `PROTON_SELECTION = "required"`. The original recovered package is preserved. Core modules, entry points, JSON inputs, requirements and tests are byte-identical to the recovered proton-selection update. No production EOS ROOT data were run locally.


## Separate display cells and minor style update

The display-only notebook reads existing count CSV/audit files and embeds the same renderer as `diagnostics/final_maps_report.py`. Four individual cells show PNG previews and relative PNG/PDF links. Requested style changes remove the proton-category subtitle and all in-panel legends, and use `Efficiency [%]` on the colorbar. Proton requirements, count calculations, hatching/empty-cell rendering and synthetic watermarks are retained. Core selection modules and numerical entry points are unchanged.

A synthetic 20,000-row tagged sample was used to execute all display notebook code cells. Four PNGs and four PDFs were produced, eight relative file links resolved, requested labels were checked in PDF text, and a PNG was visually inspected. Repeated execution created a fresh export folder and left earlier plots and both source files byte-identical. An intentionally inconsistent inclusive audit was rejected. Both notebook05 files passed notebook-format validation and all code cells parsed under Python3.9 syntax. This checks display behavior only; no CERN ROOT input or production efficiency values were available.


## Remove hatching from all four maps

Removed the sparse-cell hatch overlay from the module and embedded display-notebook renderer. Populated cells remain colored by efficiency and zero-denominator cells remain gray. Counts, intervals, proton requirements and source files are unchanged. All four synthetic map pairs were rendered with assertions that no hatch artists or legends remain and the colorbar reads `Efficiency [%]`. Sparse-cell masking and embedded-renderer parity were checked; the display notebook passes notebook-format and Python3.9 syntax checks.
