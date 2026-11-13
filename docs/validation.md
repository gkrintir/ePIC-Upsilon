# Repository validation

Run `python3 tools/run_checks.py` from the repository root after installing
`requirements-dev.txt`. Each component runs in its own process and working
directory so historical modules with the same names cannot shadow one another.

For the website, also run `node tools/check_analysis_trace.mjs` and
`node tools/check_atlas.mjs`. These check the guide routes and source links,
current rematched numerical data and figure assets, and historical-import
compatibility.

## Full-sample rematched projections

The current Atlas and analysis note use all six downstream EICrecon #3011
rematched samples with the full forward ECAL matching mode. The input SHA-256 values were
verified before analysis. Normalization uses each corrected file's actual
generated-row count. The analysis ran at framework revision `072c1fd`, and the
source manifest accompanies the current numerical exports.
[Inputs, source hashes and per-observable validation](data/paper-provenance.json)

| State | Generated rows at 9×130 | Generated rows at 9×275 |
|---|---:|---:|
| 1S | 1,007,500 | 1,137,600 |
| 2S | 1,058,300 | 1,280,400 |
| 3S | 1,041,900 | 1,254,900 |

For the proton-tagged 9×275 GeV projection at 2.5 fb⁻¹, the current results are
as follows. The weighted selected count sums the accepted simulated events.
The fitted yield comes from the simultaneous three-state mass model.

| Quantity | 1S | 2S | 3S |
|---|---:|---:|---:|
| Weighted selected events | 91.6 | 27.3 | 22.7 |
| Fitted mass yield | 91.8 (+10.2 / −9.5) | 27.1 (+6.5 / −5.7) | 22.7 (+5.8 / −5.0) |
| Fitted core width in MeV | 54.4 (+6.9 / −6.1) | 60.3 (+7.7 / −6.7) | 63.2 (+8.0 / −7.1) |
| pT² slope status | `PASS_CONDITIONAL` | `NOT_VALIDATED` | `NO_FEASIBLE_SLOPE_BINNING` |
| Adaptive W spectrum status | `NOT_VALIDATED` | `NOT_VALIDATED` | `NOT_VALIDATED` |
| Rapidity spectrum status | `NOT_VALIDATED` | `NOT_VALIDATED` | `NOT_VALIDATED` |

The mass intervals are conditional signal-only profile widths with fixed tail
shapes and relative widths. Their low-count coverage remains untested, and the
calibration-shape mismatch flags are retained. The three fitted core widths
share one common scale and fully correlated fractional width uncertainty.

The 1S maps contain 301,323 pretag and 198,265 tagged decays out of 1,126,936
eligible truth events. The tagged efficiency is 17.59%, and the conditional
proton-tag fraction is 65.80%. All 240 pT–rapidity and 288 xV–W cells are exported.
The displayed pT–rapidity range excludes 39 eligible, nine pretag and three
tagged decays. The xV–W grid contains every eligible decay.

The 1S pT² slope is `2.46 (+0.55 / −0.53) GeV⁻²` over 0–1 GeV² and passes the
conditional signal-only tests. The 2S slope has feasible bins but fails
injected-slope coverage. The separate 1S fit targeting true absolute momentum
transfer fails coverage, tail-shape and combined-closure checks. All three W
spectra fail injected-spectrum coverage and shape/combined closure. All three
rapidity spectra fail injected-spectrum coverage, with additional shape and
combined-closure failures for 1S. These diagnostics retain their recorded
statuses in the Atlas.

The separate fixed-bin W unfolding for 2S and 3S is withheld at both energies
because its highest truth-populated bin has zero reconstruction efficiency.
The adaptive-binning W diagnostics above have usable response bins and are
displayed with their failed validation statuses. At 9×130 and 10 fb⁻¹, all three
pT² slopes pass conditional validation. That exposure is ten times the 1 fb⁻¹
used for the note's early-science 9×130 mass projection.

Figure 6 uses raw producer counts from the corrected 1S files at both energies.
Truth stages 1–3 and reconstructed stages 4–10 form separate sequences. The
verified note helper splits the curve at the 3→4 boundary and leaves losses at
stages 1 and 4 undefined. Its exact output PDF and helper hashes are recorded in
the beam provenance manifest. The standard framework renderers at `8e5d380`
implement the same counter definitions. The 54.3% loss at 9×275 from stage 5 to 6 measures the
additional scattered-electron candidate requirement, including central and
backward-tagger candidates. The figure uses each producer counter before the
offline analysis windows and luminosity weighting.
[Figure 6](plots/producer-counters.pdf) · [Raw counts and adjacent losses](data/producer-counters.csv)

The two-energy kinematic panels use inclusive selections at 10 fb⁻¹ per beam.
They apply reconstructed W windows and leave generated W unrestricted. The
integrated mass-study yields additionally apply the generated-W windows.
[Beam figure provenance and selection audits](data/beams-provenance.json)

## Historical stored-matching paper replay

The [historical paper replay](../paper/README.md) starts at the stored-matching merged tuples listed in
the [CERNBox input manifest](../inputs/paper_inputs.csv). For the 9×275 GeV,
2.5 fb⁻¹ proton-tagged sample, the following quantities were reproduced from
the original inputs.

| Quantity | Result |
|---|---|
| Three-state mass fit | Expected yields 66.9, 20.8 and 17.6 events and core widths 58.4, 64.6 and 68.5 MeV for 1S, 2S and 3S. All twelve asymmetric errors agree at the paper's precision. |
| 1S efficiency | 220,772 of 1,126,936 eligible truth events pass before the proton requirement (19.59%). After proton tagging, 144,473 pass (12.82%). The two published maps contain 240 pT–rapidity and 288 xV–W cells. |
| W, rapidity and pT² | Central values and both sampling-band endpoints for all 19 displayed bins agree with the paper exports. The 1S pT² slope and profile endpoints also agree. |
| 1S momentum-transfer slope | The fit to true pT² gives `2.46 (+0.66 / −0.62) GeV⁻²` and passes the conditional validation. The distinct fit targeting true absolute momentum transfer fails tail-shape closure. |

All three W and rapidity spectra retain the overall `NOT_VALIDATED` status.
All three fail injected-spectrum coverage in W and rapidity. Additional
shape-variation checks fail for all three W spectra and for 1S rapidity.
For 2S and 3S neither slope observable has feasible binning. These outcomes
are compared alongside the numerical points. The quoted bands describe a
fixed-response, signal-only projection.

## Checks performed during historical repository assembly

The local checks used Python 3.12. Python source was also parsed against the
Python 3.9 grammar; this syntax check is not a full Python 3.9 runtime test.

| Suite | Result | Scope |
|---|---:|---|
| Repository interface | 32 passed | Workspace isolation and source guards, actual backend planning, dependency imports, plot profiles, production planning/merge bookkeeping and source collection |
| Preserved analysis backend | 32 passed | Recoil-proton selection, unchanged truth normalization, orchestration, category contracts and reoptimization selection |
| Plotting | 81 passed | Saved numerical inputs/statuses, covariance-aware display merges, unavailable cases, mass fits/intervals and reference conventions |
| Efficiency extension | 29 passed | Final selections, topology consistency, truth denominators, chunking, coordinates, maps and count intervals |
| 9×130 extension | 18 passed | Configuration/source guards, synthetic ROOT bookkeeping, support audits and unavailable inference |

The repository checker verifies **120 imported analysis files against their
original snapshot hashes** and **two amended plotting scripts against their
recorded current hashes**. Each amendment records the original baseline hash,
the current hash and the reason for the change in
`provenance/analysis_snapshot.json`. The checker also checks Python syntax,
clean notebook outputs, main documentation links and absence of detector ROOT
files or embedded archives in the source tree.

The mass/reconstruction, W and differential figures were rendered from the
supplied 9×275 proton-tagged numerical exports. PDF/PNG and CSV/JSON companions
were produced. W/rapidity rendering used explicit diagnostic opt-in because the
saved cases are NOT_VALIDATED. These statuses and unavailable higher-state pT²
cases remained in the audits; rendering did not promote them. The differential
figure and the README SVG were visually inspected.

The W notebook's display-rebin switch had a later hard-coded override; this was
removed. The active rebinning module was older than its supplied notebook/tests;
the archived update for saved-edge adjacent merges and proton-contract metadata
was recovered. Origins and before/after hashes are in
`plotting/source_inventory.json` and `plotting/SOURCE_NOTES.md`. No inference
formula was changed for repository assembly.

The standalone four-map configuration recipe and its exact plan command were
tested with explicit temporary paths. It did not read CERN event data.

## Scope of the historical replay and producer checks

- The historical event-level comparison starts with the six stored-matching merged tuples. The
  9×275 files also match the SHA-256 values in the historical paper runs.
- We ran the tuple producer on original Υ(1S) 9×275 and Υ(2S) 9×130 EICrecon inputs. All 49 branches for each set of 100 events matched the corresponding rows in the distributed tuple bit for bit. Each single-input `hadd` preserved the tree and all five histograms. [Input checksums, matched rows and paper selection](missing-production-sources.md#two-reconstructed-inputs)
- For the 2S 9×130 CERN campaign, the input list and worker records agree with all
  10,583 blocks of 100 events in the merged tuple. The comparison uses four
  electron-matching counters per input. [Campaign audit](../production/cern-2026-10/README.md)
- Both beam energies have preserved producer sources. Their raw weights differ
  from the reference for 1S at 9×130 and the excited states at 9×275. Analysis
  normalization is applied after merging. [Per-sample comparison](missing-production-sources.md)
- Notebook 07's W–xV cache and consistency checks require the original 1S tuple.
- Existing conditional or failed statistical validation, signal-only assumptions,
  response limitations and omitted experimental uncertainties are unchanged.

Keep the original analysis runs and their numerical audits. Use fresh workspace
and run names for any subsequent source or physics change.
