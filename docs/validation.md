# Repository validation

Run `python3 tools/run_checks.py` from the repository root after installing
`requirements-dev.txt`. Each component runs in its own process and working
directory so historical modules with the same names cannot shadow one another.

## Check against the early-science paper

The [paper replay](../paper/README.md) starts at the merged tuples listed in
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

## Checks performed during repository assembly

The local checks used Python 3.12. Python source was also parsed against the
Python 3.9 grammar; this syntax check is not a full Python 3.9 runtime test.

| Suite | Result | Scope |
|---|---:|---|
| Repository interface | 32 passed | Workspace isolation and source guards, actual backend planning, dependency imports, plot profiles, production planning/merge bookkeeping and source collection |
| Preserved analysis backend | 32 passed | Recoil-proton selection, unchanged truth normalization, orchestration, category contracts and reoptimization selection |
| Plotting | 81 passed | Saved numerical inputs/statuses, covariance-aware display merges, unavailable cases, mass fits/intervals and reference conventions |
| Efficiency extension | 29 passed | Final selections, topology consistency, truth denominators, chunking, coordinates, maps and count intervals |
| 9×130 extension | 18 passed | Configuration/source guards, synthetic ROOT bookkeeping, support audits and unavailable inference |

The repository checker confirms that all **122 files in the current imported
analysis baseline are byte-identical** to their supplied source snapshot. It
also checks Python syntax, clean notebook outputs, main documentation links and
absence of detector ROOT files or embedded archives in the source tree.

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

## Scope of the checks

- The event-level comparison above starts with the six merged tuples. The
  9×275 files also match the SHA-256 values in the historical paper runs.
- Both beam-specific macros compiled and loaded with ROOT ACLiC in separate clean processes under ROOT 6.40.00 and [podio v01-02](https://github.com/AIDASoft/podio/tree/d3f03813d7bf50f87c3aa4f32e619cebeb060113), using the real ObjectID and CollectionIDTable dictionaries. This verifies compilation in that environment. An original EDM4eic event file is needed to exercise the readers and compare their output.
- ROOT 6.40.00 `hadd` merged two archived 18×275 pilot tuples with 1,666 rows each. All 41 tree branches matched the concatenated input arrays, and `h_EventCount` summed to 3,332. Verifying the current 9-GeV-electron production also requires original reconstructed EDM4eic inputs and a pilot in the intended environment.
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
