# Source migration record

Imported the live source from the supplied `Ypsilon-SWAN-ProtonTagged` snapshot.
Original SHA256 hashes and source/destination paths are recorded in
`source_inventory.json`. The original snapshot remains unchanged.

## Canonical files

The active `upsilon_esr/*.py` modules were used. Nested updater payloads were not
installed over them. In particular, the active differential module retains the
later **ePIC** label in the slope inset, and the active notebook 07 retains its
latest table and `FIGURE_FONT_SCALE=1.35` typography.

`03_differential_available_cases.ipynb` becomes the canonical
`03_differential.ipynb`. It has the cleaner availability diagnostics and fresh
output behavior. The manually executed debugging copy, its stale `.py` companion,
and old updater/backups are not active entry points. Their authoritative source
snapshot and migration hashes provide the bookkeeping record.

The other canonical notebooks are 01, 02, 06 and 07. Cell outputs, execution
counts and session-specific metadata were stripped. Inline figure/layout and map
algorithms in 06/07 were retained to avoid an unverified extraction/refactor.
The standalone binning audit is also retained.

## Changes made in this consolidation

1. Replaced duplicated notebook settings and hardcoded EOS paths with one JSON
   profile loaded by `settings.py`. Source run IDs remain explicit. The example
   preserves actual state-specific Q² and electron-ID settings.
2. Added a CLI profile/check path and fresh-output handling to `render_all.py`.
   The older bin-by-bin audit CLI also reads the shared profile.
3. Fixed notebook 02's rebin switch: no later cell unconditionally overrides a
   disabled merge. Optional merges use only saved edges and preserve covariance,
   joint-toy intervals and the original validation status.
4. Added `Config.allow_diagnostic_results=False` and a strict saved-inference
   rendering gate in `reoptimization_inputs.load_case`. `NOT_VALIDATED` diagnostics
   require explicit opt-in; their status and publication flag remain unchanged.
   The opt-in is recorded in the case metadata. No inference values are changed.
5. Notebook output paths now follow the shared profile and use fresh directories;
   output links are relative to the notebook working directory where possible.
6. Added the dependencies used by notebook execution and the W–xV builder to the
   plotting requirements. No framework dependency file was changed.

The DSCB model, profile intervals, efficiencies, photon-flux conventions,
reoptimization conversion checks, reference numbers, color choices, fonts,
legends, lines, fit labels and saved slope reader remain unchanged.

## Tests recovered and adapted

Recovered the twelve original plotting test files from
`Ypsilon_ProtonTagged_Plotting.zip` under `Ypsilon_ESR_SWAN_ProtonTagged/tests/`.
This includes the previously absent `test_proton_inputs.py` fixture dependency.
The synthetic reoptimization fixture explicitly enables diagnostic rendering,
and the proton-settings fixture explicitly loads the example profile. These are
test-input adaptations to the new explicit policy, not relaxed assertions.

New profile tests cover explicit path/run selection, no fallback, missing profile,
invalid luminosity/selection settings, fresh outputs, disabled and saved-edge
merging, wrong source luminosity/category, explicit diagnostic permission, and
clean executable notebook cells.

## Deliberate limits

These are the supplied **proton-tagged, three-state** styled notebooks. An
inclusive result cannot be made tagged by a label or by disabling metadata
checks. The current layouts also cannot consume a partial one-state/all-unavailable
9×130 reoptimization. Those outcomes have a status/audit path documented in the
plotting guide. Generated PDFs/PNGs, ROOT inputs, caches, compiled Python files,
and old update installers are excluded from the source tree.

The recovered plotting archive also supplies the narrow `rebinning.py` update tested by its suite: automatic `keep_outer_bins` uses actual saved boundaries, reports invalid requested edges, and retains the proton-selection contract in merged CSV rows and merge metadata. Its source hash is recorded separately. Covariance propagation, joint-toy replay, normalization and fit formulas are unchanged.

## Verification performed

All 81 recovered plotting tests and nine new profile tests pass. Command-line preflight and actual reconstruction, W and differential rendering also passed against the supplied 9×275 proton-tagged full-state and reoptimization outputs. Diagnostic rendering was explicitly enabled for this historical-data check; all final validation statuses and unavailable states were retained. The differential PNG was visually inspected for its preserved typography, ePIC slope inset, upper-left fit legend and lack of MC curves/footer.

Notebook code cells were syntax-checked with outputs cleared. The complete 06/07 notebooks were not rerun here: notebook 07 needs access to the original external ROOT tuple. Their live layout and computation cells are retained, with shared input/output configuration substituted. Production data and smoke-test figures are external to this source repository.
