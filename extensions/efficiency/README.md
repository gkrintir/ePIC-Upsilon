# Υ efficiency diagnostics at 9×275 GeV

This is an isolated copy of the updated diagnostics package, named `Ypsilon_Efficiency_Diagnostics_ProtonTagged`. Keep the original no-tag diagnostics and notebook folders. It writes to the proton-tagged AFS/EOS framework roots shown below; the numerical selection and plotting implementation are unchanged from the tested proton-selection update.

**Implemented does not mean regenerated:** the explicit proton selection is present in this code, but neither installing the framework category nor copying this package rerenders existing 2D figures. Execute a fresh diagnostic run to create tagged maps. The full framework runner does not invoke these standalone maps. No luminosity, cross section, branching fraction, photon flux or event-weight factor enters their MC efficiencies.

The optional legacy `reference_efficiency_csv` paths in `config_9x275.json` still name the original no-tag campaign. They may be unavailable under the separate tagged output root; that is recorded as UNAVAILABLE, not a failed tag test. If you want those optional before-tag count checks, use a copied JSON configuration with the corresponding reference paths explicitly pointing to the preserved original campaign. Do not point those references at tagged CSVs.


This additive workflow investigates the lower 1S efficiency and its high-pT recovery. It reads existing merged ROOT files, writes a fresh diagnostic run, and does not change previous framework runs, normalization, extraction scripts, or notebooks 01–03. It does not rerun detector simulation.

## Run on lxplus

Place this separate package inside the proton-tagged framework directory, then run from that directory:

```bash
cd /afs/cern.ch/work/g/gkrintir/private/Ypsilon-Fit-Framework-ProtonTagged
python3 Ypsilon_Efficiency_Diagnostics_ProtonTagged/run_efficiency_diagnostics.py plan --framework-root .
python3 Ypsilon_Efficiency_Diagnostics_ProtonTagged/run_efficiency_diagnostics.py validate --framework-root .
python3 Ypsilon_Efficiency_Diagnostics_ProtonTagged/run_efficiency_diagnostics.py run \
  --framework-root . --run-id efficiency-9x275-proton-v3 --hash-inputs
```

Use a new run ID for each execution. Existing run directories are refused rather than overwritten. The workflow continues with other samples if one file is missing or fails. The final status is PARTIAL if an input fails, a saved-reference count comparison mismatches, or plotting fails. Missing optional reference CSVs are recorded as UNAVAILABLE, not claimed to pass. COMPLETE describes execution only; causal physics interpretation remains to be reviewed.

Python 3.9 or later is supported. This uses NumPy, pandas, Matplotlib and uproot, without PyROOT, SciPy or tomllib. Existing framework installations normally have these packages. If necessary, install this package's requirements in your analysis environment. The existing framework requirements.txt is untouched.

## Default inputs

All files end in `merged_Upsilon_final.root` under `/eos/cms/store/group/phys_top/gkrintir/upsilon_data/`.

| Sample key | Directory | Nominal Q² maximum [GeV²] |
|---|---|---:|
| 1s_newid | Upsilon1S/9x275/hiAcc/output_newID | 1 |
| 1s_standard | Upsilon1S/9x275/hiAcc/output | 1 |
| 2s_standard | Upsilon2S/9x275/hiAcc/output | 0.01 |
| 3s_standard | Upsilon3S/9x275/hiAcc/output | 0.01 |

These names are configuration labels, not authentication of the production revision. Edit `config_9x275.json` or pass `--config YOUR_CONFIG.json` to change inputs. `--samples 1s_newid 1s_standard` restricts a run to the two ID scenarios. Whole-file hashes are optional; ordered hashes of the five truth columns are always recorded. Equal truth hashes identify equal stored truth arrays in the same order but do not prove that ID was the only reconstruction change.

## What is counted

The base masks reproduce the supplied All States performance producer, recorded in `selection_provenance.json`. New runs additionally require an explicitly validated reconstructed proton by default.

- Eligible stored truth rows have finite true mass within ±0.2 GeV of the relevant state, nonnegative true Q² up to the configured maximum, and 20≤true W≤100 GeV.
- Selected rows additionally have 8≤reconstructed mass≤12 GeV, 20≤reconstructed W≤100 GeV, |Σ(E−pz)−18|≤1.8 GeV, and (by default) a tagged proton with topology 1/2 and positive finite momentum. These broad performance cuts differ from the narrow mass windows of the cross-section extraction.
- Each 1S sample is evaluated for its nominal Q² range, Q²≤0.01, and 0.01<Q²≤1. The boundary belongs only to the lower band, and low+high counts exactly reproduce all counts.
- The ratios use raw MC row counts, never weight branches, cross sections, branching fractions, luminosity, or fitted mass yields. Their errors are individual 68.27% Wilson intervals. Zero-success bins retain an upper bound; empty-denominator bins are undefined. The threshold of 50 eligible rows is only a visible sparsity flag.
- A rejected reconstructed-mass value can indicate a nonfinite/sentinel entry or a finite value outside the accepted interval. `failure_reasons.csv` and the mass-failure plot separate these cases; a sentinel alone does not establish which upstream reconstruction step failed.
- Neither missing producer rows nor internal tracking, cluster-association, or tagger failures can be reconstructed from these stored scalar branches. The cutflow isolates the three stored performance requirements followed by the explicit proton policy, not upstream detector algorithms.
- No absolute-value repair, renormalization, interpolation, automatic tail-bin merging, ID-difference significance, or independent-sample assumption is made.

## Outputs and reading order

Everything is saved under:

```
results/generated/runs/efficiency-9x275-proton-v3/artifacts/efficiency_diagnostics/
```

Start with `summary.md`, `diagnostic_status.json`, `sample_audits.json` and `input_comparability.json`. `plots/` contains PDF and PNG versions. The principal plots are the new-ID Q² split, matched-Q² ID comparison, Q² mixture fractions/contributions, high-pT counts, selection cutflows, pT–W/y/Q² maps, and tail-composition plots. State comparison panels retain the ID/Q² labels explicitly.

`tables/` contains:

| File | Content |
|---|---|
| efficiency.csv | Exact numerator/denominator counts and Wilson intervals versus true pT, W, and rapidity |
| q2_mixture.csv | Fractions and additive contributions of the two Q² populations; verifies the mixture identity |
| cutflow.csv | Sequential pass/fail counts by pT and Q² band |
| failure_reasons.csv | Nonfinite/sentinel values versus finite lower/upper rejections, separately for each cut and Q²/pT bin |
| nminus1.csv | Counts passing the other two cuts, isolating failure of the omitted requirement |
| maps.csv | pT versus W, rapidity and Q² for each Q² band, with efficiency and population counts |
| composition.csv | Eligible and selected W/rapidity/Q² bin fractions within fixed pT regions; not densities |
| strata.csv | Efficiency at fixed pT, W and rapidity bins, with raw counts |

The all/low/high curves overlap statistically; they are not independent measurements. Separate normalized composition histograms should always be read with their counts and audited out-of-range fractions. The matched-Q² ID plots match analysis cuts; remaining changes in production, row populations, or kinematic composition must be checked before making a causal ID claim.

The three nominal state inputs also cross-check their **pre-proton, inclusive** raw counts against `allstates-projected-L10-v1-{1s,2s,3s}/artifacts/performance/9x275/selection_efficiency.csv` where available. Normalization is immaterial to this check because only nonweight branches are used. A mismatch is recorded and makes the run PARTIAL rather than silently replacing the old result. This legacy reference check does not validate the new tagged numerator.

## SWAN

Open `04_efficiency_diagnostics.ipynb` from this extracted folder. Set FRAMEWORK_ROOT to your EOS framework copy and RUN_ID to a fresh name. The same code and JSON configuration are used. Direct ROOT access must be available in the SWAN session. Outputs go under the chosen framework root, separate from the publication figures.

To copy an lxplus diagnostic run to your existing EOS framework copy:

```bash
src=/afs/cern.ch/work/g/gkrintir/private/Ypsilon-Fit-Framework-ProtonTagged/results/generated/runs
dst=/eos/user/g/gkrintir/Ypsilon-SWAN-ProtonTagged/results/generated/runs
mkdir -p "$dst"
cp -r "$src/efficiency-9x275-proton-v3" "$dst/"
```

## Return the diagnostic results for interpretation

From the framework directory, this archives only the new numerical audits and plots, not any ROOT inputs:

```bash
zip -r efficiency-9x275-proton-v3-review.zip \
  results/generated/runs/efficiency-9x275-proton-v3/effective_config.json \
  results/generated/runs/efficiency-9x275-proton-v3/run_manifest.json \
  results/generated/runs/efficiency-9x275-proton-v3/artifacts/efficiency_diagnostics
```

## Local verification

See VALIDATION.md. Tests use synthetic inputs and independent copies of the original selection predicate. They do not establish the cause of the production-sample trends. Actual CERN EOS inputs must be run before the note is updated with a physics explanation.

## Restyle existing maps in SWAN

Open `05_final_efficiency_2d_display.ipynb`, set `RUN_ID` and `SAMPLE` to an existing completed 2D run, and select **Run All**. This standalone display notebook includes its renderer, reads only `selection_audit.json` and `tables/final_efficiency_2d.csv`, and does not execute the framework or read ROOT. No modules need to be installed or replaced in your existing notebook folder.

It shows each of the four maps in a separate cell with PNG and PDF links. Figures are written below the notebook working directory in `figures/final_efficiency_2d/<RUN_ID>/<SAMPLE>/style-v3-*/`. Each execution creates a fresh export folder. The numerical source run and earlier images remain unchanged.

The subtitle omits the proton category, the in-panel legends are removed, and the colorbar reads `Efficiency [%]`. The proton selection remains required and is checked in both CSV and audit metadata. Gray cells still have no eligible MC. All populated cells use the efficiency color scale without hatching; MC counts and intervals remain in the saved CSV. The notebook retains these explanations and writes source hashes and category metadata beside the figures.

## Four standalone final-selection maps

`run_final_efficiency_2d.py` is a separate entry point. It produces exactly four single-panel figures for one sample, each saved as PDF and PNG. It does not render the diagnostic dashboards. The default sample is `1s_newid`.

```bash
python3 Ypsilon_Efficiency_Diagnostics_ProtonTagged/run_final_efficiency_2d.py validate --framework-root .
python3 Ypsilon_Efficiency_Diagnostics_ProtonTagged/run_final_efficiency_2d.py run \
  --framework-root . --run-id efficiency-2d-9x275-proton-v3 --hash-inputs
```

The four figures are `pt_q2`, `pt_y`, `xv_q2`, and `xb_q2`. All use the full nominal final performance selection described above, including the default proton requirement, without low/high-Q² split or comparison panels. The Q² cap remains sample-specific. Use `--sample 1s_standard`, `--sample 2s_standard`, or `--sample 3s_standard` and a different run ID for another sample.

Outputs are under `results/generated/runs/efficiency-2d-9x275-proton-v3/artifacts/final_efficiency_2d/1s_newid/`. The `plots/` folder contains the four PDF/PNG pairs, `tables/final_efficiency_2d.csv` stores every bin's integer numerator/denominator and Wilson interval, and `selection_audit.json` records selection, formulas and flow counts. The new entry point refuses existing run IDs.

For SWAN, open `05_final_efficiency_2d.ipynb`. Its default framework output location is `/eos/user/g/gkrintir/Ypsilon-SWAN-ProtonTagged`. The ROOT inputs remain the configured `/eos/cms/store/group/phys_top/gkrintir/upsilon_data/...` files.

The two x maps require an event-level ROOT reread. Marginal histograms cannot reconstruct the joint x–Q² distribution. Both axes always use truth values, even for the selected numerator. xB uses Q²/(W²+Q²−mp²), and xV uses (Q²+M_state²)/(W²+Q²−mp²) with a fixed nominal state mass and the small-|t| convention. No pT/rapidity display range is imposed on the x maps. No luminosity, branching fraction, cross section or weight branch enters these efficiencies.

Logarithmic Q² display defaults to 10⁻¹⁰ GeV² through the sample cap. Values below this range, including exactly zero, are audited as outside the display and never moved into the first bin. Undefined efficiencies are gray, while genuine zero-selection bins remain efficiency zero. Sparse bins retain their intervals in the CSV and are displayed without hatching. All four plots use a fixed 0–100% color scale and no interpolation.

Optional binning overrides go inside `settings.final_2d` in a copied JSON config, with keys `pt_edges`, `y_edges`, `q2_edges`, `xb_edges`, `xv_edges`. Log-axis edges must be positive. Binning describes MC efficiency maps and does not establish feasible projected 2D measurement bins.

## Explicit proton category update

Both entry points now default to `--proton-selection required`. The reconstructed numerator must pass `reco_topology in {1,2}` **and** finite `reco_tagged_proton_p > 0`. The truth denominator, generated-event normalization and weights are unchanged. The four-map colorbar reads **Total efficiency [%]**, and the figure identifies the proton category. This is combined acceptance/reconstruction/selection efficiency relative to eligible stored truth rows, not pure geometrical acceptance or a proton efficiency conditional on acceptance.

The audited tuple producer assigns topology 1/2 to valid central/tagger lepton combinations with a proton, 3/4 to the corresponding combinations without a proton, and 0 when no valid combination exists. Topology 0/3/4 therefore requires momentum sentinel -999. Unknown, noninteger or nonfinite topology codes, or inconsistent momentum values, fail the schema check rather than being silently discarded. Validation reads both tag branches for every row.

To reproduce the old proton-inclusive comparison category, pass `--proton-selection inclusive` explicitly. That mode does not require or read either proton branch. Give it a separate fresh run ID, for example `efficiency-9x275-inclusive-recheck-v2`. Old outputs and run IDs are never overwritten. The SWAN notebooks expose `PROTON_SELECTION` explicitly and use new tagged run IDs.

The full diagnostic saves final tagged `efficiency.csv`, pre-proton `inclusive_efficiency.csv`, and `proton_selection.csv` with before/after selected raw MC yields and tag retention by pT and Q² band. Its sequential cutflow has an additional proton stage. Audits record the policy and aggregate counts before/after the proton selection. The four-map CSV also records before/after counts in every cell. No projected luminosity, cross section or event-weight factor is used in these MC yield counts.
