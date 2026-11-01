# Cross-check of the seven supplied ROOT files

## Result

The seven files are structurally valid and are the actual inputs behind most of
the retained dashboard, efficiency, gamma-p, track-correlation, reconstructed-t,
and campaign-comparison artifacts. None of those reruns calls or loads
`DVMP_JPsi_Analysis.C`.

The largest provenance correction is that the two campaign ROOT files contain
ten `TH2D` objects that were overlooked by the retained TH1-only notebook.
Every one of those ten histograms is exactly reconstructible from
`merged_elecs_test_tuples.root`, including every underflow and overflow bin.
The reconstruction rules are now encoded in
`scripts/validation/reconstruct_report_histograms_from_tuple.py`.

## Input identities and full paths used for this check

The absolute paths below are the session-local copies used for validation. The
portable identities are the filename, byte size, and SHA-256 digest.

| Role | Filename and full validation path | Bytes | SHA-256 | Contents |
|---|---|---:|---|---|
| Legacy campaign histograms | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test.root` | 70,079 | `14d350e59d63d10c08bd869e8d7e9d6b33519e92bb44fafa7c05f76a50b98c41` | 45 TH1D + 10 TH2D |
| Updated campaign histograms | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_.root` | 70,700 | `daf90faa374ce88addc660ff078afe66a3946aa07e4fc90926c9c157891f92d5` | 45 TH1D + 10 TH2D |
| Legacy 10x100 tuple | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_10t100.root` | 6,236,639 | `f714ae48f8903ffc83d34bddcf2cdcb56d409d26e84a430b2e875c35119658b8` | `UpsilonTree`, 44,018 rows |
| Updated 10x100 tuple | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_10t100_10x100.root` | 3,999,517 | `0c67638a69a81c036fc0274bc2f6e54b9e28e442e355518cdbe1edce483550aa` | `UpsilonTree`, 24,404 rows |
| Legacy 18x275 tuple | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_tuples.root` | 8,556,381 | `a9a1aa1ebdb32b861c07d612bcb4f8ffadc40a55deef0b9b358497b095777977` | `UpsilonTree`, 49,988 rows |
| Updated 18x275 tuple | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_tuples_fix2.root` | 8,504,431 | `c22591377d31526db1a6e223ef8b086cc9349313a0015cc10540fbd68eec91e3` | `UpsilonTree`, 49,988 rows |
| Expanded 18x275 tuple | `/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_tuples_fix.root` | 10,080,549 | `d3eab02fc1aa72f3467f36eed31ea2970ce2dc85d545c96cfa81fa680dea7dda` | `UpsilonTree`, 58,738 rows + 5 TH1s |

All five trees have the same 41 scalar branches: 36 doubles and five integers.
All file headers, keys, compressed payloads, and recorded end-of-file offsets
were read successfully. No truncation or decompression failure was found.

## Relationships between the files

- The campaign files have identical 55-histogram schemas but different data.
  Only two empty rapidity-bin payloads are binary-identical.
- The 18x275 tuples have identical kinematics row by row. Only `weight` and
  nonzero `weight_proj` change. The updated file is therefore a reweighted
  version, not an independent event sample.
- `merged_elecs_test_tuples_fix.root` contains every one of the 49,988 legacy
  tuple rows bit-for-bit across all 41 branches plus 8,750 inserted rows. Of
  those inserted rows, 3,958 repeat an existing event's non-weight data with
  changed weights, while 4,792 are new non-weight event records. It is the
  recovered source for the W-resolution/W-unfolding family, but it is not an
  interchangeable substitute for `fix2`.
- The updated 10x100 tuple contains 24,404 rows, all of which occur exactly in
  the 44,018-row legacy tuple when the two weight branches are ignored. It is a
  subset with a different weight convention. Raw yields from those two files
  must not be compared without normalization.
- `merged_elecs_test_.root` is the histogram companion of the high-energy
  tuples for the recovered unweighted observables. Its electron-ID histograms
  are also bin-identical to the corresponding objects stored in both
  high-energy tuple files.

## Exact artifact reruns

“Exact” below means byte-for-byte equality unless a row explicitly says
pixel-exact. Exact image reproduction is evidence of input and code provenance,
not evidence that a fit is statistically valid.

| Producer | Uploaded input(s) | Retained output checked | Result | Qualification |
|---|---|---|---|---|
| Notebook 03, cell 0 | `merged_elecs_test_tuples.root` | `Upsilon_Inclusive_Dashboard.png` | Exact | Historical high-energy dashboard input identified |
| Notebook 03, cell 1 | `merged_elecs_test_tuples_fix2.root` | `Upsilon_Inclusive_Dashboard_raw_visibleEE.png` | Exact | Updated/reweighted dashboard input identified |
| Notebook 03, cell 7 | `merged_elecs_test_10t100_10x100.root` + `merged_elecs_test_tuples_fix2.root` | `Upsilon_TrueLevel_Overlay_Dashboard.png` | Exact | Both energy inputs identified |
| Notebook 02, cell 0 | `merged_elecs_test_tuples.root` | 1S, 2S, and 3S `Kinematic_Efficiencies-Copy1.png` | All exact | Files without `-Copy1` are different stale variants |
| Notebook 06, cell 0 | `merged_elecs_test_tuples.root` | `Upsilon_Family_GammaP.png` | Exact | ROOT emits a harmless log-axis warning |
| Notebook 01, cell 1 | `merged_elecs_test_10t100_10x100.root` | `Fit_Inclusive_Upsilon_DSCB_Nice.png` | Exact after suppressing propagation of a ROOT warning | Fit reaches call limit and has an unphysical 2S yield; do not treat it as valid |
| Notebook 01, cell 2 | `merged_elecs_test_tuples_fix2.root` | `Fit_Inclusive_Upsilon_RapidityBins_raw_visible_decay.png` | Exact | Both fitted covariance matrices are non-positive-definite |
| Track-correlation script, low | `merged_elecs_test_10t100.root` | PNG/PDF and two CSVs | CSVs exact; PNG pixels exact | Stored PNG differs only in Matplotlib version metadata |
| Track-correlation script, high | `merged_elecs_test_tuples_fix2.root` | PNG/PDF and two CSVs | CSVs exact; PNG pixels exact | Stored PNG differs only in Matplotlib version metadata |
| Reconstructed-t script | `merged_elecs_test_10t100_10x100.root` + `merged_elecs_test_tuples_fix2.root` | PNG/PDF/CSV | PNG pixels and event counts exact | Weighted sums differ only at floating-point formatting level, below `9e-13` |
| W-resolution script | `merged_elecs_test_tuples_fix.root` | resolution summary and response CSVs | Exact | The added file resolves the analysis-note Figures 26-27 input; current PNG differences are Matplotlib raster drift only |
| Family gamma-p normalization | `merged_elecs_test_10t100.root` + `merged_elecs_test_tuples_fix.root` | low/high normalization CSVs | Exact | Analysis-note Figures 15 and 22-23 central inputs recovered |
| W unfolding | `merged_elecs_test_tuples_fix.root` | W response count and probability CSVs | Exact | Analysis-note Figures 24-25 response data recovered; toy errors can drift with the numerical stack |

The two campaign histogram files also regenerate all 43 nonempty archived TH1
comparison PNGs byte-for-byte when the historical raw-count labels and
Matplotlib 3.10.8 are used. The retained Notebook 07 cell has drifted: it now
normalizes both histograms, changes the legend labels, ignores TH2 objects, and
reports 45 outputs although two empty histograms are skipped. The promoted
`scripts/plotting/compare_report_campaign_histograms.py` restores raw-count
rendering, fixes the counter, and adds all ten TH2 comparisons.

## Exact tuple-to-TH2 reconstruction

All fills use unit weight. Neither `weight` nor `weight_proj` is used.

Common selections:

```python
truth = true_mass > 0
reco = (
    (reco_mass >= 8.0)
    & (reco_mass <= 12.0)
    & (abs(sigma_epz - 36.0) <= 3.6)
)
tagged = reco & (reco_t_BABE > 0)
```

| Stored TH2 key | Selection and tuple fill | Binning | Visible / flow fills | Result |
|---|---|---|---:|---|
| `GenLevel/h2_Gen_DecayLeps_Eta_vs_Eta` | truth; x=`true_decay_eMinus_eta`, y=`true_ePlus_eta` | 100x100, x/y `[-5,5]` | 48,892 / 100 | Exact |
| `GenLevel/h2_Gen_DecayLeps_pT_vs_pT` | truth; x=`true_decay_eMinus_pT`, y=`true_ePlus_pT` | 100x100, x/y `[0,12]` | 48,992 / 0 | Exact |
| `GenLevel/h2_Gen_DecayLeps_Eta_vs_pT` | truth; fill each lepton with x=pT, y=eta | 100x100, x `[0,10]`, y `[-5,5]` | 97,884 / 100 | Exact |
| `GenLevel/h2_Gen_W_vs_Rapidity` | truth; x=`true_y`, y=`true_W` | 100x100, x `[-5,5]`, y `[0,150]` | 48,992 / 0 | Exact |
| `GenLevel/h2_Gen_t_vs_W` | truth; x=`true_W`, y=`true_t` | 100x100, x `[10,150]`, y `[0,4]` | 43,870 / 5,122 | Exact |
| `h2_nTracks_Central_vs_Tagger` | all rows; x=`reco_nTracks_Central`, y=`reco_nTracks_Tagger` | 15x15, x/y `[0,15]` | 49,987 / 1 | Exact |
| `h2_eta_ePlus_eMinus` | reco; x=`reco_decay_eMinus_eta`, y=`reco_ePlus_eta` | 100x100, x/y `[-5,5]` | 17,262 / 0 | Exact |
| `h2_pT_ePlus_eMinus` | reco; x=`reco_decay_eMinus_pT`, y=`reco_ePlus_pT` | 100x100, x/y `[0,12]` | 17,262 / 0 | Exact |
| `h2_Reco_DecayLeps_Eta_vs_pT` | reco; fill each lepton with x=pT, y=eta | 100x100, x `[0,10]`, y `[-5,5]` | 34,524 / 0 | Exact |
| `h2_t_Correlation` | tagged; x=`reco_t_eX`, y=`reco_t_eXBABE` | 50x50, x/y `[0,2]` | 7,992 / 3,620 | Exact |

The generated and reconstructed eta-versus-pT maps have identical axes. The
framework therefore emits a binwise `reco/truth` efficiency candidate. The two
input maps and the ratio arithmetic are reproducible, but the precise
report-era convention for empty denominator bins was not retained, so this
derived panel is not claimed as a pixel-exact report reproduction.

## Report-map corrections enabled by the files

| Report figure | Previous status | Corrected source status |
|---|---|---|
| Original p. 6, W vs rapidity | Plot marked missing | Exact TH2 payload survives and its tuple fill is recovered exactly |
| Original p. 6, generated t vs W | Only later dashboard listed | Exact TH2 payload and tuple fill both survive |
| Original p. 7, lepton 2D maps | Later tuple dashboard only | All three relevant generator-level TH2 fills are exact |
| Original p. 9, truth/reco lepton maps | Reconstructed map marked absent | Both truth and reconstructed eta-versus-pT maps survive and are tuple-reconstructible |
| Original p. 12, W-vs-rapidity inset | Plot marked missing | Same exact `h2_Gen_W_vs_Rapidity` payload is available |
| Follow-up p. 4, raw track map | Only selected/weighted analogue listed | Exact all-event raw TH2 payload and tuple fill are available |

## Remaining limitations

- The macro-to-tuple conversion is still historical provenance, but it is not
  a runtime dependency. The supplied `UpsilonTree` branches and the encoded
  NumPy selections are sufficient for the recovered figure families.
- Running Notebook 01 cell 0 on either small campaign file does not reproduce
  the archived scaled inclusive fit. The archived fit used another historical
  state of `h_invMass_Reco_All`.
- Reproducible fit images remain statistically unsafe where ROOT reports call
  limits or non-positive-definite covariance matrices.
- Report-slide image equality and source-data equality are separate claims.
  Several source families are exact while the visible report layout is a crop,
  composite, or later styling variant.

## Automated commands

```bash
python workflow/run_analysis.py --config config/local.toml plan report_crosscheck
python workflow/run_analysis.py --config config/local.toml run report_crosscheck --hash-inputs
```

For direct use without PyROOT:

```bash
python scripts/validation/reconstruct_report_histograms_from_tuple.py \
  --tuple-root /path/to/merged_elecs_test_tuples.root \
  --reference-root /path/to/merged_elecs_test_.root \
  --out-json results/generated/tuple_th2_validation.json \
  --out-npz results/generated/reconstructed_report_th2.npz

python scripts/plotting/compare_report_campaign_histograms.py \
  --old-root /path/to/merged_elecs_test.root \
  --new-root /path/to/merged_elecs_test_.root \
  --out-dir results/generated/campaign_comparison
```
