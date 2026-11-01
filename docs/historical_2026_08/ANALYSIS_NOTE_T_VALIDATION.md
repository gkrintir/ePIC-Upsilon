# Figures 29–36 independent rerun

## Result

The promoted workflows completed successfully without changing any framework file.

- Output directory: /workspace/scratch/80bba5a372e7/tmp/corrected_note_rerun/fig29_36
- Generated files: 50 (5 benchmark products and 45 unfolding products).
- Report figure assets: 23 PNGs.
- Figures 32–36: all 20 panel PNGs are byte-for-byte identical to the assets in the attached analysis note.
- Figures 29–31: corrected reruns are intentionally not identical because the requested high-energy input contains additional events relative to the tuple used for the archived benchmark figures.

## Input resolution

There is no single low-energy ROOT file for all eight figures.

| Figures | ROOT input | SHA-256 | UpsilonTree entries | Fingerprint evidence |
|---|---|---|---:|---|
| 29–31, low-energy half | /workspace/scratch/80bba5a372e7/upload/merged_elecs_test_10t100_10x100.root | 0c67638a69a81c036fc0274bc2f6e54b9e28e442e355518cdbe1edce483550aa | 24,404 | Benchmark normalization gives entries_truth_state=24128, raw 1S cross section 103.60783084740206 fb, scale 1.0114389920424405; archived deterministic low-energy benchmark values are recovered to floating-point roundoff. |
| 33, 34, 36 | /workspace/scratch/80bba5a372e7/upload/merged_elecs_test_10t100.root | f714ae48f8903ffc83d34bddcf2cdcb56d409d26e84a430b2e875c35119658b8 | 44,018 | With 1S ±0.2 GeV, Q²≤0.001, and 20≤W≤26.5 GeV, sumw_state_raw=15165.035505339434, exactly matching all archived low-energy unfolding normalization files. |
| 29–32, 35, high-energy component | /workspace/scratch/80bba5a372e7/upload/merged_elecs_test_tuples_fix.root | d3eab02fc1aa72f3467f36eed31ea2970ce2dc85d545c96cfa81fa680dea7dda | 58,738 | With 1S ±0.2 GeV, Q²≤0.001, and 37≤W≤114 GeV, sumw_state_raw=31399.726286982, matching the archived high-energy unfolding normalization. |

The updated low tuple is required for Figures 29–31, while the legacy low tuple is required for Figures 33, 34, and 36. The framework's current single inputs.low_root wiring cannot recover both archived families without a per-producer input override.

## Code and environment

| Item | Value |
|---|---|
| Benchmark script | scripts/validation/benchmark_t_estimators.py, SHA-256 c08da94f8debe0b79f4e4faed13b609f99f32d465d568043329d0642a5a0cc0b |
| Unfolding script | scripts/unfolding/unfold_t_spectrum.py, SHA-256 771777f678c83088c1faf375d3c249256d9c0043312da8bd1036cfe6fc42b772 |
| Python | 3.12.13 |
| ROOT | 6.40.00 |
| NumPy | 2.3.5 |
| Pandas | 2.2.3 |
| Matplotlib | 3.9.4, matching the Software tag embedded in the archived PNGs |
| mplhep | 0.4.1 |
| Hash seed | PYTHONHASHSEED=0 |

Both promoted scripts are byte-identical to their archived producer versions.

The final environment prefix was:

~~~bash
PYTHONHASHSEED=0 \
MPLCONFIGDIR=/tmp/ypsilon_mpl_394 \
PYTHONPATH=/tmp/ypsilon_ref394:/tmp/ypsilon_deps
~~~

## Commands

Benchmark for Figures 29–31:

~~~bash
python scripts/validation/benchmark_t_estimators.py \
  --low-file /workspace/scratch/80bba5a372e7/upload/merged_elecs_test_10t100_10x100.root \
  --high-file /workspace/scratch/80bba5a372e7/upload/merged_elecs_test_tuples_fix.root \
  --tree UpsilonTree --state 1S --methods eX,BABE,eXBABE,proj,pT2 \
  --low-visible-xsec-fb 104.793 --high-visible-xsec-fb 327.584 \
  --lumi-invfb 10 --state-norm-mode target_visible_xsec \
  --truth-mass-window 0.20 --reco-mass-min 8 --reco-mass-max 12 \
  --q2-cut 0.001 --low-electron-energy 10 --high-electron-energy 18 \
  --low-w-min 10 --low-w-max 26.5 --high-w-min 37 --high-w-max 114 \
  --t-edges 0,0.5,1,1.5 --train-fraction 0.5 --n-iterations 4 --seed 12345 \
  --out-prefix /workspace/scratch/80bba5a372e7/tmp/corrected_note_rerun/fig29_36/benchmark/compare_t_methods_1S_twoE_separateW
~~~

All unfolding invocations used these common flags:

~~~bash
python scripts/unfolding/unfold_t_spectrum.py INPUT \
  --tree UpsilonTree --out-prefix PREFIX --state 1S --method METHOD \
  --t-edges EDGES --electron-energy EE --w-min WMIN --w-max WMAX \
  --q2-cut 0.001 --visible-xsec-fb XSEC --tuple-lumi-invfb 10 --lumi-invfb 10 \
  --weight-convention visible_decay --mass-window-sigma 3 --true-mass-window 0.2 \
  --train-fraction 0.5 --seed 12345 --n-iterations 4 --n-toys 300 \
  --zeus-csv "/workspace/scratch/80bba5a372e7/work/Ypsilon Fit/zeus_desy_11_186_fig4_exclusive_upsilon_step.csv" \
  --zeus-mode scale_to_zeus
~~~

| Figure | INPUT | PREFIX basename | METHOD | EDGES | EE | WMIN–WMAX | XSEC |
|---:|---|---|---|---|---:|---|---:|
| 32 | merged_elecs_test_tuples_fix.root | unfold_t_1S_18x275_eX_zeusBins | eX | 0,0.2,0.5,1,5 | 18 | 37–114 | 327.584 |
| 33 | merged_elecs_test_10t100.root | unfold_t_1S_10x100_eX_zeusBins | pT2 | 0,0.2,0.5,1,5 | 10 | 20–26.5 | 104.793 |
| 34 | merged_elecs_test_10t100.root | unfold_t_1S_10x100_eX_mergedLowT | eX | 0,0.5,1,5 | 10 | 20–26.5 | 104.793 |
| 35 | merged_elecs_test_tuples_fix.root | unfold_t_1S_18x275_pT2_mergedLowT | pT2 | 0,0.5,1,5 | 18 | 37–114 | 327.584 |
| 36 | merged_elecs_test_10t100.root | unfold_t_1S_10x100_pT2_mergedLowT | pT2 | 0,0.5,1,5 | 10 | 20–26.5 | 104.793 |

## Figure verification

Image metrics compare the rerun PNG directly with the corresponding analysis-note asset. AE is the number of pixels with any channel difference.

| Figure | Assets | Status | Exact comparison |
|---:|---:|---|---|
| 29 | 1 | Corrected high-energy rerun; not archival-exact | 36,689 / 7,312,500 changed pixels (0.501730%); normalized RMSE 0.0376574; NCC 0.988793 |
| 30 | 1 | Corrected high-energy rerun; not archival-exact | 120,687 / 11,250,000 changed pixels (1.072773%); normalized RMSE 0.0758628; NCC 0.948777 |
| 31 | 1 | Corrected high-energy rerun; not archival-exact | 105,657 / 11,250,000 changed pixels (0.939173%); normalized RMSE 0.0768133; NCC 0.946272 |
| 32 | 4 | Exact | 4/4 PNGs byte-identical; AE=0 for every panel |
| 33 | 4 | Exact | 4/4 PNGs byte-identical; AE=0 for every panel |
| 34 | 4 | Exact | 4/4 PNGs byte-identical; AE=0 for every panel |
| 35 | 4 | Exact | 4/4 PNGs byte-identical; AE=0 for every panel |
| 36 | 4 | Exact | 4/4 PNGs byte-identical; AE=0 for every panel |

For the 25 unfolding CSV products, 16 are byte-identical. The other 9 differ only in last-bit decimal serialization: overall maximum absolute numeric difference 3.637978807091713e-12, maximum relative difference 9.311538934866566e-16.

| Figure/family | Byte-identical CSVs | Numeric comparison for non-identical CSVs |
|---|---:|---|
| 32: 18x275 eX ZEUS bins | 1/5 | Maximum absolute difference 3.637978807091713e-12; maximum relative difference 9.311538934866566e-16 |
| 33: 10x100 filename eX, actual pT2 ZEUS bins | 5/5 | Exact |
| 34: 10x100 eX merged low-t | 5/5 | Exact |
| 35: 18x275 pT2 merged low-t | 0/5 | Maximum absolute difference 3.637978807091713e-12; maximum relative difference 9.311538934866566e-16 |
| 36: 10x100 pT2 merged low-t | 5/5 | Exact |

An isolated legacy-input probe using the 49,988-entry high-energy tuple reproduced Figures 29 and 30 byte-for-byte. This confirms that their corrected differences are caused by the explicitly requested 58,738-entry merged_elecs_test_tuples_fix.root, not by the promoted script or plotting environment.

## Caveats found

1. Figure 33 is named and captioned as an eX result, but its archived normalization and unfolding summary explicitly record method=pT2. Running pT2 reproduces all four Figure 33 panels byte-for-byte. The figure caption/filename is therefore method-mislabeled.
2. benchmark_t_estimators.py derives its closure split from Python's salted hash((sample_name, method)). Its --seed alone is insufficient for reproducibility. The corrected run sets PYTHONHASHSEED=0; the script should eventually replace the salted hash with a stable method/sample offset.
3. PyROOT emitted only its package-stage warning. All six requested workflow invocations exited with status 0, with no traceback or empty-selection failure.
