# Reproduce the 9×275 GeV paper figures

From the repository root, install the Python [analysis environment](../README.md#clone-and-install), then run:

```bash
python3 inputs/fetch_paper_inputs.py work/paper-inputs
python3 paper/replay_paper.py --inputs work/paper-inputs --output work/paper-replay
```

The [shared CERNBox folder](https://cernbox.cern.ch/s/UHPvC04xj29cxrp) holds the three Υ states at both energies. If you already have the tuples on EOS, use `--inputs /path/to/your/upsilon_data` and skip the download. The paper replay verifies and reads the three 9×275 inputs. The full input manifest also includes the three 9×130 inputs used for the note's beam comparison. Choose a fresh output directory for each run.

## What the replay does

1. Verify the raw SHA-256 values and prepare separate state tuples with the all-generated visible dielectron weights. The generation and photon-flux upper limits are 1 GeV² for 1S and 0.01 GeV² for 2S and 3S. All six reference files record `ele_id_mode=1` for reconstructed electron identification.
2. Recalculate the per-state mass histograms and calibrations with the recoil-proton tag at 2.5 fb⁻¹. Run the 2D efficiency producer and the original SWAN notebook 07. The notebook fits the combined expected mass spectrum and rereads the 1S ROOT file to compute the W–xV map. The replay checks the three yields, three widths and their twelve asymmetric errors against the archived fit at the paper's printed precision. It checks the numerator, pretag count and denominator in every efficiency cell.
3. Rerun the repository's training-sample optimization and held-out inference. For the 1S W and rapidity displays, the [paper profile](../configs/plots_9x275.paper.json) keeps the outer bins and combines the three interior fit bins. The uncertainty bands are recalculated from the joint toys. The replay compares the central values and both endpoints for all 19 displayed spectral bins, plus the slope and its two profile limits, with the reference exports.
4. Preserve the saved statistical status when comparing points. The paper's W and rapidity plots require the explicit diagnostic rendering option. They retain their failed coverage or closure checks. The 1S pT² slope passes conditional validation. [Checks by observable](../docs/validation.md)

## Where to look

Paths below are relative to `--output`.

| Output | Contents |
|---|---|
| `paper_replay.json` | Overall replay results and the path of the combined mass/efficiency PDF |
| `reconstruction/efficiency_cells.csv` | Observed denominators, pretag counts and tagged numerators for all 528 map cells |
| `spectra/w/upsilon_es_w_cross_section.pdf` | State-resolved W spectra |
| `spectra/differential/upsilon_es_differential.pdf` | 1S pT² fit with slope and all three rapidity spectra |
| `spectra/displayed_points.csv` | Rerun bin edges, central values, endpoints and validation statuses |
| `results/generated/runs/paper-reopt/artifacts/reoptimization/` | Native optimization and inference results before display merging |

The tagged mass projection gives 66.9, 20.8 and 17.6 fitted events, with core widths of 58.4, 64.6 and 68.5 MeV. The 1S map includes 144,473 tagged events out of 1,126,936 eligible truth events. This gives **12.82%**, compared with **19.59%** before proton tagging. The fitted 1S pT² slope is **2.46 (+0.66 / −0.62) GeV⁻²**.

The [numerical fingerprints](reference/checksums.json) come from the archived `allstates-proton-projected-L2p5-v2-{1s,2s,3s}` mass runs and `reopt-proton-9x275-L2p5-v1` differential run. The efficiency cell counts were checked against the paper figure. The fit fingerprint uses the paper’s printed precision; the spectrum fingerprint uses seven significant digits and includes the validation status. The map fingerprint checks the edges and exact integer counts for all three fields in every cell. The new figures are calculated from the merged tuples.
