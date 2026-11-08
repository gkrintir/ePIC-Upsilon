# Which tuple production can we trace?

The six [reference ROOT files](../inputs/paper_inputs.csv) are the inputs to the analysis-note comparison. The paper replay uses the three 9×275 files. Their checksums identify the exact data products. To trace how they were made, we compared the CERN working copies and job records with the contents of those files.

## Υ(2S) at 9×130

The [production record](../production/cern-2026-10/README.md) identifies 76 Condor jobs. They read 10,583 distinct reconstructed files and logged one 100-event tuple per input, using reconstructed electron ID and the full forward calorimeter. All 1,058,300 events are accounted for in the checksum-verified merged file.

For each input, the worker logged four counts of track–cluster matches: usable, missing, masked and ambiguous. Summing the corresponding branches in consecutive 100-row blocks reproduces **all 10,583 logged sets of counts in output-filename order**. Some inputs have identical counts. This check establishes the observed ordering of the diagnostics; the original merge command and per-output checksums would establish exact file membership independently.

The [recovered 9×130 macro](../production/macros/9x130/DVMP_Upsilon_Analysis.C) defaults to the electron-ID mode used by those jobs. The Condor worker supplied only the input and output names. Its job description requested `eic_xl:nightly`; the archived records identify the image selector but give no immutable container digest. The [9×275 macro](../production/macros/DVMP_Upsilon_Analysis.C) is a separate preserved source. The CERN copies of `ElectronFinder.C` and both style helpers match the repository copies byte-for-byte.

## What the six merged files tell us

All rows carry `ele_id_mode=1` and `femc_mode=0`. We recalculated generated photon–proton energy as `W² = (P + k − k′)²` using each beam setting and the saved scattered electron. Across 6,773,621 rows with valid generated kinematics, the largest absolute difference from the stored W was below 1.3 × 10⁻¹² GeV. The raw `weight` branch gives another, independent check:

| State | 9×130 raw weight per event | 9×275 raw weight per event |
|---|---:|---:|
| Υ(1S) | 32.7584 | 32.7584 |
| Υ(2S) | 2.9297 | 5.3133 |
| Υ(3S) | 2.3447 | 4.2969 |

For a 100-event input, the 9×130 working copy predicts the observed 2S and 3S weights but predicts 11.9685 for 1S. The 9×275 source predicts the observed 1S weight but gives 8.4780 and 6.8929 for 2S and 3S. A matching beam constant alone cannot identify the historical source revision. For the paper and note, [normalization](selections-and-normalization.md) uses each state's visible cross section and the total number of generated rows after merging.

Run the beam and weight check on the six downloaded reference files:

```bash
python3 inputs/fetch_paper_inputs.py work/paper-inputs --check-only
python3 production/cern-2026-10/check_beams.py work/paper-inputs
```

## One input identified from an older 1S pilot

The archive also contained a 100-row 1S output labelled `9x130`. Its six generated quantities (meson mass, pT and rapidity; scattered-electron energy, pT and η) match the final 1S 9×130 file exactly at zero-based entries 228300–228399. This identifies the same generated events. The pilot's stored W follows 18×275 beam settings; the final file uses 9×130. The pilot lacks the electron-ID mode branch and has zero tagged candidates in those rows, compared with seven in the final file. Use the final file for the analysis.

## Two reconstructed inputs

We reran one EICrecon file per beam energy through the [production planner and worker](../production/jobs.py), using reconstructed electron ID and the full fEMC, then merged each output with `hadd`. The runs used ROOT 6.38/04 (C++20) and [podio v01-07](https://github.com/AIDASoft/podio/tree/v01-07), matching the ROOT and podio writer versions recorded in the input files.

| Input | Rows in the distributed tuple (zero-based, inclusive) | Electron-ID tracks (usable / missing / masked-only / ambiguous) | Events with a reconstructed proton tag |
|---|---:|---:|---:|
| Υ(1S), 9×275 | 0–99 | 151 / 86 / 0 / 2 | 9 |
| Υ(2S), 9×130 | 136100–136199 | 148 / 102 / 0 / 2 | 7 |

Each input yielded 100 rows. The [comparison](../production/compare_pilot.py) used six generated observables to find a unique block of 100 consecutive events in each distributed tuple. **All 49 branches matched bit for bit for every event**, including the reconstructed kinematics, proton tag, electron ID and raw weight. Each single-file merge preserved all five histograms. The 2S block also starts where the archived job outputs predict in filename order, and its four track counts reproduce the original log.

The 9×275 input shows how tuple rows contribute to the paper. All 100 generated Υ(1S) decays enter the efficiency denominator. Of these, 13 pass the lepton selection and 8 also have a reconstructed proton tag. Applying the [paper efficiency calculation](../extensions/efficiency/diagnostics/final_maps.py) to the fresh tuple and to the matching reference rows gave identical counts in every rapidity–pT bin. The nine tags in the table count candidates before the paper selection. [Full-sample efficiency maps](https://gkrintir.github.io/ePIC-Upsilon/figures.html#efficiency)

To repeat the check, use the [one-file commands](production.md#one-file-pilot) with the matching merged reference from the [input manifest](../inputs/paper_inputs.csv). The SHA-256 digests identify the two original reconstructed inputs:

| Input | EICrecon file SHA-256 |
|---|---|
| Υ(1S), 9×275 | `f2701e1d496d931b117afd5760efe15c909061668890e2d67d6fd28d29856908` |
| Υ(2S), 9×130 | `b5f8450deba55cc4555532f11f87be1d81713ef3cefc51774f84aa7b756a01ff` |

## What remains to verify

The two pilots verify one original input per beam energy. Only the 2S 9×130 campaign also has a check covering every 100-event block. For the other four state/energy combinations, an original EICrecon input would enable the same event-level comparison. Identifying every historical input and the exact source revision used for every job would require further submission or merge records.
