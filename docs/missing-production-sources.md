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

## What remains to verify

The ZIP contains processed outputs and generator files. A producer smoke test needs an original reconstructed EDM4eic ROOT file with an `events` tree. For the first 2S 9×130 input in the private CERN list, follow the [one-file pilot](production.md#one-file-pilot). Its worker logged **148 usable, 102 missing, 0 masked and 2 ambiguous matches** over 100 events. A 1S 9×275 reconstructed input and its per-job tuple provide the corresponding direct check for the paper energy.

For the other five state/energy combinations, matching submission records, source or build hashes, and merge manifests are still needed to attribute the final files to specific jobs. The saved 3S submission settings name destinations at both energies; they alone do not identify completed jobs. For new pilots, the [production workflow](production.md) records its source snapshot and explicit merge inputs.
