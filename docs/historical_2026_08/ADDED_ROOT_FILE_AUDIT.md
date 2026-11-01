# Inventory of `merged_elecs_test_tuples_fix.root`

## Bottom line

`/workspace/scratch/80bba5a372e7/upload/merged_elecs_test_tuples_fix.root` is readable and structurally complete. It is **not** another copy of either `merged_elecs_test_tuples.root` or `merged_elecs_test_tuples_fix2.root`.

- Size: **10,080,549 bytes**
- SHA-256: **`d3eab02fc1aa72f3467f36eed31ea2970ce2dc85d545c96cfa81fa680dea7dda`**
- ROOT UUID: `2b5f0cec-56dd-11f1-ac2d-6b69b8bcbeef`
- Top-level inventory: one `TTree` and five `TH1D` objects, all at cycle 1
- `UpsilonTree`: **58,738 entries**, **41 scalar branches**
- Integrity check: uproot opened every object and read every branch through the final entry without an error. No truncation or corrupt key/basket was detected.

The important content difference is that this file contains the complete 49,988-event `merged_elecs_test_tuples.root` sample plus 8,750 additional rows inserted at the ends of the three state blocks. The original rows remain bit-for-bit identical across all 41 branches when aligned. The added rows are 3,750 1S, 3,750 2S, and 1,250 3S entries.

## File comparison

| File | Bytes | SHA-256 | Tree entries | Object inventory |
|---|---:|---|---:|---|
| `merged_elecs_test_tuples_fix.root` | 10,080,549 | `d3eab02fc1aa72f3467f36eed31ea2970ce2dc85d545c96cfa81fa680dea7dda` | 58,738 | `UpsilonTree` + 5 `TH1D` |
| `merged_elecs_test_tuples_fix2.root` | 8,504,431 | `c22591377d31526db1a6e223ef8b086cc9349313a0015cc10540fbd68eec91e3` | 49,988 | same names/classes |
| `merged_elecs_test_tuples.root` | 8,556,381 | `a9a1aa1ebdb32b861c07d612bcb4f8ffadc40a55deef0b9b358497b095777977` | 49,988 | same names/classes |
| `merged_elecs_test_10t100.root` | 6,236,639 | `f714ae48f8903ffc83d34bddcf2cdcb56d409d26e84a430b2e875c35119658b8` | 44,018 | same names/classes |
| `merged_elecs_test_10t100_10x100.root` | 3,999,517 | `0c67638a69a81c036fc0274bc2f6e54b9e28e442e355518cdbe1edce483550aa` | 24,404 | same names/classes |
| `merged_elecs_test.root` | 70,079 | `14d350e59d63d10c08bd869e8d7e9d6b33519e92bb44fafa7c05f76a50b98c41` | none | 45 `TH1D` + 10 `TH2D` in 2 directories |
| `merged_elecs_test_.root` | 70,700 | `daf90faa374ce88addc660ff078afe66a3946aa07e4fc90926c9c157891f92d5` | none | 45 `TH1D` + 10 `TH2D` in 2 directories |

The new file is 1,576,118 bytes (18.53%) larger than `fix2` and has 8,750 more entries (17.50%).

## Keys, classes, and branch schema

Top-level keys in the new file:

| Key | Class | Main content |
|---|---|---|
| `UpsilonTree` | `TTree` | 58,738 entries, 41 scalar branches |
| `h_EventCount` | `TH1D` | one bin; value 58,738 |
| `h_CutFlow` | `TH1D` | 10 bins; counts `58738, 57742, 55855, 44677, 34619, 20362, 20344, 20266, 20266, 13685` |
| `h_EleID_EpRatio` | `TH1D` | 100 bins on [0, 2], 145,960 fills |
| `h_EleID_IsoFraction` | `TH1D` | 100 bins on [0, 1.1], 145,960 fills |
| `h_EleID_CutFlow` | `TH1D` | 5 bins; counts `145960, 145960, 89362, 89337, 89337` |

The schema is identical, including branch order and types, to all four earlier tree-bearing files. There are 36 `double` branches and 5 `int32_t` branches:

```text
double: weight, weight_proj, true_mass, true_y, true_pT, true_W, true_Q2,
        true_t, true_ePlus_pT, true_ePlus_eta, true_decay_eMinus_pT,
        true_decay_eMinus_eta, true_scat_eMinus_pT, true_scat_eMinus_eta,
        true_scat_eMinus_E, reco_mass, reco_y, reco_pT, reco_W, reco_Q2,
        reco_t_eX, reco_t_BABE, reco_t_eXBABE, reco_t_n, reco_t_x,
        reco_omega, reco_t_proj, reco_ePlus_pT, reco_ePlus_eta,
        reco_decay_eMinus_pT, reco_decay_eMinus_eta, reco_scat_eMinus_pT,
        reco_scat_eMinus_eta, reco_scat_eMinus_E, reco_tagged_proton_p,
        sigma_epz

int32_t: reco_topology, true_nTracks_Central, true_nTracks_Tagger,
         reco_nTracks_Central, reco_nTracks_Tagger
```

## Exact relationship among `tuples`, `fix`, and `fix2`

The state-block layout is:

| State block | `tuples.root` / `fix2.root` | `fix.root` | Added by `fix.root` | Added weight per row |
|---|---:|---:|---:|---:|
| 1S | 19,996 | 23,746 | 3,750 | 2.620672 |
| 2S | 19,996 | 23,746 | 3,750 | 0.678240 |
| 3S block | 9,996 | 11,246 | 1,250 | 0.551432 |
| Total | 49,988 | 58,738 | 8,750 | |

All 49,988 original rows are present in `fix.root` and match `tuples.root` bit-for-bit across all 41 branches. They occur as the original 1S block, then the original 2S block after the 1S addition, then the original 3S block after the 2S addition.

The added rows require caution. Hashing the 39 non-weight branches shows:

- 4,792 added rows have new event content.
- 3,958 added rows duplicate an existing event's full non-weight record but carry the added sample's different `weight`/`weight_proj` values: 1,916 in 1S, 1,931 in 2S, and 111 in 3S.
- No two rows within the added 8,750-row block duplicate one another.

This overlap may be intentional sample combination, but it can also double-count identical events. It should be resolved before treating the raw sum of weights as a physical normalization.

`fix2.root` has the same 49,988 events and every non-weight branch is bit-identical to `tuples.root`. Its only tree changes are:

- `weight`: changed on all 49,988 entries to state-normalized values.
- `weight_proj`: changed on the 2,644 nonzero entries.

Consequently, `fix.root` is not a superset of `fix2` in the weighted sense. `fix.root` preserves the legacy weights and adds events; `fix2` keeps the legacy event set but replaces its weights. The histogram objects in `fix2` are identical to `tuples.root` except that the final two `h_CutFlow` bins are each larger by 2.

Useful normalization fingerprints for the new file are:

| Selection | Entries | Raw sum of `weight` |
|---|---:|---:|
| 1S truth window | 23,746 | 45,861.7600000000 |
| 2S truth window | 23,746 | 11,869.2000000000 |
| 3S truth window | 10,250 | 4,412.945462184874 |
| 1S, `true_Q2 <= 0.001` | 23,737 | 45,844.72248593518 |
| 1S, `true_Q2 <= 0.001`, `37 <= true_W <= 114` | 16,261 | 31,399.726286982004 |

These values exactly identify several archived outputs as having used this file.

## Code references

Direct, active canonical references:

- `notebooks/11_t_binning_optimization.ipynb`, cell 0: the active default is `~/merged_elecs_test_tuples_fix.root`. It produces the resolution-aware t-binning calibration plots, migrations, scans, and recommendation CSVs under `t_binning_resolution_aware_v4/`.
- `notebooks/09_gamma_p_cross_section_two_energies.ipynb`, cell 2: directly opens the EOS path for a Q2/weight diagnostic. This notebook is exploratory and is not the canonical cross-section extractor.

Non-active or parameterized references:

- `notebooks/02_efficiency_and_two_energy_studies.ipynb`, cell 3, mentions `fix.root` only as a commented alternative; the active high-energy input is `fix2.root`.
- `scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py`, `scripts/unfolding/unfold_gamma_p_cross_section_vs_w.py`, and `scripts/validation/validate_w_resolution.py` show `fix.root` in usage examples, but the actual input is a command-line argument.
- Archived counterparts explicitly mentioning it are `extract_sigma_gamma_p_family_fixed_swan.py` and `_v2.py`, `unfold_sigma_gp_W_swan.py`, `study_W_resolution_swan.py`, `unfold_t_distribution_swan.py`, `_v2.py`, and `_v3.py`, plus `MassFit-Copy10.ipynb` and `MassFit-Copy12.ipynb`. `MassFit-Copy2.ipynb` contains only the commented alternative.
- `config/pipeline.toml` does not bind this filename. `inputs.high_root` is blank and must be set explicitly. Therefore a pipeline run does not automatically select the new upload merely because it exists.

Archived metadata and normalization tables prove that `fix.root` was used for the `compare_10x100_18x275_truth_*`, `sigma_gp_family_18x275_fixed*`, `unfold_1S_18x275*`, high-energy t-unfolding, hybrid-t, and `Wres_1S_18x275_full*` output families. Their recorded entry counts and weight sums match the fingerprints above exactly.

## Report-figure impact

For the original 17-page **“Y meson photo-production in ePIC”** report, no requested panel in Figures 3–18 or 20–27 has a confirmed direct dependency on `merged_elecs_test_tuples_fix.root`. Those panels trace to the two histogram-only files, `merged_elecs_test_tuples.root`, `merged_elecs_test_tuples_fix2.root`, or the low-energy tuples. In particular, the current source for the rapidity fit uses `fix2`, not `fix`. The new file should therefore not be silently substituted when cross-checking those panels; its extra 8,750 rows and legacy weights would change the result.

If “report” instead means the later 65-page analysis note, the file has confirmed or strong downstream provenance for:

- Figure 15: projected gamma-p family overlay, through `sigma_gp_family_18x275_fixed*`.
- Figure 16: pT-squared / impact-parameter comparison, through the high-energy t-unfolding output whose raw normalization is exactly 31,399.726286982.
- Figure 21: low/high truth ratio; its metadata names `merged_elecs_test_tuples_fix.root` directly and records 58,738 entries.
- Figures 22–23: W-dependent efficiency and closure, through the same high-energy gamma-p family tables.
- Figures 24–25: W unfolding closure and response matrix, through `unfold_1S_18x275*`.
- Figures 26–27: W-resolution summaries and residual distributions, through `Wres_1S_18x275_full*`.

Analysis-note Figures 12–14, 17–18 use other inputs in the retained sources/metadata, chiefly `fix2.root`; Figure 20 is a photon-flux calculation rather than a tuple-derived plot.
