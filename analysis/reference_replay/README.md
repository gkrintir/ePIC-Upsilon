# Historical W replay

This directory contains actual re-executions on the recovered original ROOT
tuples. All 14 archived CSV comparisons pass, and the recovered low-curve recipe
reproduces the native note's full cross-section figure pixel for pixel. The
efficiency, closure, truth-ratio and flux-ratio panels also match their native
reference pixels exactly. It does not use the current framework's default tuple
aliases or binning.

From this directory, run:

```bash
python3 run_reference_replay.py \
  --low-root /eos/user/g/gkrintir/merged_elecs_test_10t100.root \
  --high-root /eos/user/g/gkrintir/merged_elecs_test_tuples_fix.root \
  --backend auto \
  --out-dir verified_results
```

`--backend root` requires PyROOT and runs the unchanged legacy readers.
`--backend uproot` executes the same frozen scripts with an explicit, narrow
uproot adapter for their ROOT I/O calls. The adapter leaves their physics,
selections, normalization, flux integration, histogramming, uncertainty formulas
and plotting unchanged. Unsupported calls raise errors. `auto` selects PyROOT
when installed. The delivered run used **uproot, not PyROOT**.

Python requirements for that run were NumPy2.3.5, pandas2.2.3,
Matplotlib3.10.8, uproot5.7.7 and mplhep0.4.1. Mplhep enables the legacy ROOT-style
plot appearance and was installed for the final image comparisons. The frozen scripts use `numpy.trapz`; install NumPy<2.4 when
using these historical files. Newer numerical conventions belong in the modern
framework, not in this reference replay.

Input SHA256 hashes are enforced. The original low tuple has44,018 entries;
the expanded high `fix` tuple has58,738. Updated low and `fix2` are different
samples and are deliberately not substituted. Input files are not bundled.

## What is reproduced numerically

| Stage | Exact historical settings | Evidence |
|---|---|---|
| Low1S extraction |10×100GeV,104.793fb visible,14 equal W bins20–50GeV |Archived cross-section and normalization CSVs |
| High1S/2S/3S extraction |18×275GeV,327.584/84.780/68.929fb visible,14 equal W bins20–140GeV |Archived truth-form and corrected-yield CSVs and normalizations |
| Low/high truth comparison |Three equal bins20–26.5GeV,event Q²≤0.001GeV²,state normalization before event cut |Archived comparison and metadata |
| High response and unfolding |Archived rounded edges,seed12345,50% training,4 iterations,300 toys,event Q²≤0.001GeV² |Archived response counts/probabilities,closure summary,normalization |
| High W resolution |Same archived rounded edges and event Q²≤0.001GeV² |Archived resolution summary and selected response counts |

The nominal low/high cross-section extraction has **no optional event Q² cut**;
its flux integration uses Q²max=0.001GeV². This reproduces the archived state
weight sums. Adding the event cut changes the normalization and cannot be called
the same archived calculation. Both source and analysis luminosities are10fb⁻¹.
Visible dielectron weights are divided by the corresponding branching fraction
exactly once in the cross-section conversion.

Historical uncertainties are preserved. They are not the current projected
Poisson intervals or correlation-corrected MC uncertainties. Bin-by-bin closure
uses the same truth/reconstruction sample and is algebraic; the response-matrix
closure uses the original disjoint training/test split. The latter archived
cross sections use the half-sample test yield with the original10fb⁻¹ divisor;
this intentionally reproduces that historical absolute-yield convention.

`verified_results/replay_manifest.json` records the backend, commands, code/input hashes,
independent selection/count/weight checks, output hashes and numerical results.
`verified_results/differences/` contains cell-by-cell archived-versus-rerun differences.
`verified_results/archived_vs_rerun_validation.png` displays cross sections, efficiencies,
relative numerical differences and bin-by-bin closure. The14 archived CSV
comparisons pass; response/unfolding summary and resolution CSVs are byte exact.

## Recovered native note figure and low-curve provenance

`archived_reference/upsilon_selected_overlay_1S_family.png` is the untouched
native note image. The corresponding low CSV named
`sigma_gp_1S_10x100_W20to26p5_Q2lt0p001_unfolded_withEffMC.csv` is stale: it
contains W10–26.5GeV, target304.793fb and zero reconstructed events. It cannot
produce the note's visible low curve. It is preserved as provenance and never
used as the replay's input.

The surviving low CSV was not used as evidence of its own correctness. First,
two deliberately labeled low-curve hypotheses were evaluated:

* `low_selected_candidate_*`: correct legacy tuple and104.793fb, but the stale
  CSV's W10–26.5/10-bin geometry. This produces five supported points, so it does
  not reproduce the note's four-point geometry.
* `low_visual_geometry_hypothesis_W10_25`: correct legacy tuple and104.793fb,
  event Q²≤0.001GeV²,10 bins from10–25GeV. Four supported points at19.75,21.25,
  22.75 and24.25GeV have80.841696,97.561174,112.934713 and129.054022pb. They are
  close to the native note but not pixel identical.

The successful recovered recipe is **W10–25GeV, 10 equal bins,104.793fb,
no optional event Q² cut, and flux Q²max=0.001GeV²**. The selected low points are:

| W center [GeV] | Cross section [pb] |
|---|---:|
|19.75|79.909700|
|21.25|96.436425|
|22.75|111.632727|
|24.25|127.566202|

`verified_results/recovered_note_low_10x100.csv` is the actual output of the frozen
extractor on the recovered original low tuple.
`verified_results/recovered_note_overlay_1S_family.png` combines it with the numerically
verified high family. Its 2850×2250 RGBA pixels are all identical to the native
note image. An individual vector PDF and the five-panel PDF collection are also
provided. The native image check is independent of the 14 archived CSV checks;
it does not pretend that a matching historical low CSV or command was recovered.
The recipe was recovered through the frozen producer and original inputs, then
verified against the surviving native raster.

The high curve, efficiency, closure and truth-ratio contents come from the rerun.
The flux-ratio panel rerenders an archived CSV: no eSTARlight executable was run.
`verified_results/native_image_comparison.json` records all five pixel checks and hashes.
PNG file hashes may differ because of metadata/compression even when decoded
pixels are identical. Other software/font versions can change the raster while
the CSV checks still pass.

The selected-plot renderer is also frozen and keeps its historical labels. Use
this README and manifest when distinguishing archived numerical comparisons,
the recovered recipe's exact raster reproduction, and the exploratory candidate
plots. The first two forms of evidence both pass in the delivered results.

PNG output is serialized to memory, decoded, written and fsynced to a temporary file, then atomically replaced. This fixes intermittent partial PNG writes observed during validation; no frozen plotting content changes. Every final PNG was decoded after all producers exited, and every stored output hash was independently rechecked before delivery.
