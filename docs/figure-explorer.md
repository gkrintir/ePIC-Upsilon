# Rebuild the figure explorer

The [figure explorer](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) follows the updated analysis note at revision `f7cb067`. Its default data and PDFs use the full-sample downstream EICrecon #3011 rematching. The mass, efficiency and spectral views show the proton-tagged 9×275 GeV projection at 2.5 fb⁻¹. The beam comparison shows inclusive reconstructed events at 10 fb⁻¹ per energy.

## Current rematched results

The published Atlas contains the corrected fit intervals, all 528 efficiency cells, 27 displayed spectral points, the per-state slope availability and the associated validation statuses. The provenance manifest records all six input hashes and generated-row counts, the checked source artifacts and the exported PDF hashes. [Current numerical data](data/paper.json) · [Provenance and validation](data/paper-provenance.json)

To rebuild these assets from the extracted rematched results archive, point `--results` at the directory containing `SHA256SUMS.txt` and `paper_figures/`.

```bash
python3 tools/export_rematched_atlas.py \
  --results work/rematched-results --out docs
```

The exporter checks each source artifact against the archive manifest before writing `data/paper.json`, its JavaScript companion, `data/paper-provenance.json` and the three PDF/WebP figures. The inference calculation used framework revision `072c1fd`. The published Fig. 6 PDF comes from the verified note helper at revision `f7cb067`. Its PDF and helper hashes are recorded in the beam provenance manifest. Framework revision `8e5d380` implements the same corrected counter interpretation in the standard renderers. [Beam and Fig. 6 provenance](data/beams-provenance.json)

The current note-page links use physical PDF pages. The joint mass panel is on page 27, the efficiency maps on page 46, the W spectrum on page 49 and the differential spectrum on page 50. Table 8 is on page 24. Select the corresponding updated note PDF in the browser to open these pages.

## Historical paper replay

The [historical paper replay](../paper/README.md) writes `work/paper-replay/paper.json` from the earlier stored-matching inputs. Choose that file in the explorer to inspect that run’s fit intervals, efficiency counts and spectral bins. The historical replay checks the original numerical fingerprints and writes `paper_replay.json`. The run-folder importer reads that receipt and checks the table bookkeeping.

For a completed historical run, choose **Open a run folder instead** or create its summary from the saved outputs.

```bash
python3 tools/export_explorer_paper.py --replay work/paper-replay
```

## Two-energy kinematics

Use the [analysis quick start](../README.md#analysis-quick-start) to set the corrected input directory. Run the two energies in one inclusive workspace.

```bash
python3 upsilon.py init --name inclusive-note --selection inclusive
python3 upsilon.py analysis run --workspace inclusive-note --energy all \
  --states 1S 2S 3S --lumi 10 --batch-id note-two-energy-L10 --hash-inputs
python3 tools/export_explorer_beams.py \
  --runs workspaces/inclusive-note/results/generated/runs \
  --batch-id note-two-energy-L10 --out work/explorer
```

On the beam-comparison page, choose `work/explorer/data/beams.json` to inspect the bin values from your run. Each state’s source is `artifacts/note/compare_9x130_9x275/kinematic_bins.csv` and its selection audit. These plots use generated-mass and Q² cuts, the reconstructed dielectron mass and longitudinal balance cuts, and reconstructed W in 20–70 GeV or 20–100 GeV. Generated W is unrestricted in these panels. The note’s integrated mass-study yield table additionally applies the beam-specific generated-W window. [Current beam provenance](data/beams-provenance.json)

Figures 10 and 11 use expected events per bin. Figure 12 uses yield per GeV² and retains the projective angular weight. Empty projective bins can have unavailable interval endpoints. In the updated note, Fig. 10 occupies physical pages 35–37, Fig. 11 pages 38–43 and Fig. 12 page 45.

The beam page also links the corrected Fig. 6 producer counters on physical page 22. Those values are raw 1S simulation counts. Truth stages 1–3 and reconstructed stages 4–10 form separate sequences, so the plot breaks the 3→4 connection and the loss table leaves that transition undefined. [Figure 6 PDF](plots/producer-counters.pdf) · [Counter CSV](data/producer-counters.csv)

## Open a local copy

Copy `docs/` into a new `work/explorer/` directory. Its bundled data already contains the corrected projections. To rebuild the current mass, efficiency and spectral assets from the checked archive, run the following command.

```bash
python3 tools/export_rematched_atlas.py \
  --results work/rematched-results --out work/explorer
```

Open `work/explorer/figures.html` for the tagged projections. For a new beam comparison, run the two-energy export above and open `work/explorer/beams.html`.

To display the historical reference in a separate local copy, use `tools/export_explorer_paper.py --replay work/paper-replay --out work/historical-explorer` after copying the site files there.

After editing the guide, update its JavaScript and CSS query hashes in `docs/index.html` and run both website checks.

```bash
node tools/check_analysis_trace.mjs
node tools/check_atlas.mjs
```

The guide checker verifies all five routes, their source references, local assets and version hashes. The Atlas checker verifies the current rematched data and figure assets, the 528 efficiency cells, 27 spectral rows and 3,602 beam bins, and compatibility with historical imports.
