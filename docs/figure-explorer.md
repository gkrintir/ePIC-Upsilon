# Rebuild the figure explorer

The [figure explorer](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) follows the figure numbers in analysis note revision `ad6c3b0`. Under **Paper figures**, the sample is proton tagged at 9×275 GeV and 2.5 fb⁻¹. Under **Beam comparison**, the curves show inclusive reconstructed events at 10 fb⁻¹ for each beam energy. The two views have different selection and normalization conventions.

## Paper figures

To inspect an existing run, open the figure and choose the output folder containing `paper_replay.json`. The browser reads the fit intervals and binned tables in that folder. The export below is for a standalone copy.

Follow the [paper replay](../paper/README.md) from the repository root. After the replay succeeds, export its results using the same Python environment.

```bash
python3 tools/export_explorer_paper.py --replay work/paper-replay --out work/explorer
```

Use `work/explorer/data/paper.json` with the local data control on the paper-figure page. The exporter reads the reported mass-fit intervals, efficiency counts, and displayed spectra. The status associated with each fit travels with its values.

## Two-energy kinematics

Use the [analysis quick start](../README.md#analysis-quick-start) to set the input directory. Run the two energies in one inclusive workspace.

```bash
python3 upsilon.py init --name inclusive-note --selection inclusive
python3 upsilon.py analysis run --workspace inclusive-note --energy all \
  --states 1S 2S 3S --lumi 10 --batch-id note-two-energy-L10 --hash-inputs
python3 tools/export_explorer_beams.py \
  --runs workspaces/inclusive-note/results/generated/runs \
  --batch-id note-two-energy-L10 --out work/explorer
```

On the beam-comparison page, choose `work/explorer/data/beams.json` to inspect the bin values from your run. Each state’s source is `artifacts/note/compare_9x130_9x275/kinematic_bins.csv` and its selection audit. Figures 10 and 11 use expected events per bin. Figure 12 uses yield per GeV² and retains the projective angular weight. Empty projective bins can have unavailable interval endpoints.

## Open a local copy

Copy `docs/` into a new `work/explorer/` directory before running the two exports above. Open `figures.html` for the paper figures or `beams.html` for the beam comparison. Each figure links to its analysis code. A selected note PDF opens at the physical page for that figure. Your PDF and imported numerical data are read locally in the browser.
