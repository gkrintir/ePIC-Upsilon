# Rebuild the figure explorer

The [figure explorer](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) follows the figure numbers in analysis note revision `ad6c3b0`. Under **Paper figures**, the sample is proton tagged at 9×275 GeV and 2.5 fb⁻¹. Under **Beam comparison**, the curves show inclusive reconstructed events at 10 fb⁻¹ for each beam energy. The two views have different selection and normalization conventions.

## Paper figures

The [paper replay](../paper/README.md) writes `work/paper-replay/paper.json`. Choose that file in the explorer to inspect the fit intervals, efficiency counts and spectral bins. Send the same file to a colleague to share those numbers.

For an older completed run, choose **Open a run folder instead** or create the summary from the saved outputs:

```bash
python3 tools/export_explorer_paper.py --replay work/paper-replay
```

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

Copy `docs/` into a new `work/explorer/` directory. To include your paper figures and numerical results, run:

```bash
python3 tools/export_explorer_paper.py --replay work/paper-replay --out work/explorer
```

Open `work/explorer/figures.html` for the paper figures. For the beam comparison, run the two-energy export above and open `work/explorer/beams.html`.
