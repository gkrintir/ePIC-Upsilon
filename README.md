# ePIC Υ analysis

Follow an exclusive Υ result from reconstructed ePIC events to the analysis note and early-science figure layouts. The Atlas now shows the full-sample EICrecon #3011 rematched projections for all six state and energy samples. Each stage links the calculation to its inputs and numerical results.

**New to the analysis?** [Start at the tuple producer](https://gkrintir.github.io/ePIC-Upsilon/#mass/reconstruction). If you already have the corrected merged ROOT files, [start at the inputs](https://gkrintir.github.io/ePIC-Upsilon/#mass/tuple) and follow the analysis quick start below. The separate historical replay reproduces the earlier stored-matching reference.

**Reprocessing the July samples for EICrecon #3011?** [Start an independent 1S 9×275 run](docs/condor-rematch-1s-9x275.md#start-an-independent-eos-run) or [submit the remaining five samples in parallel](docs/rematch-3011.md#submit-the-remaining-five-samples-in-parallel). To reuse finished AFS tuples, [switch the earlier campaign to EOS](docs/condor-rematch-1s-9x275.md#switch-an-existing-afs-campaign-to-eos). The walkthrough also covers [monitoring](docs/condor-rematch-1s-9x275.md#monitor-the-grouped-campaign) and [retrying failed inputs](docs/condor-rematch-1s-9x275.md#retry-failed-or-unfinished-files). The [rematching guide](docs/rematch-3011.md) covers the interactive setup, the second pilot and the analysis fits.

**Reviewing a figure?** Trace the current rematched projection back to its calculation and input tuples.

- Fig. 18 shows the [three-state mass fit (a)](https://gkrintir.github.io/ePIC-Upsilon/#mass/publication) and [tagged Υ(1S) efficiency maps (b, c)](https://gkrintir.github.io/ePIC-Upsilon/#efficiency/publication).
- [Fig. 19](https://gkrintir.github.io/ePIC-Upsilon/#w/publication) gives the cross section versus W. [Fig. 20](https://gkrintir.github.io/ePIC-Upsilon/#differential/publication) gives the pT² slope and rapidity spectra.
- The analysis note compares the two beam energies in [Figs. 10–12](https://gkrintir.github.io/ePIC-Upsilon/#beams/publication).
- [Fig. 6](docs/plots/producer-counters.pdf) shows the raw producer counters with separate truth and reconstructed sequences. Its [CSV](docs/data/producer-counters.csv) leaves the undefined transition between stages 3 and 4 blank.

## Current rematched projections

The 9×275 GeV, 2.5 fb⁻¹ proton-tagged sample contains 91.6, 27.3 and 22.7 expected selected events for 1S, 2S and 3S. The simultaneous mass fit gives 91.8, 27.1 and 22.7 events because it fits the combined calibrated line shapes. Its intervals describe conditional signal statistics. The corrected 1S tagged efficiency is 198,265 / 1,126,936 = 17.59%. [Numerical results](docs/data/paper.json)

The 1S pT² slope is `2.46 (+0.55 / −0.53) GeV⁻²` with conditional validation. The 2S slope is available as a diagnostic with failed coverage, and 3S has no feasible slope binning. The adaptive W and rapidity spectra retain their failed validation statuses. The input hashes, normalization counts, fit assumptions and per-observable checks accompany the figures. [Full provenance and validation](docs/data/paper-provenance.json)

## Clone and install

Python 3.11 or newer is recommended. Run the examples from the repository root.

```bash
git clone https://github.com/gkrintir/ePIC-Upsilon.git
cd ePIC-Upsilon
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Replay the historical paper reference

This replay uses the earlier stored-matching inputs and their original numerical fingerprints. It projects proton-tagged events at 9×275 GeV and 2.5 fb⁻¹ for all three Υ states. The historical download includes those three merged files and the three 9×130 GeV files used for the earlier beam-energy comparison. It preserves the directory structure and verifies each checksum. [Historical CERNBox inputs](https://cernbox.cern.ch/s/UHPvC04xj29cxrp)

```bash
python3 inputs/fetch_paper_inputs.py work/paper-inputs
python3 paper/replay_paper.py --inputs work/paper-inputs --output work/paper-replay
```

If you already have the files, skip the download and point `--inputs` to the directory containing `Upsilon1S`, `Upsilon2S` and `Upsilon3S`. Choose a new `--output` directory for each replay.

The replay makes all three historical paper figures and checks their numerical results against the stored-matching reference. Open `work/paper-replay/paper.json` in the [figure explorer](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) to inspect or share that historical run. The same folder contains `paper_replay.json`, which records the checks. The current rematched figures load by default in the Atlas. [Historical figure paths and expected values](paper/README.md)

## Analysis quick start

To run your own analysis, first copy the configuration file.

```bash
cp configs/local.example.json configs/local.json
```

In `configs/local.json`, set `input_base` to the absolute path of the directory containing the three state folders. Use the corrected merged files for a rematched analysis and compare their hashes with the current input manifest. The historical download above provides the stored-matching samples. [Current input manifest](docs/data/paper-provenance.json)

```bash
python3 upsilon.py init --name tagged-v1 --selection proton-tagged
python3 upsilon.py doctor --workspace tagged-v1
python3 upsilon.py analysis run --workspace tagged-v1 --energy 9x275 \
  --states 1S 2S 3S --lumi 2.5 --batch-id tagged-275-L2p5 --hash-inputs
```

Use the same batch ID to optimize the bins and fit the spectra. The `status` command shows which observables passed validation.

```bash
python3 upsilon.py reopt run --workspace tagged-v1 --energy 9x275 \
  --states 1S 2S 3S --observables pt2 abs_t_via_pt2 w y --lumi 2.5 \
  --source-batch tagged-275-L2p5 --run-id reopt-275-L2p5 --hash-inputs
python3 upsilon.py status --workspace tagged-v1 --run-id reopt-275-L2p5
```

For the note’s early-science tagged 9×130 GeV projection, use `--energy 9x130 --lumi 1` and new batch and run IDs. The full rematched adaptive-inference comparison at this energy was run at 10 fb⁻¹, where all three pT² slopes pass conditional validation. Set `--lumi 10` and use a separate run name to reproduce that exposure. [Per-run validation records](docs/data/paper-provenance.json)

For the note’s equal-luminosity comparison, create a separate workspace with `--selection inclusive` (no recoil-proton tag requirement) and use `--lumi 10` for the analysis at each energy. [Analysis commands and output paths](docs/analysis.md)

## Plot saved results

Open the [paper figures](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) or the [two-energy comparison](https://gkrintir.github.io/ePIC-Upsilon/beams.html). Each plot links back to its analysis and tuple producer. For Jupyter notebooks, display merging and diagnostic rendering, see [plotting instructions](docs/plotting.md).

## Workflow and reference

[![From reconstructed events or merged UpsilonTree files to reviewed analysis figures](docs/workflow.svg)](https://gkrintir.github.io/ePIC-Upsilon/#mass/reconstruction)

- [Produce tuples at either beam energy](docs/production.md) or [check the original 2S 9×130 campaign against its merged ROOT file](production/cern-2026-10/README.md).
- [Selections and normalization](docs/selections-and-normalization.md) define the electron identification, recoil-proton tag and state-specific Q² limits.
- [Validation](docs/validation.md) records what each result supports. [Bookkeeping](docs/bookkeeping.md) explains how to resume a run.
