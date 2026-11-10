# ePIC Υ analysis

Follow an exclusive Υ result from reconstructed ePIC events to the analysis note and early-science paper. Each stage shows the code that reads an input and the file it passes to the next stage.

**New to the analysis?** [Start at the tuple producer](https://gkrintir.github.io/ePIC-Upsilon/#mass/reconstruction). If you already have the merged ROOT files, [start at the inputs](https://gkrintir.github.io/ePIC-Upsilon/#mass/tuple) or run the two replay commands below.

**Reprocessing the July samples for EICrecon #3011?** [Submit the 1S 9×275 Condor test and full sample](docs/condor-rematch-1s-9x275.md), [monitor progress](docs/condor-rematch-1s-9x275.md#monitor-the-grouped-campaign) or [retry failed inputs](docs/condor-rematch-1s-9x275.md#retry-failed-or-unfinished-files). The [rematching guide](docs/rematch-3011.md) also covers the interactive EOS setup, the second pilot and the analysis fits.

**Reviewing a figure?** Trace a paper result back to its calculation and input tuples.

- Fig. 18 shows the [three-state mass fit (a)](https://gkrintir.github.io/ePIC-Upsilon/#mass/publication) and [tagged Υ(1S) efficiency maps (b, c)](https://gkrintir.github.io/ePIC-Upsilon/#efficiency/publication).
- [Fig. 19](https://gkrintir.github.io/ePIC-Upsilon/#w/publication) gives the cross section versus W. [Fig. 20](https://gkrintir.github.io/ePIC-Upsilon/#differential/publication) gives the pT² slope and rapidity spectra.
- The analysis note compares the two beam energies in [Figs. 10–12](https://gkrintir.github.io/ePIC-Upsilon/#beams/publication).

## Clone and install

Python 3.11 or newer is recommended. Run the examples from the repository root.

```bash
git clone https://github.com/gkrintir/ePIC-Upsilon.git
cd ePIC-Upsilon
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Reproduce the early-science paper

The paper projects proton-tagged events at 9×275 GeV and 2.5 fb⁻¹ for all three Υ states. The download includes those three merged files and the three 9×130 GeV files for the note’s beam-energy comparison. It preserves the directory structure and verifies each checksum. [CERNBox inputs](https://cernbox.cern.ch/s/UHPvC04xj29cxrp)

```bash
python3 inputs/fetch_paper_inputs.py work/paper-inputs
python3 paper/replay_paper.py --inputs work/paper-inputs --output work/paper-replay
```

If you already have the files, skip the download and point `--inputs` to the directory containing `Upsilon1S`, `Upsilon2S` and `Upsilon3S`. Choose a new `--output` directory for each replay.

The replay makes all three paper figures and checks the numerical results against the paper. Open `work/paper-replay/paper.json` in the [figure explorer](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) to inspect or share the numbers. The same folder contains `paper_replay.json`, which records the checks. [Figure paths and expected values](paper/README.md)

## Analysis quick start

To run your own analysis, first copy the configuration file.

```bash
cp configs/local.example.json configs/local.json
```

In `configs/local.json`, set `input_base` to the absolute path of the directory containing the three state folders. For the download above, use the path printed by `realpath work/paper-inputs`.

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

For the tagged 9×130 GeV projection, use `--energy 9x130 --lumi 1` and new batch and run IDs. In the tagged 1S study at this setting, every binning tested had at least one bin with fewer than five expected signal events, so the analysis could not fit a slope.

For the note’s equal-luminosity comparison, create a separate workspace with `--selection inclusive` (no recoil-proton tag requirement) and use `--lumi 10` for the analysis at each energy. [Analysis commands and output paths](docs/analysis.md)

## Plot saved results

Open the [paper figures](https://gkrintir.github.io/ePIC-Upsilon/figures.html#mass) or the [two-energy comparison](https://gkrintir.github.io/ePIC-Upsilon/beams.html). Each plot links back to its analysis and tuple producer. For Jupyter notebooks, display merging and diagnostic rendering, see [plotting instructions](docs/plotting.md).

## Workflow and reference

[![From reconstructed events or merged UpsilonTree files to reviewed analysis figures](docs/workflow.svg)](https://gkrintir.github.io/ePIC-Upsilon/#mass/reconstruction)

- [Produce tuples at either beam energy](docs/production.md) or [check the original 2S 9×130 campaign against its merged ROOT file](production/cern-2026-10/README.md).
- [Selections and normalization](docs/selections-and-normalization.md) define the electron identification, recoil-proton tag and state-specific Q² limits.
- [Validation](docs/validation.md) records what each result supports. [Bookkeeping](docs/bookkeeping.md) explains how to resume a run.
