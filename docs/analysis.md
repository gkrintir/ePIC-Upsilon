# Analysis workspaces and saved results

The repository keeps one copy of the recovered analysis source in `analysis/`.
Run it through `upsilon.py`. The launcher creates separate runtime copies for
an explicit reconstructed selection, checks their source hashes, and passes the
selection to every subprocess. It does not change an existing AFS framework or
resume its historical runs.

## One local configuration

Copy `configs/local.example.json` to `configs/local.json`. The analysis launcher
uses these keys:

```json
{
  "workspace_root": "workspaces",
  "input_base": "/eos/cms/store/group/phys_top/gkrintir/upsilon_data"
}
```

A relative workspace root is interpreted relative to this repository. The input
base must be an absolute path. You can select a separate JSON file with
`--config /path/to/local.json`, before or after the subcommand. The local file is
not a physics-parameter override: generated rates, branching fractions and Q²
settings remain those of the preserved backend.

## 1. Create an analysis workspace

```bash
python3 upsilon.py init --name tagged-v1 --selection proton-tagged
python3 upsilon.py doctor --workspace tagged-v1
```

For the inclusive reconstructed category, use a different name:

```bash
python3 upsilon.py init --name inclusive-v1 --selection inclusive
```

Here **inclusive means no explicit recoil-proton requirement**. It does not mean
inclusive production or inclusive weights. Both categories retain every
generated row for normalization. Tagged reconstruction adds
`reco_topology in {1,2}` and finite positive `reco_tagged_proton_p`.

Initialization copies source and reference files only. ROOT samples, old results,
Python caches and backup files are excluded. Existing destinations are refused.
`init --dry-run` checks the copy plan without creating directories. A failed
initialization retains an explicit incomplete marker for inspection; use a fresh
name rather than overwriting it.

Use `upsilon.py` for normal work. Direct backend invocation requires its own
explicit selection environment and is an advanced interface; the original
backend defaults to requiring a proton.

## 2. Inspect and validate the full state analysis

```bash
python3 upsilon.py analysis plan \
  --workspace tagged-v1 --energy 9x275 --states 1S 2S 3S \
  --lumi 2.5 --batch-id allstates-275-L2p5-v1

python3 upsilon.py analysis validate \
  --workspace tagged-v1 --energy 9x275 --states 1S 2S 3S \
  --lumi 2.5 --batch-id allstates-275-L2p5-v1
```

`plan` lists the exact backend commands. `validate` checks prerequisites and
input availability; it is not an assertion of statistical validation. The
launcher prints the actual command and selection, and `--dry-run` prints without
executing it.

The full workflow validates the pure-state input, prepares a separately
normalized all-row ROOT file, constructs performance and resolution diagnostics,
extracts and unfolds spectra, renders note plots and audits the normalization.
Independent tasks continue after a failure. Dependent tasks remain blocked and
are recorded in the manifest.

## 3. Execute

```bash
python3 upsilon.py analysis run \
  --workspace tagged-v1 --energy 9x275 --states 1S 2S 3S \
  --lumi 2.5 --batch-id allstates-275-L2p5-v1 --hash-inputs
```

The target luminosity is always explicit. Source tuple normalization remains
10 fb⁻¹. The default statistical mode is `projected`; `--stat-mode mc` selects the
MC-statistics mode in a separately named batch. Use `--energy 9x130` and the
intended luminosity, e.g. `--lumi 1`, for the lower-energy study. `--energy all`
uses the same supplied luminosity for both energies, so separate calls are
appropriate when the two projections have different exposures.

Current input paths are constructed as
`<input_base>/Upsilon<state>/<energy>/hiAcc/<directory>/merged_Upsilon_final.root`.
The backend selects `output_newID` for 1S and `output` for 2S/3S.

| State | Q²max [GeV²] | B_ee | σvis 9×130 [fb] | σvis 9×275 [fb] |
|---|---:|---:|---:|---:|
| 1S | 1 | 0.0238 | 119.685 | 210.136 |
| 2S | 0.01 | 0.0191 | 29.297 | 53.133 |
| 3S | 0.01 | 0.0218 | 23.447 | 42.969 |

These settings come from the preserved executable normalization contract. A
common-Q² comparison would require corresponding samples and assigned rates;
changing a plot label or note paragraph does not change this configuration.

## 4. Reoptimize using a prepared source batch

For the 9×275 study:

```bash
python3 upsilon.py reopt plan \
  --workspace tagged-v1 --energy 9x275 --states 1S 2S 3S \
  --observables pt2 abs_t_via_pt2 w y --lumi 2.5 \
  --source-batch allstates-275-L2p5-v1 --run-id reopt-275-L2p5-v1

python3 upsilon.py reopt run \
  --workspace tagged-v1 --energy 9x275 --states 1S 2S 3S \
  --observables pt2 abs_t_via_pt2 w y --lumi 2.5 \
  --source-batch allstates-275-L2p5-v1 --run-id reopt-275-L2p5-v1 \
  --hash-inputs
```

The launcher derives normalized input paths from this workspace and source
batch. It checks the recorded state, beam, normalization, Q² and category. A
PARTIAL source run is usable when input preparation succeeded. The 9×275
configuration is generated in a temporary location with explicit source runs,
then preserved by the backend as the new run's `effective_config.json`. No old
AFS input defaults are used.

For 9×130, first make a source batch with `analysis ... --energy 9x130`:

```bash
python3 upsilon.py reopt run \
  --workspace tagged-v1 --energy 9x130 --states 1S \
  --observables pt2 abs_t_via_pt2 --lumi 1 \
  --source-batch allstates-130-L1-v1 --run-id reopt-130-L1-v1 \
  --hash-inputs
```

This uses the external 9×130 extension. That extension currently supports the
proton-tagged category only. It retains the 0–1 GeV² truth fit interval for the
slope observables, with outside-range feed-in and reconstructed control bins.
It does not recover an unsupported 1–5 GeV² true-|t| tail.

`pt2` targets true pT²; `abs_t_via_pt2` targets true |t| using reconstructed pT².
W and rapidity use spectrum inference. Rapidity's primary result is inclusive ep
dσ/dy; its separate full-flux-normalized companion is a different convention.

The optimization chooses bins using training MC and the projected yield. It
then performs held-out inference and coverage, response-MC variation and
separate shape-only, tail-only and combined stress tests. The count and
migration-quality requirements are unchanged.

## 5. Inspect status before selecting figures

```bash
python3 upsilon.py status --workspace tagged-v1
python3 upsilon.py status --workspace tagged-v1 --run-id reopt-130-L1-v1
```

Outputs are located below the selected workspace:

- Full analysis: `results/generated/runs/<batch-id>-1s/artifacts/`, and similarly
  for 2S and 3S.
- Batch record: `results/generated/batches/<batch-id>/batch_manifest.json`.
- Reoptimization: `results/generated/runs/<run-id>/artifacts/reoptimization/`.

An execution status of COMPLETE is separate from physics validation. Read
`validation_status.csv` and the per-case results. `NO_FEASIBLE_SLOPE_BINNING`
means no selected binning and no slope inference. Header-only numerical tables
are expected for unavailable cases. `NOT_VALIDATED` diagnostics remain labeled
as such. `PASS_CONDITIONAL` does not include all experimental uncertainties.
Plotting must retain these statuses rather than falling back to older results.

## Continuing work and changing code

Repeat the same full-analysis command and batch ID to use its existing stage
cache. Configuration and source guards must match. For reoptimization, repeat
the exact command with `--resume`; changes to inputs, settings or code require a
fresh run ID. The launcher also checks its workspace source snapshot before
execution.

After updating repository code, initialize a new workspace name. Existing
workspaces and results remain available for read-only status and plotting. Do
not copy new Python files into an old production workspace: source hashes are
part of its run contracts. No automatic migration or relabeling of historical
results is performed.

## Legacy tasks and environment

The complete current all-state path is the main interface. Historical
`analysis/config/pipeline.toml` and `pipeline_w_validation.toml` are retained for
bookkeeping. Some historical notebook/replay assets were not present in the
supplied framework snapshot; consult the repository inventory before expecting
those particular tasks to run.

Current all-state analysis uses NumPy, SciPy, pandas, matplotlib, mplhep and
uproot. Python3.9/3.10 also need tomli. ROOT/PyROOT comes from the CERN
ROOT/LCG/SWAN environment for production macros and legacy tasks. `doctor` actually imports each required dependency in a fresh Python process,
including `scipy.linalg` and `scipy.optimize` to catch missing or unreadable
BLAS libraries. It also reports optional PyROOT import and ROOT executable
availability; missing PyROOT alone does not fail the current Python analysis
check. It does not submit jobs or read ROOT event data.
