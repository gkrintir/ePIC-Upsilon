# Four final-selection 2D efficiency maps

This is a separate calculation; the full analysis and reoptimization commands
do not generate these maps. It reads an existing merged ROOT tuple and produces
four individual figures, each in PNG and PDF, plus numerical counts and an audit.

| Output name | Horizontal coordinate | Vertical coordinate |
|---|---|---|
| `pt_q2` | True Υ pT | True Q² |
| `pt_y` | True Υ pT | True Υ rapidity |
| `xv_q2` | Derived true xV | True Q² |
| `xb_q2` | Derived true xB | True Q² |

Both numerator and denominator use truth coordinates. The W–xV map in plotting
notebook 07 is an additional calculation, checked against these four saved maps;
it does not replace one of them.

## 1. Make a local configuration with explicit inputs

Run from the repository root, in the analysis Python environment. The examples
use the `tagged-v1` workspace already created in the analysis quick start. If
`workspace_root` was changed in `configs/local.json`, replace
`workspaces/tagged-v1` below with the absolute workspace path printed at creation.

Keep this local configuration outside the copied analysis workspace's source
directories: adding files there would invalidate its source snapshot. The
following copies the preserved defaults into an ignored local directory, applies
the same `input_base` as `upsilon.py`, and removes historical optional CSV
references. It refuses to overwrite an existing local configuration.

```bash
python3 - <<'PY'
import json
from pathlib import Path
from tools.workspace import REPO, read_settings

settings = read_settings(REPO)  # configs/local.json, or local.example.json
base = Path(settings['input_base'])
config = json.loads((REPO / 'extensions/efficiency/config_9x275.json').read_text())
for sample in config['samples']:
    directory = 'output_newID' if sample['key'] == '1s_newid' else 'output'
    sample['path'] = str(base / f"Upsilon{sample['state']}" / '9x275'
                         / 'hiAcc' / directory / 'merged_Upsilon_final.root')
    sample.pop('reference_efficiency_csv', None)
config['settings']['proton_selection'] = 'required'
destination = REPO / 'workspaces/local-configs/efficiency-9x275-tagged.json'
destination.parent.mkdir(parents=True, exist_ok=True)
with destination.open('x') as stream:
    json.dump(config, stream, indent=2)
    stream.write('\n')
print(destination)
for sample in config['samples']:
    print(sample['key'], sample['path'])
PY
```

If the analysis uses an alternate local settings file, pass its path as the
second argument to `read_settings(REPO, "/absolute/path/to/local.json")` above.
Inspect the printed paths before proceeding. You may edit `path` in this copied
JSON to select another actual input; a scenario label alone does not change the
file. The extension does not expand environment variables inside JSON paths, so
keep those paths absolute.

With the supplied input base, these are the selected directories, all ending in
`merged_Upsilon_final.root`:

| Sample key | Directory below `input_base` | Q²max [GeV²] |
|---|---|---:|
| `1s_newid` | `Upsilon1S/9x275/hiAcc/output_newID` | 1 |
| `1s_standard` | `Upsilon1S/9x275/hiAcc/output` | 1 |
| `2s_standard` | `Upsilon2S/9x275/hiAcc/output` | 0.01 |
| `3s_standard` | `Upsilon3S/9x275/hiAcc/output` | 0.01 |

The `*_standard` and `1s_newid` sample keys are retained for archived run
lookup. The six [paper tuples](../inputs/paper_inputs.csv) at both energies
record `ele_id_mode=1` and `femc_mode=0`. Only three of the four rows above
belong to this set. The separate `1s_standard` comparison has its own input.
When supplied, the producer-mode branches are recorded in the input audit.

These settings describe 9×275. A 9×130 study needs its own configuration with
matching input paths, beam energies, W range and display edges; changing only
the directory name is insufficient.

## 2. Plan, validate, then run one sample

```bash
python3 extensions/efficiency/run_final_efficiency_2d.py plan \
  --framework-root workspaces/tagged-v1 \
  --config workspaces/local-configs/efficiency-9x275-tagged.json \
  --sample 1s_newid --proton-selection required

python3 extensions/efficiency/run_final_efficiency_2d.py validate \
  --framework-root workspaces/tagged-v1 \
  --config workspaces/local-configs/efficiency-9x275-tagged.json \
  --sample 1s_newid --proton-selection required

python3 extensions/efficiency/run_final_efficiency_2d.py run \
  --framework-root workspaces/tagged-v1 \
  --config workspaces/local-configs/efficiency-9x275-tagged.json \
  --sample 1s_newid --proton-selection required \
  --run-id efficiency-2d-275-1s-tagged-v1 --hash-inputs
```

`plan` prints inputs without reading ROOT. `validate` checks required branches
and scans all stored proton topology/momentum pairs for consistency. `run` reads
the tuple and writes a fresh run; existing run directories are refused. Use a
different `--sample` and fresh run ID for each additional sample. `--hash-inputs`
also hashes the complete ROOT file; ordered truth-column hashes are always saved.

The proton choice is **explicit at this standalone entry point**. Its command
line overrides the copied configuration; it does not inherit the selection of
`upsilon.py`. For the no-proton-requirement category, use
`--proton-selection inclusive`, the separate `inclusive-v1` workspace, a separate
copied configuration and fresh run ID. The current combined plotting notebooks
accept tagged maps only. The calculation adds no luminosity option because
these are unweighted MC efficiencies, independent of projected exposure.

## 3. Read the outputs and connect the plotting profile

For the command above, the directory is:

```text
workspaces/tagged-v1/results/generated/runs/efficiency-2d-275-1s-tagged-v1/
```

It contains `effective_config.json`, `run_manifest.json` and
`artifacts/final_efficiency_2d/1s_newid/`, with:

- `plots/pt_q2.{png,pdf}`, `pt_y.{png,pdf}`, `xv_q2.{png,pdf}` and
  `xb_q2.{png,pdf}`: four separate figures, with `Efficiency [%]` colorbars and
  no in-panel legend or hatching.
- `tables/final_efficiency_2d.csv`: integer denominators, counts before/after
  the proton requirement, efficiencies and Wilson intervals for every cell.
- `selection_audit.json`: actual input, cuts, tag policy, formulas, hashes and
  all out-of-range/invalid-coordinate counts.

In `configs/plots.local.json`, set `efficiency_2d_run` to
`efficiency-2d-275-1s-tagged-v1` and `map_sample` to `1s_newid`. Set `results_root`
to the absolute `results/generated/runs` directory containing this run **and**
the full-analysis/reoptimization runs named in that profile. When exporting to
SWAN, copy this exact run alongside those other named runs and use that common
EOS results root; the figure PNGs alone are not sufficient. See
[plotting and SWAN](plotting.md).

The standalone `extensions/efficiency/05_final_efficiency_2d_display.ipynb` can
also redraw the four saved tagged maps separately. It retains its original
configuration cell: explicitly set `FRAMEWORK_ROOT`, `RUN_ID` and `SAMPLE` there.
Unlike notebooks in `plotting/`, it does **not** read the shared plotting profile.
Its setup paths are historical defaults, not a recommendation to use old runs.

## What this efficiency means

The denominator is eligible **stored tuple rows**: finite true mass within
±0.2 GeV of the state mass, 0≤true Q²≤the sample cap and 20≤true W≤100 GeV.
The numerator additionally requires reconstructed mass 8–12 GeV, reconstructed
W in the same interval and `abs(sigma_epz - 2*Ee) <= 0.2*Ee` (16.2–19.8 GeV
for Ee=9 GeV). With `required`, it also requires `reco_topology in {1,2}` and
finite positive `reco_tagged_proton_p`; topology/momentum inconsistencies fail.

This is the final **performance** selection, with a broad mass window. It is
not the narrower mass selection used in differential extraction, nor a pure
geometrical acceptance or proton-only efficiency. Stored rows cannot recover
events that a producer discarded upstream. No weight branch, luminosity,
branching fraction, assigned cross section or photon flux enters these ratios.

The x coordinates are derived as
`xB = Q2 / (W**2 + Q2 - mp**2)` and
`xV = (Q2 + M_state**2) / (W**2 + Q2 - mp**2)`, using nominal state mass.
xV neglects t; it is not measured Bjorken x. Nonpositive logarithmic coordinates
and values outside display ranges are audited, not clipped into visible bins.
Empty-denominator cells are gray; populated zero-efficiency cells remain zero.
The 68.27% Wilson intervals describe finite MC counts, not projected data errors.
These display bins do not establish feasible bins for a projected 2D measurement.

The CLI and copied-configuration plan were checked during repository assembly.
The subsequent [paper replay](../paper/README.md) runs the producer on the original 9×275 1S tuple and checks the two published maps cell by cell.
