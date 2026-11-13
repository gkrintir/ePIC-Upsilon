# Plotting and SWAN

The plotting code reads named analysis products and exports PDF/PNG figures with
CSV/JSON companions. The selected run is explicit. It never searches for the
latest folder or substitutes older numerical results when a file is missing.

## Configure once

From the repository root:

```bash
cp configs/plots_9x275.example.json configs/plots.local.json
```

Edit `configs/plots.local.json` before running. The example points to the archived
9×275 proton-tagged 2.5 fb⁻¹ runs and displays every fitted bin. For the paper,
start from [`plots_9x275.paper.json`](../configs/plots_9x275.paper.json) instead.
The [paper replay](../paper/README.md) sets the paths for a fresh run.

| Setting | Meaning |
|---|---|
| `results_root` | Directory containing named run folders, ending in `results/generated/runs` |
| `output_root` | Parent for new figure directories; relative paths are resolved against the profile file |
| `full_runs` | One explicit full-framework run ID for each state, used by the mass and efficiency figures |
| `reoptimization_run` | One explicit saved inference run, used by W and pT²/rapidity figures |
| `energy`, `projection_lumi_invfb` | Requested beam and exposure; saved inference must match exactly |
| `proton_selection` | `required` for these preserved tagged notebooks |
| `q2_max`, `saved_categories` | Per-state virtuality limit and sample labels. See the verified tuple settings below. |
| `efficiency_2d_run`, `map_sample` | Run and sample key containing saved 2D efficiency counts/audits |
| `root_input_override` | Optional accessible location of the original 1S ROOT tuple for notebook 07; null uses the recorded path |
| `display_rebin` | Optional adjacent-bin display merges for 1S W or rapidity; false retains saved bins |
| `allow_diagnostic_results` | Explicit permission to render `NOT_VALIDATED` saved inference; default false |

All paths can use environment variables. Unresolved variables cause an error.
An alternate profile can be selected using `--profile` or the environment:

```bash
export UPSILON_PLOT_PROFILE="$PWD/configs/plots.local.json"
```

The six paper tuples all record reconstructed electron ID (`ele_id_mode=1`).
`saved_categories` retains the historical run labels, including “standard electron ID”
for 2S/3S, so the plotter can check the source-run contract. Existing profiles
using `categories` are also supported. Figure audits report the verified
reconstruction mode when the original tuple SHA-256 matches the input manifest.
The `output` and `output_newID` directory names locate the inputs. The event selection and flux
integration use Q²max = 1 GeV² for 1S and 0.01 GeV² for the excited states.

## Command-line plots

Use the framework analysis environment or a SWAN kernel with the packages in
`plotting/requirements.txt`. PyROOT and a LaTeX installation are not required for
these plotters. Notebook 07 additionally uses uproot to read the original tuple.

Check the mass/efficiency inputs without writing figures:

```bash
python3 plotting/render_all.py --profile configs/plots.local.json \
  --section reconstruction --check
```

Render that figure:

```bash
python3 plotting/render_all.py --profile configs/plots.local.json \
  --section reconstruction
```

Use `--section w`, `--section differential`, or `--section all` for the other
supported CLI figures. Each execution creates a fresh directory under
`output_root/<section>/render-*`, prints its path, and writes the effective plot
profile. An explicit `--out` must name a directory that does not already exist.

The supplied W/rapidity cases include `NOT_VALIDATED` statuses. Reproducing their
historical diagnostic figures therefore requires an explicit choice:

```bash
python3 plotting/render_all.py --profile configs/plots.local.json \
  --section w --allow-diagnostics
```

This permits diagnostic rendering only. Overall statuses and failed checks stay
in the numerical audit; the option cannot promote a result to a validated
measurement. A `PASS_CONDITIONAL` result still omits the experimental effects
listed in its saved scope. The original figure styling has been preserved,
including the absence of a footer, so the numerical audit must accompany
any diagnostic figure shared outside this workflow.

## Paper Figure 28

The dedicated renderer combines the corrected three-state mass fit in panel (a)
with the saved transverse-momentum spectra and slope inset in panel (b). Install
the dependencies in `plotting/requirements.txt`, then run from the repository root

```bash
python3 tools/render_paper_figure28.py \
  --results-dir paper/figure28_inputs \
  --output-dir work/figure28
```

Choose a new output directory for each run. The checked inputs are included in
`paper/figure28_inputs`, with their source-archive provenance and SHA-256 manifest.
The renderer verifies the saved mass and exponential models against their binned
predictions and uses the saved intervals at 9×275 GeV and 2.5 fb⁻¹.

The output directory contains `upsilon_es_figure28.pdf` and its PNG preview,
`note_mass_panel.pdf`, `note_differential.pdf`, and `figure28_audit.json`. The note
outputs use the same panel-drawing functions. Both paper panels carry
`26.07.1 Campaign`, with common state colors and a shared ePIC label for the
Υ(1S)/Υ(2S) slope comparison. The Υ(2S) slope retains its diagnostic status in the
audit and the figure caption.

## SWAN notebooks

Open the notebooks in `plotting/` with a fresh kernel and run all cells. Every
notebook reads the same local profile; restarting the kernel is required after
changing Python modules. The profile can also be selected explicitly by replacing
`load_profile()` with `load_profile("/absolute/path/to/profile.json")` in its setup
cell. There are no duplicated EOS paths or run IDs in the notebook code.

| Notebook | Contents | Numerical source |
|---|---|---|
| `01_reconstruction.ipynb` | Combined signal mass fit, efficiency versus W and pT | Per-state `artifacts/performance/<energy>/` CSV/JSON and run-level proton contract |
| `02_w_cross_section.ipynb` | Flux-averaged W cross sections and ZEUS reference | Saved reoptimization W cases, covariance, conversion factors and final validation status |
| `03_differential.ipynb` | pT², saved slope inset, rapidity; explicit unavailable states | Saved reoptimization pT²/y cases and the corresponding status/covariance/profile files |
| `06_mass_and_2d_efficiency.ipynb` | Mass plus y–pT and Q²–xV efficiency maps | Full-framework mass inputs and saved 2D map counts/audits |
| `07_mass_and_efficiency_w_xv.ipynb` | Latest mass table with yield/resolution intervals, y–pT and W–xV maps | Same mass/map inputs, plus the original 1S tuple to build/check W–xV |

Notebooks 06/07 retain their tested layout code and are SWAN entry points rather
than `render_all.py` sections. Notebook 07 caches W–xV counts under
`output_root/map_cache/w_xv`, checking source path, size, modification time,
settings and agreement with all four saved 2D maps. The original ROOT path must
remain accessible for that check even when a cache exists.

The mass fit is performed by the plotting package using saved signal templates;
it is not just a copy of a framework PNG. It remains a conditional signal-only
DSCB fit. Quoted resolutions are Gaussian-core widths. The shared relative state
widths and calibration assumptions are preserved. Notebook 07 profiles the same
model for its table uncertainties. Efficiency maps use raw MC counts and do not
scale with luminosity.

The numerical 2D efficiency production is kept separately under
`extensions/efficiency/`; see [the current efficiency guide](efficiency.md) for generating the
counts and selection audit before opening these notebooks.

## Display merging and unavailable results

The optimizer uses training MC at the target luminosity to choose the largest
number of bins passing its expected-yield and response criteria. The held-out
spectrum is fitted with those boundaries. The paper profile applies the
following rule to the fitted 1S spectrum for **both** observables:

```json
"display_rebin": {"1S": {"w": "keep_outer_bins", "y": "keep_outer_bins"}}
```

| 1S observable | Optimized fit edges | Paper display edges |
|---|---|---|
| W [GeV] | 20, 40, 45, 50, 60, 100 | 20, 40, 60, 100 |
| Laboratory rapidity | −1, 0.25, 0.5, 0.75, 1, 3 | −1, 0.25, 1, 3 |

This preserves the first and last saved bins and merges the interior. With
three or fewer bins it leaves the original edges unchanged. Explicit edge lists
are also accepted only when every edge is already present and the full range is
preserved. This is display aggregation, not a new response fit or optimization.
Joint toys and the full covariance determine merged sampling intervals; marginal
error-bar endpoints are never added. The source validation status is unchanged.

`NO_FEASIBLE_SLOPE_BINNING` is a completed numerical outcome. When present in a
supported three-state pT²/rapidity run, notebook 03 lists that state as unavailable
and draws no invented spectrum, slope or fit. Missing files for a supposedly
available case remain errors.

The present layouts require **all three states** and the expected observable
set. They are not a viewer for the new **1S-only 9×130** run whose differential
CSV is entirely empty. Inspect that run's `run_manifest.json`,
`artifacts/reoptimization/validation_status.csv`, `summary.md`, and the matching
`t_support/audit.json` instead. A future partial-state plot adapter needs its own
checks; no fallback to 9×275 or legacy points is used here.

The rapidity panel retains the historical quantity
`(1/Phi_total) d sigma_ep/dy`, with its existing gamma-p shorthand label. It is a
beam/selection-dependent flux-normalized ep spectrum, not a universal fixed-W
gamma-p rapidity distribution. The primary ep spectrum and the displayed
convention are retained in the points CSV and audit.

## Moving numerical outputs to SWAN

Copy the exact named runs, including JSON/CSV and covariance files. PNGs alone
cannot reproduce these figures. For example, with explicit shell variables:

```bash
source_run="$PWD/workspaces/tagged-v1/results/generated/runs/reopt-275-L2p5-v1"
runs_on_eos=/path/to/your/shared/results/generated/runs
mkdir -p "$runs_on_eos"
rsync -av --exclude='*.root' "$source_run/" "$runs_on_eos/$(basename "$source_run")/"
```

Set `results_root` to that EOS directory. Repeat for the exact full-framework and
2D runs selected in the profile. Notebook 07 uses an accessible original tuple,
which is normally already on EOS and is not part of this export.

## Checks

From the repository root:

```bash
python3 -m unittest discover -s tests -p test_plot_profiles.py -v
PYTHONPATH="$PWD/plotting" python3 -m unittest discover -s plotting/tests -v
```

The independent `plotting/audit_projected_binning.py --profile ...` utility audits
older full-framework bin-by-bin CSV bookkeeping. It does not replace the
reoptimization closure/coverage tests or certify the saved slope inference.
