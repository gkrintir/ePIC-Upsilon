# 9×130 proton-tagged support audit and reoptimization

This external add-on diagnoses the failed fixed-bin tail and runs the existing
reoptimization machinery with explicit 9×130 normalization and photon flux.
No installer is needed. Keep this directory **beside** the framework, not inside
`workflow/`, `scripts/`, or `tools/`. Existing code and results remain unchanged;
adding Python files inside those folders would affect the framework's source
contracts.

## Run the 1S check and optimization

Download `Ypsilon_9x130_Reoptimization_Update.zip` to EOS, then on lxplus:

```bash
cd /afs/cern.ch/work/g/gkrintir/private
unzip /eos/user/g/gkrintir/Ypsilon_9x130_Reoptimization_Update.zip

bash Ypsilon_9x130_Reoptimization_Update/run_1s.sh \
  Ypsilon-Fit-Framework-ProtonTagged \
  allstates-proton-9x130-projected-L2p5-20261101T164736Z
```

The script performs three steps in order:

1. Validate the normalized input, source-run provenance, proton selection and rates.
2. Audit the old `[0,0.2,0.5,1,5]` truth-|t| response before and after proton tagging.
3. Optimize and validate `pt2` and `abs_t_via_pt2`, with a new timestamped run ID.

Use the same Python environment as the framework (NumPy, SciPy, pandas,
matplotlib, uproot). This add-on does not install or change dependencies.

The default target is **1 fb⁻¹**, matching the actual `--projection-lumi-invfb 1.0`
in the failed command. The source run's `L2p5` name does not determine the
luminosity. The prepared input stays normalized to its original **10 fb⁻¹**;
reoptimization scales yields by target/source luminosity. For an intentionally
separate 2.5 fb⁻¹ projection, append `2.5` to the shell command. Each invocation
uses a fresh run ID.

A PARTIAL source workflow is acceptable if `prepare_9x130` succeeded and its
normalized input and metadata are intact. Failed unfolding does not invalidate
input preparation. Missing provenance, wrong energy/rate/branching fraction,
wrong Q² selection, or absent proton contract stops the add-on.

## What the error and optimization mean

The failed `--method pT2` script unfolds **positive `true_t`** using
**`reco_pT**2`**. Thus unsupported bin `[3]` is truth **1–5 GeV²**. The training
response has no accepted reconstructed support for that truth bin, while the
held-out truth yield is positive. Increasing projected luminosity cannot change
this response-support failure.

The audit identifies whether the observed finite-MC hole comes from absent
training truth, other reconstruction cuts, or the added proton cut. It reproduces
the legacy finite-value prefilter and random split. It also reports physical
nonnegative reconstructed-pT selection, finite-filter losses, and true-|t| tail
feed-in. Its no-proton comparison is diagnostic only. A zero in finite MC does
not by itself prove exactly zero physical acceptance.

The optimization has two distinct targets:

| Observable | Truth variable | Reconstructed variable | Fitted range |
|---|---|---|---|
| `pt2` | `true_pT**2` | `reco_pT**2` | 0–1 GeV² |
| `abs_t_via_pt2` | positive `true_t` | `reco_pT**2` | 0–1 GeV² |

**This is a restricted-range measurement. It does not recover or validate the
old 0–5 GeV² truth spectrum.** Reconstructed 1–5 GeV² is a control interval;
accepted outside-fit truth is retained as modeled feed-in and tested with
separate shape-only, tail-only and combined stresses. A truly unsupported tail
cannot be made measurable merely by merging it with an accepted bin: that would
make the result depend on the assumed within-bin spectrum. If full-range
measurement is essential, response MC and acceptance must be investigated first.

Binning is selected on training MC only. Existing expected-count, purity,
stability, effective-MC and conditioning gates are unchanged. Held-out
inference, coverage, response-MC variation and stress tests follow. At 1 fb⁻¹ a
state may have no qualifying binning; the program must report that outcome.
`COMPLETE` means execution completed, not that physics validation passed.
`PASS_CONDITIONAL` still excludes background, mass-state extraction, luminosity,
branching-ratio and flux uncertainties. Overall status is propagated to CSVs and
reports; unavailable results never fall back to older extraction.

## Outputs

The script prints both full output paths:

```text
Ypsilon-Fit-Framework-ProtonTagged/results/generated/runs/
  reopt-proton-9x130-L1p0-<UTC>-support/artifacts/t_support/
  reopt-proton-9x130-L1p0-<UTC>/artifacts/reoptimization/
```

The support directory contains `audit.json`, `response_support.csv`,
`truth_bin_cutflow.csv`, and `tail_feed_in.csv`.

The optimization directory contains `validation_status.csv`,
`differential_results.csv`, `summary.md`, `appendix_tables.tex`, per-case numerical
results and the available diagnostic figures. Inspect `validation_status.csv`
before using differential points or slope intervals in the note. Send the
support audit and optimization reports back for a production-data assessment.

This add-on does not mark the original failed task DONE, rerun independent full
workflow tasks, or replace notebook/note inputs. Those remain separate records.
After reviewing the new results, notebook inputs can be pointed at the accepted
new run explicitly.

## Manual commands and additional states

From the parent directory, define:

```bash
fw=/afs/cern.ch/work/g/gkrintir/private/Ypsilon-Fit-Framework-ProtonTagged
batch=allstates-proton-9x130-projected-L2p5-20261101T164736Z
```

Preview the configuration and sources (read-only):

```bash
python3 -B Ypsilon_9x130_Reoptimization_Update/run_reoptimization_9x130.py plan \
  --framework-root "$fw" --source-batch "$batch" \
  --states 1S --observables pt2 abs_t_via_pt2 --lumi 1
```

To run all three states once their **9×130 normalized inputs** are present:

```bash
python3 -B Ypsilon_9x130_Reoptimization_Update/run_reoptimization_9x130.py run \
  --framework-root "$fw" --source-batch "$batch" \
  --states 1S 2S 3S --observables pt2 abs_t_via_pt2 --lumi 1 \
  --run-id reopt-proton-9x130-L1-allstates-v1 --hash-inputs
```

The port also supports `--observables w y` as a separate fresh run. W uses
20–70 GeV and bin-integrated photon flux. Rapidity keeps the declared truth range
−1 to 3 and reports inclusive ep dσ/dy, with a separately labeled total-flux
companion. Feasibility and closure at this energy have to be established by that
run; configuration support is not physics validation.

To resume an interrupted optimization, repeat its exact Python `run` command
with its printed run ID and add `--resume`. Inputs, config, code and source-record
fingerprints must match. Do not use the shell helper to resume: it creates a new
run. A completed unavailable case is cached too; changes to settings require a
fresh run ID. The source batch itself is never resumed by this add-on.

## Explicit normalization

| State | Visible dielectron σ [fb] | B_ee | Q² maximum [GeV²] | All-generated yield at 1 fb⁻¹ |
|---|---:|---:|---:|---:|
| 1S | 119.685 | 0.0238 | 1 | 119.685 |
| 2S | 29.297 | 0.0191 | 0.01 | 29.297 |
| 3S | 23.447 | 0.0218 | 0.01 | 23.447 |

These are before reconstruction/fiducial selection. Independent numerical EPA
integration gives Φ=0.06880106724446 (1S) and 0.04959731276198 (2S/3S) over the
nominal W20–70 range; the integral stops at its physical support near W68.417.
This verifies numerical integration, not the generator executable provenance.
The state-dependent Q² and electron-ID settings are preserved; this is not a
matched-selection comparison across states.

## Software verification

See `VALIDATION.md`. Tests use synthetic inputs and the supplied framework code;
no claim is made that CERN ROOT production samples were run here.
