# From reconstructed events to analysis tuples

[Follow a reconstructed branch into a paper figure](https://gkrintir.github.io/ePIC-Upsilon/#w/reconstruction).

Use the existing merged samples for the analysis and plotting workflow. This guide is for **new tuple production** from downloaded ePIC reconstructed `events` trees. It does not run eSTARlight, detector simulation or EICrecon.

For the July samples affected by EICrecon #3011, follow the [ECAL rematching pilot and production instructions](rematch-3011.md).

## Production sources

Choose `--beam 9x130` or `--beam 9x275`. The planner selects the corresponding preserved C++ macro, checks its beam constants and freezes it with the electron finder and style files in the batch. The beam setting also determines the proton-candidate momentum threshold, half the proton-beam energy. The tuple stores the candidate information for downstream selections. [Source paths and SHA-256 values](../production/sources.json)

For the 2S 9×130 sample, the recovered input list and per-file worker messages account for 1,058,300 events. The track-matching counters agree with the distributed ROOT file in all 10,583 consecutive blocks of 100 events. [Original campaign and merge check](../production/cern-2026-10/README.md)

Each recovered macro sets its raw per-job event weight to the state’s visible cross section times 10 fb⁻¹ divided by the generated events in that job. Those constants reproduce the observed raw weights for some samples. The 1S 9×130 and the 2S and 3S 9×275 weights require different constants from the corresponding preserved macro. For physics yields the analysis normalizes all generated rows after merging, using the assigned rate for each state and energy. [Sample-by-sample evidence and missing records](missing-production-sources.md)

## Environment

The two pilot inputs were written with ROOT 6.38/04 and podio 1.7.0. They store collection IDs in `events___CollectionTypeInfo`. The reconstructed electron finder needs the matching podio headers and ROOT dictionaries, including the dictionary for `podio::root_utils::CollectionWriteInfo`, to resolve each track–cluster association. An older podio v01-02 installation compiled the macros but stopped when it read these files. Keep reconstructed electron ID enabled.

Use Python 3.9+ with `uproot`. The worker needs `root` on its path, and the merge needs `hadd`. Set up ROOT with C++20 and podio with its ROOT I/O libraries. The [pilot comparison](missing-production-sources.md#two-reconstructed-inputs) records the tested versions.

The archived CERN submit configuration selected `eic_xl:nightly`. Its immutable image digest is unknown. For a new batch, set `PRODUCTION_ENV_LABEL` to your tested release or image digest. Set `PRODUCTION_SETUP` to the absolute path of a noninteractive setup script that works on the execution nodes. The planner records its SHA-256. Python with `uproot` must be available in the launching environment, and the setup script configures the ROOT and `hadd` subprocesses.

The jobs use mounted shared filesystems. Workers need read access to the inputs and setup script, and write access to the work and tuple directories. The generated submit file sets `should_transfer_files=NO`. An EOS-output pilot tests the worker's write permissions at the chosen destination.

## One-file pilot

Start from one reconstructed `*.eicrecon.edm4eic.root` file. The example uses Υ(2S) at 9×130. For the paper-energy check, set `state=1S` and `beam=9x275` and use the corresponding input. Run from the repository root.

```bash
state=2S
beam=9x130
input=${PILOT_INPUT:?Set the absolute path to your reconstructed ROOT file}
batch="$PWD/workspaces/production/$state-$beam-pilot-v1"
mkdir -p workspaces/manifests
list="workspaces/manifests/$state-$beam-pilot.txt"
printf '%s\n' "$input" > "$list"

python3 production/jobs.py plan \
  --input-list "$list" --work-dir "$batch" --beam "$beam" --state "$state" \
  --electron-id reco --femc full --hash-inputs \
  --environment-label "${PRODUCTION_ENV_LABEL:?Set the tested release or digest}" \
  --setup-script "${PRODUCTION_SETUP:?Set the absolute setup script path}" --write

python3 "$batch/snapshot/jobs.py" worker "$batch/jobs/job_000000/job.json" --execute
python3 production/jobs.py merge --batch "$batch" --output "$batch/merged.root" --execute
```

Check `jobs/job_000000/result.json` and `merged.root.merge.json` in the batch directory for the final status and row count. To identify the historical 2S input from the local production ZIP, see the [CERN campaign record](../production/cern-2026-10/README.md).

If the input belongs to a distributed reference tuple, compare the event rows after merging. This locates the pilot using six generated quantities, then checks every saved branch. Add `--json` for the per-branch differences.

```bash
python3 production/compare_pilot.py "$batch/merged.root" "${REFERENCE_TUPLE:?Set the matching merged_Upsilon_final.root path}"
```

## Optional download

The recovered downloader requires `xrdfs` and `xrdcp`. Its dry run performs remote listings and updates local inventory/cache metadata, but copies no ROOT files:

```bash
python3 production/download/download_upsilon_fast.py \
  --dest "${RECO_INPUT_BASE:?Set the destination for reconstructed files}" \
  --dataset Upsilon1S/9x275/hiAcc --dry-run
```

Review the selected paths before removing `--dry-run`; `--max-files 1` is available for a transfer check. Its recovered defaults refer to `26.07.1/epic_craterlake`, `eSTARlight1.3.0-1.0`, and a `q2_0to0.01` remote path for every state. The recovered official 1S cards specify a Q² upper bound of 1 GeV², whereas 2S/3S specify 0.01 GeV². This path/card difference needs production provenance verification; a directory name does not establish the generated phase space. The downloader has not been changed to resolve that difference.

## 1. Prepare an input list for a full batch

Use one mounted reconstructed `.root` path per line, with no wildcards. Keep downloaded input files separate from tuple outputs. Every input path must identify the declared state using `1S`, `2S` or `3S`, because the preserved producer infers its raw rate from that path. Duplicated paths are rejected. Each input gets a separate ROOT output and receipt.

Export `RECO_INPUT_BASE` as the absolute path to your downloaded datasets. This example lists the reconstructed files directly under the 1S input directory.

```bash
mkdir -p workspaces/manifests
python3 - <<'PY'
from pathlib import Path
import os
source = Path(os.environ['RECO_INPUT_BASE']) / 'Upsilon1S/9x275/hiAcc'
files = sorted(p.resolve() for p in source.glob('*.eicrecon.edm4eic.root') if p.is_file())
if not files:
    raise SystemExit('No reconstructed input ROOT files found; check the source directory')
with Path('workspaces/manifests/9x275-1S-inputs.txt').open('x') as out:
    out.write('\n'.join(map(str, files)) + '\n')
print(f'Listed {len(files)} inputs; inspect this list before planning.')
PY
```

Input file identity is recorded by path, size and modification time. `--hash-inputs` additionally reads and hashes every input; this can be expensive for a large campaign.

For the full 1S 9×275 rematching campaign, [group the inputs at submission](condor-rematch-1s-9x275.md#plan-and-submit-the-full-sample). The planner option `--files-per-job N` assigns up to N consecutive input files to one Condor process. It defaults to 1 for existing callers. Outputs and completion records stay per file; the frozen `merge` command covers the whole input list.

## 2. Plan, then create a fresh batch

Run these commands from the repository root. Set the environment values to the actual tested setup; the example intentionally requires these values rather than selecting an unverified release:

```bash
export PRODUCTION_SETUP=/absolute/path/to/your/tested_epic_setup.sh
export PRODUCTION_ENV_LABEL='your exact tested release or container digest'
export PRODUCTION_WORK="$PWD/workspaces/production/9x275-1S-reco-full-v1"

python3 production/jobs.py plan \
  --input-list workspaces/manifests/9x275-1S-inputs.txt \
  --work-dir "$PRODUCTION_WORK" \
  --beam 9x275 --state 1S --electron-id reco --femc full \
  --environment-label "$PRODUCTION_ENV_LABEL" \
  --setup-script "$PRODUCTION_SETUP" --hash-inputs
```

The command above prints the plan without creating a directory. Repeat it with `--write` to create the batch. It refuses an existing directory; old jobs and results stay intact. Resource defaults are one CPU, 4096 MB memory and 8192 MB disk; adjust `--memory-mb` and `--disk-mb` after a representative job. `--job-flavour` is optional and must be supplied explicitly for the chosen CERN queue.

To write every intermediate tuple directly to EOS, add `--tuple-output-dir "$TUPLE_OUTPUT_DIR"` when planning. Set `TUPLE_OUTPUT_DIR` to a fresh absolute EOS path, separate from `PRODUCTION_WORK`. On `--write`, the planner creates that directory. Workers write one tuple per input (`job_000000.root`, `job_000001.root`, …) and record each path in its `result.json` under `$PRODUCTION_WORK/jobs/`. The merge reads the EOS tuples in manifest order. Omitting this option selects the original `jobs/job_000000/output.root` layout. The [1S 9×275 full-sample launcher](condor-rematch-1s-9x275.md#plan-and-submit-the-full-sample) selects EOS by default.

The producer arguments are:

| Planner option | C++ argument | Accepted choices |
| --- | --- | --- |
| `--data-version` | Third | `new` (default), `old` |
| `--femc` | Fourth | `full`, `remove_inner`, `remove_outer`, `keep_xpos`, `keep_yneg` |
| `--eta-split` | Fifth | Positive finite number; default 1.68 |
| `--electron-id` | Sixth, or rematching wrapper | `legacy`, `reco` or [`reco-rematch3011`](condor-rematch-1s-9x275.md), explicitly required |

`reco` selects the new reconstructed matching mode. The preserved half-plane scenarios support only `legacy`; unsupported combinations are rejected. These options do not implement a new proton-tag analysis category.

## 3. Run one job locally, then submit the remaining jobs

Use the frozen worker in the batch so later repository edits cannot alter an existing plan:

```bash
python3 "$PRODUCTION_WORK/snapshot/jobs.py" worker \
  "$PRODUCTION_WORK/jobs/job_000000/job.json"

python3 "$PRODUCTION_WORK/snapshot/jobs.py" worker \
  "$PRODUCTION_WORK/jobs/job_000000/job.json" --execute
```

The first command prints the exact ROOT invocation. The second runs it. Inspect `jobs/job_000000/root.log` and `result.json`. A successful job must return zero, contain `UpsilonTree`, preserve the input entry count, and have `h_EventCount` equal to all tuple rows. These checks validate bookkeeping, not all physical observables or input schema behavior. Review a representative tuple/cutflow before a large batch.

After a successful local job, create a submit file containing only untouched jobs:

```bash
python3 "$PRODUCTION_WORK/snapshot/jobs.py" pending \
  --batch "$PRODUCTION_WORK" --write

# Inspect pending.sub, check your queue, then submit once:
condor_submit "$PRODUCTION_WORK/pending.sub"
condor_q
```

The planner never invokes `condor_submit`. `pending` checks recorded files, not the live queue; do not submit the same batch twice. With a one-file batch there are no remaining jobs and no `pending.sub` is created. Failed/partial jobs are not silently skipped: inspect the error and create a fresh attempt with an explicit revised input list. A process killed before final bookkeeping can leave an execution lock; preserve and inspect that attempt rather than deleting its evidence to force a rerun.

## 4. Merge exactly the completed batch

For any state or beam setting, the optional [readable-input merge](condor-rematch-1s-9x275.md#merge-the-readable-inputs) can exclude confirmed corrupt original files and record them in a separate merge receipt.

A strict merge requires validated `DONE` results for every planned job. The helper checks their frozen configuration, input identities, output SHA256 and generated-row totals. It uses an explicit `hadd @filelist`, so a large number of jobs does not overflow the shell argument limit. Existing outputs and provenance files are never replaced.

For the **1S, reconstructed electron ID, full fEMC** batch above, write directly to a fresh writable EOS pilot area with the directory layout expected by the analysis front end. Keep this base separate from the current production inputs:

```bash
export PILOT_INPUT_BASE="${PILOT_INPUT_BASE:?Set a fresh absolute output directory on shared storage}"
mkdir -p "$PILOT_INPUT_BASE/Upsilon1S/9x275/hiAcc/output_newID"

python3 production/jobs.py merge \
  --batch "$PRODUCTION_WORK" \
  --output "$PILOT_INPUT_BASE/Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root" \
  --execute
```

Use the current checkout’s `production/jobs.py` to merge an existing frozen plan. A single invocation with `--execute` verifies the inputs and then runs `hadd`. Omit `--execute` for a validation-only preview, which repeats the checks if you later execute. The merger checks the final tree and `h_EventCount` against the summed input counts and writes `.merge.json` and `.inputs.txt` alongside the output. Partial outputs are retained on failure. No weights are renormalized during merge.

For large batches, use `--check-workers N` for concurrent validation and `--hadd-jobs N` for ROOT multiprocessing. Both limits default to at most four processes, bounded by visible CPUs. Choose values within your CPU allocation. Small merges use fewer processes. The input list stays in manifest order. Temporary `hadd` partials go to the final output directory unless you choose another existing directory with `--hadd-temp-dir`. Use `--check-workers 1 --hadd-jobs 1` for serial operation. The receipt records timings for both stages and the final check.

The analysis front end uses these paths beneath `input_base`. All six [merged input tuples](../inputs/paper_inputs.csv) record `ele_id_mode=1` and `femc_mode=0`, corresponding to reconstructed track–cluster matching and the full forward calorimeter. Use the following settings to match those modes in a new 9×275 pilot.

| State | Expected path beneath the pilot base | Matching production settings |
| --- | --- | --- |
| 1S | `Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root` | `--state 1S --electron-id reco --femc full` |
| 2S | `Upsilon2S/9x275/hiAcc/output/merged_Upsilon_final.root` | `--state 2S --electron-id reco --femc full` |
| 3S | `Upsilon3S/9x275/hiAcc/output/merged_Upsilon_final.root` | `--state 3S --electron-id reco --femc full` |

For 2S/3S, create separate production batches using their own input lists and the matching options above, create the corresponding EOS output directory, and run the same merge command with that batch and exact destination. Never place a reduced-fEMC or different-ID sample under one of these baseline paths. Different scenarios require a separate, explicitly configured analysis campaign.

## 5. Normalize and audit in the analysis framework

Set `input_base` in `configs/local.json` to the absolute path printed by `echo "$PILOT_INPUT_BASE"`. For example, with the default workspace directory:

```json
{
  "workspace_root": "workspaces",
  "input_base": "/path/to/upsilon-production-pilot-v1"
}
```

Create a new analysis workspace, then plan, validate and run **only 1S** for the single pilot produced above:

```bash
python3 upsilon.py init --name pilot-275-tagged-v1 --selection proton-tagged
python3 upsilon.py doctor --workspace pilot-275-tagged-v1

python3 upsilon.py analysis plan \
  --workspace pilot-275-tagged-v1 --energy 9x275 --states 1S \
  --lumi 2.5 --batch-id pilot-275-1S-L2p5-v1

python3 upsilon.py analysis validate \
  --workspace pilot-275-tagged-v1 --energy 9x275 --states 1S \
  --lumi 2.5 --batch-id pilot-275-1S-L2p5-v1

python3 upsilon.py analysis run \
  --workspace pilot-275-tagged-v1 --energy 9x275 --states 1S \
  --lumi 2.5 --batch-id pilot-275-1S-L2p5-v1 --hash-inputs

python3 upsilon.py status --workspace pilot-275-tagged-v1
```

Inspect the `plan` output and confirm the input is the pilot EOS path before running. Once all three pilot files are present and reviewed, use `--states 1S 2S 3S` and a fresh batch ID such as `pilot-275-allstates-L2p5-v1`. For an analysis without an explicit proton tag, initialize another fresh workspace with `--selection inclusive`. Neither choice reruns tuple production.

The CLI retains the current state/energy-specific visible rates, generated phase-space settings, tuple luminosity and branching fractions. `--lumi 2.5` sets the projected luminosity; it does not change the producer's raw weights or the source tuple luminosity. The framework's preparation step normalizes the complete generated sample; its audit must verify the expected total weight, row preservation, saved cross sections and flux/selection consistency. Use the normalized tuple for projected yields. The original weighted rows and job records remain available for production checks.

Keep the batch `plan.json`, individual `job.json`/`result.json`, source snapshot, logs and merge audit with the analysis provenance. This keeps the existing samples and previous tagged/untagged results independently reproducible.

## Checks shipped with this repository

```bash
python3 -m unittest discover -s tests -p test_production_plan.py -v
```

Tests cover source hashes, read-only planning, duplicate/state/energy restrictions, fresh-output protection, frozen selections, explicit execution and synthetic ROOT entry-count/merge bookkeeping. Synthetic ROOT checks need `uproot` and NumPy. No Condor job or real producer invocation is made by these tests. Run the one-file pilot in your production environment to exercise the event readers.

Interface references: [HTCondor submit description](https://htcondor.readthedocs.io/en/24.x/man-pages/condor_submit.html), [shared filesystem/file transfer behavior](https://htcondor.readthedocs.io/en/main/users-manual/file-transfer.html), [ROOT hadd indirect input lists](https://root.cern.ch/doc/v632/hadd_8cxx.html).
