# From reconstructed events to analysis tuples

[Follow a reconstructed branch into a paper figure](https://gkrintir.github.io/ePIC-Upsilon/#w/reconstruction).

Use the existing merged samples for the analysis and plotting workflow. This guide is for **new tuple production** from downloaded ePIC reconstructed `events` trees. It does not run eSTARlight, detector simulation or EICrecon.

## What is recovered and what remains missing

`production/sources.json` records the exact origin and SHA256 of every recovered file. The latest uploaded producer and electron finder are preserved byte for byte under their canonical include names. Style files come from `gkrintir/lAger_JPsi_MuMu_Dataset` at commit `64a07c2a1dfd633ff3c9ef03f5e15982e7ea2f6d`. The six cards come from `eic/Upsilon_Sim_Cards/Jul26_Campaign_Input` at commit `f62965a1612fc734dbab18e141f4809740d6b234`.

The actual CERN Condor submit description, worker/environment setup, merge script/input manifest and per-energy source revisions still need to be exported from CERN. The new `production/jobs.py` is clearly separate from those missing historical scripts. Its local tests exercise planning and synthetic ROOT bookkeeping; they do **not** establish successful compilation with the user's ROOT/podio installation.

The recovered source has these important limitations:

| Item | Preserved behavior and implication |
| --- | --- |
| Beam constants | `E_e_beam=9`, `E_p_beam=275`; the planner blocks new 9×130 production |
| Proton candidate threshold | Uses `0.5 * E_p_beam`; a 9×130 card cannot change this C++ constant |
| Raw visible rates | 327.584, 84.780, 68.929 fb for 1S, 2S, 3S, inferred from the input filename; these are historical source values, **not the current assigned analysis rates** |
| Raw luminosity | 10 fb⁻¹ per job; each job assigns `sigma * 10 / N_job` |
| Merge normalization | `hadd` sums raw weights. It does not turn per-job weights into a valid complete-sample normalization |
| Generated rows | Source fills one tuple row per input entry and writes the total to `h_EventCount`; the wrapper checks both |
| Recoil proton tagging | Producer stores tag information. The downstream tagged analysis adds its explicit selection; production mode is not the analysis category |

Do not infer which exact executable produced existing EOS files from this archive alone. Do not edit beam constants or weights silently to make a new production pass. Supply the missing source/environment revision and audit it as a new revision first.

## Environment

The planner needs Python 3.9+. Executing or merging needs `uproot` in that Python environment, plus the user's tested ROOT/podio environment. `root` and `hadd` must be available in it. The C++ source includes `podio/ObjectID.h` and `podio/CollectionIDTable.h`; reconstructed electron ID additionally uses the ROOT-I/O dictionaries and, where available, `podio/utilities/RootHelpers.h`. It loads `libpodio`, `libpodioDict` and `libpodioRootIODict`.

The exact ROOT, podio, EDM4hep/EDM4eic and ePIC environment versions used at CERN were not supplied. There is therefore no guessed CVMFS `source` command here. Provide a tested setup script and a precise environment label, preferably a release/container digest. The setup script must be readable on workers, must set the ROOT/podio environment, and must not depend on interactive shell state. Its SHA256 is recorded and checked. Python with `uproot` must already work before this script: the worker performs bookkeeping in Python and applies the setup script to its ROOT/hadd subprocess.

Jobs use a shared filesystem: the input files, fresh work directory and environment script must be readable from the execute nodes. The submit file sets `should_transfer_files=NO`; it does not transfer inputs or arrange authentication. Check CERN access on the intended nodes before a large submission. Only a limited set of environment variables is inherited, including `ROOT_INCLUDE_PATH`; other setup belongs in the explicit setup script.

## Optional download

The recovered downloader requires `xrdfs` and `xrdcp`. Its dry run performs remote listings and updates local inventory/cache metadata, but copies no ROOT files:

```bash
python3 production/download/download_upsilon_fast.py \
  --dest /eos/cms/store/group/phys_top/gkrintir/upsilon_data \
  --dataset Upsilon1S/9x275/hiAcc --dry-run
```

Review the selected paths before removing `--dry-run`; `--max-files 1` is available for a transfer check. Its recovered defaults refer to `26.07.1/epic_craterlake`, `eSTARlight1.3.0-1.0`, and a `q2_0to0.01` remote path for every state. The recovered official 1S cards specify a Q² upper bound of 1 GeV², whereas 2S/3S specify 0.01 GeV². This path/card difference needs production provenance verification; a directory name does not establish the generated phase space. The downloader has not been changed to resolve that difference.

## 1. Prepare an explicit input list

Use one mounted reconstructed `.root` path per line, with no wildcards. Keep downloaded input files separate from tuple outputs. Every input path must identify the declared state using `1S`, `2S` or `3S`, because the preserved producer infers its raw rate from that path. Duplicated paths are rejected. One input file becomes one job.

For example, create a list of the immediate ROOT files in the downloaded input directory; this does not recurse into scenario output folders:

```bash
mkdir -p workspaces/manifests
python3 - <<'PY'
from pathlib import Path
source = Path('/eos/cms/store/group/phys_top/gkrintir/upsilon_data/Upsilon1S/9x275/hiAcc')
files = sorted(p.resolve() for p in source.glob('*.root') if p.is_file())
if not files:
    raise SystemExit('No reconstructed input ROOT files found; check the source directory')
with Path('workspaces/manifests/9x275-1S-inputs.txt').open('x') as out:
    out.write('\n'.join(map(str, files)) + '\n')
print(f'Listed {len(files)} inputs; inspect this list before planning.')
PY
```

Input file identity is recorded by path, size and modification time. `--hash-inputs` additionally reads and hashes every input; this can be expensive for a large campaign.

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

The producer arguments are:

| Planner option | C++ argument | Accepted choices |
| --- | --- | --- |
| `--data-version` | Third | `new` (default), `old` |
| `--femc` | Fourth | `full`, `remove_inner`, `remove_outer`, `keep_xpos`, `keep_yneg` |
| `--eta-split` | Fifth | Positive finite number; default 1.68 |
| `--electron-id` | Sixth | `legacy` or `reco`, explicitly required |

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

All jobs declared in the batch must have validated `DONE` results. The helper checks their frozen configuration, input identities, output SHA256 and generated-row totals. It uses an explicit `hadd @filelist`, so a large number of jobs does not overflow the shell argument limit. Existing outputs and provenance files are never replaced.

For the **1S, reconstructed electron ID, full fEMC** batch above, write directly to a fresh writable EOS pilot area with the directory layout expected by the analysis front end. Keep this base separate from the current production inputs:

```bash
export PILOT_INPUT_BASE=/eos/user/g/gkrintir/upsilon-production-pilot-v1
mkdir -p "$PILOT_INPUT_BASE/Upsilon1S/9x275/hiAcc/output_newID"

python3 "$PRODUCTION_WORK/snapshot/jobs.py" merge \
  --batch "$PRODUCTION_WORK" \
  --output "$PILOT_INPUT_BASE/Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root"

python3 "$PRODUCTION_WORK/snapshot/jobs.py" merge \
  --batch "$PRODUCTION_WORK" \
  --output "$PILOT_INPUT_BASE/Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root" \
  --execute
```

The first command verifies and prints the plan. The second invokes `hadd` without `-f`, checks the merged tuple and `h_EventCount` against the summed input counts, and writes `.merge.json` and `.inputs.txt` alongside the output. Partial outputs are retained on failure. No weights are renormalized during merge.

The analysis front end uses these exact paths beneath `input_base`; capitalization matters:

| State | Expected path beneath the pilot base | Matching production settings |
| --- | --- | --- |
| 1S | `Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root` | `--state 1S --electron-id reco --femc full` |
| 2S | `Upsilon2S/9x275/hiAcc/output/merged_Upsilon_final.root` | `--state 2S --electron-id legacy --femc full` |
| 3S | `Upsilon3S/9x275/hiAcc/output/merged_Upsilon_final.root` | `--state 3S --electron-id legacy --femc full` |

For 2S/3S, create separate production batches using their own input lists and the matching options above, create the corresponding EOS output directory, and run the same merge command with that batch and exact destination. Never place a reduced-fEMC or different-ID sample under one of these baseline paths. Different scenarios require a separate, explicitly configured analysis campaign.

## 5. Normalize and audit in the analysis framework

Set `input_base` in `configs/local.json` to the pilot base used above. Preserve any intentional custom `workspace_root`; with the default local workspaces the file is:

```json
{
  "workspace_root": "workspaces",
  "input_base": "/eos/user/g/gkrintir/upsilon-production-pilot-v1"
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

The CLI retains the current state/energy-specific visible rates, generated phase-space settings, tuple luminosity and branching fractions. `--lumi 2.5` sets the projected luminosity; it does not change the producer's raw weights or the source tuple luminosity. The framework's preparation step normalizes the complete generated sample; its audit must verify the expected total weight, row preservation, saved cross sections and flux/selection consistency. Do not interpret the merged producer `weight` sum as an expected experimental yield, and do not apply the producer's historical cross sections a second time. These numerical checks still do not authenticate the producer revision or resolve the input generation provenance.

Keep the batch `plan.json`, individual `job.json`/`result.json`, source snapshot, logs and merge audit with the analysis provenance. This keeps the existing samples and previous tagged/untagged results independently reproducible.

## Checks shipped with this repository

```bash
python3 -m unittest discover -s tests -p test_production_plan.py -v
```

Tests cover source hashes, read-only planning, duplicate/state/energy restrictions, fresh-output protection, frozen selections, explicit execution and synthetic ROOT entry-count/merge bookkeeping. Synthetic ROOT checks need `uproot` and NumPy. No Condor job or real producer invocation is made by these tests. The ROOT/podio compilation and a CERN worker smoke test remain required in the user's actual environment.

Interface references: [HTCondor submit description](https://htcondor.readthedocs.io/en/24.x/man-pages/condor_submit.html), [shared filesystem/file transfer behavior](https://htcondor.readthedocs.io/en/main/users-manual/file-transfer.html), [ROOT hadd indirect input lists](https://root.cern.ch/doc/v632/hadd_8cxx.html).
