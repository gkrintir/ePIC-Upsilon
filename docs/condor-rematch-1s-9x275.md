# Condor rematching for Υ(1S) at 9×275 GeV

Run from **bare lxplus**, in an AFS checkout on `main`. Each worker reads an original EDM4eic file from EOS and writes its corrected ROOT tuple directly to EOS. The merge reads those tuples from EOS. The plan, logs and completion receipts live under `workspaces/` on AFS.

If the one-file Condor test (`condor-1S-9x275-check-v2`) already succeeded, keep its receipt. To process the full sample from the original files while an earlier campaign is active, [start an independent EOS run](#start-an-independent-eos-run). To reuse finished AFS tuples after the older campaign leaves the queue, [switch that campaign to EOS](#switch-an-existing-afs-campaign-to-eos).

## Set the paths

Run this block again from the repository root after reconnecting. The example configuration points to the July campaign; set `UPSILON_DATA` yourself if your mounted copy is elsewhere.

```bash
git pull --ff-only origin main
export UPSILON_DATA="${UPSILON_DATA:-$(python3 -c 'import json; print(json.load(open("configs/local.example.json"))["input_base"])')}"
export CORRECTED_BASE="${CORRECTED_BASE:-$UPSILON_DATA/rematch3011-v1}"
src="$UPSILON_DATA/Upsilon1S/9x275/hiAcc"
name=eSTARlight1.3.0-1.0_Upsilon1S_9x275_q2_0to0.01_hiAcc_ab_run000.0000.eicrecon.edm4eic.root
image=/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable
trial="$PWD/workspaces/condor-1S-9x275-check-v2"
full="$PWD/workspaces/condor-1S-9x275-full-grouped-eos-v1"
klist -s || kinit "$USER@CERN.CH"
head -c 4 "$src/$name" && echo
```

The last command should print `root`. Workers also need write access to `CORRECTED_BASE`. For interactive Apptainer access, see the [Kerberos cache and EOS instructions](rematch-3011.md#start-on-lxplus-and-run-one-file).

## Submit the one-file test

Run this block once. It creates its own input list. If an earlier attempt reported a missing `v1/input.txt`, use this complete block. The corrected launcher in [commit 70d8344](https://github.com/gkrintir/ePIC-Upsilon/commit/70d8344a8c3a6c6647918958985169ad6b6e6dc6) locates the frozen Python worker on AFS when Condor starts the launcher in its execution directory.

```bash
(
  set -euo pipefail
  mkdir -p "$trial"
  if [[ -e "$trial/batch" ]]; then
    echo "Batch already exists at $trial/batch. Inspect its result and queue entry first." >&2
    exit 1
  fi
  printf '%s\n' "$src/$name" > "$trial/input.txt"
  python3 production/jobs.py plan \
    --input-list "$trial/input.txt" --work-dir "$trial/batch" \
    --tuple-output-dir "$CORRECTED_BASE/Upsilon1S/9x275/hiAcc/output_newID/pilot-check-v2" \
    --state 1S --beam 9x275 --electron-id reco-rematch3011 --femc full \
    --environment-label eic_xl:26.07.1-stable --job-flavour longlunch \
    --hash-inputs --write
  grep -Fx '+AccountingGroup = "group_u_CMS.CAF.COMM"' "$trial/batch/submit.sub"
  condor_submit -maxjobs 1 -append 'getenv = False' \
    -append "+SingularityImage = \"$image\"" "$trial/batch/submit.sub" \
    | tee "$trial/submission.txt"
)
```

Condor prints a cluster ID and returns control to the shell. You can log out. A `worker --execute` command at the prompt would run the job locally.

## Check the test output

From the repository root:

```bash
batch="$PWD/workspaces/condor-1S-9x275-check-v2/batch"
condor_q "$USER"
cat "$batch/jobs/job_000000/result.json"
```

Look for `"status": "DONE"`, `"output_rows": 100` and the `"output"` path. A test planned with the command above records an EOS path. Earlier tests wrote their tuple to AFS and verified EOS reads. The recorded input hash is `f2701e1d496d931b117afd5760efe15c909061668890e2d67d6fd28d29856908`. Of the 100 rows, 23 pass the analysis cuts before the proton tag and 13 after it.

If the result file is missing or the worker failed, check the job history and logs. Some files appear only after the job starts.

```bash
cat "$batch/logs/000000.log"                # Condor events: execution, hold or exit
cat "$batch/logs/000000.err"                # Worker errors
cat "$batch/logs/000000.out"                # Worker stdout
cat "$batch/jobs/job_000000/root.log"       # ROOT output once the producer starts
```

An error about `/srv/jobs.py` identifies the old frozen launcher. Confirm the old job has left the queue, pull `main`, and create a fresh plan with the current launcher. Keep the failed directory for its logs. An existing batch continues to use its original snapshot.

## Plan and submit the full sample

The full-sample launcher uses **up to 100 files per Condor job** by default. For example, 11,376 input files require 114 Condor jobs, with 76 files in the last job. The launcher prints the actual count from your EOS inventory. In the archived 2S campaign, most Condor jobs processed 140 files; the 1S pilot and the full sample can have different runtimes.

On bare lxplus, after the one-file Condor pilot is `DONE` with 100 rows, run:

```bash
git pull --ff-only origin main && bash production/submit_1s_9x275.sh --submit
```

To see the file and job counts first, omit `--submit`. Run the same command with `--submit` afterwards to submit the frozen plan. You can set `FILES_PER_JOB=50`, for example, before the first plan; give the same value when submitting. Use `tmux` if you want to disconnect while the AFS plan is being written. After Condor returns the cluster ID, it owns the jobs.

The launcher checks the pilot receipt and writes the full input list and plan log under `workspaces/condor-1S-9x275-full-grouped-eos-v1/`. Each original input appears once in the full plan. In the default configuration:

| File | Location |
| --- | --- |
| First per-input tuple (EOS) | `$CORRECTED_BASE/Upsilon1S/9x275/hiAcc/output_newID/intermediate/job_000000.root` |
| Its receipt (AFS) | `$full/batch/jobs/job_000000/result.json` |
| Condor logs (AFS) | `$full/batch/logs/` |
| Final merged tuple (EOS) | `$CORRECTED_BASE/Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root` |

`CORRECTED_BASE` defaults to `/eos/cms/store/group/phys_top/gkrintir/upsilon_data/rematch3011-v1` with the repository configuration. To change just the per-input directory, export `TUPLE_OUTPUT_DIR` to a fresh absolute EOS path before the first plan. The planner creates this directory. If you override either variable, use the same exported values when you return to the campaign. Inspect the frozen destination or the first completed receipt on lxplus:

```bash
python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["tuple_output_dir"])' "$full/batch/plan.json"
cat "$full/batch/jobs/job_000000/result.json"
```

The planner records `files_per_job`, `condor_job_count` and the number of per-file work items (`job_count`) in `batch/plan.json`. The submission asks for `group_u_CMS.CAF.COMM`, `longlunch` and `eic_xl:26.07.1-stable`. Each Condor process runs the files in its group sequentially, using the same per-file ROOT invocation and weights as the successful pilot. Every input has its own output and receipt. The final merge follows the original manifest order.

Already submitted jobs use the output paths in their frozen plans. Choose a named run below for independent reprocessing, or the AFS-to-EOS switch to reuse completed work.

### Start an independent EOS run

After the one-file pilot succeeds, run this from the repository root on bare lxplus. `fresh-01` identifies a new full pass over the original EDM4eic files. It can start while the earlier campaign is running.

```bash
git pull --ff-only origin main && bash production/submit_1s_9x275.sh --new-run fresh-01 --submit
```

The AFS plan is `workspaces/condor-1S-9x275-full-grouped-eos-fresh-01/batch/plan.json`. Each worker writes its ROOT output directly to `$CORRECTED_BASE/Upsilon1S/9x275/hiAcc/output_newID/runs/fresh-01/intermediate/` (for example, `job_000000.root`). If you previously exported `TUPLE_OUTPUT_DIR`, unset it before planning to use this default, or point it to a fresh absolute EOS directory dedicated to this run.

Select this run explicitly for status and retries:

```bash
batch="$PWD/workspaces/condor-1S-9x275-full-grouped-eos-fresh-01/batch"
python3 production/rematch_status.py --batch "$batch" --queue-only
python3 production/rematch_status.py --batch "$batch"
```

Once this run's queue drains, preview any retry with `python3 production/retry_rematch_1s_9x275.py --batch "$batch"`. Add `--submit` to submit it. When every input is DONE, merge **inside the pinned ePIC shell** into this run's own directory:

```bash
(
  set -euo pipefail
  batch="$PWD/workspaces/condor-1S-9x275-full-grouped-eos-fresh-01/batch"
  data="${UPSILON_DATA:-$(python3 -c 'import json; print(json.load(open("configs/local.example.json"))["input_base"])')}"
  run_dir="${CORRECTED_BASE:-$data/rematch3011-v1}/Upsilon1S/9x275/hiAcc/output_newID/runs/fresh-01"
  mkdir -p "$run_dir"
  python3 "$batch/snapshot/jobs.py" merge --batch "$batch" --output "$run_dir/merged_Upsilon_final.root" --execute
)
```

For the full July inventory, check for `"status": "DONE"` and `"output_rows": 1137600` in `runs/fresh-01/merged_Upsilon_final.root.merge.json`.

### Switch an existing AFS campaign to EOS

Both earlier layouts wrote per-input tuples to AFS: the one-file campaign in `condor-1S-9x275-full-v1` and the grouped campaign in `condor-1S-9x275-full-grouped-v1`. Keep these directories, submission records and logs in place. On the original lxplus submit scheduler, check `condor_q "$USER"`. Once the old campaign and its retries have finished or been removed, run:

```bash
git pull --ff-only origin main && bash production/submit_1s_9x275.sh --after-removal --submit
```

The launcher checks the scheduler and refuses to switch while old jobs are active, including those in the removing state. Recovery checks the inventories and reconstruction settings, copies verified completed tuples to EOS, and queues groups with outstanding files. Verified files in a partially completed group are skipped. To prepare the EOS plan and copy the completed tuples now, omit `--submit`. Add it when ready to submit the remaining groups.

For a custom campaign, use `production/jobs.py plan ... --files-per-job 100`. Omitting this option keeps one input per Condor process. Grouping is a scheduler setting; the underlying per-file manifest, ROOT selection and merge interface are the same. Unhashed planning records path, size and modification time; `--hash-inputs` additionally reads every original file to record and verify its SHA-256.

## Monitor the grouped campaign

On bare lxplus, the queue-only command prints once and returns to the prompt. These commands use the default campaign. For a named run, add its explicit `--batch` argument as [shown above](#start-an-independent-eos-run).

```bash
git pull --ff-only origin main
python3 production/rematch_status.py --queue-only
python3 production/rematch_status.py
```

This reads the grouped campaign's recorded cluster IDs and worker receipts. It shows idle, running and held **Condor processes** separately from **input files** marked DONE, FAILED, untouched or started without a final receipt. One process can handle up to 100 files. The report includes hold reasons and example file names, errors and log paths. The row total comes from the DONE receipts. The final merge verifies the ROOT contents.

The full report scans per-input receipts on AFS and checks the tuple directory. To inspect a particular frozen plan, pass `--batch "$full/batch"` or the earlier AFS batch's path. For a continuously refreshed scheduler display, run:

```bash
watch -n 60 'python3 production/rematch_status.py --queue-only'
```

Ctrl+C stops the display. If a process is held, inspect its stated reason or use `condor_q -hold CLUSTER_ID` with the cluster ID from the report. For a file-only report without a scheduler query, run `python3 production/rematch_status.py --no-queue`.

## Retry failed or unfinished files

For a named run, pass its `--batch` path on both commands below.

Use the same bare lxplus submit scheduler as the initial grouped run. Inspect the reported errors first; repeated permission or input-format errors need their cause fixed. A running group can continue processing other inputs after an individual file fails. A started/partial file can still be running. Once the campaign's queue has drained (including held, transferring and removing processes), preview what would be retried:

```bash
python3 production/retry_rematch_1s_9x275.py
```

To submit that retry:

```bash
python3 production/retry_rematch_1s_9x275.py --submit
```

The retry checks completed tuples and queues only groups containing unfinished inputs. Each retry records its submission and logs in `batch/retries/attempt-…/` on AFS. It archives any partial EOS tuples under `.retries/attempt-…/` inside the tuple directory, and records that directory in `tuple-archive.txt` on AFS. Run the same status command to follow the retry. Leave submission records, worker locks and results in place. Resolve any missing cluster ID on the scheduler and finish manually launched workers before retrying.

## Check and merge the full sample

The commands below merge the default campaign. For a named run, use its [run-specific batch and merge path](#start-an-independent-eos-run).

For the original July 1S 9×275 inventory, require **all input files DONE and 1,137,600 generated rows** in the status report. Investigate a different total. An empty Condor queue alone can also mean that a launcher exited before it opened its first file.

Merge **inside the pinned ePIC shell** with access to the AFS batch and writable EOS output. Use the [lxplus container launch](rematch-3011.md#start-on-lxplus-and-run-one-file) with the same EOS credentials as the interactive test. The default output base is a separate `rematch3011-v1` directory beneath the July data area. Set `CORRECTED_BASE` to another writable absolute directory if needed. In that shell, from the repo root, run:

```bash
(
  set -euo pipefail
  batch="$PWD/workspaces/condor-1S-9x275-full-grouped-eos-v1/batch"
  data="${UPSILON_DATA:-$(python3 -c 'import json; print(json.load(open("configs/local.example.json"))["input_base"])')}"
  CORRECTED_BASE="${CORRECTED_BASE:-$data/rematch3011-v1}"
  out="$CORRECTED_BASE/Upsilon1S/9x275/hiAcc/output_newID/merged_Upsilon_final.root"
  mkdir -p "$(dirname "$out")"
  python3 "$batch/snapshot/jobs.py" merge --batch "$batch" --output "$out" --execute
  python3 - "$out.merge.json" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    result = json.load(stream)
print("Merge:", result["status"], "  Output rows:", result["output_rows"])
if result["status"] != "DONE" or result["output_rows"] != 1137600:
    raise SystemExit("Check the merged sample against the original 1S 9x275 inventory")
PY
)
```

The merger verifies each declared worker output and its checksum. It writes the input list and merge receipt alongside the ROOT file. The full merged tuple still has per-job weights. Point a new [analysis workspace](analysis.md) at `CORRECTED_BASE` to normalize the complete generated sample and rerun the fits. For the tagged 1S paper projection choose `--selection proton-tagged`, `--energy 9x275`, `--states 1S` and `--lumi 2.5`. See the [analysis commands](rematch-3011.md#recalculate-the-physics-results) when all three states are ready.
