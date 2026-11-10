# Condor rematching for Υ(1S) at 9×275 GeV

Run the test and submit the full July sample from **bare lxplus**, in an AFS checkout of this repository on `main`. Condor reads the EDM4eic files directly from EOS. Plans, worker tuples and logs go under `workspaces/` on AFS. Check your AFS quota for one small tuple and its records per input. The full production groups inputs into Condor jobs. Keep a separate writable EOS directory for the final corrected tuple.

If you already submitted **condor-1S-9x275-check-v2**, go straight to [Check the test output](#check-the-test-output). An earlier successful interactive test verifies the physics on the login node. This Condor test also checks that a batch worker can read EOS and write to AFS.

## Set the paths

Run this block again from the repository root after reconnecting. The example configuration points to the July campaign; set `UPSILON_DATA` yourself if your mounted copy is elsewhere.

```bash
git pull --ff-only origin main
export UPSILON_DATA="${UPSILON_DATA:-$(python3 -c 'import json; print(json.load(open("configs/local.example.json"))["input_base"])')}"
src="$UPSILON_DATA/Upsilon1S/9x275/hiAcc"
name=eSTARlight1.3.0-1.0_Upsilon1S_9x275_q2_0to0.01_hiAcc_ab_run000.0000.eicrecon.edm4eic.root
image=/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable
trial="$PWD/workspaces/condor-1S-9x275-check-v2"
full="$PWD/workspaces/condor-1S-9x275-full-grouped-v1"
klist -s || kinit "$USER@CERN.CH"
head -c 4 "$src/$name" && echo
```

The last command should print `root`. Resolve an EOS read error on bare lxplus before submitting. For interactive Apptainer access, see the [Kerberos cache and EOS instructions](rematch-3011.md#start-on-lxplus-and-run-one-file).

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

Look for `"status": "DONE"` and `"output_rows": 100`. The recorded input hash in `plan.json` is `f2701e1d496d931b117afd5760efe15c909061668890e2d67d6fd28d29856908`. In this file, 23 events pass the analysis cuts before the proton tag and 13 after it. Those are physics selection counts; `result.json` counts all generated rows.

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

The launcher checks the pilot receipt, lists the reconstructed hiAcc 1S 9×275 inputs directly under the EOS sample directory, and writes `input.txt` and `plan.log` in `workspaces/condor-1S-9x275-full-grouped-v1/`. You can override the repository input base with `UPSILON_DATA`. Review `input.txt` against the intended campaign inventory before interpreting the final row count. A given input, including the pilot original, appears once in the full plan. The separate test output is excluded.

The planner records `files_per_job`, `condor_job_count` and the number of per-file work items (`job_count`) in `batch/plan.json`. The submission asks for `group_u_CMS.CAF.COMM`, `longlunch` and `eic_xl:26.07.1-stable`. Each Condor process runs the files in its group sequentially, using the same per-file ROOT invocation and weights as the successful pilot. Every input has its own output and receipt. The final merge follows the original manifest order.

If you previously ran the one-job-per-file commands, a **plan only** can coexist with the new grouped plan. An earlier submission is protected even after its jobs leave the queue. The launcher never cancels jobs. For a removed one-file-per-job cluster, follow the recovery steps below. For a new grouped submission it reserves `submission.txt` before calling Condor; if the outcome is uncertain, inspect that file and the scheduler before retrying.

### Switch after removing the one-file-per-job cluster

Keep the old batch, its submission record and its logs at their original paths. From bare lxplus, using the same CERN scheduler as the original submission, run:

```bash
git pull --ff-only origin main && bash production/submit_1s_9x275.sh --after-removal --submit
```

The launcher reads the cluster ID from the earlier `submission.txt`. If the record lacks an ID, set `OLD_CLUSTER_ID` to the ID you removed. It stops if Condor still lists the old cluster, including jobs in the removing state. In that case rerun the command once the queue clears.

Recovery checks the old and new plans describe the same input inventory and reconstruction. It carries verified `DONE` outputs into the grouped plan. Interrupted inputs run again; their previous files stay in the old batch. The launcher reports how many files it reused and submits only groups with outstanding work. A group containing both completed and outstanding files skips the verified ones. Use `--after-removal` without `--submit` to preview this recovery and prepare the grouped outputs without submitting to the scheduler.

For a custom campaign, use `production/jobs.py plan ... --files-per-job 100`. Omitting this option keeps one input per Condor process. Grouping is a scheduler setting; the underlying per-file manifest, ROOT selection and merge interface are the same. Unhashed planning records path, size and modification time; `--hash-inputs` additionally reads every original file to record and verify its SHA-256.

## Monitor the grouped campaign

On bare lxplus, in the repository root, run:

```bash
git pull --ff-only origin main
python3 production/rematch_status.py
```

This reads the grouped campaign's recorded cluster IDs and worker receipts. It shows idle, running and held **Condor processes** separately from **input files** marked DONE, FAILED, untouched or started without a final receipt. One process can handle up to 100 files. The report includes hold reasons and example file names, errors and log paths. The row total comes from the DONE receipts. The final merge verifies the ROOT contents.

Run the report periodically; it scans each file's job directory on AFS. For a lightweight live view of just this campaign, run:

```bash
watch -n 60 'python3 production/rematch_status.py --queue-only'
```

Ctrl+C stops the display. If a process is held, inspect its stated reason or use `condor_q -hold CLUSTER_ID` with the cluster ID from the report. For a file-only report without a scheduler query, run `python3 production/rematch_status.py --no-queue`.

## Retry failed or unfinished files

Use the same bare lxplus submit scheduler as the initial grouped run. Inspect the reported errors first; repeated permission or input-format errors need their cause fixed. A running group can continue processing other inputs after an individual file fails. A started/partial file can still be running. Once the campaign's queue has drained (including held, transferring and removing processes), preview what would be retried:

```bash
python3 production/retry_rematch_1s_9x275.py
```

To submit that retry:

```bash
python3 production/retry_rematch_1s_9x275.py --submit
```

The retry checks the preserved outputs, archives the evidence from failed or interrupted attempts and queues only groups containing unfinished files. These groups skip verified completed inputs, including outputs carried over from a removed cluster. Each retry has a new directory under `batch/retries/attempt-…/` containing its submission receipt and Condor logs. Run the same status command to follow it. Leave the original `submission.txt`, worker locks and results in place. If a submission receipt lacks a cluster ID, first reconcile its outcome on the scheduler. Finish any manually launched workers before using the retry command.

## Check and merge the full sample

For the original July 1S 9×275 inventory, require **all input files DONE and 1,137,600 generated rows** in the status report. Investigate a different total. An empty Condor queue alone can also mean that a launcher exited before it opened its first file.

Merge **inside the pinned ePIC shell** with access to the AFS batch and writable EOS output. Use the [lxplus container launch](rematch-3011.md#start-on-lxplus-and-run-one-file) with the same EOS credentials as the interactive test. The default output base is a separate `rematch3011-v1` directory beneath the July data area. Set `CORRECTED_BASE` to another writable absolute directory if needed. In that shell, from the repo root, run:

```bash
(
  set -euo pipefail
  batch="$PWD/workspaces/condor-1S-9x275-full-grouped-v1/batch"
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
