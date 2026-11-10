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

If you previously ran the one-job-per-file commands, check `workspaces/condor-1S-9x275-full-v1/submission.txt` and `condor_q "$USER"`. A **plan only** can coexist with the new grouped plan. The launcher refuses to submit over a recorded submission or a matching live cluster. If the old batch is running, keep it separate and reconcile its completed files before any replacement. The launcher never cancels jobs. In the new grouped directory it reserves `submission.txt` before submission to prevent a second launch; if submission fails or its outcome is uncertain, inspect that file and the scheduler first.

For a custom campaign, use `production/jobs.py plan ... --files-per-job 100`. Omitting this option keeps one input per Condor process. Grouping is a scheduler setting; the underlying per-file manifest, ROOT selection and merge interface are the same. Unhashed planning records path, size and modification time; `--hash-inputs` additionally reads every original file to record and verify its SHA-256.

## Check and merge the full sample

Use `condor_q "$USER"` on lxplus to follow the submitted cluster. To count completed input files across all groups, run this from the repository root:

```bash
batch="$PWD/workspaces/condor-1S-9x275-full-grouped-v1/batch"
python3 - "$batch/plan.json" <<'PY'
import json, sys
from collections import Counter
from pathlib import Path
plan = json.loads(Path(sys.argv[1]).read_text())
counts = Counter()
rows = 0
for job in plan["jobs"]:
    receipt = Path(job).with_name("result.json")
    if receipt.is_file():
        try:
            result = json.loads(receipt.read_text())
        except json.JSONDecodeError:
            counts["unreadable receipt; recheck"] += 1
            continue
        counts[result["status"]] += 1
        if result["status"] == "DONE":
            rows += result["output_rows"]
    else:
        counts["no receipt yet"] += 1
print("Condor jobs:", plan["condor_job_count"], "  Input files:", len(plan["jobs"]), "  File results:", dict(counts), "  DONE rows:", rows)
PY
```

Every input file must be `DONE`. For a match to the original sample, the row total should be **1,137,600**. Investigate missing/failed jobs or a different row total before using physics yields. For example, the first group writes `$batch/logs/group_000000.{out,err,log}`. Its first input writes `$batch/jobs/job_000000/root.log` and `result.json`. If Condor restarts a group, its worker verifies completed outputs before skipping them. It preserves partial or failed attempts for inspection and reports a nonzero exit if any input fails.

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
