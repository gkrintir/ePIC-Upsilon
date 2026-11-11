#!/usr/bin/env bash
# Plan July hiAcc rematching campaigns on lxplus; submit only with --submit.
set -euo pipefail
die() { echo "$*" >&2; exit 1; }
submit=false
after_removal=false
run_label=
remaining=false
state=1S
beam=9x275
explicit_cell=false
while (( $# )); do
  case "$1" in
    --submit) submit=true ;;
    --after-removal) after_removal=true ;;
    --remaining) remaining=true ;;
    --state)
      [[ ${2:-} =~ ^[123]S$ ]] || die 'Use --state 1S, 2S, or 3S.'
      state=$2; explicit_cell=true; shift ;;
    --beam)
      [[ ${2:-} == 9x130 || ${2:-} == 9x275 ]] || die 'Use --beam 9x130 or 9x275.'
      beam=$2; explicit_cell=true; shift ;;
    --new-run)
      [[ -z "$run_label" && ${2:-} =~ ^[A-Za-z0-9][A-Za-z0-9_-]{0,47}$ && "$2" != v1 ]] ||
        die 'Supply one new run label (1–48 letters, digits, underscores or hyphens, starting with a letter or digit; v1 is reserved).'
      run_label=$2
      shift ;;
    *) die "Usage: $0 [--remaining | --state STATE --beam BEAM] [--after-removal | --new-run LABEL] [--submit]" ;;
  esac
  shift
done
if [[ -n "$run_label" ]] && "$after_removal"; then die 'Choose either --new-run or --after-removal.'; fi
cd "$(dirname "$0")/.."
repo=$(pwd -P)
data="${UPSILON_DATA:-$(python3 -c 'import json; print(json.load(open("configs/local.example.json"))["input_base"])')}"
corrected="${CORRECTED_BASE:-$data/rematch3011-v1}"
group_size=${FILES_PER_JOB:-100}
[[ "$group_size" =~ ^[1-9][0-9]*$ ]] || die 'FILES_PER_JOB must be a positive integer.'

if "$remaining"; then
  "$explicit_cell" && die 'Choose either --remaining or a specific state and beam.'
  [[ -n "$run_label" ]] || die '--remaining requires --new-run LABEL.'
  [[ -z ${TUPLE_OUTPUT_DIR:-} ]] || die 'Unset TUPLE_OUTPUT_DIR for --remaining. Use CORRECTED_BASE for a common EOS base.'
  options=(--new-run "$run_label")
  if "$submit"; then options+=(--submit); fi
  launch_cell() {
    local state=$1 beam=$2 full ids
    full="$repo/workspaces/condor-$state-$beam-full-grouped-eos-$run_label"
    if [[ -e "$full/submission.txt" ]]; then
      # Only a confirmed scheduler acknowledgement belongs in the "already submitted" count.
      ids=$(PYTHONPATH="$repo/production" python3 - "$full/batch" "$state" "$beam" "$corrected" "$run_label" "$group_size" <<'PY'
import json, sys
from pathlib import Path
from rematch_condor import recorded_clusters
batch = Path(sys.argv[1])
try:
    plan = json.loads((batch / "plan.json").read_text())
    if (plan["state"], plan["beam"], plan["work_dir"]) != (sys.argv[2], sys.argv[3], str(batch)):
        raise ValueError(f"Submission record and batch disagree: {batch}")
    subdir = "output_newID" if plan["state"] == "1S" else "output"
    tuples = Path(sys.argv[4]) / ("Upsilon" + plan["state"]) / plan["beam"] / "hiAcc" / subdir / "runs" / sys.argv[5] / "intermediate"
    if plan["tuple_output_dir"] != str(tuples.resolve()) or plan["files_per_job"] != int(sys.argv[6]):
        raise ValueError(f"Requested output base or grouping differs from the submitted plan: {batch}")
    print(", ".join(map(str, recorded_clusters(batch))))
except (OSError, ValueError, KeyError) as exc:
    sys.exit(str(exc))
PY
      ) || return 1
      [[ -n "$ids" ]] || { echo "Unconfirmed submission under $full; check Condor." >&2; return 1; }
      echo "$state $beam already submitted: cluster(s) $ids."
      return 0
    fi
    bash "$repo/production/submit_rematch.sh" --state "$state" --beam "$beam" "${options[@]}"
  }
  pids=() cells=() failed=()
  for state in 1S 2S 3S; do
    for beam in 9x130 9x275; do
      if [[ $state == 1S && $beam == 9x275 ]]; then continue; fi
      cells+=("$state/$beam")
      launch_cell "$state" "$beam" &
      pids+=("$!")
    done
  done
  for i in "${!pids[@]}"; do
    if ! wait "${pids[i]}"; then failed+=("${cells[i]}"); fi
  done
  if (( ${#failed[@]} )); then die "Check these campaigns: ${failed[*]}. Reuse the same label when resuming."; fi
  if "$submit"; then
    echo "All five campaigns have been submitted for $run_label."
  else
    echo "To submit the five campaigns: FILES_PER_JOB=${FILES_PER_JOB:-100} bash production/submit_rematch.sh --remaining --new-run $run_label --submit"
  fi
  exit 0
fi

if [[ $state/$beam != 1S/9x275 && -z "$run_label" ]]; then
  die 'Use --new-run LABEL for this state and beam. AFS migration is specific to 1S 9x275.'
fi
single="$repo/workspaces/condor-1S-9x275-full-v1"
grouped="$repo/workspaces/condor-1S-9x275-full-grouped-v1"
full="$repo/workspaces/condor-$state-$beam-full-grouped-eos-${run_label:-v1}"
trial="$repo/workspaces/condor-1S-9x275-check-v2"
batch="$full/batch"
subdir=output
if [[ $state == 1S ]]; then subdir=output_newID; fi
output_dir="$corrected/Upsilon$state/$beam/hiAcc/$subdir"
predecessors=("$single" "$grouped")
if [[ -n "$run_label" ]]; then
  output_dir="$output_dir/runs/$run_label"
  predecessors=()
fi
tuple_dir="${TUPLE_OUTPUT_DIR:-$output_dir/intermediate}"
[[ "$tuple_dir" == /* ]] || die 'The tuple output directory must be an absolute mounted path.'
[[ ! -e "$full/submission.txt" ]] || die "Submission already attempted; check $full/submission.txt and condor_q."

has_logs() { [[ -d "$1/batch/logs" && -n "$(find "$1/batch/logs" -maxdepth 1 -type f -print -quit)" ]]; }
has_logs "$full" && die "Found prior EOS-campaign Condor logs under $batch/logs. Check the scheduler before submitting."
old="$single"
old_clusters=()
unknown=()
shopt -s nullglob
for previous in "${predecessors[@]}"; do
  [[ ! -d "$previous/recovery.lock" ]] || die "Earlier recovery is in progress or interrupted at $previous/recovery.lock. Inspect it before switching."
  records=("$previous"/batch/retries/attempt-*/submission.txt)
  [[ ! -e "$previous/submission.txt" ]] || records+=("$previous/submission.txt")
  if (( ${#records[@]} )) || has_logs "$previous" || [[ -e "$previous/batch/after-removal.sub" ]]; then
    "$after_removal" || die "Found an earlier campaign at $previous. Once its jobs have left the queue, use --after-removal to reuse its finished tuples on EOS."
    old="$previous"  # the old grouped campaign already carried over the one-file campaign
    if (( ${#records[@]} == 0 )) && has_logs "$previous"; then unknown+=("$previous/batch/logs"); fi
    for record in "${records[@]}"; do
      ids=$(sed -nE 's/.*submitted to cluster ([0-9]+)\.[[:space:]]*$/\1/p' "$record" | sort -u)
      if [[ -z "$ids" ]]; then
        unknown+=("$record")
      else
        while read -r id; do old_clusters+=("$id"); done <<< "$ids"
      fi
    done
  fi
done
if "$after_removal"; then
  [[ -f "$old/batch/plan.json" && -f "$old/input.txt" ]] || die "Cannot locate the old plan and inventory under $old."
  if (( ${#old_clusters[@]} == 0 && ${#unknown[@]} == 0 )) && [[ ! -e "$old/batch/after-removal.sub" ]]; then
    unknown+=("$old/submission.txt")
  fi
  if (( ${#unknown[@]} )); then
    (( ${#unknown[@]} == 1 )) && [[ "${OLD_CLUSTER_ID:-}" =~ ^[1-9][0-9]*$ ]] || die "Missing cluster ID in ${unknown[*]}. Reconcile these records with the submit scheduler. For one missing ID set OLD_CLUSTER_ID."
    old_clusters+=("$OLD_CLUSTER_ID")
  fi
  echo "Reusing verified outputs from $old/batch."
fi

check_queue() {
  local active id campaign
  active=$(condor_q "$(id -un)" -af ClusterId Iwd Cmd) || die 'Could not check the Condor queue.'
  for campaign in "$full" "${predecessors[@]}"; do
    if [[ "$active" == *"/workspaces/${campaign##*/}/"* ]]; then
      die "Condor still lists jobs for $campaign. Check that campaign's queue before submitting."
    fi
  done
  if "$after_removal"; then
    for id in "${old_clusters[@]}"; do
      if grep -Eq "^[[:space:]]*${id}([[:space:]]|$)" <<< "$active"; then
        die "Old cluster $id still has queue entries; retry after removal finishes."
      fi
    done
  fi
}
if "$submit" || "$after_removal"; then check_queue; fi

plan_campaign() (
  python3 - "$trial/batch/jobs/job_000000/result.json" <<'PY'
import json, sys
from pathlib import Path
path = Path(sys.argv[1])
if not path.is_file():
    raise SystemExit(f"Missing successful Condor pilot receipt: {path}")
r = json.loads(path.read_text())
if (r.get("status"), r.get("state"), r.get("beam"), r.get("output_rows"), r.get("electron_id")) != ("DONE", "1S", "9x275", 100, "reco-rematch3011"):
    raise SystemExit("The expected 100-row 1S 9x275 Condor pilot must be DONE")
PY
  src="$data/Upsilon$state/$beam/hiAcc"
  [[ -d "$src" && -r "$src" ]] || die "Cannot read the input directory: $src"
  pattern="eSTARlight1.3.0-1.0_Upsilon${state}_${beam}_*_hiAcc_ab_run*.eicrecon.edm4eic.root"
  if [[ $state/$beam == 1S/9x275 ]]; then
    name=eSTARlight1.3.0-1.0_Upsilon1S_9x275_q2_0to0.01_hiAcc_ab_run000.0000.eicrecon.edm4eic.root
    head -c 4 "$src/$name" | grep -qx root || die "Cannot read the expected ROOT file under $src"
  fi
  mkdir -p "$full"
  mkdir "$full/planning.lock" || die "Another planner may be active. Inspect $full/planning.lock."
  inventory="$full/planning.lock/input.txt"
  trap 'rm -f -- "$inventory"; rmdir "$full/planning.lock"' EXIT
  [[ ! -e "$batch" ]] || die "An existing batch is at $batch; inspect it before retrying."
  find "$src" -maxdepth 1 -type f -name "$pattern" -print \
    | LC_ALL=C sort > "$inventory"
  n=$(wc -l < "$inventory")
  (( n > 1 )) || die "Found only $n file(s); check the EOS inventory."
  if [[ $state/$beam == 1S/9x275 ]]; then
    grep -Fxq "$src/$name" "$inventory" || die 'The input list is missing the validated pilot file.'
  fi
  if [[ -z "$run_label" && -f "$old/input.txt" ]]; then
    cmp -s "$old/input.txt" "$inventory" || die "The previous inventory at $old/input.txt differs from the current EOS list."
  fi
  if [[ -e "$full/input.txt" ]]; then
    cmp -s "$full/input.txt" "$inventory" || die "Inventory differs from $full/input.txt; inspect before retrying."
  else
    mv "$inventory" "$full/input.txt"
  fi
  echo "$state $beam: planning $n input files with up to $group_size per Condor job. Details: $full/plan.log"
  python3 production/jobs.py plan --input-list "$full/input.txt" --work-dir "$batch" \
    --state "$state" --beam "$beam" --electron-id reco-rematch3011 --femc full \
    --environment-label eic_xl:26.07.1-stable --job-flavour longlunch \
    --files-per-job "$group_size" --tuple-output-dir "$tuple_dir" --write >> "$full/plan.log"
)

# A plan made without --submit can be submitted later, with the same group size.
if [[ ! -f "$batch/submit.sub" ]]; then plan_campaign; fi

n=$(wc -l < "$full/input.txt")
groups=$(python3 - "$batch/plan.json" "$group_size" "$n" "$tuple_dir" "$state" "$beam" <<'PY'
import json, sys
from pathlib import Path
with open(sys.argv[1]) as f:
    plan = json.load(f)
width, count = map(int, sys.argv[2:4])
if (plan["state"], plan["beam"], plan["electron_id"], plan["femc"]) != (sys.argv[5], sys.argv[6], "reco-rematch3011", "full"):
    raise SystemExit("The existing plan uses a different state, beam or reconstruction.")
if plan.get("tuple_output_dir") != str(Path(sys.argv[4]).resolve()):
    raise SystemExit(f"The frozen plan writes tuples to {plan.get('tuple_output_dir')}. Check TUPLE_OUTPUT_DIR or CORRECTED_BASE.")
groups = (count + width - 1) // width
if (plan["files_per_job"], plan["job_count"], plan["condor_job_count"]) != (width, count, groups):
    raise SystemExit("The existing plan's grouping or input count differs from this request.")
print(groups)
PY
)
echo "$state $beam: $n files -> $groups Condor jobs (maximum $group_size files/job, group_u_CMS.CAF.COMM)"
echo "Per-file ROOT tuples: $tuple_dir/job_XXXXXX.root"
if [[ -n "$run_label" ]]; then
  echo "Monitor run $run_label: python3 production/rematch_status.py --batch '$batch' --queue-only"
fi
submit_file="$batch/submit.sub"
limit=$groups
if "$after_removal"; then
  # Make validation, carry-over and the guarded submission one recovery attempt.
  mkdir "$full/recovery.lock" || die "Another recovery may be in progress; inspect $full/recovery.lock."
  trap 'rmdir "$full/recovery.lock"' EXIT
  [[ ! -e "$full/submission.txt" ]] || die "Grouped submission already attempted; inspect $full/submission.txt."
  check_queue
  python3 production/reuse_removed_outputs.py "$old/batch" "$batch"
  submit_file="$batch/after-removal.sub"
  if [[ ! -s "$submit_file" ]]; then
    echo 'All files already have verified outputs in the grouped batch. Proceed to the full merge.'
    exit 0
  fi
  limit=$(grep -Ec '^[0-9]+$' "$submit_file")
  (( limit > 0 && limit <= groups )) || die 'Invalid recovery queue size.'
  echo "$limit Condor groups still have work."
fi
grep -Fxq '+AccountingGroup = "group_u_CMS.CAF.COMM"' "$submit_file" || die 'Accounting group mismatch.'
[[ $(grep -Ec '^[0-9]+$' "$submit_file") -eq $limit ]] || die 'Queue and plan sizes differ.'
if "$submit"; then
  check_queue
  ( set -C; : > "$full/submission.txt" ) || die 'Another submission attempt already owns this batch.'
  condor_submit -maxjobs "$limit" -append 'getenv = False' \
    -append '+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable"' \
    "$submit_file" | tee -a "$full/submission.txt"
  ack=$(grep -Ei '[1-9][0-9]* job\(s\) submitted to cluster [1-9][0-9]*\.' "$full/submission.txt") ||
    die "No confirmed cluster ID for $state $beam. Check the scheduler and $full/submission.txt."
  echo "$state $beam: $ack"
else
  options=(--state "$state" --beam "$beam")
  if "$after_removal"; then options+=(--after-removal); fi
  if [[ -n "$run_label" ]]; then options+=(--new-run "$run_label"); fi
  echo "Plan ready. To submit, run: FILES_PER_JOB=$group_size bash production/submit_rematch.sh ${options[*]} --submit"
fi
