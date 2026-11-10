#!/usr/bin/env bash
# Plan the July 1S 9x275 full sample on lxplus; submit only with --submit.
set -euo pipefail
die() { echo "$*" >&2; exit 1; }
submit=false
after_removal=false
for option in "$@"; do
  case "$option" in
    --submit) submit=true ;;
    --after-removal) after_removal=true ;;
    *) die "Usage: $0 [--after-removal] [--submit] (FILES_PER_JOB defaults to 100)" ;;
  esac
done
cd "$(dirname "$0")/.."
repo=$(pwd -P)
full="$repo/workspaces/condor-1S-9x275-full-grouped-v1"
old="$repo/workspaces/condor-1S-9x275-full-v1"
trial="$repo/workspaces/condor-1S-9x275-check-v2"
batch="$full/batch"
group_size=${FILES_PER_JOB:-100}
[[ "$group_size" =~ ^[1-9][0-9]*$ ]] || die 'FILES_PER_JOB must be a positive integer.'
if ! "$after_removal" && [[ -e "$old/submission.txt" ]]; then
  die "The earlier one-file-per-job batch has a submission record. After removing that cluster, use --after-removal to verify and resume."
fi
[[ ! -e "$full/submission.txt" ]] || die "Submission already attempted; check $full/submission.txt and condor_q."
for previous in "$old" "$full"; do
  if "$after_removal" && [[ "$previous" == "$old" ]]; then continue; fi
  if [[ -d "$previous/batch/logs" && -n "$(find "$previous/batch/logs" -maxdepth 1 -type f -print -quit)" ]]; then
    die "Found previous Condor logs in $previous/batch/logs. Check the cluster before submitting."
  fi
done

if "$after_removal"; then
  [[ -f "$old/batch/plan.json" && -f "$old/input.txt" ]] || die "Cannot locate the old plan and inventory under $old."
  old_clusters=""
  if [[ -f "$old/submission.txt" ]]; then
    old_clusters=$(sed -nE 's/.*submitted to cluster ([0-9]+)\.[[:space:]]*$/\1/p' "$old/submission.txt" | sort -u)
  fi
  if [[ -z "$old_clusters" ]]; then
    [[ "${OLD_CLUSTER_ID:-}" =~ ^[0-9]+$ ]] || die "No cluster ID in $old/submission.txt. Set OLD_CLUSTER_ID to the ID you removed."
    old_clusters=$OLD_CLUSTER_ID
  fi
fi

check_queue() {
  local active id
  active=$(condor_q "$(id -un)" -af ClusterId Iwd Cmd) || die 'Could not check the Condor queue.'
  if [[ "$active" == *"/workspaces/condor-1S-9x275-full-"* ]]; then
    die 'A 1S 9x275 full-sample batch is still in the queue. If already removed, retry after Condor finishes cleanup.'
  fi
  if "$after_removal"; then
    for id in $old_clusters; do
      if grep -Eq "^[[:space:]]*${id}([[:space:]]|$)" <<< "$active"; then
        die "Old cluster $id still has queue entries; retry after removal finishes."
      fi
    done
  fi
}
if "$submit" || "$after_removal"; then check_queue; fi

# A plan made without --submit can be submitted later, with the same group size.
if [[ ! -f "$batch/submit.sub" ]]; then
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
  data="${UPSILON_DATA:-$(python3 -c 'import json; print(json.load(open("configs/local.example.json"))["input_base"])')}"
  src="$data/Upsilon1S/9x275/hiAcc"
  name=eSTARlight1.3.0-1.0_Upsilon1S_9x275_q2_0to0.01_hiAcc_ab_run000.0000.eicrecon.edm4eic.root
  head -c 4 "$src/$name" | grep -qx root || die "Cannot read the expected ROOT file under $src"
  mkdir -p "$repo/workspaces"
  mkdir "$full" || die "An incomplete plan may exist at $full; inspect it before retrying."
  find "$src" -maxdepth 1 -type f -name 'eSTARlight1.3.0-1.0_Upsilon1S_9x275_*_hiAcc_ab_run*.eicrecon.edm4eic.root' -print \
    | LC_ALL=C sort > "$full/input.txt"
  n=$(wc -l < "$full/input.txt")
  (( n > 1 )) || die "Found only $n file(s); check the EOS inventory."
  grep -Fxq "$src/$name" "$full/input.txt" || die 'The input list is missing the validated pilot file.'
  if [[ -f "$old/input.txt" ]]; then
    cmp -s "$old/input.txt" "$full/input.txt" || die "The previous inventory and the current EOS list differ. Inspect $full/input.txt."
  fi
  echo "Planning $n input files with up to $group_size per Condor job. Details: $full/plan.log"
  python3 production/jobs.py plan --input-list "$full/input.txt" --work-dir "$batch" \
    --state 1S --beam 9x275 --electron-id reco-rematch3011 --femc full \
    --environment-label eic_xl:26.07.1-stable --job-flavour longlunch \
    --files-per-job "$group_size" --write > "$full/plan.log"
fi

n=$(wc -l < "$full/input.txt")
groups=$(python3 - "$batch/plan.json" "$group_size" "$n" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    plan = json.load(f)
width, count = map(int, sys.argv[2:])
groups = (count + width - 1) // width
if (plan["files_per_job"], plan["job_count"], plan["condor_job_count"]) != (width, count, groups):
    raise SystemExit("The existing plan's grouping or input count differs from this request.")
print(groups)
PY
)
echo "$n files -> $groups Condor jobs (maximum $group_size files/job, group_u_CMS.CAF.COMM)"
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
  if "$after_removal"; then check_queue; fi
  ( set -C; : > "$full/submission.txt" ) || die 'Another submission attempt already owns this batch.'
  condor_submit -maxjobs "$limit" -append 'getenv = False' \
    -append '+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable"' \
    "$submit_file" | tee -a "$full/submission.txt"
else
  options=()
  if "$after_removal"; then options+=(--after-removal); fi
  echo "Plan ready. To submit, run: FILES_PER_JOB=$group_size bash production/submit_1s_9x275.sh ${options[*]} --submit"
fi
