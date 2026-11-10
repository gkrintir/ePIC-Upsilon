#!/usr/bin/env bash
# Plan the July 1S 9x275 full sample on lxplus; submit only with --submit.
set -euo pipefail
case "${1:-}" in
  ''|--submit) ;;
  *) echo "Usage: $0 [--submit]   (set FILES_PER_JOB before planning; default 100)" >&2; exit 2 ;;
esac
[[ $# -le 1 ]] || { echo 'Expected at most --submit.' >&2; exit 2; }
die() { echo "$*" >&2; exit 1; }
cd "$(dirname "$0")/.."
repo=$(pwd -P)
full="$repo/workspaces/condor-1S-9x275-full-grouped-v1"
old="$repo/workspaces/condor-1S-9x275-full-v1"
trial="$repo/workspaces/condor-1S-9x275-check-v2"
batch="$full/batch"
group_size=${FILES_PER_JOB:-100}
[[ "$group_size" =~ ^[1-9][0-9]*$ ]] || die 'FILES_PER_JOB must be a positive integer.'
[[ ! -e "$old/submission.txt" ]] || die "The earlier one-file-per-job batch has a submission record. Check $old/submission.txt and condor_q before starting another campaign."
[[ ! -e "$full/submission.txt" ]] || die "Submission already attempted; check $full/submission.txt and condor_q."
for previous in "$old" "$full"; do
  if [[ -d "$previous/batch/logs" && -n "$(find "$previous/batch/logs" -maxdepth 1 -type f -print -quit)" ]]; then
    die "Found previous Condor logs in $previous/batch/logs. Check the cluster before submitting."
  fi
done

if [[ "${1:-}" == --submit ]]; then
  active=$(condor_q "$(id -un)" -af Iwd Cmd) || die 'Could not check the Condor queue.'
  if [[ "$active" == *"/workspaces/condor-1S-9x275-full-"* ]]; then
    die 'A 1S 9x275 full-sample batch is already in the queue. Do not submit a second copy.'
  fi
fi

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
grep -Fxq '+AccountingGroup = "group_u_CMS.CAF.COMM"' "$batch/submit.sub" || die 'Accounting group mismatch.'
[[ $(grep -Ec '^[0-9]+$' "$batch/submit.sub") -eq $groups ]] || die 'Queue and plan sizes differ.'
echo "$n files -> $groups Condor jobs (maximum $group_size files/job, group_u_CMS.CAF.COMM)"
if [[ "${1:-}" == --submit ]]; then
  ( set -C; : > "$full/submission.txt" ) || die 'Another submission attempt already owns this batch.'
  condor_submit -maxjobs "$groups" -append 'getenv = False' \
    -append '+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable"' \
    "$batch/submit.sub" | tee -a "$full/submission.txt"
else
  echo "Plan ready. To submit, run: FILES_PER_JOB=$group_size bash production/submit_1s_9x275.sh --submit"
fi
