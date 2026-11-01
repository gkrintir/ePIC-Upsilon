#!/usr/bin/env bash
# Keep this package outside the framework. Only fresh generated output is written.
set -euo pipefail
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "Usage: bash run_1s.sh FRAMEWORK_ROOT SOURCE_BATCH [TARGET_LUMI_FB=1]" >&2
  exit 2
fi
package_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
framework_dir=$(cd -- "$1" && pwd)
source_batch=$2
target_lumi=${3:-1.0}
run_id="reopt-proton-9x130-L${target_lumi//./p}-$(date -u +%Y%m%dT%H%M%SZ)"
audit_dir="$framework_dir/results/generated/runs/${run_id}-support/artifacts/t_support"

python3 -B "$package_dir/run_reoptimization_9x130.py" validate \
  --framework-root "$framework_dir" --source-batch "$source_batch" \
  --states 1S --observables pt2 abs_t_via_pt2 --lumi "$target_lumi"

python3 -B "$package_dir/audit_t_support.py" \
  --framework-root "$framework_dir" \
  --input "$framework_dir/results/generated/runs/${source_batch}-1s/artifacts/inputs/9x130_normalized.root" \
  --state 1S --source-lumi 10 --projection-lumi "$target_lumi" \
  --out-dir "$audit_dir"

python3 -B "$package_dir/run_reoptimization_9x130.py" run \
  --framework-root "$framework_dir" --source-batch "$source_batch" \
  --states 1S --observables pt2 abs_t_via_pt2 --lumi "$target_lumi" \
  --run-id "$run_id" --hash-inputs

printf '\nSupport audit: %s\nReoptimization: %s\n' "$audit_dir" \
  "$framework_dir/results/generated/runs/$run_id/artifacts/reoptimization"
