#!/usr/bin/env bash
# Read-only: run both audits even if one reports a failed consistency check.
set -u
if [ "$#" -lt 2 ]; then
  echo "Usage: bash Ypsilon_Normalization_Audit/run_two_energies.sh RUN_9x130 RUN_9x275 [OUTPUT_DIR] [SAMPLE_SUBDIR]" >&2
  echo "SAMPLE_SUBDIR defaults to output; use output_newID, output_keep_yneg, etc. only for matching runs." >&2
  exit 2
fi
audit_scripts=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
audit_output=${3:-results/normalization-audit}
audit_sample=${4:-output}
audit_data=/eos/cms/store/group/phys_top/gkrintir/upsilon_data/Upsilon1S
audit_python=${PYTHON:-python3}
audit_status=0
# This wrapper is for reruns prepared with the verified July 1S baseline.
# Audit historical outputs with the direct CLI and their recorded assignments.
audit_source='Published July 1S cards f62965a1612fc734dbab18e141f4809740d6b234; reference eSTARlight 939b11a24499398392d959db81c7502aeec91046; completed 100000-event runs; B_ee=0.0238' 
"$audit_python" "$audit_scripts/audit_upsilon_normalization.py" \
  --merged "$audit_data/9x130/hiAcc/$audit_sample/merged_Upsilon_final.root" \
  --run-dir "$1" --electron-energy 9 --proton-energy 130 \
  --assigned-visible-xsec-fb 119.685 --lumi-invfb 10 --br-ee 0.0238 \
  --generator-visible-xsec-fb 119.685 --generator-source "$audit_source" \
  --w-edges 20,25,30,35,40,45,50,55,60,65,70 \
  --out-dir "$audit_output/9x130" --hash-inputs || audit_status=$?
audit_reference=()
if [ -f "$audit_output/9x130/audit.json" ]; then
  audit_reference=(--reference-report "$audit_output/9x130/audit.json")
fi
"$audit_python" "$audit_scripts/audit_upsilon_normalization.py" \
  --merged "$audit_data/9x275/hiAcc/$audit_sample/merged_Upsilon_final.root" \
  --run-dir "$2" --electron-energy 9 --proton-energy 275 \
  --assigned-visible-xsec-fb 210.136 --lumi-invfb 10 --br-ee 0.0238 \
  --generator-visible-xsec-fb 210.136 --generator-source "$audit_source" \
  --w-edges 20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100 \
  --out-dir "$audit_output/9x275" --hash-inputs "${audit_reference[@]}" || audit_status=$?
echo "Reports: $audit_output/9x130/REPORT.md and $audit_output/9x275/REPORT.md"
echo "Checked July reference assignments: 119.685 fb (9x130), 210.136 fb (9x275), B_ee=0.0238. Use corrected framework run IDs."
exit "$audit_status"
