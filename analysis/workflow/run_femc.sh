#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."

action="${1:-compare}"       # compare, note, or all
stat_mode="${STAT_MODE:-projected}"
projection_lumi="${PROJECTION_LUMI:-10}"
batch="$(date -u +%Y%m%dT%H%M%SZ)"

base="/eos/cms/store/group/phys_top/gkrintir/upsilon_data/Upsilon1S"
comparison_dir="results/generated/runs/femc-comparisons-${batch}/artifacts/beam_comparison"

case "$action" in
    compare|note|all) ;;
    *) echo "Usage: bash workflow/run_femc.sh {compare|note|all}" >&2; exit 2 ;;
esac

case "$stat_mode" in
    mc|projected) ;;
    *) echo "STAT_MODE must be mc or projected" >&2; exit 2 ;;
esac

input_path() {
    local beam="$1" scenario="$2" directory="output"
    if [[ "$scenario" != "full" ]]; then
        directory="output_${scenario}"
    fi
    printf '%s/%s/hiAcc/%s/merged_Upsilon_final.root\n' \
        "$base" "$beam" "$directory"
}

# Check all six inputs before starting.
for beam in 9x130 9x275; do
    for scenario in full remove_outer remove_inner; do
        input="$(input_path "$beam" "$scenario")"
        if [[ ! -r "$input" ]]; then
            echo "Missing or unreadable input: $input" >&2
            exit 1
        fi
    done
done

compare_pair() {
    local ref_beam="$1" ref_scenario="$2"
    local cand_beam="$3" cand_scenario="$4"
    local output="${comparison_dir}/${ref_beam}-${ref_scenario}_vs_${cand_beam}-${cand_scenario}"

    mkdir -p "$output"

    python3 scripts/validation/compare_beam_samples.py \
        --reference-file "$(input_path "$ref_beam" "$ref_scenario")" \
        --candidate-file "$(input_path "$cand_beam" "$cand_scenario")" \
        --reference-label "${ref_beam} GeV, fEMC ${ref_scenario}" \
        --candidate-label "${cand_beam} GeV, fEMC ${cand_scenario}" \
        --reference-electron-energy 9 \
        --candidate-electron-energy 9 \
        --tree UpsilonTree \
        --state 1S \
        --out-dir "$output"
}

if [[ "$action" == compare || "$action" == all ]]; then
    for beam in 9x130 9x275; do
        compare_pair "$beam" full "$beam" remove_outer
        compare_pair "$beam" full "$beam" remove_inner
        compare_pair "$beam" remove_outer "$beam" remove_inner
    done

    for scenario in full remove_outer remove_inner; do
        compare_pair 9x130 "$scenario" 9x275 "$scenario"
    done

    echo "Comparison output: $comparison_dir"
fi

if [[ "$action" == note || "$action" == all ]]; then
    for beam in 9x130 9x275; do
        for scenario in full remove_outer remove_inner; do
            run_id="femc-${beam}-${scenario}-${stat_mode}-L${projection_lumi}-${batch}"

            python3 workflow/run_statistics.py \
                --set "inputs.root_${beam}=$(input_path "$beam" "$scenario")" \
                --set "analysis.stat_mode=${stat_mode}" \
                --set "analysis.projection_lumi_invfb=${projection_lumi}" \
                --set analysis.count_interval=median_garwood \
                run "note_${beam}" \
                --run-id "$run_id" \
                --hash-inputs

            echo "Note output: results/generated/runs/${run_id}/artifacts/note/${beam}/"
        done
    done
fi
