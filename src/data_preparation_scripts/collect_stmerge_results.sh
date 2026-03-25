#!/bin/bash
# Checks that all StringTie merge + SQANTI3 results exist, then collects the
# classification and junctions files into a flat directory structure for plotting:
#
#   stringtie_merge_results/
#   ├── pacbio/{tool}/{B100K0,B0K100}_STMERGE_{classification.txt,junctions.txt}
#   └── ont/{tool}/{B100K0,B0K100}_STMERGE_{classification.txt,junctions.txt}
#
# Usage (on cluster):
#   bash collect_stmerge_results.sh [--dry-run]

set -euo pipefail

PROJ_ROOT="/storage/gge/home_members/fabianje/repos/join_and_call_paper"
DATA_OUT="${PROJ_ROOT}/data/output"
CONDITIONS=(B100K0 B0K100)
COLLECT_DIR="${DATA_OUT}/stringtie_merge_results"
DRY_RUN=false

if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "=== DRY RUN (check only, no copies) ==="
fi

declare -A TOOL_MAP=(
    ["isoseq/flair_ar_sr/stringtie_merge"]="pacbio/flair_ar_sr"
    ["isoseq/isoquant/stringtie_merge"]="pacbio/isoquant"
    ["isoseq/bambu/stringtie_merge"]="pacbio/bambu"
    ["isoseq/talon/stringtie_merge"]="pacbio/talon"
    ["isoseq/isoseq_pipeline/stringtie_merge"]="pacbio/isoseq_pipeline"
    ["isoseq/mandalorion/stringtie_merge"]="pacbio/mandalorion"
    ["ont/flair_ar_sr/stringtie_merge"]="ont/flair_ar_sr"
    ["ont/isoquant/stringtie_merge"]="ont/isoquant"
    ["ont/bambu/stringtie_merge"]="ont/bambu"
    ["ont/talon/stringtie_merge"]="ont/talon"
)

all_ok=true
total=0
found=0

echo ""
echo "=== Checking StringTie merge SQANTI3 results ==="
echo ""

for result_name in "${!TOOL_MAP[@]}"; do
    dest_subdir="${TOOL_MAP[$result_name]}"
    stmerge_data="${DATA_OUT}/${result_name}_data"
    stmerge_link="${stmerge_data}/stringtie_merge_condition"

    echo "--- ${dest_subdir} ---"

    if [[ ! -d "$stmerge_link" ]]; then
        echo "  MISSING: stringtie_merge_condition dir not found at ${stmerge_link}"
        all_ok=false
        continue
    fi

    for cond in "${CONDITIONS[@]}"; do
        total=$((total + 1))
        sq3_dir="${stmerge_link}/${cond}_STMERGE/${cond}_STMERGE_sq3"
        class_file="${sq3_dir}/${cond}_STMERGE_classification.txt"
        junc_file="${sq3_dir}/${cond}_STMERGE_junctions.txt"

        if [[ -f "$class_file" ]]; then
            found=$((found + 1))
            echo "  OK: ${cond}_STMERGE_classification.txt ($(wc -l < "$class_file") lines)"
        else
            echo "  MISSING: ${class_file}"
            all_ok=false
        fi
        if [[ -f "$junc_file" ]]; then
            echo "  OK: ${cond}_STMERGE_junctions.txt ($(wc -l < "$junc_file") lines)"
        else
            echo "  MISSING: ${junc_file}"
            all_ok=false
        fi
    done
done

echo ""
echo "=== Summary: ${found}/${total} classification files found ==="

if ! $all_ok; then
    echo "WARNING: Some results are missing. Check the output above."
    if ! $DRY_RUN; then
        echo "Proceeding with available results only."
    fi
fi

if $DRY_RUN; then
    echo ""
    echo "Dry run complete. Re-run without --dry-run to collect files."
    exit 0
fi

echo ""
echo "=== Collecting results into ${COLLECT_DIR} ==="

for result_name in "${!TOOL_MAP[@]}"; do
    dest_subdir="${TOOL_MAP[$result_name]}"
    stmerge_link="${DATA_OUT}/${result_name}_data/stringtie_merge_condition"
    dest_dir="${COLLECT_DIR}/${dest_subdir}"

    mkdir -p "$dest_dir"

    for cond in "${CONDITIONS[@]}"; do
        sq3_dir="${stmerge_link}/${cond}_STMERGE/${cond}_STMERGE_sq3"
        for suffix in classification.txt junctions.txt; do
            src_file="${sq3_dir}/${cond}_STMERGE_${suffix}"
            if [[ -f "$src_file" ]]; then
                cp "$src_file" "$dest_dir/"
                echo "  Copied: ${dest_subdir}/${cond}_STMERGE_${suffix}"
            fi
        done
    done
done

echo ""
echo "=== Done. Results collected at: ${COLLECT_DIR} ==="
echo "Download with:"
echo "  rsync -avz fabianje@garnatxa.uv.es:${COLLECT_DIR}/ /local/path/stringtie_merge_results/"
