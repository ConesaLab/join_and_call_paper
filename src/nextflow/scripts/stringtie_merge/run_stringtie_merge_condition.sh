#!/bin/bash

# Runs StringTie merge per condition (tissue), grouping individual sample GTFs
# by the condition column (column 1) of the metadata_samples file.
# Mirrors the logic of run_tama_condition.sh but uses StringTie merge instead.

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --joblog) joblog="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

WD=$(realpath "$WD")
SCRIPT_DIR=$(dirname "$(realpath "$0")")

if [ ! -d "$WD" ]; then
    mkdir -p "$WD"
fi

metadata_stmerge="${WD}/metadata_stmerge_condition.tsv"
echo -e "# name\tgtf" > "$metadata_stmerge"

array_config="${WD}/stmerge_condition_config.tsv"
> "$array_config"

readarray -t conditions < <(cut -f 1 "$metadata_samples" | grep -v '^#' | sort -u)
nConditions=0
for cond in "${conditions[@]}"; do
    ((++nConditions))
    gtf_fofn="${WD}/${cond}_gtf.fofn"
    > "$gtf_fofn"
    while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq_file sample_gtf count; do
        if [ "$sample_cond" == "$cond" ]; then
            expr_gtf="${sample_gtf%.gtf}.expr.gtf"
            if [ -f "$expr_gtf" ]; then
                echo "$expr_gtf" >> "$gtf_fofn"
            elif [ -f "$sample_gtf" ]; then
                echo "WARNING: .expr.gtf not found for $sample_id, falling back to .gtf" >&2
                echo "$sample_gtf" >> "$gtf_fofn"
            else
                echo "ERROR: neither .expr.gtf nor .gtf found for $sample_id: $sample_gtf" >&2
                exit 1
            fi
        fi
    done < <(grep -v '^#' "$metadata_samples")

    echo -e "${cond}\t${gtf_fofn}" >> "$array_config"

    result_gtf="${WD}/${cond}_STMERGE/${cond}_STMERGE.gtf"
    echo -e "${cond}_STMERGE\t${result_gtf}" >> "$metadata_stmerge"
done

jobid=$(sbatch --wait --array=1-$nConditions "$SCRIPT_DIR/execute_stringtie_merge.sbatch" "$WD" "$array_config" | awk '{print $NF}')
echo -e "STRINGTIE_MERGE_CONDITION\t${jobid}" >> "$joblog"
