#!/bin/bash

# Runs SQANTI3 QC on StringTie merge output GTFs (no count files).
# Simplified version of run_sqanti3_qc.sh that reads the simpler
# metadata_stmerge_condition.tsv format (name, gtf) and uses the
# nocount SQANTI3 sbatch variant.

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_stmerge) metadata_stmerge="$2"; shift ;;
        --sqanti_path) sqanti_path="$2"; shift ;;
        --joblog) joblog="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ ! -d "$WD" ]; then
    mkdir -p "$WD"
fi

WD=$(realpath "$WD")

SCRIPT_DIR=$(dirname "$(realpath "$0")")
SQANTI_SCRIPT_DIR=$(realpath "$SCRIPT_DIR/../sqanti3")

sq3_inputs="${WD}/sq3_stmerge_inputs.tsv"
> "$sq3_inputs"

metadata_sq3="${WD}/metadata_sq3_stmerge.tsv"
echo -e "# name\tgtf\tsq3dir" > "$metadata_sq3"

while IFS=$'\t' read -r name gtf; do
    sq3dir="${gtf%.gtf}_sq3"
    echo -e "$name\t$gtf\t$sq3dir" >> "$metadata_sq3"
    echo "$gtf" >> "$sq3_inputs"
done < <(grep -v '^#' "$metadata_stmerge")

nSamples=$(wc -l < "$sq3_inputs")

if [ "$nSamples" -gt 0 ]; then
    jobid=$(sbatch --wait --array=1-$nSamples "$SQANTI_SCRIPT_DIR/execute_sqanti3_qc_nocount.sbatch" "$sq3_inputs" "$genome" "$annotation" "$sqanti_path" | awk '{print $NF}')
    echo -e "SQANTI3_STMERGE\t${jobid}" >> "$joblog"
fi
