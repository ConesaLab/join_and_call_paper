#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_merged) metadata_merged="$2"; shift ;;
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

SCRIPT_DIR=$(dirname "$(realpath $0)")

sq3_inputs="${WD}/sq3_inputs.tsv"
> $sq3_inputs

metadata_sq3="${WD}/metadata_sq3.tsv"
echo -e "# name\tsample_id\tpool\tbam_file\tfastq_file\tgtf\tcount\tmerge\tsq3dir" > $metadata_sq3

while IFS=$'\t' read -r line; do
    # need to switch to semicolon-delimited because multiple whitespace (e.g. tab) delimiters are treated as a single delimiter rather than empty fields
    line="${line//$'\t'/';'}"
    IFS=';' read -r sample_cond sample_id pool bam_file fastq_file gtf count merge <<< "$line"

    sq3dir="${gtf%.gtf}_sq3"
    echo -e "$sample_cond\t$sample_id\t$pool\t$bam_file\t$fastq_file\t$gtf\t$count\t$merge\t$sq3dir" >> $metadata_sq3
    echo -e "$gtf\t$count" >> $sq3_inputs
done < <(grep -v '^#' "$metadata_merged")

nSamples=$(wc -l < $sq3_inputs)

if [ $nSamples -gt 0 ]; then
    jobid=$(sbatch --wait --array=1-$nSamples $SCRIPT_DIR/execute_sqanti3_qc.sbatch $sq3_inputs $genome $annotation $sqanti_path $SCRIPT_DIR | awk '{print $NF}')
    echo -e "SQANTI3\t${jobid}" >> $joblog
fi
