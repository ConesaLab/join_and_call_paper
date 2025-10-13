#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --metadata_ind) metadata_ind="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
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

sq3_inputs="${WD}/sq3_filter_inputs.tsv"
> $sq3_inputs

metadata_filter_ind="${WD}/metadata_filter_ind.tsv"
echo -e "# name\tsample_id\tpool\tbam_file\tfastq\tsr\tgtf\tcount" > $metadata_filter_ind
nindSamples=0
# parse paths to gtf files of all metadata files
while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq sr gtf count sq3dir; do
    echo -e "$gtf\t$sq3dir" >> $sq3_inputs
    output_prefix=$(basename "$gtf" .gtf)
    output_location=$(dirname "$gtf")/"${output_prefix}_sq3_filter"
    gtf=${output_location}/$(basename "$gtf" .gtf).filtered.gtf
    echo -e "$sample_cond\t$sample_id\t$pool\t$bam_file\t$fastq\t$sr\t$gtf\t$count" >> $metadata_filter_ind
    nindSamples=$(($nindSamples + 1))
done < <(grep -v '^#' "$metadata_ind")

metadata_filtered_concat="${WD}/metadata_filter_concat.tsv"
echo -e "# name\tsample_id\tpool\tbam_file\tfastq\tsr\tgtf\tcount" > $metadata_filtered_concat
while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq sr gtf count sq3dir; do
    echo -e "$gtf\t$sq3dir" >> $sq3_inputs
    output_prefix=$(basename "$gtf" .gtf)
    output_location=$(dirname "$gtf")/"${output_prefix}_sq3_filter"
    gtf=${output_location}/$(basename "$gtf" .gtf).filtered.gtf
    echo -e "$sample_cond\t$sample_id\t$pool\t$bam_file\t$fastq\t$sr\t$gtf\t$count" >> $metadata_filtered_concat
done < <(grep -v '^#' "$metadata_concat")


nSamples=$(wc -l < $sq3_inputs)

if [ $nSamples -gt 0 ]; then
    jobid=$(sbatch --wait --parsable --array=1-$nindSamples $SCRIPT_DIR/execute_sqanti3_filter.sbatch $sq3_inputs $sqanti_path $SCRIPT_DIR)
    echo -e "ISOSEQ_filter_IND\t$jobid" >> "$joblog"
fi

if [ $nSamples -gt 0 ]; then
    jobid=$(sbatch --wait --parsable --array=$(($nindSamples + 1))-$nSamples $SCRIPT_DIR/execute_sqanti3_filter.sbatch $sq3_inputs $sqanti_path $SCRIPT_DIR)
    echo -e "ISOSEQ_filter_CONCAT\t$jobid" >> "$joblog"
fi