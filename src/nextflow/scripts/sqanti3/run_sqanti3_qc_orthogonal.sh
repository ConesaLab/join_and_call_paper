#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_ind) metadata_ind="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
        --cage) cage="$2"; shift ;;
        --polyA) polyA="$2"; shift ;;
        --sqanti_path) sqanti_path="$2"; shift ;;
        --force_id_ignore) force_id_ignore="$2"; shift ;;
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

metadata_orth_ind="${WD}/metadata_orth_ind.tsv"
echo -e "# name\tsample_id\tpool\tbam_file\tfastq\tsr\tgtf\tcount\tsq3dir" > $metadata_orth_ind

# parse paths to gtf files of all metadata files
nindSamples=0
while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq sr gtf count; do
    echo -e "$gtf\t$count\t$sr" >> $sq3_inputs
    sq3dir="${gtf%.gtf}_sq3_orth"
    echo -e "$sample_cond\t$sample_id\t$pool\t$bam_file\t$fastq\t$sr\t$gtf\t$count\t$sq3dir" >> $metadata_orth_ind
    nindSamples=$(($nindSamples + 1))
done < <(grep -v '^#' "$metadata_ind")

metadata_orth_concat="${WD}/metadata_orth_concat.tsv"
echo -e "# name\tsample_id\tpool\tbam_file\tfastq\tsr\tgtf\tcount\tsq3dir" > $metadata_orth_concat
while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq sr gtf count; do
    echo -e "$gtf\t$count\t$sr" >> $sq3_inputs
    sq3dir="${gtf%.gtf}_sq3_orth"
    echo -e "$sample_cond\t$sample_id\t$pool\t$bam_file\t$fastq\t$sr\t$gtf\t$count\t$sq3dir" >> $metadata_orth_concat
done < <(grep -v '^#' "$metadata_concat")

nSamples=$(wc -l < $sq3_inputs)

if [ $nSamples -gt 0 ]; then
    jobid=$(sbatch --wait --parsable --array=1-$nindSamples $SCRIPT_DIR/execute_sqanti3_qc_orthogonal.sbatch $sq3_inputs $genome $annotation $sqanti_path $force_id_ignore $cage $polyA)
    echo -e "ISOSEQ_SQANTI_IND\t$jobid" >> "$joblog"
fi

if [ $nSamples -gt 0 ]; then
    jobid=$(sbatch --wait --parsable --array=$(($nindSamples + 1))-$nSamples $SCRIPT_DIR/execute_sqanti3_qc_orthogonal.sbatch $sq3_inputs $genome $annotation $sqanti_path $force_id_ignore $cage $polyA)
    echo -e "ISOSEQ_SQANTI_CONCAT\t$jobid" >> "$joblog"
fi