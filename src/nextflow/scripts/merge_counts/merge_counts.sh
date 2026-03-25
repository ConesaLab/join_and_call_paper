#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --metadata_ind) metadata_ind="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
        --metadata_tama) metadata_tama="$2"; shift ;;
        --metadata_tama_full) metadata_tama_full="$2"; shift ;;
        --joblog) joblog="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# determine location of script
SCRIPT_DIR=$(dirname "$(realpath $0)")

if [ ! -d "$WD" ]; then
    mkdir -p "$WD"
fi

WD=$(realpath "$WD")

# cd "$WD"

metadata_merged="${WD}/metadata_merged.tsv"
echo -e "# name\tsample_id\tpool\tbam_file\tfastq_file\tgtf\tcount\tmerge" > $metadata_merged

# parse paths to gtf files of all metadata files
while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq_file gtf count; do
    echo -e "$sample_cond\t$sample_id\t$pool\t$bam_file\t$fastq_file\t$gtf\t$count\t" >> $metadata_merged
done < <(grep -v '^#' "$metadata_ind")

while IFS=$'\t' read -r cond gtf count; do
    echo -e "$cond\tconcat\t\t\t\t$gtf\t$count\t" >> $metadata_merged
done < <(grep -v '^#' "$metadata_concat")

for metadata_file in $metadata_tama $metadata_tama_full; do
    while IFS=$'\t' read -r name gtf merge count; do
        echo -e "$name\tTAMA\t\t\t\t$gtf\t$count\t$merge" >> $metadata_merged
    done < <(grep -v '^#' "$metadata_file")
done


# perform twice: once for ind, once for concat;

brain_quant_fofn="${WD}/B100K0_quant.fofn"
> $brain_quant_fofn
kidney_quant_fofn="${WD}/B0K100_quant.fofn"
> $kidney_quant_fofn
ind_quant_fofn="${WD}/ind_quant.fofn"
> $ind_quant_fofn
concat_quant_fofn="${WD}/concat_quant.fofn"
> $concat_quant_fofn

brain_merge=""
kidney_merge=""
ind_merge=""
concat_merge=""

# read relevant information from metadata
while IFS=$'\t' read -r line; do
    # need to switch to semicolon-delimited because multiple whitespace (e.g. tab) delimiters are treated as a single delimiter rather than empty fields
    line="${line//$'\t'/';'}"
    IFS=';' read -r sample_cond sample_id pool bam_file fastq_file gtf count merge <<< "$line"

    if [ "$sample_id" == "concat" ]; then
        echo -e "$sample_cond\t$count" >> "$concat_quant_fofn"
    elif [ "$sample_id" != "TAMA" ]; then
        echo -e "$sample_id\t$count" >> "$ind_quant_fofn"
        if [[ $sample_id == B* ]]; then
            echo -e "$sample_id\t$count" >> "$brain_quant_fofn"
        elif [[ $sample_id == K* ]]; then
            echo -e "$sample_id\t$count" >> "$kidney_quant_fofn"
        fi
    else # sample_id == "TAMA"
        if [ "$sample_cond" == "B100K0_TAMA" ]; then
            brain_merge=$merge
            brain_out_mat=$count
        elif [ "$sample_cond" == "B0K100_TAMA" ]; then
            kidney_merge=$merge
            kidney_out_mat=$count
        elif [ "$sample_cond" == "ind_TAMA" ]; then
            ind_merge=$merge
            ind_out_mat=$count
        elif [ "$sample_cond" == "concat_TAMA" ]; then
            concat_merge=$merge
            concat_out_mat=$count
        fi
    fi
done < <(grep -v '^#' "$metadata_merged")

# Sort quant_fofn by the second column: B0K100/K31 < B100K0/B31 --> Kidney first
sort -k2,2 $ind_quant_fofn -o $ind_quant_fofn
sort -k2,2 $concat_quant_fofn -o $concat_quant_fofn

> $brain_out_mat
> $kidney_out_mat
> $ind_out_mat
> $concat_out_mat

merge_array_config="${WD}/merge_counts_array_config.tsv"
echo -e "${brain_quant_fofn}\t${brain_merge}\t${brain_out_mat}" > $merge_array_config
echo -e "${kidney_quant_fofn}\t${kidney_merge}\t${kidney_out_mat}" >> $merge_array_config
echo -e "${ind_quant_fofn}\t${ind_merge}\t${ind_out_mat}" >> $merge_array_config
echo -e "${concat_quant_fofn}\t${concat_merge}\t${concat_out_mat}" >> $merge_array_config

jobid=$(sbatch --wait --array=1-4 "${SCRIPT_DIR}/merge_counts.sbatch" $merge_array_config $SCRIPT_DIR | awk '{print $NF}')
echo -e "MERGE_COUNTS\t${jobid}" >> $joblog