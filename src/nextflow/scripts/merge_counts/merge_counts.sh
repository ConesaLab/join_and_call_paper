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


# Discover conditions dynamically from individual metadata (column 1)
readarray -t conditions < <(cut -f 1 "$metadata_ind" | grep -v '^#' | sort -u)

# Create per-condition quant fofn files and associative arrays for TAMA merge/count
declare -A cond_quant_fofn
declare -A cond_merge
declare -A cond_out_mat
for cond in "${conditions[@]}"; do
    cond_quant_fofn[$cond]="${WD}/${cond}_quant.fofn"
    > "${cond_quant_fofn[$cond]}"
done

ind_quant_fofn="${WD}/ind_quant.fofn"
> "$ind_quant_fofn"
concat_quant_fofn="${WD}/concat_quant.fofn"
> "$concat_quant_fofn"
ind_merge=""
ind_out_mat=""
concat_merge=""
concat_out_mat=""

# Route samples from merged metadata into the appropriate fofn files
while IFS=$'\t' read -r line; do
    line="${line//$'\t'/';'}"
    IFS=';' read -r sample_cond sample_id pool bam_file fastq_file gtf count merge <<< "$line"

    if [ "$sample_id" == "concat" ]; then
        echo -e "$sample_cond\t$count" >> "$concat_quant_fofn"
    elif [ "$sample_id" != "TAMA" ]; then
        echo -e "$sample_id\t$count" >> "$ind_quant_fofn"
        for cond in "${conditions[@]}"; do
            if [ "$sample_cond" == "$cond" ]; then
                echo -e "$sample_id\t$count" >> "${cond_quant_fofn[$cond]}"
            fi
        done
    else
        if [ "$sample_cond" == "ind_TAMA" ]; then
            ind_merge=$merge
            ind_out_mat=$count
        elif [ "$sample_cond" == "concat_TAMA" ]; then
            concat_merge=$merge
            concat_out_mat=$count
        else
            for cond in "${conditions[@]}"; do
                if [ "$sample_cond" == "${cond}_TAMA" ]; then
                    cond_merge[$cond]=$merge
                    cond_out_mat[$cond]=$count
                fi
            done
        fi
    fi
done < <(grep -v '^#' "$metadata_merged")

sort -k2,2 "$ind_quant_fofn" -o "$ind_quant_fofn"
sort -k2,2 "$concat_quant_fofn" -o "$concat_quant_fofn"

# Build the merge_counts array config: one line per condition + ind + concat
> "$ind_out_mat"
> "$concat_out_mat"

merge_array_config="${WD}/merge_counts_array_config.tsv"
> "$merge_array_config"

for cond in "${conditions[@]}"; do
    if [ -z "${cond_out_mat[$cond]}" ] || [ -z "${cond_merge[$cond]}" ]; then
        echo "WARNING: No TAMA condition entry found for condition '$cond' -- skipping per-condition merge"
        continue
    fi
    > "${cond_out_mat[$cond]}"
    echo -e "${cond_quant_fofn[$cond]}\t${cond_merge[$cond]}\t${cond_out_mat[$cond]}" >> "$merge_array_config"
done
echo -e "${ind_quant_fofn}\t${ind_merge}\t${ind_out_mat}" >> "$merge_array_config"
echo -e "${concat_quant_fofn}\t${concat_merge}\t${concat_out_mat}" >> "$merge_array_config"

nTasks=$(wc -l < "$merge_array_config")
jobid=$(sbatch --wait --array=1-$nTasks "${SCRIPT_DIR}/merge_counts.sbatch" "$merge_array_config" "$SCRIPT_DIR" | awk '{print $NF}')
echo -e "MERGE_COUNTS\t${jobid}" >> $joblog