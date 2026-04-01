#!/bin/bash

module load R/4.2.1
module load anaconda
source activate SQANTI3.env

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --empty_report) empty_report="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --metadata_merged) metadata_merged="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ ! -d "$WD" ]; then
    mkdir -p "$WD"
fi

WD=$(realpath "$WD")

if [ -d "$empty_report" ]; then
    cp -a "$empty_report/." "$WD/"
fi

# determine location of script
SCRIPT_DIR=$(dirname "$(realpath $0)")

data_dir="${WD}/data"

# Discover conditions dynamically from merged metadata
# Individual samples have sample_id != "concat" and != "TAMA"; strip _TAMA suffix to find base condition
readarray -t conditions < <(grep -v '^#' "$metadata_merged" | while IFS=$'\t' read -r line; do
    line="${line//$'\t'/';'}"
    IFS=';' read -r cond sample_id rest <<< "$line"
    if [[ "$sample_id" != "concat" && "$sample_id" != "TAMA" ]]; then
        echo "$cond"
    fi
done | sort -u)

# Create per-condition classification and junction directories
for cond in "${conditions[@]}"; do
    mkdir -p "${data_dir}/classification_${cond}"
    mkdir -p "${data_dir}/junctions_${cond}"
done

# Route SQANTI3 output files into per-condition directories
while IFS=$'\t' read -r line; do
    line="${line//$'\t'/';'}"
    IFS=';' read -r cond sample_id pool bam_file fastq_file gtf count merge sq3dir <<< "$line"

    if [[ "$cond" != "ind_TAMA" && "$cond" != "concat_TAMA" ]]; then
        base_cond="${cond%_TAMA}"
        class_target="${data_dir}/classification_${base_cond}"
        junc_target="${data_dir}/junctions_${base_cond}"

        id=0
        if [ "$sample_id" == "concat" ]; then
            id=1
        elif [ "$sample_id" == "TAMA" ]; then
            id=2
        else
            id=3
        fi

        cp "${sq3dir}"/*_classification.txt "${class_target}/${id}_${cond}_${sample_id}_classification.txt"
        cp "${sq3dir}"/*_junctions.txt "${junc_target}/${id}_${cond}_${sample_id}_junctions.txt"
    else
        cp "${sq3dir}"/*_classification.txt "${data_dir}/${cond}_classification.txt"
        cp "$gtf" "$data_dir"
        cp "$count" "$data_dir"
    fi
done < <(grep -v '^#' "$metadata_merged")

# Build per-condition .fofn files and a conditions_config.tsv for the R script
conditions_config="${data_dir}/conditions_config.tsv"
> "$conditions_config"

for cond in "${conditions[@]}"; do
    for dtype in "classification" "junctions"; do
        folder="${data_dir}/${dtype}_${cond}"
        fofn="${folder}.fofn"
        > "$fofn"
        for file in "$folder"/*; do
            [ -e "$file" ] && realpath "$file" >> "$fofn"
        done
        sort "$fofn" -o "$fofn"
    done
    echo -e "${cond}\t${data_dir}/classification_${cond}.fofn\t${data_dir}/junctions_${cond}.fofn" >> "$conditions_config"
done

Rscript $SCRIPT_DIR/generate_report_rdata.R "${WD}/helper_functions" \
    "${conditions_config}" \
    "${data_dir}/ind_TAMA.counts.tsv" \
    "${data_dir}/concat_TAMA.counts.tsv" \
    "${data_dir}/ind_TAMA_classification.txt" \
    "${data_dir}/concat_TAMA_classification.txt" \
    $genome \
    "${data_dir}/ind_TAMA.gtf" \
    "${data_dir}/concat_TAMA.gtf" \
    "${WD}/"