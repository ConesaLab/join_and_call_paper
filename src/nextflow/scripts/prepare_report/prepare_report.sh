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
class_brain="${data_dir}/classification_brain"
class_kidney="${data_dir}/classification_kidney"
junc_brain="${data_dir}/junctions_brain"
junc_kidney="${data_dir}/junctions_kidney"

for dir in "$class_brain" "$class_kidney" "$junc_brain" "$junc_kidney"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
done

# read relevant information from metadata
while IFS=$'\t' read -r line; do
    # need to switch to semicolon-delimited because multiple whitespace (e.g. tab) delimiters are treated as a single delimiter rather than empty fields
    line="${line//$'\t'/';'}"
    IFS=';' read -r cond sample_id pool bam_file fastq_file gtf count merge sq3dir <<< "$line"

    if [[ "$cond" != "ind_TAMA" && "$cond" != "concat_TAMA" ]]; then
        if [[ "$cond" == "B100K0" || "$cond" == "B100K0_TAMA" ]]; then
            class_target="$class_brain"
            junc_target="$junc_brain"
        elif [[ "$cond" == "B0K100" || "$cond" == "B0K100_TAMA" ]]; then
            class_target="$class_kidney"
            junc_target="$junc_kidney"
        fi
        id=0
        # define order of class and junc files for .fofn
        if [ "$sample_id" == "concat" ]; then
            id=1
        elif [ "$sample_id" == "TAMA" ]; then
            id=2
        else # for individual samples
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

# while IFS=$'\t' read -r quant_fofn merge out_mat; do
#     cp "$out_mat" "$data_dir"
# done < "$merge_counts_array_config"

for folder in "$class_brain" "$class_kidney" "$junc_brain" "$junc_kidney"; do
    fofn="${folder}.fofn"
    # Empty the .fofn file if it already exists
    > "$fofn"
    # List all files in the directory write their absolute paths to the .fofn file
    for file in "$folder"/*; do
        realpath "$file" >> "$fofn"
    done
    # Sort the file paths in the .fofn file alphabetically
    sort "$fofn" -o "$fofn"
done

Rscript $SCRIPT_DIR/generate_report_rdata.R "${WD}/helper_functions" \
    "${class_brain}.fofn" \
    "${class_kidney}.fofn" \
    "${junc_brain}.fofn" \
    "${junc_kidney}.fofn" \
    "${data_dir}/ind_TAMA.counts.tsv" \
    "${data_dir}/concat_TAMA.counts.tsv" \
    "${data_dir}/ind_TAMA_classification.txt" \
    "${data_dir}/concat_TAMA_classification.txt" \
    $genome \
    "${data_dir}/ind_TAMA.gtf" \
    "${data_dir}/concat_TAMA.gtf" \
    "${WD}/"
    