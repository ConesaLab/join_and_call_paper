#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --tama_path) tama_path="$2"; shift ;;
        --joblog) joblog="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

WD=$(realpath "$WD")

# determine location of script
SCRIPT_DIR=$(dirname "$(realpath $0)")

if [ ! -d "$WD" ]; then
    mkdir -p "$WD"
fi
# cd "$WD"

metadata_tama="${WD}/metadata_tama_condition.tsv"
echo -e "# name\tgtf\tmerge\tcount" > $metadata_tama

array_config="${WD}/tama_condition_config.tsv"
> $array_config

# create gtf_fofn for each condition
readarray -t conditions < <(cut -f 1 "$metadata_samples" | grep -v '^#' | sort -u)
nConditions=0
for cond in "${conditions[@]}"; do
    ((++nConditions))
    gtf_fofn="${WD}/${cond}_gtf.fofn"
    > "$gtf_fofn"
    while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq_file sample_gtf count; do
        if [ "$sample_cond" == "$cond" ]; then
            echo "$sample_gtf" >> "$gtf_fofn"
        fi
    done < <(grep -v '^#' "$metadata_samples")

    echo -e "${cond}\t${gtf_fofn}" >> $array_config
    
    result_file="${WD}/${cond}_TAMA/${cond}_TAMA"
    count_file="${result_file}.counts.tsv"
    echo -e "${cond}_TAMA\t${result_file}.gtf\t${result_file}_merge.txt\t${count_file}" >> $metadata_tama
done

# run tama as array job for each condition
jobid=$(sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_tama.sbatch $WD $array_config $tama_path | awk '{print $NF}')
echo -e "TAMA_CONDITION\t${jobid}" >> $joblog
