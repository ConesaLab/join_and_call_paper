#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
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

metadata_tama="${WD}/metadata_tama_full.tsv"
echo -e "# name\tgtf\tmerge\tcount" > $metadata_tama

array_config="${WD}/tama_full_config.tsv"
> $array_config

# create gtf_fofn for individual samples
ind_gtf_fofn="${WD}/ind_gtf.fofn"
> "$ind_gtf_fofn"
while IFS=$'\t' read -r sample_cond sample_id pool bam_file fastq_file sample_gtf count; do
    echo "$sample_gtf" >> "$ind_gtf_fofn"
done < <(grep -v '^#' "$metadata_samples")
echo -e "ind\t${ind_gtf_fofn}" >> $array_config

result_file="${WD}/ind_TAMA/ind_TAMA"
count_file="${result_file}.counts.tsv"
echo -e "ind_TAMA\t${result_file}.gtf\t${result_file}_merge.txt\t${count_file}" >> $metadata_tama

# create gtf_fofn for concat
concat_gtf_fofn="${WD}/concat_gtf.fofn"
> "$concat_gtf_fofn"
while IFS=$'\t' read -r sample_cond sample_gtf count; do
    echo "$sample_gtf" >> "$concat_gtf_fofn"
done < <(grep -v '^#' "$metadata_concat")
echo -e "concat\t${concat_gtf_fofn}" >> $array_config

result_file="${WD}/concat_TAMA/concat_TAMA"
count_file="${result_file}.counts.tsv"
echo -e "concat_TAMA\t${result_file}.gtf\t${result_file}_merge.txt\t${count_file}" >> $metadata_tama

# run tama as array job for individuals and concat
jobid=$(sbatch --wait --array=1-2 $SCRIPT_DIR/execute_tama.sbatch $WD $array_config $tama_path | awk '{print $NF}')
echo -e "TAMA_FULL\t${jobid}" >> $joblog
