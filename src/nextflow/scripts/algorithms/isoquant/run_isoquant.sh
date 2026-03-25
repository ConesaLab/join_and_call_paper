#!/bin/bash

fl_data=false

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
		--data_type) data_type="$2"; shift ;;
		--fl_data) fl_data="$2"; shift ;;
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

out_iso_detect="$WD/isoform_detection"

# define output files for nextflow
isoquant_out_metadata_ind="${WD}/isoquant_out_metadata_ind.tsv"
isoquant_out_metadata_concat="${WD}/isoquant_out_metadata_concat.tsv"

> $isoquant_out_metadata_ind
> $isoquant_out_metadata_concat

sort_script=$(realpath "$SCRIPT_DIR/../../util/sort_gtf.sh")

# ######### Individual Samples #########

# number of samples without header line
nSamples=$(wc -l < $metadata_samples)
nSamples=$((nSamples - 1)) 

# run isoquant on individual samples, wait for completion
if [ $nSamples -ge 1 ]; then
	echo "Starting isoquant job..." 
	echo "SCRIPT_DIR: ${SCRIPT_DIR}"
	echo "nSamples: ${nSamples}"

    jobid=$(sbatch --wait --array=1-$nSamples $SCRIPT_DIR/execute_isoquant.sbatch $WD $genome $annotation $metadata_samples $out_iso_detect $data_type $fl_data | awk '{print $NF}')
	echo -e "ISOQUANT_IND\t${jobid}" >> $joblog
fi

# create output metadata for individual samples
echo "writing isoquant ind metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $isoquant_out_metadata_ind
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# change counts to sqanti format
	isoquant_counts="${out_iso_detect}/${cond}/${sample_id}/OUT/OUT.transcript_model_counts.tsv"
	sq3_counts="${out_iso_detect}/${cond}/${sample_id}/OUT/OUT.transcript_model_counts_sq3.tsv"
	echo -e "pbid\tcount_fl" > $sq3_counts
	tail -n +2 "$isoquant_counts" | head -n -3 | awk -F'\t' 'BEGIN {OFS="\t"} {print $1, int($2)}' >> "$sq3_counts"

	# rename gtf
	isoquant_gtf="${out_iso_detect}/${cond}/${sample_id}/OUT/OUT.transcript_models.gtf"
	new_isoquant_gtf="${out_iso_detect}/${cond}/${sample_id}/${sample_id}.isoquant.gtf"
	
	echo "Running $sort_script to sort $isoquant_gtf into $new_isoquant_gtf ..."
	bash $sort_script -i $isoquant_gtf -o $new_isoquant_gtf

	echo -e "${line}\t${new_isoquant_gtf}\t${sq3_counts}" >> $isoquant_out_metadata_ind
done < "$metadata_samples"

# ######### Concat Samples #########

# number of conditions without header line
nConditions=$(wc -l < $metadata_concat)
nConditions=$((nConditions - 1)) 

# run isoquant on concat samples, wait for completion
if [ $nConditions -ge 1 ]; then
	echo "Starting isoquant job..." 
	echo "SCRIPT_DIR: ${SCRIPT_DIR}"
	echo "nConditions: ${nConditions}"

    # jobid=$(sbatch --qos medium -t 1-00:00:00 --cpus-per-task 80 --mem 240gb --wait --array=1-$nConditions $SCRIPT_DIR/execute_isoquant.sbatch $WD $genome $annotation $metadata_concat $out_iso_detect $data_type $fl_data | awk '{print $NF}')
	jobid=$(sbatch --qos medium -t 7-00:00:00 --wait --array=1-$nConditions $SCRIPT_DIR/execute_isoquant.sbatch $WD $genome $annotation $metadata_concat $out_iso_detect $data_type $fl_data | awk '{print $NF}')
	echo -e "ISOQUANT_CONCAT\t${jobid}" >> $joblog
		
	# sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_isoquant.sbatch $WD $genome $annotation $metadata_concat $out_iso_detect
	jobid=$(sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_isoquant_concat_quant.sbatch $WD $genome $annotation $metadata_samples $metadata_concat $out_iso_detect $data_type $fl_data | awk '{print $NF}')
	echo -e "ISOQUANT_CONCAT_QUANT\t${jobid}" >> $joblog
fi

# create output metadata file for concat
echo "writing isoquant concat metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "# cond\tgtf\tcount" >> $isoquant_out_metadata_concat
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# rename gtf
	isoquant_gtf="${out_iso_detect}/${cond}/${sample_id}/OUT/OUT.transcript_models.gtf"
	new_isoquant_gtf="${out_iso_detect}/${cond}/${sample_id}/${cond}.isoquant.gtf"
	
	echo "Running $sort_script to sort $isoquant_gtf into $new_isoquant_gtf ..."
	bash $sort_script -i $isoquant_gtf -o $new_isoquant_gtf

	# adjust format of counts for sqanti format
	isoquant_counts="${out_iso_detect}/${cond}/${cond}/OUT/OUT.transcript_grouped_counts.tsv"
	sq3_counts="${out_iso_detect}/${cond}/${cond}/OUT/OUT.transcript_grouped_counts_sq3.tsv"
	

	# Read the first line of isoquant_counts
	header=$(head -n 1 "$isoquant_counts")
	
	# Split the header into an array by tab
	IFS=$'\t' read -r -a fields <<< "$header"
	
	# Replace the first field with "pbid"
	fields[0]="superPBID"
	
	# Process each following field
	# for ((i = 1; i < ${#fields[@]}; i++)); do
	# 	field=${fields[i]}
	# 	fields[i]="${field:0:2}${field:3:1}"
	# done
	
	# Join the fields back into a single line separated by tabs
	new_header=$(IFS=$'\t'; echo "${fields[*]}")
	
	# Print the new header into sq3_counts
	echo -e "$new_header" >> "$sq3_counts"

	tail -n +2 "$isoquant_counts" | head -n -3 | awk -F'\t' 'BEGIN {OFS="\t"} {for (i=2; i<=NF; i++) $i=int($i); print}' >> "$sq3_counts"

	echo -e "${cond}\t${new_isoquant_gtf}\t${sq3_counts}" >> $isoquant_out_metadata_concat
done < "$metadata_concat"
