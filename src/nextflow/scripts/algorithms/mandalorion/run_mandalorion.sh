#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
		--mandalorion_location) mandalorion_location="$2"; shift ;;
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
mandalorion_out_metadata_ind="${WD}/mandalorion_out_metadata_ind.tsv"
mandalorion_out_metadata_concat="${WD}/mandalorion_out_metadata_concat.tsv"

> $mandalorion_out_metadata_ind
> $mandalorion_out_metadata_concat


# ######### Individual Samples #########

# number of samples without header line
nSamples=$(wc -l < $metadata_samples)
nSamples=$((nSamples - 1)) 

# run mandalorion on individual samples, wait for completion
if [ $nSamples -ge 1 ]; then
	echo "Starting mandalorion job..." 
	echo "SCRIPT_DIR: ${SCRIPT_DIR}"
	echo "nSamples: ${nSamples}"

    jobid=$(sbatch --wait --array=1-$nSamples $SCRIPT_DIR/execute_mandalorion.sbatch $WD $genome $annotation $metadata_samples $out_iso_detect $mandalorion_location | awk '{print $NF}')
	echo -e "MANDALORION_IND\t${jobid}" >> $joblog
fi


# create output metadata for individual samples
echo "writing mandalorion ind metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $mandalorion_out_metadata_ind
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# change counts to sqanti format
	mandalorion_counts="${out_iso_detect}/${cond}/${sample_id}/Isoforms.filtered.clean.quant"
	sq3_counts="${out_iso_detect}/${cond}/${sample_id}/Isoforms.filtered.clean_sq3.tsv"
	echo -e "pbid\tcount_fl" > $sq3_counts
	awk -F'\t' '{print $1 "\t" $3}' "$mandalorion_counts" | tail -n +2 >> "$sq3_counts"

	mandalorion_gtf="${out_iso_detect}/${cond}/${sample_id}/Isoforms.filtered.clean.gtf"
	new_mandalorion_gtf="${out_iso_detect}/${cond}/${sample_id}/${sample_id}_Isoforms.filtered.clean.gtf"

	mv "$mandalorion_gtf" "$new_mandalorion_gtf"

	echo -e "${line}\t${new_mandalorion_gtf}\t${sq3_counts}" >> $mandalorion_out_metadata_ind
done < "$metadata_samples"

# ######### Concat Samples #########

# number of conditions without header line
nConditions=$(wc -l < $metadata_concat)
nConditions=$((nConditions - 1))

# run mandalorion on concat samples, wait for completion
if [ $nConditions -ge 1 ]; then
	echo "Starting mandalorion job..." 
	echo "SCRIPT_DIR: ${SCRIPT_DIR}"
	echo "nConditions: ${nConditions}"

    jobid=$(sbatch --mem 128gb --wait --array=1-$nConditions $SCRIPT_DIR/execute_mandalorion.sbatch $WD $genome $annotation $metadata_concat $out_iso_detect $mandalorion_location | awk '{print $NF}')
	echo -e "MANDALORION_CONCAT\t${jobid}" >> $joblog
		
	# concat quantification
	# sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_mandalorion.sbatch $WD $genome $annotation $metadata_concat $out_iso_detect
	# jobid=$(sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_mandalorion_concat_quant.sbatch $WD $genome $annotation $metadata_samples $metadata_concat $out_iso_detect $data_type $mandalorion_location | awk '{print $NF}')
	# echo -e "MANDALORION_CONCAT_QUANT\t${jobid}" >> $joblog
fi

# create output metadata file for concat
echo "writing mandalorion concat metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "# cond\tgtf\tcount" >> $mandalorion_out_metadata_concat
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	mandalorion_gtf="${out_iso_detect}/${cond}/${sample_id}/Isoforms.filtered.clean.gtf"
	new_mandalorion_gtf="${out_iso_detect}/${cond}/${sample_id}/${cond}_Isoforms.filtered.clean.gtf"

	mv "$mandalorion_gtf" "$new_mandalorion_gtf"

	# adjust format of counts for sqanti format
	mandalorion_counts="${out_iso_detect}/${cond}/${sample_id}/Isoforms.filtered.clean.quant"
	sq3_counts="${out_iso_detect}/${cond}/${sample_id}/Isoforms.filtered.clean.quant_sq3.tsv"
	
	echo -e "pbid\tcount_fl" > $sq3_counts
	
	# tail -n +2 "$mandalorion_counts" | awk -F'\t' 'BEGIN {OFS="\t"} {for (i=2; i<=NF; i++) $i=int($i); print}' >> "$sq3_counts"
	awk -F'\t' '{print $1 "\t" $3}' "$mandalorion_counts" | tail -n +2 >> "$sq3_counts"

	echo -e "${cond}\t${new_mandalorion_gtf}\t${sq3_counts}" >> $mandalorion_out_metadata_concat
done < "$metadata_concat"
