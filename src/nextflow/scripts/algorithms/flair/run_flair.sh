#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
		--use_sr) use_sr="$2"; shift ;;
		--stringent) stringent="$2"; shift ;;
		--sr_junctions) sr_junctions="$2"; shift ;;
		--sr_junctions_concat) sr_junctions_concat="$2"; shift ;;
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
out_quant="$WD/quantification"

# define output files for nextflow
flair_out_metadata_ind="${WD}/flair_out_metadata_ind.tsv"
flair_out_metadata_concat="${WD}/flair_out_metadata_concat.tsv"

> $flair_out_metadata_ind
> $flair_out_metadata_concat


# ######### Individual Samples #########

# number of samples without header line
nSamples=$(wc -l < $metadata_samples)
nSamples=$((nSamples - 1)) 

# run flair on individual samples, wait for completion
if [ $nSamples -ge 1 ]; then
	echo "Starting flair job..." 
	echo "SCRIPT_DIR: ${SCRIPT_DIR}"
	echo "nSamples: ${nSamples}"


    jobid=$(sbatch --wait --array=1-$nSamples $SCRIPT_DIR/execute_flair_ind.sbatch $WD $genome $annotation $metadata_samples $out_iso_detect $out_quant $use_sr $sr_junctions $sr_junctions_concat $stringent | awk '{print $NF}')
	echo -e "FLAIR_IND\t${jobid}" >> $joblog
fi

# create output metadata for individual samples
echo "writing flair ind metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $flair_out_metadata_ind
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# change counts to sqanti format
	flair_counts="${out_quant}/${cond}/${sample_id}.counts.tsv"
	sq3_counts="${out_quant}/${cond}/${sample_id}_sq3.counts.tsv"
	echo -e "pbid\tcount_fl" > "$sq3_counts"

	# Adjust the following line to properly format the first column before saving to sq3_counts
	tail -n +2 "$flair_counts" | awk -F'\t' 'BEGIN {OFS="\t"}
	{
		if (index($1, "|") > 0) {
			n = split($1, arr, "_");
			$1 = arr[1];
		} else {
			n = split($1, arr, "_");
			$1 = arr[1];
			for (i = 2; i < n; i++) {
				$1 = $1 "_" arr[i];
			}
		}
		print $1, int($2)
	}' >> "$sq3_counts"

	echo -e "${line}\t${out_iso_detect}/flair_collapse/${cond}/${sample_id}.isoforms.gtf\t${sq3_counts}" >> $flair_out_metadata_ind
done < "$metadata_samples"

# ######### Concat Samples #########

# Create array to store concatenated bed files
declare -a concat_bed_files

echo "[FLAIR] Concatenate BED files"

while IFS= read -r line; do
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	cond_dir="${out_iso_detect}/flair_align_correct/${cond}"
	fastq_fofn="${cond_dir}/fastq.fofn"
	concatenated_bed="${cond_dir}/split_files/${cond}_all_corrected.bed"

	# if using short reads, use separate corrected files for concat
	if [ "$use_sr" = true ]; then
		current_bed="${cond_dir}/${sample_id}_concat_all_corrected.bed"
	else
		current_bed="${cond_dir}/${sample_id}_all_corrected.bed"
	fi

	directory="${cond_dir}/split_files"
 
	echo "[FLAIR] $sample_id, $cond"

	if [ ! -d "$directory" ]; then
		mkdir -p "$directory"
		echo "[FLAIR] Directory created: $directory"
	fi

	# Add concatenated bed file to array
	if [[ ! " ${concat_bed_files[*]} " =~ " ${concatenated_bed} " ]]; then
		concat_bed_files+=("$concatenated_bed")
		echo -n "" > $concatenated_bed
		echo -n "" > $fastq_fofn
		echo -n "" > "${cond_dir}/split_files/${cond}_split_by_chrom_id.tsv"
	fi

	echo "$fastq_file" >> "$fastq_fofn"
	
	# Concatenate the bed file for the current sample into the condition's bed file
	echo "[FLAIR] $current_bed TO $concatenated_bed"
	if [[ -f "$current_bed" ]]; then
		cat "$current_bed" >> "$concatenated_bed"
	else
		echo "[FLAIR] ERROR: Concatenated BED file does not exist"
	fi

	echo "--------------------------------------------------------------------------"
done < "$metadata_samples"

# Split the concatenated bed files by chromosome (recommended by flair for large files)
for concatenated_bed in "${concat_bed_files[@]}"; do
	cond=$(basename "$concatenated_bed" | cut -d_ -f1)
	cond_dir="${out_iso_detect}/flair_align_correct/${cond}"
	fastq_fofn="${cond_dir}/fastq.fofn"
	samples_id_file="${out_iso_detect}/flair_align_correct/${cond}/split_files/${cond}_split_by_chrom_id.tsv"
	
	# Split by chromosome
	echo "[FLAIR] Split BED files by chrom: $cond"
	cut -f1 "$concatenated_bed" | sort | uniq | while read -r chrom; do
		output_split_bed="${cond_dir}/split_files/${cond}_${chrom}_all_corrected.bed"
		awk -v chrom="$chrom" -F'\t' '$1 == chrom' "$concatenated_bed" > "$output_split_bed"

		# when using --annotation_reliant generate, we also need to separate the reference gtf by chromosome
		output_split_gtf="${cond_dir}/split_files/$(basename "$annotation" .gtf)_${chrom}.gtf"
		awk -v chrom="$chrom" -F'\t' '$1 == chrom' "$annotation" > "$output_split_gtf"

		echo -e "${output_split_bed}\t${output_split_gtf}" >> "$samples_id_file"
	done

	# Run FLAIR on concat samples and wait for completion
	# This is done sequentially for the two conditions, could be improved by starting both conditions in an array and THEN starting the array for the chromosome files
	nChromSamples=$(wc -l < $samples_id_file)
	echo "nChromSamples: ${nChromSamples}"
	
	jobid=$(sbatch --wait --array=1-$nChromSamples $SCRIPT_DIR/execute_flair_concat.sbatch "$WD" "$genome" "$samples_id_file" "$out_iso_detect" "$cond" "$fastq_fofn" "$stringent" | awk '{print $NF}')
	echo -e "FLAIR_CONCAT\t$jobid" >> "$joblog"
done

# perform quantification for concat files
sort_script=$(realpath "$SCRIPT_DIR/../../util/sort_gtf.sh")
nConditions=$(awk 'NR > 1 && !/^#/ {print $1}' "$metadata_samples" | sort -u | wc -l)
jobid=$(sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_flair_concat_quant.sbatch "$WD" "$genome" "$metadata_samples" "$out_iso_detect" "$out_quant" "$sort_script" | awk '{print $NF}')
echo -e "FLAIR_CONCAT_QUANT\t${jobid}" >> $joblog



# create output metadata file for concat
echo "writing flair concat metadata output"
first_line=true
while IFS= read -r cond; do
	if $first_line && [[ "$cond" == \#* ]]; then
		echo -e "# cond\tgtf\tcount" >> $flair_out_metadata_concat
	fi
	first_line=false
	if [[ "$cond" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# sum up counts for sqanti format
	flair_counts="${out_quant}/${cond}/${cond}_ConcatReads.counts.tsv"
	sq3_counts="${out_quant}/${cond}/${cond}_ConcatReads_sq3.counts.tsv"
	
	# Adjust formatting of flair counts file for multiple samples to be compatible with SQANTI3 FL counts
	{
		read -r header
		# Process the header
		new_header="superPBID"
		IFS=$'\t' read -ra fields <<< "$header"
		for field in "${fields[@]:1}"; do
			new_header+="\t${field%%_*}"
		done
		echo -e "$new_header" > "$sq3_counts"
		
		# Process the remaining lines
		while IFS= read -r line; do
			IFS=$'\t' read -ra fields <<< "$line"
			# if [[ "${fields[0]}" == *"|"* ]]; then
			# 	new_line="${fields[0]%%_*}"
			# else
			# 	new_line="${fields[0]%_*}"
			# fi
			new_line="${fields[0]}"
			for field in "${fields[@]:1}"; do
				new_line+="\t$(printf "%.0f" "$field")"
			done
			echo -e "$new_line" >> "$sq3_counts"
		done
	} < "$flair_counts"

	orig_gtf="${out_iso_detect}/flair_collapse/${cond}/${cond}.join.isoforms.gtf"
	renamed_gtf="${out_iso_detect}/flair_collapse/${cond}/${cond}.join.isoforms_renamed.gtf"
	
	$SCRIPT_DIR/rename_gtf.sh "$orig_gtf" > "$renamed_gtf"
	
	echo -e "${cond}\t${renamed_gtf}\t${sq3_counts}" >> $flair_out_metadata_concat
done < <(awk '{print $1}' "$metadata_samples" | sort -u)
