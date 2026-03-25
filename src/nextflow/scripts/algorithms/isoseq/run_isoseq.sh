#!/bin/bash

# Input
echo $@
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
		--joblog) joblog="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

module load samtools
module load anaconda
source activate isoseq.4.0

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
isoseq_out_metadata_ind="${WD}/isoseq_out_metadata_ind.tsv"
isoseq_out_metadata_concat="${WD}/isoseq_out_metadata_concat.tsv"

> $isoseq_out_metadata_ind
> $isoseq_out_metadata_concat

genome_name=$(basename "${genome%.*}")

for dir in $out_iso_detect $out_quant; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
done
head $metadata_concat
nSamples=$(( $(wc -l < "$metadata_samples") - 1 )) # calculate number of samples excluding header line
nConditions=$(( $(wc -l < "$metadata_concat") - 1 )) # calculate number of conditions excluding header line


# cd $WD

echo "[Iso-Seq] metadata_samples: $metadata_samples"
echo "[Iso-Seq] metadata_concat: $metadata_concat"
echo "[Iso-Seq] genome_name: $genome_name"
echo "[Iso-Seq] out_iso_detect: $out_iso_detect"
echo "[Iso-Seq] out_quant: $out_quant"
echo "[Iso-Seq] nSamples: $nSamples"
echo "[Iso-Seq] nConditions: $nConditions"


# call execute_isoseq for each sample
if [ $nSamples -ge 1 ]; then
    jobid=$(sbatch --parsable --mem 45gb --wait --array=1-$nSamples $SCRIPT_DIR/execute_isoseq.sbatch $WD $genome $metadata_samples $out_iso_detect $out_quant)
	echo -e "ISOSEQ_IND\t${jobid}" >> $joblog
fi

# run isoseq for each concatenated file (per condition)
if [ $nConditions -ge 1 ]; then
    # assign extra resources for concat runs
    jobid=$(sbatch --parsable --mem 45gb --wait --array=1-$nConditions $SCRIPT_DIR/execute_isoseq.sbatch $WD $genome $metadata_concat $out_iso_detect $out_quant)
    echo -e "ISOSEQ_CONCAT\t$jobid" >> "$joblog"
fi

# create output metadata for individual samples
echo "writing Iso-Seq ind metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $isoseq_out_metadata_ind
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fq sr <<< "$line"

	# add sqanti3 fl counts to metadata
	sq3_counts="${out_quant}/isoseq/${cond}/${sample_id}/${sample_id}.flnc_count.txt"
	
	# add isoseq gtf to metadata
	isoseq_gtf="${out_iso_detect}/isoseq/${cond}/${sample_id}/collapse/${sample_id}.gtf"

    echo -e "${line}\t${isoseq_gtf}\t${sq3_counts}" >> $isoseq_out_metadata_ind
done < "$metadata_samples"

# create output metadata for concat samples
echo "writing Iso-Seq concat metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $isoseq_out_metadata_concat
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fq sr <<< "$line"

	# add sqanti3 fl counts to metadata
	sq3_counts="${out_quant}/isoseq/${cond}/${sample_id}/${sample_id}.flnc_count.txt"
	
	# add isoseq gtf to metadata
	isoseq_gtf="${out_iso_detect}/isoseq/${cond}/${sample_id}/collapse/${sample_id}.gtf"

	echo -e "${line}\t${isoseq_gtf}\t${sq3_counts}" >> $isoseq_out_metadata_concat
done < "$metadata_concat"