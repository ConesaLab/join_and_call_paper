#!/bin/bash

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
        --joblog) joblog="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

module load samtools
module load anaconda
source activate talon

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
talon_out_metadata_ind="${WD}/talon_out_metadata_ind.tsv"
talon_out_metadata_concat="${WD}/talon_out_metadata_concat.tsv"

> $talon_out_metadata_ind
> $talon_out_metadata_concat

genome_name=$(basename "${genome%.*}")
annotation_name=$(basename "${annotation%.*}")

for dir in $out_iso_detect $out_quant; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
done

nSamples=$(( $(wc -l < "$metadata_samples") - 1 )) # calculate number of samples excluding header line
nConditions=$(( $(wc -l < "$metadata_concat") - 1 )) # calculate number of conditions excluding header line

annotation_talon="$WD/annotation/$(basename "${annotation%.gtf}_talon.gtf")"

empty_database="${out_iso_detect}/${genome_name}_empty"
empty_database_file="${empty_database}.db"

# cd $WD

echo "[TALON] metadata_samples: $metadata_samples"
echo "[TALON] metadata_concat: $metadata_concat"
echo "[TALON] genome_name: $genome_name"
echo "[TALON] annotation_name: $annotation_name"
echo "[TALON] out_iso_detect: $out_iso_detect"
echo "[TALON] out_quant: $out_quant"
echo "[TALON] nSamples: $nSamples"
echo "[TALON] nConditions: $nConditions"
echo "[TALON] annotation_talon: $annotation_talon"
echo "[TALON] empty_database: $empty_database"
echo "[TALON] empty_database_file: $empty_database_file"

sort_script=$(realpath "$SCRIPT_DIR/../../util/sort_gtf.sh")

# convert BAM to SAM
for metadata_file in $metadata_samples $metadata_concat; do
    grep -v "^#" $metadata_file | while read -r condition sample pool bam fastq; do
        samfile="${bam%.bam}.sam"
        if [ ! -f "$samfile" ]; then
            echo "Converting BAM to SAM: ${bam}" 
            samtools view -h $bam > $samfile
        fi
    done
done

# reformat and clean gtf file
bash $SCRIPT_DIR/reformat_annotation.sh -a $annotation --delete_tmp --annotation_talon $annotation_talon

# create talon database from gtf
if [ ! -f "$empty_database_file" ]; then
    talon_initialize_database \
        --f ${annotation_talon} \
        --g ${genome_name} \
        --a ${annotation_name} \
        --o ${empty_database}
fi

# call execute_talon for each sample
if [ $nSamples -ge 1 ]; then
    jobid=$(sbatch --wait --array=1-$nSamples $SCRIPT_DIR/execute_talon.sbatch $WD $genome_name $genome $annotation_name $annotation_talon $metadata_samples $empty_database_file $out_iso_detect $out_quant | awk '{print $NF}')
	echo -e "TALON_IND\t${jobid}" >> $joblog
fi

# run talon for each concatenated file (per condition)
if [ $nConditions -ge 1 ]; then
    # assign extra resources for concat runs
    jobid=$(sbatch --qos medium -t 7-00:00:00 --mem 500gb --wait --array=1-$nConditions $SCRIPT_DIR/execute_talon.sbatch $WD $genome_name $genome $annotation_name $annotation_talon $metadata_concat $empty_database_file $out_iso_detect $out_quant | awk '{print $NF}')
    echo -e "TALON_CONCAT\t${jobid}" >> $joblog
fi

# create output metadata for individual samples
echo "writing TALON ind metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $talon_out_metadata_ind
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# change counts to sqanti format
	talon_counts="${out_quant}/${cond}/${sample_id}/${sample_id}_talon_abundance_filtered.tsv"
	sq3_counts="${out_quant}/${cond}/${sample_id}/${cond}_sq3.counts.tsv"
	echo -e "pbid\tcount_fl" > $sq3_counts
	tail -n +2 "$talon_counts" | awk -F'\t' 'BEGIN {OFS="\t"} {print $4, $NF}' >> "$sq3_counts"
	
	# rename talon gtf
	talon_gtf="${out_iso_detect}/talon/${cond}/${sample_id}/${sample_id}_talon.gtf"
	new_talon_gtf="${out_iso_detect}/talon/${cond}/${sample_id}/${sample_id}.talon.sorted.gtf"

	echo "Running $sort_script to sort $talon_gtf into $new_talon_gtf ..."
	bash $sort_script -i $talon_gtf -o $new_talon_gtf

	echo -e "${line}\t${new_talon_gtf}\t${sq3_counts}" >> $talon_out_metadata_ind
done < "$metadata_samples"

# create output metadata for concat samples
echo "writing TALON concat metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "# cond\tgtf\tcount" >> $talon_out_metadata_concat
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi
    
	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# sum up counts for sqanti format
	talon_counts="${out_quant}/${cond}/${sample_id}/${sample_id}_talon_abundance_filtered.tsv"
	sq3_counts="${out_quant}/${cond}/${sample_id}/${cond}_sq3.counts.tsv"
	echo -e "pbid\tcount_fl" > $sq3_counts
	tail -n +2 "$talon_counts" | awk -F'\t' 'BEGIN {OFS="\t"} {print $4, $NF}' >> "$sq3_counts"

	# rename talon gtf
	talon_gtf="${out_iso_detect}/talon/${cond}/${sample_id}/${sample_id}_talon.gtf"
	new_talon_gtf="${out_iso_detect}/talon/${cond}/${sample_id}/${cond}.talon.sorted.gtf"
	
	echo "Running $sort_script to sort $talon_gtf into $new_talon_gtf ..."
	bash $sort_script -i $talon_gtf -o $new_talon_gtf

	echo -e "${cond}\t${new_talon_gtf}\t${sq3_counts}" >> $talon_out_metadata_concat
done < "$metadata_concat"
