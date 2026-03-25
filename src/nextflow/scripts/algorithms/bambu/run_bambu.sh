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
bambu_out_metadata_ind="${WD}/bambu_out_metadata_ind.tsv"
bambu_out_metadata_concat="${WD}/bambu_out_metadata_concat.tsv"

> $bambu_out_metadata_ind
> $bambu_out_metadata_concat


for dir in $out_iso_detect $out_quant; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi
done

nSamples=$(( $(wc -l < "$metadata_samples") - 1 )) # calculate number of samples excluding header line
nConditions=$(( $(wc -l < "$metadata_concat") - 1 )) # calculate number of conditions excluding header line

# cd $WD

echo "[bambu] metadata_samples: $metadata_samples"
echo "[bambu] metadata_concat: $metadata_concat"
echo "[bambu] genome: $genome"
echo "[bambu] annotation: $annotation"
echo "[bambu] out_iso_detect: $out_iso_detect"
echo "[bambu] out_quant: $out_quant"
echo "[bambu] nSamples: $nSamples"
echo "[bambu] nConditions: $nConditions"




# call execute_bambu for each sample
if [ $nSamples -ge 1 ]; then
    jobid=$(sbatch --wait --array=1-$nSamples --export=SCRIPT_DIR=$SCRIPT_DIR,REFORMAT_SCRIPT="${SCRIPT_DIR}/../../sqanti3/rever_9th_gtf.py" $SCRIPT_DIR/execute_bambu.sbatch $genome $annotation $metadata_samples $out_iso_detect $out_quant | awk '{print $NF}')
	echo -e "bambu_IND\t${jobid}" >> $joblog
fi

# run bambu for each concatenated file (per condition)
if [ $nConditions -ge 1 ]; then
    # assign extra resources for concat runs
    jobid=$(sbatch --qos short -t 1-00:00:00 --mem 80gb --wait --array=1-$nConditions --export=SCRIPT_DIR=$SCRIPT_DIR,REFORMAT_SCRIPT="${SCRIPT_DIR}/../../sqanti3/rever_9th_gtf.py" $SCRIPT_DIR/execute_bambu.sbatch $genome $annotation $metadata_concat $out_iso_detect $out_quant | awk '{print $NF}')
    echo -e "bambu_CONCAT\t${jobid}" >> $joblog
	jobid=$(sbatch --qos short -t 1-00:00:00 --mem 32gb --wait --array=1-$nConditions --export=SCRIPT_DIR=$SCRIPT_DIR $SCRIPT_DIR/execute_bambu_quant.sbatch $genome $metadata_concat $metadata_samples $out_iso_detect $out_quant | awk '{print $NF}')
    echo -e "bambu_quant_CONCAT\t${jobid}" >> $joblog
fi


# create output metadata for individual samples
echo "writing bambu ind metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "$line\tgtf\tcount" >> $bambu_out_metadata_ind
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# get counts and gtf files
	bambu_counts="${out_quant}/${cond}/${sample_id}/${sample_id}_counts.tsv"	
	bambu_gtf="${out_iso_detect}/bambu/${cond}/${sample_id}/${sample_id}.gtf"

	echo -e "${line}\t${bambu_gtf}\t${bambu_counts}" >> $bambu_out_metadata_ind
done < "$metadata_samples"

# create output metadata for concat samples
echo "writing bambu concat metadata output"
first_line=true
while IFS= read -r line; do
	if $first_line && [[ "$line" == \#* ]]; then
		echo -e "# cond\tgtf\tcount" >> $bambu_out_metadata_concat
	fi
	first_line=false
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi
    
	# Split line into variables
	read -r cond sample_id pool bam_file fastq_file <<< "$line"

	# get counts and gtf files
	bambu_counts="${out_quant}/${cond}/${sample_id}/${sample_id}_counts.tsv"
	bambu_counts_final="${out_quant}/${cond}/${sample_id}/${cond}_counts.tsv"
	bambu_gtf="${out_iso_detect}/bambu/${cond}/${sample_id}/${sample_id}.gtf"
	bambu_gtf_final="${out_iso_detect}/bambu/${cond}/${sample_id}/${cond}.gtf"
	mv ${bambu_gtf} ${bambu_gtf_final}
	mv ${bambu_counts} ${bambu_counts_final}


	echo -e "${cond}\t${bambu_gtf_final}\t${bambu_counts_final}" >> $bambu_out_metadata_concat
done < "$metadata_concat"
