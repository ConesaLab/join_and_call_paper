#!/bin/bash

# ./run_talon_experiment.sh --wd /home/fabianje/repos/documenting_NIH/fabian/data/output/ont_subset/talon_experiment --genome /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa --annotation /home/fabianje/repos/documenting_NIH/fabian/data/output/ont/talon_run_full_1/talon/annotation/mm39.ncbiRefSeq_SIRV_talon.gtf --metadata_samples /storage/gge/Fabian/nih/data/metadata/ont_samples_brain_subset.tsv --metadata_concat /storage/gge/Fabian/nih/data/metadata/ont_concat_samples_brain_subset.tsv

# Input
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --wd) WD="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --annotation) annotation="$2"; shift ;;
        --metadata_samples) metadata_samples="$2"; shift ;;
        --metadata_concat) metadata_concat="$2"; shift ;;
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

annotation_talon="$annotation"

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
    sbatch --wait $SCRIPT_DIR/talon_ind_experiment.sbatch $WD $genome_name $genome $annotation_name $annotation_talon $metadata_samples $empty_database_file $out_iso_detect $out_quant B100K0
fi

# run talon for each concatenated file (per condition)
if [ $nConditions -ge 1 ]; then
    # assign extra resources for concat runs
    sbatch --wait $SCRIPT_DIR/talon_join_experiment.sbatch $WD $genome_name $genome $annotation_name $annotation_talon $metadata_concat $empty_database_file $out_iso_detect $out_quant
    # sbatch --wait --array=1-$nConditions $SCRIPT_DIR/execute_talon.sbatch $WD $genome_name $genome $annotation_name $annotation_talon $metadata_concat $empty_database_file $out_iso_detect $out_quant
fi
