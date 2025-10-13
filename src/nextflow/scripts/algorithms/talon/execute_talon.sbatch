#!/bin/bash
#SBATCH --job-name=talon
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=70gb # Per processor memory
#SBATCH -t 2-00:00:00     # Walltime
#SBATCH -o log/talon_%A_%a.out

module load anaconda

# Input
WD="$1"
genome_name="$2"
genome="$3"
annotation_name="$4"
annotation="$5"
metadata="$6"
empty_database="$7"
out_iso_detect="$8"
out_quant="$9"
line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# Change working directory
# cd $WD

# Get sample metadata
read -r cond sample_id pool bam_file fastq_file <<< "$(head -n $line_number_metadata $metadata | tail -n 1)"
sam_file="${bam_file%.bam}.sam"

out_transcriptclean="${out_iso_detect}/transcript_clean/${cond}/${sample_id}"
out_label_reads="${out_iso_detect}/talon_label_reads/${cond}/${sample_id}/${sample_id}"
out_talon="${out_iso_detect}/talon/${cond}/${sample_id}/${sample_id}"
filtered="${out_talon}_filtered.csv"
out_quant="${out_quant}/${cond}/${sample_id}/${sample_id}"

labeled_sam_file="${out_label_reads}_labeled.sam"

# Print Input
echo "Sample id: $sample_id"
echo "Condition: $cond"
echo "Pool: $pool"
echo "FASTQ reads: $fastq_file"
echo "SAM reads: $sam_file"
echo "Genome: $genome"
echo "Annotation: $annotation"
echo "Out dir transcriptclean: $out_transcriptclean"
echo "Out dir label_reads: $out_label_reads"
echo "Out dir talon: $out_talon"
echo "Out dir quant: $out_quant"
echo "N cores: $n_cores"
echo ""

# Check if output directories exist
directories=("$out_iso_detect" $(dirname "$out_quant") "$out_transcriptclean" $(dirname "$out_label_reads") $(dirname "$out_talon"))

for directory in "${directories[@]}"; do
	if [ ! -d "$directory" ]; then
		mkdir -p "$directory"
		echo "Directory created: $directory"
	else
		echo "Directory already exists: $directory"
	fi
done

# copy empty database
talon_database="${out_talon}.db"
cp -f $empty_database $talon_database

# write config file
if [ ! -d "$out_iso_detect/configs" ]; then
    mkdir -p "$out_iso_detect/configs"
fi

config_file="${out_iso_detect}/configs/${sample_id}_sample_config.csv"
echo -n "${cond},${sample_id},X,${labeled_sam_file}" > $config_file

source activate transcriptclean

echo "running TranscriptClean ..."

# run TranscriptClean
transcriptclean \
	--sam $sam_file \
	--genome $genome \
	--tmpDir "${out_transcriptclean}/tmp" \
	--deleteTmp \
	--threads $n_cores \
	--outprefix $out_transcriptclean

conda deactivate
source activate talon

echo "running talon_label_reads ..."

# run TALON steps
talon_label_reads \
	--f "${out_transcriptclean}/TC_clean.sam" \
	--g $genome \
	--tmpDir "${out_label_reads}/tmp" \
	--deleteTmp \
	--t $n_cores \
	--o $out_label_reads

echo "running talon ..."
talon \
	--f $config_file \
	--db $talon_database \
	--build $genome_name \
	--tmpDir "$(dirname $out_talon)/tmp" \
	--threads $n_cores \
	--o $out_talon

echo "running talon_filter_transcripts ..."
talon_filter_transcripts \
	--db $talon_database \
	-a $annotation_name \
	--o $filtered

echo "running talon_create_GTF ..."
talon_create_GTF \
	--db $talon_database \
	-b $genome_name \
	-a $annotation_name \
	--whitelist $filtered \
	--o $out_talon

echo "running talon_abundance ..."
talon_abundance \
	--db $talon_database \
	-b $genome_name \
	-a $annotation_name \
	--whitelist $filtered \
	--o $out_quant
