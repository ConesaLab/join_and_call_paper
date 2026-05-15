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
module load samtools

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
max_reads="${10:-0}"
talon_main_threads="${11:-1}"
cap_blas="${12:-0}"
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
echo "talon_max_reads (0 = full SAM): ${max_reads:-0}"
echo "talon_main_threads (--threads for main talon): ${talon_main_threads}"
echo "cap_blas (1 = export OMP_/OpenBLAS/MKL=1 before talon): ${cap_blas}"
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

# Convert BAM to SAM if SAM doesn't exist
if [ ! -f "$sam_file" ]; then
	sam_file="${out_transcriptclean}/${sample_id}.sam"
	echo "Converting BAM to SAM: $bam_file -> $sam_file"
	samtools view -h "$bam_file" > "$sam_file"
fi

# Optional smoke test: first N primary alignments only (header + N data lines)
if [[ "${max_reads}" =~ ^[0-9]+$ ]] && [ "${max_reads}" -gt 0 ]; then
	subset_sam="${out_transcriptclean}/${sample_id}.maxreads${max_reads}.sam"
	if [ -f "$subset_sam" ]; then
		echo "Reusing subset SAM (${max_reads} reads): $subset_sam"
	else
		echo "Subsetting SAM to first ${max_reads} alignments -> $subset_sam"
		{
			samtools view -H "$sam_file"
			samtools view "$sam_file" | head -n "${max_reads}"
		} > "$subset_sam" || { echo "ERROR: SAM subset failed"; exit 1; }
	fi
	sam_file="$subset_sam"
fi

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

# Fail fast if main talon crashes (avoids misleading filter/GTF/abundance errors on empty DB).
set -e

echo "running talon_label_reads ..."
echo "talon_label_reads threads: 1 (avoid multi-thread heap crashes; TranscriptClean uses n_cores=$n_cores; main talon uses talon_main_threads=$talon_main_threads)"

# run TALON steps
talon_label_reads \
	--f "${out_transcriptclean}/TC_clean.sam" \
	--g $genome \
	--tmpDir "${out_label_reads}/tmp" \
	--deleteTmp \
	--t 1 \
	--o $out_label_reads

echo "running talon ..."
echo "talon main --threads: ${talon_main_threads} (Slurm CPUs: ${n_cores}); cap_blas=${cap_blas}"

if [ "${cap_blas}" = "1" ]; then
	export OMP_NUM_THREADS=1
	export OPENBLAS_NUM_THREADS=1
	export MKL_NUM_THREADS=1
	export NUMEXPR_NUM_THREADS=1
	echo "Set OMP_/OpenBLAS/MKL/NUMEXPR thread caps to 1 for main talon only"
fi

talon \
	--f $config_file \
	--db $talon_database \
	--build $genome_name \
	--tmpDir "$(dirname $out_talon)/tmp" \
	--threads "${talon_main_threads}" \
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
