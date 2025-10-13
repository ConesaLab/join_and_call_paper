#!/bin/bash
#SBATCH --job-name=map_short_reads
#SBATCH --output=log/map_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=14:00:00
#SBATCH --array=0-0

source activate mapper

pwd;date

DEFAULT_PAIRS="pairs.tsv"
DEFAULT_TARGET="/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped"

pairs="${1:-$DEFAULT_PAIRS}"
target="${2:-$DEFAULT_TARGET}"

#Path to index
index_path="/storage/gge/Alejandro/genome_index/mouse_nih"

# Pair ends sufix
r1="_1.fq.gz"
r2="_2.fq.gz"

# Extract the line corresponding to the SLURM_ARRAY_TASK_ID
line=$(awk "NR == ${SLURM_ARRAY_TASK_ID} { print \$0 }" $pairs)

# Separate the fields into $condition and $file_basenames assuming tab-separated columns
condition=$(echo "$line" | awk -F'\t' '{print $1}')
file_basenames=$(echo "$line" | awk -F'\t' '{print $2}')

# Initialize empty arrays for modified filenames
files_r1_array=()
files_r2_array=()

# Split file_basenames by commas and append "_1.fq.gz" and "_2.fq.gz" to each
IFS=',' read -r -a paths <<< "$file_basenames"
for path in "${paths[@]}"; do
    files_r1_array+=("${path}${r1}")
    files_r2_array+=("${path}${r2}")
done

# Join the arrays into comma-separated strings
files_r1=$(IFS=,; echo "${files_r1_array[*]}")
files_r2=$(IFS=,; echo "${files_r2_array[*]}")

mkdir $target/$condition

echo $files_r1
echo $files_r2

# Map
STAR --runThreadN 8 \
	--genomeDir $index_path \
	--readFilesIn $files_r1 $files_r2 \
	--outFileNamePrefix $target/$condition/$condition \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterType BySJout \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNoverLmax 0.04 \
	--outFilterMismatchNmax 999 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--sjdbScore 1 \
	--genomeLoad NoSharedMemory \
	--outSAMtype BAM SortedByCoordinate \
	--twopassMode Basic \
	--readFilesCommand zcat \
	--quantMode GeneCounts
