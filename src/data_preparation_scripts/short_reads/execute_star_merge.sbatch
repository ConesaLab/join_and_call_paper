#!/bin/bash
#SBATCH --job-name=map_short_reads_merge
#SBATCH --output=log/map_merge_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=14:00:00
#SBATCH --array=0-0

source activate mapper

pwd;date

# Check if required arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 pairs1.tsv pairs2.tsv genome_index"
    exit 1
fi

pairs_file1="$1"
pairs_file2="$2"
genome_index="$3"

# Get the directory where pairs files are located
base_dir1=$(dirname "$pairs_file1")
base_dir2=$(dirname "$pairs_file2")

# Extract the lines corresponding to the SLURM_ARRAY_TASK_ID
line1=$(awk "NR == ${SLURM_ARRAY_TASK_ID} { print \$0 }" $pairs_file1)
line2=$(awk "NR == ${SLURM_ARRAY_TASK_ID} { print \$0 }" $pairs_file2)

# Separate the fields into IDs assuming tab-separated columns
id1_1=$(echo "$line1" | awk -F'\t' '{print $1}')
id1_2=$(echo "$line1" | awk -F'\t' '{print $2}')
id2_1=$(echo "$line2" | awk -F'\t' '{print $1}')
id2_2=$(echo "$line2" | awk -F'\t' '{print $2}')

# Construct full file paths
file1_1="${base_dir1}/${id1_1}.fastq.gz"
file1_2="${base_dir1}/${id1_2}.fastq.gz"
file2_1="${base_dir2}/${id2_1}.fastq.gz"
file2_2="${base_dir2}/${id2_2}.fastq.gz"

# Create comma-separated lists of files
files1="${file1_1},${file2_1}"
files2="${file1_2},${file2_2}"

# Create output directory based on input file names
output_dir="$base_dir1"
output_prefix="${id1_1}_${id1_2}_${id2_1}_${id2_2}"

# Create output directory
mkdir -p "$output_dir/$output_prefix"

# Map
STAR --runThreadN 8 \
    --genomeDir $genome_index \
    --readFilesIn $files1 $files2 \
    --outFileNamePrefix $output_dir/$output_prefix/$output_prefix \
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