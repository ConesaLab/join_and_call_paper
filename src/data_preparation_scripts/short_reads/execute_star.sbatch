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

# Check if required arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 pairs.tsv genome_index"
    exit 1
fi

pairs_file="$1"
genome_index="$2"

# Get the directory where pairs.tsv is located
base_dir=$(dirname "$pairs_file")

# Extract the line corresponding to the SLURM_ARRAY_TASK_ID
line=$(awk "NR == ${SLURM_ARRAY_TASK_ID} { print \$0 }" $pairs_file)

# Separate the fields into $id1 and $id2 assuming tab-separated columns
id1=$(echo "$line" | awk -F'\t' '{print $1}')
id2=$(echo "$line" | awk -F'\t' '{print $2}')

# Construct full file paths
file1="${base_dir}/${id1}.fastq.gz"
file2="${base_dir}/${id2}.fastq.gz"

# Create output directory based on input file names
output_dir="$base_dir"
output_prefix="${id1}_${id2}"

# Map
STAR --runThreadN 8 \
    --genomeDir $genome_index \
    --readFilesIn $file1 $file2 \
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
