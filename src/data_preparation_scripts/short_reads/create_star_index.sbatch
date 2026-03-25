#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=log/star_index_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --qos=short
#SBATCH --time=24:00:00

source activate mapper

pwd;date

# Check if required arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 genome.fa annotation.gtf output_dir"
    exit 1
fi

genome_fa="$1"
annotation_gtf="$2"
output_dir="$3"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Generate STAR index
STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir "$output_dir" \
    --genomeFastaFiles "$genome_fa" \
    --sjdbGTFfile "$annotation_gtf" \
    --sjdbOverhang 100