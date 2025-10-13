#!/bin/bash
#SBATCH --job-name=map_fastq
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8             # Adjust CPU count as needed
#SBATCH --mem-per-cpu=100gb            # Adjust memory as needed
#SBATCH -t 1-00:00:00                 # Walltime
#SBATCH -o log/map_fastq_%A.out

# Load required modules
module load samtools
module load anaconda

source activate minimap2

# Get the full path to the FASTQ file
FASTQ_FILE="$1"

# Get the directory containing the FASTQ file; this will be our output directory.
OUTPUT_DIR=$(dirname "$FASTQ_FILE")

# Extract the base name from the FASTQ file (without extension) to name the BAM file accordingly.
base_name=$(basename "$FASTQ_FILE")
base_name="${base_name%.*}"

# Define output file names: unsorted and final sorted BAM file.
UNSORTED_BAM="${OUTPUT_DIR}/${base_name}.bam"
FINAL_BAM="${OUTPUT_DIR}/${base_name}.sorted.bam"

# Define the reference genome file (update this path as needed)
REFERENCE="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Map FASTQ to BAM using minimap2
minimap2 -ax splice:hq -uf --MD -t 8 "${REFERENCE}" "${FASTQ_FILE}" | samtools view -b -o "${UNSORTED_BAM}"

# Sort and index the BAM file using samtools
samtools sort "${UNSORTED_BAM}" -o "${FINAL_BAM}"
samtools index "${FINAL_BAM}"

# Optionally remove the unsorted BAM file to save space
rm -f "${UNSORTED_BAM}"