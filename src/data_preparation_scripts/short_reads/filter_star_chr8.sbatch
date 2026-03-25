#!/bin/bash
#SBATCH --job-name=filter_chr8
#SBATCH --output=log/filter_chr8_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --array=0-0

module load samtools

pwd;date

# Check if required argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 paths.txt"
    exit 1
fi

paths_file="$1"

# Get the line corresponding to the SLURM_ARRAY_TASK_ID
input_path=$(awk "NR == ${SLURM_ARRAY_TASK_ID} { print \$0 }" $paths_file)

# Extract the directory name (last component of the path)
dir_name=$(basename "$input_path")

# Create chr8 subdirectory
chr8_dir="${input_path}/chr8"
mkdir -p "$chr8_dir"

# Find and process BAM file
bam_file="${input_path}/${dir_name}Aligned.sortedByCoord.out.bam"
if [ -f "$bam_file" ]; then
    # Index the input BAM file if index doesn't exist
    if [ ! -f "${bam_file}.bai" ]; then
        samtools index "$bam_file"
    fi
    
    output_bam="${chr8_dir}/${dir_name}Aligned.sortedByCoord.out_chr8.bam"
    samtools view -b "$bam_file" chr8 > "$output_bam"
    samtools index "$output_bam"
else
    echo "BAM file not found: $bam_file"
fi

# Find and process SJ file
sj_file="${input_path}/${dir_name}SJ.out.tab"
if [ -f "$sj_file" ]; then
    output_sj="${chr8_dir}/${dir_name}SJ.out_chr8.tab"
    awk '$1 == "chr8"' "$sj_file" > "$output_sj"
else
    echo "SJ file not found: $sj_file"
fi

echo "Processing completed for $dir_name"
