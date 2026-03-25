#!/bin/bash
#SBATCH --job-name=prepare_isoseq
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=10gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/prepare_isoseq_%A.out

module load samtools

# Define the base directory
base_dir="/storage/gge/nih/PacBio_IsoSeq/merged_fl_reads"

# Loop through all subdirectories and process .bam files
find "$base_dir" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
    bam_file="$dir/$(basename "$dir").concat.bam"
    if [ -f "$bam_file" ]; then
        # Create a sorted .bam file and an index (.bai) for the .bam file
        sorted_bam_file="${bam_file%.bam}.sorted.bam"
        # samtools sort "$bam_file" -o "$sorted_bam_file"
        # samtools index "$sorted_bam_file"
        
        # Convert the .bam file to .sam with the same name
        sam_file="${sorted_bam_file%.bam}.sam"
        samtools view -h -o "$sam_file" "$sorted_bam_file"
    fi
done
