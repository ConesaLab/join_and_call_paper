#!/bin/bash


bam_location="/storage/gge/nih/Nanopore/data/bams"
file_pattern="B3[1-5]_primary_aln_sorted.bam"
target_chr="SIRV5"
target_location="/storage/gge/Fabian/nih/data/ont/SIRV5"

for bam_file in "$bam_location"/$file_pattern; do
    # Check if file exists (in case no files match the pattern)
    if [[ -f "$bam_file" ]]; then
        target_file="${target_location}/$(basename "${bam_file}" .bam)_${target_chr}.bam"
        echo sbatch split_chr_from_bam.sbatch "$bam_file" "$target_chr" "$target_file"
        sbatch split_chr_from_bam.sbatch "$bam_file" "$target_chr" "$target_file"
    else
        echo "No matching files found for pattern: $file_pattern"
    fi
done
