#!/bin/bash
#SBATCH --job-name=concat_samples_isoseq
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=40gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o concat_samples_isoseq_%A_%a.out

condition=$1

fastq_files=(/storage/gge/nih/PacBio_IsoSeq/merged_reads/mouse_${condition}/*.fastq)
destination_dir="/home/fabianje/repos/documenting_NIH/fabian/data/isoseq/fastq/$condition"
mkdir -p "$destination_dir"
cp "${fastq_files[@]}" "$destination_dir"
cat "${fastq_files[@]}" > "$destination_dir/${condition}_concat.fastq"
