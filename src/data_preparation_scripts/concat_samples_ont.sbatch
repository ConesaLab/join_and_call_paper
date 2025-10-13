#!/bin/bash
#SBATCH --job-name=concat_samples
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=40gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o concat_samples_%A_%a.out

fastq_dir=$1
bam_dir=$2
fastq_target_dir=$3
bam_target_dir=$4

cat $fastq_dir/B3*.fastq > $fastq_target_dir/B100K0_concat.fastq
cat $fastq_dir/K3*.fastq > $fastq_target_dir/B0K100_concat.fastq

module load samtools

samtools merge -o "${bam_target_dir}/B100K0_concat_aln_sorted.bam" "${bam_dir}"/B3*.bam
samtools merge -o "${bam_target_dir}/B0K100_concat_aln_sorted.bam" "${bam_dir}"/K3*.bam

