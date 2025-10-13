#!/bin/bash
#SBATCH --job-name=split_bam
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=10gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/split_bam_%A.out

module load samtools

in_bam=$1
target_chr=$2
output=$3

samtools view $in_bam $target_chr -b > $output

samtools index $output
