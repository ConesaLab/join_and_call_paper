#!/bin/bash
#SBATCH --job-name=bam2gff
#SBATCH --output=logs/bam2gff_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=0:05:00
#SBATCH --array=11-384



source activate test_is3

# Create an array of files
readarray myarray < mapped_reads.fofn

# Read the file corresponding to the task
bam=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $bam

# Get file name wo estension
reads_dir=$(dirname $bam)
reads_file=$(basename $bam)

filename="${reads_file%.*}"

echo $filename


spliced_bam2gff -t 1000000 -M ${bam} > ${reads_dir}/${filename}.gff