#!/bin/bash
#SBATCH --job-name=mapQC_ont
#SBATCH --output=../analysis/logs/mapQC_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
bam_dir="/storage/gge/nih/Nanopore/bams/"
qc_dir="/storage/gge/nih/Nanopore/analysis/flagstat/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

samtools flagstat ${bam_dir}/${filename}_primary_aln_sorted.bam > ${qc_dir}/${filename}_primary_aln_flagstat.txt
