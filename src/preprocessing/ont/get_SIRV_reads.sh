#!/bin/bash
#SBATCH --job-name=SIRV_counts_nanopore
#SBATCH --output=../analysis/logs/SIRV_counts_nanopore_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15gb
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

dir="/storage/gge/nih/Nanopore/analysis/SIRV_counts/"

bam_dir="/storage/gge/nih/Nanopore/data/bams/"
gff_dir="/storage/gge/nih/Nanopore/data/gffs/"
gtf_dir="/storage/gge/nih/Nanopore/data/gtfs/"
sqanti_dir="/storage/gge/nih/Nanopore/analysis/run_SQANTI/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename


cat ${sqanti_dir}/${filename}_SIRV/${filename}_SIRV_classification.txt | cut -f 8 | sort | uniq -c > ${dir}/${filename}_SIRV_counts.txt
