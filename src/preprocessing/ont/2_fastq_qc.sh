#!/bin/bash
#SBATCH --job-name=fastqc_ont
#SBATCH --output=../analysis/logs/fastqc_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-24
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env
module load fastqc

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

dir="./storage/gge/nih/Nanopore/analysis/fastq_qc/"

file_dir=$(dirname $file)
ori_name=$(basename $file)

filename="${ori_name%.fastq}"

echo $filename

# Run fastqc

fastqc -t 8 -o ../analysis/fastqc/ -f fastq $file

