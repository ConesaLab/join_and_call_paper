#!/bin/bash
#SBATCH --job-name=fastqc_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/fastqc_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/fastqc_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env
module load fastqc

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

analysis_dir="/storage/gge/Fabian/ont_r10_sy5y/analysis"
mkdir -p ${analysis_dir}/fastqc

fastqc -t 8 -o ${analysis_dir}/fastqc/ -f fastq $file
