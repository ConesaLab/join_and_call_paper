#!/bin/bash
#SBATCH --job-name=fastqc_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/fastqc_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/fastqc_%A_%a.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=fjetzinger@biobam.com
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

source config.sh

mkdir -p "${fastqc_dir}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs_for_map.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "FASTQ: ${file}"

fastqc -t 8 -o "${fastqc_dir}/" -f fastq "${file}"
