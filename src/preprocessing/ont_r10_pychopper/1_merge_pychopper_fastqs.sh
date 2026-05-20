#!/bin/bash
#SBATCH --job-name=merge_pychopper_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/merge_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/merge_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-3

source ~/.bashrc

source config.sh

mkdir -p "${fastq_out}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

ori_name=$(basename "$file")
filename="${ori_name%.fastq}"

fl="${fastq_out}/${filename}_full_length.fastq"
rescued="${fastq_out}/${filename}_rescued.fastq"
out="${fastq_out}/${filename}_for_map.fastq"

echo "Merging ${fl} + ${rescued} -> ${out}"

cat "${fl}" "${rescued}" > "${out}"

echo "Done: ${filename}"
