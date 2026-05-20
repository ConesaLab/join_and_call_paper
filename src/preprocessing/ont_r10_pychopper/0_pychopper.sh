#!/bin/bash
#SBATCH --job-name=pychopper_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/pychopper_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/pychopper_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3
#
# Pychopper orientation/trimming for ONT R10 SY5Y (alternative to Dorado-only + minimap2 -ub).
# Design reference: src/preprocessing/ont/3_pychopper.sh (mouse PCS111; here LSK114 for SQK-LSK114).
# Input: Dorado-trimmed FASTQs in fastq/ (same as default ont_r10 branch).
# Kit: -k LSK114 (not the full product string "ONT SQK-LSK114").

source ~/.bashrc

module load samtools
conda deactivate
conda activate pychopper

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${_script_dir}/config.sh"

mkdir -p "${fastq_out}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < "${_script_dir}/list_fastqs.fof"

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "Input: ${file}"

ori_name=$(basename "$file")
filename="${ori_name%.fastq}"

echo "Sample: ${filename} | Kit: ${PYCHOPPER_KIT}"

pychopper -k "${PYCHOPPER_KIT}" \
  -r "${fastq_out}/${filename}_pychopper_report.pdf" \
  -u "${fastq_out}/${filename}_unclassified.fastq" \
  -w "${fastq_out}/${filename}_rescued.fastq" \
  -t "${SLURM_CPUS_PER_TASK:-8}" \
  "${fastq_in}/${filename}.fastq" \
  "${fastq_out}/${filename}_full_length.fastq"

echo "Done: ${filename}"
