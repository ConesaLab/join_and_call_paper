#!/bin/bash
#SBATCH --job-name=pychopper_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/pychopper_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/pychopper_%A_%a.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=fjetzinger@biobam.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --qos=long-mem
#SBATCH --time=7-00:00:00
#SBATCH --array=0-3
#
# Submit from this directory: cd .../ont_r10_pychopper && sbatch 0_pychopper.sh
# Kit: -k LSK114. Uses -B batching to limit peak RAM on full SRR FASTQs.

source ~/.bashrc

module load samtools
conda deactivate
conda activate pychopper

source config.sh

mkdir -p "${fastq_out}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="
echo "mem=${SLURM_MEM_PER_NODE:-?} batch=${PYCHOPPER_BATCH_SIZE}"

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "Input: ${file}"

ori_name=$(basename "$file")
filename="${ori_name%.fastq}"

echo "Sample: ${filename} | Kit: ${PYCHOPPER_KIT}"

pychopper -k "${PYCHOPPER_KIT}" \
  -B "${PYCHOPPER_BATCH_SIZE}" \
  -r "${fastq_out}/${filename}_pychopper_report.pdf" \
  -u "${fastq_out}/${filename}_unclassified.fastq" \
  -w "${fastq_out}/${filename}_rescued.fastq" \
  -t 4 \
  "${fastq_in}/${filename}.fastq" \
  "${fastq_out}/${filename}_full_length.fastq"

echo "Done: ${filename}"
