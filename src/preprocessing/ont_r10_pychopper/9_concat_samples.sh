#!/bin/bash
#SBATCH --job-name=concat_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/concat_samples_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/concat_samples_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=fjetzinger@biobam.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=24:00:00

source ~/.bashrc

module load samtools

source config.sh

mkdir -p "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="

echo "Concatenating Pychopper mapping FASTQs..."
cat "${fastq_out}/SRR31732191_for_map.fastq" \
  "${fastq_out}/SRR31732198_for_map.fastq" \
  "${fastq_out}/SRR31732209_for_map.fastq" \
  "${fastq_out}/SRR31732210_for_map.fastq" \
  > "${fastq_out}/SY5Y_concat_pychopper.fastq"

echo "Merging BAMs..."
samtools merge -o "${bam_dir}/SY5Y_concat_aln_sorted_pychopper.bam" \
  "${bam_dir}/SRR31732191_primary_aln_sorted.bam" \
  "${bam_dir}/SRR31732198_primary_aln_sorted.bam" \
  "${bam_dir}/SRR31732209_primary_aln_sorted.bam" \
  "${bam_dir}/SRR31732210_primary_aln_sorted.bam"

samtools index "${bam_dir}/SY5Y_concat_aln_sorted_pychopper.bam"

echo "Done. Outputs:"
echo "  ${fastq_out}/SY5Y_concat_pychopper.fastq"
echo "  ${bam_dir}/SY5Y_concat_aln_sorted_pychopper.bam"
