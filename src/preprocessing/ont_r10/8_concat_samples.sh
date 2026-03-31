#!/bin/bash
#SBATCH --job-name=concat_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/concat_samples_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/concat_samples_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=24:00:00

source ~/.bashrc

module load samtools

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
fastq_dir="${base_dir}/fastq"
bam_dir="${base_dir}/bam"

# Concatenate all 4 replicates into single FASTQ for Join&Call
echo "Concatenating FASTQs..."
cat ${fastq_dir}/SRR31732191.fastq \
    ${fastq_dir}/SRR31732198.fastq \
    ${fastq_dir}/SRR31732209.fastq \
    ${fastq_dir}/SRR31732210.fastq \
    > ${fastq_dir}/SY5Y_concat.fastq

echo "Merging BAMs..."
samtools merge -o ${bam_dir}/SY5Y_concat_aln_sorted.bam \
    ${bam_dir}/SRR31732191_primary_aln_sorted.bam \
    ${bam_dir}/SRR31732198_primary_aln_sorted.bam \
    ${bam_dir}/SRR31732209_primary_aln_sorted.bam \
    ${bam_dir}/SRR31732210_primary_aln_sorted.bam

samtools index ${bam_dir}/SY5Y_concat_aln_sorted.bam

echo "Done. Outputs:"
echo "  ${fastq_dir}/SY5Y_concat.fastq"
echo "  ${bam_dir}/SY5Y_concat_aln_sorted.bam"
