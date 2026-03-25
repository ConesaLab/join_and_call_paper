#!/bin/bash
#SBATCH --job-name=mapQC_ont_r10
#SBATCH --output=../analysis/logs/mapQC_%A_%a.out
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

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

ori_name=$(basename $file)
filename="${ori_name%.fastq}"

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
bam_dir="${base_dir}/bam"
qc_dir="${base_dir}/analysis/flagstat"

mkdir -p ${qc_dir}

echo $filename

samtools flagstat ${bam_dir}/${filename}_primary_aln_sorted.bam > ${qc_dir}/${filename}_primary_aln_flagstat.txt
