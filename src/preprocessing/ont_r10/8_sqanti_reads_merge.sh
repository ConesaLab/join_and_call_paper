#!/bin/bash
#SBATCH --job-name=sqanti_reads_merge_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_reads_merge_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_reads_merge_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --qos=short
#SBATCH --time=24:00:00

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

sqanti_dir="/home/cmonzo/software/SQANTI3-5.2"

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
ref_dir="${base_dir}/ref"
sqanti_out_dir="${base_dir}/analysis/sqanti_reads"

ref_annotation="${ref_dir}/gencode.v49_SIRV.gtf"
assembly="${ref_dir}/GRCh38_SIRV.fa"

export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/

python3 ${sqanti_dir}/sqanti_reads.py \
    --genome ${assembly} \
    --annotation ${ref_annotation} \
    -de design_ont_r10.csv \
    -i ${sqanti_out_dir} \
    -f "condition" \
    -p "ONT_R10" \
    -d ${sqanti_out_dir}/out \
    --force_id_ignore \
    -t 12 -n 1 --skip_hash
