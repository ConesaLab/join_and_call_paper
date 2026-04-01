#!/bin/bash
#SBATCH --job-name=sqanti_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
ref_dir="${base_dir}/ref"
gff_dir="${base_dir}/gff"
qc_dir="${base_dir}/analysis"

ref_annotation="${ref_dir}/gencode.v49_SIRV.gtf"
assembly="${ref_dir}/GRCh38_SIRV.fa"

isoforms_gff="${gff_dir}/${filename}_primary_aln.gff"
isoforms_gtf="${gff_dir}/${filename}_primary_aln.gtf"

# Convert GFF to GTF -- SQANTI3_dev's gtfToGenePred chokes on GFF input
gffread ${isoforms_gff} -T -o ${isoforms_gtf}

# Clean previous output to avoid stale cached files
rm -rf ${qc_dir}/run_SQANTI/${filename}
mkdir -p ${qc_dir}/run_SQANTI/${filename}

sqanti_dir="$HOME/tools/SQANTI3_dev"

python3 ${sqanti_dir}/sqanti3_qc.py \
    --min_ref_len 0 --skipORF \
    --dir "${qc_dir}/run_SQANTI/${filename}" \
    --output "${filename}" \
    ${isoforms_gtf} ${ref_annotation} ${assembly}
