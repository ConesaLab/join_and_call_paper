#!/bin/bash
#SBATCH --job-name=sqanti_SIRV_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_SIRV_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_SIRV_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --qos=short
#SBATCH --time=10:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

# SIRV GFF subsets are extracted during 3_map.sh:
#   head -n1 ${gff} > ${sirv_gff}; grep SIRV ${gff} >> ${sirv_gff}

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

isoforms_gff="${gff_dir}/${filename}_SIRV.gff"
isoforms_gtf="${gff_dir}/${filename}_SIRV.gtf"

gffread ${isoforms_gff} -T -o ${isoforms_gtf}

rm -rf ${qc_dir}/run_SQANTI/${filename}_SIRV
mkdir -p ${qc_dir}/run_SQANTI/${filename}_SIRV

sqanti_dir="$HOME/tools/SQANTI3_dev"

python3 ${sqanti_dir}/sqanti3_qc.py \
    --min_ref_len 0 --skipORF \
    --dir "${qc_dir}/run_SQANTI/${filename}_SIRV" \
    --output "${filename}_SIRV" \
    ${isoforms_gtf} ${ref_annotation} ${assembly}
