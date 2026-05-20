#!/bin/bash
#SBATCH --job-name=sqanti_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

source config.sh

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs_for_map.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

ori_name=$(basename "$file")
filename="${ori_name%_for_map.fastq}"

isoforms_gff="${gff_dir}/${filename}_primary_aln.gff"
isoforms_gtf="${gff_dir}/${filename}_primary_aln.gtf"

gffread "${isoforms_gff}" -T -o "${isoforms_gtf}"

rm -rf "${sqanti_dir}/${filename}"
mkdir -p "${sqanti_dir}/${filename}"

sqanti_tool="$HOME/tools/SQANTI3_dev"

python3 "${sqanti_tool}/sqanti3_qc.py" \
  --min_ref_len 0 --skipORF \
  --dir "${sqanti_dir}/${filename}" \
  --output "${filename}" \
  "${isoforms_gtf}" "${ref_annotation}" "${assembly}"
