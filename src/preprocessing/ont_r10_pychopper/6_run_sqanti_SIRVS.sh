#!/bin/bash
#SBATCH --job-name=sqanti_SIRV_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_SIRV_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_SIRV_%A_%a.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=fjetzinger@biobam.com
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

source config.sh

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs_for_map.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

ori_name=$(basename "$file")
filename="${ori_name%_for_map.fastq}"

isoforms_gff="${gff_dir}/${filename}_SIRV.gff"
isoforms_gtf="${gff_dir}/${filename}_SIRV.gtf"

gffread "${isoforms_gff}" -T -o "${isoforms_gtf}"

rm -rf "${sqanti_dir}/${filename}_SIRV"
mkdir -p "${sqanti_dir}/${filename}_SIRV"

sqanti_tool="$HOME/tools/SQANTI3_dev"

python3 "${sqanti_tool}/sqanti3_qc.py" \
  --min_ref_len 0 --skipORF \
  --dir "${sqanti_dir}/${filename}_SIRV" \
  --output "${filename}_SIRV" \
  "${isoforms_gtf}" "${ref_annotation}" "${assembly}"
