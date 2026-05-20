#!/bin/bash
#SBATCH --job-name=sqanti_reads_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_reads_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_reads_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

source config.sh

sqanti_tool="$HOME/tools/SQANTI3_dev"

mkdir -p "${sqanti_reads_dir}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs_for_map.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

ori_name=$(basename "$file")
filename="${ori_name%_for_map.fastq}"

bam="${bam_dir}/${filename}_primary_aln_sorted.bam"

gff_path="${sqanti_reads_dir}/${filename}.gff"
gtf_path="${sqanti_reads_dir}/${filename}.gtf"

spliced_bam2gff -t 1000000 -M "${bam}" > "${gff_path}"

gffread "${gff_path}" -T -o "${gtf_path}"

rm -f "${gff_path}"

python3 "${sqanti_tool}/sqanti3_qc.py" \
  --skipORF --dir "${sqanti_reads_dir}/${filename}" \
  --output "${filename}_reads" --min_ref_len "0" \
  --report "skip" --force_id_ignore \
  "${gtf_path}" "${ref_annotation}" "${assembly}"
