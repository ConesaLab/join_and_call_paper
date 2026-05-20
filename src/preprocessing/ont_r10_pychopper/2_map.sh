#!/bin/bash
#SBATCH --job-name=map_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/map_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/map_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3
#
# Map Pychopper-oriented reads with minimap2 -uf (forward strand; cf. ont_r10/2_map.sh -ub).
# Design reference: src/preprocessing/ont/4_map.sh (-uf after Pychopper).

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${_script_dir}/config.sh"

mkdir -p "${bam_dir}" "${gff_dir}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < "${_script_dir}/list_fastqs_for_map.fof"

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "Input: ${file}"

ori_name=$(basename "$file")
# SRR31732191_for_map.fastq -> SRR31732191
filename="${ori_name%_for_map.fastq}"

echo "Mapping ${filename}"

# Pychopper LSK114 orients reads; -uf forces forward strand (vs ont_r10 -ub for unoriented Dorado reads).
# Keep R10 junction sensitivity flags from ont_r10/2_map.sh.
minimap2 -ax splice -uf --secondary=no --MD -t 8 -k 14 -w 4 \
  "${assembly}" "${file}" > "${bam_dir}/${filename}_splice.sam"

samtools view -bS -F0x900 "${bam_dir}/${filename}_splice.sam" \
  | samtools sort -o "${bam_dir}/${filename}_primary_aln_sorted.bam"

samtools index "${bam_dir}/${filename}_primary_aln_sorted.bam"

spliced_bam2gff -t 1000000 -M "${bam_dir}/${filename}_primary_aln_sorted.bam" \
  > "${gff_dir}/${filename}_primary_aln.gff"

head -n1 "${gff_dir}/${filename}_primary_aln.gff" > "${gff_dir}/${filename}_SIRV.gff"
grep SIRV "${gff_dir}/${filename}_primary_aln.gff" >> "${gff_dir}/${filename}_SIRV.gff"

echo "Done: ${filename}"
