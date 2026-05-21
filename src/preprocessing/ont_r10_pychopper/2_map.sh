#!/bin/bash
#SBATCH --job-name=map_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/map_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/map_%A_%a.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=fjetzinger@biobam.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

source config.sh

mkdir -p "${bam_dir}" "${gff_dir}" "${logs_dir}"

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs_for_map.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo "Input: ${file}"

ori_name=$(basename "$file")
filename="${ori_name%_for_map.fastq}"

echo "Mapping ${filename}"

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
