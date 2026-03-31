#!/bin/bash
#SBATCH --job-name=map_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/map_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/map_%A_%a.log
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

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
bam_dir="${base_dir}/bam"
gff_dir="${base_dir}/gff"

assembly="${base_dir}/ref/GRCh38_SIRV.fa"

mkdir -p ${bam_dir} ${gff_dir}

echo "Mapping ${filename}"

# Parameters from source paper; adapters already trimmed by Dorado, no pychopper
# -ub: both strands (reads not pre-oriented)
# --secondary=no: no secondary alignments
# -k 14 -w 4: sensitive junction detection
# --MD: mismatch info for downstream tools
minimap2 -ax splice -ub --secondary=no --MD -t 8 -k 14 -w 4 \
    ${assembly} ${file_dir}/${filename}.fastq > ${bam_dir}/${filename}_splice.sam

# Filter supplementary alignments, convert to sorted BAM
samtools view -bS -F0x900 ${bam_dir}/${filename}_splice.sam \
    | samtools sort -o ${bam_dir}/${filename}_primary_aln_sorted.bam

samtools index ${bam_dir}/${filename}_primary_aln_sorted.bam

spliced_bam2gff -t 1000000 -M ${bam_dir}/${filename}_primary_aln_sorted.bam \
    > ${gff_dir}/${filename}_primary_aln.gff

# Extract SIRV-only GFF subset for separate SQANTI analysis
head -n1 ${gff_dir}/${filename}_primary_aln.gff > ${gff_dir}/${filename}_SIRV.gff
grep SIRV ${gff_dir}/${filename}_primary_aln.gff >> ${gff_dir}/${filename}_SIRV.gff

echo "Done: ${filename}"
