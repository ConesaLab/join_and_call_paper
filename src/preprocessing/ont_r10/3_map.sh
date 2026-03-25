#!/bin/bash
#SBATCH --job-name=map_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/map_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/map_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --qos=short
#SBATCH --time=12:00:00
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

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
pychopper_dir="${base_dir}/fastq/pychopper"
bam_dir="${base_dir}/bam"
gff_dir="${base_dir}/gff"

assembly="${base_dir}/ref/GRCh38_SIRV.fa"

mkdir -p ${bam_dir} ${gff_dir}

echo "Mapping ${filename}"

# Same minimap2 settings as R9 pipeline (src/preprocessing/ont/4_map.sh)
minimap2 -ax splice -uf --MD -t 4 ${assembly} ${pychopper_dir}/${filename}_fl_rescued.fastq \
    > ${bam_dir}/${filename}_splice.sam

samtools view -bS -F0x900 ${bam_dir}/${filename}_splice.sam \
    | samtools sort -o ${bam_dir}/${filename}_primary_aln_sorted.bam

samtools index ${bam_dir}/${filename}_primary_aln_sorted.bam

spliced_bam2gff -t 1000000 -M ${bam_dir}/${filename}_primary_aln_sorted.bam \
    > ${gff_dir}/${filename}_primary_aln.gff

# Extract SIRV-only GFF subset for separate SQANTI analysis
head -n1 ${gff_dir}/${filename}_primary_aln.gff > ${gff_dir}/${filename}_SIRV.gff
grep SIRV ${gff_dir}/${filename}_primary_aln.gff >> ${gff_dir}/${filename}_SIRV.gff

echo "Done: ${filename}"
