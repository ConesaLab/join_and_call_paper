#!/bin/bash
#SBATCH --job-name=map_ont
#SBATCH --output=../analysis/logs/map_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --qos=short
#SBATCH --time=12:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
bam_dir="/storage/gge/nih/Nanopore/data/bams/"
gff_dir="/storage/gge/nih/Nanopore/data/gffs/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"

assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

#minimap2 -ax splice -uf --MD -t 4 /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa <(cat ../fastqs/merged_fastqs/B15_2_full_length.fastq ../fastqs/merged_fastqs/B15_2_rescued.fastq ../fastqs/merged_fastqs/B15_2_unclassified.fastq) > ../bams/test_spliceUF_B15_2.sam

# run minimap
minimap2 -ax splice -uf --MD -t 4 ${assembly} ${file_dir}/${filename}.fastq > ${bam_dir}/${filename}_splice.sam

# SAM to GTF
samtools view -bS -F0x900 ${bam_dir}/${filename}_splice.sam | samtools sort -o ${bam_dir}/${filename}_primary_aln_sorted.bam

spliced_bam2gff -t 1000000 -M ${bam_dir}/${filename}_primary_aln_sorted.bam > ${gff_dir}/${filename}_primary_aln.gff
