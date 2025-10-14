#!/bin/bash
#SBATCH --job-name=UJC_nanopore
#SBATCH --output=../analysis/logs/UJC_nanopore_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15gb
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

# Script adapted from Luis Ferrandez's work

module load bedtools

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

dir="/storage/gge/nih/Nanopore/analysis/UJC"

bam_dir="/storage/gge/nih/Nanopore/data/bams/"
gff_dir="/storage/gge/nih/Nanopore/data/gffs/"
gtf_dir="/storage/gge/nih/Nanopore/data/gtfs/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

# Clean GTF. First, modify the primary_aln.gff by adding gene_name=NA to fit GTFTOOLS (works OK now)
awk -F'\t' -v OFS="\t" '{$9=$9" gene_id \"NA\";"};1' ${gff_dir}/${filename}_primary_aln.gff > ${gtf_dir}/${filename}_primary_aln_mod.gtf

# GTFtools to get introns per transcript
gtftools -i ${dir}/${filename}_introns.bed ${gtf_dir}/${filename}_primary_aln_mod.gtf

# BEDtools groupby to obtain Unique Junction Chains (UJC). Note intronic start coord is +1 since input is a BED file
awk -F'\t' -v OFS="\t" '{print $5,"chr"$1,$4,$2+1"_"$3}' ${dir}/${filename}_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{print $1,$2"_"$3"_"$4}' > ${dir}/${filename}_Nanopore_UJC.txt

# After finishing the job, run
# cat *UJC.txt > MAS-seq_all_UJC_B31.txt
# In prev. dir, get flnc ids that have 0 as min_cov (any SJ is not supported by Illumina)
# find -name "*b31*classification.txt" -exec cat {} \; | awk -F'\t' '$19=="0" {print $1}' > UJC_B31_all_SJ_covered/b31_flnc_ids_mas-seq_min_cov_eq_0.txt
# Then filter out previous UJC if flnc had 0 min_cov
# grep -w -v -Ff b31_flnc_ids_mas-seq_min_cov_eq_0.txt MAS-seq_all_UJC_B31.txt > MAS-seq_UJC_min_cov_gt_0_B31.txt
