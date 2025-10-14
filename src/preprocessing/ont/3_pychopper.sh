#!/bin/bash
#SBATCH --job-name=pychopper_ont
#SBATCH --output=pychopper_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --array=0-168
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate pychopper

# Create array of files
readarray myarray < list_fastqs.txt

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

#dir="./storage/gge/nih/Nanopore/analysis/fastqs/merged_fastqs/"

file_dir=$(dirname $file)
#bam_dir="/storage/gge/nih/Nanopore/bams/"
#gff_dir="/storage/gge/nih/Nanopore/gffs/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

#ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"

#assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

pychopper -k PCS111 -r ${file_dir}/${filename}_pychopper_report.pdf -u ${file_dir}/${filename}_unclassified.fastq -w ${file_dir}/${filename}_rescued.fastq -t 4 ${file_dir}/${filename}.fastq ${file_dir}/${filename}_full_length.fastq
