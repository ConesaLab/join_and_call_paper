#!/bin/bash
#SBATCH --job-name=map_nih
#SBATCH --output=map_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=25gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-95 # 191

source activate test_is3
module load samtools

# genome 
# Define reference annotation, reference assembly and FASTQ for isoforms
ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

# Create an array of files
readarray myarray < flnc.fofn

# Read the file corresponding to the task
bam=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $bam


# Get file name wo estension
reads_dir=$(dirname $bam)
reads_file=$(basename $bam)

filename="${reads_file%.*.*}"

echo $filename

# Convert to bam to fastq
bamtools convert -format fastq -in $bam > ${reads_dir}/${filename}.flnc.fastq

# map using minimap2
minimap2 -ax splice:hq -uf --MD -t 8 ${assembly} ${reads_dir}/${filename}.flnc.fastq > ${reads_dir}/${filename}.sam

# SAM to GTF
samtools view -bS -F0x900 ${reads_dir}/${filename}.sam | samtools sort -o ${reads_dir}/${filename}_primary_aln_sorted.bam