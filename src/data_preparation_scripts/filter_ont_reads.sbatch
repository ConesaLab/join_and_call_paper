#!/bin/bash
#SBATCH --job-name=filter_fastq
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o filter_fastq_%A_%a.out

task_number="$SLURM_ARRAY_TASK_ID"

module load fastp
module load anaconda
conda deactivate
conda activate SQANTI3.env

dir="/home/fabianje/repos/documenting_NIH/fabian/data/output/ont/fastq"

# List all .fastq files, sort them, and pick the one corresponding to the task number
fastq_files=($(ls $dir/merged/*.fastq | sort))
selected_fastq=${fastq_files[$task_number-1]}

base_filename=$(basename $selected_fastq)

# Define the output filenames
filtered_fastq="$dir/length_filtered/${base_filename%.fastq}_300bp.fastq"

# Filter reads under 300bp using fastp
# disable quality filtering
# fastp -i "$selected_fastq" -o "$filtered_fastq" --length_required 300 -Q -A 


# create subset of 10000 reads
subset_file="$dir/length_filtered/subset/${base_filename%.fastq}_300bp_subset.fastq"
# head -n 40000 "$filtered_fastq" > $subset_file

# map
filename="${base_filename%.fastq}"
filename_subset="$(basename $subset_file .fastq)"
bam_dir="$(dirname "$dir")/bam/length_filtered"

assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

# run minimap
# minimap2 -ax splice -uf --MD -t 4 ${assembly} $filtered_fastq > ${bam_dir}/${filename}_splice.sam
# minimap2 -ax splice -uf --MD -t 4 ${assembly} $subset_file > ${bam_dir}/subset/${filename_subset}_splice.sam

# sam to bam
samtools view -bS -F0x900 ${bam_dir}/${filename}_splice.sam | samtools sort -o ${bam_dir}/${filename}_primary_aln_sorted.bam
samtools view -bS -F0x900 ${bam_dir}/subset/${filename_subset}_splice.sam | samtools sort -o ${bam_dir}/subset/${filename_subset}_primary_aln_sorted.bam
