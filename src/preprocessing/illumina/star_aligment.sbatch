#!/bin/bash
#SBATCH --job-name=map_short_reads
#SBATCH --output=map_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=14:00:00
#SBATCH --array=0-25

source activate mapper

pwd;date

#Path to index
index_path="/storage/gge/Alejandro/genome_index/mouse_nih"

# Pair ends sufix
r1="_1.fq.gz"
r2="_2.fq.gz"

# Create an array of files
readarray -t myarray < pairs.fofn

# Read the file corresponding to the task
fastq=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $fastq

aligments_dir="/storage/gge/nih/Illumina_short_reads/short_reads/mapped"

# create a folder for the pair end reads
reads_dir=$(dirname $fastq)
reads_file=$(basename $fastq)
readarray -d "_" -t sample_arr <<< $reads_file

mkdir $aligments_dir/${sample_arr[0]}

# Map
STAR --runThreadN 8 \
	--genomeDir $index_path \
	--readFilesIn ${fastq}${r1} ${fastq}${r2} \
	--outFileNamePrefix $aligments_dir/${sample_arr[0]}/$reads_file \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterType BySJout \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNoverLmax 0.04 \
	--outFilterMismatchNmax 999 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--sjdbScore 1 \
	--genomeLoad NoSharedMemory \
	--outSAMtype BAM SortedByCoordinate \
	--twopassMode Basic \
	--readFilesCommand zcat
