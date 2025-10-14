#!/bin/bash
#SBATCH --job-name=lima_nih
#SBATCH --output=lima_bam_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-7

source activate test_is3

pwd;date

# primers fasta
primers=/storage/gge/nih/PacBio_IsoSeq/Iso-Seq_12_Barcoded_cDNA_Primers.fasta

# Create an array of files
readarray -t myarray < hifi.fofn

# Read the file corresponding to the task
bam=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $bam

# create a folder for the demultiplexed reads
reads_dir=$(dirname $bam)

mkdir $reads_dir/demultiplexed_bam


# Get file name wo estension

reads_file=$(basename $bam)

filename="${reads_file%%.*}"

lima $bam $primers $reads_dir/demultiplexed_bam/$filename.fl.bam --isoseq --peek-guess --dump-clips -j 8
