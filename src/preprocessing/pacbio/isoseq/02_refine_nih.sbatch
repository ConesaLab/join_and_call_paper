#!/bin/bash
#SBATCH --job-name=refine_nih
#SBATCH --output=refine_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --qos=short
#SBATCH --time=5:00:00
#SBATCH --array=0-95

source activate test_is3

# primers fasta
primers=/storage/gge/nih/PacBio_IsoSeq/Iso-Seq_12_Barcoded_cDNA_Primers.fasta

# Create an array of files
readarray myarray < fl.fofn

# Read the file corresponding to the task
bam=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $bam

# create a folder for the demultiplexed reads
reads_dir=$(dirname $bam)

mkdir $reads_dir/flnc


# Get file name wo xstension

reads_file=$(basename $bam)

filename="${reads_file%.*}"

isoseq3 refine $bam $primers $reads_dir/flnc/$filename.flnc.bam --require-polya -j 8
