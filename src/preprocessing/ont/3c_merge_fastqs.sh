#!/bin/bash
#SBATCH --job-name=merge_fastqs
#SBATCH --output=../analysis/logs/merge_fastqs_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

# Create array of files
readarray myarray < concatenations.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

eval $file
