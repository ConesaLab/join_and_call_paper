#!/bin/bash
#SBATCH --job-name=flair_junc
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH --array 1-12     # array
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/flair_junc_%A_%a.out

module load anaconda
source activate flair_basic_conda_env

flair_location="$1"
config_file="$2"
out_location="$3"
joblog="$4"

export PATH="${flair_location}/bin:$PATH"

line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1))

read -r cond sample bam <<< "$(head -n $line_number_metadata $config_file | tail -n 1)"

echo -e "${SLURM_JOB_ID}\t${SLURM_ARRAY_TASK_ID}\t${sample}" >> $joblog

echo junctions_from_sam -s "${bam}" -n "${out_location}/${sample}"
junctions_from_sam -s "${bam}" -n "${out_location}/${sample}"
