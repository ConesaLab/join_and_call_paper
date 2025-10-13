#!/bin/bash
#SBATCH --job-name=merge_counts
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=10gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/merge_counts_%A_%a.out

# Input
array_config="$1"
SCRIPT_DIR="$2"
n_cores="$SLURM_CPUS_PER_TASK"
task_number="$SLURM_ARRAY_TASK_ID"

# read corresponding line of array config
IFS=$'\t' read -r quant_fofn merge out_mat <<< "$(sed -n "${task_number}p" $array_config)"

module load R/4.2.1
module load anaconda
#source activate SQANTI3.env

# Run count table merge
Rscript "${SCRIPT_DIR}/merge_quantification.R" \
        "$quant_fofn" \
        "$merge" \
        "$out_mat"

