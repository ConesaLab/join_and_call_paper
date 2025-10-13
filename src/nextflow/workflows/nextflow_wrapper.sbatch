#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --qos long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH -t 14-00:00:00     # Walltime
#SBATCH -o log/nextflow_%A.out

# use this to run a longer nextflow workflow
# example: sbatch nextflow_wrapper.sbatch --data ont_subset --algorithm flair --result_name ont_subset/flair/test1

module load anaconda
conda activate nextflow

echo "nextflow main_workflow.nf $@"

nextflow main_workflow.nf "$@"
