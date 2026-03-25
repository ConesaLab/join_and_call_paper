#!/bin/bash
#SBATCH --job-name=nf_stmerge
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=5gb
#SBATCH -t 1-00:00:00
#SBATCH -o log/nf_stmerge_%A.out

module load anaconda
conda activate nextflow

echo "nextflow stringtie_merge_workflow.nf $@"

nextflow stringtie_merge_workflow.nf "$@"
