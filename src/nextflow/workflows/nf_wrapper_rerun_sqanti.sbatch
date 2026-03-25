#!/bin/bash
#SBATCH --job-name=nf_sqanti
#SBATCH --qos long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH -t 14-00:00:00     # Walltime
#SBATCH -o log/nf_sqanti_%A.out

# use this to run a longer nextflow workflow
# example: sbatch rerun_sqanti_wrapper.sbatch --metadata_merged /storage/gge/Alejandro/for_fabian/paper_stuff/isoseq_bambu_run1_data/merge_counts/metadata_merged.tsv --result_name alejandro/isoseq/bambu/run1

module load anaconda
conda activate nextflow

echo "nextflow rerun_sqanti.nf $@"

nextflow rerun_sqanti.nf "$@"
