#!/bin/bash
#SBATCH --job-name=nf_mando
#SBATCH --qos long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH -t 14-00:00:00     # Walltime
#SBATCH -o log/nf_mando_%A.out

# use this to run a longer nextflow workflow
# example: sbatch rerun_mando_jc.sbatch \
#           --metadata_ind /home/fabianje/repos/documenting_NIH/fabian/data/output/isoseq_trimmed/mandalorion/run1_data/mandalorion/mandalorion_out_metadata_ind.tsv
#           --metadata_concat /home/fabianje/repos/documenting_NIH/fabian/data/output/isoseq_trimmed/mandalorion/run1_data/mandalorion/mandalorion_out_metadata_concat.tsv
#           --metadata_tama_condition /home/fabianje/repos/documenting_NIH/fabian/data/output/isoseq_trimmed/mandalorion/run1_data/tama_condition/metadata_tama_condition.tsv
#           --result_name isoseq_trimmed/mandalorion/run1


# params.metadata_ind = ""
# params.metadata_concat = ""
# params.metadata_tama_condition = ""
module load anaconda
conda activate nextflow

echo "nextflow rerun_mando_jc.nf $@"

nextflow rerun_mando_jc.nf "$@"
