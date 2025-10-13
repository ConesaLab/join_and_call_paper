#!/bin/bash
#SBATCH --job-name=workflow_isoquant
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=10gb # Per processor memory
#SBATCH -t 7-00:00:00     # Walltime
#SBATCH -o log/workflow_isoquant_%A.out

# paths
work_dir="/home/fabianje/repos/documenting_NIH/fabian/data/output/ont/isoquant_run_full_1"
metadata_samples="/home/fabianje/repos/documenting_NIH/fabian/data/metadata/ont_samples.tsv"
metadata_concat="/home/fabianje/repos/documenting_NIH/fabian/data/metadata/ont_concat_samples.tsv"
empty_report="/home/fabianje/repos/documenting_NIH/fabian/reports/empty_report"
report_dir="/home/fabianje/repos/documenting_NIH/fabian/reports/$(basename $work_dir)"
src_dir="/home/fabianje/repos/documenting_NIH/fabian/src/nextflow/scripts"

# files
genome="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
annotation_db="/home/fabianje/repos/documenting_NIH/fabian/data/metadata/isoquant/mm39.ncbiRefSeq_SIRV.db"

# tools
tama_location="/home/fabianje/tools/tama"
# use dev version of SQANTI3 cloned from github to get the intergenic/genic intron bugfix
sqanti_location="/home/fabianje/tools/SQANTI3_dev"

# workflow

echo "Starting workflow...."

# isoquant
isoquant_wd="${work_dir}/isoquant"
isoquant_metadata_ind="${isoquant_wd}/isoquant_out_metadata_ind.tsv"
isoquant_metadata_concat="${isoquant_wd}/isoquant_out_metadata_concat.tsv"

echo "Running isoquant..."

bash "${src_dir}/algorithms/isoquant/run_isoquant.sh" \
  --wd "${isoquant_wd}" \
  --genome "${genome}" \
  --annotation "${annotation_db}" \
  --metadata_samples "${metadata_samples}" \
  --metadata_concat "${metadata_concat}" \
  --data_type "nanopore"

echo "finished run_isoquant.sh"  

# tama_condition
tama_condition_wd="${work_dir}/tama_condition"
metadata_tama_condition="${tama_condition_wd}/metadata_tama_condition.tsv"

echo "Running TAMA for conditions..."

bash "${src_dir}/tama/run_tama_condition.sh" \
    --wd "$tama_condition_wd" \
    --metadata_samples "${isoquant_metadata_ind}" \
    --tama_path "${tama_location}"

# tama_full
tama_full_wd="${work_dir}/tama_full"
metadata_tama_full="${tama_full_wd}/metadata_tama_full.tsv"

echo "Running TAMA for ind and concat..."

bash "${src_dir}/tama/run_tama_full.sh" \
    --wd "$tama_full_wd" \
    --metadata_samples "${isoquant_metadata_ind}" \
    --metadata_concat "${isoquant_metadata_concat}" \
    --tama_path "${tama_location}"

# sqanti3
sqanti3_wd="${work_dir}/sq3"
metadata_merged="${sqanti3_wd}/metadata_merged.tsv"

echo "Running SQANTI3..."

bash "${src_dir}/sqanti3/run_sqanti3_qc.sh" \
    --wd "$sqanti3_wd" \
    --genome "${genome}" \
    --annotation "${annotation}" \
    --metadata_ind "${isoquant_metadata_ind}" \
    --metadata_concat "${isoquant_metadata_concat}" \
    --metadata_tama "${metadata_tama_condition}" \
    --metadata_tama_full "${metadata_tama_full}" \
    --sqanti_path "${sqanti_location}" \
    --force_id_ignore true

# merge counts
merge_counts_wd="${work_dir}/merge_counts"
merge_counts_array_config="${merge_counts_wd}/merge_counts_array_config.tsv"

echo "Running merge counts..."

bash "${src_dir}/merge_counts/merge_counts.sh" \
    --wd "${merge_counts_wd}" \
    --metadata_merged "${metadata_merged}"
    
echo "Preparing report directory..."

bash "${src_dir}/prepare_report/prepare_report.sh" \
    --wd "${report_dir}" \
    --empty_report "${empty_report}" \
    --metadata_merged "${metadata_merged}" \
    --merge_counts_array_config "${merge_counts_array_config}"
    
echo "Finished!"
