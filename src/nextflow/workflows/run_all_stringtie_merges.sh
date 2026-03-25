#!/bin/bash
#SBATCH --job-name=all_stmerge
#SBATCH --qos long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 15-00:00:00
#SBATCH -o log/all_stmerge_%A.out

# Submits the StringTie merge + SQANTI3 workflow for all 10 data/algorithm
# combinations sequentially (each waits for the previous to finish).

PB_BASE="/home/fabianje/repos/documenting_NIH/fabian/data/output/isoseq"
ONT_BASE="/home/fabianje/repos/documenting_NIH/fabian/data/output/ont"
# SLURM copies batch scripts to its spool dir, so realpath "$0" resolves
# incorrectly. Use SLURM_SUBMIT_DIR (where sbatch was invoked) instead.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$(realpath "$0")")}"
CORRECTED_META="$(realpath "${SUBMIT_DIR}/../scripts/stringtie_merge/metadata")"

# --- PacBio ---

echo "=== PacBio FLAIR ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${CORRECTED_META}/pb_flair_metadata_ind.tsv" \
    --result_name isoseq/flair_ar_sr/stringtie_merge

echo "=== PacBio IsoQuant ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${PB_BASE}/isoquant/run3_data/isoquant/isoquant_out_metadata_ind.tsv" \
    --result_name isoseq/isoquant/stringtie_merge

echo "=== PacBio Bambu ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${PB_BASE}/bambu/run3_data/bambu/bambu_out_metadata_ind.tsv" \
    --result_name isoseq/bambu/stringtie_merge

echo "=== PacBio TALON ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${PB_BASE}/talon/run3_data/talon/talon_out_metadata_ind.tsv" \
    --result_name isoseq/talon/stringtie_merge

echo "=== PacBio IsoSeq (Alejandro) ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${CORRECTED_META}/pb_isoseq_metadata_ind.tsv" \
    --result_name isoseq/isoseq_pipeline/stringtie_merge

echo "=== PacBio Mandalorion (Alejandro) ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${CORRECTED_META}/pb_mandalorion_metadata_ind.tsv" \
    --result_name isoseq/mandalorion/stringtie_merge

# --- ONT ---

echo "=== ONT FLAIR ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${CORRECTED_META}/ont_flair_metadata_ind.tsv" \
    --result_name ont/flair_ar_sr/stringtie_merge

echo "=== ONT IsoQuant ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${ONT_BASE}/isoquant/run2_data/isoquant/isoquant_out_metadata_ind.tsv" \
    --result_name ont/isoquant/stringtie_merge

echo "=== ONT Bambu ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${ONT_BASE}/bambu/run2_data/bambu/bambu_out_metadata_ind.tsv" \
    --result_name ont/bambu/stringtie_merge

echo "=== ONT TALON ==="
sbatch --wait nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "${ONT_BASE}/talon/run1_data/talon/talon_out_metadata_ind.tsv" \
    --result_name ont/talon/stringtie_merge

echo "=== All StringTie merge runs submitted ==="
