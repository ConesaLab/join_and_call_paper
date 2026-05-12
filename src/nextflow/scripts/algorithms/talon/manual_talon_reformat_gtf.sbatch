#!/bin/bash
#SBATCH --job-name=talon_reformat_gtf
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=32gb
#SBATCH -t 1-00:00:00     # Walltime (1 day = 24h)
#SBATCH -o /storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_reformat_gtf/log/%x_%j.out
#SBATCH -e /storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_reformat_gtf/log/%x_%j.err

set -euo pipefail

REPO="/home/fabianje/repos/join_and_call_paper"
REF_DIR="/storage/gge/Fabian/ont_r10_sy5y/ref"
WORKDIR="${REF_DIR}/talon_manual_reformat_gtf"
OUT_GTF="${WORKDIR}/annotation/gencode.v49_SIRV_fixed_talon.gtf"
DB_PREFIX="${WORKDIR}/GRCh38_SIRV_empty"
REFORMAT="${REPO}/src/nextflow/scripts/algorithms/talon/reformat_annotation.sh"

mkdir -p "${WORKDIR}/annotation" "${WORKDIR}/log"

module load anaconda
source activate talon

rm -f "${OUT_GTF}"

# Omit --delete_tmp so annotation/tmp intermediates remain for validation.
bash "${REFORMAT}" \
  -a "${REF_DIR}/gencode.v49_SIRV_fixed.gtf" \
  -o "${OUT_GTF}" \
  --skip_gene_id_gene_name

# rm -f "${DB_PREFIX}.db" "${DB_PREFIX}.db-wal" "${DB_PREFIX}.db-shm"

# talon_initialize_database \
#   --f "${OUT_GTF}" \
#   --g GRCh38_SIRV \
#   --a gencode.v49_SIRV_fixed \
#   --o "${DB_PREFIX}"

echo "Done. GTF: ${OUT_GTF}"
