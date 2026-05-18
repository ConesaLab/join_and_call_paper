#!/bin/bash
#SBATCH --job-name=ont_r10_bam_readlen
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/bam_readlen_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/bam_readlen_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#
# Writes lengths_seq_only/ont/<sample>.bam.readlen.txt — one integer per line (SEQ length
# for each primary-mapped alignment, -F 2308). Same format as src/util/get_read_lengths.sbatch
# for mouse NIH .bam.readlen.txt files (no subsampling).
#
# Submit from repo: sbatch src/preprocessing/ont_r10/10_get_read_lengths_bam.sh
# Requires: list_fastqs.fof and *_primary_aln_sorted.bam under base_dir/bam/.

set -euo pipefail
LC_ALL=C

source ~/.bashrc
module load samtools

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="/storage/gge/Fabian/ont_r10_sy5y"
BAMDIR="${BASE_DIR}/bam"
OUT_DIR="${BASE_DIR}/analysis/read_qc/lengths_seq_only/ont"
FOF="${SCRIPT_DIR}/list_fastqs.fof"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "=== Job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="
echo "OUT_DIR=${OUT_DIR}"

mkdir -p "${OUT_DIR}"
mkdir -p "${BASE_DIR}/analysis/logs"

write_bam_seq_lengths() {
  local bam="$1"
  local out="$2"
  if [[ ! -e "$bam" ]]; then
    echo "WARN: BAM not found: $bam" >&2
    : >"$out"
    return
  fi
  samtools view -@ "$THREADS" -F 2308 -- "$bam" |
    awk '{s=$10; if (s == "*" || s == "") { print 0 } else { print length(s) } }' >"$out"
}

if [[ ! -e "$FOF" ]]; then
  echo "ERROR: missing ${FOF}" >&2
  exit 1
fi

mapfile -t FASTQ_PATHS < "$FOF"

for fq in "${FASTQ_PATHS[@]}"; do
  [[ -z "${fq// }" ]] && continue
  ori_name="$(basename "$fq")"
  sample="$ori_name"
  [[ "$sample" == *.gz ]] && sample="${sample%.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"
  bam="${BAMDIR}/${sample}_primary_aln_sorted.bam"
  out="${OUT_DIR}/${sample}.bam.readlen.txt"
  echo "[${sample}] ${bam} -> ${out}"
  write_bam_seq_lengths "$bam" "$out"
done

# Optional: merged BAM (uncomment if needed)
# write_bam_seq_lengths "${BAMDIR}/SY5Y_concat_aln_sorted.bam" "${OUT_DIR}/SY5Y_concat.bam.readlen.txt"

echo "Done. Outputs in: ${OUT_DIR}"
