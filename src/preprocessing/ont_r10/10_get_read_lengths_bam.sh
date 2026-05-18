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
# Submit (Slurm copies this script to spool; see 9_count_reads_joint.sh header).
#   cd /path/to/join_and_call_paper && sbatch src/preprocessing/ont_r10/10_get_read_lengths_bam.sh
# or set LIST_FASTQS_FOF, or place list_fastqs.fof under BASE_DIR.
#
# Requires: list_fastqs.fof and *_primary_aln_sorted.bam under BASE_DIR/bam/.
#
# Strict mode after sourcing bashrc: with nounset enabled first, /etc/bashrc can fail on
# unset BASHRCSOURCED when sourcing ~/.bashrc.

LC_ALL=C

source ~/.bashrc
module load samtools

set -euo pipefail

BASE_DIR="/storage/gge/Fabian/ont_r10_sy5y"
BAMDIR="${BASE_DIR}/bam"
OUT_DIR="${BASE_DIR}/analysis/read_qc/lengths_seq_only/ont"

FOF=""
if [[ -n "${LIST_FASTQS_FOF:-}" && -f "${LIST_FASTQS_FOF}" ]]; then
  FOF="${LIST_FASTQS_FOF}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  for c in "${SLURM_SUBMIT_DIR}/src/preprocessing/ont_r10/list_fastqs.fof" \
    "${SLURM_SUBMIT_DIR}/list_fastqs.fof"; do
    if [[ -f "$c" ]]; then
      FOF="$c"
      break
    fi
  done
fi
if [[ -z "${FOF}" && -f "${BASE_DIR}/list_fastqs.fof" ]]; then
  FOF="${BASE_DIR}/list_fastqs.fof"
fi
if [[ -z "${FOF}" ]]; then
  _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  if [[ -f "${_script_dir}/list_fastqs.fof" ]]; then
    FOF="${_script_dir}/list_fastqs.fof"
  fi
fi

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

if [[ -z "${FOF}" || ! -f "${FOF}" ]]; then
  echo "ERROR: could not find list_fastqs.fof (Slurm spool breaks paths next to this script)." >&2
  echo "  Try: cd <repo> && sbatch src/preprocessing/ont_r10/10_get_read_lengths_bam.sh" >&2
  echo "  Or:  export LIST_FASTQS_FOF=/path/to/list_fastqs.fof" >&2
  echo "  Or:  cp list_fastqs.fof ${BASE_DIR}/" >&2
  exit 1
fi
echo "Using list_fastqs.fof: ${FOF}"

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
