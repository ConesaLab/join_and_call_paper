#!/bin/bash
#SBATCH --job-name=ont_r10_pychopper_bam_readlen
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/bam_readlen_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/bam_readlen_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#
# Submit: cd .../ont_r10_pychopper && sbatch 11_get_read_lengths_bam.sh

LC_ALL=C

source ~/.bashrc
module load samtools

set -euo pipefail

BASE_DIR="/storage/gge/Fabian/ont_r10_sy5y"
BAMDIR="${BASE_DIR}/bam_pychopper"
OUT_DIR="${BASE_DIR}/analysis/read_qc_pychopper/lengths_seq_only/ont"

FOF=""
if [[ -n "${LIST_FASTQS_FOF:-}" && -f "${LIST_FASTQS_FOF}" ]]; then
  FOF="${LIST_FASTQS_FOF}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  for c in "${SLURM_SUBMIT_DIR}/list_fastqs_for_map.fof"; do
    if [[ -f "$c" ]]; then
      FOF="$c"
      break
    fi
  done
fi
if [[ -z "${FOF}" && -f list_fastqs_for_map.fof ]]; then
  FOF="list_fastqs_for_map.fof"
fi

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "=== Job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="
echo "OUT_DIR=${OUT_DIR}"

mkdir -p "${OUT_DIR}"
mkdir -p "${BASE_DIR}/analysis/logs_pychopper"

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
  echo "ERROR: could not find list_fastqs_for_map.fof (run sbatch from ont_r10_pychopper/)." >&2
  exit 1
fi
echo "Using list_fastqs_for_map.fof: ${FOF}"

mapfile -t FASTQ_PATHS < "$FOF"

for fq in "${FASTQ_PATHS[@]}"; do
  [[ -z "${fq// }" ]] && continue
  ori_name="$(basename "$fq")"
  sample="${ori_name%_for_map.fastq}"
  bam="${BAMDIR}/${sample}_primary_aln_sorted.bam"
  out="${OUT_DIR}/${sample}.bam.readlen.txt"
  echo "[${sample}] ${bam} -> ${out}"
  write_bam_seq_lengths "$bam" "$out"
done

echo "Done. Outputs in: ${OUT_DIR}"
