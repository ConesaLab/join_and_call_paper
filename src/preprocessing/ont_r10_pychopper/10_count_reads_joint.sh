#!/bin/bash
#SBATCH --job-name=ont_r10_pychopper_read_counts
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/count_reads_joint_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/count_reads_joint_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#
# Same TSV schema as ont_r10/9_count_reads_joint.sh; writes to read_qc_pychopper/.
#
# Submit:
#   cd <repo>/src/preprocessing/ont_r10_pychopper && sbatch 10_count_reads_joint.sh
# or: export LIST_FASTQS_FOF=/path/to/list_fastqs_for_map.fof

LC_ALL=C

source ~/.bashrc
module load samtools

set -euo pipefail

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${_script_dir}/config.sh"

OUT_TSV="${read_qc_dir}/read_numbers_joint.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

FOF=""
if [[ -n "${LIST_FASTQS_FOF:-}" && -f "${LIST_FASTQS_FOF}" ]]; then
  FOF="${LIST_FASTQS_FOF}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  for c in "${SLURM_SUBMIT_DIR}/src/preprocessing/ont_r10_pychopper/list_fastqs_for_map.fof" \
    "${SLURM_SUBMIT_DIR}/list_fastqs_for_map.fof"; do
    if [[ -f "$c" ]]; then
      FOF="$c"
      break
    fi
  done
fi
if [[ -z "${FOF}" && -f "${_script_dir}/list_fastqs_for_map.fof" ]]; then
  FOF="${_script_dir}/list_fastqs_for_map.fof"
fi

echo "=== Job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="
echo "OUT_TSV=${OUT_TSV}"

mkdir -p "${read_qc_dir}" "${logs_dir}"

count_fastq_reads() {
  local fq="$1"
  if [[ -z "${fq:-}" || ! -e "$fq" ]]; then
    echo 0
    return
  fi
  if [[ "$fq" =~ \.gz$ ]]; then
    gzip -cd -- "$fq" | awk 'END{print NR/4}'
  else
    awk 'END{print NR/4}' "$fq"
  fi
}

count_primary_mapped_bam() {
  local bam="$1"
  if [[ -z "${bam:-}" || ! -e "$bam" ]]; then
    echo 0
    return
  fi
  samtools view -@ "$THREADS" -c -F 2308 -- "$bam"
}

if [[ -z "${FOF}" || ! -f "${FOF}" ]]; then
  echo "ERROR: could not find list_fastqs_for_map.fof" >&2
  exit 1
fi
echo "Using: ${FOF}"

mapfile -t FASTQ_PATHS < "$FOF"

{
  printf '%s\t%s\t%s\t%s\t%s\n' sample ont_prim_aln ont_fastq pb_prim_aln pb_fastq
  for fq in "${FASTQ_PATHS[@]}"; do
    [[ -z "${fq// }" ]] && continue
    ori_name="$(basename "$fq")"
    sample="${ori_name%_for_map.fastq}"
    bam="${bam_dir}/${sample}_primary_aln_sorted.bam"
    ont_fq="$(count_fastq_reads "$fq")"
    ont_pr="$(count_primary_mapped_bam "$bam")"
    printf '%s\t%s\t%s\t0\t0\n' "$sample" "$ont_pr" "$ont_fq"
  done
} > "$OUT_TSV"

echo "Wrote: ${OUT_TSV}"
wc -l "$OUT_TSV"
