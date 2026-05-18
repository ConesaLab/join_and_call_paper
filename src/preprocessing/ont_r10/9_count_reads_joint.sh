#!/bin/bash
#SBATCH --job-name=ont_r10_read_counts
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/count_reads_joint_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/count_reads_joint_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#
# Writes read_numbers_joint.tsv in the same schema as src/util/count_reads.sbatch
# (mouse NIH): sample, ont_prim_aln, ont_fastq, pb_prim_aln, pb_fastq.
# PacBio columns are 0 (ONT-only SY5Y dataset).
#
# Submit (Slurm copies this script to spool; do not rely on the script file path):
#   cd /path/to/join_and_call_paper && sbatch src/preprocessing/ont_r10/9_count_reads_joint.sh
# or: cd .../join_and_call_paper/src/preprocessing/ont_r10 && sbatch 9_count_reads_joint.sh
# or: export LIST_FASTQS_FOF=/abs/path/list_fastqs.fof && sbatch ...
# or: copy list_fastqs.fof to ${BASE_DIR}/list_fastqs.fof (see BASE_DIR below).
#
# Requires: list_fastqs.fof and data under BASE_DIR.
#
# Strict mode after sourcing bashrc: with nounset enabled first, /etc/bashrc can fail on
# unset BASHRCSOURCED when sourcing ~/.bashrc.

LC_ALL=C

source ~/.bashrc
module load samtools

set -euo pipefail

BASE_DIR="/storage/gge/Fabian/ont_r10_sy5y"
BAMDIR="${BASE_DIR}/bam"
OUT_DIR="${BASE_DIR}/analysis/read_qc"
OUT_TSV="${OUT_DIR}/read_numbers_joint.tsv"

# list_fastqs.fof: Slurm runs a *copy* of this script under /var/spool/slurm/..., so
# dirname(BASH_SOURCE) is not the repo. Resolve FOF in order:
#   LIST_FASTQS_FOF -> SLURM_SUBMIT_DIR paths -> BASE_DIR -> script dir (local runs).
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
echo "OUT_TSV=${OUT_TSV}"

mkdir -p "${OUT_DIR}"
mkdir -p "${BASE_DIR}/analysis/logs"

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

# primary mapped = exclude unmapped(4), secondary(256), supplementary(2048) => 2308
count_primary_mapped_bam() {
  local bam="$1"
  if [[ -z "${bam:-}" || ! -e "$bam" ]]; then
    echo 0
    return
  fi
  samtools view -@ "$THREADS" -c -F 2308 -- "$bam"
}

if [[ -z "${FOF}" || ! -f "${FOF}" ]]; then
  echo "ERROR: could not find list_fastqs.fof (Slurm spool breaks paths next to this script)." >&2
  echo "  Try: cd <repo> && sbatch src/preprocessing/ont_r10/9_count_reads_joint.sh" >&2
  echo "  Or:  export LIST_FASTQS_FOF=/path/to/list_fastqs.fof" >&2
  echo "  Or:  cp list_fastqs.fof ${BASE_DIR}/" >&2
  exit 1
fi
echo "Using list_fastqs.fof: ${FOF}"

mapfile -t FASTQ_PATHS < "$FOF"

{
  printf '%s\t%s\t%s\t%s\t%s\n' sample ont_prim_aln ont_fastq pb_prim_aln pb_fastq
  for fq in "${FASTQ_PATHS[@]}"; do
    [[ -z "${fq// }" ]] && continue
    ori_name="$(basename "$fq")"
    sample="$ori_name"
    [[ "$sample" == *.gz ]] && sample="${sample%.gz}"
    sample="${sample%.fastq}"
    sample="${sample%.fq}"
    bam="${BAMDIR}/${sample}_primary_aln_sorted.bam"
    ont_fq="$(count_fastq_reads "$fq")"
    ont_pr="$(count_primary_mapped_bam "$bam")"
    printf '%s\t%s\t%s\t0\t0\n' "$sample" "$ont_pr" "$ont_fq"
  done
} > "$OUT_TSV"

echo "Wrote: ${OUT_TSV}"
wc -l "$OUT_TSV"
