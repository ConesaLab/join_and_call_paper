#!/bin/bash
#SBATCH --job-name=seq_len_only
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 1-00:00:00
#SBATCH -o log/seq_len_only_%A.out

set -euo pipefail
LC_ALL=C

module load samtools

# ---- inputs ----
ONT_TSV="/storage/gge/Fabian/nih/data/metadata/ont_samples.tsv"
PB_TSV="/storage/gge/Fabian/nih/data/metadata/isoseq_fl_samples.tsv"

# ---- output root ----
OUT_DIR="./lengths_seq_only"
mkdir -p "${OUT_DIR}/ont" "${OUT_DIR}/pb"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---- helpers ----

# One number per FASTQ read: sequence length.
write_fastq_lengths () {
  local fq="$1" out="$2"
  if [[ ! -e "$fq" ]]; then
    echo "WARN: FASTQ not found: $fq" >&2
    : > "$out"
    return
  fi
  if [[ "$fq" =~ \.gz$ ]]; then
    gzip -cd -- "$fq" | awk 'NR%4==2{print length($0)}' > "$out"
  else
    awk 'NR%4==2{print length($0)}' "$fq" > "$out"
  fi
}

# One number per BAM alignment: SEQ field length (primary mapped only).
# If SEQ is "*" (rare), print 0.
write_bam_seq_lengths () {
  local bam="$1" out="$2"
  if [[ ! -e "$bam" ]]; then
    echo "WARN: BAM not found: $bam" >&2
    : > "$out"
    return
  fi
  # primary mapped only: exclude unmapped(4), secondary(256), supplementary(2048) => -F 2308
  samtools view -@ "$THREADS" -F 2308 -- "$bam" \
  | awk '{s=$10; if(s=="*"||s==""){print 0}else{print length(s)}}' > "$out"
}

process_table () {
  local tsv="$1" subdir="$2"
  [[ -e "$tsv" ]] || { echo "WARN: Missing table $tsv" >&2; return; }

  # columns: condition  sample  pool  bam  fastq
  while IFS=$'\t' read -r condition sample pool bam fastq; do
    [[ -z "${sample:-}" || "$sample" == "sample" ]] && continue

    bam_out="${OUT_DIR}/${subdir}/${sample}.bam.readlen.txt"
    fq_out="${OUT_DIR}/${subdir}/${sample}.fastq.readlen.txt"

    echo "[${subdir}] $sample -> $bam_out ; $fq_out"

    write_bam_seq_lengths   "$bam"  "$bam_out"
    write_fastq_lengths     "$fastq" "$fq_out"

  done < <(awk -F'\t' 'NR==1{next} NF>=5' "$tsv")
}

# ---- run ----
process_table "$ONT_TSV" "ont"
process_table "$PB_TSV"  "pb"

echo "Done. Outputs in: $OUT_DIR"
