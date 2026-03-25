#!/bin/bash
#SBATCH --job-name=count_reads_joint
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=10gb
#SBATCH -t 1-00:00:00
#SBATCH -o log/count_reads_joint_%A.out

set -euo pipefail
LC_ALL=C

module load samtools

# ---- inputs ----
ONT_TSV="/storage/gge/Fabian/nih/data/metadata/ont_samples.tsv"
PB_TSV="/storage/gge/Fabian/nih/data/metadata/isoseq_fl_samples.tsv"

# ---- output ----
OUT_TSV="./read_numbers_joint.tsv"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---- helpers ----
count_fastq_reads () {
  local fq="$1"
  if [[ -z "${fq:-}" || ! -e "$fq" ]]; then
    echo 0; return
  fi
  if [[ "$fq" =~ \.gz$ ]]; then
    gzip -cd -- "$fq" | awk 'END{print NR/4}'
  else
    awk 'END{print NR/4}' "$fq"
  fi
}

count_primary_mapped_bam () {
  local bam="$1"
  if [[ -z "${bam:-}" || ! -e "$bam" ]]; then
    echo 0; return
  fi
  # primary mapped = exclude unmapped(4), secondary(256), supplementary(2048) => 2308
  samtools view -@ "$THREADS" -c -F 2308 -- "$bam"
}

# ---- storage ----
declare -A ONT_PRIM ONT_FQ PB_PRIM PB_FQ SEEN

# ---- read ONT (no pipelines into while!) ----
if [[ -e "$ONT_TSV" ]]; then
  while IFS=$'\t' read -r condition sample pool bam fastq; do
    [[ -z "${sample:-}" || "$sample" == "sample" ]] && continue  # skip header/blank
    ONT_PRIM["$sample"]="$(count_primary_mapped_bam "$bam")"
    ONT_FQ["$sample"]="$(count_fastq_reads "$fastq")"
    SEEN["$sample"]=1
  done < <(awk -F'\t' 'NR==1{next} NF>=5' "$ONT_TSV")
else
  echo "WARN: ONT metadata not found at $ONT_TSV" >&2
fi

# ---- read PacBio ----
if [[ -e "$PB_TSV" ]]; then
  while IFS=$'\t' read -r condition sample pool bam fastq; do
    [[ -z "${sample:-}" || "$sample" == "sample" ]] && continue
    PB_PRIM["$sample"]="$(count_primary_mapped_bam "$bam")"
    PB_FQ["$sample"]="$(count_fastq_reads "$fastq")"
    SEEN["$sample"]=1
  done < <(awk -F'\t' 'NR==1{next} NF>=5' "$PB_TSV")
else
  echo "WARN: PacBio metadata not found at $PB_TSV" >&2
fi

# ---- write output ----
echo -e "sample\tont_prim_aln\tont_fastq\tpb_prim_aln\tpb_fastq" > "$OUT_TSV"

# iterate all samples (sorted, stable)
for sample in $(printf "%s\n" "${!SEEN[@]}" | sort); do
  printf "%s\t%s\t%s\t%s\t%s\n" \
    "$sample" \
    "${ONT_PRIM[$sample]:-0}" \
    "${ONT_FQ[$sample]:-0}" \
    "${PB_PRIM[$sample]:-0}" \
    "${PB_FQ[$sample]:-0}" >> "$OUT_TSV"
done

echo "Wrote: $OUT_TSV"
