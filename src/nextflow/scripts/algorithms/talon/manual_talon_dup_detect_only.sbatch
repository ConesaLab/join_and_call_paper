#!/bin/bash
#SBATCH --job-name=talon_dup_detect
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=16gb
#SBATCH -t 06:00:00     # Walltime (single awk pass over ~800M GTF; increase if needed)
#SBATCH -o /storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_dup_detect/log/%x_%j.out
#SBATCH -e /storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_dup_detect/log/%x_%j.err

set -euo pipefail

# Edit if your formatted GTF lives elsewhere (output of the "formatted" awk in reformat_annotation.sh).
# Override at submit time, e.g.: sbatch --export=ALL,FORMATTED=/path/to/*_formatted.gtf,OUT=/path/out.txt ...
FORMATTED="${FORMATTED:-/storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_reformat_gtf/annotation/tmp/gencode.v49_SIRV_fixed_exons_fixedsirv_reformatted_formatted.gtf}"
OUT="${OUT:-/storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_dup_detect/duplicate_genes.txt}"

mkdir -p "$(dirname "$OUT")" "/storage/gge/Fabian/ont_r10_sy5y/ref/talon_manual_dup_detect/log"

test -r "$FORMATTED" || { echo "Missing or unreadable: $FORMATTED" >&2; exit 1; }

echo "INPUT:  $FORMATTED"
echo "OUTPUT: $OUT"

# Same logic as reformat_annotation.sh (gene_name present on >=2 distinct seqnames, $3==gene only).
time awk 'BEGIN {FS=OFS="\t"}
$3 == "gene" {
  gene_name = ""
  if (match($0, /gene_name "[^"]+"/)) {
    gene_name = substr($0, RSTART + 10, RLENGTH - 11)
    gsub(/^"/, "", gene_name)
  }
  if (gene_name != "") seen[gene_name SUBSEP $1] = 1
}
END {
  for (k in seen) {
    split(k, a, SUBSEP)
    g = a[1]
    nchr[g]++
  }
  for (g in nchr) {
    if (nchr[g] > 1) print g
  }
}' "$FORMATTED" > "$OUT"

echo "line count:"
wc -l "$OUT"
echo "first 40 lines:"
head -40 "$OUT"
