#!/bin/bash
#SBATCH --job-name=sqanti_reads_ont_r10
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_reads_%A_%a.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs/sqanti_reads_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --qos=medium
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3
#
# Optional: write BAM->GFF and gffread GTF to a separate directory so existing
#   ${base_dir}/analysis/sqanti_reads/<SRR>.{gff,gtf} are not overwritten, and keep the .gff.
#   Example (submit from repo with list_fastqs.fof):
#     sbatch --export=ALL,SQANTI_READS_INTERMEDIATE_OUTDIR=/storage/gge/Fabian/ont_r10_sy5y/analysis/sqanti_reads_bam2gff_probe \
#       src/preprocessing/ont_r10/6_sqanti_reads.sh
# SQANTI3 still writes under analysis/sqanti_reads/<SRR>/ unless you change that separately.

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

echo "=== Job ${SLURM_JOB_ID} | Task ${SLURM_ARRAY_TASK_ID:-N/A} | $(hostname) | $(date) ==="

sqanti_dir="$HOME/tools/SQANTI3_dev"

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
ref_dir="${base_dir}/ref"
bam_dir="${base_dir}/bam"
outdir="${base_dir}/analysis/sqanti_reads"

ref_annotation="${ref_dir}/gencode.v49_SIRV.gtf"
assembly="${ref_dir}/GRCh38_SIRV.fa"

mkdir -p $outdir

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

ori_name=$(basename $file)
filename="${ori_name%.fastq}"
bam="${bam_dir}/${filename}_primary_aln_sorted.bam"

echo $filename

# BAM -> GFF -> GTF: default paths under outdir; optional alternate dir (no overwrite of originals)
intermediate_dir="${SQANTI_READS_INTERMEDIATE_OUTDIR:-}"
if [[ -n "${intermediate_dir}" ]]; then
  mkdir -p "${intermediate_dir}"
  gff_path="${intermediate_dir}/${filename}.gff"
  gtf_path="${intermediate_dir}/${filename}.gtf"
  echo "Using intermediate dir (GFF/GTF not written to main sqanti_reads root): ${intermediate_dir}"
else
  gff_path="${outdir}/${filename}.gff"
  gtf_path="${outdir}/${filename}.gtf"
fi

# Convert BAM to GFF, then GFF to GTF for SQANTI reads input
spliced_bam2gff -t 1000000 -M "${bam}" > "${gff_path}"

gffread "${gff_path}" -T -o "${gtf_path}"

# Remove GFF only when using default layout (saves space). Always keep .gff when using a separate intermediate dir.
if [[ -z "${intermediate_dir}" ]]; then
  rm -f "${gff_path}"
else
  echo "Kept intermediate GFF: ${gff_path}"
fi

python3 ${sqanti_dir}/sqanti3_qc.py \
    --skipORF --dir "${outdir}/${filename}" \
    --output "${filename}_reads" --min_ref_len "0" \
    --report "skip" --force_id_ignore \
    "${gtf_path}" ${ref_annotation} ${assembly}
