#!/bin/bash
#SBATCH --job-name=stmerge
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem=10gb
#SBATCH -t 1-00:00:00
#SBATCH -o log/stmerge_%A_%a.out

module load anaconda
source activate stringtie

WD=$1
array_config=$2
line_number="$SLURM_ARRAY_TASK_ID"

cd "$WD"

IFS=$'\t' read -r name gtf_fofn <<< "$(sed -n "${line_number}p" "$array_config")"

out_dir="${WD}/${name}_STMERGE"
mkdir -p "$out_dir"

merged_gtf="${out_dir}/${name}_STMERGE.gtf"

# No -G: merge ONLY the per-sample reconstructed transcripts, do NOT
# include the reference annotation in the merged output.
# Disable ALL filters to preserve every input transcript:
#   -F 0 -T 0 -c 0: no FPKM/TPM/coverage thresholds
#   -f 0:            no minimum isoform fraction (default 0.01 drops rare isoforms)
#   -m 0:            no minimum transcript length (default 50)
#   -i:              keep transcripts with retained introns
stringtie --merge \
    -o "$merged_gtf" \
    -F 0 -T 0 -c 0 \
    -f 0 -m 0 -i \
    "$gtf_fofn"
