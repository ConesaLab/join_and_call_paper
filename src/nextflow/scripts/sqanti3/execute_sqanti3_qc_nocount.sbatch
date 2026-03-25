#!/bin/bash
#SBATCH --job-name=sq3_nc
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mem=20gb
#SBATCH -t 1-00:00:00
#SBATCH -o log/sq3_nc_%A_%a.out

# Simplified SQANTI3 QC execution without count-based GTF filtering.
# Input file has one GTF path per line (no count column).
# Skips filter_gtf_count.sbatch and runs SQANTI3 directly on the GTF.

module load anaconda
source activate SQANTI3.env

sq3_input=$1
genome=$2
annotation=$3
sqanti_location=$4
n_cores="$SLURM_CPUS_PER_TASK"
task_number="$SLURM_ARRAY_TASK_ID"

gtf=$(sed -n "${task_number}p" "$sq3_input")
output_prefix=$(basename "$gtf" .gtf)
output_location=$(dirname "$gtf")/"${output_prefix}_sq3"

echo "gtf: ${gtf}"
echo "output_prefix: ${output_prefix}"
echo "output_location: ${output_location}"

python "$sqanti_location/sqanti3_qc.py" \
    "$gtf" "$annotation" "$genome" \
    -o "${output_prefix}" \
    -d "${output_location}" \
    --skipORF \
    -t "$n_cores"
