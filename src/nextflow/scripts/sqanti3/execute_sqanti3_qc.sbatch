#!/bin/bash
#SBATCH --job-name=sq3
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/sq3_%A_%a.out

module load anaconda
source activate SQANTI3.env

# Input
sq3_input=$1 # Working directory
genome=$2 # Reference genome
annotation=$3 # Reference annotation
sqanti_location=$4
SCRIPT_DIR=$5
n_cores="$SLURM_CPUS_PER_TASK"
task_number="$SLURM_ARRAY_TASK_ID"

# read corresponding line of config
read -r gtf count <<< $(sed -n "${task_number}p" $sq3_input)
output_prefix=$(basename "$gtf" .gtf)
output_location=$(dirname "$gtf")/"${output_prefix}_sq3"

echo "gtf: ${gtf}"
echo "count: ${count}"

sq3_gtf="${gtf%.*}.expr.gtf"
sbatch --wait $SCRIPT_DIR/filter_gtf_count.sbatch $count $gtf

if [[ -z "$count" ]]; then
    # If count is empty
    echo "Running SQANTI3 QC WITHOUT count table"
    python $sqanti_location/sqanti3_qc.py \
        ${sq3_gtf} $annotation $genome \
        -o ${output_prefix} \
        -d ${output_location} \
        --skipORF \
        -t $n_cores
else
    echo "Running SQANTI3 QC with count table"
    python $sqanti_location/sqanti3_qc.py \
        ${sq3_gtf} $annotation $genome \
        -fl ${count} \
        -o ${output_prefix} \
        -d ${output_location} \
        --skipORF \
        -t $n_cores
fi
