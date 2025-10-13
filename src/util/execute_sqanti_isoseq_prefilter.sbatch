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

# sbatch --array=1-12 execute_sqanti_isoseq_prefilter.sbatch isoseq_prefilter.tsv /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf /home/fabianje/tools/SQANTI3_dev /home/fabianje/repos/documenting_NIH/fabian/src/nextflow/scripts/sqanti3 /home/fabianje/repos/documenting_NIH/fabian/src/util/isoseq_prefilter_sq3

# Input
sq3_input=$1 # Working directory
genome=$2 # Reference genome
annotation=$3 # Reference annotation
sqanti_location=$4
output=$5
n_cores="$SLURM_CPUS_PER_TASK"
task_number="$SLURM_ARRAY_TASK_ID"

# read corresponding line of config
read -r tissue sample gtf count <<< $(sed -n "${task_number}p" $sq3_input)
output_prefix=$(basename "$gtf" .gtf)
output_location=$output/"${output_prefix}_sq3"

echo "gtf: ${gtf}"
echo "count: ${count}"

mkdir -p $output_location

sq3_gtf="${output_location}/$(basename "$gtf" .gtf).expr.gtf"

sbatch --wait ./isoseq_filter_gtf_count.sbatch $count $gtf $sq3_gtf

echo "Running SQANTI3 QC with count table"
python $sqanti_location/sqanti3_qc.py \
    ${sq3_gtf} $annotation $genome \
    -fl ${count} \
    -o ${output_prefix} \
    -d ${output_location} \
    --skipORF \
    -t $n_cores

