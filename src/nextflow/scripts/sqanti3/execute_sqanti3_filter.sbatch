#!/bin/bash
#SBATCH --job-name=sq3_filter
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/sq3_filter_%A_%a.out

module load anaconda
source activate SQANTI3.env

# Input
sq3_input=$1 # Working directory
sqanti_location=$2
src=$3

n_cores="$SLURM_CPUS_PER_TASK"
task_number="$SLURM_ARRAY_TASK_ID"

# read corresponding line of config
read -r gtf sq3dir <<< $(sed -n "${task_number}p" $sq3_input)
output_prefix=$(basename "$gtf" .gtf)
output_location=$(dirname "$gtf")/"${output_prefix}_sq3_filter"

# Create TP and TN sets
python3 ${src}/define_ml_sets.py ${sq3dir}/${output_prefix}_classification.txt \
	--output_dir ${output_location}

# Run the filter
python3 ${sqanti_location}/sqanti3_filter.py ml \
	-d ${output_location} \
	--isoforms ${sq3dir}/${output_prefix}_corrected.fasta \
	--gtf ${gtf} \
	-o ${output_prefix} \
	--skip_report \
	${sq3dir}/${output_prefix}_classification.txt \
	-p ${output_location}/TP_list.txt \
	-n ${output_location}/TN_list.txt \
	-r ${output_location}/exclusion_list.txt
	
python3 ${src}/rever_9th_gtf.py ${output_location}/$(basename "$gtf" .gtf).filtered.gtf
mv ${output_location}/$(basename "$gtf" .gtf).filtered.gtf.tmp ${output_location}/$(basename "$gtf" .gtf).filtered.gtf