#!/bin/bash
#SBATCH --job-name=flair_concat
#SBATCH --qos medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=40gb # Per processor memory
#SBATCH -t 7-00:00:00     # Walltime
#SBATCH -o log/flair_concat_%A_%a.out

module load anaconda

# source activate flair
# source activate flair_basic_conda_env
source activate flair_conda_env

# Input
WD="$1"
genome="$2"
samples_id_file="$3"
out_iso_detect="$4"
cond="$5"
fastq_fofn="$6"
stringent="$7"
line_number="$SLURM_ARRAY_TASK_ID"
n_cores="$SLURM_CPUS_PER_TASK"

cd $WD

# Prepare reads fastq files
fastq_list=$(cat $fastq_fofn | tr '\n' ',')
fastq_list=${fastq_list%,}

# Get BED file
read bed_file chrom_annotation < <(awk -F'\t' "NR==$line_number {print \$1, \$2}" "$samples_id_file")

out_dir="${out_iso_detect}/flair_collapse/${cond}/split_by_chrom"
if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
	echo "[FLAIR] Directory created: ${out_dir}"
fi

# Isoform detection
echo "[FLAIR] Run FLAIR collapse ${SLURM_ARRAY_TASK_ID}: ${cond}, ${bed_file}"
echo $genome
echo $chrom_annotation
echo $bed_file
echo $fastq_list
echo "${out_dir}/${cond}_split_${SLURM_ARRAY_TASK_ID}"
echo $n_cores

if [ "$stringent" = true ]; then
	echo "[FLAIR] Running FLAIR collapse stringent"
	flair collapse \
		-g $genome \
		--gtf $chrom_annotation \
		-q $bed_file \
		-r $fastq_list \
		--output "${out_dir}/${cond}_split_${SLURM_ARRAY_TASK_ID}" \
		--check_splice \
		--stringent \
		--threads $n_cores \
		--generate_map \
		--annotation_reliant generate
else
	echo "[FLAIR] Running FLAIR collapse"
	flair collapse \
		-g $genome \
		--gtf $chrom_annotation \
		-q $bed_file \
		-r $fastq_list \
		--output "${out_dir}/${cond}_split_${SLURM_ARRAY_TASK_ID}" \
		--check_splice \
		--threads $n_cores \
		--generate_map
fi

