#!/bin/bash
#SBATCH --job-name=unzip_fastq
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o unzip_fastq_%A_%a.out

task_number="$SLURM_ARRAY_TASK_ID"

fastq_path="/home/fabianje/repos/documenting_NIH/fabian/data/ont/fastq"

# Unzip .fastq.gz files

# Create directory for unzipped fastq files if it doesn't exist
mkdir -p $fastq_path/full_length
mkdir -p $fastq_path/rescued

# Read samples.tsv and unzip each .fastq.gz file
# Skip the header line and read only the $task_number line
line=$(sed -n "$((task_number+1))p" samples.tsv)

# Extract the last column (sample file name)
sample_file=$(echo "$line" | awk '{print $NF}')
sample_rescued="${sample_file%_full_length.fastq.gz}_rescued.fastq.gz"

# Put the basename of sample_file without gz into a variable
unzipped_file_basename=$(basename "${sample_file%.gz}")
unzipped_file_basename_rescued=$(basename "${sample_rescued%.gz}")
echo "Unzipping ${sample_file} to: ${unzipped_file_basename}"
echo "Unzipping ${sample_rescued} to: ${unzipped_file_basename_rescued}"

# Unzip the .fastq.gz file into the "fastq" directory
gunzip -c "$sample_file" > "$fastq_path/full_length/${unzipped_file_basename}"
gunzip -c "$sample_rescued" > "$fastq_path/rescued/${unzipped_file_basename_rescued}"

cat "$fastq_path/full_length/${unzipped_file_basename}" "$fastq_path/rescued/${unzipped_file_basename_rescued}" > "$fastq_path/merged/${unzipped_file_basename%_full_length.fastq}_fl_r.fastq"

