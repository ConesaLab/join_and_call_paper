#!/bin/bash
#SBATCH --job-name=isoseq
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=15gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/isoseq_%A_%a.out

module load anaconda
source activate isoseq.4.0

# Input
WD="$1"
genome="$2"
metadata="$3"
out_iso_detect="$4"
out_quant="$5"
line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# Change working directory
# cd $WD

# Get sample metadata
read -r cond sample_id pool bam_file fq sr <<< "$(head -n $line_number_metadata $metadata | tail -n 1)"


out_isoseq="${out_iso_detect}/isoseq/${cond}/${sample_id}/"
out_quant="${out_quant}/isoseq/${cond}/${sample_id}/${sample_id}"


# Print Input
echo "Sample id: $sample_id"
echo "Condition: $cond"
echo "Pool: $pool"
echo "BAM reads: $bam_file"
echo "Genome: $genome"
echo "Out dir isoseq: $out_isoseq"
echo "Out dir quant: $out_quant"
echo "N cores: $n_cores"
echo ""

# Check if output directories exist
directories=("$out_iso_detect" $(dirname "$out_quant") $(dirname "$out_isoseq"))

for directory in "${directories[@]}"; do
	if [ ! -d "$directory" ]; then
		mkdir -p "$directory"
		echo "Directory created: $directory"
	else
		echo "Directory already exists: $directory"
	fi
done

# Run the steps in the isoseq pipeline

echo "Running clustering step"
mkdir -p ${out_isoseq}cluster
isoseq cluster2 -j $SLURM_CPUS_PER_TASK $bam_file ${out_isoseq}cluster/${sample_id}_clustered.bam

echo "Running mapping step"
mkdir -p ${out_isoseq}map
pbmm2 align --preset ISOSEQ --sort -j $SLURM_CPUS_PER_TASK $genome ${out_isoseq}cluster/${sample_id}_clustered.bam ${out_isoseq}map/${sample_id}_aligned.bam

echo "Running collapse step"
mkdir -p ${out_isoseq}collapse
isoseq collapse --do-not-collapse-extra-5exons -j $SLURM_CPUS_PER_TASK --log-level INFO --log-file ${out_isoseq}collapse/${sample_id}_info.log ${out_isoseq}map/${sample_id}_aligned.bam $bam_file ${out_isoseq}collapse/${sample_id}.gff
mv ${out_isoseq}collapse/${sample_id}.gff ${out_isoseq}collapse/${sample_id}.gtf

echo "Moving quantification files"

# Count the number of columns in the file
column_count=$(head -1 "${out_isoseq}collapse/${sample_id}.flnc_count.txt" | awk -F',' '{print NF}')

# Apply transformations based on the column count
if [ "$column_count" -eq 2 ]; then
  # Use pbid instead of superPBID for files with 2 columns (call&join)
  {
    echo -e "pbid\tcount_fl"
    tail -n +2 "${out_isoseq}collapse/${sample_id}.flnc_count.txt" | sed $'s/,/\t/g'
  } > "${out_quant}.flnc_count.txt"
else
  # Use superPBID for files with more than 2 columns
  cat "${out_isoseq}collapse/${sample_id}.flnc_count.txt" | sed $'s/,/\t/g' | sed 's/id/superPBID/g' > "${out_quant}.flnc_count.txt"
fi