#!/bin/bash
#SBATCH --job-name=mandalorion
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=48gb # Per processor memory
#SBATCH -t 7-00:00:00     # Walltime
#SBATCH -o log/mandalorion_%A_%a.out

module load anaconda
source activate mandalorion

# Input
WD="$1"
genome="$2"
annotation="$3"
samples_metadata="$4"
out_iso_detect="$5"
mandalorion_location="$6"
line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# Change working directory
cd $WD

# Get sample metadata
read -r cond sample_id pool bam_file fastq_file <<< "$(head -n $line_number_metadata $samples_metadata | tail -n 1)"
out_cond="${out_iso_detect}/${cond}"

# Print Input
echo "[mandalorion] Sample id: $sample_id"
echo "[mandalorion] Condition: $cond"
echo "[mandalorion] Pool: $pool"
echo "[mandalorion] FASTQ reads: $fastq_file"
echo "[mandalorion] BAM reads: $bam_file"
echo "[mandalorion] Genome: $genome"
echo "[mandalorion] Annotation: $annotation"
echo "[mandalorion] Out dir: $out_cond"
echo "[mandalorion] N cores: $n_cores"
echo ""

if [ ! -d "$out_cond" ]; then
	mkdir -p "$out_cond"
	echo "[mandalorion] Directory created: $out_cond"
else
	echo "[mandalorion] Directory already exists: $out_cond"
fi

# Run mandalorion for long-read isoform detection
echo "[mandalorion] Running mandalorion isoform detection: $sample_id"

python $mandalorion_location/Mando.py -p "${out_cond}/${sample_id}" -g $annotation -G $genome -f $fastq_file -t $n_cores -W '[SIRV,"basic"]'
