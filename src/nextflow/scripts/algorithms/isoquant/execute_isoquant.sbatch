#!/bin/bash
#SBATCH --job-name=isoquant
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=80gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/isoquant_%A_%a.out

module load anaconda
source activate isoquant

# Input
WD="$1"
genome="$2"
annotation="$3"
samples_metadata="$4"
out_iso_detect="$5"
data_type="$6" # options: nanopore, pacbio_ccs
fl_data="$7"
line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# if fl_data was not provided, set it to false
fl_data="${fl_data:-false}"


# Change working directory
cd $WD

# Get sample metadata
read -r cond sample_id pool bam_file fastq_file <<< "$(head -n $line_number_metadata $samples_metadata | tail -n 1)"
out_cond="${out_iso_detect}/${cond}"

# Print Input
echo "[isoquant] Sample id: $sample_id"
echo "[isoquant] Condition: $cond"
echo "[isoquant] Pool: $pool"
echo "[isoquant] FASTQ reads: $fastq_file"
echo "[isoquant] BAM reads: $bam_file"
echo "[isoquant] Genome: $genome"
echo "[isoquant] Annotation: $annotation"
echo "[isoquant] Out dir: $out_cond"
echo "[isoquant] N cores: $n_cores"
echo ""

if [ ! -d "$out_cond" ]; then
	mkdir -p "$out_cond"
	echo "[isoquant] Directory created: $out_cond"
else
	echo "[isoquant] Directory already exists: $out_cond"
fi

# Run isoquant for long-read isoform detection
echo "[isoquant] Running isoquant isoform detection: $sample_id"

if [ $fl_data = true ]; then
    isoquant.py -t $n_cores --reference $genome --genedb $annotation --bam $bam_file --data_type $data_type --fl_data -o "${out_cond}/${sample_id}"
else
    isoquant.py -t $n_cores --reference $genome --genedb $annotation --bam $bam_file --data_type $data_type -o "${out_cond}/${sample_id}"
fi
