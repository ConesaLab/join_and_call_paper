#!/bin/bash
#SBATCH --job-name=isoquant_old
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=40gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/isoquant_old_%A_%a.out

module load anaconda
source activate isoquant_old

# Input
WD="$1"
bam_file="$2"
sample_id="$3"
genome="$4"
annotation="$5"
out_iso_detect="$6"
data_type="$7" # options: nanopore, pacbio_ccs
fl_data="$8"
prefix="$9"
n_cores="$SLURM_CPUS_PER_TASK"

# if fl_data was not provided, set it to false
fl_data="${fl_data:-false}"

# Change working directory
cd $WD

# Print Input
echo "[isoquant] Sample id: $sample_id"
echo "[isoquant] BAM reads: $bam_file"
echo "[isoquant] Genome: $genome"
echo "[isoquant] Annotation: $annotation"
echo "[isoquant] Out dir: $out_iso_detect"
echo "[isoquant] N cores: $n_cores"
echo ""

if [ ! -d "$out_iso_detect" ]; then
	mkdir -p "$out_iso_detect"
	echo "[isoquant] Directory created: $out_iso_detect"
else
	echo "[isoquant] Directory already exists: $out_iso_detect"
fi

# Run isoquant for long-read isoform detection
echo "[isoquant] Running isoquant isoform detection: $sample_id"

if [ $fl_data = true ]; then
    isoquant.py -t $n_cores --reference $genome --genedb $annotation --bam $bam_file --data_type $data_type -p $prefix --fl_data -o "${out_iso_detect}/${sample_id}"
else
    isoquant.py -t $n_cores --reference $genome --genedb $annotation --bam $bam_file --data_type $data_type -p $prefix -o "${out_iso_detect}/${sample_id}"
fi

