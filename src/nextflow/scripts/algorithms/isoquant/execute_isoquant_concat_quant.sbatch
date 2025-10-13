#!/bin/bash
#SBATCH --job-name=isoquant
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=40gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/isoquant_%A_%a.out

module load anaconda
source activate isoquant

# Input
WD="$1"
genome="$2"
annotation="$3"
samples_metadata="$4"
concat_metadata="$5"
out_quant="$6"
data_type="$7" # options: nanopore, pacbio_ccs
fl_data="$8"
line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# if fl_data was not provided, set it to false
fl_data="${fl_data:-false}"

read -r cond concat_sample_id pool bam_file fastq_file <<< "$(head -n $line_number_metadata $concat_metadata | tail -n 1)"

# Change working directory
cd $WD

bam_files=""
sample_ids=""
# Get sample metadata
while IFS= read -r line; do
    [[ "$line" == \#* ]] && continue
    read -r curr_cond sample_id pool bam_file fastq_file <<< "$line"
    if [ "$cond" = "$curr_cond" ]; then
        bam_files="$bam_files $bam_file"
        sample_ids="$sample_ids $sample_id"
    fi
done < "$samples_metadata"

out_cond="${out_quant}/${cond}"

isoquant_gtf="${out_cond}/${concat_sample_id}/OUT/OUT.transcript_models.gtf"

# Print Input
echo "[isoquant] Condition: $cond"
echo "[isoquant] BAM reads: $bam_files"
echo "[isoquant] Genome: $genome"
echo "[isoquant] Annotation: $isoquant_gtf"
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
echo "[isoquant] Running isoquant for concat quantification: $cond"

if [ $fl_data = true ]; then
    isoquant.py -t $n_cores --reference $genome --genedb $isoquant_gtf --bam $bam_files --labels $sample_ids --data_type $data_type --fl_data --no_model_construction -o "${out_cond}/${cond}"
else
    isoquant.py -t $n_cores --reference $genome --genedb $isoquant_gtf --bam $bam_files --labels $sample_ids --data_type $data_type --no_model_construction -o "${out_cond}/${cond}"
fi
