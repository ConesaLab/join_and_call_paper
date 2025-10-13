#!/bin/bash
#SBATCH --job-name=tusco_exec
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=4gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/tusco_exec_%A_%a.out

module load anaconda
source activate tusco_env

# constants
sirv_gtf="/storage/gge/genomes/mouse_ref_NIH/reference_genome/SIRV/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf"

# Get the directory where this script is located
SCRIPT_DIR="$(pwd)"

# Set absolute paths for TUSCO files
tusco_mouse="${SCRIPT_DIR}/tusco_mouse.tsv"
tusco_brain="${SCRIPT_DIR}/tusco_mmu_brain.tsv"
tusco_kidney="${SCRIPT_DIR}/tusco_mmu_kidney.tsv"

# Check if config file exists
if [ ! -f "tusco_config.csv" ]; then
    echo "Error: tusco_config.csv not found"
    exit 1
fi

# Get the line corresponding to this array task ID (adding 1 to skip header)
line_number=$((SLURM_ARRAY_TASK_ID + 1))
line=$(sed -n "${line_number}p" tusco_config.csv)

# Read the fields from the line
IFS=, read -r classification_file data_type tool_name tissue_name sample_name <<< "$line"

echo "Processing array task $SLURM_ARRAY_TASK_ID: $sample_name ($tissue_name)"

# Determine which tissue-specific TUSCO file to use
if [[ "$tissue_name" == "B100K0" ]]; then
    tusco_tissue="$tusco_brain"
    tissue_type="brain"
elif [[ "$tissue_name" == "B0K100" ]]; then
    tusco_tissue="$tusco_kidney"
    tissue_type="kidney"
else
    echo "Error: Unknown tissue type: $tissue_name"
    exit 1
fi

# Create output directory structure for both full and tissue-specific analyses
for analysis_type in "full" "$tissue_type"; do
    out_dir="$SCRIPT_DIR/tusco_results/tusco_${analysis_type}/${data_type}/${tool_name}/${tissue_name}"
    mkdir -p "$out_dir"
done

# Run TUSCO with mouse reference
echo "Processing $sample_name with mouse reference..."
Rscript $SCRIPT_DIR/tusco_ind.R \
    "$classification_file" \
    "$sirv_gtf" \
    "$tusco_mouse" \
    "$SCRIPT_DIR/tusco_results/tusco_full" \
    "$data_type" \
    "$tool_name" \
    "$tissue_name" \
    "$sample_name"

# Run TUSCO with tissue-specific reference
echo "Processing $sample_name with ${tissue_type} reference..."
Rscript $SCRIPT_DIR/tusco_ind.R \
    "$classification_file" \
    "$sirv_gtf" \
    "$tusco_tissue" \
    "$SCRIPT_DIR/tusco_results/tusco_${tissue_type}" \
    "$data_type" \
    "$tool_name" \
    "$tissue_name" \
    "$sample_name"

echo "Completed processing array task $SLURM_ARRAY_TASK_ID: $data_type $tool_name $tissue_name $sample_name"

