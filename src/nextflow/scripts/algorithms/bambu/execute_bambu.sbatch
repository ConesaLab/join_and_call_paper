#!/bin/bash
#SBATCH --job-name=bambu
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=30gb # Per processor memory
#SBATCH -t 12:00:00     # Walltime
#SBATCH -o log/bambu_%A_%a.out

module load anaconda
source activate bambu_env


# Input
genome="$1"
annotation="$2"
metadata="$3"
out_iso_detect="$4"
out_quant="$5"


line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# Change working directory
# cd $WD

# Get sample metadata
read -r cond sample_id pool bam_file fastq_file <<< "$(head -n $line_number_metadata $metadata | tail -n 1)"

out_quant="${out_quant}/${cond}/${sample_id}/${sample_id}"
out_bambu="${out_iso_detect}/bambu/${cond}/${sample_id}/${sample_id}"


# Print Input
echo "Sample id: $sample_id"
echo "Condition: $cond"
echo "Pool: $pool"
echo "Genome: $genome"
echo "Annotation: $annotation"
echo "Out dir bambu: $out_bambu"
echo "Out dir quant: $out_quant"
echo "N cores: $n_cores"
echo ""

# Check if output directories exist
directories=("$out_iso_detect" $(dirname "$out_quant")  $(dirname "$out_bambu"))

for directory in "${directories[@]}"; do
	if [ ! -d "$directory" ]; then
		mkdir -p "$directory"
		echo "Directory created: $directory"
	else
		echo "Directory already exists: $directory"
	fi
done

# Run bambu
Rscript $SCRIPT_DIR/bambu.R $genome $annotation $bam_file $out_bambu $n_cores
echo "bambu finished writing the outputs"

# Move counts to final location
echo "moving counts"
mv ${out_bambu}_counts.tsv ${out_quant}_counts.tsv

# reorder gtf
echo "formating bambu gtf"
# echo "Original bambu:"
# head ${out_bambu}.gtf
python $SCRIPT_DIR/fix_bambu_gtf.py ${out_bambu}.gtf ${out_bambu}.gtf.tmp
# echo "fixed genes"
# head ${out_bambu}.gtf.tmp
gffread ${out_bambu}.gtf.tmp -T -o ${out_bambu}.gtf
# echo "sorted gtf"
# head ${out_bambu}.gtf
python $REFORMAT_SCRIPT ${out_bambu}.gtf
# echo "ordered gene_id"
# head ${out_bambu}.gtf.tmp
mv ${out_bambu}.gtf.tmp ${out_bambu}.gtf
