#!/bin/bash
#SBATCH --job-name=flair_ind
#SBATCH --qos medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 7-00:00:00     # Walltime
#SBATCH -o log/flair_ind_%A_%a.out

module load anaconda

# source activate flair
# source activate flair_basic_conda_env
source activate flair_conda_env

# Input
WD="$1"
genome="$2"
annotation="$3"
samples_metadata="$4"
out_iso_detect="$5"
out_quant="$6"
use_sr="$7"
short_read_junctions="$8"
short_read_junctions_concat="$9"
stringent="${10}"
line_number_metadata=$((SLURM_ARRAY_TASK_ID + 1)) # Add 1 because of header
n_cores="$SLURM_CPUS_PER_TASK"

# Change working directory
cd $WD

# Get sample metadata
read -r cond sample_id pool bam_file fastq_file <<< "$(head -n $line_number_metadata $samples_metadata | tail -n 1)"
out_align="${out_iso_detect}/flair_align_correct/${cond}"
out_collapse="${out_iso_detect}/flair_collapse/${cond}"
out_quant="${out_quant}/${cond}"
read_manifest="$out_quant"/"$sample_id".read_manifest.tsv

sr_junctions=""

if [ "$use_sr" = true ]; then
    while IFS=$'\t' read -r sr_cond sr_sample sr_bam; do
        if [ "$sample_id" = "$sr_sample" ]; then
            sr_junctions="$sr_bam"
            break
        fi
    done < <(tail -n +2 "$short_read_junctions")
fi

sr_junctions_concat=""

if [ "$use_sr" = true ]; then
    while IFS=$'\t' read -r sr_cond sr_sample sr_bam; do
        if [ "$cond" = "$sr_cond" ]; then
            sr_junctions_concat="$sr_bam"
            break
        fi
    done < <(tail -n +2 "$short_read_junctions_concat")
fi

# Print Input
echo "[FLAIR] Sample id: $sample_id"
echo "[FLAIR] Condition: $cond"
echo "[FLAIR] Pool: $pool"
echo "[FLAIR] FASTQ reads: $fastq_file"
echo "[FLAIR] BAM reads: $bam_file"
echo "[FLAIR] Genome: $genome"
echo "[FLAIR] Annotation: $annotation"
echo "[FLAIR] Out dir align: $out_align"
echo "[FLAIR] Out dir collapse: $out_collapse"
echo "[FLAIR] Out dir quant: $out_quant"
echo "[FLAIR] N cores: $n_cores"
echo "[FLAIR] SR junctions: $sr_junctions"
echo "[FLAIR] SR junctions concat: $sr_junctions_concat"
echo ""

# Check if output directories exist
directories=("$out_align" "$out_collapse" "$out_quant")

for directory in "${directories[@]}"; do
	if [ ! -d "$directory" ]; then
		mkdir -p "$directory"
		echo "[FLAIR] Directory created: $directory"
	else
		echo "[FLAIR] Directory already exists: $directory"
	fi
done

# Run FLAIR for long-read isoform detection
echo "[FLAIR] Running FLAIR isoform detection: $sample_id"

# NOTE: NOT RUNNING FLAIR ALIGN

echo "[FLAIR] Convert BAM file to BED12"
bam2Bed12 $bam_file > "$out_align"/"$sample_id".aln.bed


if [ "$use_sr" = true ]; then
	echo "[FLAIR] Running FLAIR correct with short reads"
	flair correct \
		-q "$out_align"/"$sample_id".aln.bed \
		-f $annotation \
		-g $genome \
		--output "$out_align"/"$sample_id" \
		--threads $n_cores \
		--shortread $sr_junctions

	# concat
	flair correct \
		-q "$out_align"/"$sample_id".aln.bed \
		-f $annotation \
		-g $genome \
		--output "$out_align"/"$sample_id"_concat \
		--threads $n_cores \
		--shortread $sr_junctions_concat
else
	echo "[FLAIR] Running FLAIR correct"
	flair correct \
		-q "$out_align"/"$sample_id".aln.bed \
		-f $annotation \
		-g $genome \
		--output "$out_align"/"$sample_id" \
		--threads $n_cores
fi

if [ "$stringent" = true ]; then
	echo "[FLAIR] Running FLAIR collapse stringent"
	flair collapse \
		-g $genome \
		--gtf $annotation \
		-q "$out_align"/"$sample_id"_all_corrected.bed \
		-r $fastq_file \
		--output "$out_collapse"/"$sample_id" \
		--stringent \
		--check_splice \
		--threads $n_cores \
		--generate_map \
		--annotation_reliant generate
		
else
	echo "[FLAIR] Running FLAIR collapse"
	flair collapse \
		-g $genome \
		--gtf $annotation \
		-q "$out_align"/"$sample_id"_all_corrected.bed \
		-r $fastq_file \
		--output "$out_collapse"/"$sample_id" \
		--check_splice \
		--threads $n_cores \
		--generate_map
fi

# Run FLAIR for quantification
echo "[FLAIR] Running FLAIR quantification"

# Write FLAIR read manifest
# flair read manifest does not allow underscores in the first 3 fields
sample_id_no_underscores=$(echo $sample_id | tr -d '_')
cond_no_underscores=$(echo $cond | tr -d '_')
pool_no_underscores=$(echo $pool | tr -d '_')
printf "%s\t%s\t%s\t%s" $sample_id_no_underscores $cond_no_underscores $pool_no_underscores $fastq_file > $read_manifest

flair quantify \
	-i "$out_collapse"/"$sample_id".isoforms.fa \
	--reads_manifest "$read_manifest" \
	--threads $n_cores \
	--output "$out_quant"/"$sample_id" \
	--tpm \
	--check_splice \
	--stringent \
	--isoform_bed "$out_collapse"/"$sample_id".isoforms.bed
