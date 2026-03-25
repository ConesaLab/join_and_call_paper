#!/bin/bash
#SBATCH --job-name=flair_concat_quant
#SBATCH --qos medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=120gb # Per processor memory
#SBATCH -t 7-00:00:00     # Walltime
#SBATCH -o log/flair_concat_quant_%A_%a.out

module load anaconda

# source activate flair
# source activate flair_basic_conda_env
source activate flair_conda_env

# Input
WD="$1"
genome="$2"
samples_metadata="$3"
out_iso_detect="$4"
out_quant="$5"
sort_script="$6"
n_cores="$SLURM_CPUS_PER_TASK"
task_id="$SLURM_ARRAY_TASK_ID"

# select the condition specified by the task id
cond=$(awk -v task_id="$task_id" 'NR > 1 && !/^#/ { if (!seen[$1]++) entries[++count] = $1 } END { print entries[task_id] }' "$samples_metadata")

echo $cond

out_quant="${out_quant}/${cond}"
isoforms_gtf="${out_iso_detect}/flair_collapse/${cond}/${cond}.join.isoforms.unsorted.gtf"
isoforms_sorted_gtf="${out_iso_detect}/flair_collapse/${cond}/${cond}.join.isoforms.gtf"
isoforms_fa="${out_iso_detect}/flair_collapse/${cond}/${cond}.join.isoforms.fa"
isoforms_bed="${out_iso_detect}/flair_collapse/${cond}/${cond}.join.isoforms.bed"

cd $WD

# Concatenate files splitted by chrom
cat $(ls "$out_iso_detect"/flair_collapse/"$cond"/split_by_chrom/*.isoforms.gtf) > $isoforms_gtf
cat $(ls "$out_iso_detect"/flair_collapse/"$cond"/split_by_chrom/*.isoforms.fa) > $isoforms_fa
cat $(ls "$out_iso_detect"/flair_collapse/"$cond"/split_by_chrom/*.isoforms.bed) > $isoforms_bed

bash $sort_script -i $isoforms_gtf -o $isoforms_sorted_gtf

# Write read manifest for FLAIR quantification
read_manifest="${out_quant}/${cond}.read_manifest.tsv"
echo -n "" > $read_manifest

while IFS= read -r line; do
	if [[ "$line" == \#* ]]; then
		continue  # Skip lines that start with #
	fi

	# Split line into variables
	read -r curr_cond sample_id pool bam_file fastq_file <<< "$line"

	if [ "$cond" = "$curr_cond" ]; then
		# flair read manifest does not allow underscores in the first 3 fields
		sample_id_no_underscores=$(echo $sample_id | tr -d '_')
		cond_no_underscores=$(echo $cond | tr -d '_')
		pool_no_underscores=$(echo $pool | tr -d '_')
		printf "%s\t%s\t%s\t%s\n" $sample_id_no_underscores $cond_no_underscores $pool_no_underscores $fastq_file >> $read_manifest
	fi

done < "$samples_metadata"

echo "[FLAIR] isoforms_fa: $isoforms_fa"
echo "[FLAIR] read_manifest: $read_manifest"
echo "[FLAIR] n_cores: $n_cores"
echo "[FLAIR] out_quant: $out_quant"
echo "[FLAIR] cond: $cond"
echo "[FLAIR] isoforms_bed: $isoforms_bed"

# Run FLAIR quantification
flair quantify \
	-i $isoforms_fa \
	--reads_manifest $read_manifest \
	--threads $n_cores \
	--output "$out_quant"/"$cond"_ConcatReads \
	--tpm \
	--isoform_bed $isoforms_bed \
	--temp_dir "./tmp"

	# running with check_splice and/or stringent here causes an error
	# --check_splice \
	# --stringent \

echo "[FLAIR] End"
