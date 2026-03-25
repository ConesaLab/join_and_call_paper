#!/bin/bash
#SBATCH --job-name=tama_merge
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=10gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/tama_merge_%A_%a.out

module load anaconda
source activate tama

# Input
WD=$1
array_config=$2
tama_path=$3

tama_format_converter_path="${tama_path}/tama_go/format_converter"
n_cores="$SLURM_CPUS_PER_TASK"
line_number="$SLURM_ARRAY_TASK_ID"

# Change working directory
cd $WD

IFS=$'\t' read -r name gtf_fofn <<< "$(sed -n "${line_number}p" $array_config)"

out_path="${WD}/${name}_TAMA"
bed_dir="${out_path}/beds"
merged_bed="${out_path}/${name}_TAMA"
bed_fofn="${merged_bed}.bed.fofn"

# create directory
for dir in "$out_path" "$bed_dir"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "Directory created: $dir"
    else
        echo "Directory already exists: $dir"
    fi
done

if [ -f $bed_fofn ]; then
	rm $bed_fofn
fi

# GTF to BED12
while IFS= read -r line; do
	tama_bed="${bed_dir}/$(basename "${line%.gtf}").tama.bed"
	# echo "$tama_bed" | awk -F/ '{split($NF, a, "."); gsub(".gtf", "", a[1]); print $0 "\tno_cap\t1,1,1\t" a[1]}' >> $bed_fofn
	# echo "$tama_bed" | awk -F/ '{split($NF, a, "."); gsub(".gtf", "", a[1]); print $0 "\tcapped\t1,1,1\t" a[1]}' >> $bed_fofn
	echo "$tama_bed" | awk -F/ '{split($NF, a, "."); gsub(".gtf", "", a[1]); split(a[1], b, "_"); print $0 "\tcapped\t1,1,1\t" b[1]}' >> $bed_fofn
	echo "$line $tama_bed"

	gtf_filtered="${line/.gtf/.filtered.gtf}"
	grep -E 'gene_id.*transcript_id|transcript_id.*gene_id' "$line" > "$gtf_filtered"

	# python ${tama_format_converter_path}/tama_format_gtf_to_bed12_ncbi.py \
	python ${tama_format_converter_path}/tama_format_gff_to_bed12_cupcake.py \
		$gtf_filtered \
		$tama_bed
done < $gtf_fofn

# Run TAMA merge
python $tama_path/tama_merge.py \
	-f $bed_fofn \
	-p $merged_bed \
	-a 50 \
	-z 50 \
	-m 0 \
	-d merge_dup

# BED to GTF
python $tama_format_converter_path/tama_convert_bed_gtf_ensembl_no_cds.py \
	$merged_bed.bed \
	$merged_bed.gtf
