#!/bin/bash
#SBATCH --job-name=sqanti3_nih
#SBATCH --output=logs/sqanti3_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15gb
#SBATCH --qos=short
#SBATCH --time=03:00:00
#SBATCH --array=0-0

pwd; date
# genome 
# Define reference annotation, reference assembly and FASTQ for isoforms, cage and polyA
ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
sr_dir="/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/"
cage_data="/storage/gge/genomes/mouse_ref_NIH/reference_genome/lft_mm39_CAGE.bed"
polyA="/home/apadepe/lr_pipelines/SQANTI3/data/polyA_motifs/mouse_and_human.polyA_motif.txt"

# Remember that runs from 1 to 4 belong to the P1 and samples from 5-8 belong
# to P2. Use the appropiate bc-sample correspondance file.
bc_samples="bc_sample_P1.tsv"

# Create an array of files
readarray myarray < P1_gff.fofn

# Read the file corresponding to the task
gff=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $gff

# Get the path to the gtf file
reads_dir=$(dirname $gff)


# Get file name wo estension

reads_file=$(basename $gff)

filename="${reads_file%.*}"

echo $filename

# Get the matching  short reads sample
# Extract the barcode from the name
readarray -t -d "." arr <<< $gff
bc=$(echo ${arr[2]} | cut -f1 -d"_")

# Get the corresponding sample
sample=$(grep $bc $bc_samples | cut -f1)

echo $sample

# Get the corresponding sr sample
sj=${sr_dir}${sample}*SJ.out.tab
bam=${sr_dir}${sample}*.sortedByCoord.out.bam

# Write bam dir to file
echo $bam > tmp_$SLURM_ARRAY_TASK_ID.fofn

echo $sj
echo $bam

# Run SQANTI3
source activate sq3

python3 /home/apadepe/lr_pipelines/SQANTI3/sqanti3_qc.py \
	--skipORF --dir "${reads_dir}/run_SQANTI/${filename}" \
	--output "${filename}" -c $sj --min_ref_len "0" \
	--SR_bam "tmp_${SLURM_ARRAY_TASK_ID}.fofn" \
	--CAGE_peak "$cage_data" --polyA_motif_list "$polyA" \
	--report "both" \
	${gff} ${ref_annotation} ${assembly}

pwd;date
