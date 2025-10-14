#!/bin/bash
#SBATCH --job-name=sqanti_ont
#SBATCH --output=logs/sqanti_ont_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --array=0-9

sqanti_dir="/home/apadepe/lr_pipelines/SQANTI3"
ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
outdir="/storage/gge/Alejandro/nih/nanopore/reads"

mkdir -p $outdir
source activate test_is3

# Create an array of files
readarray myarray < ont_mapped_reads.fofn

# Read the file corresponding to the task
bam=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $bam

# Get file name wo estension
sample=$(basename $bam _primary_aln_sorted.bam)


echo $sample


# spliced_bam2gff -t 1000000 -M ${bam} > ${outdir}/${sample}.gff

# gffread ${outdir}/${sample}.gff -T -o ${outdir}/${sample}.gtf

# rm ${outdir}/${sample}.gff

conda deactivate
conda activate SQANTI3.env

python3 ${sqanti_dir}/sqanti3_qc.py \
	--skipORF --dir "${outdir}/${sample}" \
	--output "${sample}" --min_ref_len "0" \
	--report "skip" --force_id_ignore \
	${outdir}/${sample}.gtf ${ref_annotation} ${assembly}

