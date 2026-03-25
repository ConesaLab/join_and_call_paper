#!/bin/bash
#SBATCH --job-name=sqanti_reads_ont_r10
#SBATCH --output=../analysis/logs/sqanti_reads_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

sqanti_dir="/home/cmonzo/software/SQANTI3-5.2"

base_dir="/storage/gge/Fabian/ont_r10_sy5y"
ref_dir="${base_dir}/ref"
bam_dir="${base_dir}/bam"
outdir="${base_dir}/analysis/sqanti_reads"

ref_annotation="${ref_dir}/gencode.v49_SIRV.gtf"
assembly="${ref_dir}/GRCh38_SIRV.fa"

mkdir -p $outdir

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

ori_name=$(basename $file)
filename="${ori_name%.fastq}"
bam="${bam_dir}/${filename}_primary_aln_sorted.bam"

echo $filename

# Convert BAM to GFF, then GFF to GTF for SQANTI reads input
spliced_bam2gff -t 1000000 -M ${bam} > ${outdir}/${filename}.gff

gffread ${outdir}/${filename}.gff -T -o ${outdir}/${filename}.gtf

rm ${outdir}/${filename}.gff

export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/

python3 ${sqanti_dir}/sqanti3_qc.py \
    --skipORF --dir "${outdir}/${filename}" \
    --output "${filename}" --min_ref_len "0" \
    --report "skip" --force_id_ignore \
    ${outdir}/${filename}.gtf ${ref_annotation} ${assembly}
