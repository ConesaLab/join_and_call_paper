#!/bin/bash
#SBATCH --job-name=get_UJC_Iso-seq
#SBATCH --time 10:00:00
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --output=logs/get_UJC.out
#SBATCH --error=logs/get_UJC.err
#SBATCH --array=0-23

# Script adapated from Nathalia doucmenting_nih to run in paralel
# This script takes the .gtf from SQANTI files from to create UJC per sample
# UJC are generated for chr1-19, X and Y (default in gtftools) - discarding noncanonical chr. See parameter -c (https://www.genemine.org/gtftools.php) to include other contigs (chrM, unlocalized, random contigs...).
# gtftools installation https://bioconda.github.io/recipes/gtftools/README.html

#module purge && module load anaconda
source activate gtftools 
module load bedtools

# outdir
outdir="/storage/gge/Alejandro/nih/iso_seq/combined_UJC"
mkdir $outdir

readarray -t  myarray < isoseq_gtf.fofn
gtf_file=${myarray[$SLURM_ARRAY_TASK_ID]}

# Get prefix name
name=`basename ${gtf_file} | sed 's/.gtf//g'`

# Get introns
gtftools -i ${outdir}/${name}_introns.bed ${gtf_file}

# Get UJC 
awk -F'\t' -v OFS="\t" '{print $5,"chr"$1,$4,$2+1"_"$3}' ${outdir}/${name}_introns.bed | bedtools groupby -g 1 -c 2,3,4 -o distinct,distinct,collapse | sed 's/,/_/g' | awk -F'\t' -v OFS="\t" '{print $1,$2"_"$3"_"$4}' > ${outdir}/${name}_UJC.txt

# remove intermidiate files
rm ${outdir}/${name}_introns.bed
rm ${gtf_file}.ensembl