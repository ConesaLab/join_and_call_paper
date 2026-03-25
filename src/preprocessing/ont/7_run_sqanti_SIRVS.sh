#!/bin/bash
#SBATCH --job-name=sqanti_SIRV_ont
#SBATCH --output=../analysis/logs/sqanti_SIRV_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --qos=short
#SBATCH --time=10:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

# SIRV subsetted GTF files extracted like this:
#(flair) cmonzo@master:/storage/gge/nih/Nanopore/gffs
#$ for file in $(ls *_primary_aln.gff); do head -n1 $file > ${file%_primary_aln.gff}_SIRV.gff; grep SIRV $file >> ${file%_primary_aln.gff}_SIRV.gff; done

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
qc_dir="/storage/gge/nih/Nanopore/analysis/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"


# SQ input
isoforms_gff="/storage/gge/nih/Nanopore/gffs/${filename}_SIRV.gff"

export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/

python3 /home/cmonzo/software/SQANTI3-5.2/sqanti3_qc.py --min_ref_len 0 --skipORF --dir "${qc_dir}/run_SQANTI/${filename}_SIRV" --output "${filename}_SIRV" ${isoforms_gff} ${ref_annotation} ${assembly}
