#!/bin/bash
#SBATCH --job-name=sqanti_ont
#SBATCH --output=../analysis/logs/sqanti_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45gb
#SBATCH --qos=short
#SBATCH --time=10:00:00
#SBATCH --array=0-23
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

# Create array of files
readarray myarray < list_fastqs.fof

# Read the file corresponding to the task
file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
bam_dir="/storage/gge/nih/Nanopore/bams/"
qc_dir="/storage/gge/nih/Nanopore/analysis/"
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"


# SQ input
isoforms_gff="/storage/gge/nih/Nanopore/gffs/${filename}_primary_aln.gff"

export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/

python3 /home/cmonzo/software/SQANTI3-5.2/sqanti3_qc.py --min_ref_len 0 --skipORF --dir "${qc_dir}/run_SQANTI/${filename}" --output "${filename}" ${isoforms_gff} ${ref_annotation} ${assembly}
