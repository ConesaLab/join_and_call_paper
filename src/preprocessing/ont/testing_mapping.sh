#!/bin/bash
#SBATCH --job-name=sqanti_ont
#SBATCH --output=../analysis/logs/sqanti_%A_%a.out 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --qos=short
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

#minimap2 -ax splice -uf --MD -t 4 /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa <(cat ../fastqs/merged_fastqs/B15_2_full_length.fastq ../fastqs/merged_fastqs/B15_2_rescued.fastq ../fastqs/merged_fastqs/B15_2_unclassified.fastq) > ../bams/test_spliceUF_B15_2.sam

#minimap2 -ax splice --MD -t 4 /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa <(cat ../fastqs/merged_fastqs/B15_2_full_length.fastq ../fastqs/merged_fastqs/B15_2_rescued.fastq ../fastqs/merged_fastqs/B15_2_unclassified.fastq) > ../bams/test_splice_B15_2.sam

#samtools view -bS -F0x900 ../bams/test_spliceUF_B15_2.sam | samtools sort -o ../bams/test_spliceUF_B15_2_primary_aln_sorted.bam

#spliced_bam2gff -t 1000000 -M ../bams/test_spliceUF_B15_2_primary_aln_sorted.bam > ../bams/test_spliceUF_B15_2_primary_aln.gff

#samtools view -bS -F0x900 ../bams/test_splice_B15_2.sam | samtools sort -o ../bams/test_splice_B15_2_primary_aln_sorted.bam

#spliced_bam2gff -t 1000000 -M ../bams/test_splice_B15_2_primary_aln_sorted.bam > ../bams/test_splice_B15_2_primary_aln.gff

# Remove the unstranded from the minimap output without UF
#(2023_carol) (2023_carol) cmonzo@master:/storage/gge/nih/Nanopore/bams
#$ awk '$7 != "." || NR < 3' test_splice_B15_2_primary_aln.gff > test_splice_B15_2_primary_aln_corrected.gff 

export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/home/cmonzo/software/cDNA_Cupcake/

#python3 /home/cmonzo/software/SQANTI3-5.2/sqanti3_qc.py --min_ref_len 0 --skipORF --dir ../bams/test_run_SQANTI/spliceUF --output spliceUF ../bams/test_spliceUF_B15_2_primary_aln.gff /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa

#python3 /home/cmonzo/software/SQANTI3-5.2/sqanti3_qc.py --min_ref_len 0 --skipORF --dir ../bams/test_run_SQANTI/splice --output splice ../bams/test_splice_B15_2_primary_aln_corrected.gff /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa

#minimap2 -ax splice -uf --MD -t 4 /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa <(cat ../fastqs/merged_fastqs/B15_2_full_length.fastq ../fastqs/merged_fastqs/B15_2_rescued.fastq) > ../bams/test_spliceUF_onlyFL_B15_2.sam

#samtools view -bS -F0x900 ../bams/test_spliceUF_onlyFL_B15_2.sam | samtools sort -o ../bams/test_spliceUF_B15_2_primary_aln_sorted_onlyFL.bam

#spliced_bam2gff -t 1000000 -M ../bams/test_spliceUF_B15_2_primary_aln_sorted_onlyFL.bam > ../bams/test_spliceUF_B15_2_primary_aln_onlyFL.gff

python3 /home/cmonzo/software/SQANTI3-5.2/sqanti3_qc.py --min_ref_len 0 --skipORF --dir ../analysis/run_SQANTI/BK20_80_4 --output BK20_80_4 ../gffs/BK20_80_4_primary_aln.gff /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa
