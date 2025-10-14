#!/bin/bash
#SBATCH --job-name=merge_fastqs
#SBATCH --output=merge_fastqs.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=25gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #Send e-mails
#SBATCH --mail-user=carolina.monzo@csic.es

source ~/.bashrc

# Command to make the merge_fastqs.txt
#(2023_carol) cmonzo@master:/storage/gge/nih/Nanopore/scripts
#$ for file in $(ls ../fastq_pass/); do echo "zcat ../fastq_pass/${file}/*.fastq.gz > ../fastq_pass/merged_fastqs/${file}_merged.fastq.gz"; done > merge_fastqs.txt

zcat ../fastq_pass/barcode01/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode01_merged.fastq.gz
zcat ../fastq_pass/barcode02/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode02_merged.fastq.gz
zcat ../fastq_pass/barcode03/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode03_merged.fastq.gz
zcat ../fastq_pass/barcode04/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode04_merged.fastq.gz
zcat ../fastq_pass/barcode05/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode05_merged.fastq.gz
zcat ../fastq_pass/barcode06/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode06_merged.fastq.gz
zcat ../fastq_pass/barcode07/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode07_merged.fastq.gz
zcat ../fastq_pass/barcode08/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode08_merged.fastq.gz
zcat ../fastq_pass/barcode09/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode09_merged.fastq.gz
zcat ../fastq_pass/barcode10/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode10_merged.fastq.gz
zcat ../fastq_pass/barcode11/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode11_merged.fastq.gz
zcat ../fastq_pass/barcode12/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode12_merged.fastq.gz
zcat ../fastq_pass/barcode13/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode13_merged.fastq.gz
zcat ../fastq_pass/barcode14/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode14_merged.fastq.gz
zcat ../fastq_pass/barcode15/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode15_merged.fastq.gz
zcat ../fastq_pass/barcode16/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode16_merged.fastq.gz
zcat ../fastq_pass/barcode17/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode17_merged.fastq.gz
zcat ../fastq_pass/barcode18/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode18_merged.fastq.gz
zcat ../fastq_pass/barcode19/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode19_merged.fastq.gz
zcat ../fastq_pass/barcode20/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode20_merged.fastq.gz
zcat ../fastq_pass/barcode21/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode21_merged.fastq.gz
zcat ../fastq_pass/barcode22/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode22_merged.fastq.gz
zcat ../fastq_pass/barcode23/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode23_merged.fastq.gz
zcat ../fastq_pass/barcode24/*.fastq.gz > ../fastq_pass/merged_fastqs/barcode24_merged.fastq.gz

exit
