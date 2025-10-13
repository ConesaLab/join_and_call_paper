#!/bin/bash
#SBATCH --job-name=old_joins
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/old_joins_%A.out

module load samtools

ont_bams="/storage/gge/nih/Nanopore/data/bams"
ont_fastqs="/storage/gge/nih/Nanopore/data/fastqs/merged_fastqs"
isoseq_trimmed="/storage/gge/nih/PacBio_IsoSeq/merged_reads/old_mouse"
isoseq="/storage/gge/nih/PacBio_IsoSeq/merged_fl_reads"

ont_target="/storage/gge/Fabian/nih/data/ont/old_mice"
ont_target_partial="$ont_target/partial_joins"
isoseq_trimmed_target="/storage/gge/Fabian/nih/data/isoseq_trimmed/old_mice"
isoseq_trimmed_target_partial="$isoseq_trimmed_target/partial_joins"
isoseq_target="/storage/gge/Fabian/nih/data/isoseq/old_mice"
isoseq_target_partial="$isoseq_target/partial_joins"

mkdir -p $ont_target_partial/2samples
mkdir -p $ont_target_partial/3samples
mkdir -p $isoseq_trimmed_target_partial/2samples
mkdir -p $isoseq_trimmed_target_partial/3samples
mkdir -p $isoseq_target_partial/2samples
mkdir -p $isoseq_target_partial/3samples

# ONT
## bams
# samtools merge -o $ont_target/B15_primary_aln_sorted.bam $ont_bams/B151_primary_aln_sorted.bam $ont_bams/B152_primary_aln_sorted.bam $ont_bams/B153_primary_aln_sorted.bam $ont_bams/B154_primary_aln_sorted.bam
# samtools merge -o $ont_target_partial/3samples/B151152153_primary_aln_sorted.bam $ont_bams/B151_primary_aln_sorted.bam $ont_bams/B152_primary_aln_sorted.bam $ont_bams/B153_primary_aln_sorted.bam
# samtools merge -o $ont_target_partial/2samples/B151152_primary_aln_sorted.bam $ont_bams/B151_primary_aln_sorted.bam $ont_bams/B152_primary_aln_sorted.bam

## fastqs
# cat $ont_fastqs/B151.fastq.gz $ont_fastqs/B152.fastq.gz $ont_fastqs/B153.fastq.gz $ont_fastqs/B154.fastq.gz > $ont_target/B15.fastq.gz
# cat $ont_fastqs/B151.fastq.gz $ont_fastqs/B152.fastq.gz $ont_fastqs/B153.fastq.gz > $ont_target_partial/3samples/B151152153.fastq.gz
# cat $ont_fastqs/B151.fastq.gz $ont_fastqs/B152.fastq.gz > $ont_target_partial/2samples/B151152.fastq.gz

# ISOSEQ_TRIMMED
## bams
# samtools merge -o $isoseq_trimmed_target/B15.aln.sorted.bam $isoseq_trimmed/B151.aln.sorted.bam $isoseq_trimmed/B152.aln.sorted.bam $isoseq_trimmed/B153.aln.sorted.bam $isoseq_trimmed/B154.aln.sorted.bam
# samtools merge -o $isoseq_trimmed_target_partial/3samples/B151152153.aln.sorted.bam $isoseq_trimmed/B151.aln.sorted.bam $isoseq_trimmed/B152.aln.sorted.bam $isoseq_trimmed/B153.aln.sorted.bam
# samtools merge -o $isoseq_trimmed_target_partial/2samples/B151152.aln.sorted.bam $isoseq_trimmed/B151.aln.sorted.bam $isoseq_trimmed/B152.aln.sorted.bam

## fastqs
cat $isoseq_trimmed/B151.concat.fastq $isoseq_trimmed/B152.concat.fastq $isoseq_trimmed/B153.concat.fastq $isoseq_trimmed/B154.concat.fastq | gzip > $isoseq_trimmed_target/B15.concat.fastq.gz
cat $isoseq_trimmed/B151.concat.fastq $isoseq_trimmed/B152.concat.fastq $isoseq_trimmed/B153.concat.fastq | gzip > $isoseq_trimmed_target_partial/3samples/B151152153.concat.fastq.gz
cat $isoseq_trimmed/B151.concat.fastq $isoseq_trimmed/B152.concat.fastq | gzip > $isoseq_trimmed_target_partial/2samples/B151152.concat.fastq.gz

# ISOSEQ
## bams
samtools merge -o $isoseq_target/B15.concat.sorted.bam $isoseq/B151/B151.concat.sorted.bam $isoseq/B152/B152.concat.sorted.bam $isoseq/B153/B153.concat.sorted.bam $isoseq/B154/B154.concat.sorted.bam
samtools merge -o $isoseq_target_partial/3samples/B151152153.concat.sorted.bam $isoseq/B151/B151.concat.sorted.bam $isoseq/B152/B152.concat.sorted.bam $isoseq/B153/B153.concat.sorted.bam
samtools merge -o $isoseq_target_partial/2samples/B151152.concat.sorted.bam $isoseq/B151/B151.concat.sorted.bam $isoseq/B152/B152.concat.sorted.bam

## fastqs
# cat $isoseq/B151/B151.concat.fastq $isoseq/B152/B152.concat.fastq $isoseq/B153/B153.concat.fastq $isoseq/B154/B154.concat.fastq | gzip > $isoseq_target/B15.concat.fastq.gz
# cat $isoseq/B151/B151.concat.fastq $isoseq/B152/B152.concat.fastq $isoseq/B153/B153.concat.fastq | gzip > $isoseq_target_partial/3samples/B151152153.concat.fastq.gz
# cat $isoseq/B151/B151.concat.fastq $isoseq/B152/B152.concat.fastq | gzip > $isoseq_target_partial/2samples/B151152.concat.fastq.gz
