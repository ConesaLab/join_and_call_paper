#!/bin/bash
#SBATCH --job-name=partial_joins
#SBATCH --qos short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/partial_joins_%A.out

module load samtools

ont_bams="/storage/gge/nih/Nanopore/data/bams"
ont_fastqs="/storage/gge/nih/Nanopore/data/fastqs/merged_fastqs"
isoseq_trimmed_brain="/storage/gge/nih/PacBio_IsoSeq/merged_reads/mouse_B100K0"
isoseq_trimmed_kidney="/storage/gge/nih/PacBio_IsoSeq/merged_reads/mouse_B0K100"
isoseq="/storage/gge/nih/PacBio_IsoSeq/merged_fl_reads"

ont_target="/storage/gge/Fabian/nih/data/ont/partial_joins/"
isoseq_trimmed_target="/storage/gge/Fabian/nih/data/isoseq_trimmed/partial_joins/"
isoseq_target="/storage/gge/Fabian/nih/data/isoseq/partial_joins/"

mkdir -p $ont_target/2samples
mkdir -p $ont_target/3samples
mkdir -p $ont_target/4samples
mkdir -p $isoseq_trimmed_target/2samples
mkdir -p $isoseq_trimmed_target/3samples
mkdir -p $isoseq_trimmed_target/4samples
mkdir -p $isoseq_target/2samples
mkdir -p $isoseq_target/3samples
mkdir -p $isoseq_target/4samples


# MERGE ONT
# merge ont brain bams
# samtools merge -o $ont_target/4samples/B31323334_primary_aln_sorted.bam $ont_bams/B31_primary_aln_sorted.bam $ont_bams/B32_primary_aln_sorted.bam $ont_bams/B33_primary_aln_sorted.bam $ont_bams/B34_primary_aln_sorted.bam
# samtools merge -o $ont_target/3samples/B313233_primary_aln_sorted.bam $ont_bams/B31_primary_aln_sorted.bam $ont_bams/B32_primary_aln_sorted.bam $ont_bams/B33_primary_aln_sorted.bam
# samtools merge -o $ont_target/2samples/B3132_primary_aln_sorted.bam $ont_bams/B31_primary_aln_sorted.bam $ont_bams/B32_primary_aln_sorted.bam

# # merge ont kidney bams
# samtools merge -o $ont_target/4samples/K31323334_primary_aln_sorted.bam $ont_bams/K31_primary_aln_sorted.bam $ont_bams/K32_primary_aln_sorted.bam $ont_bams/K33_primary_aln_sorted.bam $ont_bams/K34_primary_aln_sorted.bam
# samtools merge -o $ont_target/3samples/K313233_primary_aln_sorted.bam $ont_bams/K31_primary_aln_sorted.bam $ont_bams/K32_primary_aln_sorted.bam $ont_bams/K33_primary_aln_sorted.bam
# samtools merge -o $ont_target/2samples/K3132_primary_aln_sorted.bam $ont_bams/K31_primary_aln_sorted.bam $ont_bams/K32_primary_aln_sorted.bam

# # merge ont brain fastqs
# cat $ont_fastqs/B31.fastq.gz $ont_fastqs/B32.fastq.gz $ont_fastqs/B33.fastq.gz $ont_fastqs/B34.fastq.gz > $ont_target/4samples/B31323334.fastq.gz
# cat $ont_fastqs/B31.fastq.gz $ont_fastqs/B32.fastq.gz $ont_fastqs/B33.fastq.gz > $ont_target/3samples/B313233.fastq.gz
# cat $ont_fastqs/B31.fastq.gz $ont_fastqs/B32.fastq.gz > $ont_target/2samples/B3132.fastq.gz

# # merge ont kidney fastqs
# cat $ont_fastqs/K31.fastq.gz $ont_fastqs/K32.fastq.gz $ont_fastqs/K33.fastq.gz $ont_fastqs/K34.fastq.gz > $ont_target/4samples/K31323334.fastq.gz
# cat $ont_fastqs/K31.fastq.gz $ont_fastqs/K32.fastq.gz $ont_fastqs/K33.fastq.gz > $ont_target/3samples/K313233.fastq.gz
# cat $ont_fastqs/K31.fastq.gz $ont_fastqs/K32.fastq.gz > $ont_target/2samples/K3132.fastq.gz


# MERGE ISOSEQ_TRIMMED
# merge isoseq brain bams
samtools merge -o $isoseq_trimmed_target/4samples/B31323334.aln.sorted.bam $isoseq_trimmed_brain/B31.aln.sorted.bam $isoseq_trimmed_brain/B32.aln.sorted.bam $isoseq_trimmed_brain/B33.aln.sorted.bam $isoseq_trimmed_brain/B34.aln.sorted.bam
samtools merge -o $isoseq_trimmed_target/3samples/B313233.aln.sorted.bam $isoseq_trimmed_brain/B31.aln.sorted.bam $isoseq_trimmed_brain/B32.aln.sorted.bam $isoseq_trimmed_brain/B33.aln.sorted.bam
samtools merge -o $isoseq_trimmed_target/2samples/B3132.aln.sorted.bam $isoseq_trimmed_brain/B31.aln.sorted.bam $isoseq_trimmed_brain/B32.aln.sorted.bam

# merge isoseq kidney bams
samtools merge -o $isoseq_trimmed_target/4samples/K31323334.aln.sorted.bam $isoseq_trimmed_kidney/K31.aln.sorted.bam $isoseq_trimmed_kidney/K32.aln.sorted.bam $isoseq_trimmed_kidney/K33.aln.sorted.bam $isoseq_trimmed_kidney/K34.aln.sorted.bam
samtools merge -o $isoseq_trimmed_target/3samples/K313233.aln.sorted.bam $isoseq_trimmed_kidney/K31.aln.sorted.bam $isoseq_trimmed_kidney/K32.aln.sorted.bam $isoseq_trimmed_kidney/K33.aln.sorted.bam
samtools merge -o $isoseq_trimmed_target/2samples/K3132.aln.sorted.bam $isoseq_trimmed_kidney/K31.aln.sorted.bam $isoseq_trimmed_kidney/K32.aln.sorted.bam

# merge isoseq brain fastqs
cat $isoseq_trimmed_brain/B31.concat.fastq $isoseq_trimmed_brain/B32.concat.fastq $isoseq_trimmed_brain/B33.concat.fastq $isoseq_trimmed_brain/B34.concat.fastq | gzip > $isoseq_trimmed_target/4samples/B31323334.concat.fastq.gz
cat $isoseq_trimmed_brain/B31.concat.fastq $isoseq_trimmed_brain/B32.concat.fastq $isoseq_trimmed_brain/B33.concat.fastq | gzip > $isoseq_trimmed_target/3samples/B313233.concat.fastq.gz
cat $isoseq_trimmed_brain/B31.concat.fastq $isoseq_trimmed_brain/B32.concat.fastq | gzip > $isoseq_trimmed_target/2samples/B3132.concat.fastq.gz

# merge isoseq kidney fastqs
cat $isoseq_trimmed_kidney/K31.concat.fastq $isoseq_trimmed_kidney/K32.concat.fastq $isoseq_trimmed_kidney/K33.concat.fastq $isoseq_trimmed_kidney/K34.concat.fastq | gzip > $isoseq_trimmed_target/4samples/K31323334.concat.fastq.gz
cat $isoseq_trimmed_kidney/K31.concat.fastq $isoseq_trimmed_kidney/K32.concat.fastq $isoseq_trimmed_kidney/K33.concat.fastq | gzip > $isoseq_trimmed_target/3samples/K313233.concat.fastq.gz
cat $isoseq_trimmed_kidney/K31.concat.fastq $isoseq_trimmed_kidney/K32.concat.fastq | gzip > $isoseq_trimmed_target/2samples/K3132.concat.fastq.gz

# MERGE ISOSEQ
# merge isoseq brain bams
samtools merge -o $isoseq_target/4samples/B31323334.aln.sorted.bam $isoseq/B31/B31.aln.sorted.bam $isoseq/B32/B32.aln.sorted.bam $isoseq/B33/B33.aln.sorted.bam $isoseq/B34/B34.aln.sorted.bam
samtools merge -o $isoseq_target/3samples/B313233.aln.sorted.bam $isoseq/B31/B31.aln.sorted.bam $isoseq/B32/B32.aln.sorted.bam $isoseq/B33/B33.aln.sorted.bam
samtools merge -o $isoseq_target/2samples/B3132.aln.sorted.bam $isoseq/B31/B31.aln.sorted.bam $isoseq/B32/B32.aln.sorted.bam

# merge isoseq kidney bams
samtools merge -o $isoseq_target/4samples/K31323334.aln.sorted $isoseq/K31/K31.aln.sorted.bam $isoseq/K32/K32.aln.sorted.bam $isoseq/K33/K33.aln.sorted.bam $isoseq/K34/K34.aln.sorted.bam
samtools merge -o $isoseq_target/3samples/K313233.aln.sorted.b $isoseq/K31/K31.aln.sorted.bam $isoseq/K32/K32.aln.sorted.bam $isoseq/K33/K33.aln.sorted.bam
samtools merge -o $isoseq_target/2samples/K3132.aln.sorted.bam $isoseq/K31/K31.aln.sorted.bam $isoseq/K32/K32.aln.sorted.bam

# merge isoseq brain fastqs
cat $isoseq/B31/B31.concat.fastq $isoseq/B32/B32.concat.fastq $isoseq/B33/B33.concat.fastq $isoseq/B34/B34.concat.fastq | gzip > $isoseq_target/4samples/B31323334.concat.fastq.gz
cat $isoseq/B31/B31.concat.fastq $isoseq/B32/B32.concat.fastq $isoseq/B33/B33.concat.fastq | gzip > $isoseq_target/3samples/B313233.concat.fastq.gz
cat $isoseq/B31/B31.concat.fastq $isoseq/B32/B32.concat.fastq | gzip > $isoseq_target/2samples/B3132.concat.fastq.gz

# merge isoseq kidney fastqs
cat $isoseq/K31/K31.concat.fastq $isoseq/K32/K32.concat.fastq $isoseq/K33/K33.concat.fastq $isoseq/K34/K34.concat.fastq | gzip > $isoseq_target/4samples/K31323334.concat.fastq.gz
cat $isoseq/K31/K31.concat.fastq $isoseq/K32/K32.concat.fastq $isoseq/K33/K33.concat.fastq | gzip > $isoseq_target/3samples/K313233.concat.fastq.gz
cat $isoseq/K31/K31.concat.fastq $isoseq/K32/K32.concat.fastq | gzip > $isoseq_target/2samples/K3132.concat.fastq.gz
