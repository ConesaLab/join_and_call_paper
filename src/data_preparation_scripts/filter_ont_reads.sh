#!/bin/bash
# Count the number of .fastq files in the merged directory
num_fastq_files=$(find /home/fabianje/repos/documenting_NIH/fabian/data/ont/fastq/merged -type f -name "*.fastq" | wc -l)

# Submit the sbatch array job, specifying the range from 0 to one less than the number of .fastq files
sbatch --array=1-$num_fastq_files filter_ont_reads.sbatch
