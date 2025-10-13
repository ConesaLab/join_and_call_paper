#!/bin/bash
# Count the number of lines in samples.tsv and subtract 1 for the header
num_samples=$(($(wc -l < samples.tsv) - 1))

# Submit the sbatch array job, specifying the range from 1 to the number of samples
sbatch --wait --array=1-$num_samples unzip_fastq.sbatch

