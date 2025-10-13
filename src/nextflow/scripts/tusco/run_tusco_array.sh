#!/bin/bash

# Check if config file exists
if [ ! -f "tusco_config.csv" ]; then
    echo "Error: tusco_config.csv not found"
    exit 1
fi

# Calculate number of lines (minus header) in config file
array_size=$(($(wc -l < tusco_config.csv) - 1))

if [ $array_size -le 0 ]; then
    echo "Error: No data lines found in tusco_config.csv"
    exit 1
fi

echo "Found $array_size samples to process"

# Submit the array job with the calculated array size
job_id=$(sbatch --array=1-$array_size execute_tusco.sbatch | grep -o '[0-9]\+')

echo "Submitted array job $job_id with $array_size tasks" 