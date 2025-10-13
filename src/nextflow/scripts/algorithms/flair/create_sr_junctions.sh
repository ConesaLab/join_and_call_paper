#!/bin/bash

flair_location="/home/fabianje/tools/flair"
config_file="/storage/gge/Fabian/nih/data/metadata/flair/flair_sr_junc_config.tsv"
out_location="/storage/gge/Fabian/nih/data/metadata/flair"
joblog="${out_location}/joblog_flair_sr_junc.tsv"

nSamples=12

sbatch --array=1-$nSamples create_sr_junctions.sbatch $flair_location $config_file $out_location $joblog
