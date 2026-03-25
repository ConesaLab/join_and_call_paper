#!/bin/bash

ont_samples="/storage/gge/Fabian/nih/data/metadata/ont_samples.tsv"
ont_subset_target="/storage/gge/Fabian/nih/data/ont/subset"

sbatch create_data_subset.sbatch $ont_samples $ont_subset_target
