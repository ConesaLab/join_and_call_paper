#!/bin/bash

source="/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped"
target="/storage/gge/Fabian/nih/data/star_sj_subset"

sbatch create_star_sj_subset.sbatch $source $target