#!/usr/bin/env bash

DEFAULT_DIR="."

dir="${1:-$DEFAULT_DIR}"

# essentials
conda env create -f $dir/SQANTI3.env.yaml
conda env create -f $dir/tama.yaml
conda env create -f $dir/tusco_env.yaml

# isoform identification tools
conda env create -f $dir/isoquant.yaml
# conda env create -f $dir/flair_conda_env.yaml
# conda env create -f $dir/talon.yaml && conda env create -f $dir/transcriptclean.yaml
# conda env create -f $dir/mandalorion.yaml
# conda env create -f $dir/bambu.yaml && conda activate bambu_env && Rscript bambu_setup.R
