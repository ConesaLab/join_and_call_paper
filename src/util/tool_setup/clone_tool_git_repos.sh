#!/bin/bash

DEFAULT_TOOLS_DIR="$HOME/tools"

tools_dir="${1:-$DEFAULT_TOOLS_DIR}"

original_wd=$(pwd)

mkdir -p $tools_dir

git clone https://github.com/ConesaLab/SQANTI3.git "$tools_dir/SQANTI3_dev"
git -C "$tools_dir/SQANTI3_dev" checkout 5f182ae54b3decf65f18bd8099cdd22e0914ad6d # commit of 5.3.6 release

git clone https://github.com/GenomeRIK/tama.git "$tools_dir/tama"
git -C "$tools_dir/tama" checkout 2fa3c308282190c413e9bf0e0b49e63086eef7d4 # latest commit as of 2025-02-05 (commit from 2023-04-12)

# for TranscriptClean (needed by TALON)
# git clone https://github.com/mortazavilab/TranscriptClean "$tools_dir/TranscriptClean"
# git -C "$tools_dir/TranscriptClean" checkout 56fe8f25fe3ededc1f8f048b7cc00139a55ef15c # latest commit as of 2025-02-05 (commit from 2024-01-31)

## TranscriptClean setup
# conda activate transcriptclean
# cd $tools_dir/TranscriptClean
# pip install -e .
# cd $original_wd
# conda deactivate

# for Mandalorion
# git clone https://github.com/christopher-vollmers/Mandalorion "$tools_dir/Mandalorion"
# git -C "$tools_dir/Mandalorion" checkout 6b8ae0558cc7d52014d7df32af1f1357337157aa # latest commit as of 2025-02-05 (commit from 2024-02-19)

## Mandalorion setup
# cd $tools_dir/Mandalorion
# git clone https://github.com/lh3/minimap2 "$tools_dir/Mandalorion/minimap2"
# git -C "$tools_dir/Mandalorion/minimap2" checkout 8170693de39b667d11c8931d343c94a23a7690d2 # latest release (2.28-r1209) as of 2025-02-05
# cd minimap2 && make
# cd $tools_dir/Mandalorion

# wget https://github.com/yangao07/abPOA/releases/download/v1.4.1/abPOA-v1.4.1.tar.gz
# tar -zxvf abPOA-v1.4.1.tar.gz
# cd abPOA-v1.4.1; make
# cd $tools_dir/Mandalorion
# rm abPOA-v1.4.1.tar.gz
# cd $original_wd
