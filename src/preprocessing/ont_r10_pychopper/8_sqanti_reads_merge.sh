#!/bin/bash
#SBATCH --job-name=sqanti_reads_merge_ont_r10_pychopper
#SBATCH --output=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_reads_merge_%j.log
#SBATCH --error=/storage/gge/Fabian/ont_r10_sy5y/analysis/logs_pychopper/sqanti_reads_merge_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb
#SBATCH --qos=medium
#SBATCH --time=7-00:00:00

source ~/.bashrc

module load samtools
conda deactivate
conda activate SQANTI3.env

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "${_script_dir}/config.sh"

sqanti_tool="$HOME/tools/SQANTI3_dev"
design_csv="${_script_dir}/design_ont_r10_pychopper.csv"

echo "=== Job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="

python3 "${sqanti_tool}/sqanti_reads.py" \
  --genome "${assembly}" \
  --annotation "${ref_annotation}" \
  -de "${design_csv}" \
  -i "${sqanti_reads_dir}" \
  -f "condition" \
  -p "ONT_R10_pychopper" \
  -d "${sqanti_reads_dir}/out" \
  --force_id_ignore \
  -t 12 -n 1
