#!/bin/bash
#SBATCH --job-name=create_gffutils_db
#SBATCH --output=logs/create_gffutils_db_%j.log
#SBATCH --error=logs/create_gffutils_db_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --qos=short
#SBATCH --time=5:00:00

if [ $# -ne 2 ]; then
    echo "Usage: sbatch create_gffutils_db.sh <input.gtf> <output.db>"
    exit 1
fi

INPUT_GTF="$1"
OUTPUT_DB="$2"

source ~/.bashrc
set -euo pipefail
conda deactivate
conda activate isoquant

echo "=== Job ${SLURM_JOB_ID} | $(hostname) | $(date) ==="
echo "Input GTF: ${INPUT_GTF}"
echo "Output DB: ${OUTPUT_DB}"

python -c "
import gffutils
print('Creating gffutils database...')
gffutils.create_db(
    '${INPUT_GTF}',
    '${OUTPUT_DB}',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True
)
print('Done')
"
