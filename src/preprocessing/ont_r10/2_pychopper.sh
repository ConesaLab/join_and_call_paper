#!/bin/bash
#SBATCH --job-name=pychopper_ont_r10
#SBATCH --output=../analysis/logs/pychopper_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=15gb
#SBATCH --qos=short
#SBATCH --time=24:00:00
#SBATCH --array=0-3

source ~/.bashrc

module load samtools
conda deactivate
conda activate pychopper

readarray myarray < list_fastqs.fof

file=${myarray[$SLURM_ARRAY_TASK_ID]}

echo $file

file_dir=$(dirname $file)
ori_name=$(basename $file)
filename="${ori_name%.fastq}"

echo $filename

pychopper_dir="/storage/gge/Fabian/ont_r10_sy5y/fastq/pychopper"
mkdir -p ${pychopper_dir}

pychopper -k LSK114 \
    -r ${pychopper_dir}/${filename}_pychopper_report.pdf \
    -u ${pychopper_dir}/${filename}_unclassified.fastq \
    -w ${pychopper_dir}/${filename}_rescued.fastq \
    -t 4 \
    ${file_dir}/${filename}.fastq \
    ${pychopper_dir}/${filename}_full_length.fastq

# Merge full_length + rescued into single file for mapping
cat ${pychopper_dir}/${filename}_full_length.fastq \
    ${pychopper_dir}/${filename}_rescued.fastq \
    > ${pychopper_dir}/${filename}_fl_rescued.fastq

echo "Done: ${filename}"
echo "  Full-length: ${pychopper_dir}/${filename}_full_length.fastq"
echo "  Rescued: ${pychopper_dir}/${filename}_rescued.fastq"
echo "  Merged (for mapping): ${pychopper_dir}/${filename}_fl_rescued.fastq"
