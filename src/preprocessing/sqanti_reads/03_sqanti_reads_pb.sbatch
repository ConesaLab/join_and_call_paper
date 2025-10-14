#!/bin/bash
#SBATCH --job-name=sqanti_reads_pb   	# Job name (showed with squeue)
#SBATCH --output=sqanti_reads_pb_%j.out  	# Standard output and error log
#SBATCH --ntasks=1			# Required only 1 task 
#SBATCH --cpus-per-task=1               # Required only 1 cpu
#SBATCH --mem=20G			# Required 10GB of memory
#SBATCH --time=08:00:00			# Required 5 minutes of execution.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

sqanti_dir="/home/apadepe/lr_pipelines/SQANTI3"
ref_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
assembly="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
sqanti_out_dir="/storage/gge/Alejandro/nih/iso_seq/SQANTI_reads"


python3 ${sqanti_dir}/sqanti_reads.py --genome $assembly \
    --annotation ${ref_annotation} -de design_pb.csv -i $sqanti_out_dir \
    -f "condition" -p "PacBio" -d ${sqanti_out_dir}/out --skip_hash