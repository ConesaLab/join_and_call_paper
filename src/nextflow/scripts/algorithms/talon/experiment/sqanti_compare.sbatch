#!/bin/bash
#SBATCH --job-name=sq3
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 01:00:00     # Walltime
#SBATCH -o log/sq3_%A_%a.out

# sbatch ./sqanti_compare.sbatch /home/fabianje/repos/documenting_NIH/fabian/data/output/ont_subset/talon_experiment/isoform_detection/talon/B100K0/indiv/indiv_talon.sorted.gtf /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa /home/fabianje/repos/documenting_NIH/fabian/data/output/ont_subset/talon_experiment/isoform_detection/talon/B100K0/Bconcat/Bconcat_talon.sorted.gtf /home/fabianje/tools/SQANTI3_dev
# sbatch ./sqanti_compare.sbatch /home/fabianje/repos/documenting_NIH/fabian/data/output/ont_subset/talon_experiment/isoform_detection/talon/B100K0/Bconcat/Bconcat_talon.sorted.gtf /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa /home/fabianje/repos/documenting_NIH/fabian/data/output/ont_subset/talon_experiment/isoform_detection/talon/B100K0/indiv/indiv_talon.sorted.gtf /home/fabianje/tools/SQANTI3_dev

module load anaconda
source activate SQANTI3.env

export PYTHONPATH=$PYTHONPATH:$HOME/tools/cDNA_Cupcake
export PYTHONPATH=$PYTHONPATH:$HOME/tools/cDNA_Cupcake/sequence/

# Input
gtf=$1 # Working directory
genome=$2 # Reference genome
annotation=$3 # Reference annotation
sqanti_location=$4
n_cores="$SLURM_CPUS_PER_TASK"

# read corresponding line of config
output_prefix=$(basename "$gtf" .gtf)
output_location=$(dirname "$gtf")/"${output_prefix}_sq3"

python $sqanti_location/sqanti3_qc.py \
    $gtf $annotation $genome \
    -o ${output_prefix} \
    -d ${output_location} \
    --skipORF \
    -t $n_cores