#!/bin/bash
#SBATCH --job-name=sq3
#SBATCH --qos short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=20gb # Per processor memory
#SBATCH -t 1-00:00:00     # Walltime
#SBATCH -o log/sq3_%A_%a.out

module load anaconda
source activate SQANTI3.env

export PYTHONPATH=$PYTHONPATH:$HOME/tools/cDNA_Cupcake
export PYTHONPATH=$PYTHONPATH:$HOME/tools/cDNA_Cupcake/sequence/

# Input
sq3_input=$1 # Working directory
genome=$2 # Reference genome
annotation=$3 # Reference annotation
sqanti_location=$4
cage=$5
polyA_bed=$6
n_cores="$SLURM_CPUS_PER_TASK"
task_number="$SLURM_ARRAY_TASK_ID"

polyA="${sqanti_location}/data/polyA_motifs/mouse_and_human.polyA_motif.txt"

# read corresponding line of config
read -r gtf count sr <<< $(sed -n "${task_number}p" $sq3_input)
output_prefix=$(basename "$gtf" .gtf)
output_location=$(dirname "$gtf")/"${output_prefix}_sq3_orth"

# Get the corresponding sr sample
sj=${sr}/*SJ.out.tab
sj=$(echo $sj | sed 's/ /,/g')
bam=${sr}/*.sortedByCoord.out.bam


# Write bam dir to file
echo $bam | sed 's/ /\n/g' > tmp_$SLURM_ARRAY_TASK_ID.fofn


python $sqanti_location/sqanti3_qc.py \
    ${gtf} $annotation $genome \
    -fl ${count} \
    -o ${output_prefix} \
    -d ${output_location} \
    --skipORF \
    -t $n_cores \
    -c ${sj} \
    --CAGE_peak ${cage} \
    --polyA_peak ${polyA_bed} \
    --polyA_motif_list $polyA \
    --report "skip" \
    --SR_bam tmp_$SLURM_ARRAY_TASK_ID.fofn 
