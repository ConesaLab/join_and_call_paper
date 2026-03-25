#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --qos long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH -t 15-00:00:00     # Walltime
#SBATCH -o log/old_nextflow_%A.out

# ONT
echo "sbatch --wait nextflow_wrapper.sbatch --data ont_old --algorithm flair --use_sr true --sr_config star_old --stringent true --result_name ont/flair/old_mice/full"
sbatch --wait nextflow_wrapper.sbatch --data ont_old --algorithm flair --use_sr true --sr_config star_old --stringent true --result_name ont/flair/old_mice/full
echo "sbatch --wait nextflow_wrapper.sbatch --data ont_old_p3 --algorithm flair --use_sr true --sr_config star_old_p3 --stringent true --result_name ont/flair/old_mice/partial_join_3sample"
sbatch --wait nextflow_wrapper.sbatch --data ont_old_p3 --algorithm flair --use_sr true --sr_config star_old_p3 --stringent true --result_name ont/flair/old_mice/partial_join_3sample
echo "sbatch --wait nextflow_wrapper.sbatch --data ont_old_p2 --algorithm flair --use_sr true --sr_config star_old_p2 --stringent true --result_name ont/flair/old_mice/partial_join_2sample"
sbatch --wait nextflow_wrapper.sbatch --data ont_old_p2 --algorithm flair --use_sr true --sr_config star_old_p2 --stringent true --result_name ont/flair/old_mice/partial_join_2sample

# IsoSeq
echo "sbatch --wait nextflow_wrapper.sbatch --data isoseq_old --algorithm flair --use_sr true --sr_config star_old --stringent true --result_name isoseq/flair/old_mice/full"
sbatch --wait nextflow_wrapper.sbatch --data isoseq_old --algorithm flair --use_sr true --sr_config star_old --stringent true --result_name isoseq/flair/old_mice/full
echo "sbatch --wait nextflow_wrapper.sbatch --data isoseq_old_p3 --algorithm flair --use_sr true --sr_config star_old_p3 --stringent true --result_name isoseq/flair/old_mice/partial_join_3sample"
sbatch --wait nextflow_wrapper.sbatch --data isoseq_old_p3 --algorithm flair --use_sr true --sr_config star_old_p3 --stringent true --result_name isoseq/flair/old_mice/partial_join_3sample
echo "sbatch --wait nextflow_wrapper.sbatch --data isoseq_old_p2 --algorithm flair --use_sr true --sr_config star_old_p2 --stringent true --result_name isoseq/flair/old_mice/partial_join_2sample"
sbatch --wait nextflow_wrapper.sbatch --data isoseq_old_p2 --algorithm flair --use_sr true --sr_config star_old_p2 --stringent true --result_name isoseq/flair/old_mice/partial_join_2sample
