#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH --qos long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH -t 15-00:00:00     # Walltime
#SBATCH -o log/pj_nextflow_%A.out

# ONT
echo "sbatch --wait nextflow_wrapper.sbatch --data ont_p4 --algorithm flair --use_sr true --sr_config star_p4 --stringent true --result_name ont/flair/partial_join_4sample"
sbatch --wait nextflow_wrapper.sbatch --data ont_p4 --algorithm flair --use_sr true --sr_config star_p4 --stringent true --result_name ont/flair/partial_join_4sample
echo "sbatch --wait nextflow_wrapper.sbatch --data ont_p3 --algorithm flair --use_sr true --sr_config star_p3 --stringent true --result_name ont/flair/partial_join_3sample"
sbatch --wait nextflow_wrapper.sbatch --data ont_p3 --algorithm flair --use_sr true --sr_config star_p3 --stringent true --result_name ont/flair/partial_join_3sample
echo "sbatch --wait nextflow_wrapper.sbatch --data ont_p2 --algorithm flair --use_sr true --sr_config star_p2 --stringent true --result_name ont/flair/partial_join_2sample"
sbatch --wait nextflow_wrapper.sbatch --data ont_p2 --algorithm flair --use_sr true --sr_config star_p2 --stringent true --result_name ont/flair/partial_join_2sample

# IsoSeq
echo "sbatch --wait nextflow_wrapper.sbatch --data isoseq_p4 --algorithm flair --use_sr true --sr_config star_p4 --stringent true --result_name isoseq/flair/partial_join_4sample"
sbatch --wait nextflow_wrapper.sbatch --data isoseq_p4 --algorithm flair --use_sr true --sr_config star_p4 --stringent true --result_name isoseq/flair/partial_join_4sample
echo "sbatch --wait nextflow_wrapper.sbatch --data isoseq_p3 --algorithm flair --use_sr true --sr_config star_p3 --stringent true --result_name isoseq/flair/partial_join_3sample"
sbatch --wait nextflow_wrapper.sbatch --data isoseq_p3 --algorithm flair --use_sr true --sr_config star_p3 --stringent true --result_name isoseq/flair/partial_join_3sample
echo "sbatch --wait nextflow_wrapper.sbatch --data isoseq_p2 --algorithm flair --use_sr true --sr_config star_p2 --stringent true --result_name isoseq/flair/partial_join_2sample"
sbatch --wait nextflow_wrapper.sbatch --data isoseq_p2 --algorithm flair --use_sr true --sr_config star_p2 --stringent true --result_name isoseq/flair/partial_join_2sample
