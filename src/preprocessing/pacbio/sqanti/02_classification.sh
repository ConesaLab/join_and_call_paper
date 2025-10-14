#!/bin/bash
#SBATCH --job-name=classification_iso-seq
#SBATCH --time 3:00:00
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --output=logs/classification_iso-seq_%A_%a.out
#SBATCH --error=logs/classification_iso-seq_%A_%a.err

# Script to generate the concatenated classification files and genepred files
# from the different iso-seq runs


# Pool1
pool="Pool1"
bc_file="bc_sample_P1.tsv"

# Isoseq base dir
isoseq_dir="/storage/gge/nih/PacBio_IsoSeq"
# To  make classification files
out_dir="/storage/gge/Alejandro/nih/iso_seq/combined_class"
out_dir_junc="/storage/gge/Alejandro/nih/iso_seq/combined_junc"
out_dir_gp="/storage/gge/Alejandro/nih/iso_seq/combined_anotation"
mkdir $out_dir
mkdir $out_dir_gp
mkdir $out_dir_junc

while read -r line
do
    readarray -t -d $'\t' line_list <<< $line
    bc=$(basename ${line_list[1]}) # remove thr trailing "\n"
    sample=$(basename ${line_list[0]})
    # get the path to the classification files
    ls ${isoseq_dir}/$pool/RUN*/Cell*/demultiplexed_bam/flnc/run_SQANTI/*${bc}*/*classification.txt > ${sample}_classification.fofn
    ls ${isoseq_dir}/$pool/RUN*/Cell*/demultiplexed_bam/flnc/run_SQANTI/*${bc}*/*corrected.genePred > ${sample}_genepred.fofn
    ls ${isoseq_dir}/$pool/RUN*/Cell*/demultiplexed_bam/flnc/run_SQANTI/*${bc}*/*junctions.txt > ${sample}_junctions.fofn

    # initiate the classification file with the header
    cp classification_header.txt ${out_dir}/${sample}_classification.txt

    tail -n +2 $(cat ${sample}_classification.fofn) >> ${out_dir}/${sample}_classification.txt

    # Delete space lines or empty
    sed -i '/classification.txt/d' ${out_dir}/${sample}_classification.txt
    sed -i '/^[[:space:]]*$/d' ${out_dir}/${sample}_classification.txt

    # write the genepred files
    cat $(cat ${sample}_genepred.fofn) >> ${out_dir_gp}/${sample}.genepred
    sort -k2,2 -k4,4n -k5,5n ${out_dir_gp}/${sample}.genepred > ${out_dir_gp}/${sample}_sorted.genepred
    ./genePredToGtf file ${out_dir_gp}/${sample}_sorted.genepred ${out_dir_gp}/${sample}.gtf

    cp $(head -n 1 ${sample}_junctions.fofn) ${out_dir}/${sample}_junctions.txt
    tail -n +2 ${sample}_junctions.fofn | xargs -I {} tail -n +2 {} >> ${out_dir}/${sample}_junctions.txt
    

done < $bc_file


# Pool2
pool="Pool2"
bc_file="bc_sample_P2.tsv"

while read -r line
do
    readarray -t -d $'\t' line_list <<< $line
    bc=$(basename ${line_list[1]}) # remove thr trailing "\n"
    sample=$(basename ${line_list[0]})
    # get the path to the classification files
    ls ${isoseq_dir}/$pool/RUN*/Cell*/demultiplexed_bam/flnc/run_SQANTI/*${bc}*/*classification.txt > ${sample}_classification.fofn
    ls ${isoseq_dir}/$pool/RUN*/Cell*/demultiplexed_bam/flnc/run_SQANTI/*${bc}*/*corrected.genePred > ${sample}_genepred.fofn
    ls ${isoseq_dir}/$pool/RUN*/Cell*/demultiplexed_bam/flnc/run_SQANTI/*${bc}*/*junctions.txt > ${sample}_junctions.fofn

    # initiate the classification file with the header
    cp classification_header.txt ${out_dir}/${sample}_classification.txt

    tail -n +2 $(cat ${sample}_classification.fofn) >> ${out_dir}/${sample}_classification.txt

    # Delete space lines or empty
    sed -i '/classification.txt/d' ${out_dir}/${sample}_classification.txt
    sed -i '/^[[:space:]]*$/d' ${out_dir}/${sample}_classification.txt

    # write the genepred files
    cat $(cat ${sample}_genepred.fofn) >> ${out_dir_gp}/${sample}.genepred
    sort -k2,2 -k4,4n -k5,5n ${out_dir_gp}/${sample}.genepred > ${out_dir_gp}/${sample}_sorted.genepred
    ./genePredToGtf file ${out_dir_gp}/${sample}_sorted.genepred ${out_dir_gp}/${sample}.gtf

    cp $(head -n 1 ${sample}_junctions.fofn) ${out_dir}/${sample}_junctions.txt
    tail -n +2 ${sample}_junctions.fofn | xargs -I {} tail -n +2 {} >> ${out_dir}/${sample}_junctions.txt
    
done < $bc_file