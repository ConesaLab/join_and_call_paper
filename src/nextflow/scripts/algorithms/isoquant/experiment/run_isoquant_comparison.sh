#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath $0)")

joblog="${SCRIPT_DIR}/output/joblog_subset.tsv"
> $joblog

# data
genome="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
annotation_db="/storage/gge/Fabian/nih/data/metadata/isoquant/mm39.ncbiRefSeq_SIRV.db"

# ONT
# bam_file="/storage/gge/nih/Nanopore/data/bams/B31_primary_aln_sorted.bam"
bam_file="/storage/gge/Fabian/nih/data/ont_pilot/bam/length_filtered/subset/B3_1_fl_r_300bp_subset_primary_aln_sorted.bam"

jobid=$(sbatch $SCRIPT_DIR/execute_isoquant.sbatch $SCRIPT_DIR $bam_file B31_ONT_SUBSET_NEW $genome $annotation_db $SCRIPT_DIR/output_subset nanopore false ONT_SUBSET_NEW | awk '{print $NF}')
echo -e "IQ_ONT_SUBSET_NEW\t${jobid}" >> $joblog
jobid=$(sbatch $SCRIPT_DIR/execute_isoquant_old.sbatch $SCRIPT_DIR $bam_file B31_ONT_SUBSET_OLD $genome $annotation_db $SCRIPT_DIR/output_subset nanopore false ONT_SUBSET_OLD | awk '{print $NF}')
echo -e "IQ_ONT_SUBSET_OLD\t${jobid}" >> $joblog

# PacBio
# bam_file="/storage/gge/nih/PacBio_IsoSeq/merged_reads/mouse_B100K0/B31.aln.sorted.bam"

# jobid=$(sbatch $SCRIPT_DIR/execute_isoquant.sbatch $SCRIPT_DIR $bam_file B31_PB_NEW $genome $annotation_db $SCRIPT_DIR/output pacbio_ccs true PB_NEW | awk '{print $NF}')
# echo -e "IQ_PB_NEW\t${jobid}" >> $joblog
# jobid=$(sbatch $SCRIPT_DIR/execute_isoquant_old.sbatch $SCRIPT_DIR $bam_file B31_PB_OLD $genome $annotation_db $SCRIPT_DIR/output pacbio_ccs true PB_OLD | awk '{print $NF}')
# echo -e "IQ_PB_OLD\t${jobid}" >> $joblog
