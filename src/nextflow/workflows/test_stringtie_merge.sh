#!/bin/bash

# Test StringTie merge pipeline on a single dataset:
# PacBio IsoQuant (5 brain + 5 kidney samples).
# Run this from the workflows/ directory on garnatxa.

PB_BASE="/home/fabianje/repos/documenting_NIH/fabian/data/output/isoseq"

METADATA_IND="${PB_BASE}/isoquant/run3_data/isoquant/isoquant_out_metadata_ind.tsv"

if [ ! -f "$METADATA_IND" ]; then
    echo "ERROR: metadata_ind not found at: $METADATA_IND"
    echo "Run check_inputs_stringtie_merge.sh first to diagnose."
    exit 1
fi

echo "Testing StringTie merge on PacBio IsoQuant data..."
echo "  metadata_ind: $METADATA_IND"
echo "  result_name:  isoseq/isoquant/stringtie_merge_test"
echo ""

sbatch nf_wrapper_stringtie_merge.sbatch \
    --metadata_ind "$METADATA_IND" \
    --result_name isoseq/isoquant/stringtie_merge_test
