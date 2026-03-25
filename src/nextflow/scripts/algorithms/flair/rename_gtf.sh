#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.gtf>"
    exit 1
fi

input_gtf="$1"

# Process the GTF file
awk -F'\t' '
BEGIN { OFS="\t" }
{    
    if ($9 ~ /transcript_id/) {
        
        # Extract gene_id and transcript_id
        match($9, /gene_id "([^"]+)"/, gene)
        match($9, /transcript_id "([^"]+)"/, trans)
                
        # Create new attribute string
        new_attr = $9
        sub(/transcript_id "[^"]+"/, "transcript_id \"" trans[1] "_" gene[1] "\"", new_attr)
        
        $9 = new_attr
    }
    print
}' "$input_gtf"
