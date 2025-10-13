#!/bin/bash

gtf_file="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"

awk -F '\t' '
  # Only process lines where the first column starts with "SIRV"
  $1 ~ /^SIRV/ {
    # Use a regular expression to extract the transcript_id value
    match($9, /transcript_id "([^"]+)"/, arr)
    if (arr[1] != "") {
      # Print the extracted transcript_id
      print arr[1]
    }
  }
' "$gtf_file" | sort -u
