#!/bin/bash
# reformat annotation for TALON

# module load anaconda
# source activate talon

# example: ./reformat_annotation.sh -a /storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf -o ~/repos/documenting_NIH/fabian/data

# Default values
delete_tmp=false

# Parse command line options
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -a|--annotation)
      annotation="$2"
      shift 2
      ;;
    -o|--annotation_talon)
      annotation_talon="$2"
      shift 2
      ;;
    -d|--delete_tmp)
      delete_tmp=true
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Check if mandatory options are provided
if [ -z "$annotation" ] || [ -z "$annotation_talon" ]; then
  echo "Missing mandatory option: -a or --annotation (reference annotation as .gtf) and/or --annotation_talon (output file for TALON formatted annotation)"
  exit 1
fi

tmp=$(dirname "${annotation_talon}")/tmp

# create directories if they do not exist
if [ ! -d "$tmp" ]; then
    echo "Creating directory: $tmp"
    mkdir -p "$tmp"
fi

annotation_basename=$(basename "${annotation%.gtf}")
annotation_exons="$tmp/${annotation_basename}_exons.gtf"
annotation_fixedsirv="${tmp}/$(basename "${annotation_exons%.gtf}_fixedsirv.gtf")"
annotation_reformatted="${annotation_fixedsirv%.gtf}_reformatted.gtf"
formatted_annotation="${tmp}/$(basename "${annotation_reformatted%.gtf}_formatted.gtf")"
annotation_deduplicated="${tmp}/$(basename "${annotation_reformatted%.gtf}_deduplicated.gtf")"
annotation_deduplicated2="${tmp}/$(basename "${annotation_reformatted%.gtf}_deduplicated2.gtf")"
annotation_expanded="${tmp}/$(basename "${annotation_reformatted%.gtf}_expanded.gtf")"
final_gtf="$annotation_talon"

if [ ! -f "$final_gtf" ]; then
  # remove all gene and transcript entries; keep only exon and other entries
  echo "removing all 'gene' and 'transcript' entries."
  awk '$3 != "gene" && $3 != "transcript"' $annotation > $annotation_exons

  echo "Processing SIRV entries for strand-specific gene separation"
  awk 'BEGIN {FS=OFS="\t"} 
  /^#/{next} 
  $2=="LexogenSIRVData" {
    gene_id=$9; 
    match(gene_id, /gene_id "([^"]+)"/, arr); 
    gene_id=arr[1]; 
    if ($7=="+") {
      plus_genes[gene_id]=1; 
    } else if ($7=="-") {
      minus_genes[gene_id]=1;
    }
  } 
  END {
    for (gene in plus_genes) {
      if (gene in minus_genes) {
        print gene, gene "_plus", gene "_minus";
      }
    }
  }' $annotation_exons > "${tmp}/genes_to_process.txt"

  # Create a copy of the exons file to work on
  cp $annotation_exons "${annotation_exons}.tmp"

  # Process each gene to adjust its gene_id based on strand
  while read gene plus minus; do
    awk -v gene="$gene" -v plus="$plus" -v minus="$minus" 'BEGIN {FS=OFS="\t"} 
    {
      if ($9 ~ "gene_id \"" gene "\"") {
        if ($7=="+") {
          gsub("gene_id \"" gene "\"", "gene_id \"" plus "\"", $9);
        } else if ($7=="-") {
          gsub("gene_id \"" gene "\"", "gene_id \"" minus "\"", $9);
        }
      }
    } {print}' "${annotation_exons}.tmp" > "${annotation_exons}.tmp2"
    # Use the updated file for the next iteration
    mv "${annotation_exons}.tmp2" "${annotation_exons}.tmp"
  done < "${tmp}/genes_to_process.txt"

  # Move the final updated file to the desired location
  mv "${annotation_exons}.tmp" $annotation_fixedsirv

  # run talon_reformat_gtf to create transcript and gene entries
  echo "running talon_reformat_gtf"
  talon_reformat_gtf -gtf $annotation_fixedsirv


  # fix formatting (the gene_id sometimes contains the suffix string "gene_name"; following the gene_id, the space after the semicolon before the gene_name is sometimes missing)
  awk '{ gsub("gene_id \"[^\"]*gene_name", "gene_id \""); gsub(";gene_name", "; gene_name"); print }' "$annotation_reformatted" > "$formatted_annotation"


  # determine all duplicate genes
  dup_genes="${tmp}/duplicate_genes.txt"
  echo "determining duplicate gene_id entries on multiple chromosomes"
  awk 'BEGIN {OFS="\t"} $3 == "gene" { gene_name = ""; match($0, /gene_name "[^"]+"/); if (RSTART) gene_name = substr($0, RSTART + 10, RLENGTH - 11); if (gene_name != "") genes[gene_name]++; } END { for (gene in genes) { if (genes[gene] > 1) { gsub(/^"/, "", gene); print gene; } } }' $formatted_annotation > $dup_genes

  tmp_in="${formatted_annotation%.gtf}_in.gtf"
  tmp_out="${formatted_annotation%.gtf}_out.gtf"
  cp "$formatted_annotation" "$tmp_in"

  # rename duplicate genes
  while read -r gene_id; do
    echo "renaming entries with duplicate gene_id ${gene_id}"
    awk -v gene_id="$gene_id" 'BEGIN {OFS="\t"} { chr = $1; gsub(gene_id, gene_id "-" chr, $0); } { print }' $tmp_in > $tmp_out
    cp $tmp_out $tmp_in
  done < $dup_genes

  cp $tmp_out $annotation_deduplicated

  # for gene entries, copy gene_name content into gene_id as well
  awk -F'\t' '
    $3 == "gene" {
      if (match($9, /gene_name "([^"]*)"/, gene_name_match)) {
        gene_name_content = gene_name_match[1];
        gsub(/gene_id "[^"]*"/, "gene_id \"" gene_name_content "\"", $9);
      }
    }
    { print }
  ' OFS='\t' "$annotation_deduplicated" > "$annotation_deduplicated2"
  
  # add gene_name and exon_number entries to SIRV data
  awk -F'\t' 'BEGIN {exon_number=0; last_gene_id=""} 
  {
    if($2=="LexogenSIRVData") {
      split($9, a, "; ?");
      for(i in a) {
        if(a[i] ~ /^gene_id/) {
          gene_id=gensub(/.*"([^";]+)"/, "\\1", "g", a[i])
          $9=$9" gene_name \""gene_id"\";"
          # if($3=="exon") {
          #   if(gene_id!=last_gene_id) {
          #     exon_number=1
          #     last_gene_id=gene_id
          #   } else {
          #     exon_number++
          #   }
          #   $9=$9" exon_number \""exon_number"\";"
          # }
        }
      }
    }
    print $0
  }' OFS='\t' "$annotation_deduplicated2" > "$annotation_expanded"

  cp "$annotation_expanded" "$final_gtf"

  echo "finished reformatting gtf, result in ${final_gtf}"

  # remove temp directory
  if [ -d "$tmp" ] && [ "$delete_tmp" == true ]; then
      echo "Deleting the tmp directory and all its contents: $tmp"
      rm -rf "$tmp"
  fi

fi