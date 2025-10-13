#!/bin/bash

# Join Chromosome Results Script
# This script combines chromosome-specific FLAIR collapse results into a single output

set -e  # Exit on any error

# Input parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --chromosomes) chromosomes="$2"; shift ;;
        --output_prefix) output_prefix="$2"; shift ;;
        --chrom_dir) chrom_dir="$2"; shift ;;
        --collapse_dir) collapse_dir="$2"; shift ;;
        --sort_script) sort_script="$2"; shift ;;
        --rename_script) rename_script="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate required parameters
if [[ -z "$chromosomes" || -z "$output_prefix" || -z "$chrom_dir" || -z "$collapse_dir" ]]; then
    echo "Usage: $0 --chromosomes <chromosome_list> --output_prefix <prefix> --chrom_dir <chromosome_directory> --collapse_dir <collapse_directory> [--sort_script <sort_script_path>] [--rename_script <rename_script_path>]"
    echo ""
    echo "Required parameters:"
    echo "  --chromosomes: Space-separated list of chromosome names"
    echo "  --output_prefix: Output file prefix"
    echo "  --chrom_dir: Directory containing chromosome-specific results"
    echo "  --collapse_dir: Directory for final combined results"
    echo ""
    echo "Optional parameters:"
    echo "  --sort_script: Path to GTF sorting script"
    echo "  --rename_script: Path to GTF renaming script"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$collapse_dir"

echo "=== Joining Chromosome Results ==="
echo "Chromosomes: $chromosomes"
echo "Output prefix: $output_prefix"
echo "Chromosome directory: $chrom_dir"
echo "Collapse directory: $collapse_dir"
echo ""

# Combine all chromosome GTF files
combined_gtf="${collapse_dir}/${output_prefix}_combined.isoforms.gtf"
> "$combined_gtf"  # Create empty file

echo "Combining chromosome GTF files..."

for chrom in $chromosomes; do
    chrom_gtf="${chrom_dir}/${output_prefix}_${chrom}.isoforms.gtf"
    if [[ -f "$chrom_gtf" ]]; then
        echo "Adding chromosome $chrom results..."
        cat "$chrom_gtf" >> "$combined_gtf"
    else
        echo "Warning: No results found for chromosome $chrom at $chrom_gtf"
    fi
done

echo "Combined GTF file created: $combined_gtf"
echo ""

# Sort the combined GTF file
sorted_gtf="${collapse_dir}/${output_prefix}.isoforms.gtf"
echo "Sorting combined GTF file..."

if [[ -n "$sort_script" && -f "$sort_script" ]]; then
    echo "Using custom sort script: $sort_script"
    "$sort_script" -i "$combined_gtf" -o "$sorted_gtf"
else
    echo "Using simple sort (chromosome, then position)"
    sort -k1,1 -k4,4n "$combined_gtf" > "$sorted_gtf"
fi

echo "Sorted GTF file: $sorted_gtf"
echo ""

# Rename GTF file for consistency
final_gtf="${collapse_dir}/${output_prefix}.isoforms_renamed.gtf"
echo "Renaming GTF file..."

if [[ -n "$rename_script" && -f "$rename_script" ]]; then
    echo "Using custom rename script: $rename_script"
    "$rename_script" "$sorted_gtf" > "$final_gtf"
else
    echo "Using simple copy (no renaming)"
    cp "$sorted_gtf" "$final_gtf"
fi

echo "Final GTF file: $final_gtf"
echo ""

# Combine FASTA files
combined_fasta="${collapse_dir}/${output_prefix}.isoforms.fa"
> "$combined_fasta"  # Create empty file

echo "Combining chromosome FASTA files..."

for chrom in $chromosomes; do
    chrom_fasta="${chrom_dir}/${output_prefix}_${chrom}.isoforms.fa"
    if [[ -f "$chrom_fasta" ]]; then
        echo "Adding chromosome $chrom FASTA..."
        cat "$chrom_fasta" >> "$combined_fasta"
    else
        echo "Warning: No FASTA results found for chromosome $chrom at $chrom_fasta"
    fi
done

echo "Combined FASTA file: $combined_fasta"
echo ""

# Combine BED files
combined_bed="${collapse_dir}/${output_prefix}.isoforms.bed"
> "$combined_bed"  # Create empty file

echo "Combining chromosome BED files..."

for chrom in $chromosomes; do
    chrom_bed="${chrom_dir}/${output_prefix}_${chrom}.isoforms.bed"
    if [[ -f "$chrom_bed" ]]; then
        echo "Adding chromosome $chrom BED..."
        cat "$chrom_bed" >> "$combined_bed"
    else
        echo "Warning: No BED results found for chromosome $chrom at $chrom_bed"
    fi
done

echo "Combined BED file: $combined_bed"
echo ""

echo "=== Chromosome Results Joined Successfully ==="
echo "Final outputs:"
echo "  GTF file: $final_gtf"
echo "  FASTA file: $combined_fasta"
echo "  BED file: $combined_bed"
echo "" 