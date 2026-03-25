#!/bin/bash

# FLAIR Collapse for Single Chromosome
# This script runs FLAIR collapse on a single chromosome
# Can be run in parallel for multiple chromosomes

set -e  # Exit on any error

# Input parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --chrom_bed) chrom_bed="$2"; shift ;;
        --chrom_gtf) chrom_gtf="$2"; shift ;;
        --genome) genome="$2"; shift ;;
        --reads) reads="$2"; shift ;;
        --output_prefix) output_prefix="$2"; shift ;;
        --chrom) chrom="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --output_dir) output_dir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate required parameters
if [[ -z "$chrom_bed" || -z "$chrom_gtf" || -z "$genome" || -z "$reads" || -z "$output_prefix" || -z "$chrom" ]]; then
    echo "Usage: $0 --chrom_bed <chromosome.bed> --chrom_gtf <chromosome.gtf> --genome <genome.fa> --reads <reads.fastq> --output_prefix <prefix> --chrom <chromosome_name> [--threads <num_threads>] [--output_dir <output_directory>]"
    echo ""
    echo "Required parameters:"
    echo "  --chrom_bed: BED file for this chromosome"
    echo "  --chrom_gtf: GTF annotation file for this chromosome"
    echo "  --genome: Genome FASTA file"
    echo "  --reads: Input reads FASTQ file"
    echo "  --output_prefix: Output file prefix"
    echo "  --chrom: Chromosome name"
    echo ""
    echo "Optional parameters:"
    echo "  --threads: Number of threads (default: 4)"
    echo "  --output_dir: Output directory (default: current directory)"
    exit 1
fi

# Set default values
if [[ -z "$threads" ]]; then
    threads=4
fi

if [[ -z "$output_dir" ]]; then
    output_dir="."
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

echo "=== FLAIR Collapse for Chromosome $chrom ==="
echo "Chromosome BED: $chrom_bed"
echo "Chromosome GTF: $chrom_gtf"
echo "Genome: $genome"
echo "Reads: $reads"
echo "Output prefix: $output_prefix"
echo "Threads: $threads"
echo "Output directory: $output_dir"
echo ""

# Check if input files exist
if [[ ! -f "$chrom_bed" ]]; then
    echo "Error: Chromosome BED file not found: $chrom_bed"
    exit 1
fi

if [[ ! -f "$chrom_gtf" ]]; then
    echo "Error: Chromosome GTF file not found: $chrom_gtf"
    exit 1
fi

if [[ ! -f "$genome" ]]; then
    echo "Error: Genome file not found: $genome"
    exit 1
fi

if [[ ! -f "$reads" ]]; then
    echo "Error: Reads file not found: $reads"
    exit 1
fi

# Run FLAIR collapse for this chromosome
echo "Running FLAIR collapse for chromosome $chrom..."

flair collapse \
    -g "$genome" \
    --gtf "$chrom_gtf" \
    -q "$chrom_bed" \
    -r "$reads" \
    --output "${output_dir}/${output_prefix}_${chrom}" \
    --check_splice \
    --stringent \
    --threads "$threads" \
    --generate_map \
    --annotation_reliant generate

echo "Collapse completed for chromosome $chrom"
echo "Output files:"
echo "  GTF: ${output_dir}/${output_prefix}_${chrom}.isoforms.gtf"
echo "  FASTA: ${output_dir}/${output_prefix}_${chrom}.isoforms.fa"
echo "  BED: ${output_dir}/${output_prefix}_${chrom}.isoforms.bed"
echo "" 