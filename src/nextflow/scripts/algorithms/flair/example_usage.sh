#!/bin/bash

# Example usage of the simplified FLAIR pipeline scripts
# This script demonstrates how to run the FLAIR pipeline with the new simplified scripts

echo "=== FLAIR Pipeline Example Usage ==="
echo ""

# Example parameters - modify these for your specific use case
WORKING_DIR="/path/to/working/directory"
GENOME="/path/to/genome.fa"
ANNOTATION="/path/to/annotation.gtf"
READS="/path/to/reads.fastq"
SR_JUNCTIONS="/path/to/sr_junctions.bed"
OUTPUT_PREFIX="sample_name"
THREADS=8
MAX_PARALLEL=4

echo "Example parameters:"
echo "  Working directory: $WORKING_DIR"
echo "  Genome: $GENOME"
echo "  Annotation: $ANNOTATION"
echo "  Reads: $READS"
echo "  Short read junctions: $SR_JUNCTIONS"
echo "  Output prefix: $OUTPUT_PREFIX"
echo "  Threads per job: $THREADS"
echo "  Max parallel jobs: $MAX_PARALLEL"
echo ""

# Example 1: Run the complete pipeline using the main script
echo "=== Example 1: Complete Pipeline ==="
echo "Running the complete FLAIR pipeline..."
echo ""

# Uncomment the following lines to run the actual pipeline:
# ./run_flair_simplified_v2.sh \
#     --wd "$WORKING_DIR" \
#     --genome "$GENOME" \
#     --annotation "$ANNOTATION" \
#     --reads "$READS" \
#     --sr_junctions "$SR_JUNCTIONS" \
#     --output_prefix "$OUTPUT_PREFIX" \
#     --threads "$THREADS" \
#     --max_parallel "$MAX_PARALLEL"

echo ""

# Example 2: Run individual steps manually
echo "=== Example 2: Manual Step-by-Step Execution ==="
echo ""

# Step 1: FLAIR Align (uncomment when ready)
echo "Step 1: FLAIR Align"
echo "flair align -g $GENOME -r $READS -o $WORKING_DIR/align/$OUTPUT_PREFIX -t $THREADS"
echo ""

# Step 2: Convert BAM to BED12
echo "Step 2: Convert BAM to BED12"
echo "bam2Bed12 $WORKING_DIR/align/$OUTPUT_PREFIX.bam > $WORKING_DIR/align/$OUTPUT_PREFIX.bed"
echo ""

# Step 3: FLAIR Correct
echo "Step 3: FLAIR Correct"
echo "flair correct -q $WORKING_DIR/align/$OUTPUT_PREFIX.bed -f $ANNOTATION -g $GENOME --output $WORKING_DIR/correct/$OUTPUT_PREFIX --threads $THREADS --shortread $SR_JUNCTIONS"
echo ""

# Step 4: Split by chromosome and run collapse
echo "Step 4: Split by chromosome and run collapse"
echo "# Get chromosomes from corrected BED file"
echo "chromosomes=\$(cut -f1 $WORKING_DIR/correct/${OUTPUT_PREFIX}_all_corrected.bed | sort | uniq)"
echo ""

echo "# For each chromosome, run collapse in parallel:"
echo "for chrom in \$chromosomes; do"
echo "    ./run_collapse_chromosome.sh \\"
echo "        --chrom_bed \$WORKING_DIR/chromosomes/${OUTPUT_PREFIX}_\${chrom}_corrected.bed \\"
echo "        --chrom_gtf \$WORKING_DIR/chromosomes/\$(basename $ANNOTATION .gtf)_\${chrom}.gtf \\"
echo "        --genome $GENOME \\"
echo "        --reads $READS \\"
echo "        --output_prefix $OUTPUT_PREFIX \\"
echo "        --chrom \$chrom \\"
echo "        --threads $THREADS \\"
echo "        --output_dir \$WORKING_DIR/chromosomes &"
echo "done"
echo "wait"
echo ""

# Step 5: Join chromosome results
echo "Step 5: Join chromosome results"
echo "./join_chromosome_results.sh \\"
echo "    --chromosomes \"\$chromosomes\" \\"
echo "    --output_prefix $OUTPUT_PREFIX \\"
echo "    --chrom_dir \$WORKING_DIR/chromosomes \\"
echo "    --collapse_dir \$WORKING_DIR/collapse"
echo ""

# Step 6: FLAIR Quantify
echo "Step 6: FLAIR Quantify"
echo "flair quantify \\"
echo "    -i \$WORKING_DIR/collapse/${OUTPUT_PREFIX}.isoforms.fa \\"
echo "    --reads_manifest \$WORKING_DIR/quant/${OUTPUT_PREFIX}.read_manifest.tsv \\"
echo "    --threads $THREADS \\"
echo "    --output \$WORKING_DIR/quant/$OUTPUT_PREFIX \\"
echo "    --tpm \\"
echo "    --check_splice \\"
echo "    --stringent \\"
echo "    --isoform_bed \$WORKING_DIR/collapse/${OUTPUT_PREFIX}.isoforms.bed"
echo ""

echo "=== Example 3: Using with SLURM ==="
echo ""

# Example SLURM job script
cat << 'EOF'
#!/bin/bash
#SBATCH --job-name=flair_pipeline
#SBATCH --qos medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH -t 7-00:00:00
#SBATCH -o log/flair_%j.out

# Load required modules
module load anaconda
source activate flair_conda_env

# Run the pipeline
./run_flair_simplified_v2.sh \
    --wd "/path/to/working/directory" \
    --genome "/path/to/genome.fa" \
    --annotation "/path/to/annotation.gtf" \
    --reads "/path/to/reads.fastq" \
    --sr_junctions "/path/to/sr_junctions.bed" \
    --output_prefix "sample_name" \
    --threads 8 \
    --max_parallel 4
EOF

echo ""
echo "=== Notes ==="
echo "1. Make sure FLAIR is installed and activated in your conda environment"
echo "2. Ensure all input files exist and are accessible"
echo "3. The BAM file from FLAIR align should be placed in the align directory"
echo "4. Adjust thread counts and memory based on your system resources"
echo "5. The pipeline automatically handles chromosome splitting and parallel execution"
echo "6. All outputs are organized in subdirectories: align/, correct/, collapse/, quant/, chromosomes/"
echo "" 