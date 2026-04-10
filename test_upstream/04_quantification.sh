#!/bin/bash

# [[test_upstream/04_quantification.sh]]
# TEST VERSION: Quantifies reads per gene and biotype using featureCounts.
# Logic: Runs featureCounts against Chr21 annotations for GeneID and Biotype.
# Divergence from Prod: Uses Chr21 GTF and outputs to counts_test directory.

set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
BAM_DIR="$BASE_DIR/../_data/bam_test"
GENOME_DIR="$BASE_DIR/../_data/genome_test"
COUNTS_DIR="$BASE_DIR/../_data/counts_test"
UTILS_DIR="$BASE_DIR/utils"
GTF_FILE="$GENOME_DIR/chr21_annotations.gtf"

mkdir -p "$COUNTS_DIR"

if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file not found at $GTF_FILE. Run 01_genome_prep.sh first."
    exit 1
fi

# Find all BAM files
files=$(find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam")
if [ -z "$files" ]; then
    echo "Error: No BAM files found in $BAM_DIR."
    exit 1
fi

echo "=== [TEST] Starting Quantification with featureCounts (Chr21) ==="

# 1. Gene Counts (for Differential Expression)
# -p: paired-end
# -T: threads
# -t: feature type (exon)
# -g: attribute type (gene_id)
echo "-----------------------------------------------------"
echo "1. Counting Reads per Gene (gene_id)..."
featureCounts -p -T 4 -t exon -g gene_id -a "$GTF_FILE" -o "$COUNTS_DIR/gene_counts.txt" $files

# 2. Biotype Counts (for QC)
echo "-----------------------------------------------------"
echo "2. Counting Reads per Biotype (gene_biotype)..."
featureCounts -p -T 4 -t exon -g gene_biotype -a "$GTF_FILE" -o "$COUNTS_DIR/biotype_counts.txt" $files

# 3. Format Biotype counts for MultiQC
echo "-----------------------------------------------------"
echo "3. Formatting biotype results for MultiQC..."

# Pre-check for pandas dependency
if ! python3 -c "import pandas" &> /dev/null; then
    echo "Warning: Python 'pandas' library not found."
    echo "Skipping biotype formatting. To enable this, run: pip install pandas"
else
    python3 "$UTILS_DIR/biotype_to_multiqc.py" "$COUNTS_DIR/biotype_counts.txt" "$COUNTS_DIR/biotype_counts_mqc.txt"
fi

echo "=== [TEST] Quantification Complete ==="
