#!/bin/bash

# [[test_upstream/05_multiqc.sh]]
# TEST VERSION: Aggregates QC reports.
# Logic: Runs MultiQC on bam_test, qc_test, and counts_test directories.
# Logic: Runs featureCounts against Chr21 annotations for GeneID and Biotype.
# Divergence from Prod: Uses test-specific paths for results.

set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
BAM_DIR="$BASE_DIR/../_data/bam_test"
GENOME_DIR="$BASE_DIR/../_data/genome_test"
COUNTS_DIR="$BASE_DIR/../_data/counts_test"
UTILS_DIR="$BASE_DIR/utils"
GTF_FILE="$GENOME_DIR/chr21_annotations.gtf"

mkdir -p "$COUNTS_DIR"

# Find BAM files
files=$(find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam")
if [ -z "$files" ]; then
    echo "Error: No BAM files found in $BAM_DIR."
    exit 1
fi

echo "=== [TEST] Starting Quantification (Chr21) ==="

# 1. Gene Counts
echo "Counting Reads per Gene..."
featureCounts -p -T 4 -t exon -g gene_id -a "$GTF_FILE" -o "$COUNTS_DIR/gene_counts.txt" $files

# 2. Biotype Counts
echo "Counting Reads per Biotype..."
featureCounts -p -T 4 -t exon -g gene_biotype -a "$GTF_FILE" -o "$COUNTS_DIR/biotype_counts.txt" $files

# 3. Format Biotype for MultiQC
if python3 -c "import pandas" &> /dev/null; then
    python3 "$UTILS_DIR/biotype_to_multiqc.py" "$COUNTS_DIR/biotype_counts.txt" "$COUNTS_DIR/biotype_counts_mqc.txt"
fi

echo "=== [TEST] Quantification Complete ==="
