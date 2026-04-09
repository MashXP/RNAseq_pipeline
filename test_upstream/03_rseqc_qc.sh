#!/bin/bash

# test_upstream/03_rseqc_qc.sh
# TEST VERSION: Performs full QC on test BAM files (aligned with production).

set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
BAM_DIR="$BASE_DIR/../_data/bam_test"
GENOME_DIR="$BASE_DIR/../_data/genome_test"
QC_DIR="$BASE_DIR/../_data/qc_test"
BED_FILE="$GENOME_DIR/human_eb113_chr21.bed"

mkdir -p "$QC_DIR"

# Generate BED from GTF if missing
if [ ! -f "$BED_FILE" ]; then
    echo "Creating BED file for RSeQC (Chr21)..."
    python3 "$BASE_DIR/utils/gtf2bed.py" "$GENOME_DIR/chr21_annotations.gtf" > "$BED_FILE"
fi

echo "=== [TEST] Starting RSeQC Analysis ==="

find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam" | while read -r bam_file
do
    sample_name=$(basename "$(dirname "$bam_file")")
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample_name"
    
    SAMPLE_QC_DIR="$QC_DIR/$sample_name"
    mkdir -p "$SAMPLE_QC_DIR"

    # Index BAM if index doesn't exist (Aligned with production)
    if [ ! -f "${bam_file}.bai" ]; then
        echo "Indexing BAM..."
        samtools index "$bam_file"
    fi

    # 1. Gene Body Coverage
    echo "Running geneBody_coverage.py..."
    geneBody_coverage.py -i "$bam_file" -r "$BED_FILE" -o "$SAMPLE_QC_DIR/${sample_name}"

    # 2. Junction Annotation
    echo "Running junction_annotation.py..."
    junction_annotation.py -i "$bam_file" -r "$BED_FILE" -o "$SAMPLE_QC_DIR/${sample_name}" 2> "$SAMPLE_QC_DIR/${sample_name}.junction_annotation.log"

    # 3. Read Distribution
    echo "Running read_distribution.py..."
    read_distribution.py -i "$bam_file" -r "$BED_FILE" > "$SAMPLE_QC_DIR/${sample_name}.read_distribution.txt"

    echo "Done with QC for $sample_name"
done

echo "=== [TEST] QC Complete ==="
