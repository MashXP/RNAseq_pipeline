#!/bin/bash

# 03_rseqc_qc.sh
# Runs RSeQC suite on aligned BAM files.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
BAM_DIR="$BASE_DIR/../_data/bam"
GENOME_DIR="$BASE_DIR/../_data/genome"
QC_DIR="$BASE_DIR/../_data/rseqc"
BED_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.113.bed"

mkdir -p "$QC_DIR"

if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found at $BED_FILE. Run 01_genome_prep.sh first."
    exit 1
fi

echo "=== Starting RSeQC Analysis ==="

# Find all BAM files in subdirectories
find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam" | while read -r bam_file
do
    sample_name=$(basename "$(dirname "$bam_file")")
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample_name"
    
    SAMPLE_QC_DIR="$QC_DIR/$sample_name"
    mkdir -p "$SAMPLE_QC_DIR"

    # Index BAM if index doesn't exist
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

echo "=== RSeQC Analysis Complete ==="
