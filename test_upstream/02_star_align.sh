#!/bin/bash

# test_upstream/02_star_align.sh
# TEST VERSION: Performs alignment using the Chr21 index.

set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
FASTQ_DIR="$BASE_DIR/../_data/fastq_trimmed_test"
INDEX_DIR="$BASE_DIR/../_data/index_test"
ALIGN_DIR="$BASE_DIR/../_data/bam_test"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../MiaPAca-2_Sample_Data_Table.csv"

echo "=== [TEST] Starting STAR Alignment (Chr21) ==="
mkdir -p "$ALIGN_DIR"

# Parse CSV and loop through samples
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2
do
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample"
    
    # Check if FASTQ files exist
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "Error: FASTQ files not found in $FASTQ_DIR. Skipping."
        continue
    fi

    OUT_DIR="$ALIGN_DIR/$sample"
    mkdir -p "$OUT_DIR"

    # Run STAR
    echo "Running STAR..."
    STAR --genomeDir "$INDEX_DIR" \
         --readFilesIn "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix "$OUT_DIR/${sample}_" \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 4 \
         --quantMode GeneCounts

    echo "Done with $sample"
done

echo "=== [TEST] STAR Alignment Complete ==="
