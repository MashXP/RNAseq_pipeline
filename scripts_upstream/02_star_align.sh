#!/bin/bash

# [[scripts_upstream/02_star_align.sh]]
# Performs paired-end alignment using STAR for all samples in the CSV table.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
FASTQ_DIR="$BASE_DIR/../_data/fastq_trimmed"
INDEX_DIR="$BASE_DIR/../_data/index"
ALIGN_DIR="$BASE_DIR/../_data/bam"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../MiaPAca-2_Sample_Data_Table.csv"

echo "=== Starting STAR Alignment ==="
mkdir -p "$ALIGN_DIR"

# Parse CSV and loop through samples
# format: sample_name r1_file r2_file
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2
do
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample"
    
    # Check if FASTQ files exist
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "Error: FASTQ files not found in $FASTQ_DIR:"
        echo "  R1: $r1"
        echo "  R2: $r2"
        echo "Skipping sample $sample."
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
         --runThreadN 32 \
         --quantMode GeneCounts || { echo "ERROR: STAR failed for $sample. Aborting."; exit 1; }
    
    echo "Done with $sample"
    # Report Mapping Rate
    map_rate=$(grep "Uniquely mapped reads %" "$OUT_DIR/${sample}_Log.final.out" | cut -f2)
    echo "STATUS: Uniquely Mapped Reads for $sample: $map_rate"
done

echo "=== STAR Alignment Complete ==="
