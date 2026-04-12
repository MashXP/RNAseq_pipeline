#!/bin/bash

# [[scripts_upstream/03_rseqc_qc.sh]]
# Runs RSeQC suite on aligned BAM files.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
BAM_DIR="$BASE_DIR/../_data/bam"
GENOME_DIR="$BASE_DIR/../_data/genome"
QC_DIR="$BASE_DIR/../_data/rseqc"
CSV_FILE="$BASE_DIR/../drPhuong_Sample_Data_Table.csv"

mkdir -p "$QC_DIR"

echo "=== Starting RSeQC Analysis ==="

# Find all BAM files in subdirectories
find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam" | while read -r bam_file
do
    sample_name=$(basename "$(dirname "$bam_file")")
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample_name"
    
    SAMPLE_QC_DIR="$QC_DIR/$sample_name"
    mkdir -p "$SAMPLE_QC_DIR"

    # Determine species
    sample_species=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | grep "^$sample_name " | awk '{print $4}')
    if [ "$sample_species" == "Human" ]; then
        BED_FILE="$GENOME_DIR/Homo_sapiens.GRCh38.113.bed"
    else
        BED_FILE="$GENOME_DIR/Canis_lupus_familiaris.ROS_Cfam_1.0.113.bed"
    fi

    if [ ! -f "$BED_FILE" ]; then
        echo "Error: BED file not found at $BED_FILE for $sample_name ($sample_species)."
        continue
    fi

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
