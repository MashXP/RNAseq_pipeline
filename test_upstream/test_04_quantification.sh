#!/bin/bash

# [[test_upstream/test_04_quantification.sh]]
# Quantifies reads per gene and biotype using featureCounts.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
# --- DEVIATION: Output to test-specific directories
BAM_DIR="$BASE_DIR/../_data/bam_test"
GENOME_DIR="$BASE_DIR/../_data/genome_test"
COUNTS_DIR="$BASE_DIR/../_data/counts_test"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../drPhuong_Sample_Data_Table.csv"

mkdir -p "$COUNTS_DIR"

# --- Resource Auto-detection ---
# --- DEVIATION: Fallback changed from 16 to 10 threads to match local machine specs
THREADS=${SLURM_CPUS_PER_TASK:-10}
echo "Resources: $THREADS threads."
# -------------------------------

echo ""
echo "================================================================================"
echo "   QUANTIFICATION: Starting featureCounts (per Group)"
echo "================================================================================"

groups=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk '{print $5}' | sort -u)

for group in $groups; do
    echo "-----------------------------------------------------"
    echo "Processing Group: $group"

    # Get species for this group
    species=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk -v grp="$group" '$5 == grp {print $4}' | head -n 1)
    
    if [ "$species" == "Human" ]; then
        # --- DEVIATION: Select test subset GTF instead of primary assembly
        GTF_FILE="$GENOME_DIR/Human/chr21_annotations.gtf"
    else
        # --- DEVIATION: Select test subset GTF instead of primary assembly
        GTF_FILE="$GENOME_DIR/Dog/chr38_annotations.gtf"
    fi

    if [ ! -f "$GTF_FILE" ]; then
        echo "Error: GTF file not found at $GTF_FILE."
        exit 1
    fi

    # Get all sample names for this group
    group_samples=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk -v grp="$group" '$5 == grp {print $1}')
    
    files=""
    for s in $group_samples; do
        if [ -f "$BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam" ]; then
            files="$files $BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam"
        fi
    done

    if [ -z "$files" ]; then
        echo "Warning: No BAM files found for group $group."
        continue
    fi

    echo "1. Counting Reads per Gene (gene_id)..."
    featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_id -a "$GTF_FILE" -o "$COUNTS_DIR/gene_counts_${group}.txt" $files

    echo "2. Counting Reads per Biotype (gene_biotype)..."
    featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_biotype -a "$GTF_FILE" -o "$COUNTS_DIR/biotype_counts_${group}.txt" $files

    if python3 -c "import pandas" &> /dev/null; then
        python3 "$UTILS_DIR/biotype_to_multiqc.py" "$COUNTS_DIR/biotype_counts_${group}.txt" "$COUNTS_DIR/biotype_counts_mqc_${group}.txt"
    fi 
done

echo "=== Quantification Complete ==="
