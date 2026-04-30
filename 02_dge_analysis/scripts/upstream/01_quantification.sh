#!/bin/bash

# [[02_dge_analysis/scripts/upstream/04_quantification.sh]]
# Quantifies reads per gene and biotype using featureCounts.

# Exit on error
set -e

# Directories
# --- Project Root Discovery ---
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}"
# ------------------------------
BAM_DIR="$PROJECT_ROOT/_data/bam"
GENOME_DIR="$PROJECT_ROOT/_data/genome"
FEATURECOUNTS_DIR="$PROJECT_ROOT/_data/featurecounts"
UTILS_DIR="$PROJECT_ROOT/utils"
CSV_FILE="$PROJECT_ROOT/_data/_metadata/drPhuong_Sample_Data_Table.csv"

mkdir -p "$FEATURECOUNTS_DIR"

# --- Resource Auto-detection ---
THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "Resources: $THREADS threads."
# -------------------------------

echo ""
echo "================================================================================"
echo "   QUANTIFICATION: Starting featureCounts (per Species)"
echo "================================================================================"

species_list=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk '{print $4}' | sort -u)

for sp in $species_list; do
    echo "-----------------------------------------------------"
    echo "Processing Species: $sp"
    
    # Use lowercase for paths
    sp_lower=$(echo "$sp" | tr '[:upper:]' '[:lower:]')
    mkdir -p "$FEATURECOUNTS_DIR/$sp_lower"

    if [ "$sp" == "Human" ]; then
        GTF_FILE="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf"
    else
        GTF_FILE="$GENOME_DIR/Canine/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"
    fi

    if [ ! -f "$GTF_FILE" ]; then
        echo "Error: GTF file not found at $GTF_FILE."
        exit 1
    fi

    # Get all sample names for this species
    sp_samples=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | awk -v species="$sp" '$4 == species {print $1}')
    
    files=""
    for s in $sp_samples; do
        if [ -f "$BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam" ]; then
            files="$files $BAM_DIR/$s/${s}_Aligned.sortedByCoord.out.bam"
        fi
    done

    if [ -z "$files" ]; then
        echo "Warning: No BAM files found for species $sp."
        continue
    fi

    echo "1. Counting Reads per Gene (gene_id)..."
    featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_id -a "$GTF_FILE" -o "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_counts.tsv" $files

    echo "2. Counting Reads per Biotype (gene_biotype)..."
    # Output to raw file first, then process with biotype_to_multiqc.py
    featureCounts -p -s 2 -T "$THREADS" -t exon -g gene_biotype -a "$GTF_FILE" -o "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts_raw.tsv" $files

    if python3 -c "import pandas" &> /dev/null; then
        echo "3. Processing Biotypes..."
        python3 "$UTILS_DIR/biotype_to_multiqc.py" \
            "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts_raw.tsv" \
            "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts.tsv"
    else
        echo "Warning: pandas not found. Skipping biotype processing. Raw file saved as ${sp_lower}_featurecounts_biotype_counts_raw.tsv"
        mv "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts_raw.tsv" "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts.tsv"
        mv "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts_raw.tsv.summary" "$FEATURECOUNTS_DIR/$sp_lower/${sp_lower}_featurecounts_biotype_counts.tsv.summary"
    fi 
done

echo "=== Quantification Complete ==="
