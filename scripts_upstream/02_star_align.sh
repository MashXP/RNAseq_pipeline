#!/bin/bash

# [[scripts_upstream/02_star_align.sh]]
# Performs paired-end alignment using STAR for all samples in the CSV table.

# Exit on error
set -e
set -o pipefail

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
FASTQ_DIR="$BASE_DIR/../_data/fastq_trimmed"
INDEX_DIR="$BASE_DIR/../_data/index"
ALIGN_DIR="$BASE_DIR/../_data/bam"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../drPhuong_Sample_Data_Table.csv"

echo ""
echo "================================================================================"
echo "   ALIGNMENT: Starting STAR Alignment"
echo "================================================================================"
mkdir -p "$ALIGN_DIR"

# --- Resource Auto-detection ---
# Default to 32 threads or respect SLURM
THREADS=${SLURM_CPUS_PER_TASK:-32}

# Detect RAM and set limits
if [ -n "$SLURM_MEM_PER_NODE" ]; then
    MEM_GB=$(($SLURM_MEM_PER_NODE / 1024))
else
    # Fallback to system check or default to 128
    MEM_GB=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo 128)
fi

# Allocate RAM for STAR sorting (set to ~65-70% of total RAM to avoid swap)
# Value must be in bytes. 1GB approx = 1,000,000,000 bytes
SORT_RAM_BYTES=$(($MEM_GB * 700000000))

echo "Resources: $THREADS threads, Sort RAM: $SORT_RAM_BYTES bytes (~${MEM_GB}G total)."
# -------------------------------

# Parse CSV and loop through samples
# format: sample_name r1_file r2_file
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2 species group
do
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample"
    
    # Check if FASTQ files exist AND are valid GZIPs
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "Error: FASTQ files not found in $FASTQ_DIR:"
        echo "  R1: $r1"
        echo "  R2: $r2"
        echo "Skipping sample $sample."
        continue
    fi

    # Quick integrity check
    if ! gzip -t "$FASTQ_DIR/$r1" 2>/dev/null || ! gzip -t "$FASTQ_DIR/$r2" 2>/dev/null; then
        echo "--------------------------------------------------------------------------------"
        echo "CRITICAL ERROR: Trimmed FASTQ files for $sample are CORRUPTED (unexpected EOF)."
        echo "ACTION REQUIRED: Delete the corrupted files in $FASTQ_DIR and re-run trimming."
        echo "--------------------------------------------------------------------------------"
        exit 1
    fi

    OUT_DIR="$ALIGN_DIR/$sample"
    mkdir -p "$OUT_DIR"

    # Skip if BAM already exists to save time/resources
    BAM_OUT="$OUT_DIR/${sample}_Aligned.sortedByCoord.out.bam"
    if [ -f "$BAM_OUT" ]; then
        echo "BAM for $sample already exists ($BAM_OUT). Skipping."
        continue
    fi

    # Run STAR
    echo "Running STAR for $species..."
    STAR --genomeDir "$INDEX_DIR/$species" \
         --readFilesIn "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix "$OUT_DIR/${sample}_" \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN "$THREADS" \
         --limitBAMsortRAM "$SORT_RAM_BYTES" \
         --quantMode GeneCounts || { echo "ERROR: STAR failed for $sample. Aborting."; exit 1; }
    
    # Verification of BAM generation
    BAM_OUT="$OUT_DIR/${sample}_Aligned.sortedByCoord.out.bam"
    if [ ! -f "$BAM_OUT" ]; then
        echo "ERROR: STAR failed to generate BAM file for $sample: $BAM_OUT"
        exit 1
    fi
    
    # Check if BAM is non-empty
    if [ ! -s "$BAM_OUT" ]; then
        echo "ERROR: Generated BAM file for $sample is EMPTY."
        exit 1
    fi

    echo "Done with $sample"
    # Report Mapping Rate
    map_rate=$(grep "Uniquely mapped reads %" "$OUT_DIR/${sample}_Log.final.out" | cut -f2)
    echo "STATUS: Uniquely Mapped Reads for $sample: $map_rate"
done

echo "=== STAR Alignment Complete ==="
