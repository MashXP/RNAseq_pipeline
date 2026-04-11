#!/bin/bash

# cleanup.sh
# Utility to wipe partial or corrupted results from the MiaPAca-2 pipeline.

# Exit on error
set -e

BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/_data"

echo "=== MiaPAca-2 Pipeline Cleanup Utility ==="
echo "Target Directory: $DATA_DIR"
echo ""

# Function to safely delete a directory if it exists
clean_dir() {
    local dir_name=$1
    local full_path="$DATA_DIR/$dir_name"
    if [ -d "$full_path" ]; then
        echo "Found $dir_name. Cleaning..."
        rm -rf "$full_path"
        echo "Done."
    else
        echo "Skipping $dir_name (not found)."
    fi
}

usage() {
    echo "Usage: bash cleanup.sh {trimmed|bam|counts|all}"
    echo "  trimmed : Wipes the fastq_trimmed directory (useful if files are corrupted)"
    echo "  bam     : Wipes the BAM alignment directory"
    echo "  counts  : Wipes the gene quantification results"
    echo "  all     : Wipes everything listed above"
    exit 1
}

if [ "$#" -eq 0 ]; then
    usage
fi

case $1 in
    trimmed)
        clean_dir "fastq_trimmed"
        ;;
    bam)
        clean_dir "bam"
        ;;
    counts)
        clean_dir "counts"
        ;;
    all)
        echo "WARNING: This will wipe ALL upstream results. FASTQ raw files will be kept."
        read -p "Are you sure? [y/N]: " confirm
        if [[ "$confirm" =~ ^[Yy]$ ]]; then
            clean_dir "fastq_trimmed"
            clean_dir "bam"
            clean_dir "counts"
        else
            echo "Abort."
            exit 1
        fi
        ;;
    *)
        usage
        ;;
esac

echo ""
echo "Cleanup complete. You can now re-run './run up' or your sbatch script."
