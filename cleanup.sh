#!/bin/bash

# cleanup.sh
# Utility to wipe partial or corrupted results from the RNA-seq pipeline.

# Exit on error
set -e

BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/_data"

echo "=== RNA-seq Pipeline Cleanup Utility ==="
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
    echo "Usage: bash cleanup.sh {trimmed|bam|counts|test|all}"
    echo "  trimmed : Wipes the fastq_trimmed directory (useful if files are corrupted)"
    echo "  bam     : Wipes the BAM alignment directory"
    echo "  counts  : Wipes the gene quantification results"
    echo "  test    : Wipes all test-related output directories"
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
    test)
        echo "Cleaning all test output directories..."
        clean_dir "bam_test"
        clean_dir "counts_test"
        clean_dir "fastq_trimmed_test"
        clean_dir "genome_test"
        clean_dir "index_test"
        clean_dir "multiqc_test"
        clean_dir "qc_test"
        if [ -d "$BASE_DIR/results_test" ]; then
            echo "Found results_test. Cleaning..."
            rm -rf "$BASE_DIR/results_test"
            echo "Done."
        else
            echo "Skipping results_test (not found)."
        fi
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
