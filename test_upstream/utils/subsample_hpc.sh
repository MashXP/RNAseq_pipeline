#!/bin/bash

# subsample_hpc.sh
# Run this on the HPC inside the directory containing 'data/drPhuong'
# It creates a flat folder 'test_subset' with 100k reads per sample.

# Set the base directory (where data/ resides)
DATA_ROOT="/mnt/22T/vuphong/data/drPhuong"
OUTPUT_DIR="$HOME/test_subset"

if [ ! -d "$DATA_ROOT" ]; then
    echo "Error: Directory $DATA_ROOT not found."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "=== Starting Subsampling for Verification ==="
echo "Output directory: $OUTPUT_DIR"

find "$DATA_ROOT" -name "*.fastq.gz" | while read -r f; do
    fname=$(basename "$f")
    echo "Processing $fname..."
    
    # Extract first 400,000 lines (100,000 reads)
    # Using zcat and gzip.
    zcat "$f" | head -n 400000 | gzip > "$OUTPUT_DIR/$fname"
done

echo "=== Subsampling Complete ==="
echo "You can now download the '$OUTPUT_DIR' folder to your local machine."
echo "Suggested download command (from local):"
echo "  scp -r phongdinh@worker-01:/mnt/22T/vuphong/test_subset /path/to/MiaPAca-2_pipeline/_data/fastq"
