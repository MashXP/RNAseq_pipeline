#!/bin/bash

# [[01_core_alignment/scripts/00_download_data.sh]]
# This script downloads raw FASTQ files and metadata from Google Drive using rclone
# and performs extensive data integrity checks.

set -e
set -o pipefail

# Directories
FASTQ_DIR="/mnt/22T/vuphong/cancer_rna_drPhuong_data"

# Google Drive Info
GDRIVE_FOLDER_ID="1zcwHtgM3Ds3XPDezQvBW1oSI5Ek3wX7N"
# Default to "gdrive" if no argument is provided
REMOTE_NAME=${1:-"dinh.chanphngg@gmail.com"}

echo "=== Step 1: Downloading from Google Drive ==="
echo "Target Folder ID: $GDRIVE_FOLDER_ID"
echo "Destination: $FASTQ_DIR"
echo "Rclone Remote: $REMOTE_NAME"

mkdir -p "$FASTQ_DIR"

# Download using rclone. The --checksum flag ensures the transfer uses MD5 hashes to verify.
rclone copy "${REMOTE_NAME}:Hong P Nguyen" "$FASTQ_DIR" \
    --drive-shared-with-me \
    --progress \
    --checksum \
    --transfers 8 \
    --checkers 8

echo ""
echo "=== Step 2: Rclone Server-Side Integrity Check ==="
echo "Comparing destination against Google Drive source..."
rclone check "${REMOTE_NAME}:Hong P Nguyen" "$FASTQ_DIR" \
    --drive-shared-with-me \
    --one-way

echo ""
echo "=== Step 3: Local MD5 Integrity Check ==="
echo "Scanning for md5sums.txt, md5.txt, or *.md5 files..."
MD5_FILES=$(find "$FASTQ_DIR" -type f \( -name "*.md5" -o -name "md5sums.txt" -o -name "md5sum.txt" \))

if [ -n "$MD5_FILES" ]; then
    for md5file in $MD5_FILES; do
        target_dir=$(dirname "$md5file")
        echo "Verifying checksums in ${target_dir} using $(basename "$md5file")..."
        # Run md5sum within the directory of the md5 file
        if (cd "$target_dir" && md5sum -c "$(basename "$md5file")" --quiet); then
             echo " [OK] MD5 integrity check passed for $(basename "$md5file")"
        else
             echo " [ERROR] MD5 integrity check FAILED for $(basename "$md5file")!"
             echo "         Rerun this script to retry failed downloads."
             exit 1
        fi
    done
else
    echo "No MD5 manifest file found. Relying on rclone checks."
fi

echo ""
echo "=== Step 4: GZIP Integrity Check for FASTQ Files ==="
FASTQ_GZ_FILES=$(find "$FASTQ_DIR" -type f -name "*.fastq.gz" -o -name "*.fq.gz")

if [ -n "$FASTQ_GZ_FILES" ]; then
    echo "Testing gzip validity of downloaded .gz files..."
    for f in $FASTQ_GZ_FILES; do
        if gzip -t "$f" 2>/dev/null; then
             echo " [OK] gzip integrity passed: $(basename "$f")"
        else
             echo " [ERROR] Corrupted file detected: $(basename "$f")"
             echo "         Please delete $(basename "$f") and re-run this script."
             exit 1
        fi
    done
else
    echo "No .fastq.gz files found to run gzip -t check."
fi

echo ""
echo "=== Download and Integrity Checks Complete! ==="
