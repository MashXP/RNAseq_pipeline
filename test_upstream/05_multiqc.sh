#!/bin/bash

# test_upstream/05_multiqc.sh
# TEST VERSION: Aggregates QC results from test directories.

set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/../_data"
BAM_DIR="$DATA_DIR/bam_test"
QC_DIR="$DATA_DIR/qc_test"
COUNTS_DIR="$DATA_DIR/counts_test"
OUT_DIR="$DATA_DIR/multiqc_test"

mkdir -p "$OUT_DIR"

echo "=== [TEST] Running MultiQC ==="

multiqc "$BAM_DIR" "$QC_DIR" "$COUNTS_DIR" -o "$OUT_DIR" -f

echo "=== [TEST] MultiQC Complete. Report available in $OUT_DIR ==="
