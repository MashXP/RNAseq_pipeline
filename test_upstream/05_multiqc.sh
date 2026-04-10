#!/bin/bash

# test_upstream/05_multiqc.sh
# TEST VERSION: Aggregates all reports into a final MultiQC summary for test data.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/../_data"
MULTIQC_DIR="$DATA_DIR/multiqc_test"
BAM_DIR="$DATA_DIR/bam_test"
QC_DIR="$DATA_DIR/qc_test"
COUNTS_DIR="$DATA_DIR/counts_test"

echo "=== [TEST] Running MultiQC ==="
mkdir -p "$MULTIQC_DIR"

# Run MultiQC on specific test result directories to catch:
# 1. FastQC (if run)
# 2. STAR alignment logs
# 3. RSeQC logs
# 4. featureCounts summaries
# 5. Biotype custom content

multiqc "$BAM_DIR" "$QC_DIR" "$COUNTS_DIR" -o "$MULTIQC_DIR" -f -n "multiqc_test_report"

echo "=== [TEST] MultiQC Complete ==="
echo "Report generated at: $MULTIQC_DIR/multiqc_test_report.html"
