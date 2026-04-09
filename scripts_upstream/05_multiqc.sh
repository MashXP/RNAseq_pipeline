#!/bin/bash

# 05_multiqc.sh
# Aggregates all reports into a final MultiQC summary.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/../_data"
MULTIQC_DIR="$DATA_DIR/multiqc"

echo "=== Running MultiQC ==="
mkdir -p "$MULTIQC_DIR"

# Run MultiQC on the entire _data directory to catch:
# 1. FastQC (if run)
# 2. STAR alignment logs
# 3. RSeQC logs
# 4. featureCounts summaries
# 5. Biotype custom content

multiqc "$DATA_DIR" -o "$MULTIQC_DIR" -n "multiqc_summary_report"

echo "=== MultiQC Complete ==="
echo "Report generated at: $MULTIQC_DIR/multiqc_summary_report.html"
