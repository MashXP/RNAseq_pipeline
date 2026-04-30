#!/bin/bash

# [[scripts_upstream/05_multiqc.sh]]
# Aggregates all reports into a final MultiQC summary.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/../_data"
MULTIQC_DIR="$DATA_DIR/multiqc"
BAM_DIR="$DATA_DIR/bam"
QC_DIR="$DATA_DIR/qc"
QC_RAW_DIR="$DATA_DIR/qc_raw"
COUNTS_DIR="$DATA_DIR/featurecounts"

echo ""
echo "================================================================================"
echo "   SUMMARY: Running MultiQC"
echo "================================================================================"
mkdir -p "$MULTIQC_DIR"

# Run MultiQC on specific result directories to catch:
# 1. FastQC (raw and trimmed if run)
# 2. STAR alignment logs
# 3. RSeQC logs
# 4. featureCounts summaries
# 5. Biotype custom content

multiqc "$BAM_DIR" "$QC_DIR" "$QC_RAW_DIR" "$COUNTS_DIR" -o "$MULTIQC_DIR" -f -n "rna_seq_pipeline_summary" -c "$BASE_DIR/../multiqc_config.yaml"

echo "=== MultiQC Complete ==="
echo "Report generated at: $MULTIQC_DIR/rna_seq_pipeline_summary.html"
