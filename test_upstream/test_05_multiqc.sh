#!/bin/bash
 
# [[test_upstream/test_05_multiqc.sh]]
# Aggregates all reports into a final MultiQC summary.
 
# Exit on error
set -e
 
# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/../_data"
# --- DEVIATION: Points to test output directory
MULTIQC_DIR="$DATA_DIR/multiqc_test"
# --- DEVIATION: Points to test BAM directory
BAM_DIR="$DATA_DIR/bam_test"
# --- DEVIATION: Points to test QC directory
QC_DIR="$DATA_DIR/qc_test"
QC_RAW_DIR="$DATA_DIR/qc_raw_test"
# --- DEVIATION: Points to test Counts directory
COUNTS_DIR="$DATA_DIR/counts_test"

echo ""
echo "================================================================================"
echo "   SUMMARY: Running MultiQC"
echo "================================================================================"
mkdir -p "$MULTIQC_DIR"

# Run MultiQC on specific result directories to catch:
# 1. FastQC (raw and trimmed if run)
# 2. STAR alignment logs
# 3. Picard logs (Alignment QC)
# 4. featureCounts summaries
# 5. Biotype custom content

# --- DEVIATION: Report filename includes _test suffix
multiqc "$BAM_DIR" "$QC_DIR" "$QC_RAW_DIR" "$COUNTS_DIR" -o "$MULTIQC_DIR" -f -n "rna_seq_pipeline_summary_test" 
echo "=== MultiQC Complete ==="
# --- DEVIATION: Report filename includes _test suffix
echo "Report generated at: $MULTIQC_DIR/rna_seq_pipeline_summary_test.html"
