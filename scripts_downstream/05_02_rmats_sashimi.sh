#!/bin/bash

# [[scripts_downstream/05_02_rmats_sashimi.sh]]
# Goal: Generate Sashimi plots for specific rMATS events.
#
# Usage: bash 07_rmats_sashimi.sh <cell_line> <comparison> <as_type> <event_id>
# Example: bash 07_rmats_sashimi.sh CNK89 Kromastat_6nM_vs_Romidepsin_6nM SE 2318

set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <cell_line> <comparison> <as_type> <event_id>"
    exit 1
fi

CELL_LINE=$1
COMP=$2
AS_TYPE=$3
EVENT_ID=$4

# --- Paths ---
BASE_DIR=$(dirname "$(realpath "$0")")
PIPELINE_ROOT="$BASE_DIR/.."
BAM_LISTS_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/bam_lists"
RESULTS_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/results"
GENOME_DIR="$PIPELINE_ROOT/_data/genome"

# Determine Species/GTF
if [[ "$CELL_LINE" == "UL1" || "$CELL_LINE" == "CNK89" ]]; then
    GTF="$GENOME_DIR/Canine/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"
else
    GTF="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf"
fi

# Input files
B1="$BAM_LISTS_DIR/$CELL_LINE/$COMP/b1.txt"
B2="$BAM_LISTS_DIR/$CELL_LINE/$COMP/b2.txt"
EVENT_FILE="$RESULTS_DIR/$CELL_LINE/$COMP/raw/${AS_TYPE}.MATS.JCEC.txt"

# Output
DOWNSTREAM_RESULTS="$PIPELINE_ROOT/results/alt_splicing/sashimi"
OUT_DIR="$DOWNSTREAM_RESULTS/$CELL_LINE/$COMP/${AS_TYPE}_${EVENT_ID}"
mkdir -p "$OUT_DIR"

# Labels
GROUP1=$(echo "$COMP" | cut -d'_' -f1,2)
GROUP2=$(echo "$COMP" | cut -d'_' -f4,5)

echo ">>> Generating Sashimi Plot for $CELL_LINE | $COMP | $AS_TYPE | ID: $EVENT_ID"

# Run rmats2sashimiplot
# Note: It requires coordinate-sorted BAM files with .bai indexes.
# The BAM lists (b1.txt, b2.txt) contain paths to BAMs.
# We need to convert them to comma-separated for rmats2sashimiplot.

BAMS1=$(cat "$B1" | tr '\n' ',' | sed 's/,$//')
BAMS2=$(cat "$B2" | tr '\n' ',' | sed 's/,$//')

rmats2sashimiplot \
    --b1 "$BAMS1" \
    --b2 "$BAMS2" \
    --event-type "$AS_TYPE" \
    -e "$EVENT_FILE" \
    --l1 "$GROUP1" \
    --l2 "$GROUP2" \
    --exon_s 1 \
    --intron_s 1 \
    -o "$OUT_DIR" \
    || {
        echo "ERROR: rmats2sashimiplot failed for $COMP | $AS_TYPE | ID: $EVENT_ID"
        exit 1
    }

echo "[OK] Sashimi plot generated in: $OUT_DIR"
