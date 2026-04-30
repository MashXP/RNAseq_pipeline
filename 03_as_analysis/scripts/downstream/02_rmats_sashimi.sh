#!/bin/bash

# [[03_as_analysis/scripts/downstream/02_rmats_sashimi.sh]]
# Goal: Generate Sashimi plots for specific rMATS events.

set -e

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <cell_line> <comparison> <as_type> <event_id>"
    exit 1
fi

CELL_LINE=$1
COMP=$2
AS_TYPE=$3
EVENT_ID=$4

# --- Project Root Discovery ---
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}"
# ------------------------------

BAM_LISTS_DIR="$PROJECT_ROOT/_data/rmats/bam_lists"
RESULTS_DIR="$PROJECT_ROOT/_data/rmats/results"
GENOME_DIR="$PROJECT_ROOT/_data/genome"

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
OUT_DIR="$PROJECT_ROOT/results/03_as_analysis/sashimi/$CELL_LINE/$COMP/${AS_TYPE}_${EVENT_ID}"
mkdir -p "$OUT_DIR"

# Labels
GROUP1=$(echo "$COMP" | cut -d'_' -f1,2)
GROUP2=$(echo "$COMP" | cut -d'_' -f4,5)

echo ">>> Generating Sashimi Plot for $CELL_LINE | $COMP | $AS_TYPE | ID: $EVENT_ID"

# Check if inputs exist
if [ ! -f "$B1" ] || [ ! -f "$B2" ] || [ ! -f "$EVENT_FILE" ]; then
    echo "ERROR: Missing input files for Sashimi plot."
    echo "B1: $B1"
    echo "B2: $B2"
    echo "Event: $EVENT_FILE"
    exit 1
fi

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
