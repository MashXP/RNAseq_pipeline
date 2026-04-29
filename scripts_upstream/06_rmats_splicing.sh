#!/bin/bash

# [[scripts_upstream/06_rmats_splicing.sh]]
# rMATS-turbo Alternative Splicing Pipeline
#
# Parallelization: Picard-style internal background job management.
# This script requests a large block of resources and manages concurrent 
# comparisons internally to stay within a predictable resource footprint.

set -e
set -o pipefail

# ──────────────────────────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────────────────────────

BASE_DIR=$(dirname "$(realpath "$0")")
PIPELINE_ROOT="$BASE_DIR/.."
UTILS_DIR="$BASE_DIR/utils"

GENOME_DIR="$PIPELINE_ROOT/_data/genome"
GTF_HUMAN="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf"
GTF_CANINE="$GENOME_DIR/Canine/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"

BAM_LISTS_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/bam_lists"
RESULTS_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/results"
TMP_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/tmp"
LOG_DIR="$PIPELINE_ROOT/_data/logs/rmats_array"

mkdir -p "$RESULTS_DIR" "$TMP_DIR" "$LOG_DIR"

# ──────────────────────────────────────────────────────────────────────────────
# PARAMETERS & RESOURCES
# ──────────────────────────────────────────────────────────────────────────────

READ_LENGTH=150
LIB_TYPE="fr-secondstrand"
READ_TYPE="paired"
CSTAT=0.05
ANCHOR_LENGTH=1

# --- Resource Auto-detection (Picard-style) ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_TASKS=4  # Number of comparisons to run concurrently
THREADS_PER_TASK=$((TOTAL_CPUS / PARALLEL_TASKS))
[ "$THREADS_PER_TASK" -lt 1 ] && THREADS_PER_TASK=1

echo "========================================================================"
echo "   rMATS-turbo Alternative Splicing Pipeline"
echo "========================================================================"
echo "Resources : $TOTAL_CPUS total CPUs"
echo "Parallel  : $PARALLEL_TASKS concurrent comparisons ($THREADS_PER_TASK threads each)"
echo "Read type : $READ_TYPE, length = $READ_LENGTH bp"
echo "Library   : $LIB_TYPE"
echo "------------------------------------------------------------------------"

# ──────────────────────────────────────────────────────────────────────────────
# PRE-FLIGHT
# ──────────────────────────────────────────────────────────────────────────────

if ! command -v rmats.py >/dev/null 2>&1; then
    if ! (command -v python3 >/dev/null 2>&1 && python3 -c "import rmats" 2>/dev/null); then
        echo "[ERROR] rmats.py not found. Is the conda environment active?"
        exit 1
    fi
fi

# ──────────────────────────────────────────────────────────────────────────────
# STEP 0 — Generate b1.txt / b2.txt BAM lists
# ──────────────────────────────────────────────────────────────────────────────

echo ">>> STEP 0: Generating BAM list files from metadata CSV..."
python3 "$UTILS_DIR/generate_bam_lists.py"

# ──────────────────────────────────────────────────────────────────────────────
# TASK FUNCTION
# ──────────────────────────────────────────────────────────────────────────────

# Function to run a full rMATS analysis for one comparison
run_comparison() {
    local CELL_LINE=$1
    local COMP=$2
    local SPECIES=$3
    local GTF=$4

    local COMP_BAM_DIR="$BAM_LISTS_DIR/$CELL_LINE/$COMP"
    local B1="$COMP_BAM_DIR/b1.txt"
    local B2="$COMP_BAM_DIR/b2.txt"
    local TMP_OUT="$TMP_DIR/$CELL_LINE/$COMP"
    local POST_OUT="$RESULTS_DIR/$CELL_LINE/$COMP"

    if [ ! -f "$B1" ] || [ ! -f "$B2" ]; then
        echo "[WARN] Skipping $CELL_LINE / $COMP — BAM lists missing."
        return 0
    fi

    mkdir -p "$TMP_OUT" "$POST_OUT"

    echo "[$(date +%T)] [$CELL_LINE|$COMP] Starting Prep..."
    rmats.py \
        --gtf "$GTF" --b1 "$B1" --b2 "$B2" --od "$TMP_OUT" --tmp "$TMP_OUT/prep_tmp" \
        -t "$READ_TYPE" --libType "$LIB_TYPE" --readLength "$READ_LENGTH" \
        --variable-read-length --anchorLength "$ANCHOR_LENGTH" \
        --nthread "$THREADS_PER_TASK" --novelSS --cstat "$CSTAT" --task prep \
        > /dev/null 2>&1

    echo "[$(date +%T)] [$CELL_LINE|$COMP] Starting Post..."
    rmats.py \
        --gtf "$GTF" --b1 "$B1" --b2 "$B2" --od "$POST_OUT" --tmp "$TMP_OUT/prep_tmp" \
        -t "$READ_TYPE" --libType "$LIB_TYPE" --readLength "$READ_LENGTH" \
        --variable-read-length --anchorLength "$ANCHOR_LENGTH" \
        --nthread "$THREADS_PER_TASK" --novelSS --cstat "$CSTAT" --task post \
        > /dev/null 2>&1

    # Significance Filtering
    for AS_TYPE in SE A5SS A3SS MXE RI; do
        local JCEC="$POST_OUT/${AS_TYPE}.MATS.JCEC.txt"
        local FILTERED="$POST_OUT/${AS_TYPE}.MATS.JCEC.filtered.txt"
        if [ -f "$JCEC" ]; then
            # Determine column indices (MXE has 2 extra exon columns)
            if [ "$AS_TYPE" == "MXE" ]; then
                IJC1=15; SJC1=16; IJC2=17; SJC2=18; FDR_COL=22; DPSI_COL=25
            else
                IJC1=13; SJC1=14; IJC2=15; SJC2=16; FDR_COL=20; DPSI_COL=23
            fi

            awk -F'\t' -v i1=$IJC1 -v s1=$SJC1 -v i2=$IJC2 -v s2=$SJC2 -v fc=$FDR_COL -v dc=$DPSI_COL \
            'NR==1 { print > "'"$FILTERED"'"; next }
            {
                n1=split($i1, ijc1, ","); split($s1, sjc1, ",")
                sum1=0; for(i=1;i<=n1;i++) sum1+=(ijc1[i]+sjc1[i])
                avg1=sum1/n1
                n2=split($i2, ijc2, ","); split($s2, sjc2, ",")
                sum2=0; for(i=1;i<=n2;i++) sum2+=(ijc2[i]+sjc2[i])
                avg2=sum2/n2
                fdr=$fc+0; dpsi=$dc+0
                # Filtering: FDR <= 0.05, |DeltaPSI| >= 0.1, Average Counts >= 10
                if (fdr <= 0.05 && (dpsi >= 0.1 || dpsi <= -0.1) && avg1 >= 10 && avg2 >= 10)
                    print >> "'"$FILTERED"'"
            }' "$JCEC"
        fi
    done
    echo "[$(date +%T)] [$CELL_LINE|$COMP] Complete."

    # --- Clean up & Organize ---
    mkdir -p "$POST_OUT/filtered" "$POST_OUT/raw" "$POST_OUT/intermediate"
    mv "$POST_OUT"/*.filtered.txt "$POST_OUT/filtered/" 2>/dev/null || true
    mv "$POST_OUT"/*.MATS.JCEC.txt "$POST_OUT/raw/" 2>/dev/null || true
    mv "$POST_OUT"/summary.txt "$POST_OUT/raw/" 2>/dev/null || true
    mv "$POST_OUT"/*.MATS.JC.txt "$POST_OUT/intermediate/" 2>/dev/null || true
    mv "$POST_OUT"/*.raw.input.*.txt "$POST_OUT/intermediate/" 2>/dev/null || true
    mv "$POST_OUT"/fromGTF.*.txt "$POST_OUT/intermediate/" 2>/dev/null || true
    [ -d "$POST_OUT/tmp" ] && mv "$POST_OUT/tmp" "$POST_OUT/intermediate/"
}

export -f run_comparison
export READ_TYPE LIB_TYPE READ_LENGTH ANCHOR_LENGTH CSTAT THREADS_PER_TASK BAM_LISTS_DIR TMP_DIR RESULTS_DIR

# ──────────────────────────────────────────────────────────────────────────────
# MAIN LOOP (Background Management)
# ──────────────────────────────────────────────────────────────────────────────

CELL_LINES=("UL1" "CNK89" "H9" "SUPM2")
COMPARISONS=(
    "Romidepsin_6nM_vs_DMSO_Romidepsin"
    "Kromastat_6nM_vs_DMSO_Kromastat"
    "Kromastat_6nM_vs_Romidepsin_6nM"
)

echo ">>> Launching comparisons..."

for CL in "${CELL_LINES[@]}"; do
    # Determine species/GTF
    if [[ "$CL" == "UL1" || "$CL" == "CNK89" ]]; then
        SPECIES="Canine"; GTF="$GTF_CANINE"
    else
        SPECIES="Human"; GTF="$GTF_HUMAN"
    fi

    for COMP in "${COMPARISONS[@]}"; do
        # Manage parallel slots
        while [ $(jobs -r | wc -l) -ge "$PARALLEL_TASKS" ]; do
            sleep 10
        done

        echo ">>> Submitting: $CL | $COMP"
        COMP_LOG="$LOG_DIR/${CL}_${COMP}.log"
        
        run_comparison "$CL" "$COMP" "$SPECIES" "$GTF" > "$COMP_LOG" 2>&1 &
    done
done

echo "Waiting for all rMATS jobs to complete..."
wait

echo "========================================================================"
echo "   rMATS-turbo pipeline COMPLETE"
echo "========================================================================"
