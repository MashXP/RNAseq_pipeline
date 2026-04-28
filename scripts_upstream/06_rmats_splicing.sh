#!/bin/bash

# [[scripts_upstream/06_rmats_splicing.sh]]
# rMATS-turbo Alternative Splicing Pipeline
#
# Workflow:
#   Step 0  — Generate b1.txt / b2.txt BAM lists from metadata CSV
#   Step 1  — rMATS Prep step per comparison (per-sample splicing graph)
#   Step 2  — rMATS Post step per comparison (event detection + statistics)
#
# Post-hoc significance filters applied to output tables (not rMATS CLI args):
#   FDR <= 0.01  |  |ΔPSI| >= 0.05  |  avg read coverage >= 10
#
# Library type: fr-secondstrand (dUTP stranded, standard for Illumina TruSeq)
# Read type:    paired-end
# Novel splice sites: enabled (--novelSS)

set -e
set -o pipefail

# ──────────────────────────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────────────────────────

BASE_DIR=$(dirname "$(realpath "$0")")
PIPELINE_ROOT="$BASE_DIR/.."
UTILS_DIR="$BASE_DIR/utils"

# Genome / annotation
GENOME_DIR="$PIPELINE_ROOT/_data/genome"

# GTF files (unzipped by 01_genome_prep.sh)
GTF_HUMAN="$GENOME_DIR/Human/Homo_sapiens.GRCh38.113.gtf"
GTF_CANINE="$GENOME_DIR/Canine/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf"

# rMATS-turbo output under _data/rMATS-turbo/
BAM_LISTS_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/bam_lists"
RESULTS_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/results"
TMP_DIR="$PIPELINE_ROOT/_data/rMATS-turbo/tmp"


# ──────────────────────────────────────────────────────────────────────────────
# PARAMETERS
# ──────────────────────────────────────────────────────────────────────────────

READ_LENGTH=150                 # Set to your actual sequenced read length
LIB_TYPE="fr-secondstrand"      # dUTP stranded (Illumina TruSeq Stranded)
READ_TYPE="paired"              # Paired-end
CSTAT=0.05                      # Null hypothesis cutoff for |ΔPSI| (>= 0.05)
ANCHOR_LENGTH=1                 # Minimum overhang for splice junction reads

# Resource detection
THREADS=${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 8)}

# Cell line → species mapping
declare -A CELL_LINE_SPECIES
CELL_LINE_SPECIES["UL1"]="Canine"
CELL_LINE_SPECIES["CNK89"]="Canine"
CELL_LINE_SPECIES["H9"]="Human"
CELL_LINE_SPECIES["SUPM2"]="Human"

# The 3 contrasts per cell line (matches generate_bam_lists.py)
COMPARISONS=(
    "Romidepsin_6nM_vs_DMSO_Romidepsin"
    "Kromastat_6nM_vs_DMSO_Kromastat"
    "Kromastat_6nM_vs_Romidepsin_6nM"
)

# Build a flat list of all (CellLine, Comparison) pairs for parallelization
TASK_LIST=()
CELL_LINES=("UL1" "CNK89" "H9" "SUPM2")
for CL in "${CELL_LINES[@]}"; do
    for CP in "${COMPARISONS[@]}"; do
        TASK_LIST+=("$CL|$CP")
    done
done

# If TASK_INDEX is provided (e.g. from Slurm Array), only run that specific task
# Otherwise, loop through all (original behavior)
if [ -n "$1" ]; then
    TASK_INDEX=$1
elif [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    TASK_INDEX=$SLURM_ARRAY_TASK_ID
fi

if [ -n "$TASK_INDEX" ]; then
    IFS='|' read -r TARGET_CELL_LINE TARGET_COMP <<< "${TASK_LIST[$TASK_INDEX]}"
    echo "[INFO] Running as PARALLEL TASK #$TASK_INDEX: $TARGET_CELL_LINE / $TARGET_COMP"
fi

# ──────────────────────────────────────────────────────────────────────────────
# LOGGING
# ──────────────────────────────────────────────────────────────────────────────

LOG_DIR="$PIPELINE_ROOT/_data/logs"
mkdir -p "$LOG_DIR"

# Use Task ID in log name if running in parallel to avoid collisions
LOG_SUFFIX="${TASK_INDEX:-$(date +%H%M%S)}"
LOG_FILE="$LOG_DIR/rmats_${LOG_SUFFIX}.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "Log file: $LOG_FILE"

# ──────────────────────────────────────────────────────────────────────────────
# PRE-FLIGHT CHECKS
# ──────────────────────────────────────────────────────────────────────────────

echo ""
echo "========================================================================"
echo "   rMATS-turbo Alternative Splicing Pipeline"
echo "========================================================================"
echo "Resources : $THREADS threads"
echo "Read type : $READ_TYPE, length = $READ_LENGTH bp"
echo "Library   : $LIB_TYPE"
echo "CSTAT     : $CSTAT  (|ΔΨPSI| null hypothesis cutoff)"
echo "Post-hoc  : FDR <= 0.01  |  |ΔPSI| >= 0.05  |  avg reads >= 10"
echo "------------------------------------------------------------------------"

command -v rmats.py >/dev/null 2>&1 || \
    command -v python3 >/dev/null 2>&1 && python3 -c "import rmats" 2>/dev/null || \
    { echo "[ERROR] rmats.py not found. Is the conda environment active?"; exit 1; }

if [ ! -f "$GTF_HUMAN" ]; then
    echo "[ERROR] Human GTF not found: $GTF_HUMAN"
    echo "        Run scripts_upstream/01_genome_prep.sh first."
    exit 1
fi

if [ ! -f "$GTF_CANINE" ]; then
    echo "[ERROR] Canine GTF not found: $GTF_CANINE"
    echo "        Run scripts_upstream/01_genome_prep.sh first."
    exit 1
fi

mkdir -p "$RESULTS_DIR" "$TMP_DIR"

# ──────────────────────────────────────────────────────────────────────────────
# STEP 0 — Generate b1.txt / b2.txt BAM lists
# ──────────────────────────────────────────────────────────────────────────────

echo ""
echo "========================================================================"
echo "   STEP 0 — Generating BAM list files from metadata CSV"
echo "========================================================================"

python3 "$UTILS_DIR/generate_bam_lists.py"
echo "[INFO] BAM lists ready under: $BAM_LISTS_DIR"

# ──────────────────────────────────────────────────────────────────────────────
# STEP 1 — rMATS Prep step (per comparison)
# ──────────────────────────────────────────────────────────────────────────────

echo ""
echo "========================================================================"
echo "   STEP 1 — rMATS Prep (per-sample splicing graph extraction)"
echo "========================================================================"

for CELL_LINE in "${!CELL_LINE_SPECIES[@]}"; do
    # Skip if we are targeting a different cell line in parallel mode
    if [ -n "$TARGET_CELL_LINE" ] && [ "$CELL_LINE" != "$TARGET_CELL_LINE" ]; then continue; fi

    SPECIES="${CELL_LINE_SPECIES[$CELL_LINE]}"

    # Select species-correct GTF
    if [ "$SPECIES" = "Human" ]; then
        GTF="$GTF_HUMAN"
    else
        GTF="$GTF_CANINE"
    fi

    for COMP in "${COMPARISONS[@]}"; do
        # Skip if we are targeting a different comparison in parallel mode
        if [ -n "$TARGET_COMP" ] && [ "$COMP" != "$TARGET_COMP" ]; then continue; fi
        COMP_BAM_DIR="$BAM_LISTS_DIR/$CELL_LINE/$COMP"
        B1="$COMP_BAM_DIR/b1.txt"
        B2="$COMP_BAM_DIR/b2.txt"

        # Skip if BAM lists are missing (e.g., head-to-head if desired to skip)
        if [ ! -f "$B1" ] || [ ! -f "$B2" ]; then
            echo "[WARN] BAM list not found for $CELL_LINE / $COMP — skipping."
            continue
        fi

        TMP_OUT="$TMP_DIR/$CELL_LINE/$COMP"
        mkdir -p "$TMP_OUT"

        echo "------------------------------------------------------------------------"
        echo "[PREP] $CELL_LINE | $COMP | Species: $SPECIES"
        echo "       GTF   : $GTF"
        echo "       b1    : $B1"
        echo "       b2    : $B2"
        echo "       tmp   : $TMP_OUT"

        rmats.py \
            --gtf     "$GTF" \
            --b1      "$B1" \
            --b2      "$B2" \
            --od      "$TMP_OUT" \
            --tmp     "$TMP_OUT/prep_tmp" \
            -t        "$READ_TYPE" \
            --libType "$LIB_TYPE" \
            --readLength "$READ_LENGTH" \
            --variable-read-length \
            --anchorLength "$ANCHOR_LENGTH" \
            --nthread "$THREADS" \
            --novelSS \
            --cstat   "$CSTAT" \
            --task    prep \
        || { echo "[ERROR] rMATS prep failed for $CELL_LINE / $COMP. Aborting."; exit 1; }

        echo "[DONE ] Prep complete: $CELL_LINE / $COMP"
    done
done

# ──────────────────────────────────────────────────────────────────────────────
# STEP 2 — rMATS Post step (per comparison)
# ──────────────────────────────────────────────────────────────────────────────

echo ""
echo "========================================================================"
echo "   STEP 2 — rMATS Post (event detection + differential splicing stats)"
echo "========================================================================"

for CELL_LINE in "${!CELL_LINE_SPECIES[@]}"; do
    # Skip if we are targeting a different cell line in parallel mode
    if [ -n "$TARGET_CELL_LINE" ] && [ "$CELL_LINE" != "$TARGET_CELL_LINE" ]; then continue; fi

    SPECIES="${CELL_LINE_SPECIES[$CELL_LINE]}"

    if [ "$SPECIES" = "Human" ]; then
        GTF="$GTF_HUMAN"
    else
        GTF="$GTF_CANINE"
    fi

    for COMP in "${COMPARISONS[@]}"; do
        # Skip if we are targeting a different comparison in parallel mode
        if [ -n "$TARGET_COMP" ] && [ "$COMP" != "$TARGET_COMP" ]; then continue; fi
        COMP_BAM_DIR="$BAM_LISTS_DIR/$CELL_LINE/$COMP"
        B1="$COMP_BAM_DIR/b1.txt"
        B2="$COMP_BAM_DIR/b2.txt"

        if [ ! -f "$B1" ] || [ ! -f "$B2" ]; then
            echo "[WARN] Skipping post for $CELL_LINE / $COMP — no BAM list."
            continue
        fi

        TMP_OUT="$TMP_DIR/$CELL_LINE/$COMP"
        POST_OUT="$RESULTS_DIR/$CELL_LINE/$COMP"
        mkdir -p "$POST_OUT"

        echo "------------------------------------------------------------------------"
        echo "[POST] $CELL_LINE | $COMP | Species: $SPECIES"
        echo "       GTF    : $GTF"
        echo "       tmp    : $TMP_OUT"
        echo "       output : $POST_OUT"

        rmats.py \
            --gtf     "$GTF" \
            --b1      "$B1" \
            --b2      "$B2" \
            --od      "$POST_OUT" \
            --tmp     "$TMP_OUT/prep_tmp" \
            -t        "$READ_TYPE" \
            --libType "$LIB_TYPE" \
            --readLength "$READ_LENGTH" \
            --variable-read-length \
            --anchorLength "$ANCHOR_LENGTH" \
            --nthread "$THREADS" \
            --novelSS \
            --cstat   "$CSTAT" \
            --task    post \
        || { echo "[ERROR] rMATS post failed for $CELL_LINE / $COMP. Aborting."; exit 1; }

        echo "[DONE ] Post complete: $CELL_LINE / $COMP"

        # ── Post-hoc significance filter on the JCEC output (SE events as example)
        # Filters: FDR <= 0.01 | |ΔPSI| >= 0.05 | avg(IJC+SJC) >= 10 per group
        echo "[INFO ] Applying post-hoc filters to JCEC output files..."
        for AS_TYPE in SE A5SS A3SS MXE RI; do
            JCEC="$POST_OUT/${AS_TYPE}.MATS.JCEC.txt"
            FILTERED="$POST_OUT/${AS_TYPE}.MATS.JCEC.filtered.txt"
            if [ -f "$JCEC" ]; then
                # awk filter:
                #   col 19 = FDR, col 22 = IncLevelDifference (ΔPSI)
                #   cols 7,8 = IJC_SAMPLE_1, SJC_SAMPLE_1 (comma-separated replicates)
                #   cols 9,10 = IJC_SAMPLE_2, SJC_SAMPLE_2
                awk -F'\t' 'NR==1 { print > "'"$FILTERED"'"; next }
                {
                    # Average read count group 1: mean of (IJC+SJC) across replicates
                    n1=split($7, ijc1, ","); split($8, sjc1, ",")
                    sum1=0; for(i=1;i<=n1;i++) sum1+=(ijc1[i]+sjc1[i])
                    avg1=sum1/n1

                    # Average read count group 2
                    n2=split($9, ijc2, ","); split($10, sjc2, ",")
                    sum2=0; for(i=1;i<=n2;i++) sum2+=(ijc2[i]+sjc2[i])
                    avg2=sum2/n2

                    fdr=$19+0; dpsi=$22+0
                    if (fdr <= 0.01 && (dpsi >= 0.05 || dpsi <= -0.05) && avg1 >= 10 && avg2 >= 10)
                        print >> "'"$FILTERED"'"
                }' "$JCEC"
                echo "        ${AS_TYPE}: $(wc -l < "$FILTERED") significant events → $FILTERED"
            fi
        done
    done
done

echo ""
echo "========================================================================"
echo "   rMATS-turbo pipeline COMPLETE"
echo "========================================================================"
echo "  Results  : $RESULTS_DIR"
echo "  Temp     : $TMP_DIR"
echo "  Log      : $LOG_FILE"
echo "========================================================================"
