# [[scripts_upstream/03_alignment_qc.sh]]
# Runs Picard CollectRnaSeqMetrics on aligned BAM files.

# Exit on error
set -e
set -o pipefail

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
BAM_DIR="$BASE_DIR/../_data/bam"
GENOME_DIR="$BASE_DIR/../_data/genome"
QC_DIR="$BASE_DIR/../_data/qc"
LOG_DIR="$BASE_DIR/../_data/logs/qc"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../drPhuong_Sample_Data_Table.csv"

mkdir -p "$QC_DIR"

# --- Resource Auto-detection ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-32}
PARALLEL_SAMPLES=4
CPUS_PER_SAMPLE=$((TOTAL_CPUS / PARALLEL_SAMPLES))
[ "$CPUS_PER_SAMPLE" -lt 1 ] && CPUS_PER_SAMPLE=1

# Memory management for Picard (Java)
JAVA_MEM_PER_SAMPLE=$((128 / PARALLEL_SAMPLES - 4))
[ "$JAVA_MEM_PER_SAMPLE" -lt 4 ] && JAVA_MEM_PER_SAMPLE=4

echo "=== Starting Alignment QC (Picard) ==="
echo "Resources: $PARALLEL_SAMPLES concurrent samples, ${JAVA_MEM_PER_SAMPLE}G RAM per Picard job."

# Function to run QC for a single sample
run_qc() {
    local bam_file=$1
    local sample_name=$(basename "$(dirname "$bam_file")")
    local SAMPLE_QC_DIR="$QC_DIR/$sample_name"
    local OUTPUT_METRICS="$SAMPLE_QC_DIR/${sample_name}.picard_rnaseq_metrics"

    if [ -f "$OUTPUT_METRICS" ]; then
        echo "[Sample: $sample_name] Picard metrics already exist. Skipping."
        return 0
    fi

    mkdir -p "$SAMPLE_QC_DIR"

    echo "[Sample: $sample_name] Detecting species..."
    local sample_species=$(python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | grep "^$sample_name " | awk '{print $4}')
    if [ -z "$sample_species" ]; then
        echo "Error: Could not determine species for $sample_name."
        return 1
    fi

    # Select Reference files
    if [ "$sample_species" == "Human" ]; then
        REFFLAT="$GENOME_DIR/Homo_sapiens.GRCh38.113.refFlat"
    else
        REFFLAT="$GENOME_DIR/Canis_lupus_familiaris.ROS_Cfam_1.0.113.refFlat"
    fi

    if [ ! -f "$REFFLAT" ]; then
        echo "Error: RefFlat file not found at $REFFLAT for $sample_name."
        return 1
    fi

    # Picard CollectRnaSeqMetrics
    echo "[Sample: $sample_name] Running Picard CollectRnaSeqMetrics..."
    picard -Xmx${JAVA_MEM_PER_SAMPLE}g CollectRnaSeqMetrics \
        I="$bam_file" \
        O="$SAMPLE_QC_DIR/${sample_name}.metrics.tmp" \
        REF_FLAT="$REFFLAT" \
        STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
        RIBOSOMAL_INTERVALS=null \
        CHART_OUTPUT="$SAMPLE_QC_DIR/${sample_name}.coverage.pdf.tmp"

    mv "$SAMPLE_QC_DIR/${sample_name}.metrics.tmp" "$OUTPUT_METRICS"
    mv "$SAMPLE_QC_DIR/${sample_name}.coverage.pdf.tmp" "$SAMPLE_QC_DIR/${sample_name}.picard_coverage.pdf"

    echo "[Sample: $sample_name] Done."
}

export -f run_qc
export QC_DIR BAM_DIR GENOME_DIR UTILS_DIR CSV_FILE JAVA_MEM_PER_SAMPLE

# Simple background job management
BAM_LIST=$(find "$BAM_DIR" -maxdepth 2 -name "*_Aligned.sortedByCoord.out.bam")
SAMPLE_COUNT=$(echo "$BAM_LIST" | wc -l)
CURRENT_SAMPLE=0

for bam in $BAM_LIST; do
    sample_name=$(basename "$(dirname "$bam")")
    CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
    
    while [ $(jobs -r | wc -l) -ge "$PARALLEL_SAMPLES" ]; do
        sleep 5
    done

    mkdir -p "$LOG_DIR"
    echo ">>> Submitting [$CURRENT_SAMPLE/$SAMPLE_COUNT]: $sample_name"
    SAMPLE_LOG="$LOG_DIR/${sample_name}.picard_qc.log"
    
    bash -c 'run_qc "'"$bam"'"' > "$SAMPLE_LOG" 2>&1 &
done

echo "Waiting for all Picard jobs to complete..."
wait

echo "=== Alignment QC Complete ==="
