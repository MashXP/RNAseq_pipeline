#!/bin/bash
#SBATCH --job-name=QC_Regen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=hpc_qc_regen_%j.log
#SBATCH --error=hpc_qc_regen_%j.err

# [[_hpc/hpc_qc_only.sbatch]]
# Temporary script to (re)generate FastQC and MultiQC reports.
# Assumes alignment, quantification, and Picard QC are already complete.

# 1. Initialize Micromamba (Safe method for batch shells)
eval "$(micromamba shell hook --shell bash)"
micromamba activate cancer_rnaseq

# 2. Setup Directories
# Assumes submission from project root
BASE_DIR=$(pwd)
FASTQ_DIR="$BASE_DIR/_data/fastq"
RAW_QC_DIR="$BASE_DIR/_data/qc_raw"
CSV_FILE="$BASE_DIR/drPhuong_Sample_Data_Table.csv"
UTILS_DIR="$BASE_DIR/scripts_upstream/utils"

mkdir -p "$RAW_QC_DIR"

echo "================================================================================"
echo "   QC REGENERATION: Starting FastQC and MultiQC"
echo "================================================================================"
echo "Project Root: $BASE_DIR"

# --- Resource Auto-detection ---
THREADS=${SLURM_CPUS_PER_TASK:-16}
echo "Using $THREADS threads for FastQC."

# 3. Run FastQC on Raw Reads
echo ""
echo ">>> Step 1: Running FastQC on Raw Reads..."
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2 species group
do
    echo "-----------------------------------------------------"
    echo "Processing Sample: $sample"
    
    # Check if raw FASTQ exists
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "WARNING: FASTQ files not found in $FASTQ_DIR for $sample. Skipping."
        continue
    fi

    # Check if reports already exist
    if [ ! -f "$RAW_QC_DIR/${r1%.fastq.gz}_fastqc.html" ] || [ ! -f "$RAW_QC_DIR/${r2%.fastq.gz}_fastqc.html" ]; then
        echo "Generating FastQC reports..."
        fastqc -t "$THREADS" "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" -o "$RAW_QC_DIR"
    else
        echo "FastQC results for $sample already exist. (Delete $RAW_QC_DIR to force regeneration)."
    fi
done

# 4. Run MultiQC
echo ""
echo ">>> Step 2: Running MultiQC Summary..."
# Calling the upstream script ensures all paths (BAM, QC, Counts) are aggregated correctly
bash scripts_upstream/05_multiqc.sh

echo ""
echo "================================================================================"
echo "   QC REGENERATION COMPLETE"
echo "================================================================================"
echo "MultiQC Report: _data/multiqc/rna_seq_pipeline_summary.html"
