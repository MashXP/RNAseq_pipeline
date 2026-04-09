#!/bin/bash

# [[test_upstream/01_genome_prep.sh]]
# TEST VERSION: Prepares the Chromosome 21 test environment.
# Logic: 
#   1. Verifies Chr21 FASTA/GTF presence.
#   2. Trims subsampled reads using Trimmomatic (matching production parameters).
# Divergence from Prod: Uses Chr21 only; outputs to _test directories.

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
THREADS=4 # Uniform for test; Production uses 8
GENOME_DIR="$BASE_DIR/../_data/genome_test"
INDEX_DIR="$BASE_DIR/../_data/index_test"

echo "=== [TEST] Step 1: Checking Chr21 Reference Files ==="
if [ ! -f "$GENOME_DIR/chr21.fa" ] || [ ! -f "$GENOME_DIR/chr21_annotations.gtf" ]; then
    echo "Error: Test genome files not found in $GENOME_DIR. Please run ./verify_pipeline.sh first."
    exit 1
fi

echo "Test genome files found."

echo "=== [TEST] Step 2: Checking STAR Index ==="
if [ -z "$(ls -A "$INDEX_DIR")" ]; then
    echo "STAR index not found in $INDEX_DIR. Building it now..."
    STAR --runThreadN 4 \
         --runMode genomeGenerate \
         --genomeDir "$INDEX_DIR" \
         --genomeFastaFiles "$GENOME_DIR/chr21.fa" \
         --sjdbGTFfile "$GENOME_DIR/chr21_annotations.gtf" \
         --sjdbOverhang 99 \
         --genomeSAindexNbases 11
else
    echo "STAR index already exists. Skipping build."
fi

echo "[TEST] Step 3: Trimming Reads (Trimmomatic) ==="
TRIM_DIR="$BASE_DIR/../_data/fastq_trimmed_test"
FASTQ_DIR="$BASE_DIR/../_data/fastq"
CSV_FILE="$BASE_DIR/../MiaPAca-2_Sample_Data_Table.csv"
UTILS_DIR="$BASE_DIR/utils"

mkdir -p "$TRIM_DIR"

# Locate adapters
ADAPTERS="/opt/conda/share/trimmomatic/adapters/TruSeq3-PE.fa"
if [ ! -f "$ADAPTERS" ]; then
    ADAPTERS=$(find / -name "TruSeq3-PE.fa" 2>/dev/null | head -n 1)
fi

python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2
do
    echo "Trimming Test Sample: $sample"
    if [ ! -f "$FASTQ_DIR/$r1" ]; then continue; fi

    trimmomatic PE -threads 4 \
        "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
        "$TRIM_DIR/$r1" "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" \
        "$TRIM_DIR/$r2" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    rm "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz"
done

echo "[TEST] Upstream prep complete."
