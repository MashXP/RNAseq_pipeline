#!/bin/bash

# [[test_upstream/01_genome_prep.sh]]
# TEST VERSION: Prepares the Chromosome 21 test environment.
# Note: This is an automated "mini" version of the production prep script.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
GENOME_DIR="$BASE_DIR/../_data/genome_test"
INDEX_DIR="$BASE_DIR/../_data/index_test"
FASTQ_DIR="$BASE_DIR/../_data/fastq"
TRIM_DIR="$BASE_DIR/../_data/fastq_trimmed_test"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../MiaPAca-2_Sample_Data_Table.csv"

# 1. Ensembl 113 URLs (Chr21 subset)
FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"

# 2. Downloading Genome Files
echo "=== Step 1: Downloading Genome Files (Chr21) ==="
mkdir -p "$GENOME_DIR"

if [ ! -f "$GENOME_DIR/chr21.fa.gz" ]; then
    echo "Downloading Chr21 FASTA..."
    wget -O "$GENOME_DIR/chr21.fa.gz" "$FASTA_URL"
else
    echo "Chr21 FASTA already exists. Skipping."
fi

if [ ! -f "$GENOME_DIR/human_eb113.gtf.gz" ]; then
    echo "Downloading Human GTF (Ensembl 113)..."
    wget -O "$GENOME_DIR/human_eb113.gtf.gz" "$GTF_URL"
else
    echo "Human GTF already exists. Skipping."
fi

# 3. Extracting and Filtering Files
echo "=== Step 2: Extracting and Filtering Files ==="
if [ ! -f "$GENOME_DIR/chr21.fa" ]; then
    echo "Extracting FASTA..."
    gunzip -c "$GENOME_DIR/chr21.fa.gz" > "$GENOME_DIR/chr21.fa"
fi

if [ ! -f "$GENOME_DIR/chr21_annotations.gtf" ]; then
    echo "Extracting and filtering GTF for Chromosome 21..."
    # Filter for Chr21 to speed up test indexing
    gunzip -c "$GENOME_DIR/human_eb113.gtf.gz" | grep "^21\|^#" > "$GENOME_DIR/chr21_annotations.gtf"
fi

# 4. Generating BED12 for RSeQC
echo "=== Step 3: Generating BED12 for RSeQC ==="
if [ ! -f "$GENOME_DIR/chr21.bed" ]; then
    python3 "$UTILS_DIR/gtf2bed.py" "$GENOME_DIR/chr21_annotations.gtf" > "$GENOME_DIR/chr21.bed"
fi

# 5. Creating STAR Index
echo "=== Step 4: Creating STAR Index (Chr21) ==="
mkdir -p "$INDEX_DIR"
if [ -z "$(ls -A "$INDEX_DIR")" ]; then
    echo "Running STAR genomeGenerate..."
    STAR --runThreadN 4 \
         --runMode genomeGenerate \
         --genomeDir "$INDEX_DIR" \
         --genomeFastaFiles "$GENOME_DIR/chr21.fa" \
         --sjdbGTFfile "$GENOME_DIR/chr21_annotations.gtf" \
         --sjdbOverhang 99 \
         --genomeSAindexNbases 11
else
    echo "STAR index already exists. Skipping."
fi

# 6. Trimming with Trimmomatic
echo "=== Step 5: Trimming Reads ==="
mkdir -p "$TRIM_DIR"

# 6.1 Defensive Check: Verify FASTQ presence before starting
echo "Checking for required FASTQ files..."
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2
do
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "--------------------------------------------------------------------------------"
        echo "CRITICAL ERROR: Required test FASTQ files are missing for $sample."
        echo "ACTION REQUIRED: Please run 'subsample_hpc.sh' on the HPC and download"
        echo "                 the results into $FASTQ_DIR before proceeding."
        echo "--------------------------------------------------------------------------------"
        exit 1
    fi
done || exit 1

# Locate adapters — check active conda env first, then search conda envs
ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"
if [ ! -f "$ADAPTERS" ]; then
    ADAPTERS=$(find "$CONDA_PREFIX/share" -name "TruSeq3-PE.fa" 2>/dev/null | head -n 1)
fi

if [ -z "$ADAPTERS" ]; then
    echo "ERROR: Could not find TruSeq3-PE.fa adapter file. Cannot proceed with trimming."
    echo "Ensure Trimmomatic adapters are installed (check conda: trimmomatic/adapters/)."
    exit 1
fi

python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2
do
    echo "-----------------------------------------------------"
    echo "Trimming Sample: $sample"
    
    # Check if raw FASTQ exists
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "Error: FASTQ files not found in $FASTQ_DIR:"
        echo "  R1: $r1"
        echo "  R2: $r2"
        echo "Skipping sample $sample."
        continue
    fi

    if [ -f "$TRIM_DIR/$r1" ] && [ -f "$TRIM_DIR/$r2" ]; then
        echo "Trimmed files already exist for $sample. Skipping."
        continue
    fi

    echo "Running Trimmomatic..."
    trimmomatic PE -threads 4 \
        "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
        "$TRIM_DIR/$r1" "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" \
        "$TRIM_DIR/$r2" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || { echo "ERROR: Trimmomatic failed for $sample. Aborting."; exit 1; }
    
    rm "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz"
done

echo "[TEST] Upstream preprocessing complete."
