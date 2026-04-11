#!/bin/bash

# [[scripts_upstream/01_genome_prep.sh]]
# This script downloads the human genome (GRCh38) and creates a STAR index.

# Exit on error
set -e

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
GENOME_DIR="$BASE_DIR/../_data/genome"
INDEX_DIR="$BASE_DIR/../_data/index"
FASTQ_DIR="$BASE_DIR/../_data/fastq"
TRIM_DIR="$BASE_DIR/../_data/fastq_trimmed"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../MiaPAca-2_Sample_Data_Table.csv"

# 1. Ensembl 113 URLs (based on CSV info)
FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"

# 2. Downloading Genome Files
echo "=== Step 1: Downloading Genome Files ==="
mkdir -p "$GENOME_DIR"

if [ ! -f "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]; then
    echo "Downloading FASTA..."
    wget -P "$GENOME_DIR" "$FASTA_URL"
else
    echo "FASTA already exists. Skipping."
fi

if [ ! -f "$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf.gz" ]; then
    echo "Downloading GTF..."
    wget -P "$GENOME_DIR" "$GTF_URL"
else
    echo "GTF already exists. Skipping."
fi

# 3. Extracting Files
echo "=== Step 2: Extracting Files ==="
# Decompressing (required for STAR and gtf2bed)
if [ ! -f "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    gunzip -c "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" > "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fi

if [ ! -f "$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf" ]; then
    gunzip -c "$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf.gz" > "$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf"
fi

# 4. Generating BED12 for RSeQC
echo "=== Step 3: Generating BED12 for RSeQC ==="
if [ ! -f "$GENOME_DIR/Homo_sapiens.GRCh38.113.bed" ]; then
    python3 "$UTILS_DIR/gtf2bed.py" "$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf" > "$GENOME_DIR/Homo_sapiens.GRCh38.113.bed"
fi

# 5. Creating STAR Index
echo "=== Step 4: Creating STAR Index ==="
mkdir -p "$INDEX_DIR"
if [ -z "$(ls -A "$INDEX_DIR")" ]; then
    echo "Running STAR genomeGenerate..."
    STAR --runThreadN 32 \
         --runMode genomeGenerate \
         --genomeDir "$INDEX_DIR" \
         --genomeFastaFiles "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
         --sjdbGTFfile "$GENOME_DIR/Homo_sapiens.GRCh38.113.gtf" \
         --sjdbOverhang 99
else
    echo "STAR index already exists. Skipping."
fi

# 6. Trimming with Trimmomatic
echo "=== Step 5: Trimming Reads ==="
mkdir -p "$TRIM_DIR"

# --- Resource Auto-detection ---
# Default to 32 threads or respect SLURM
THREADS=${SLURM_CPUS_PER_TASK:-32}

# Detect RAM and set Java heap size (~90% of available)
if [ -n "$SLURM_MEM_PER_NODE" ]; then
    MEM_GB=$(($SLURM_MEM_PER_NODE / 1024))
else
    # Fallback to system check or default to 128
    MEM_GB=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}' || echo 128)
fi

# Allocate ~90% of RAM to Java, leaving some for OS/threading overhead
JAVA_MEM=$(($MEM_GB - 12))
[ "$JAVA_MEM" -lt 4 ] && JAVA_MEM=4 # Safety minimum

echo "Resources: $THREADS threads, allocated ${JAVA_MEM}G Java Heap."
# -------------------------------

# 6.1 Defensive Check: Verify FASTQ presence before starting
echo "Checking for required FASTQ files..."
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2
do
    if [ ! -f "$FASTQ_DIR/$r1" ] || [ ! -f "$FASTQ_DIR/$r2" ]; then
        echo "--------------------------------------------------------------------------------"
        echo "CRITICAL ERROR: Required FASTQ files are missing for $sample."
        echo "ACTION REQUIRED: Please ensure all FASTQ files listed in the CSV are present"
        echo "                 in $FASTQ_DIR before proceeding."
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

    echo "Running Trimmomatic (Memory: ${JAVA_MEM}G, Threads: ${THREADS})..."
    trimmomatic -Xmx${JAVA_MEM}g PE -threads "$THREADS" \
        "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
        "$TRIM_DIR/$r1" "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" \
        "$TRIM_DIR/$r2" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || { echo "ERROR: Trimmomatic failed for $sample. Aborting."; exit 1; }

    rm "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz"
done

echo "Upstream preprocessing complete."
