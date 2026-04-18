#!/bin/bash

# [[scripts_upstream/01_genome_prep.sh]]
# This script downloads the human genome (GRCh38) and creates a STAR index.

# Exit on error
set -e
set -o pipefail # Catch errors in piped loops

# Directories
BASE_DIR=$(dirname "$(realpath "$0")")
GENOME_DIR="$BASE_DIR/../_data/genome"
INDEX_DIR="$BASE_DIR/../_data/index"
FASTQ_DIR="$BASE_DIR/../_data/fastq"
TRIM_DIR="$BASE_DIR/../_data/fastq_trimmed"
UTILS_DIR="$BASE_DIR/utils"
CSV_FILE="$BASE_DIR/../drPhuong_Sample_Data_Table.csv"

# 1. Ensembl 113 URLs
declare -A FASTA_URLS
declare -A GTF_URLS
FASTA_URLS[Human]="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URLS[Human]="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"
FASTA_URLS[Dog]="https://ftp.ensembl.org/pub/release-113/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz"
GTF_URLS[Dog]="https://ftp.ensembl.org/pub/release-113/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf.gz"

# 2. Downloading Genome Files
echo "================================================================================"
echo "   GENOME PREPARATION: Downloading, Extracting, and Indexing"
echo "================================================================================"
mkdir -p "$GENOME_DIR"

for species in Human Dog; do
    echo "Processing Genome for: $species"
    
    # Organize genome files into species subdirectories
    SPECIES_GENOME_DIR="$GENOME_DIR/$species"
    mkdir -p "$SPECIES_GENOME_DIR"
    
    fasta_url=${FASTA_URLS[$species]}
    gtf_url=${GTF_URLS[$species]}
    fasta_gz=$(basename "$fasta_url")
    gtf_gz=$(basename "$gtf_url")
    fasta_unzipped="${fasta_gz%.gz}"
    gtf_unzipped="${gtf_gz%.gz}"
    
    # Download
    if [ ! -f "$SPECIES_GENOME_DIR/$fasta_gz" ]; then
        echo "Downloading FASTA for $species..."
        wget -P "$SPECIES_GENOME_DIR" "$fasta_url"
    fi
    if [ ! -f "$SPECIES_GENOME_DIR/$gtf_gz" ]; then
        echo "Downloading GTF for $species..."
        wget -P "$SPECIES_GENOME_DIR" "$gtf_url"
    fi
    
    # Extract
    if [ ! -f "$SPECIES_GENOME_DIR/$fasta_unzipped" ]; then gunzip -c "$SPECIES_GENOME_DIR/$fasta_gz" > "$SPECIES_GENOME_DIR/$fasta_unzipped"; fi
    if [ ! -f "$SPECIES_GENOME_DIR/$gtf_unzipped" ]; then gunzip -c "$SPECIES_GENOME_DIR/$gtf_gz" > "$SPECIES_GENOME_DIR/$gtf_unzipped"; fi
    
    # Picard RefFlat (required for CollectRnaSeqMetrics)
    refflat_file="${gtf_unzipped%.gtf}.refFlat"
    if [ ! -f "$SPECIES_GENOME_DIR/$refflat_file" ]; then
        echo "Generating Picard refFlat for $species..."
        gtfToGenePred -genePredExt -geneNameAsName2 "$SPECIES_GENOME_DIR/$gtf_unzipped" "$SPECIES_GENOME_DIR/refflat.tmp.txt"
        awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' "$SPECIES_GENOME_DIR/refflat.tmp.txt" > "$SPECIES_GENOME_DIR/$refflat_file"
        rm -f "$SPECIES_GENOME_DIR/refflat.tmp.txt"
    fi
    
    # STAR Index
    species_index="$INDEX_DIR/$species"
    mkdir -p "$species_index"
    if [ -z "$(ls -A "$species_index")" ]; then
        echo "Running STAR genomeGenerate for $species..."
        STAR --runThreadN 32 \
             --runMode genomeGenerate \
             --genomeDir "$species_index" \
             --genomeFastaFiles "$SPECIES_GENOME_DIR/$fasta_unzipped" \
             --sjdbGTFfile "$SPECIES_GENOME_DIR/$gtf_unzipped" \
             --sjdbOverhang 99
    else
        echo "STAR index for $species already exists. Skipping."
    fi
done

# 3. FastQC on Raw Reads
echo ""
echo "================================================================================"
echo "   QUALITY CONTROL: Running FastQC on Raw Reads"
echo "================================================================================"
RAW_QC_DIR="$BASE_DIR/../_data/qc_raw"
mkdir -p "$RAW_QC_DIR"

# --- Resource Auto-detection ---
# Default to 32 threads or respect SLURM
THREADS=${SLURM_CPUS_PER_TASK:-32}

python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2 species group
do
    echo "Running FastQC for Sample: $sample"
    # Skip if both FastQC reports exist
    if [ ! -f "$RAW_QC_DIR/${r1%.fastq.gz}_fastqc.html" ] || [ ! -f "$RAW_QC_DIR/${r2%.fastq.gz}_fastqc.html" ]; then
        fastqc -t "$THREADS" "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" -o "$RAW_QC_DIR"
    else
        echo "FastQC results for $sample already exist. Skipping."
    fi
done

# 4. Trimming with Trimmomatic
echo ""
echo "================================================================================"
echo "   PREPROCESSING: Trimming Reads with Trimmomatic"
echo "================================================================================"
mkdir -p "$TRIM_DIR"

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
python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2 species group
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

python3 "$UTILS_DIR/parse_samples.py" "$CSV_FILE" | while read -r sample r1 r2 species group
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

    # Integrity Check: skip only if files exist AND are valid GZIPs
    if [ -f "$TRIM_DIR/$r1" ] && [ -f "$TRIM_DIR/$r2" ]; then
        echo "Checking integrity of existing trimmed files for $sample..."
        if gzip -t "$TRIM_DIR/$r1" 2>/dev/null && gzip -t "$TRIM_DIR/$r2" 2>/dev/null; then
            echo "Trimmed files are valid. Skipping."
            continue
        else
            echo "WARNING: Corrupted trimmed files detected (from a previous failed run). Deleting and re-trimming..."
            rm -f "$TRIM_DIR/$r1" "$TRIM_DIR/$r2"
        fi
    fi

    echo "Running Trimmomatic (Memory: ${JAVA_MEM}G, Threads: ${THREADS})..."
    trimmomatic -Xmx${JAVA_MEM}g PE -threads "$THREADS" \
        "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" \
        "$TRIM_DIR/$r1" "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" \
        "$TRIM_DIR/$r2" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 || { echo "ERROR: Trimmomatic failed for $sample. Aborting."; exit 1; }

    rm -f "$TRIM_DIR/${r1%.fastq.gz}_unpaired.fastq.gz" "$TRIM_DIR/${r2%.fastq.gz}_unpaired.fastq.gz"
done

echo "Upstream preprocessing complete."
