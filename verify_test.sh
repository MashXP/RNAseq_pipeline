#!/bin/bash

# verify_pipeline.sh
# This script sets up a mini-RNA-seq environment using Human Chromosome 21.
# Run this locally to verify the pipeline logic.

set -e

# 1. Directories
BASE_DIR=$(dirname "$(realpath "$0")")
DATA_DIR="$BASE_DIR/_data"
GENOME_DIR="$DATA_DIR/genome_test"
INDEX_DIR="$DATA_DIR/index_test"
FASTQ_DIR="$DATA_DIR/fastq"
BAM_DIR="$DATA_DIR/bam_test"

mkdir -p "$GENOME_DIR" "$INDEX_DIR" "$FASTQ_DIR" "$BAM_DIR"

echo "=== Step 1: Downloading Chromosome 21 (GRCh38) ==="

# Ensembl 113 URLs
FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"

if [ ! -f "$GENOME_DIR/chr21.fa.gz" ]; then
    wget -O "$GENOME_DIR/chr21.fa.gz" "$FASTA_URL"
fi

if [ ! -f "$GENOME_DIR/human_eb113.gtf.gz" ]; then
    wget -O "$GENOME_DIR/human_eb113.gtf.gz" "$GTF_URL"
fi

echo "=== Step 2: Extracting Reference Files ==="
if [ ! -f "$GENOME_DIR/chr21.fa" ]; then
    gunzip -c "$GENOME_DIR/chr21.fa.gz" > "$GENOME_DIR/chr21.fa"
fi

if [ ! -f "$GENOME_DIR/human_eb113.gtf" ]; then
    gunzip -c "$GENOME_DIR/human_eb113.gtf.gz" > "$GENOME_DIR/human_eb113.gtf"
fi

# Filter GTF for Chr21 to speed up indexing and quantify
if [ ! -f "$GENOME_DIR/chr21_annotations.gtf" ]; then
    echo "Filtering GTF for Chr21..."
    grep "^21\|^#" "$GENOME_DIR/human_eb113.gtf" > "$GENOME_DIR/chr21_annotations.gtf"
fi

echo "=== Step 3: Generating Mini STAR Index (~1.5GB RAM) ==="
if [ -z "$(ls -A "$INDEX_DIR")" ]; then
    STAR --runThreadN 4 \
         --runMode genomeGenerate \
         --genomeDir "$INDEX_DIR" \
         --genomeFastaFiles "$GENOME_DIR/chr21.fa" \
         --sjdbGTFfile "$GENOME_DIR/chr21_annotations.gtf" \
         --sjdbOverhang 99 \
         --genomeSAindexNbases 11
else
    echo "Small index already exists. Skipping."
fi

echo "=== Verification Setup Ready ==="
echo ""
echo "Next Steps:"
echo "1. Run subsample_hpc.sh on the HPC."
echo "2. Download 'test_subset/*.fastq.gz' into '$FASTQ_DIR'."
echo "3. Run the pipeline targeting this test index."
echo ""
echo "Note: You can run the alignment with:"
echo "  STAR --genomeDir $INDEX_DIR --readFilesIn ... --outSAMtype BAM SortedByCoordinate"
