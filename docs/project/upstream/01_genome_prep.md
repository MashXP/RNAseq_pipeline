# Upstream Dissection: 01_genome_prep.sh

This script is the "janitor" and the "librarian" of the pipeline. It handles the raw data cleaning and builds the genomic references for both Human and Canine models.

---

## 1. Safety & Directory Setup
```bash
set -e
set -o pipefail
BASE_DIR=$(dirname "$(realpath "$0")")
...
```
- **The Job**: Initializes the environment with safety flags.
- **The Reasoning**: 
    - `set -e`: Stops everything if a single task fails. This prevents "cascading errors" where the pipeline keeps running on bad data.
    - `set -o pipefail`: Ensures that if a command results in an error *inside* a pipe, the script catches it.
    - `BASE_DIR`: Uses absolute paths to avoid "File Not Found" errors if you run the script from a different folder.

---

## 2. Genomic Source Definitions
```bash
FASTA_URLS[Human]="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/..."
GTF_URLS[Human]="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/..."
```
- **The Job**: Maps species names (Human/Canine) to specific Ensembl 113 download links.
- **The Reasoning**: Centralizes your external dependencies. If you ever need to upgrade to Ensembl 114, you only change these lines.

---

## 3. The Species Loop: Mapping & Preparation
```bash
for species in Human Canine; do
    ...
done
```
- **The Job**: This is the core logic that ensures one pipeline can handle both species' genomes in a single execution.
- **The Reasoning**: This supports the mentor's vision of a cross-species study by treating "Species" as a standard variable, not a special case.

### 3.1 Downloading & Extraction
```bash
if [ ! -f "$SPECIES_GENOME_DIR/$fasta_gz" ]; then
    wget -P "$SPECIES_GENOME_DIR" "$fasta_url"
fi
```
- **The Job**: Only downloads and unzips if files are missing.
- **The Reasoning**: Genome files are huge (>1GB). This makes the script "idempotent"—if it crashes, you can restart it and it won't waste time redownloading data it already has.

### 3.2 Picard refFlat Transformation
```bash
gtfToGenePred -genePredExt -geneNameAsName2 "$SPECIES_GENOME_DIR/$gtf_unzipped" "$SPECIES_GENOME_DIR/refflat.tmp.txt"
awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' ... > "$SPECIES_GENOME_DIR/$refflat_file"
```
- **The Job**: Converts modern annotation into a legacy 11-column format required by Picard.
- **The Reasoning**: This enables high-rigor checks (Exon vs Intron distribution) that standard counting tools ignore.

### 3.3 STAR Indexing
```bash
if [ -z "$(ls -A "$species_index")" ]; then
    STAR --runMode genomeGenerate --genomeDir "$species_index" ...
fi
```
- **The Job**: Builds the searchable "index" for the genome.
- **Why it matters**: It converts 3 billion letters of DNA into a binary database that STAR can search in seconds. This uses high-performance settings (`--sjdbOverhang 99`) for standard sequencing.

---

## 4. Pre-Trimming Quality Control (FastQC)
```bash
fastqc -t "$THREADS" "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2" -o "$RAW_QC_DIR"
```
- **The Job**: Generates a report on the raw, untouched sequencing data.
- **The Reasoning**: We run this *before* trimming to see if there were issues at the sequencing facility that we need to be aware of.

---

## 5. Trimming (Trimmomatic) logic

### 5.1 Resource Math
```bash
JAVA_MEM=$(($MEM_GB - 12))
```
- **The Job**: Automatically calculates how much RAM to give to the Java Virtual Machine.
- **The Reasoning**: Trimmomatic is a Java tool. We allocate most of your RAM to it, but leave 12GB for the system to prevent the computer from freezing.

### 5.2 The Trimming Operation
```bash
trimmomatic PE -threads "$THREADS" ... \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
- **The Job**: Cleans the reads by removing adapters and low-quality sequences.
- **Parameters**:
    - `ILLUMINACLIP`: Snips off the synthetic Illumina DNA.
    - `SLIDINGWINDOW:4:15`: Cleans "blurry" reads where quality drops below 97% accuracy.
    - `MINLEN:36`: Discards reads too short to be unique (to avoid false mapping).

> [!NOTE]
> **Optional Step**: If your raw FASTQ data is already clean (e.g., provided as pre-processed by a sequencing facility), you can skip this tool and place your files directly into `_data/fastq_trimmed/`.
