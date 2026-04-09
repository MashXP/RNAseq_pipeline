# MiaPAca-2 Cancer RNA-seq Pipeline

This directory contains the automated end-to-end pipeline for analyzing MiaPAca-2 pancreatic cancer RNA-seq data (human GRCh38).

## Research Purpose
The objective of this study is to investigate the dose-dependent transcriptomic response of **MiaPAca-2 pancreatic cancer cells** to **Kromastat** treatment (0.5nM, 1nM, and 5nM) using high-throughput RNA-seq. The pipeline is designed to:
1.  **Identify Differentially Expressed Genes (DEGs)** for each dosage relative to a DMSO (NG) baseline.
2.  **Compare overlapping drug effects** across concentrations using multi-way Venn Diagrams.
3.  **Perform Functional Enrichment (GSEA)** to identify activated and suppressed signaling pathways and biological processes driven by Kromastat.

## Directory Structure

```text
MiaPAca-2_pipeline/
├── _data/                   # Upstream data files (git ignored)
│   ├── fastq/               # Raw .fastq.gz files
│   ├── fastq_trimmed/       # Adapter/Quality trimmed reads (Trimmomatic)
│   ├── genome/              # Reference genome
│   ├── index/               # STAR index
│   ├── bam/                 # Aligned BAMs
│   ├── counts/              # Gene counts
│   ├── rseqc/               # QC results
│   └── multiqc/             # Aggregated report
├── scripts_upstream/        # Upstream automation (BASH)
│   ├── 01_genome_prep.sh    # Genome prep + Trimmomatic Trimming
│   ├── 02_star_align.sh     # STAR mapping (uses trimmed reads)
│   ├── 03_rseqc_qc.sh
│   ├── 04_quantification.sh
│   ├── 05_multiqc.sh
│   └── utils/               # Format conversion & parsers
├── scripts_downstream/      # Downstream analysis (R)
│   ├── 01_data_prep.R
│   ├── 02_deseq2_dge.R
│   ├── 03_enrichment.R
│   ├── 04_visualization.R
├── results/                 # Final Output
│   ├── figures/             # PCA, Volcano, Heatmap
│   └── tables/              # DGE results, GO/KEGG tables
├── MiaPAca-2_Sample_Data_Table.csv
├── scripts_upstream_docs.md   # Upstream pipeline documentation
├── scripts_downstream_docs.md # Downstream analysis documentation
├── setup_env.sh               # Environment setup script
├── run                        # Master runner (Upstream + Downstream)
└── README.md                  # This file
```

## Prerequisites & Environment Setup

It is highly recommended to use **Mamba** (or Conda) to manage the pipeline's dependencies.

1. **Initialize the Environment**:
   ```bash
   bash setup_env.sh
   ```
2. **Activate the Environment**:
   ```bash
   mamba activate cancer_rnaseq
   ```

## Quick Start

### 1. Upstream (Alignment & Quantification)
Place raw FASTQ files in `_data/fastq/`. You can run only the upstream pipeline using:
```bash
./run up
```

### 2. Downstream (Differential Expression)
Run only the downstream R analysis using:
```bash
./run down
```

### 3. Full Pipeline
Run everything from start to finish:
```bash
./run all
```

## Verification & Testing (Local Run)

For rapid verification of the pipeline logic on local hardware, a **Chromosome 21 Test Suite** is provided. This allows you to run the entire pipeline in minutes instead of hours.

### 1. Overview
The test suite utilizes a subsampled dataset and a single chromosome (Chr21) to minimize RAM and storage requirements.

| Mode | Target | STAR RAM | Time |
| :--- | :--- | :--- | :--- |
| **Production** | Full GRCh38 | ~30GB | 8-10 Hours |
| **Test** | [[Chr21 Only]] | ~1.5GB | ~10 Minutes |

### 2. Usage
```bash
./test all      # Run entire verification (Upstream -> Downstream)
./test up       # Run BASH upstream only
./test down     # Run R analysis only
```

For more details on test thresholds and script logic, refer to the internal documentation headers in [[test_upstream/]] and [[test_downstream/]].

## Requirements
- **Bioinformatics**: STAR, samtools, featureCounts (subread), RSeQC, MultiQC, Trimmomatic.
- **Python**: `pandas` (required for `biotype_to_multiqc.py`). 
- **R Libraries**: `DESeq2`, `tidyverse`, `clusterProfiler`, `org.Hs.eg.db`, `EnhancedVolcano`, `pheatmap`, `ggVennDiagram`.

For technical details, see the documentation for [upstream scripts](scripts_upstream_docs.md) and [downstream analysis](scripts_downstream_docs.md).
