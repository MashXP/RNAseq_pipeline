# MiaPAca-2 Cancer RNA-seq Pipeline

This directory contains the automated end-to-end pipeline for analyzing MiaPAca-2 pancreatic cancer RNA-seq data (human GRCh38).

## Research Purpose
The objective of this study is to investigate the dose-dependent transcriptomic response of **MiaPAca-2 pancreatic cancer cells** to **Kromastat** treatment (0.5nM, 1nM, and 5nM) using high-throughput RNA-seq. The pipeline is designed to:
1.  **Identify Differentially Expressed Genes (DEGs)** for each dosage relative to a DMSO (NG) baseline.
2.  **Compare overlapping drug effects** across concentrations using multi-way Venn Diagrams.
3.  **Perform Functional Enrichment (GSEA & ORA)** to identify activated and suppressed signaling pathways, Hallmark gene sets, and biological processes driven by Kromastat.

## Directory Structure

```text
MiaPAca-2_pipeline/
├── _data/                     # Upstream data files (git ignored)
│   ├── fastq/                 # Raw .fastq.gz files
│   ├── fastq_trimmed/         # Adapter/Quality trimmed reads (Trimmomatic)
│   ├── genome/                # Reference genome
│   ├── index/                 # STAR index
│   ├── bam/                   # Aligned BAMs
│   ├── counts/                # Production gene counts
│   ├── counts_test/           # Test gene counts (Chr21 subsampled)
│   ├── rseqc/                 # QC results
│   └── multiqc/               # Aggregated report
├── scripts_upstream/          # Upstream automation (BASH)
│   ├── 01_genome_prep.sh      # Genome prep + Trimmomatic trimming
│   ├── 02_star_align.sh       # STAR mapping (uses trimmed reads)
│   ├── 03_rseqc_qc.sh
│   ├── 04_quantification.sh
│   ├── 05_multiqc.sh
│   └── utils/                 # Format conversion & parsers
├── scripts_downstream/        # Production downstream analysis (R)
│   ├── 01_data_prep.R
│   ├── 02_deseq2_dge.R
│   ├── 03_enrichment.R
│   └── 04_visualization.R
├── test_upstream/             # Test upstream scripts (Chr21)
│   ├── 01_genome_prep.sh
│   ├── 02_star_align.sh
│   ├── 03_rseqc_qc.sh
│   ├── 04_quantification.sh
│   ├── 05_multiqc.sh
│   └── utils/
├── test_downstream/           # Test downstream analysis (R, relaxed thresholds)
│   ├── 01_data_prep.R
│   ├── 02_deseq2_dge.R
│   ├── 03_enrichment.R
│   └── 04_visualization.R
├── results/                   # Production output (git ignored)
│   ├── figures/               # PCA, Volcano, Heatmaps, NES, GSEA Dotplots
│   └── tables/                # DGE results, GO/KEGG/Hallmark tables
├── results_test/              # Test output (git ignored)
│   ├── figures/
│   └── tables/
├── MiaPAca-2_Sample_Data_Table.csv
├── scripts_upstream_docs.md   # Upstream pipeline documentation
├── scripts_downstream_docs.md # Downstream analysis documentation
├── environment.yml            # Conda/Mamba environment definition
├── setup_env.sh               # Environment setup script
├── subsample_hpc.sh           # HPC subsampling utility for test data
├── run                        # Master runner (Production)
├── test                       # Master runner (Test/Chr21)
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
The test suite utilizes a subsampled dataset and a single chromosome (Chr21) to minimize RAM and storage requirements. Test scripts use relaxed statistical thresholds to ensure meaningful output from the limited gene set.

| Mode | Target | STAR RAM | Time |
| :--- | :--- | :--- | :--- |
| **Production** | Full GRCh38 | ~30GB | 8-10 Hours |
| **Test** | Chr21 Only | ~1.5GB | ~10 Minutes |

### 2. Usage
```bash
./test all      # Run entire verification (Upstream -> Downstream)
./test up       # Run BASH upstream only
./test down     # Run R analysis only
```

For more details on test thresholds and script logic, refer to the internal documentation headers in `test_upstream/` and `test_downstream/`.

## Downstream Outputs

The visualization script (`04_visualization.R`) generates the following figures:

| Figure | Description |
| :--- | :--- |
| PCA Plot | Sample clustering by drug dose |
| Volcano Plots | Individual + combined per-dose DEG plots |
| Venn Diagram | DEG overlap: Kromastat concentrations vs DMSO |
| Pathway Heatmap | Top enriched GO pathways with gene expression |
| Top Variable Heatmap | Top 30 (test) / 50 (prod) most variable genes |
| Hallmark NES Barplots | Individual + combined Hallmark enrichment scores |
| GSEA Dotplots | Individual + combined GO GSEA enrichment |

## Requirements

### Bioinformatics Tools
STAR, samtools, featureCounts (subread), RSeQC, MultiQC, Trimmomatic.

### Python
`pandas` (required for `biotype_to_multiqc.py`).

### R Libraries
`DESeq2`, `tidyverse`, `ggplot2`, `clusterProfiler`, `org.Hs.eg.db`, `EnhancedVolcano`, `ComplexHeatmap`, `circlize`, `ggVennDiagram`, `patchwork`, `msigdbr`, `enrichplot`, `apeglm`.

> All R/Python dependencies are managed via `environment.yml`. Run `bash setup_env.sh` to install everything.

For technical details, see the documentation for [upstream scripts](scripts_upstream_docs.md) and [downstream analysis](scripts_downstream_docs.md).
