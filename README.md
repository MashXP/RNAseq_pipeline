# RNA-seq Pipeline

This directory contains the automated end-to-end pipeline for analyzing cancer RNA-seq data, now organized into a modular, domain-specific architecture.

## Research Purpose
The objective of this pipeline is to investigate transcriptomic responses (e.g., drug treatments) using high-throughput RNA-seq. The pipeline is designed to:
1.  **Identify Differentially Expressed Genes (DEGs)** and Alternative Splicing events.
2.  **Compare overlapping effects** across treatments and species.
3.  **Perform Functional Enrichment (GSEA & ORA)** to identify activated and suppressed signaling pathways.

## Modular Architecture

The repository is organized into three primary analysis modules, segregating Core Alignment, DGE, and Alternative Splicing (AS).

```text
RNAseq_pipeline/
├── 01_core_alignment/     # Module: Genome Indexing, STAR Mapping, MultiQC
├── 02_dge_analysis/       # Module: featureCounts, DESeq2, Enrichment, DGE Plots
├── 03_as_analysis/        # Module: rMATS-turbo (Upstream & Downstream)
├── utils/                 # Centralized Python/Bash/R utilities
├── _data/                 # Centralized Data Storage (git ignored)
│   └── _metadata/
│       └── drPhuong_Sample_Data_Table.csv # STUDY METADATA
├── results/               # Modular results (git ignored)
│   ├── 01_core_alignment/
│   ├── 02_dge_analysis/
│   └── 03_as_analysis/
├── docs/                  # Documentation
│   └── test_deviation.md  # Historical log of test-suite overrides
├── environment.yml        # Conda/Mamba environment definition
└── run                    # Master Modular Runner
```


## Quick Start

The master `run` script now targets modules directly.

### 1. Execute Specific Module
```bash
./run core    # Alignment & QC
./run dge     # Differential Expression analysis
./run as      # Alternative Splicing analysis
```

### 2. Execute Full Pipeline
```bash
./run all     # Runs Core -> DGE -> AS in sequence
```

## Pipeline Architecture

### 01_core_alignment (Scripts)
- `01_genome_prep.sh`: Indexing & Trimming (STAR, Trimmomatic)
- `02_star_align.sh`: Mapping (Alignment) (STAR)
- `03_alignment_qc.sh`: Health Check (Picard, QC Stats)
- `04_multiqc.sh`: Final Reporting (MultiQC)

### 02_dge_analysis (Scripts)
- `upstream/01_quantification.sh`: Counting Genes (featureCounts)
- `downstream/01_data_prep.R`: Clean & organize metadata
- `downstream/02_deseq2_dge.R`: DGE Engine (DESeq2)
- `downstream/03_enrichment.R`: Functional analysis (GSEA, ORA)
- `downstream/04_pca.R`: Sample clustering / PCA
- `downstream/05_volcano.R`: DEG significance (Volcano)
- `downstream/06_venn.R`: Multi-contrast overlap (Venn)
- `downstream/07_heatmap_pathway.R`: Top enriched pathways gene expression
- `downstream/08_heatmap_variable.R`: Top 50 most variable genes
- `downstream/09-13_...`: Modular plotting suite (NES, Dotplots, UpSet, Correlation)
- `downstream/14-16_...`: Comparative engine (Alluvial, Ortholog Overlaps)

### 03_as_analysis (Scripts)
- `upstream/01_rmats_splicing.sh`: Splicing detection (rMATS-turbo)
- `downstream/01_rmats_summary.R`: Splicing event summaries

## Documentation

Full, line-by-line documentation for every script is located within the `docs/` subfolder of each module.

*(Note: The legacy verification test suite has been removed for redesign; see [docs/test_deviation.md](docs/test_deviation.md) for details on previous implementation logic).*

---
*This pipeline builds upon methodological patterns established in the [YeastAnalysis](https://github.com/MashXP/YeastAnalysis) project.*
