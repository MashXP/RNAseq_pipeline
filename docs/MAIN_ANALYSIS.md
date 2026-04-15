# RNA-seq Pipeline Analysis: Detailed Report

## 1. Technical Validation (Upstream)
- **Mapping Quality:** **EXCELLENT**. Unique alignment rates range from **92.8%** (Sample 30) to **96.6%** (Sample 73).
    - *Source: `_data/multiqc/rna_seq_pipeline_summary_data/star_summary_table.txt`*
- **Library Strandedness:** Verified as **Reverse-Stranded (ISR)**.
    - *Source: Comparison of `star_gene_counts_Reverse_Stranded.txt` (~22M/sample) vs `star_gene_counts_Same_Stranded.txt` (~0.8M/sample).*
- **Intergenic Rate:** Sample 30 (Romi) shows **17.8% intergenic bases** vs 15.1% in control (Sample 27).
    - *Source: `_data/multiqc/rna_seq_pipeline_summary_data/multiqc_picard_RnaSeqMetrics.txt`, Column: `PCT_INTERGENIC_BASES`.*
- **Raw QC:** FastQC initial check has been **integrated** into `01_genome_prep.sh` for future runs.
    - *Source: `scripts_upstream/01_genome_prep.sh` updated to run FastQC on raw reads.*

## 2. Biological Interpretation: Romidepsin (6nM)

### A. Human Response (SUPM2 / H9)
- **Primary Mechanism:** Potent **Anti-Proliferative** arrest. 
- **Suppressed:** DNA Replication (padj = 7.08e-14), Ribosome Biogenesis (padj = 2.86e-08).
    - *Source: `results/SUPM2/tables/03_go_ora_Romi_6nM_vs_DMSO_Romi.csv`.*
- **Activated:** Vesicle/Vacuole Organization (NES = 1.75, padj = 8.15e-07).
    - *Source: `results/SUPM2/tables/03_gsea_activated_Romi_6nM_vs_DMSO_Romi.csv`.*
- **H9 Specific:** Strong **B-cell activation** (padj = 1.81e-08).
    - *Source: `results/H9/tables/03_go_ora_Romi_6nM_vs_DMSO_Romi.csv`.*

### B. Dog Response (UL1 / CNK89)
- **Primary Mechanism:** **Structural & Immunomodulatory** remodeling.
- **Activated:** Cell Adhesion (padj = 0.012), Integrin Signaling (padj = 1.39e-11).
    - *Source: `results/UL1/tables/03_go_ora_Romi_6nM_vs_DMSO_Romi.csv` and `03_kegg_ora_...csv`.*
- **Significance:** In UL1, **10,773 / 14,320** genes (~75%) are significantly differentially expressed (padj < 0.05).
    - *Source: `results/UL1/tables/02_dge_Romi_6nM_vs_DMSO_Romi.csv`.*

### C. Conserved Signatures (Cross-Species)
- **Cytoskeletal Remodeling:** All lines show impact on **Cilium Organization** (SUPM2 padj = 5.28e-12; H9 padj = 3.90e-08; CNK89 padj = 0.006).
    - *Source: Cross-comparison of `03_go_ora_...` tables across all result directories.*

## 3. Drug Comparison: Romidepsin vs Kromastat
- **Superiority:** Romidepsin significantly outperforms Kromastat in suppressing **E2F Targets** (NES = -1.66) and **G2M Checkpoint** (NES = -1.70).
    - *Source: `results/UL1/tables/03_gsea_hallmark_Romi_6nM_vs_Kromastat_6nM.csv`.*
- **Inflammation:** Romidepsin provides stronger suppression of TNFA Signaling via NFKB (NES = -1.52) compared to Kromastat in Dog UL1.
    - *Source: `results/UL1/tables/03_gsea_hallmark_Romi_6nM_vs_Kromastat_6nM.csv`.*

## 4. Conclusion
Data confirms Romidepsin as more potent than Kromastat. Response is species-divergent: Human lines show direct cell cycle arrest, while Dog lines prioritize structural/adhesion remodeling and immune activation.
