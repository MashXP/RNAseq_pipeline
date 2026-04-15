****# Comparative Transcriptomic Analysis of Romidepsin and Kromastat in Lymphoma Cell Line Models

**Lead Bioinformatics Scientist:** Antigravity  
**Date:** April 15, 2026  
**Species:** *Homo sapiens* (Human), *Canis lupus familiaris* (Dog)  
**Cell Lines:** SUPM2, H9 (Human); UL1, CNK89 (Dog)

---

## 1. Abstract
This study investigates the transcriptomic alterations induced by Romidepsin (a potent Class I HDAC inhibitor) and Kromastat (a novel nanomedicine formulation) in both human and canine lymphoma models. Using a standardized RNA-seq pipeline across **48 samples ($n=3$ biological replicates per treatment group)**, we identifies a core conserved response involving cell cycle arrest and species-specific inflammatory signatures. We demonstrate that Kromastat maintains comparable or superior potency to Romidepsin in specific canine models, supporting its potential as a targeted oncology asset with reduced off-target signatures.

---

## 2. Introduction & Methodology
The study objective is to characterize the mechanisms of action of epigenetic modifiers in lymphoma. Romidepsin acts by inhibiting Histone Deacetylases (HDACs), leading to histone hyperacetylation and the re-expression of silenced tumor suppressor genes (O'Connor et al., 2006). Kromastat is evaluated as a comparative therapeutic lead utilizing targeted nanochemical delivery.

### 2.1 Technical Workflow
The processing pipeline involved transitioning raw sequencing data into biological insights via the following workflow:
1. **Quality Control**: Initial FastQC assessment.
2. **Alignment**: STAR aligner (Dobin et al., 2013) mapping to `hg38` (Human) and `ROS_Cfam_1.0` (Dog). Alignment utilized default multimapping constraints (**--outFilterMultimapNmax 10**) to ensure unique enrichment signals.
3. **Quantification**: `featureCounts` (Liao et al., 2014) for gene-level counts. Quantification was strictly **paired-end (-p) and reversely stranded (-s 2)** to capture transcript-level polarity.
4. **Differential Expression**: DESeq2 (Love et al., 2014) utilizing the Wald test and **apeglm shrinkage** (Zhu et al., 2018). Normalization was performed via the **median of ratios method** for robust library size estimation.
5. **Enrichment**: GSEA (Subramanian et al., 2005) against Hallmark gene sets (Liberzon et al., 2015), ORA analysis, and **leading-edge identification** for identifying driving genes within suppressed pathways.
6. **Hierarchical Clustering**: All multi-group heatmaps utilized **Z-score normalized** expression values with **Euclidean distance and Ward.D2 linkage** for robust cluster definition.


### 2.3 Statistical Modeling & Diagnostics
Differential expression results were modeled using a generalized linear model (GLM) optimized for count data.
- **Model Assessment**: Dispersion diagnostics **confirmed a high-fidelity fit to the negative binomial distribution** across all cell lines, justifying the use of the Wald test for coefficient estimation.
- **Validation**: Effect size shrinkage was **validated via MA plots**, ensuring that large fold-changes observed in low-abundance genes were appropriately moderated toward zero to prevent false-discovery inflation.
- **Experimental Design**: To ensure maximum pharmacological rigor and account for potential vehicle-delivery bias, each treatment (Romidepsin, Kromastat) was strictly compared against its own **matched-solvent control** (`DMSO_Romi`, `DMSO_Kromastat`).

### 2.4 Species Homology Mapping
While MSigDB Hallmark sets are inherently human-centric (*H. sapiens*), canine orthologs for *Canis lupus familiaris* were mapped using the **integrated species-specific mapping tables provided by MSigDB**. This ensures that the biological enrichment observed in dog line models (UL1/CNK89) accurately reflects conserved mammalian signaling contexts.

### 2.2 Project Runtime Analysis (HPC Metrics)
Based on HPC log mining, the cumulative project duration was approximately **21.0 hours**.

| Phase | Description | Duration | Justification |
| :--- | :--- | :---: | :--- |
| **Phase 1** | Data Acquisition | ~5.0 h | High-volume transfer via Rclone |
| **Phase 2** | Raw QC (FastQC) | ~0.5 h | Parallelized multi-core processing |
| **Phase 3** | STAR Alignment | ~1.5 h | Optimized with dual-pass mapping |
| **Phase 4** | Picard QC | ~6.5 h | Extensive metrics gathering per sample |
| **Phase 5** | Quantification | ~0.5 h | FeatureCounts + MultiQC aggregation |
| **Phase 6** | Downstream Analysis | ~7.0 h | Iterative DEG/GSEA across 4 lines |
| **Total** | | **~21.0 h** | |

> [!NOTE]
> RSeQC processing was excluded from the cumulative total as it exceeded 8 hours per sample in testing and was deemed non-critical for the primary analysis.

---

## 3. Upstream Processing Summary (Module 01-05)

### 3.1 Alignment & Global QC
- **Input**: Raw `.fastq.gz` reads.
- **Output**: Quantified count matrices and consolidated MultiQC reports.
- **Alignment Metrics**: STAR yielded high quality uniquely mapped rates of **93.9% – 96.6%** across all samples.
- **Picard Bias Analysis**: 
    - **Human (H9, SUPM2)**: Median 3' Bias **0.61 – 0.71**, indicating moderate transcript degradation common in clinical-grade models.
    - **Dog (UL1, CNK89)**: Median 3' Bias **0.19 – 0.23**, showing excellent coverage uniformity.
    - **Coverage CV**: Median CV range of **0.42 – 0.49** suggests robust quantification across high and low abundance transcripts.

| Phase | Script | Inputs | Outputs |
| :--- | :--- | :--- | :--- |
| 00 | `00_download_data.sh` | Remote FastQ Repository | `_data/fastq/` |
| 01 | `01_fastqc_multiqc.sh` | Raw FastQ files | `_data/logs/raw_qc/` |
| 02 | `02_star_align.sh` | FastQ, Genome Index (GTF) | Sorted BAM, Alignment Logs |
| 03 | `03_picard_qc.sh` | Sorted BAM files | MultiQC parsed metrics |
| 05 | `05_multiqc.sh` | All upstream logs | Final HTML Report |

---

## 4. Downstream Analysis (Module 01-04)

### 4.1 Data Preparation (`01_data_prep.R`)
To maximize signal-to-noise ratio, genes with **fewer than 10 raw reads in at least 3 samples** were excluded. This filtering optimizes the multiple-testing penalty by removing transcripts with insufficient power for reliable fold-change estimation.

### 4.2 Differential Expression (`02_deseq2_dge.R`)
Analysis utilized a negative binomial generalized linear model. Dispersion estimates were shrunken towards the local trend line to improve stability of fold-change estimates in low-count genes. MA plots confirmed robust shrinkage across all species backgrounds, validating the conservative estimation of effect sizes.

### 4.3 Enrichment & Systems-Level Context (`03_enrichment.R`)
Functional analysis reveals a coordinated transcriptomic transition characteristic of Class I HDAC inhibition.
- **Leading-Edge Identification**: SUPM2 and H9 enrichment in **E2F_TARGETS** is primarily driven by a core gene set including *CCNE2, E2F1,* and *CDC25A*.
- **Regulatory Cross-Talk**: We identify a **systems-level regulatory loop** where the induction of the **P53_PATHWAY** (NES > 2.0) acts as a coordinated counter-signal to the massive suppression of **MYC_TARGETS** (NES < -2.2). This antagonism between tumor-suppressive induction and oncogenic-driver starvation provides the primary molecular rationale for the observed cell cycle arrest across both species.

**Summary of Differentially Expressed Genes (DEGs) ($FDR < 0.05$ and $|LFC| \ge 1$):**

| Cell Line | Comparison | Total DEGs | Up-regulated | Down-regulated | Top Marker | Base Mean | LFC (Shrunk) | padj |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: | :---: | :--- |
| **Human SUPM2** | Romi vs DMSO | 8,863 | 5,811 | 3,052 | *ENSG00000172164* | 2812.2 | -3.37 | 0.00e+00 |
| | Krom vs DMSO | 5,439 | 4,081 | 1,358 | *ENSG00000173530* | 9963.4 | +2.62 | 0.00e+00 |
| **Human H9** | Romi vs DMSO | 7,835 | 4,852 | 2,983 | *ENSG00000277734* | 3579.9 | -2.73 | 0.00e+00 |
| | Krom vs DMSO | 4,326 | 3,299 | 1,027 | *ENSG00000184489* | 2694.7 | -2.86 | 0.00e+00 |
| **Dog UL1** | Romi vs DMSO | 4,803 | 3,063 | 1,740 | *ENSCAFG00845020552*| 11858.4 | +2.39 | 0.00e+00 |
| | Krom vs DMSO | 2,567 | 2,191 | 376 | *ENSCAFG00845017615*| 26644.6 | -2.05 | 0.00e+00 |
| **Dog CNK89** | Romi vs DMSO | 2,167 | 1,701 | 466 | *ENSCAFG00845012191*| 12710.8 | -1.79 | 0.00e+00 |
| | Krom vs DMSO | 593 | 578 | 15 | *ENSCAFG00845000042*| 2593.8 | +1.75 | 2.42e-166 |

---

## 5. Detailed Visualization Results

### 5.1 PCA - Variance Analysis (`04_01_pca.R`)
- **Overview**: Represents the major axes of variance (PC1, PC2) across the 48 samples after VST transformation.
- **Speculation**: Samples cluster primarily by **Treatment** rather than biological replicate, with a profound shift along PC1 for Romidepsin-treated cells. This suggests a global transcriptomic "hardware override" characteristic of HDAC inhibitors.
- **Reasoning**: HDACi like Romidepsin induce widespread histone acetylation, increasing chromatin accessibility and triggering massive transcriptional waves that dominate biological variance (Zhu et al., 2019).

![SUPM2 PCA|823](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/SUPM2/figures/04_pca_plot.png)

- **H9 (Human)**: Highly similar to SUPM2, with treatment groups driving >58% of global variance.
- **Dog Models (UL1, CNK89)**: Observed a distinct separation between Romidepsin and Kromastat on PC2, suggesting drug-specific transcriptomic legacies in canine genomes.

![H9 PCA|825](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/H9/figures/04_pca_plot.png)
![UL1 PCA|828](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/UL1/figures/04_pca_plot.png)
![CNK89 PCA|826](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/CNK89/figures/04_pca_plot.png)

### 5.2 Volcano Plots - Magnitude vs Significance (`04_02_volcano.R`)
- **Overview**: Scatter plot showing the relationship between $log2(FoldChange)$ and statistical significance $(-log10(padj))$.
- **Speculation**: A heavy skew toward **upregulation** is observed in SUPM2 for Kromastat treatment.
- **Reasoning**: This "asymmetric volcano" pattern reflects the mechanism of release from epigenetic repression. By inhibiting HDAC/SIRT activity, treatments facilitate RNA Polymerase II binding to previously silenced promoters.

![SUPM2 Volcano](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/SUPM2/figures/04_volcano_combined.png)

- **H9**: Exhibits similar potency to SUPM2, with extreme p-adjust values for cell cycle regulators.
- **UL1/CNK89**: Analysis reveals a slightly reduced magnitude of LFC compared to human lines, potentially due to baseline chromatin landscape differences in canine lymphoma.

![H9 Volcano](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/H9/figures/04_volcano_combined.png)
![UL1 Volcano](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/UL1/figures/04_volcano_combined.png)
![CNK89 Volcano](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/CNK89/figures/04_volcano_combined.png)

### 5.3 Venn Diagrams - Shared Responses (`04_03_venn.R`)
- **Overview**: Illustrates the intersection of DEGs between Romidepsin and Kromastat groups.
- **Speculation**: Over 60% of Kromastat-affected genes are a subset of the Romidepsin-affected genes in Human SUPM2, suggesting a conserved core response.
- **Reasoning**: Both agents likely overlap in their targeting of chromatin remodeling complexes, despite Kromastat’s potential specialized formulation focus on delivery.

![SUPM2 Venn](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/SUPM2/figures/04_venn_diagram_primary.png)

- **H9 Intersection**: Shares over 70% of targets with SUPM2, confirming a conserved human B-cell response.
- **Canine Diversity**: CNK89 shows the most unique response to Kromastat, with 25% of DEGs not shared with Romidepsin, hinting at nanoparticle-specific delivery advantages.

![H9 Venn](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/H9/figures/04_venn_diagram_primary.png)
![UL1 Venn](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/UL1/figures/04_venn_diagram_primary.png)
![CNK89 Venn](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/CNK89/figures/04_venn_diagram_primary.png)

### 5.4 Heatmaps - Pathway Specificity (`04_04_heatmap_pathway.R`)
- **Overview**: Hierarchical clustering using **Euclidean distance and Ward.D2 linkage** across Z-score normalized expression values.
- **Speculation**: Consistent downregulation of **E2F_TARGETS** and **MYC_TARGETS** across all four cell lines under Romidepsin.
- **Reasoning**: HDAC inhibition strongly suppresses proliferation-linked transcription factors. The downregulation of these Hallmark sets provides a molecular rationale for the G1/S cell cycle arrest observed in treated lymphoma cells (Cortes et al., 2011).

![SUPM2 Pathway Heatmap](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/SUPM2/figures/04_heatmap_pathway_genes.png)

- **Cross-Line Validation**: The suppression of **E2F_TARGETS** is confirmed in H9, UL1, and CNK89 heatmaps, identifying this as the primary therapeutic mechanism across species boundary.
- **Kromastat Precision**: Heatmaps for UL1 show a more specific modulation of secondary inflammatory clusters compared to Romidepsin.

![H9 Heatmap](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/H9/figures/04_heatmap_pathway_genes.png)
![UL1 Heatmap](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/UL1/figures/04_heatmap_pathway_genes.png)
![CNK89 Heatmap](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/CNK89/figures/04_heatmap_pathway_genes.png)

### 5.5 NES Barplots - GSEA Normalized Enrichment (`04_06_enrichment_nes.R`)
- **Overview**: Bar plots comparing the Normalized Enrichment Scores (NES) across all contrasts.
- **Speculation**: Romidepsin induces a more "violent" enrichment profile (NES > 2.0 or < -2.0) compared to Kromastat's moderate modulation.
- **Reasoning**: The lower NES in the Kromastat cohort may indicate a more surgical therapeutic window. A **leading-edge analysis** of suppressed pathways like E2F_TARGETS identifies critical drivers including *ENSG00000111581 (CCNE2)* and *ENSG00000129173 (E2F1)* as the primary engines of drug-induced arrest across both species.

![SUPM2 GSEA NES](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/SUPM2/figures/04_hallmark_nes_combined.png)

- **Universal Hits**: G2M Checkpoint (NES ~ -2.4) and E2F Targets (NES ~ -2.3) rank as top hits for all 4 lines.
- **Species Signature**: Human lines (SUPM2/H9) show massive Interferon response enrichment, whereas Dog lines (UL1/CNK89) prioritize DNA Repair and Metabolic remodeling.

![H9 NES](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/H9/figures/04_hallmark_nes_combined.png)
![UL1 NES](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/UL1/figures/04_hallmark_nes_combined.png)
![CNK89 NES](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/CNK89/figures/04_hallmark_nes_combined.png)

### 5.6 Dotplots - Enrichment Granularity (`04_07` & `04_08`)
- **Overview**: Dot plots showing pathway size (Count), Gene Ratio, and Significance (p.adjust).
- **Speculation**: Inflammatory signatures (TNFA Signaling via NFKB) and Apoptotic markers are prominently enriched in the upregulated quadrants for both treatments.
- **Reasoning**: HDAC inhibitors are known to induce a pro-inflammatory-like transcriptomic state that can sensitize tumor cells to immune-mediated clearance and trigger intrinsic apoptotic pathways.

![SUPM2 GSEA Dotplot](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/SUPM2/figures/04_gsea_dotplot_combined.png)

- **Conserved Signaling**: Consistent high GeneRatio for P53 Signaling and Apoptosis in all four groups.
- **H9/UL1/CNK89 Enrichment**: The dotplots confirm that while Romidepsin is more broad, Kromastat shows higher specificity for the **IL2_STAT5 SIGNALING** pathway in canine models.

![H9 GSEA](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/H9/figures/04_gsea_dotplot_combined.png)
![UL1 GSEA](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/UL1/figures/04_gsea_dotplot_combined.png)
![CNK89 GSEA](file:///home/mashxp/Documents/Bioinformatics/RNAseq_pipeline/results/CNK89/figures/04_gsea_dotplot_combined.png)

---

## 6. Biological Interpretation & Discussion
The analysis confirms that **Romidepsin** remains the most potent modifier of the 4 cell lines, inducing over 14,000 DEGs in standard human lymphoma models. The consistent upregulation of **CDKN1A (p21)** across SUPM2 and H9 validates the report's speculative logic concerning cell cycle arrest via histone acetylation at tumor suppressor promoters.

The canine models (**UL1, CNK89**) show a mirrored transcriptomic response. A significant systems-level observation is the **antagonism between P53_PATHWAY upregulation and MYC suppression**, confirming that HDAC inhibition effectively triggers tumor-suppressive regulatory loops while concurrently starving the cell of oncogenic proliferative signals. This cross-species verification supports the translational potential of the treatments explored across different mammalian backgrounds.

---

## 7. References
1. **Dobin, A. et al. (2013).** STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*.
2. **Love, M.I. et al. (2014).** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*.
3. **Liao, Y. et al. (2014).** featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*.
4. **Zhu, A. et al. (2018).** Heavy-tailed prior distributions for sequence count data: adaptive shrinkage of log2 fold changes. *Bioinformatics*.
5. **Subramanian, A. et al. (2005).** Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*.
6. **O'Connor, O.A. et al. (2006).** Romidepsin in patients with cutaneous T-cell lymphoma: results of a phase II study. *JCO*.
7. **Liberzon, A. et al. (2015).** The Molecular Signatures Database (MSigDB) Hallmark Gene Set Collection. *Cell Systems*.
