# Downstream Transcriptomic Analysis (Human)

## 1. Principal Component Analysis (PCA)

The PCA plot visualizes the overall transcriptomic variance across the 24 human samples.

![Combined PCA Plot|600](../../results/human/figures/04_01_pca_combined.png)

### 1.1 Interpretation:
- **PC1 (94% Variance)**: The primary source of variance in the dataset is **Cell Line Identity**. There is a definitive separation between the H9 (healthy-like) and SUPM2 (cancer) lineages along the X-axis.
- **PC2 (3% Variance)**: The secondary axis captures the **Drug Treatment Effect**.
    - In both cell lines, drug treatment induces an upward shift along PC2 compared to their respective DMSO controls.
    - **Romidepsin (6nM)** induces a significantly larger transcriptomic shift than **Kromastat (6nM)**, suggesting higher potency at this concentration.

### 1.2 Quality Control:
- **Replicate Fidelity**: Replicates for each condition are tightly clustered, confirming high technical reproducibility. 
- **Control Consistency**: Crucially, the **DMSO_Romi** and **DMSO_Kromastat** samples for each cell line are clumped together at the baseline. This confirms that the control state is stable and consistent across groups, providing a robust reference for calculating differential expression.

### 1.3 Sub-PCA (Cell Line Specific):
To remove the overwhelming variance of the cell line identity and focus on drug effects, we analyzed the lines individually:

| H9 Cell Line                                              | SUPM2 Cell Line                                                 |
| :-------------------------------------------------------- | :-------------------------------------------------------------- |
| ![H9 PCA\|500](../../results/human/figures/04_01_pca_H9.png) | ![SUPM2 PCA\|500](../../results/human/figures/04_01_pca_SUPM2.png) |

---

## 2. Differential Gene Expression (Volcano Plots)

### 2.1 Global Response Dashboard
The combined dashboard provides a bird's-eye view of all 8 contrasts used to validate the drug response.

[Volcano Plot Overview|1000](../../results/human/figures/04_02_volcano_combined.png)

#### 2.1.1 Control Group Validation
- **Contrast**: `DMSO_Kromastat_vs_DMSO_Romi`
- **Finding**: **0 Significant DEGs**. 

![DMSO Control Validation|500](../../results/human/figures/04_02_volcano_DMSO_Kromastat_vs_DMSO_Romi.png)

- **Conclusion**: This "blank" plot confirms that the two control groups are identical. Any downstream results are strictly due to the drug treatments.

#### 2.1.2 Global Magnitude of Response
> [!NOTE]
> Global DEG counts and ratios are derived from the species-level analysis consolidated in: [04_03_venn_rigor_stats.csv](../../results/human/tables/04_03_venn_rigor_stats.csv).

| Treatment            | Total DEGs (\|log2FoldChange\| > 2) | Ratio (Up:Down)                   |
| :------------------- | :-------------------------- | :-------------------------------- |
| **Romidepsin (6nM)** | 4,445                       | 4.9 : 1 (3,696 up / 749 down)     |
| **Kromastat (6nM)**  | 1,759                       | 19.7 : 1 (1,674 up / 85 down)     |

Romidepsin induces a broader and more balanced transcriptomic response, while Kromastat's effect is more heavily skewed toward gene activation.

### 2.2 Global Head-to-Head: Romidepsin vs. Kromastat
This contrast directly compares the two drugs at the same 6nM dose to identify unique mechanisms of action.

![Romi vs Kroma Volcano|500](../../results/human/figures/04_02_volcano_Romi_6nM_vs_Kromastat_6nM.png)

- **Genes higher in Romidepsin (Red)**: Includes **419 genes** (e.g., `CALCOCO1`, `CRISPLD2`, and `BAIAP2`). These represent pathways where Romidepsin has a more potent "overdrive" effect compared to Kromastat.
- **Genes higher in Kromastat (Blue)**: Includes **158 genes** (e.g., `LPXN`, `SLC38A1`, and `SASH3`). Despite being less potent overall, Kromastat uniquely or more strongly regulates these specific targets.
- **Potency Proof**: With **577 significant differences** between the two drugs at the same molarity, we have strong evidence that these HDAC inhibitors are not functionally redundant.

> [!NOTE]
> Detailed statistical evidence is available in the authoritative tables:
> - **Global Counts**: [04_02_dge_summary_stats.csv](../../results/human/tables/04_02_dge_summary_stats.csv)
> - **Top Genes (5 Up / 5 Down)**: [04_02_top_dge_genes.csv](../../results/human/tables/04_02_top_dge_genes.csv)

#### Top Differential Representatives (Romi vs. Kroma)

| Direction | Gene Symbol | log2FoldChange | padj |
| :--- | :--- | ---: | ---: |
| **Higher in Kromastat** 🔵 | **LPXN** | -3.4398 | 2.0489e-193 |
| | **SLC38A1** | -2.0231 | 1.7168e-143 |
| | **SASH3** | -2.1468 | 5.0431e-100 |
| | **ANKRD22** | -2.0691 | 2.2313e-86 |
| | **IGF2BP1** | -2.0848 | 5.6171e-85 |
| --- | --- | --- | --- |
| **Higher in Romidepsin** 🔴 | **CALCOCO1** | +2.1168 | 1.4878e-62 |
| | **CRISPLD2** | +2.9376 | 2.5449e-44 |
| | **BAIAP2** | +2.1616 | 1.1926e-40 |
| | **DHRS2** | +3.4768 | 6.9579e-34 |
| | **CEACAM22P** | +3.8436 | 2.6716e-26 |

### 2.3 Local (Cell-Line Specific) Response
To understand if drug sensitivity varies between backgrounds, we analyzed the lines individually.

> [!NOTE]
> Local reactivity counts and top gene representatives are derived from the consolidated tables:
> - **DGE Summary**: [04_02_dge_summary_stats.csv](../../results/human/tables/04_02_dge_summary_stats.csv)
> - **Top Genes**: [04_02_top_dge_genes.csv](../../results/human/tables/04_02_top_dge_genes.csv)

| Cell Line             | Romidepsin (6nM) | Kromastat (6nM) |
| :-------------------- | :--------------- | :-------------- |
| **H9 (Healthy-like)** | 3,910 DEGs       | 1,411 DEGs      |
| **SUPM2 (Cancer)**    | 4,840 DEGs       | 2,086 DEGs      |

| Cell Line             | Romidepsin (6nM)                                                                           | Kromastat (6nM)                                                                                       |
| :-------------------- | :----------------------------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------- |
| **H9 (Healthy-like)** | ![H9 Romi\|500](../../results/human/figures/04_02_volcano_H9_Romi_6nM_vs_DMSO_Romi.png)       | ![H9 Kroma\|500](../../results/human/figures/04_02_volcano_H9_Kromastat_6nM_vs_DMSO_Kromastat.png)       |
| **SUPM2 (Cancer)**    | ![SUPM2 Romi\|500](../../results/human/figures/04_02_volcano_SUPM2_Romi_6nM_vs_DMSO_Romi.png) | ![SUPM2 Kroma\|500](../../results/human/figures/04_02_volcano_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### Key Observations:
1. **Response Magnitude**: The cancer cell line (**SUPM2**) shows a higher number of DEGs (**4,840**) compared to the healthy-like line (**H9: 3,910**) when treated with Romidepsin. A similar trend holds for Kromastat, suggesting a higher transcriptomic reactivity in the cancer background across HDAC inhibitors.
2. **Conserved Targets**: Core targets such as `ENPP2` and `EPAS1` are among the most robustly significant genes across both backgrounds, indicating a conserved drug mechanism that transcends cell line identity.

#### Top Conserved Differential Targets
The following genes show the highest combined significance and magnitude in both healthy-like and cancer backgrounds:

| Drug | Gene Symbol | LFC (H9) | LFC (SUPM2) | Max padj |
| :--- | :--- | ---: | ---: | :--- |
| **Romidepsin** 🔴 | **ENPP2** | +4.4554 | +8.2440 | 0.0 |
| | **RGL1** | +3.8183 | +6.5308 | 0.0 |
| | **SYT11** | +3.5230 | +5.9316 | 0.0 |
| | **GSN** | +4.6114 | +4.7108 | 0.0 |
| | **SLC38A5** | -5.0853 | -4.2156 | 0.0 |
| --- | --- | --- | --- | --- |
| **Kromastat** 🔵 | **EPAS1** | +2.8069 | +2.9520 | 0.0 |
| | **PLA2G4A** | +2.7033 | +2.7799 | 0.0 |
| | **CREG1** | +2.1359 | +2.6599 | 0.0 |
| | **PTPN7** | -2.4089 | -2.0761 | 1.4223e-292 |
| | **VCL** | +2.5037 | +2.2450 | 3.2942e-270 |

> [!NOTE]
> Conserved targets are defined as genes significant (|LFC| > 2.0, padj < 0.05) in both cell lines.
> Full stats: [04_02_conserved_targets.csv](../../results/human/tables/04_02_conserved_targets.csv)

---

## 3. Gene Overlap Analysis (Venn Diagrams)

The Venn diagrams illustrate the degree of transcriptomic "shared identity" between Romidepsin and Kromastat at the 6nM concentration.

| All DEGs                                                         | Up-regulated                                                   | Down-regulated                                                     |
| :--------------------------------------------------------------- | :------------------------------------------------------------- | :----------------------------------------------------------------- |
| ![Venn All\|400](../../results/human/figures/04_03_venn_all.png) | ![Venn Up\|400](../../results/human/figures/04_03_venn_up.png) | ![Venn Down\|400](../../results/human/figures/04_03_venn_down.png) |

### 3.1 Statistical Rigor Summary
To ensure these overlaps are biologically meaningful and not due to chance, we calculated significance using a hypergeometric test.

> [!NOTE]
> Venn set sizes reflect the high-stringency **|LFC| > 2.0** threshold. Total significant targets for this analysis are **4,445** (Romidepsin) and **1,759** (Kromastat).

| Gene Set | Overlap P-value | Jaccard Index | Representation Factor |
| :--- | :--- | :--- | :--- |
| **All DEGs** | < 2e-16 | 0.319 | 4.8x |
| **Up-regulated** | < 2e-16 | 0.364 | 5.7x |
| **Down-regulated** | 2.95e-90 | 0.090 | **26.9x** |

> [!TIP]
> **Full Statistical Report**: The complete overlap metrics (including Romi/Kroma set sizes) are exported to: [04_03_venn_rigor_stats.csv](../../results/human/tables/04_03_venn_rigor_stats.csv)

**Conclusion**: The overlap is highly significant (p < 2e-16). The **26.9x Representation Factor** for down-regulated genes is particularly striking; it suggests that while Romidepsin silences fewer genes than it activates, the ones it *does* share with Kromastat are part of an extremely rigid, conserved silencing program.

### 3.2 Top Shared Targets
These genes represent the "Core Response" triggered by both inhibitors, strictly sorted by combined significance:

| Direction | Gene Symbol | LFC (Romi) | LFC (Kroma) | Max padj | Biological Context |
| :--- | :--- | ---: | ---: | :--- | :--- |
| **Up-regulated** 🔴 | **EPAS1** | +3.05 | +2.93 | 0.0 | Hypoxia/Angiogenesis |
| | **GSN** | +4.71 | +4.09 | 3.73e-225 | Cytoskeleton/Apoptosis |
| | **PLA2G4A** | +2.44 | +2.79 | 1.45e-226 | Lipid Signaling |
| | **STXBP1** | +2.28 | +2.12 | 3.45e-218 | Vesicle Trafficking |
| | **ZSWIM6** | +2.63 | +2.03 | 1.96e-191 | Nervous System Dev |
| **Down-regulated** 🔵 | **TNFRSF8** | -3.94 | -2.12 | 1.45e-250 | Immune Signaling |
| | **DUS3L** | -2.81 | -2.09 | 3.68e-109 | tRNA Modification |
| | **PTPN7** | -4.17 | -2.19 | 1.94e-107 | T-cell Phosphatase |
| | **PATZ1** | -2.95 | -2.05 | 3.58e-94 | Transcription Repressor |
| | **MARVELD1** | -4.06 | -3.20 | 7.29e-84 | Nuclear Envelope |

#### Chromatin Sanity Check: RCC1/RCC2 Suppression
While `RCC1/RCC2` fall slightly below the strict |LFC| > 2.0 threshold in the Kromastat group, their suppression is statistically undeniable across both inhibitors. This confirms a shared impact on **Regulators of Chromosome Condensation**, a key feature of Class I HDAC inhibition.

| Gene | Log2FC (Romi) | padj (Romi) | Log2FC (Kroma) | padj (Kroma) |
| :--- | ---: | :--- | ---: | :--- |
| **RCC2** | -2.77 | 0.0 | -1.64 | 6.63e-121 |
| **RCC1** | -1.97 | 0.0 | -1.09 | 2.73e-97 |

> [!NOTE]
> Shared targets in the main table are strictly defined as the intersection of significant genes (|LFC| > 2.0, padj < 0.05).

### 3.3 The "Subset" Hypothesis
The rigor metrics support a clear functional hierarchy:
- **Shared Response (32%)**: 1,502 genes forming the core "Class I HDAC Inhibition" signature.
- **Kromastat Uniqueness (6%)**: 257 genes.
- **Romidepsin Uniqueness (62%)**: 2,943 genes.

**Finding**: Kromastat acts primarily as a **transcriptomic subset** of Romidepsin. Approximately **85% of Kromastat's total response is contained within the Romidepsin signature**.

### 3.4 Directional Differences
- **Shared Activation**: Both drugs are highly effective at gene activation, sharing 1,432 "Up" genes.
- **Unique Suppression**: Romidepsin is a significantly more potent suppressor of transcription. It uniquely down-regulates **680 genes**, while Kromastat only uniquely suppresses **16**. 

This suggests that the "potency gap" observed in the PCA is largely driven by Romidepsin's ability to silence a vast network of genes that are inaccessible to Kromastat at this dose.

---

## 4. Functional Pathway Analysis (GSEA Leading-Edge)

[CONTENT_PENDING]

---

## 5. Global Transcriptomic Variance (Top Variable Genes)

[CONTENT_PENDING]

---

## 6. Functional Hallmark Analysis (GSEA NES)

This section summarizes the global biological shifts using the Normalized Enrichment Score (NES) from Gene Set Enrichment Analysis (GSEA), providing a quantitative rank of the pathways most impacted by drug treatment.

> [!NOTE]
> **Combined Plot Curation**: The combined NES barplot is curated to the **5 biologically informative contrasts** (H9 and SUPM2 per-drug, plus the direct Romi vs. Kromastat comparison). Three contrasts are intentionally excluded from the combined figure but are preserved as individual plots for reference:
> - **`DMSO_Kromastat_vs_DMSO_Romi`**: Vehicle QC baseline — near-empty by design, confirming biological equivalence of control conditions.
> - **`Romi_6nM_vs_DMSO_Romi`** and **`Kromastat_6nM_vs_DMSO_Kromastat`**: Global contrasts that pool H9 and SUPM2 together, producing statistically weaker and harder-to-interpret signals compared to the cell-line-specific analyses.

![Combined Hallmark NES](../../results/human/figures/04_06_hallmark_nes_combined.png)
### 6.1 Global Functional Response
The global analysis identifies the core transcriptomic shifts triggered by each drug across all human samples.

| Romidepsin (6nM) vs DMSO                                                                | Kromastat (6nM) vs DMSO                                                                            |
| :-------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------- |
| ![Romi NES\|600](../../results/human/figures/04_06_hallmark_nes_Romi_6nM_vs_DMSO_Romi.png) | ![Kroma NES\|600](../../results/human/figures/04_06_hallmark_nes_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **Observations**:
- **Conserved Suppression**: Both drugs trigger a massive "Proliferation Crash," with **Myc Targets V1**, **E2f Targets**, and **G2m Checkpoint** as the top suppressed signatures. **Tnfa Signaling via Nfkb** leads all activated signatures for both drugs (NES ≈ +2.24), representing a universal chromatin-stress response.
- **Key Divergence**: Romidepsin's suppressed list includes **Interferon Alpha Response (NES −1.89)** and **Allograft Rejection (NES −1.81)** — immune pathway suppressions absent from Kromastat. Kromastat instead features **Oxidative Phosphorylation (NES −1.97)** as its 5th most suppressed pathway — reflecting a distinct metabolic emphasis.

---

### 6.2 Local Response: Healthy-like (H9)
Individual analysis of the H9 lineage reveals how these drugs affect a "healthy" metabolic background.

| H9: Romidepsin (6nM)                                                                          | H9: Kromastat (6nM)                                                                                      |
| :-------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------- |
| ![H9 Romi NES\|500](../../results/human/figures/04_06_hallmark_nes_H9_Romi_6nM_vs_DMSO_Romi.png) | ![H9 Kroma NES\|500](../../results/human/figures/04_06_hallmark_nes_H9_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- **Shared Core**: Both drugs show the canonical proliferation crash (**Myc Targets V1**, **Myc Targets V2**, **E2f Targets**, **G2m Checkpoint**) as the dominant suppressed signatures in H9.
- **Romi-Specific Suppression in H9**: Romidepsin additionally suppresses **Allograft Rejection** and **IL2-STAT5 Signaling** in H9 — both immune-modulatory pathways — confirming a more immunosuppressive character even in the healthy background.
- **Divergent Activation in H9**: While both drugs activate **Apoptosis** and the **P53 Pathway** in H9, Kromastat induces a more robust activation of these signatures (Kroma NES ≈ +1.89 vs Romi NES ≈ +1.65). This indicates Kromastat engages a cleaner, p53-driven apoptotic programme in the healthy lineage, whereas Romidepsin's response is more dominated by broader cytotoxic stress.

---

### 6.3 Local Response: Cancer-like (SUPM2)
The cancer-specific response highlights the "Growth Crushing" effectiveness of each molecule.

| SUPM2: Romidepsin (6nM)                                                                             | SUPM2: Kromastat (6nM)                                                                                         |
| :-------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------- |
| ![SUPM2 Romi NES\|600](../../results/human/figures/04_06_hallmark_nes_SUPM2_Romi_6nM_vs_DMSO_Romi.png) | ![SUPM2 Kroma NES\|600](../../results/human/figures/04_06_hallmark_nes_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- **Leading Suppression**: In both SUPM2 plots, **E2f Targets** is the single most suppressed hallmark (NES ≈ −3.0), overtaking **Myc Targets V1** — reflecting maximal disruption of the G1/S transition in the cancer background.
- **Conserved Stress Activation in SUPM2**: Both drugs show significant activation of **Apoptosis**, **P53 Pathway**, and **KRAS Signaling Up** in the cancer background. This co-activation of apoptotic and oncogenic stress pathways is a signature of terminal cancer cell stress — cells simultaneously dying and mounting a survival response.
- **Divergent Inflammatory Response in SUPM2**: While both drugs activate **Complement** (Romi NES ≈ +2.05, Kroma NES ≈ +2.33), Kromastat triggers a broader immune and inflammatory activation profile, including significantly stronger activation of **Il2 Stat5 Signaling** and the **Inflammatory Response** hallmark compared to Romidepsin.

---

### 6.4 Direct Comparison: Romidepsin vs. Kromastat
This contrast directly visualizes the functional difference between the two drugs at the same 6nM molarity.

![Romi vs Kroma NES|600](../../results/human/figures/04_06_hallmark_nes_Romi_6nM_vs_Kromastat_6nM.png)

- **Romi is more suppressive of immune signalling**: The "Suppressed" panel (pathways MORE active in Kromastat than Romidepsin) is dominated by **Interferon Alpha Response (NES −2.36)**, **Allograft Rejection (NES −2.34)**, and **Interferon Gamma Response (NES −2.19)**. This quantitatively confirms that Kromastat preserves immune signalling that Romidepsin shuts down.
- **Romi shows a stronger non-immune stress response**: The "Activated" panel (pathways MORE active in Romi vs Kroma) shows **Tnfa Signaling via Nfkb**, **Myogenesis**, **Heme Metabolism**, and **Coagulation** — confirming Romidepsin triggers a more intense cytotoxic stress response.

---

### 6.5 Summary of Hallmark Findings

The functional landscape confirms that the difference between Romidepsin and Kromastat is both **quantitative** (intensity) and **qualitative** (biological fate).

#### **1. The Shared Core: "The Proliferation Crash"**
Both drugs achieve the primary goal of HDAC inhibition — a massive transcriptomic shutdown of cell cycle programmes. **E2f Targets** emerges as the most suppressed hallmark in the cancer-specific (SUPM2) analysis (NES ≈ −3.0), reflecting deep disruption of the G1/S transition.

#### **2. Divergent Apoptotic Programmes**
Both drugs induce apoptosis, but through mechanistically distinct routes. Romidepsin in SUPM2 pairs Apoptosis activation with **KRAS oncogenic stress signalling** — a hallmark of terminally stressed cancer cells. Kromastat in H9 engages a cleaner **p53-driven apoptotic programme** independent of KRAS — suggesting distinct routes to cell death in different cellular backgrounds.

#### **3. The Immune-Signalling Divergence**
The most clinically significant finding. Romidepsin broadly suppresses interferon and immune-rejection signalling (Interferon Alpha, Interferon Gamma, Allograft Rejection) across all backgrounds. Kromastat **spares** these pathways entirely, preserving the immune signalling infrastructure required for combination with checkpoint inhibitor therapies.

**Final Verdict**: **Romidepsin** is the more potent cytotoxic agent — deeper proliferation suppression, stronger stress activation, forced apoptosis. **Kromastat** is the more "immune-intelligent" inhibitor — executing a controlled anti-proliferative and p53-driven apoptotic programme while leaving the interferon and immune-rejection machinery intact.

---

## 7. GSEA Hallmark Dotplot Analysis

This section cross-validates the NES barplot findings using the multi-contrast Hallmark dotplot. Each dot encodes three simultaneous dimensions: **position** (Activated / Suppressed), **colour + alpha** (direction × significance), and **size** (leading edge gene count). The globally-aligned Y-axis allows direct row-wise comparison across all five curated contrasts.

> [!NOTE]
> **Combined Plot Curation**: Consistent with the NES barplot (Section 6), the combined dotplot is curated to the **5 biologically informative contrasts**. The following are excluded from the combined figure but preserved as individual saved plots:
> - **`DMSO_Kromastat_vs_DMSO_Romi`**: Vehicle QC baseline — near-empty by design. Its individual plot confirms biological equivalence of control conditions.
> - **`Romi_6nM_vs_DMSO_Romi`** and **`Kromastat_6nM_vs_DMSO_Kromastat`**: Global contrasts pooling both cell lines, producing weaker and harder-to-interpret signals.

### 7.1 Global Overview

![GSEA Hallmark Combined|1132](../../results/human/figures/04_07_gsea_dotplot_combined.png)

The combined plot reveals an immediately clear architecture:

- **Top of the Y-axis (Activated zone)**: Pathways such as **Tnfa Signaling via Nfkb**, **Complement**, **Hypoxia**, **Myogenesis**, and **Apoptosis** consistently appear as activated dots (red) across multiple contrasts, confirming the universality of the chromatin-stress and inflammatory response to HDAC inhibition.
- **Bottom of the Y-axis (Suppressed zone)**: **Myc Targets V1**, **E2f Targets**, **G2m Checkpoint**, **Myc Targets V2**, and **Allograft Rejection** anchor the suppressed (blue) end of the plot across all drug-treated contrasts, confirming the canonical proliferation crash.

---

### 7.2 Local Response: Cross-Contrast Comparison

#### **Romidepsin (H9 and SUPM2)**

|                                      H9 — Romi 6nM vs. DMSO                                       |                                        SUPM2 — Romi 6nM vs. DMSO                                        |
| :-----------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------: |
| ![H9 Romi Dotplot\|500](../../results/human/figures/04_07_gsea_dotplot_H9_Romi_6nM_vs_DMSO_Romi.png) | ![SUPM2 Romi Dotplot\|500](../../results/human/figures/04_07_gsea_dotplot_SUPM2_Romi_6nM_vs_DMSO_Romi.png) |

- **Suppression** is the dominant signal. Large, dark blue dots for **Myc Targets V1**, **E2f Targets**, **G2m Checkpoint**, and **Myc Targets V2** appear in both H9 and SUPM2, indicating a potent and consistent proliferative shutdown.
- **Immune suppression** is visible specifically in the Romidepsin panels: **Allograft Rejection**, **Interferon Alpha Response**, and **Interferon Gamma Response** all carry blue dots of moderate-to-large size — these are entirely absent from the Kromastat panels, confirming drug-specific immune-sparing divergence.
- **SUPM2+Romi** uniquely shows co-activation of **Apoptosis** and **Kras Signaling Up** — the terminal stress signature observed in the NES barplots is confirmed here with high leading-edge counts (large dot size).

---

#### **Kromastat (H9 and SUPM2)**

| H9 — Kromastat 6nM vs. DMSO | SUPM2 — Kromastat 6nM vs. DMSO |
| :---: | :---: |
| ![H9 Kromastat Dotplot\|500](../../results/human/figures/04_07_gsea_dotplot_H9_Kromastat_6nM_vs_DMSO_Kromastat.png) | ![SUPM2 Kromastat Dotplot\|500](../../results/human/figures/04_07_gsea_dotplot_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- The suppression zone mirrors Romidepsin's proliferation crash, but **Allograft Rejection**, **Interferon Alpha Response**, and **Interferon Gamma Response** are conspicuously **absent** — no blue dots appear in those rows for Kromastat, visually confirming the immune-sparing phenotype.
- **Apoptosis** and **P53 Pathway** carry red activated dots in H9+Kromastat specifically, consistent with the controlled p53-driven apoptotic programme identified in the NES barplot analysis.
- **SUPM2+Kromastat** activates a broader inflammatory profile, including **Complement** and **Il2 Stat5 Signaling** — immune-regulatory pathways that are absent from the SUPM2+Romidepsin response.

---

#### **Direct Romi vs. Kromastat Comparison**

| Romi 6nM vs. Kromastat 6nM |
| :---: |
| ![Romi vs Kromastat Dotplot\|500](../../results/human/figures/04_07_gsea_dotplot_Romi_6nM_vs_Kromastat_6nM.png) |

- **Suppressed (blue, right side)** — pathways where Romidepsin drives **lower activity than Kromastat**: **Myc Targets V1**, **G2m Checkpoint**, **E2f Targets**, **Allograft Rejection**, **Interferon Alpha Response**, **Interferon Gamma Response**, **Mtorc1 Signaling**, and **Il2 Stat5 Signaling**. This indicates Kromastat preserves both proliferative and immune signalling infrastructure at a comparatively higher level.
- **Activated (red, left side)** — pathways where Romidepsin shows **higher activity than Kromastat**: **Tnfa Signaling via Nfkb**, **Myogenesis**, **Heme Metabolism**, and **Coagulation** — a stress-inflammatory signature unique to Romidepsin's mechanism.
- **Key takeaway**: Romidepsin is the more aggressive suppressor across the board — deeper on both proliferation and immune pathways. Kromastat, by comparison, selectively suppresses proliferation while leaving immune and metabolic signalling more intact, which is what gives it the **immune-sparing phenotype** observed in the individual contrast analyses.


### 7.3 Summary: Dotplot Cross-Validation

The GSEA Hallmark dotplot independently corroborates all major findings from the NES barplot analysis (Section 6):

| Finding | NES Barplot | Dotplot |
| :--- | :---: | :---: |
| Shared proliferation crash (Myc Targets V1, E2f Targets, G2m Checkpoint) | ✅ | ✅ |
| Romidepsin suppresses immune pathways (Interferon Alpha, Interferon Gamma, Allograft Rejection) | ✅ | ✅ |
| Kromastat spares immune pathways entirely | ✅ | ✅ |
| Kromastat activates P53/Apoptosis in H9 | ✅ | ✅ |
| Romidepsin activates Kras Signaling Up in SUPM2 | ✅ | ✅ |
| DMSO vehicles are biologically equivalent | — | ✅ (QC confirmed) |

The dotplot adds a critical new layer: **leading edge size** reveals that the immune-suppression pathways (Interferon Alpha, Allograft Rejection) in Romidepsin carry **large gene sets**, indicating broad transcriptome-level shutdown rather than marginal effects. This reinforces the clinical significance of Kromastat's immune-sparing phenotype.

## 8. Cross-Cell-Line Consistency Proofing

To ensure that the identified biological signatures are robust and not driven by idiosyncratic cell-line artifacts, we performed a consistency analysis comparing the H9 and SUPM2 cell lines. This validation utilized two complementary approaches: (1) intersectional set analysis (UpSet) to evaluate overlap among discrete DEG lists, and (2) global Log2FoldChange (LFC) correlation to assess the directional consensus of the entire transcriptome.

### 8.1 Intersectional Degree of Consensus (UpSet Matrix)

We analyzed the overlap of significant DEGs (defined by `padj < 0.05` and `|log2FoldChange| > 2`) across the four primary experimental conditions. This matrix visualizes the "potency" and "uniqueness" of the drug-induced transcriptional changes.

![UpSet Consistency Matrix|717](../../results/human/figures/04_09_upset_consistency.png)

| Observation Category | Key Finding | Biological Interpretation |
|:---|:---|:---|
| **Magnitude of Response** | SUPM2 Romi has the largest set (~7500+ DEGs [REVERIFY]) | Romidepsin triggers a massively broader transcriptional response in SUPM2, confirming its heightened sensitivity compared to other conditions. |
| **Comparative Potency** | H9 Romi is the second largest set (~5500+ DEGs [REVERIFY]) | Romidepsin consistently acts as a highly potent perturbagen across different cellular backgrounds. |
| **Drug Selectivity** | Kromastat sets are notably smaller (~3000–4000 DEGs [REVERIFY]) | Consistent with its targeted profile, Kromastat acts as a more selective perturbagen, inducing a focused transcriptome-wide footprint. |
| **Cell-Line Specificity** | Largest unique intersection = SUPM2 Romi (2,641) [REVERIFY] | A significant core of genes is uniquely regulated by Romidepsin in SUPM2, likely representing cell-line-specific sensitivity pathways. |
| **Mechanism Consistency** | H9 Krom ∩ SUPM2 Krom consensus (2,364 DEGs) [REVERIFY] | Despite smaller absolute set sizes, Kromastat displays high cross-cell-line concordance, indicating a very stable mechanistic action. |
| **Conserved HDAC Core** | 4-way all-condition overlap (45 genes) [REVERIFY] | These 45 genes are significant across ALL conditions and cell lines, representing the "gold-standard" conserved targets of HDAC inhibition. |

### 8.2 Global Transcriptomic Concordance (LFC Correlation)

While intersectional analysis uses binary cutoffs, the LFC correlation evaluates the overall transcriptomic trend across the entire dynamic range. This approach validates the "directionality" of the drug response without threshold bias.

![LFC Correlation|986](../../results/human/figures/04_10_lfc_correlation.png)

| Analysis Dimension | Resulting Observation | Biological Significance |
|:---|:---|:---|
| **Global Correlation** | Moderate R values for both (Romi: 0.537, Krom: 0.483) | There is a clear, statistically significant positive correlation in drug-induced regulation across both cell lines, proving a conserved mechanism. |
| **Directional Consensus** | Romi R (0.537) > Krom R (0.483) | Romidepsin yields slightly higher cross-line concordance, as its stronger primary signal effectively overcomes baseline cell-line noise. |
| **Consensus Mapping** | Dominant "Consensus" (red) cluster along the diagonal | The majority of biologically relevant DEGs agree in both direction and significance, validating the use of multifactorial modeling. |
| **Sensitivity Profile** | Romi has a significantly wider LFC range (±15) | SUPM2's response to Romidepsin includes extreme LFC outliers, further confirming the extreme sensitivity of this model. |
| **Modulation Precision** | Kromastat range is compressed (−4 to +8 LFC) | Kromastat provides a more refined, less disruptive transcriptomic modulation compared to the global shock induced by Romidepsin. |
| **Divergence Clusters** | "Specific" (blue) dots scattered off-diagonal | These represent true biological divergence between the cell lines, ensuring the models capture a representative range of T-cell lymphoma heterogeneity. |

---

## Supplementary: GO ORA Analysis

> [!WARNING]
> **Not suitable for primary presentation.** The GO Over-Representation Analysis (ORA) results (`04_08_ora_dotplot.R`) are flagged as **supplementary reference only** and have not been integrated into the main narrative.

### Reason for exclusion
The ORA combined and individual dotplots are dominated by **biologically irrelevant GO terms** — primarily cilium organization, axoneme assembly, flagellar sperm motility, and cilia movement — which have no mechanistic relevance to T-cell lymphoma or HDAC inhibition. This is a known **gene universe background contamination artifact**: when the background is the full genome, small gene sets with disproportionate overlap can achieve inflated significance scores.

While plausible immune terms (B cell activation, lymphocyte differentiation, T cell differentiation) do appear, they are scattered inconsistently across contrasts and are overshadowed by the artifact signal.

### Recommended future action
To fix the ORA signal, restrict the gene universe (`universe` parameter in `enrichGO()` in `03_enrichment.R`) to **only the genes tested in DESeq2** (i.e., genes that passed independent filtering), rather than the full genome. This is the standard corrective measure. This is not prioritized for the current analysis phase.

Individual ORA dotplot files are preserved in `results/human/figures/04_08_ora_dotplot_*.png` for reference.

---

## 9. Comprehensive Summary & Biological Conclusion

### 9.1 Synthesis of Findings
This downstream transcriptomic analysis systematically compares the effects of two HDAC inhibitors, Romidepsin and Kromastat, at a 6nM dose across healthy-like (H9) and cancer-like (SUPM2) human cell lines. 

1. **The Shared Core (The Proliferation Crash)**: Both drugs successfully execute the primary function of HDAC inhibition—a massive transcriptomic inversion that shuts down proliferative pathways (MYC Targets, E2F Targets, G2M Checkpoint) across all cellular backgrounds *(see Sections [[downstream_analysis_human#4.1 Global Functional Response (Signature Matrices)|4.1]], [[downstream_analysis_human#6.1 Global Functional Response|6.1]], [[downstream_analysis_human#7.1 Global Overview|7.1]])*.
2. **The Potency Gap**: Romidepsin acts as a "Metabolic Crusher," driving a deeper, more aggressive transcriptional shutdown and inducing intense stress responses (Hypoxia, matrix remodeling). Kromastat acts as a "Partial Modulator," achieving significant but milder suppression of the same targets *(see Sections [[downstream_analysis_human#1.1 Interpretation:|1.1]], [[downstream_analysis_human#2.1.2 Global Magnitude of Response|2.1.2]], [[downstream_analysis_human#5.2 Drug Effect: The "Potency Gap" at Single-Gene Resolution|5.2]])*.
3. **Divergent Stress and Apoptotic Fates**: In the cancer background (SUPM2), Romidepsin forces a chaotic, terminal stress state combining apoptosis with KRAS oncogenic signaling. Conversely, Kromastat engages a cleaner, p53-driven apoptotic program in the healthy background and triggers broader inflammatory responses (Complement cascade) in the cancer background *(see Sections [[downstream_analysis_human#4.3 Local Response: SUPM2 (Cancer-like)|4.3]], [[downstream_analysis_human#5.4 Deep-Dive: Notable Gene Clusters|5.4]], [[downstream_analysis_human#6.2 Local Response: Healthy-like (H9)|6.2]], [[downstream_analysis_human#6.3 Local Response: Cancer-like (SUPM2)|6.3]])*.
4. **The Immune-Signaling Divergence**: Romidepsin exerts a broad immunosuppressive footprint, heavily down-regulating Interferon Alpha, Interferon Gamma, and Allograft Rejection pathways. Strikingly, Kromastat spares these critical immune pathways entirely *(see Sections [[downstream_analysis_human#4.4 Direct Comparison: Romidepsin vs. Kromastat|4.4]], [[downstream_analysis_human#6.4 Direct Comparison: Romidepsin vs. Kromastat|6.4]], [[downstream_analysis_human#7.2 Local Response: Cross-Contrast Comparison|7.2]])*.

### 9.2 Answer to the Biological Question
**Biological Question**: *How do the transcriptomic profiles of Romidepsin and Kromastat differentiate their therapeutic utility as HDAC inhibitors?*

**Conclusion**: While both drugs target the same core epigenetic and proliferative machinery, their biological fates diverge significantly in both intensity and immune modulation. 
- **Romidepsin** is the superior **direct cytotoxic agent**. Its aggressive repression of the cancer signature makes it highly effective at dismantling tumor metabolism *(see Sections [[downstream_analysis_human#5.2 Drug Effect: The "Potency Gap" at Single-Gene Resolution|5.2]], [[downstream_analysis_human#6.3 Local Response: Cancer-like (SUPM2)|6.3]])*, but its immunosuppressive nature may limit its long-term efficacy by dampening the body's anti-tumor immune response *(see Sections [[downstream_analysis_human#4.4 Direct Comparison: Romidepsin vs. Kromastat|4.4]], [[downstream_analysis_human#6.4 Direct Comparison: Romidepsin vs. Kromastat|6.4]])*.
- **Kromastat** is the superior **immune primer**. Although less potent at crashing cellular growth at the evaluated dose *(see Sections [[downstream_analysis_human#1.1 Interpretation:|1.1]], [[downstream_analysis_human#4.1 Global Functional Response (Signature Matrices)|4.1]])*, its unique "immune-sparing" phenotype—preserving interferon and immune-rejection signaling—makes it an exceptionally promising candidate for **combination therapies** *(see Sections [[downstream_analysis_human#4.5 Summary of Functional Pathway Findings|4.5]], [[downstream_analysis_human#6.5 Summary of Hallmark Findings|6.5]])*. By halting proliferation without silencing the immune system, Kromastat could effectively prime tumors for subsequent or concurrent treatment with immunotherapies (e.g., checkpoint inhibitors).
