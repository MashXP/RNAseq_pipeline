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

| Treatment            | Total DEGs (\|log2FC\| > 1) | Ratio (Up:Down)                   |
| :------------------- | :-------------------------- | :-------------------------------- |
| **Romidepsin (6nM)** | 9,594                       | 2.5 : 1 (6,867 up / 2,727 down)   |
| **Kromastat (6nM)**  | 5,131                       | 4.6 : 1 (4,222 up / 909 down)     |

Romidepsin induces a broader and more balanced transcriptomic response, while Kromastat's effect is more heavily skewed toward gene activation.

### 2.2 Global Head-to-Head: Romidepsin vs. Kromastat
This contrast directly compares the two drugs at the same 6nM dose to identify unique mechanisms of action.

![Romi vs Kroma Volcano|500](../../results/human/figures/04_02_volcano_Romi_6nM_vs_Kromastat_6nM.png)

- **Genes higher in Romidepsin (Red)**: Includes `SPAG9`, `GM2A`, and `SLC25A4`. These represent pathways where Romidepsin has a more potent "overdrive" effect compared to Kromastat.
- **Genes higher in Kromastat (Blue)**: Includes `STAT6`, `HDGF`, and `TNFRSF8`. Despite being less potent overall, Kromastat uniquely or more strongly regulates these specific targets.
- **Potency Proof**: With **3,084 significant differences** between the two drugs at the same molarity, we have strong evidence that these HDAC inhibitors are not functionally redundant.

### 2.3 Local (Cell-Line Specific) Response
To understand if drug sensitivity varies between backgrounds, we analyzed the lines individually.

> [!NOTE]
> Local reactivity counts are derived from individual cell-line DGE tables (e.g., `02_dge_H9_*.csv`) located in: [results/human/tables/](../../results/human/tables/).

| Cell Line             | Romidepsin (6nM) | Kromastat (6nM) |
| :-------------------- | :--------------- | :-------------- |
| **H9 (Healthy-like)** | 8,258 DEGs       | 4,499 DEGs      |
| **SUPM2 (Cancer)**    | 9,561 DEGs       | 5,742 DEGs      |

| Cell Line             | Romidepsin (6nM)                                                                           | Kromastat (6nM)                                                                                       |
| :-------------------- | :----------------------------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------- |
| **H9 (Healthy-like)** | ![H9 Romi\|500](../../results/human/figures/04_02_volcano_H9_Romi_6nM_vs_DMSO_Romi.png)       | ![H9 Kroma\|500](../../results/human/figures/04_02_volcano_H9_Kromastat_6nM_vs_DMSO_Kromastat.png)       |
| **SUPM2 (Cancer)**    | ![SUPM2 Romi\|500](../../results/human/figures/04_02_volcano_SUPM2_Romi_6nM_vs_DMSO_Romi.png) | ![SUPM2 Kroma\|500](../../results/human/figures/04_02_volcano_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### Key Observations:
1. **Response Magnitude**: The cancer cell line (**SUPM2**) shows a higher number of DEGs (**9,561**) compared to the healthy-like line (**H9: 8,258**) when treated with Romidepsin. A similar trend holds for Kromastat, suggesting a higher transcriptomic reactivity in the cancer background across HDAC inhibitors.
2. **Conserved Targets**: Despite the difference in scale, core targets such as `PMEL` and `CHGB` are highly significant across both lines, indicating a conserved drug mechanism that transcends the cell line identity.

---

## 3. Gene Overlap Analysis (Venn Diagrams)

The Venn diagrams illustrate the degree of transcriptomic "shared identity" between Romidepsin and Kromastat at the 6nM concentration.

| All DEGs                                                         | Up-regulated                                                   | Down-regulated                                                     |
| :--------------------------------------------------------------- | :------------------------------------------------------------- | :----------------------------------------------------------------- |
| ![Venn All\|400](../../results/human/figures/04_03_venn_all.png) | ![Venn Up\|400](../../results/human/figures/04_03_venn_up.png) | ![Venn Down\|400](../../results/human/figures/04_03_venn_down.png) |

### 3.1 Statistical Rigor Summary
To ensure these overlaps are biologically meaningful and not due to chance, we calculated significance using a hypergeometric test.

> [!NOTE]
> Venn set sizes below are from the `04_03_venn_rigor_stats.csv` snapshot. Current DGE tables yield **9,593** (Romi) and **5,697** (Kroma) at the same thresholds; the proportional conclusions are unchanged.

| Gene Set | Overlap P-value | Jaccard Index | Representation Factor |
| :--- | :--- | :--- | :--- |
| **All DEGs** | < 2.2e-16 | 0.415 | 2.2x |
| **Up-regulated** | < 2.2e-16 | 0.479 | 3.1x |
| **Down-regulated** | < 2.2e-16 | 0.248 | **7.2x** |

> [!TIP]
> **Full Statistical Report**: The complete overlap metrics (including Romi/Kroma set sizes) are exported to: [04_03_venn_rigor_stats.csv](../../results/human/tables/04_03_venn_rigor_stats.csv)

**Conclusion**: The overlap is highly significant (p < 2.2e-16). The **7.2x Representation Factor** for down-regulated genes is particularly striking; it suggests that while Romidepsin silences more genes, the ones it *does* share with Kromastat are part of an extremely rigid, conserved silencing program.

### 3.2 Top Shared Targets
These genes represent the "Core Response" triggered by both inhibitors:

| Direction | Top Shared Genes (padj ~ 0) | Biological Context |
| :--- | :--- | :--- |
| **Up-regulated** | `EPAS1`, `GM2A`, `GSN`, `ABCD1` | Angiogenesis, Lipid Transport |
| **Down-regulated** | `RCC1`, `RCC2`, `HDGF`, `TNFRSF8` | **Chromatin Condensation**, Growth Factors |

The shared suppression of **RCC1/RCC2** (Regulators of Chromosome Condensation) is a perfect "sanity check" for an HDAC inhibitor study, as these drugs directly impact chromatin structure.

### 3.3 The "Subset" Hypothesis
The rigor metrics support a clear functional hierarchy:
- **Shared Response (42%)**: 4,322 genes forming the core "Class I HDAC Inhibition" signature.
- **Kromastat Uniqueness (8%)**: Only 809 genes.
- **Romidepsin Uniqueness (51%)**: 5,272 genes.

**Finding**: Kromastat acts primarily as a **transcriptomic subset** of Romidepsin. Approximately **85% of Kromastat's total response is contained within the Romidepsin signature**.

### 3.4 Directional Differences
- **Shared Activation**: Both drugs are highly effective at gene activation, sharing 3,593 "Up" genes.
- **Unique Suppression**: Romidepsin is a significantly more potent suppressor of transcription. It uniquely down-regulates **2,005 genes**, while Kromastat only uniquely suppresses **187**. 

This suggests that the "potency gap" observed in the PCA is largely driven by Romidepsin's ability to silence a vast network of genes that are inaccessible to Kromastat at this dose.

---

## 4. Functional Pathway Analysis (GSEA Leading-Edge)

Pathway analysis connects gene lists to biological phenomena. Both drugs induce a massive suppression of the cell cycle machinery.

### 4.1 Global Functional Response (Signature Matrices)

| Romidepsin vs DMSO                                                                             | Kromastat vs DMSO                                                                                         |
| :--------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------- |
| ![Romi_Heatmap\|825](../../results/human/figures/04_04_heatmap_pathway_Romi_6nM_vs_DMSO_Romi.png) | ![Kroma_Heatmap\|825](../../results/human/figures/04_04_heatmap_pathway_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **Romidepsin (Global)**
| Pathway            |   D_R    | R_6nM |   D_K    |        K_6nM         |
| :----------------- | :------: | :---: | :------: | :------------------: |
| **MYC Targets V1** | light 🟥 |  🟦   | light 🟥 | generally lighter 🟦 |
| **E2F Targets**    |    //    |  🟦   |    //    |          //          |
| **G2M Checkpoint** |    //    |  🟦   |    //    |          //          |
| **MYC Targets V2** |    //    |  🟦   |    //    |          //          |
| **TNFa via NFkB**  | light 🟦 |  🟥   | light 🟦 | generally lighter 🟥 |
| **Myogenesis**     |    //    |  🟥   |    //    |          //          |

#### **Kromastat (Global)**
| Pathway            |   D_R    | R_6nM |   D_K    |        K_6nM         |
| :----------------- | :------: | :---: | :------: | :------------------: |
| **MYC Targets V1** | light 🟥 |  🟦   | light 🟥 | generally lighter 🟦 |
| **MYC Targets V2** |    //    |  🟦   |    //    |          //          |
| **E2F Targets**    |    //    |  🟦   |    //    |          //          |
| **G2M Checkpoint** |    //    |  🟦   |    //    |          //          |
| **TNFa via NFkB**  | light 🟦 |  🟥   | light 🟦 | generally lighter 🟥 |
| **EMT**            |    //    |  🟥   |    //    |          //          |

#### **Interpreting "The Big Flip"**
The matrices reveal a striking **Reciprocal Signature** between the control state and the drug state:
1. **The Baseline Baseline**: In both control groups (`D_R`, `D_K`), the cells exhibit a mild "Light Red" signature in metabolic pathways (MYC) and a "Light Blue" signature in stress pathways (TNFa). This represents a basal state of active proliferation and low stress.
2. **The Full Flip (Romidepsin)**: Romidepsin treatment completely inverts this biology, turning proliferative pathways into **Deep Blue** (Absolute Suppression) and stress pathways into **Deep Red** (Maximum Activation).
3. **The Partial Potency (Kromastat)**: While Kromastat attempts the same "flip," its signature in the `K_6nM` column is **generally lighter**. This provides visual confirmation of the "Potency Gap" identified in the PCA—Kromastat hits the same targets but with less repressive force at the same molarity.

---

### 4.2 Local Response: Healthy-like (H9)

| H9: Romidepsin (6nM)                                                                         | H9: Kromastat (6nM)                                                                                     |
| :------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| ![H9 Romi\|825](../../results/human/figures/04_04_heatmap_pathway_H9_Romi_6nM_vs_DMSO_Romi.png) | ![H9 Kroma\|825](../../results/human/figures/04_04_heatmap_pathway_H9_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **H9 Local Matrix (Healthy-like)**
##### R vs D

| Pathway            |   D_R    | R_6nM |   D_K    |                 K_6nM                  |
| :----------------- | :------: | :---: | :------: | :------------------------------------: |
| **MYC Targets V1** | light 🟥 |  🟦   | light 🟥 | generally lighter 🟦, with a bit of 🟥 |
| **MYC Targets V2** |    //    |  🟦   |    //    |                   //                   |
| **E2f Targets**    |    //    |  🟦   |    //    |                   //                   |
| **G2M Checkpoint** |    //    |  🟦   |    //    |                   //                   |
| **Allograft Rej.** |    //    |  🟦   |    //    |             more light 🟥              |
| **Myogenesis**     | light 🟦 |  🟥   | light 🟦 |                   //                   |
##### K vs D

| Pathway            |   D_R    | R_6nM |   D_K    |        K_6nM         |
| :----------------- | :------: | :---: | :------: | :------------------: |
| **MYC Targets V1** | light 🟥 |  🟦   | light 🟥 | generally lighter 🟦 |
| **MYC Targets V2** |    //    |  🟦   |    //    |          //          |
| **E2f Targets**    |    //    |  🟦   |    //    |          //          |
| **G2M Checkpoint** |    //    |  🟦   |    //    |          //          |
| **TNFa via NFkB**  | light 🟦 |  🟥   | light 🟦 |          🟥          |
| **Myogenesis**     |    //    |  🟥   |    //    |          🟥          |

#### **Analysis: The H9 Specificity**
The hand-verified H9 matrices reveal three critical biological nuances:
1. **Potency Gradation**: The presence of "a bit of red" in the Kromastat (`K_6nM`) column for **MYC Targets V1** (despite the overall blue shift) provides granular proof of Kromastat's weaker repressive force compared to Romidepsin's absolute suppression (`🟦`).
2. **The Immune Toggle**: The most significant divergence is seen in **Allograft Rejection**. While Romidepsin forces a total "Flip" to deep blue (immunosuppression), Kromastat actually increases the baseline activation ("more light red"). This suggests a potential immuno-stimulatory or stress-inducing effect unique to Kromastat in healthy-like cells.
3. **Conserved Differentiation**: Despite the potency gap and immune differences, both drugs successfully flip the **Myogenesis** signature from light blue to red, confirming that this differentiation pathway is a robust, conserved response to HDAC inhibition in the H9 background.

---

### 4.3 Local Response: SUPM2 (Cancer-like)

| SUPM2: Romidepsin (6nM)                                                                            | SUPM2: Kromastat (6nM)                                                                                        |
| :------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------ |
| ![SUPM2 Romi\|825](../../results/human/figures/04_04_heatmap_pathway_SUPM2_Romi_6nM_vs_DMSO_Romi.png) | ![SUPM2 Kroma\|825](../../results/human/figures/04_04_heatmap_pathway_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **SUPM2 Local Matrix (Cancer-like)**
##### R vs D

| Pathway            |   D_R    | R_6nM |   D_K    |       K_6nM        |
| :----------------- | :------: | :---: | :------: | :----------------: |
| **E2F Targets**    | light 🟥 |  🟦   | light 🟥 | generally light 🟦 |
| **MYC Targets V1** |    //    |  🟦   |    //    |         //         |
| **G2M Checkpoint** |    //    |  🟦   |    //    |         //         |
| **TNFa via NFkB**  | light 🟦 |  🟥   | light 🟦 |      light 🟥      |
| **MYC Targets V2** | light 🟥 |  🟦   | light 🟥 |      light 🟦      |
| **Hypoxia**        | light 🟦 |  🟥   | light 🟦 |      light 🟥      |

##### K vs D

| Pathway            |   D_R    | R_6nM |   D_K    |        K_6nM         |
| :----------------- | :------: | :---: | :------: | :------------------: |
| **E2F Targets**    | light 🟥 |  🟦   | light 🟥 | generally lighter 🟦 |
| **MYC Targets V1** |    //    |  🟦   |    //    |          //          |
| **G2M Checkpoint** |    //    |  🟦   |    //    |          //          |
| **MYC Targets V2** |    //    |  🟦   |    //    |          //          |
| **TNFa via NFkB**  | light 🟦 |  🟥   | light 🟦 |          🟥          |
| **Complement**     |    //    |  🟥   |    //    |          🟥          |

#### **Analysis: The Cancer Response**
1. **Conserved Stress Signature**: In the cancer background, both drugs converge on the activation of **TNFa Signaling** and the **Complement** cascade. This suggests that HDAC inhibition in cancer cells triggers a universal "inflammatory stress" state regardless of the specific inhibitor used.
2. **Potency Across the Axis**: The verified matrices provide the most granular evidence of the "Potency Gap." While Romidepsin achieves absolute suppression (**Deep Blue**) of E2F and MYC targets, Kromastat only reaches a **Light Blue** state. This suggests that at 6nM, Kromastat is a "Sub-optimal" inhibitor of the cancer cell cycle compared to Romidepsin.
3. **Hypoxia Priming**: Romidepsin forces a significantly stronger **Hypoxia** signature than Kromastat. This "metabolic suffocating" effect may explain why Romidepsin is more effective at inducing cell death in cancer lineages.

---

### 4.4 Direct Comparison: Romidepsin vs. Kromastat
This heatmap visualizes the functional difference between the two drugs at the same 6nM concentration.

![Romi vs Kroma Heatmap|825](../../results/human/figures/04_04_heatmap_pathway_Romi_6nM_vs_Kromastat_6nM.png)

| Pathway                       |             D_R             | R_6nM |                      D_K                       |                       K_6nM                        |
| :---------------------------- | :-------------------------: | :---: | :--------------------------------------------: | :------------------------------------------------: |
| **MYC Targets V1**            |          light 🟥           |  🟦   |                    light 🟥                    | generally very pale and balanced between 🟥 and 🟦 |
| **G2M Checkpoint**            |             //              |  🟦   |                       //                       |                         //                         |
| **E2F Targets**               |             //              |  🟦   |                       //                       |                         //                         |
| **MYC Targets V2**            |             //              |  🟦   |                       //                       |                         //                         |
| **Interferon Alpha Response** | light 🟥 with patches of 🟦 |  🟦   | mostly light 🟥 with small patches of light 🟦 |                  slightly more 🟥                  |
| **Allograft Rej.**            |             //              |  🟦   |                       //                       |                         //                         |

#### **Analysis: The Direct Divergence**
1. **The Metabolic Gap**: The contrast between Romidepsin’s absolute suppression (`🟦`) and Kromastat’s "pale/balanced" signature for **MYC and E2F targets** provides the ultimate explanation for the potency differences seen throughout the pipeline. Romidepsin acts as a total metabolic silencer at 6nM, whereas Kromastat acts as a partial modulator.
2. **Immune Signaling Preservation**: The most functionally significant divergence is seen in the **Interferon Alpha Response**. Romidepsin suppresses this vital immune-signaling pathway, while Kromastat maintains or slightly increases its baseline activation.
3. **Clinical Implication**: This suggests that Kromastat, despite its lower "growth-crushing" potency, may be a superior candidate for combination therapies with **Immunotherapy**, as it lacks the broad immunosuppressive footprint (Allograft Rejection/Interferon suppression) of Romidepsin.

---

### 4.5 Summary of Functional Pathway Findings

The functional landscape across H9 and SUPM2 reveals that the difference between Romidepsin and Kromastat is both **quantitative** (potency) and **qualitative** (biological fate).

#### **1. The Shared Core: "The Big Flip"**
Both drugs successfully achieve the primary goal of HDAC inhibition: a massive **transcriptomic inversion**. They flip the basal "Light Red" proliferative signatures (MYC, E2F, G2M) into a suppressed state. This confirms that the core anticancer mechanism is conserved across both molecules.

#### **2. The Potency Gradient (Intensity)**
The matrices visually prove that Romidepsin is a **"Metabolic Crusher,"** inducing a "Deep Blue" wall of suppression. Kromastat, by contrast, acts as a **"Metabolic Modulator,"** resulting in a "Pale/Light Blue" signature. It hits the same targets but with significantly less repressive force at the 6nM concentration.

#### **3. The Functional Divergence**
*   **Immune Sparing**: Romidepsin broadly suppresses immune-related pathways (**Interferon Alpha, Allograft Rejection**). Kromastat uniquely **preserves or activates** these pathways.
*   **Stress Fates**: In cancer cells (SUPM2), Romidepsin pushes the cells toward a **Hypoxia** stress state, while Kromastat triggers a **Complement** immune-cascade.
*   **Healthy-Line Safety**: In H9 cells, Kromastat's weaker repression suggests it may be less "toxic" to healthy cell metabolism than Romidepsin.

#### **Final Verdict**
**Romidepsin** is the superior "Growth Crusher" (Metabolic Overdrive), but **Kromastat** is the superior "Immune Primer." Kromastat's ability to inhibit growth without shutting down interferon signaling makes it a high-potential candidate for combination with **Immunotherapies**.

---

## 5. Global Transcriptomic Variance (Top Variable Genes)

This analysis focuses on the top 50 genes with the highest variance across the entire human dataset, providing a high-resolution view of the biological drivers behind the PCA clusters. Cluster assignments are cross-verified against the DGE tables (`baseMean` expression).

![Top Variable Genes Heatmap|825](../../results/human/figures/04_05_heatmap_top_variable.png)

### 5.1 Biological Drivers of PC1 (Cell Line Identity)
The most dominant feature is the massive vertical divide between **H9 (Green)** and **SUPM2 (Purple)**, confirming that the primary source of variance is biological lineage rather than drug treatment. The heatmap resolves into two distinct identity clusters:

**Top Cluster → SUPM2 (Cancer) Identity Genes** *(Blue in H9, Red in SUPM2 controls)*
These genes are highly expressed in the cancer background and silenced in the healthy-like baseline:

| Gene | Role | SUPM2 baseMean | Drug Effect (SUPM2+Romi log2FC) |
| :--- | :--- | ---: | ---: |
| **ALK** | Anaplastic Lymphoma Kinase — key lymphoma oncogene | 10,932 | −1.58 (padj ≈ 0) |
| **SLC26A4** | Ion transporter, metabolic reprogramming | 10,287 | −2.90 (padj ≈ 0) |
| **TNC** | Extracellular matrix / tumour microenvironment | 11,380 | −2.70 (padj ≈ 0) |
| **ABCC9** | ATP-sensitive K⁺ channel — tumour stress resistance | 9,896 | −5.07 (padj ≈ 0) |
| **MYBPC2** | Cytoskeletal remodelling | 10,925 | −3.81 (padj ≈ 0) |
| **SERPINA1** | Protease inhibitor — immune evasion | 63,712 | not significantly changed |

> **Oncogenic Validation**: The high expression of **ALK** (Anaplastic Lymphoma Kinase) in the SUPM2 baseline provides a critical internal validation of its lymphoma-derived, cancer-like phenotype.

**Bottom Cluster → H9 (Healthy-like) Identity Genes** *(Red in H9, Blue in SUPM2 controls)*
These genes define the healthy-cell homeostatic programme and are almost silent in SUPM2:

| Gene | Role | SUPM2 baseMean (near-zero) |
| :--- | :--- | ---: |
| **JUP** | Junction plakoglobin — epithelial adhesion & signalling | 11 |
| **COTL1** | Actin cytoskeleton organisation | 40 |
| **GLUL** | Glutamine synthetase — metabolic homeostasis | 10 |
| **L1CAM** | Neural cell adhesion molecule | 12 |
| **TCF7** | Wnt signalling transcription factor | — |

### 5.2 Drug Effect: The "Potency Gap" at Single-Gene Resolution
The variable gene heatmap provides the most granular confirmation of the potency gap identified throughout the pipeline. Focusing on the top SUPM2-identity markers:

- **Romidepsin (R_6nM)**: Forces `ABCC9` (log2FC = **−5.07**) and `MYBPC2` (log2FC = **−3.81**) into **Deep Blue (absolute suppression)**. The SUPM2 cancer signature is comprehensively dismantled.
- **Kromastat (K_6nM)**: Achieves only a **partial reduction** of the same markers (log2FC = −2.46 and −2.09 respectively, SUPM2-specific). The cancer signature fades but is not eliminated at 6nM.

This single-gene level view confirms: Romidepsin is a **"Cancer Signature Crusher"** while Kromastat remains a **"Partial Modulator"** at the same dose.

### 5.3 Technical Quality
- **Replicate Uniformity**: The vertical consistency within each 3-replicate block is near-perfect across all 24 samples, confirming high library prep quality and that all observed patterns are biologically robust rather than technical artefacts.

### 5.4 Deep-Dive: Notable Gene Clusters

#### Mid-Cluster (SERPINB1, TIMP1, CA2, IRF4, G0S2, ADAM19, SERPINB4) — SUPM2 Stress Activation by Romidepsin

This cluster occupies a visually striking zone in the heatmap: **deep red specifically in SUPM2+Romi**, while remaining pale in all H9 columns and under Kromastat. The baseMean values confirm these are near-silent in H9 but massively expressed in SUPM2, and are **further induced by Romidepsin**:

| Gene | H9 baseMean | SUPM2 baseMean | SUPM2+Romi log2FC | SUPM2+Romi padj | Biological Role |
| :--- | ---: | ---: | ---: | :--- | :--- |
| **SERPINB1** | 87 | 38,342 | +5.10 | ≈ 0 | Neutrophil elastase inhibitor — inflammation / immune evasion |
| **TIMP1** | 324 | 110,218 | +3.86 | ≈ 0 | Matrix metalloproteinase inhibitor — ECM stress response |
| **CA2** | 35 | 26,634 | +2.43 | ≈ 0 | Carbonic anhydrase — pH regulation under metabolic stress |
| **IRF4** | 302 | 44,807 | +1.10 | ≈ 0 | Lymphoma survival factor — constitutively elevated in SUPM2 |
| **G0S2** | 8 | 11,988 | −0.22 | ≈ 0 | G0/G1 switch gene — mild suppression suggesting cell cycle arrest |
| **ADAM19** | 7 | 13,634 | not sig | — | SUPM2 identity marker, not drug-responsive |
| **SERPINB4** | 6 | 7,441 | not sig | — | SUPM2 identity marker, not drug-responsive |

> **Interpretation**: The "deep red" in this cluster for Romidepsin/SUPM2 is **not baseline expression** — it is a **drug-induced stress activation** unique to the cancer background. The massive induction of TIMP1 (+3.86) and SERPINB1 (+5.10) by Romidepsin in SUPM2 suggests the cancer cells are mounting a desperate **matrix remodelling and protease-inhibitor stress response** as Romidepsin dismantles their proliferative programme. IRF4's continued elevation is a known survival mechanism in lymphoma — a potential Romidepsin resistance pathway.

---

#### Bottom Cluster (SPOCK2, RCSD1, GAS7, PDE3B, RASSF2, GPR174, TRBV13, MARCKS, IGFBP2, ACTG2) — H9-Specific Kromastat Suppression

This cluster is the mirror image: **strongly expressed in H9 (Green), near-zero in SUPM2 (Purple)**. All genes are significantly suppressed by Kromastat in H9, but nearly absent in SUPM2 (so they cannot be suppressed there):

| Gene | H9 baseMean | SUPM2 baseMean | H9+Kroma log2FC | H9+Kroma padj | Biological Role |
| :--- | ---: | ---: | ---: | :--- | :--- |
| **RASSF2** | 20,244 | 35 | −1.11 | ≈ 0 | Ras-association tumour suppressor |
| **SPOCK2** | 20,417 | 28 | −0.75 | ≈ 0 | Proteoglycan — ECM organisation |
| **ACTG2** | 28,504 | 19 | −0.34 | ≈ 0 | Actin gamma 2 — cytoskeletal integrity |
| **IGFBP2** | 17,638 | 11 | −0.97 | ≈ 0 | IGF binding protein — growth/metabolic signalling |
| **RCSD1** | 14,079 | 8 | −1.46 | ≈ 0 | Actin-associated — cytoskeletal remodelling |
| **TRBV13** | 8,559 | 4 | −1.10 | ≈ 0 | T-cell receptor beta — immune cell identity marker |
| **GPR174** | 8,345 | 6 | −0.95 | ≈ 0 | G-protein coupled receptor — immune cell homing |
| **MARCKS** | 8,827 | 10 | −0.66 | ≈ 0 | PKC substrate — cytoskeletal/membrane signalling |
| **PDE3B** | 8,542 | 5 | −0.60 | ≈ 0 | Phosphodiesterase — cAMP/lipid metabolism |
| **GAS7** | 8,567 | 5 | −0.50 | ≈ 0 | Growth-arrest specific — neuronal/cytoskeletal function |

> **Interpretation**: This cluster constitutes the **healthy-cell homeostatic programme** of H9 — a set of genes governing cytoskeletal integrity, immune identity (TRBV13, GPR174), and growth regulation. The consistent suppression of all these genes by **Kromastat in H9** reveals that at 6nM, Kromastat is actively **dismantling the healthy-cell maintenance programme**, raising a meaningful safety concern relative to Romidepsin, which shows no significant effect on these H9-specific genes. The near-zero baseMean in SUPM2 means this cluster is effectively absent from the cancer background — these genes are simply not part of the cancer cell's identity.

> **Note**: `LOC101927690` visible in the heatmap could not be resolved to a gene symbol in the DGE tables, consistent with its status as an unannotated or low-confidence locus.

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
- **Kroma-Specific Activation in H9**: Kromastat uniquely activates **Apoptosis** and **P53 Pathway** in H9 (both present in its top 10 activated panel; absent from H9+Romi). This indicates Kromastat induces a **controlled, p53-driven apoptotic programme** in the healthy lineage — a qualitatively distinct fate from Romidepsin's stress-dominated response.

---

### 6.3 Local Response: Cancer-like (SUPM2)
The cancer-specific response highlights the "Growth Crushing" effectiveness of each molecule.

| SUPM2: Romidepsin (6nM)                                                                             | SUPM2: Kromastat (6nM)                                                                                         |
| :-------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------- |
| ![SUPM2 Romi NES\|600](../../results/human/figures/04_06_hallmark_nes_SUPM2_Romi_6nM_vs_DMSO_Romi.png) | ![SUPM2 Kroma NES\|600](../../results/human/figures/04_06_hallmark_nes_SUPM2_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- **Leading Suppression**: In both SUPM2 plots, **E2f Targets** is the single most suppressed hallmark (NES ≈ −3.0), overtaking **Myc Targets V1** — reflecting maximal disruption of the G1/S transition in the cancer background.
- **Romi-Specific Activation in SUPM2**: Romidepsin uniquely activates **Apoptosis**, **P53 Pathway**, and **KRAS Signaling Up** in SUPM2. The co-activation of apoptotic and oncogenic stress pathways is a signature of **terminal cancer cell stress** — cells simultaneously dying and fighting back.
- **Kroma-Specific Activation in SUPM2**: Kromastat instead activates **Complement**, **UV Response Down**, **IL2-STAT5 Signaling**, and **Inflammatory Response** — a broader immune and inflammatory activation profile entirely absent from Romidepsin's SUPM2 response.

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

We analyzed the overlap of significant DEGs (defined by `padj < 0.05` and `|LFC| > 1`) across the four primary experimental conditions. This matrix visualizes the "potency" and "uniqueness" of the drug-induced transcriptional changes.

![UpSet Consistency Matrix|717](../../results/human/figures/04_09_upset_consistency.png)

| Observation Category | Key Finding | Biological Interpretation |
|:---|:---|:---|
| **Magnitude of Response** | SUPM2 Romi has the largest set (~7500+ DEGs) | Romidepsin triggers a massively broader transcriptional response in SUPM2, confirming its heightened sensitivity compared to other conditions. |
| **Comparative Potency** | H9 Romi is the second largest set (~5500+ DEGs) | Romidepsin consistently acts as a highly potent perturbagen across different cellular backgrounds. |
| **Drug Selectivity** | Kromastat sets are notably smaller (~3000–4000 DEGs) | Consistent with its targeted profile, Kromastat acts as a more selective perturbagen, inducing a focused transcriptome-wide footprint. |
| **Cell-Line Specificity** | Largest unique intersection = SUPM2 Romi (2,641) | A significant core of genes is uniquely regulated by Romidepsin in SUPM2, likely representing cell-line-specific sensitivity pathways. |
| **Mechanism Consistency** | H9 Krom ∩ SUPM2 Krom consensus (2,364 DEGs) | Despite smaller absolute set sizes, Kromastat displays high cross-cell-line concordance, indicating a very stable mechanistic action. |
| **Conserved HDAC Core** | 4-way all-condition overlap (45 genes) | These 45 genes are significant across ALL conditions and cell lines, representing the "gold-standard" conserved targets of HDAC inhibition. |

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
