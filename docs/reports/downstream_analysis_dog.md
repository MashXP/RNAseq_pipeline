# Downstream Transcriptomic Analysis (Dog)

## 1. Principal Component Analysis (PCA)

The PCA plot visualizes the overall transcriptomic variance across the dog samples.

![Combined PCA Plot|600](../../results/dog/figures/04_01_pca_combined.png)

### 1.1 Interpretation:
- **PC1 (96% Variance)**: The dominant source of variance is **Cell Line Identity**. CNK89 (NK-cell cancer) and UL1 (T-ALL cancer) lineages are sharply bifurcated along the primary axis, indicating significant baseline transcriptomic differences between these dog lymphoid models.
- **PC2 (3% Variance)**: This axis captures the **Drug Treatment Effect**.
    - Both dog cell lines respond to treatment with a consistent vertical shift along PC2.
    - **Romidepsin (6nM)** appears to induce a more pronounced shift than **Kromastat (6nM)**, reflecting the high potency of both inhibitors in the canine context.

### 1.2 Quality Control:
- **Replicate Fidelity**: Replicates for all conditions cluster with high precision, demonstrating excellent technical consistency in the dog dataset.
- **Control Consistency**: The **DMSO_Romi** and **DMSO_Kromastat** controls for each cell line are perfectly overlaid for each group (indicating excellent vehicle stability), though they are not located at the origin (0,0). In the combined plot, controls are positioned at approximately PC2 = -9 (UL1) and PC2 = -3 (CNK89).

### 1.3 Sub-PCA (Cell Line Specific):
Individual PCA analysis for each cell line confirms that when the variance of the other lineage is removed, the drug effect becomes the primary driver of variance (>90%).

| CNK89 (NK-cell Cancer) | UL1 (T-ALL Cancer) |
| :---: | :---: |
| ![CNK89 PCA](../../results/dog/figures/04_01_pca_CNK89.png) | ![UL1 PCA](../../results/dog/figures/04_01_pca_UL1.png) |

## 2. Differential Gene Expression (Volcano Plots)

### 2.1 Global Response Dashboard
The combined dashboard provides a bird's-eye view of all primary contrasts used to validate the drug response in dog lymphoid models.

[Volcano Plot Overview|1000](../../results/dog/figures/04_02_volcano_combined.png)

#### 2.1.1 Control Group Validation
- **Contrast**: `DMSO_Kromastat_vs_DMSO_Romi`
- **Finding**: **0 Significant DEGs** (padj < 0.05).

![DMSO Control Validation|500](../../results/dog/figures/04_02_volcano_DMSO_Kromastat_vs_DMSO_Romi.png)

- **Conclusion**: The "blank" plot confirms that the two control groups are identical. This baseline stability ensures that all observed differential expression is directly attributable to the drug treatments.

#### 2.1.2 Global Magnitude of Response
> [!WARNING]
> **Statistical Caveat**: Global DEG counts are derived from a species-level analysis that pools CNK89 and UL1. Because 96% of the variance is driven by biological lineage (Section 1), pooling these vastly different populations violates basic DGE assumptions and introduces significant confounding. These global metrics are used only for high-level potency screening; the local (cell-line specific) results in Section 2.3 are the definitive findings.

| Treatment            | Total DEGs (\|log2FoldChange\| > 2) | Ratio (Up:Down)               |
| :------------------- | :-------------------------- | :---------------------------- |
| **Romidepsin (6nM)** | 3,494 [REVERIFY]            | 3.4 : 1 (2,705 up / 789 down) |
| **Kromastat (6nM)**  | 1,327 [REVERIFY]            | 22.3 : 1 (1,270 up / 57 down) |

Both inhibitors show a significant bias toward gene activation in dog cells. However, Kromastat's response is exceptionally skewed, with almost no gene suppression observed at this concentration.

### 2.2 Global Head-to-Head: Romidepsin vs. Kromastat
This contrast directly compares the two drugs at the same 6nM dose.

![Romi vs Kroma Volcano|500](../../results/dog/figures/04_02_volcano_Romi_6nM_vs_Kromastat_6nM.png)

- **Genes higher in Romidepsin (Red)**: Includes 394 targets where Romidepsin has a more potent "overdrive" effect.
- **Genes higher in Kromastat (Blue)**: Includes 290 targets such as `RXRA`, `RASA3`, and `TMEM71`. These represent genes that are significantly better preserved or more strongly activated by Kromastat compared to Romidepsin.
- **Conclusion**: With **684 significant differences** between the two drugs at the same molarity, the dog data confirms that these inhibitors are not functionally redundant.

### 2.3 Local (Cell-Line Specific) Response
Analysis of the lineages individually highlights significant differences in transcriptomic reactivity.

> [!NOTE]
> Local reactivity counts are derived from individual cell-line DGE tables (e.g., `02_dge_UL1_*.csv`) located in: [results/dog/tables/](../../results/dog/tables/).

| Cell Line             | Romidepsin (6nM) | Kromastat (6nM) |
| :-------------------- | :--------------- | :-------------- |
| **UL1 (T-ALL Cancer)** | 5,045 DEGs [REVERIFY] | 2,659 DEGs [REVERIFY] |
| **CNK89 (NK Cancer)** | 2,245 DEGs [REVERIFY] | 589 DEGs [REVERIFY]   |

| Cell Line             | Romidepsin (6nM)                                                                           | Kromastat (6nM)                                                                                       |
| :-------------------- | :----------------------------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------- |
| **UL1 (T-ALL Cancer)** | ![UL1 Romi\|500](../../results/dog/figures/04_02_volcano_UL1_Romi_6nM_vs_DMSO_Romi.png)       | ![UL1 Kroma\|500](../../results/dog/figures/04_02_volcano_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png)       |
| **CNK89 (NK Cancer)** | ![CNK89 Romi\|500](../../results/dog/figures/04_02_volcano_CNK89_Romi_6nM_vs_DMSO_Romi.png) | ![CNK89 Kroma\|500](../../results/dog/figures/04_02_volcano_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### Key Observations:
1. **Reactivity**: The T-cell leukemia line (**UL1**) shows significantly higher transcriptomic reactivity than the NK-cell cancer line (**CNK89**) across both treatments. Specifically, UL1 responds with **5,045 DEGs** to Romidepsin compared to **2,245** in CNK89, and maintains higher reactivity to Kromastat (**2,659** vs **589**).
2. **Conserved Targets**: Highly significant targets such as `RXRA` and `TACC3` are consistently modulated across both dog cell lines, indicating a conserved drug mechanism that transcends the cell line identity.
## 3. Gene Overlap Analysis (Venn Diagrams)

The Venn diagrams illustrate the degree of transcriptomic "shared identity" between Romidepsin and Kromastat at the 6nM concentration.

| All DEGs                                                       | Up-regulated                                                 | Down-regulated                                                   |
| :------------------------------------------------------------- | :----------------------------------------------------------- | :--------------------------------------------------------------- |
| ![Venn All\|400](../../results/dog/figures/04_03_venn_all.png) | ![Venn Up\|400](../../results/dog/figures/04_03_venn_up.png) | ![Venn Down\|400](../../results/dog/figures/04_03_venn_down.png) |

### 3.1 Statistical Rigor Summary
To ensure these overlaps are biologically meaningful and not due to chance, we calculated significance using a hypergeometric test.

> [!NOTE]
> Venn set sizes below are from the `04_03_venn_rigor_stats.csv` snapshot. Current DGE tables yield **3,494 [REVERIFY]** (Romi) and **1,327 [REVERIFY]** (Kroma) at the same thresholds.

| Gene Set           | Overlap P-value | Jaccard Index | Representation Factor |
| :----------------- | :-------------- | :------------ | :-------------------- |
| **All DEGs**       | < 2.2e-16       | 0.311         | 4.0x                  |
| **Up-regulated**   | < 2.2e-16       | 0.382         | 5.2x                  |
| **Down-regulated** | 6.87e-49        | 0.056         | **16.4x**             |

> [!TIP]
> **Full Statistical Report**: The complete overlap metrics (including Romi/Kroma set sizes) are exported to: [04_03_venn_rigor_stats.csv](../../results/dog/tables/04_03_venn_rigor_stats.csv)

**Conclusion**: The overlap is highly significant across all sets. The **16.4x Representation Factor** for down-regulated genes is statistically striking, though it is derived from a minimal absolute pool of genes (57), indicating a nascent rather than robust silencing program for Kromastat at this dose.

### 3.2 Top Shared Targets
These genes represent the "Core Response" triggered by both inhibitors:

| Direction          | Top Shared Genes (padj ~ 0)                                            | Biological Context                              |
| :----------------- | :--------------------------------------------------------------------- | :---------------------------------------------- |
| **Up-regulated**   | `TMEM241`, `DNAJC12`, `HOOK2`, `H2AC9`, `ELAPOR1`                      | Epigenetics, Proteostasis, Lysosomal Function   |
| **Down-regulated** | `GAS7`, `CTSS`, `CRACR2A`, `FERMT3`, `IRAK1`, `E2F8`, `ITGB2`, `CRLF2` | **Cell Cycle**, Immune Signaling, Cell Adhesion |

The shared suppression of **E2F8** and **GAS7** appears to confirm a conserved anti-proliferative mechanism; however, this "global" signal is almost entirely driven by the high-reactivity **UL1** line. In reality, Kromastat only suppresses **16 genes total** in the CNK89 NK cancer line (Section 2.3), suggesting that the anti-proliferative program is largely resisted by the malignant background at 6nM.

### 3.3 The "Subset" Hypothesis
The rigor metrics support a clear functional hierarchy:
- **Shared Response (31%) [REVERIFY]**: 1,144 genes [REVERIFY] forming the core HDAC inhibition signature.
- **Kromastat Uniqueness (5%) [REVERIFY]**: Only 183 genes [REVERIFY].
- **Romidepsin Uniqueness (64%) [REVERIFY]**: 2,350 genes [REVERIFY].

**Finding**: Kromastat acts primarily as a **transcriptomic subset** of Romidepsin in dog cell lines. Approximately **86% of Kromastat's total response is contained within the Romidepsin signature**. While this confirms that Kromastat rarely triggers novel pathways, its **differentiation** arises from:
1. **The 183 Unique Targets**: These Kromastat-specific genes (including immune receptors like `CD86`) may drive the specialized cytokine signatures observed later.
2. **Lineage-Specific Immune Preservation**: By *not* suppressing the ~2,000 genes that Romidepsin does (including key immune regulators), Kromastat preserves the cellular infrastructure in the T-cell lineage (**UL1**) necessary for its role as a **Superior Immune Activator**. Its value lies in what it **selectively spares** in the T-cell background compared to the broader suppression of Romidepsin.

## 4. Functional Pathway Analysis (GSEA Leading-Edge)

Pathway analysis connects gene lists to biological phenomena. Both drugs induce a massive suppression of the cell cycle machinery in the canine context.

### 4.1 Global Functional Response (Signature Matrices)

| Romidepsin vs DMSO                                                                             | Kromastat vs DMSO                                                                                         |
| :--------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------- |
| ![Romi_Heatmap\|825](../../results/dog/figures/04_04_heatmap_pathway_Romi_6nM_vs_DMSO_Romi.png) | ![Kroma_Heatmap\|825](../../results/dog/figures/04_04_heatmap_pathway_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **Romidepsin (Global)**
| Pathway             |   D_R    |        R_6nM         |   D_K    |             K_6nM             |
| :------------------ | :------: | :------------------: | :------: | :---------------------------: |
| **G2M Checkpoint**  | light 🟥 | 🟦 (stronger in UL1) | light 🟥 | light (UL1)🟦/light 🟥(CNK89) |
| **E2F Targets**     |    //    |          //          |    //    |              //               |
| **MYC Targets V1**  |    //    |          //          |    //    |              //               |
| **MYC Targets V2**  |    //    |          //          |    //    |              //               |
| **Mitotic Spindle** |    //    |          //          |    //    |              //               |
| **Allograft Rej.**  |    //    |          //          |    //    |              //               |

#### **Kromastat (Global)**
| Pathway             |   D_R    | R_6nM |   D_K    |             K_6nM             |
| :------------------ | :------: | :---: | :------: | :---------------------------: |
| **E2F Targets**     | light 🟥 |  🟦   | light 🟥 | light (UL1)🟦/light 🟥(CNK89) |
| **G2M Checkpoint**  |    //    |  🟦   |    //    |              //               |
| **MYC Targets V1**  |    //    |  🟦   |    //    |              //               |
| **MYC Targets V2**  |    //    |  🟦   |    //    |              //               |
| **TNFa via NFkB**   | light 🟦 |  🟥   | light 🟦 |      🟥 (strong in UL1)       |
| **Mitotic Spindle** | light 🟥 |  🟦   | light 🟥 | light (UL1)🟦/light 🟥(CNK89) |


#### **Interpreting "The Big Flip"**
The matrices reveal a striking **Reciprocal Signature** between the control state and the drug state, with notable cell-line-specific nuances:
1. **The Proliferative Baseline**: In both control groups (`D_R`, `D_K`), the cells exhibit a mild "Light Red" signature in cell cycle and metabolic pathways (G2M, E2F, MYC). This represents a basal state of active proliferation across both dog lineages.
2. **The Uneven Flip (Romidepsin)**: Romidepsin forces a massive inversion, turning proliferative pathways into **Deep Blue** (Absolute Suppression). However, this suppression is visually much stronger in the T-ALL leukemia line (**UL1**) than in the NK-cell cancer line (**CNK89**), suggesting the NK cancer line retains some proliferative resilience at 6nM.
3. **The Partial Potency & Stress Activation (Kromastat)**: Kromastat attempts the same "flip," but its signature in the `K_6nM` column is **generally lighter** (light blue or even white/neutral), providing visual confirmation of the "Potency Gap." Furthermore, Kromastat uniquely and strongly activates **TNFa Signaling via NFkB** in the UL1 T-ALL background, indicating a distinct stress or inflammatory response not shared by Romidepsin.

---

### 4.2 Local Response: T-ALL Leukemia (UL1)

| UL1: Romidepsin (6nM)                                                                         | UL1: Kromastat (6nM)                                                                                     |
| :------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| ![UL1 Romi\|825](../../results/dog/figures/04_04_heatmap_pathway_UL1_Romi_6nM_vs_DMSO_Romi.png) | ![UL1 Kroma\|825](../../results/dog/figures/04_04_heatmap_pathway_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **UL1 Local Matrix (T-ALL Leukemia)**
##### R vs D

| Pathway             |               D_R                |       R_6nM        |               D_K                |                  K_6nM                  |
| :------------------ | :------------------------------: | :----------------: | :------------------------------: | :-------------------------------------: |
| **E2F Targets**     |             light 🟥             | 🟦 (Strong in UL1) |             light 🟥             | light 🟦 in UL1, light 🟥 in the other. |
| **G2M Checkpoint**  |                //                |         //         |                //                |                   //                    |
| **MYC Targets V1**  |                //                |         //         |                //                |                   //                    |
| **MYC Targets V2**  |                //                |         //         |                //                |                   //                    |
| **Allograft Rej.**  | light 🟥, a bit darker 🟥 in UL1 |         //         | light 🟥, a bit darker 🟥 in UL1 |         Mix of light 🟦 and 🟥          |
| **Mitotic Spindle** |             light 🟥             |         //         |             light 🟥             | light 🟦 in UL1, light 🟥 in the other. |

##### K vs D

| Pathway            |   D_R    | R_6nM |   D_K    |        K_6nM         |
| :----------------- | :------: | :---: | :------: | :------------------: |
| **E2F Targets**    | light 🟥 |  🟦   | light 🟥 |          🟦          |
| **G2M Checkpoint** |    //    |  🟦   |    //    |          🟦          |
| **MYC Targets V1** |    //    |  🟦   |    //    |          🟦          |
| **MYC Targets V2** |    //    |  🟦   |    //    |          🟦          |
| **TNFa via NFkB**  | light 🟦 |  🟥   | light 🟦 |          🟥          |
| **Hypoxia**        | light 🟦 |  🟥   | light 🟦 |          🟥          |

#### **Analysis: UL1 Lineage Response**
Referencing the UL1 local matrices provides direct visual evidence for these biological shifts:
1. **Differential Anti-Proliferative Potency:** As shown in the UL1 matrix for E2F Targets, Romidepsin completely inverts the baseline (`light 🟥`) to absolute suppression (`🟦`). Kromastat achieves a noticeably weaker suppression (`light 🟦`), confirming lower potency at 6nM.
2. **Baseline Immune Pathway Modulation:** The `D_R` and `D_K` columns highlight a `darker 🟥` baseline for Allograft Rejection in UL1. Romidepsin represses this signature. Kromastat, however, results in a "Mix of `light 🟦` and `🟥`" in the `K_6nM` column, demonstrating partial preservation of baseline immune signaling in this T-cell background.
3. **Targeted Stress Response Induction:** In the K vs D matrix, the `K_6nM` column exclusively shows activation (`🟥`) for TNFα via NF-κB and Hypoxia. This red signature is absent in the Romidepsin profile, identifying Kromastat as a distinct stress-modulator in the UL1 leukemia lineage.

---

### 4.3 Local Response: NK-Cell Cancer (CNK89)

| CNK89: Romidepsin (6nM)                                                                            | CNK89: Kromastat (6nM)                                                                                        |
| :------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------ |
| ![CNK89 Romi\|825](../../results/dog/figures/04_04_heatmap_pathway_CNK89_Romi_6nM_vs_DMSO_Romi.png) | ![CNK89 Kroma\|825](../../results/dog/figures/04_04_heatmap_pathway_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **CNK89 Local Matrix (NK-Cell Cancer)**
##### R vs D

| Pathway              |   D_R    | R_6nM |   D_K    |  K_6nM   |
| :------------------- | :------: | :---: | :------: | :------: |
| **G2M Checkpoint**   | light 🟥 |  🟦   | light 🟥 | light 🟦 |
| **E2F Targets**      |    //    |  🟦   |    //    | light 🟦 |
| **MYC Targets V1**   |    //    |  🟦   |    //    | light 🟦 |
| **MYC Targets V2**   |    //    |  🟦   |    //    | light 🟦 |
| **Mitotic Spindle**  |    //    |  🟦   |    //    | light 🟦 |
| **Mtorc1 Signaling** |    //    | 🟦 (Subtle) |    //    |    //    |

##### K vs D

| Pathway              |   D_R    |                   R_6nM                   |   D_K    |                      K_6nM                       |
| :------------------- | :------: | :---------------------------------------: | :------: | :----------------------------------------------: |
| **E2F Targets**      | light 🟥 |           Dark 🟦 mostly in UL1           | light 🟥 | mostly light 🟦 in UL1, with patches of light 🟥 |
| **G2M Checkpoint**   |    //    |                    //                     |    //    |                        //                        |
| **MYC Targets V1**   |    //    |                    //                     |    //    |                        //                        |
| **TNFa via NFkB**    | light 🟦 | light 🟦 in UL1, but DARK 🟥 in the other | light 🟦 |            // (slightly more 🟥 here)            |
| **Mtorc1 Signaling** | light 🟥 |    🟦 in UL1 (except for IDH1 DARK 🟥)    | light 🟥 |                        //                        |
| **MYC Targets V2**   | light 🟥 |                    //                     | light 🟥 |                        //                        |

#### **Analysis: CNK89 Lineage Response**
The CNK89 local matrices illustrate a more resistant malignant phenotype:
1. **Transcriptomic Rigidity:** In the R vs D matrix, Romidepsin achieves a clear shift to suppression (`🟦`) across G2M and E2F targets. However, the K vs D matrix shows Kromastat fails to replicate this; the `K_6nM` column remains predominantly `light 🟦` or retains `light 🟥` patches, indicating resistance to systematic cell cycle silencing in the NK-cell background.
2. **Metabolic Pathway Alterations:** While visually subtle in the heatmap, the CNK89 NK-cancer line shows early signs of mTORC1 Signaling suppression under Romidepsin (confirmed by GSEA in Section 6.3). This suggests that while Romidepsin can modulate this metabolic node, the NK background retains higher resilience compared to the T-cell leukemia lineage.

---

### 4.4 Direct Comparison: Romidepsin vs. Kromastat
This heatmap visualizes the functional difference between the two drugs at the same 6nM concentration.

![Romi vs Kroma Heatmap|825](../../results/dog/figures/04_04_heatmap_pathway_Romi_6nM_vs_Kromastat_6nM.png)

| Pathway                    |   D_R    |                  R_6nM                   |   D_K    |              K_6nM               |
| :------------------------- | :------: | :--------------------------------------: | :------: | :------------------------------: |
| **MYC Targets V1**         | light 🟥 | darker 🟦 on UL1, paler on the other end | light 🟥 |        light 🟦/light 🟥         |
| **G2M Checkpoint**         |    //    |                    //                    |    //    |                //                |
| **E2F Targets**            |    //    |                    //                    |    //    |                //                |
| **Allograft Rejection**    |    //    |                    //                    |    //    | // (but dark 🟥 at CD86 of UL1)  |
| **MYC Targets V1**         |    //    |                    //                    |    //    | // (but dark 🟥 at DUSP2 of UL1) |
| **Interferon Alpha Resp.** |    //    |   // (but HERC6 in CNK89 is darker 🟦)   |    //    |                //                |

#### **Analysis: Direct Drug Divergence**
The direct comparison matrix (Romidepsin vs. Kromastat) highlights key gene-level deviations:
1. **Lineage-Dependent Efficacy:** In the `R_6nM` column for MYC Targets V1, the matrix explicitly notes the suppression is "darker 🟦 on UL1, paler on the other end" (CNK89). This visually confirms that Romidepsin's efficacy is restricted in the NK-cell cancer background.
2. **Immune Activation vs. Suppression:** The `K_6nM` column reveals isolated immune priming, specifically noted as "dark 🟥 at CD86" (Allograft Rejection) and "dark 🟥 at DUSP2" (MYC Targets V2) in the UL1 T-ALL background. Conversely, the `R_6nM` column demonstrates targeted immune suppression, marked by a "darker 🟦" state at HERC6 (Interferon Alpha Response).

---

### 4.5 Summary of Functional Pathway Findings

By tracking the matrix shifts across UL1 and CNK89, the data indicates both quantitative and qualitative differences between the two drugs at 6nM.

#### **1. Conserved Anti-Proliferative Effects**
Both compounds successfully shift basal proliferative signatures (E2F, G2M, MYC targets) from red to blue, confirming a shared core mechanism of cell cycle suppression.

#### **2. Lineage-Dependent Potency Gap**
Romidepsin exhibits a deep blue, uniform suppression, particularly in the UL1 line. Kromastat shows a lighter blue signature overall and fails to systematically silence the cell cycle machinery (`light 🟥` retention) of the CNK89 NK cancer line at this concentration.

#### **3. Divergent Immune Modulation**
Romidepsin broadly represses immune-related signaling (shifting Allograft Rejection and Interferon Alpha responses to blue). Kromastat demonstrates a divergent functional fate by selectively preserving or activating specific immune nodes in the T-cell leukemia lineage, resulting in distinct dark red hotspots like CD86. 

#### Final Verdict
Romidepsin operates as a direct cytotoxic agent with a broad immunosuppressive footprint. Kromastat functions as an immune-modulating agent, executing a partial anti-proliferative program while preserving the specific red biological signatures (such as CD86 expression) necessary for potential synergy with immunotherapies, specifically in a T-cell malignancy background.

---


## 5. Global Transcriptomic Variance (Top Variable Genes)

> [!IMPORTANT]
> **Gene Selection Bias**: Because 96% of the dataset's variance is driven purely by Cell Line Identity (Section 1), the "Top 50 Most Variable Genes" are predominantly lineage markers. Analyzing drug mechanisms based on these genes highlights **lineage disruption** rather than overall drug efficacy, which is better evaluated through treatment-specific contrasts.

![Top Variable Genes Heatmap|825](../../results/dog/figures/04_05_heatmap_top_variable.png)

### 5.1 Biological Drivers of PC1 (Cell Line Identity)
The most dominant feature is the vertical divide between **UL1 (Green)** and **CNK89 (Purple)**, confirming that the primary source of variance is biological lineage rather than drug treatment. The heatmap resolves into two distinct identity clusters:

**Top Cluster → UL1 (T-ALL Leukemia) Identity Genes** *(Red in UL1, Blue in CNK89)*
These genes are highly expressed in the T-ALL leukemia background and near-silent in the NK cancer line:

| Gene | Role | UL1 baseMean | CNK89 baseMean |
| :--- | :--- | ---: | ---: |
| **VWA1** | von Willebrand factor A domain — extracellular matrix integrity | 23,508 | 285 |
| **LITAF** | LPS-induced TNF-alpha factor — inflammatory signalling | 21,723 | 42 |
| **ATP9A** | Phospholipid flippase — membrane homeostasis | 15,404 | 11 |
| **CDH17** | Cadherin 17 — cell-cell adhesion | 12,076 | 76 |
| **HORMAD1** | HORMA domain protein — DNA damage checkpoint | 9,968 | 7 |
| **ARHGAP18** | Rho GTPase activating protein — cytoskeletal regulation | 10,307 | 11 |
| **MAF** | Transcription factor — T-cell immune differentiation | 9,576 | 57 |

> **Lineage Validation**: The high constitutive expression of T-cell and structural regulators (**LITAF**, **MAF**, **ARHGAP18**) in UL1 relative to CNK89 confirms the T-lymphocyte lineage of this model — genes governing T-cell signalling and cytoskeletal integrity dominate its baseline transcriptome.

**Bottom Cluster → CNK89 (NK-cell Cancer) Identity Genes** *(Red in CNK89, Blue in UL1)*
These genes are highly expressed in the NK cancer background and near-absent in the T-ALL leukemia line:

| Gene | Role | CNK89 baseMean | UL1 baseMean |
| :--- | :--- | ---: | ---: |
| **LCK** | Lymphocyte-specific kinase — T-cell receptor signalling | 13,675 | 206 |
| **BRINP3** | BMP/Retinoic acid signalling — cancer identity | 14,199 | 52 |
| **TRAT1** | T-cell receptor-associated transmembrane adapter | 18,352 | 54 |
| **PTK2B** | Protein tyrosine kinase — T-cell activation | 13,707 | 65 |
| **CD96** | T-cell surface glycoprotein — immune checkpoint | 6,002 | 8 |
| **CD101** | Immunoglobulin superfamily — T-cell co-inhibitory receptor | 12,278 | 92 |
| **CALR** | Calreticulin — ER stress / cancer immunogenicity | 36,636 | 496 |
| **LTB** | Lymphotoxin-beta — lymphoid tumour microenvironment | 15,846 | 19 |
| **THEMIS** | Thymocyte selection — T-cell cancer identity | 29,625 | 33 |
| **NFKB1** | NF-κB subunit — cancer survival / immune evasion | 10,155 | 45 |
| **RASSF2** | Ras-association tumour suppressor | 26,166 | 95 |
| **CAMK4** | Calcium/calmodulin kinase — T-cell signalling | 19,782 | 58 |

> **Oncogenic Validation**: The top of this cluster is anchored by **LCK**, **TRAT1**, **CD96**, and **CD101** — all core components of T-cell receptor signalling and immune checkpoint biology. This pattern of co-expression confirms CNK89's lymphoid tumour-derived identity. The elevated **NFKB1** and **CALR** further indicate active cancer survival and ER stress programmes constitutively engaged in the NK-cancer line.

### 5.2 Drug Effect: The "Potency Gap" at Single-Gene Resolution
The variable gene heatmap provides granular confirmation of the potency gap identified throughout the pipeline, observable at the level of specific lineage identity markers:

- **Romidepsin (R_6nM)**: Strongly suppresses CNK89 identity genes — **TRAT1** (log2FC = **−1.29**, padj ≈ 0) and **BRINP3** (log2FC = **−1.22**, padj ≈ 0) are driven toward UL1-level expression. In UL1, **MAF** is significantly repressed (log2FC = **−1.44**, padj ≈ 0), indicating a broad, cross-lineage suppressive force.
- **Kromastat (K_6nM)**: Achieves only marginal changes at the same CNK89 identity markers — **TRAT1** (log2FC = −0.19) and **BRINP3** (log2FC = −0.18) show only mild suppression, indicating the NK cancer signature is largely preserved at 6nM. In UL1, **CDH17** is significantly activated (log2FC = **+1.01**, padj ≈ 0), suggesting Kromastat reinforces structural adhesion programmes in the T-ALL leukemia background rather than dismantling them.

This single-gene level view confirms: Romidepsin acts as a broad **"Lineage Disruptor"** — aggressively dismantling the identity programmes of both lineages — while Kromastat exerts **selective and partial modulation**, largely preserving the NK cancer signature at the evaluated dose.


### 5.3 Technical Quality
- **Replicate Uniformity**: The vertical consistency within each 3-replicate block is near-perfect across all samples, confirming high library prep quality and that all observed patterns are biologically robust rather than technical artefacts.

---

## 6. Functional Hallmark Analysis (GSEA NES)

This section summarizes the global biological shifts using the Normalized Enrichment Score (NES) from Gene Set Enrichment Analysis (GSEA), providing a quantitative rank of the pathways most impacted by drug treatment in the canine context.

> [!NOTE]
> **Combined Plot Curation**: The combined NES barplot is curated to the **5 biologically informative contrasts** (UL1 and CNK89 per-drug, plus the direct Romi vs. Kromastat comparison). Three contrasts are intentionally excluded from the combined figure but preserved as individual plots:
> - **`DMSO_Kromastat_vs_DMSO_Romi`**: Vehicle QC baseline — near-empty by design.
> - **`Romi_6nM_vs_DMSO_Romi`** and **`Kromastat_6nM_vs_DMSO_Kromastat`**: Global contrasts pooling both cell lines, producing weaker signals.

![Combined Hallmark NES](../../results/dog/figures/04_06_hallmark_nes_combined.png)

### 6.1 Global Functional Response

| Romidepsin (6nM) vs DMSO | Kromastat (6nM) vs DMSO |
| :--- | :--- |
| ![Romi NES\|600](../../results/dog/figures/04_06_hallmark_nes_Romi_6nM_vs_DMSO_Romi.png) | ![Kroma NES\|600](../../results/dog/figures/04_06_hallmark_nes_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **Observations**:
> [!WARNING]
> **GSEA Statistical Bias**: The "Global" functional response is derived by pooling UL1 and CNK89. Because the baseline profiles are 96% divergent, combining them into a single ranked list for GSEA introduces significant mathematical distortion. These scores are a high-level abstraction and should be interpreted with caution compared to the local (per-cell-line) results.

- **Conserved Suppression**: Both drugs trigger a massive proliferation crash. **G2M Checkpoint** (Romi NES −3.18 / Kroma NES −3.23), **E2F Targets** (Romi −3.06 / Kroma −3.30), and **MYC Targets V1** (Romi −2.75 / Kroma −2.93) are the top suppressed signatures. **TNFα Signaling via NF-κB** leads activated signatures for Romidepsin (NES +1.65).
- **Key Divergence**: Romidepsin suppresses **Allograft Rejection (NES −1.96)** and **mTORC1 Signaling (NES −1.95)** — absent from Kromastat's suppressed list. Kromastat uniquely activates **TNFα Signaling via NF-κB (NES +2.30)**, **Epithelial-Mesenchymal Transition (NES +1.97)**, **IL2-STAT5 Signaling (NES +1.89)**, and **Inflammatory Response (NES +1.72)** — a broad immune activation signature entirely absent from Romidepsin's global profile.

---

### 6.2 Local Response: T-ALL Leukemia (UL1)

| UL1: Romidepsin (6nM) | UL1: Kromastat (6nM) |
| :--- | :--- |
| ![UL1 Romi NES\|500](../../results/dog/figures/04_06_hallmark_nes_UL1_Romi_6nM_vs_DMSO_Romi.png) | ![UL1 Kroma NES\|500](../../results/dog/figures/04_06_hallmark_nes_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- **Shared Core**: Both drugs show the canonical proliferation crash (**E2F Targets**, **G2M Checkpoint**, **MYC Targets V1/V2**) as the dominant suppressed signatures in UL1.
- **Romi-Specific Suppression in UL1**: Romidepsin additionally suppresses **Allograft Rejection (NES −2.28)** and **Inflammatory Response (NES −1.56)** — confirming a broader immunosuppressive character in the T-cell leukemia background.
- **Kroma-Specific Activation in UL1**: Kromastat uniquely activates **TNFα Signaling via NF-κB (NES +2.27)** and **Hypoxia (NES +2.14)**, indicating a distinct oxidative and inflammatory stress programme in the T-cell leukemia lineage.

---

### 6.3 Local Response: NK-Cell Cancer (CNK89)

| CNK89: Romidepsin (6nM) | CNK89: Kromastat (6nM) |
| :--- | :--- |
| ![CNK89 Romi NES\|600](../../results/dog/figures/04_06_hallmark_nes_CNK89_Romi_6nM_vs_DMSO_Romi.png) | ![CNK89 Kroma NES\|600](../../results/dog/figures/04_06_hallmark_nes_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- **Leading Suppression**: **G2M Checkpoint** is the single most suppressed hallmark in both CNK89 plots (Romi NES −2.65, Kroma NES −2.75), indicating strong disruption of mitotic progression in the NK-cell cancer background.
- **Romi-Specific Activation in CNK89**: Romidepsin activates **TNFα Signaling via NF-κB (NES +1.80)** and **Xenobiotic Metabolism (NES +1.65)** — a stress and detoxification signature consistent with cellular crisis.
- **Kroma-Specific Activation in CNK89**: Kromastat drives a dramatically stronger immune activation: **TNFα Signaling via NF-κB (NES +2.21)**, **IL2-STAT5 Signaling (NES +2.06)**, **IL6-JAK-STAT3 Signaling (NES +2.05)**, **Coagulation (NES +2.01)**, and **Inflammatory Response (NES +2.01)**. This broad cytokine activation is entirely absent from Romidepsin's CNK89 response, suggesting Kromastat actively rewires immune signalling in the NK-cell cancer background.

---

### 6.4 Direct Comparison: Romidepsin vs. Kromastat

![Romi vs Kroma NES|600](../../results/dog/figures/04_06_hallmark_nes_Romi_6nM_vs_Kromastat_6nM.png)

- **Romi is more suppressive across all axes**: The panel of pathways **more active in Kromastat than Romidepsin** is dominated by **MYC Targets V1 (NES −2.73)**, **G2M Checkpoint (NES −2.62)**, **Allograft Rejection (NES −2.41)**, **Interferon Alpha Response (NES −2.02)**, **KRAS Signaling Up (NES −1.87)**, and **Interferon Gamma Response (NES −1.77)**. This confirms Kromastat preserves both proliferative and immune signalling infrastructure at a comparatively higher level.
- **Romi activates a narrower stress footprint**: The pathways more active in Romidepsin are anchored by **Heme Metabolism (NES +1.71)** and **Protein Secretion (NES +1.52)** — a narrower, metabolic signature unique to Romidepsin's mechanism.

---

### 6.5 Summary of Hallmark Findings

#### **1. The Shared Core: "The Proliferation Crash"**
Both drugs achieve a massive transcriptomic shutdown of cell cycle programmes. **E2F Targets** and **G2M Checkpoint** consistently emerge as the top suppressed hallmarks (NES < −2.5) across all conditions, reflecting deep disruption of mitotic and S-phase entry.

#### **2. Divergent Immune Fates**
Romidepsin broadly suppresses immune signalling (Allograft Rejection, Interferon Alpha/Gamma), confirming an immunosuppressive footprint. Kromastat executes the opposite — it spares and amplifies immune-related pathways, driving a strong **IL2-STAT5, IL6-JAK-STAT3, and Inflammatory Response** signature in CNK89.

#### **3. Lineage-Specific Stress Response**
In the T-ALL leukemia background, Kromastat triggers a **Hypoxia and TNFα** stress signature absent from Romidepsin, underscoring that Kromastat engages qualitatively different stress pathways in leukemia cells.

**Final Verdict**: **Romidepsin** is the superior direct cytotoxic agent — deeper proliferation suppression, broader immune silencing. **Kromastat** is the superior **immune activator** — executing an anti-proliferative programme while amplifying cytokine signalling, suggesting high potential for synergy with immunotherapies in the canine lymphoma context.

---


## 7. GSEA Hallmark Dotplot Analysis

This section cross-validates the NES barplot findings using the multi-contrast Hallmark dotplot. Each dot encodes three simultaneous dimensions: **position** (Activated / Suppressed), **colour + alpha** (direction × significance), and **size** (leading edge gene count). The globally-aligned Y-axis allows direct row-wise comparison across all five curated contrasts.

> [!NOTE]
> **Combined Plot Curation**: Consistent with the NES barplot (Section 6), the combined dotplot is curated to the **5 biologically informative contrasts**. The DMSO vehicle baseline and global pooled contrasts are excluded from the combined figure but preserved as individual saved plots.

### 7.1 Global Overview

![GSEA Hallmark Combined|1132](../../results/dog/figures/04_07_gsea_dotplot_combined.png)

The combined plot reveals a clear architecture:

- **Top of the Y-axis (Activated zone)**: **TNFα Signaling via NF-κB**, **Inflammatory Response**, **Hypoxia**, **Myogenesis**, and **IL2-STAT5 Signaling** consistently appear as activated (red) dots across multiple contrasts, confirming the universality of the chromatin-stress and immune response to HDAC inhibition in canine cells.
- **Bottom of the Y-axis (Suppressed zone)**: **E2F Targets**, **G2M Checkpoint**, **MYC Targets V1/V2**, and **Mitotic Spindle** anchor the suppressed (blue) end across all drug-treated contrasts, confirming the canonical proliferation crash.

---

### 7.2 Local Response: Cross-Contrast Comparison

#### **Romidepsin (UL1 and CNK89)**

| UL1 — Romi 6nM vs. DMSO | CNK89 — Romi 6nM vs. DMSO |
| :---: | :---: |
| ![UL1 Romi Dotplot\|500](../../results/dog/figures/04_07_gsea_dotplot_UL1_Romi_6nM_vs_DMSO_Romi.png) | ![CNK89 Romi Dotplot\|500](../../results/dog/figures/04_07_gsea_dotplot_CNK89_Romi_6nM_vs_DMSO_Romi.png) |

- **Suppression** is the dominant signal. Large, dark blue dots for **E2F Targets**, **G2M Checkpoint**, **MYC Targets V1/V2**, and **Mitotic Spindle** appear in both UL1 and CNK89, indicating potent and consistent proliferative shutdown.
- **Immune suppression** is specific to Romidepsin: **Allograft Rejection** carries a blue dot in UL1 — entirely absent from Kromastat panels — visually confirming drug-specific immune suppression.
- **UL1+Romidepsin** uniquely shows suppression of **Inflammatory Response**, while **CNK89+Romidepsin** retains activation of **TNFα Signaling via NF-κB** — a pattern of leukemia-specific stress retained even under Romidepsin's aggressive suppression.

---

#### **Kromastat (UL1 and CNK89)**

| UL1 — Kromastat 6nM vs. DMSO | CNK89 — Kromastat 6nM vs. DMSO |
| :---: | :---: |
| ![UL1 Kromastat Dotplot\|500](../../results/dog/figures/04_07_gsea_dotplot_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png) | ![CNK89 Kromastat Dotplot\|500](../../results/dog/figures/04_07_gsea_dotplot_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png) |

- The suppression zone mirrors Romidepsin's proliferation crash, but **Allograft Rejection** and **Interferon Alpha/Gamma Response** are conspicuously **absent** as blue dots — visually confirming the immune-sparing phenotype of Kromastat.
- **TNFα Signaling via NF-κB** carries large red activated dots in both UL1 and CNK89 Kromastat panels. However, **Hypoxia** and **IL2-STAT5 Signaling** are predominantly UL1-specific stress responses and are not significantly activated in the CNK89 NK-cell cancer background at this dose.
- **CNK89+Kromastat** uniquely activates a broad inflammatory cluster including **IL6-JAK-STAT3 Signaling** and **Coagulation** — large-dot activation not present in any other contrast.

---

#### **Direct Romi vs. Kromastat Comparison**

| Romi 6nM vs. Kromastat 6nM |
| :---: |
| ![Romi vs Kromastat Dotplot\|500](../../results/dog/figures/04_07_gsea_dotplot_Romi_6nM_vs_Kromastat_6nM.png) |

- **Suppressed (blue)** — pathways where Romidepsin drives **lower activity than Kromastat**: **MYC Targets V1**, **G2M Checkpoint**, **E2F Targets**, **Allograft Rejection**, **Interferon Alpha Response**, and **Interferon Gamma Response**. This confirms Kromastat preserves both proliferative and immune signalling at a higher relative level.
- **Activated (red)** — pathways where Romidepsin shows **higher activity than Kromastat**: **Heme Metabolism** and **Protein Secretion** — a narrow metabolic stress signature unique to Romidepsin's mechanism.

---

### 7.3 Summary: Dotplot Cross-Validation

The GSEA Hallmark dotplot independently corroborates all major findings from the NES barplot analysis (Section 6):

| Finding | NES Barplot | Dotplot |
| :--- | :---: | :---: |
| Shared proliferation crash (E2F Targets, G2M Checkpoint, MYC Targets V1/V2) | ✅ | ✅ |
| Romidepsin suppresses immune pathways (Allograft Rejection, Interferon Alpha/Gamma) | ✅ | ✅ |
| Kromastat spares immune pathways and activates TNFα / IL2-STAT5 | ✅ | ✅ |
| Kromastat activates Hypoxia in UL1 (leukemia stress) | ✅ | ✅ |
| Kromastat activates IL6-JAK-STAT3 / Coagulation specifically in CNK89 | ✅ | ✅ |
| DMSO vehicles are biologically equivalent | — | ✅ (Verified in full QC plots) |

The dotplot adds a critical new layer: **leading edge size** reveals that the immune-activation pathways in Kromastat (TNFα, IL2-STAT5, IL6-JAK-STAT3) carry **large gene sets**, indicating broad transcriptome-level immune engagement rather than marginal effects. This reinforces the clinical significance of Kromastat's immune-amplifying phenotype in the canine lymphoma context.

---

## 8. Cross-Cell-Line Consistency Proofing

To ensure that the identified biological signatures are robust and not driven by idiosyncratic cell-line artefacts, we performed a consistency analysis comparing the UL1 and CNK89 cell lines. This validation utilized two complementary approaches: (1) intersectional set analysis (UpSet) to evaluate overlap among discrete DEG lists, and (2) global Log2FoldChange (LFC) correlation to assess the directional consensus of the entire transcriptome.

### 8.1 Intersectional Degree of Consensus (UpSet Matrix)

We analyzed the overlap of significant DEGs (defined by `padj < 0.05` and `\|log2FoldChange\| > 2`) across the four primary experimental conditions.

![UpSet Consistency Matrix|717](../../results/dog/figures/04_09_upset_consistency.png)

| Observation Category | Key Finding | Biological Interpretation |
|:---|:---|:---|
| **Magnitude of Response** | UL1 Romi has the largest set (~5,000 DEGs [REVERIFY]) | Romidepsin triggers a broader transcriptional response in the T-ALL background, consistent with its heightened sensitivity profile in the dog dataset. |
| **Comparative Potency** | CNK89 Romi is the second largest set (~2,245 DEGs [REVERIFY]) | Romidepsin consistently acts as a highly potent perturbagen, though the NK cancer line shows reduced sensitivity compared to UL1. |
| **Drug Selectivity** | Kromastat sets are notably smaller (UL1: ~2,659 [REVERIFY] / CNK89: ~589 DEGs [REVERIFY]) | Consistent with its targeted profile, Kromastat acts as a more selective perturbagen with a compressed transcriptomic footprint, especially in the NK-cell cancer background. |
| **Cancer-Line Resistance** | CNK89 Kroma is the smallest set (~589 DEGs [REVERIFY]) | The extremely low DEG count for CNK89+Kromastat indicates significant transcriptomic resistance of the NK-cell cancer line at 6nM — a key finding for dose optimisation. |
| **Mechanism Consistency** | UL1 Krom ∩ CNK89 Krom consensus | Despite smaller absolute set sizes, Kromastat displays cross-cell-line concordance, indicating a stable mechanistic core shared across lineages. |
| **Conserved HDAC Core** | 4-way all-condition overlap | Genes significant across ALL conditions represent the "gold-standard" conserved targets of HDAC inhibition in the canine context. |

### 8.2 Global Transcriptomic Concordance (LFC Correlation)

While intersectional analysis uses binary cutoffs, the LFC correlation evaluates the overall transcriptomic trend across the entire dynamic range, validating the "directionality" of the drug response without threshold bias.

![LFC Correlation|986](../../results/dog/figures/04_10_lfc_correlation.png)

| Analysis Dimension | Resulting Observation | Biological Significance |
|:---|:---|:---|
| **Global Correlation** | Positive R values (~0.42) for both Romi and Kroma | Indicates partial transcriptomic concordance across lineages. However, an R value of 0.42 (R² ≈ 0.17) means 83% of the response is divergent, highlighting significant cell-line specific drug dynamics. |
| **Directional Consensus** | Dominant "Consensus" cluster along the diagonal | The majority of biologically relevant DEGs agree in both direction and significance, though the low correlation coefficient suggests the "Consensus" is limited to core HDAC targets. |
| **Sensitivity Profile** | Romi has a wider LFC range | UL1's response to Romidepsin includes stronger LFC values, confirming the heightened sensitivity of the T-ALL leukemia model compared to CNK89. |
| **Modulation Precision** | Kromastat range is compressed | Kromastat provides a more refined, less disruptive transcriptomic modulation — particularly in the NK cancer background — compared to Romidepsin. |
| **Divergence Clusters** | "Specific" dots scattered off-diagonal | These represent true biological divergence between the cell lines, ensuring the models capture a representative range of canine lymphoid heterogeneity. |

---

## Supplementary: GO ORA Analysis

> [!WARNING]
> **Not suitable for primary presentation.** The GO Over-Representation Analysis (ORA) results (`04_08_ora_dotplot.R`) are flagged as **supplementary reference only** and have not been integrated into the main narrative.

### Reason for Exclusion
The ORA combined and individual dotplots are dominated by **biologically irrelevant GO terms** — primarily housekeeping cellular processes and non-lymphoid developmental pathways — which have limited mechanistic relevance to canine lymphoma or HDAC inhibition. This is a known **gene universe background contamination artefact**: when the background is the full genome, small gene sets with disproportionate overlap can achieve inflated significance scores.

### Recommended Future Action
To correct the ORA signal, restrict the gene universe (`universe` parameter in `enrichGO()`) to **only the genes tested in DESeq2** (i.e., genes that passed independent filtering), rather than the full genome. This is the standard corrective measure and is not prioritized for the current analysis phase.

Individual ORA dotplot files are preserved in `results/dog/figures/04_08_ora_dotplot_*.png` for reference.

---

## 9. Comprehensive Summary & Biological Conclusion

### 9.1 Synthesis of Findings
This downstream transcriptomic analysis systematically compares the effects of two HDAC inhibitors, Romidepsin and Kromastat, at a 6nM dose across T-ALL leukemia (UL1) and NK-cell cancer (CNK89) canine lymphoid cell lines.

1. **The Shared Core (The Proliferation Crash)**: Both drugs successfully execute the primary function of HDAC inhibition — a massive transcriptomic inversion shutting down proliferative pathways (MYC Targets, E2F Targets, G2M Checkpoint) across all cellular backgrounds *(see Sections [[downstream_analysis_dog#4.1 Global Functional Response (Signature Matrices)|4.1]], [[downstream_analysis_dog#6.1 Global Functional Response|6.1]], [[downstream_analysis_dog#7.1 Global Overview|7.1]])*.
2. **The Potency Gap**: Romidepsin acts as a "Metabolic Crusher," driving a deeper, more aggressive transcriptional shutdown — particularly in UL1 (**5,045 DEGs**) — while CNK89 shows reduced sensitivity (**2,245 DEGs**). Kromastat acts as a "Partial Modulator," achieving significant but milder suppression, with CNK89's response (**589 DEGs**) representing a striking case of transcriptomic resistance *(see Sections [[downstream_analysis_dog#1.1 Interpretation:|1.1]], [[downstream_analysis_dog#2.1.2 Global Magnitude of Response|2.1.2]], [[downstream_analysis_dog#5.2 Drug Effect: The "Potency Gap" at Single-Gene Resolution|5.2]])*.
3. **Divergent Immune Modulation**: Romidepsin broadly represses immune-related signalling (Allograft Rejection, Interferon Alpha/Gamma) across both lineages. Conversely, Kromastat actively amplifies immune-receptor pathways — including IL2-STAT5, IL6-JAK-STAT3, and Inflammatory Response — particularly in the CNK89 NK-cell cancer background, a phenotype entirely absent from Romidepsin *(see Sections [[downstream_analysis_dog#4.4 Direct Comparison: Romidepsin vs. Kromastat|4.4]], [[downstream_analysis_dog#6.4 Direct Comparison: Romidepsin vs. Kromastat|6.4]], [[downstream_analysis_dog#7.2 Local Response: Cross-Contrast Comparison|7.2]])*.
4. **The Subset Relationship**: Approximately **86% of Kromastat's total transcriptomic response is contained within the Romidepsin signature**, confirming that Kromastat operates as a transcriptomic subset of Romidepsin in canine cells — rarely engaging pathways not already captured by Romidepsin *(see Section [[downstream_analysis_dog#3.3 The "Subset" Hypothesis|3.3]])*.

### 9.2 Answer to the Biological Question
**Biological Question**: *How do the transcriptomic profiles of Romidepsin and Kromastat differentiate their therapeutic utility as HDAC inhibitors in a canine lymphoma model?*

**Conclusion**: While both drugs target the same core epigenetic and proliferative machinery, their biological fates diverge significantly in both intensity and immune modulation.
- **Romidepsin** is the superior **direct cytotoxic agent**. Its aggressive repression of both the cancer signature *(see Section [[downstream_analysis_dog#5.2 Drug Effect: The "Potency Gap" at Single-Gene Resolution|5.2]])* and immune signalling *(see Sections [[downstream_analysis_dog#4.4 Direct Comparison: Romidepsin vs. Kromastat|4.4]], [[downstream_analysis_dog#6.4 Direct Comparison: Romidepsin vs. Kromastat|6.4]])* makes it highly effective at dismantling tumor metabolism — but its immunosuppressive nature may limit its long-term efficacy by dampening the anti-tumor immune response.
- **Kromastat** is a potential **immune-modulating agent**, but faces a significant **Efficacy Paradox**. While it executes a partial anti-proliferative program *(see Sections [[downstream_analysis_dog#1.1 Interpretation:|1.1]], [[downstream_analysis_dog#4.1 Global Functional Response (Signature Matrices)|4.1]])* without broad immunosuppression, its robust "immune activator" signature (IFN, IL2-STAT5) is currently lineage-dependent (specifically in the T-cell lineage) *(see Sections [[downstream_analysis_dog#4.5 Summary of Functional Pathway Findings|4.5]], [[downstream_analysis_dog#6.5 Summary of Hallmark Findings|6.5]])*. In the malignant NK-cell CNK89 background, Kromastat's response is severely compressed (~589 DEGs), indicating resistance at 6nM *(see Section [[downstream_analysis_dog#2.3 Local (Cell-Line Specific) Response|2.3]])*. Its utility as a cancer therapeutic likely depends on dose-escalation or specific combination strategies to overcome this lineage-specific resistance while leveraging its immune-amplifying properties in compatible cellular backgrounds.
