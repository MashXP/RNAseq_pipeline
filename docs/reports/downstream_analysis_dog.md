# Downstream Transcriptomic Analysis (Dog)

## 1. Principal Component Analysis (PCA)

The PCA plot visualizes the overall transcriptomic variance across the dog samples.

![Combined PCA Plot|600](../../results/dog/figures/04_01_pca_combined.png)

### 1.1 Interpretation:
- **PC1 (96% Variance)**: The dominant source of variance is **Cell Line Identity**. CNK89 (cancer) and UL1 (Healthy-like) lineages are sharply bifurcated along the primary axis, indicating significant baseline transcriptomic differences between these dog lymphoid models.
- **PC2 (3% Variance)**: This axis captures the **Drug Treatment Effect**.
    - Both dog cell lines respond to treatment with a consistent vertical shift along PC2.
    - **Romidepsin (6nM)** appears to induce a more pronounced shift than **Kromastat (6nM)**, reflecting the high potency of both inhibitors in the canine context.

### 1.2 Quality Control:
- **Replicate Fidelity**: Replicates for all conditions cluster with high precision, demonstrating excellent technical consistency in the dog dataset.
- **Control Consistency**: The **DMSO_Romi** and **DMSO_Kromastat** controls for each cell line are perfectly overlaid at the origin for each group. This confirms the stability of the vehicle control across experimental arms.

### 1.3 Sub-PCA (Cell Line Specific):
Individual PCA analysis for each cell line confirms that when the variance of the other lineage is removed, the drug effect becomes the primary driver of variance (>90%).

| CNK89 (cancer) | UL1 (Healthy-like) |
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
> [!NOTE]
> Global DEG counts and ratios are derived from the species-level analysis consolidated in: [04_03_venn_rigor_stats.csv](../../results/dog/tables/04_03_venn_rigor_stats.csv).

| Treatment            | Total DEGs (\|log2FC\| > 1) | Ratio (Up:Down)               |
| :------------------- | :-------------------------- | :---------------------------- |
| **Romidepsin (6nM)** | 3,494                       | 3.4 : 1 (2,705 up / 789 down) |
| **Kromastat (6nM)**  | 1,327                       | 22.3 : 1 (1,270 up / 57 down) |

Both inhibitors show a significant bias toward gene activation in dog cells. However, Kromastat's response is exceptionally skewed, with almost no gene suppression observed at this concentration.

### 2.2 Global Head-to-Head: Romidepsin vs. Kromastat
This contrast directly compares the two drugs at the same 6nM dose.

![Romi vs Kroma Volcano|500](../../results/dog/figures/04_02_volcano_Romi_6nM_vs_Kromastat_6nM.png)

- **Genes higher in Romidepsin (Red)**: Includes 352 targets where Romidepsin has a more potent "overdrive" effect.
- **Genes higher in Kromastat (Blue)**: Includes 291 targets such as `RXRA`, `RASA3`, and `TMEM71`. These represent genes that are significantly better preserved or more strongly activated by Kromastat compared to Romidepsin.
- **Conclusion**: With **643 significant differences** between the two drugs at the same molarity, the dog data confirms that these inhibitors are not functionally redundant.

### 2.3 Local (Cell-Line Specific) Response
Analysis of the lineages individually highlights significant differences in transcriptomic reactivity.

> [!NOTE]
> Local reactivity counts are derived from individual cell-line DGE tables (e.g., `02_dge_UL1_*.csv`) located in: [results/dog/tables/](../../results/dog/tables/).

| Cell Line             | Romidepsin (6nM) | Kromastat (6nM) |
| :-------------------- | :--------------- | :-------------- |
| **UL1 (Healthy-like)** | 5,045 DEGs   | 2,647 DEGs     |
| **CNK89 (Cancer)**    | 2,245 DEGs   | 590 DEGs       |

| Cell Line             | Romidepsin (6nM)                                                                           | Kromastat (6nM)                                                                                       |
| :-------------------- | :----------------------------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------- |
| **UL1 (Healthy-like)** | ![UL1 Romi\|500](../../results/dog/figures/04_02_volcano_UL1_Romi_6nM_vs_DMSO_Romi.png)       | ![UL1 Kroma\|500](../../results/dog/figures/04_02_volcano_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png)       |
| **CNK89 (Cancer)**    | ![CNK89 Romi\|500](../../results/dog/figures/04_02_volcano_CNK89_Romi_6nM_vs_DMSO_Romi.png) | ![CNK89 Kroma\|500](../../results/dog/figures/04_02_volcano_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### Key Observations:
1. **Reactivity**: The healthy-like line (**UL1**) shows significantly higher transcriptomic reactivity than the cancer line (**CNK89**) across both treatments. Specifically, UL1 responds with **5,045 DEGs** to Romidepsin compared to **2,245** in CNK89, and maintains higher reactivity to Kromastat (**2,647** vs **590**).
2. **Conserved Targets**: Highly significant targets such as `RXRA` and `TACC3` are consistently modulated across both dog cell lines, indicating a conserved drug mechanism that transcends the cell line identity.
## 3. Gene Overlap Analysis (Venn Diagrams)

The Venn diagrams illustrate the degree of transcriptomic "shared identity" between Romidepsin and Kromastat at the 6nM concentration.

| All DEGs                                                       | Up-regulated                                                 | Down-regulated                                                   |
| :------------------------------------------------------------- | :----------------------------------------------------------- | :--------------------------------------------------------------- |
| ![Venn All\|400](../../results/dog/figures/04_03_venn_all.png) | ![Venn Up\|400](../../results/dog/figures/04_03_venn_up.png) | ![Venn Down\|400](../../results/dog/figures/04_03_venn_down.png) |

### 3.1 Statistical Rigor Summary
To ensure these overlaps are biologically meaningful and not due to chance, we calculated significance using a hypergeometric test.

> [!NOTE]
> Venn set sizes below are from the `04_03_venn_rigor_stats.csv` snapshot. Current DGE tables yield **3,494** (Romi) and **1,327** (Kroma) at the same thresholds.

| Gene Set           | Overlap P-value | Jaccard Index | Representation Factor |
| :----------------- | :-------------- | :------------ | :-------------------- |
| **All DEGs**       | < 2.2e-16       | 0.311         | 4.0x                  |
| **Up-regulated**   | < 2.2e-16       | 0.382         | 5.2x                  |
| **Down-regulated** | 6.87e-49        | 0.056         | **16.4x**             |

> [!TIP]
> **Full Statistical Report**: The complete overlap metrics (including Romi/Kroma set sizes) are exported to: [04_03_venn_rigor_stats.csv](../../results/dog/tables/04_03_venn_rigor_stats.csv)

**Conclusion**: The overlap is highly significant across all sets. The **16.4x Representation Factor** for down-regulated genes is particularly striking, indicating a highly conserved silencing program despite the small number of genes suppressed by Kromastat.

### 3.2 Top Shared Targets
These genes represent the "Core Response" triggered by both inhibitors:

| Direction          | Top Shared Genes (padj ~ 0)                                            | Biological Context                              |
| :----------------- | :--------------------------------------------------------------------- | :---------------------------------------------- |
| **Up-regulated**   | `TMEM241`, `DNAJC12`, `HOOK2`, `H2AC9`, `ELAPOR1`                      | Epigenetics, Proteostasis, Lysosomal Function   |
| **Down-regulated** | `GAS7`, `CTSS`, `CRACR2A`, `FERMT3`, `IRAK1`, `E2F8`, `ITGB2`, `CRLF2` | **Cell Cycle**, Immune Signaling, Cell Adhesion |

The shared suppression of **E2F8** and **GAS7** confirms a conserved anti-proliferative mechanism, while the suppression of immune regulators like **IRAK1** and **CRLF2** highlights a potential immunosuppressive profile in dog lymphoid cells.

### 3.3 The "Subset" Hypothesis
The rigor metrics support a clear functional hierarchy:
- **Shared Response (31%)**: 1,144 genes forming the core HDAC inhibition signature.
- **Kromastat Uniqueness (5%)**: Only 183 genes.
- **Romidepsin Uniqueness (64%)**: 2,350 genes.

**Finding**: Kromastat acts primarily as a **transcriptomic subset** of Romidepsin in dog cell lines. Approximately **86% of Kromastat's total response is contained within the Romidepsin signature**, demonstrating that while selective, Kromastat rarely triggers pathways that are not already captured by Romidepsin.

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
2. **The Uneven Flip (Romidepsin)**: Romidepsin forces a massive inversion, turning proliferative pathways into **Deep Blue** (Absolute Suppression). However, this suppression is visually much stronger in the healthy-like **UL1** line than in the cancer line **CNK89**, suggesting the cancer line retains some proliferative resilience at 6nM.
3. **The Partial Potency & Stress Activation (Kromastat)**: Kromastat attempts the same "flip," but its signature in the `K_6nM` column is **generally lighter** (light blue or even white/neutral), providing visual confirmation of the "Potency Gap." Furthermore, Kromastat uniquely and strongly activates **TNFa Signaling via NFkB** in the UL1 background, indicating a distinct stress or inflammatory response not shared by Romidepsin.

---

### 4.2 Local Response: Healthy-like (UL1)

| UL1: Romidepsin (6nM)                                                                         | UL1: Kromastat (6nM)                                                                                     |
| :------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------ |
| ![UL1 Romi\|825](../../results/dog/figures/04_04_heatmap_pathway_UL1_Romi_6nM_vs_DMSO_Romi.png) | ![UL1 Kroma\|825](../../results/dog/figures/04_04_heatmap_pathway_UL1_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **UL1 Local Matrix (Healthy-like)**
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
2. **Baseline Immune Pathway Modulation:** The `D_R` and `D_K` columns highlight a `darker 🟥` baseline for Allograft Rejection in UL1. Romidepsin represses this signature. Kromastat, however, results in a "Mix of `light 🟦` and `🟥`" in the `K_6nM` column, demonstrating partial preservation of baseline immune signaling.
3. **Targeted Stress Response Induction:** In the K vs D matrix, the `K_6nM` column exclusively shows activation (`🟥`) for TNFα via NF-κB and Hypoxia. This red signature is absent in the Romidepsin profile, identifying Kromastat as a distinct stress-modulator.

---

### 4.3 Local Response: CNK89 (Cancer-like)

| CNK89: Romidepsin (6nM)                                                                            | CNK89: Kromastat (6nM)                                                                                        |
| :------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------ |
| ![CNK89 Romi\|825](../../results/dog/figures/04_04_heatmap_pathway_CNK89_Romi_6nM_vs_DMSO_Romi.png) | ![CNK89 Kroma\|825](../../results/dog/figures/04_04_heatmap_pathway_CNK89_Kromastat_6nM_vs_DMSO_Kromastat.png) |

#### **CNK89 Local Matrix (Cancer-like)**
##### R vs D

| Pathway              |   D_R    | R_6nM |   D_K    |  K_6nM   |
| :------------------- | :------: | :---: | :------: | :------: |
| **G2M Checkpoint**   | light 🟥 |  🟦   | light 🟥 | light 🟦 |
| **E2F Targets**      |    //    |  🟦   |    //    | light 🟦 |
| **MYC Targets V1**   |    //    |  🟦   |    //    | light 🟦 |
| **MYC Targets V2**   |    //    |  🟦   |    //    | light 🟦 |
| **Mitotic Spindle**  |    //    |  🟦   |    //    | light 🟦 |
| **Mtorc1 Signaling** |    //    |  //   |    //    |    //    |

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
1. **Transcriptomic Rigidity:** In the R vs D matrix, Romidepsin achieves a clear shift to suppression (`🟦`) across G2M and E2F targets. However, the K vs D matrix shows Kromastat fails to replicate this; the `K_6nM` column remains predominantly `light 🟦` or retains `light 🟥` patches, indicating resistance to systematic cell cycle silencing.
2. **Metabolic Pathway Alterations:** The R vs D matrix for CNK89 shows mTORC1 Signaling retains baseline characteristics (`//`) rather than shifting to blue under Romidepsin. This indicates pathway resilience and an inability of Romidepsin to fully suppress this metabolic node in the cancer line.

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
1. **Lineage-Dependent Efficacy:** In the `R_6nM` column for MYC Targets V1, the matrix explicitly notes the suppression is "darker 🟦 on UL1, paler on the other end" (CNK89). This visually confirms that Romidepsin's efficacy is restricted in the cancer background.
2. **Immune Activation vs. Suppression:** The `K_6nM` column reveals isolated immune priming, specifically noted as "dark 🟥 at CD86" (Allograft Rejection) and "dark 🟥 at DUSP2" (MYC Targets V2). Conversely, the `R_6nM` column demonstrates targeted immune suppression, marked by a "darker 🟦" state at HERC6 (Interferon Alpha Response).

---

### 4.5 Summary of Functional Pathway Findings

By tracking the matrix shifts across UL1 and CNK89, the data indicates both quantitative and qualitative differences between the two drugs at 6nM.

#### **1. Conserved Anti-Proliferative Effects**
Both compounds successfully shift basal proliferative signatures (E2F, G2M, MYC targets) from red to blue, confirming a shared core mechanism of cell cycle suppression.

#### **2. Lineage-Dependent Potency Gap**
Romidepsin exhibits a deep blue, uniform suppression, particularly in the UL1 line. Kromastat shows a lighter blue signature overall and fails to systematically silence the cell cycle machinery (`light 🟥` retention) of the CNK89 cancer line at this concentration.

#### **3. Divergent Immune Modulation**
Romidepsin broadly represses immune-related signaling (shifting Allograft Rejection and Interferon Alpha responses to blue). Kromastat demonstrates a divergent functional fate by selectively preserving or activating specific immune nodes, resulting in distinct dark red hotspots like CD86. 

#### Final Verdict
Romidepsin operates as a direct cytotoxic agent with a broad immunosuppressive footprint. Kromastat functions as an immune-modulating agent, executing a partial anti-proliferative program while preserving the specific red biological signatures (such as CD86 expression) necessary for potential synergy with immunotherapies.