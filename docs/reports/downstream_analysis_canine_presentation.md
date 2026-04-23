---
marp: true
theme: default
paginate: true
backgroundColor: #fff
style: |
  section {
    font-family: 'Inter', sans-serif;
  }
  h1 {
    color: #0f172a;
    border-bottom: 2px solid #38bdf8;
    padding-bottom: 10px;
  }
  h2 {
    color: #1e293b;
    border-bottom: 1px solid #94a3b8;
  }
  h3 {
    color: #475569;
  }
  blockquote {
    background: #f8fafc;
    border-left: 10px solid #38bdf8;
    margin: 1.5em 10px;
    padding: 0.5em 10px;
    font-style: italic;
  }
  .pending-tag {
    background: #fef3c7;
    color: #92400e;
    padding: 5px 15px;
    border-radius: 20px;
    font-size: 0.7em;
    font-weight: bold;
    border: 1px solid #f59e0b;
    display: inline-block;
    margin-bottom: 10px;
  }
  .verified-tag {
    background: #dcfce7;
    color: #166534;
    padding: 5px 15px;
    border-radius: 20px;
    font-size: 0.7em;
    font-weight: bold;
    border: 1px solid #22c55e;
    display: inline-block;
    margin-bottom: 10px;
  }
  footer {
    color: #94a3b8;
  }
---

# Downstream Transcriptomic Analysis (Canine)
## Verified Audit: Romidepsin vs. Kromastat

**Status**: Technical Verification Phase (LFC > 2.0)
**Objective**: Dissecting functional conservation and divergence in canine lymphoid models.

---

## 1. Principal Component Analysis (PCA)
### Global Transcriptomic Variance

![bg right:60% contain](../../results/canine/figures/04_01_pca_combined.png)

- <div class="verified-tag">AUDITED</div>
- **PC1 (96%)**: Lineage bifurcation (UL1 vs CNK89).
- **PC2 (3%)**: Consistent Drug Effect.
- High technical fidelity across technical replicates.

---

## 1.3 Sub-PCA: Lineage-Specific Response
### Removing Background Variance

| CNK89 (NK-cell Cancer) | UL1 (T-ALL Cancer) |
| :---: | :---: |
| ![height:400px](../../results/canine/figures/04_01_pca_CNK89.png) | ![height:400px](../../results/canine/figures/04_01_pca_UL1.png) |

- Drug effect becomes primary driver (**>90% variance**) in local analysis.

---

## 2.1 Global Response Magnitude
### Verified Activation vs. Suppression

<p align="center">
  <img src="../../results/canine/figures/04_02_volcano_combined.png" width="700px"/>
</p>

| Treatment | Total DEGs (|LFC| > 2) | Ratio (Up:Down) |
| :--- | :---: | :---: |
| **Romidepsin (6nM)** | **1,474** | 11.2 : 1 |
| **Kromastat (6nM)** | **147** | **146.0 : 1** |

---

## 2.2 Global Head-to-Head
### Romidepsin vs. Kromastat (6nM)

![bg left:60% contain](../../results/canine/figures/04_02_volcano_Romi_6nM_vs_Kromastat_6nM.png)

- <div class="verified-tag">AUDITED</div>
- **10 Significant Differences** identified at 6nM.
- **Top Representatives**:
  - Higher in Romi: `ENPP2`, `ELOVL2`
  - Higher in Kroma: `CRLF2`, `ARK2C`
- **Verdict**: Inhibitors are functionally distinct in the canine context.

---

## 2.3 Local Response: UL1 vs CNK89
### The Reactivity Gap

| UL1 (T-ALL Cancer) | CNK89 (NK-cell Cancer) |
| :---: | :---: |
| ![height:350px](../../results/canine/figures/04_02_volcano_UL1_Romi_6nM_vs_DMSO_Romi.png) | ![height:350px](../../results/canine/figures/04_02_volcano_CNK89_Romi_6nM_vs_DMSO_Romi.png) |
| **2,254 DEGs** (Romi) | **1,015 DEGs** (Romi) |

- **UL1** is significantly more reactive to HDAC inhibition than **CNK89**.

---

## 3. Gene Overlap Analysis
### Verified Core Identity

| All DEGs | Up-regulated | Down-regulated |
| :---: | :---: | :---: |
| ![height:300px](../../results/canine/figures/04_03_venn_all.png) | ![height:300px](../../results/canine/figures/04_03_venn_up.png) | ![height:300px](../../results/canine/figures/04_03_venn_down.png) |

- **94%** of Kromastat's activation is a subset of Romidepsin.
- **10.4x** Representation Factor confirms shared activation is highly non-random.

---

## 4.1 Functional Response Heatmaps
### Signature Matrix Preview

![bg right:60% contain](../../results/canine/figures/04_04_heatmap_pathway_Romi_6nM_vs_DMSO_Romi.png)

- <div class="pending-tag">AUDIT PENDING</div>
- Preliminary scan indicates broad proliferation crash.
- Formal pathway re-verification in progress.

---

## 5. Top Variable Genes
### Lineage vs. Drug Variance
---
![bg left:60% contain](../../results/canine/figures/04_05_heatmap_top_variable.png)

- <div class="pending-tag">AUDIT PENDING</div>
- Dominant divide between UL1 and CNK89 preserved at 6nM.
- Audit focusing on the degree of lineage-signature preservation under Kromastat.

---

## 6. GSEA Hallmark NES (Combined)
### Quantitative Rankings
---
![bg right:60% contain](../../results/canine/figures/04_06_hallmark_nes_combined.png)

- <div class="pending-tag">AUDIT PENDING</div>
- Rankings provide evidence for the "Immune Sparing" phenotype.
- Waiting for final leading-edge verification.

---

## 8.1 UpSet Consistency Matrix
### Intersectional Degree of Consensus
---
![bg center:90% contain](../../results/canine/figures/04_09_upset_consistency.png)

---

## 8.2 LFC Correlation
### Directional Transcriptomic Consensus
---
![bg center:85% contain](../../results/canine/figures/04_10_lfc_correlation.png)

---

# 9. Alluvial Analysis
## Cross-Species Phenotypic Conservation

- <div class="verified-tag">NEW: COMPARATIVE ENGINE</div>
- **Goal**: Trace the "fate" of Hallmark pathways across species and cell lines.
- **Conservation**: Identifying core HDAC inhibition signatures that transcend biological boundaries.

> [!NOTE]
> Alluvial flows visualize how pathway status (Activated/Suppressed/NS) shifts between conditions.

---

## 9.1 Romidepsin Flow (Human vs. Canine)
<!-- ### Aggressive Cancer Bridge (SUPM2 vs. CNK89) -->

![bg right:60% contain](../../results/comparative/figures/04_11_alluvial_comp_Romi_SUPM2_vs_CNK89.png)

- **High Signal Conservation**: Core pathways (E2F, G2M, MYC) remain "Locked" in Suppressed status across species.
- **Divergent Nodes**: Identifying species-specific immune activation patterns.

---

## 9.2 Master Romidepsin Panorama
### Healthy -> Cancer | Human -> Canine
---
![bg center:90% contain](../../results/comparative/figures/04_11_alluvial_master_romidepsin.png)

---

## 9.3 Master Drug-Species Bridge
### Aggressive Cancer Lines Comparison
---
![bg center:90% contain](../../results/comparative/figures/04_11_alluvial_master_bridge_aggressive.png)

---

# 10. Comprehensive Summary
### Biological Audit Conclusion

- <div class="verified-tag">DATA VERIFIED</div>
- **Shared Core**: 94% Activation Subset (Conserved HDAC response).
- **Potency Gap**: Romidepsin is a ~10x broader modulator by DEG count.
- **Suppression Program**:
  - **Romidepsin**: Aggressive gene silencing (121 genes).
  - **Kromastat**: Minimal silencing (1 gene) at high stringency.

**Action**: Finalize Section 4 & 5 interpretations to confirm specialized immune preservation in CNK89 line.
