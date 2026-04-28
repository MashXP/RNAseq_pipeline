# Introduction
**Alternative Splicing (AS)**, let's think of a gene not as a single recipe, but as a modular blueprint.

We have **exons** and **introns** -> many ways of removing introns or not removing them at all -> Generate **isoforms**

The isoforms have different functions

whether **Romi** or **Kromi** change _how_ the gene is put together
- Romi and Kromi are **HDAC inhibitors**, they modify how tightly DNA is packed, which can directly interfere with the speed and accuracy of the splicing machinery ⚠️

Gene in Apoptosis (2 isoform): 
- 1: promote survival
- 2: death -> drug make effect shift here maybe???
without changing the overall amount of the gene's mRNA

It no longer matter of **Activation** or **Suppression**
- If it's activated, we don't know if:
	- **Functional domain?** if lost, protein is useless
	- **Changes localization?** if it goes into nucleus instead of cytoplas?
	- **Switches function?** apoptotic genes are dual purpose. It can kill or promote cell activity

|**Category**|**DGE (Quantity)**|**AS (Quality)**|**Biological "Story"**|
|---|---|---|---|
|**I: Quantity Shift** 📈|Significant|Not Significant|The cell makes more or less of the _same_ protein.|
|**II: Quality Shift** 🛠️|Not Significant|Significant|The cell makes the _same amount_ of protein, but it's a _different version_.|
|**III: The Double Hit** 💥|Significant|Significant|The cell changes both how much is made and which version is produced.|
# Implementation
Focusing _only_ on genes that don't change in expression (Category 2) could indeed be seen as "cherry-picking" if we didn't acknowledge the broader context.

In a robust analysis, we don't choose one over the other. Instead, we use **Alternative Splicing (AS)** to add a second layer of evidence. The goal is to be more thorough, not less.

- **Category 3 (The Double Hit):** If a gene is downregulated _and_ starts producing a broken isoform, you have a "slam dunk" case for why that pathway is failing.
- **Category 2 (The Hidden Hit):** If a gene looks fine on the satellite map but is "broken" on the street level, you’ve discovered a mechanism that a standard DGE analysis would have missed.

## Metrics
FDR < 0.05  💚
- **False Discovery Rate:** expect no more than 5% of them to be false positives (noise). It’s our "confidence filter."
**Inclusion Level Difference (Δψ)**.
- **ψ (Psi)(Percent Spliced In) = 0.8:** 80% of the transcripts include that exon.💚
	- ψ is the probability that a specific exon is included in the final mRNA transcript.
		- **High ψ (e.g., 0.9 or 90%):** The exon is "spliced in" or included. The protein is likely built in its full, standard form. 
		- **Low ψ (e.g., 0.1 or 10%):** The exon is mostly skipped and removed along with the introns.
- **Δψ (Delta Psi) = 0.3:** The drug caused a 30% shift in how that exon is used compared to the control.💚
	- rMATS-turbo calculates **Δψ (Delta Psi)**, which is simply the difference: $\psi_{Treatment} ​ −\psi_{Control}$​.
		- If Δψ < 0, the drug is causing that exon to be skipped.
		- If Δψ > 0, the drug is forcing that exon to be included.

## **Functional Domain Mapping**.
To deduce _exactly_ **where the drug is hitting**, we perform a process called **Functional Domain Mapping**.

We don't just guess; we use specialized bioinformatic tools and databases to "zoom in" on the protein:
- **Genomic Coordinates**: rMATS-turbo tells us the exact start and end of the affected exon on the chromosome.
- **Protein Translation**: We map those coordinates to the corresponding amino acids in the protein sequence.
- **Domain Lookup**: We cross-reference that sequence with databases like **UniProt** or **Pfam**. 

- These databases act like an encyclopedia of protein "parts," telling us if a specific region is a binding site, a structural hinge, or a catalytic core.

## Data Preparation
Input: 
	**BAM files (Binary Alignment Map)**: contain **junction reads**
	**GTF files (Gene Transfer Format):** 
generates intermediate `.rmats` files that store these junction counts without performing any statistical comparisons yet

it can also detect **novel** splicing events ⚠️
	(Is it actually novel or just false positive? )

### Handling Replicates (The "Group" Approach)

In rMATS, you don't actually do individual 1-vs-1 comparisons (Replicate 1 vs Replicate 1, 1 vs 2, etc.). Instead, you provide **all triplicates** for a group at once, separated by commas.
- **Group 1 (Treatment):** Path/to/Kroma_Rep1.bam, Path/to/Kroma_Rep2.bam, Path/to/Kroma_Rep3.bam
- **Group 2 (Control):** Path/to/Control_Rep1.bam, Path/to/Control_Rep2.bam, Path/to/Control_Rep3.bam

### Kromastat vs. Romidepsin (The Head-to-Head)
- **Drug vs. Control:** Tells you what each drug does to the cell.
- **Drug A vs. Drug B:** Tells you the **difference in mechanism**. For example, if Romidepsin causes a specific splicing error that Kromastat doesn't, that might explain why Romidepsin is more effective in certain cell lines.

| **Feature**                | **Species Pooling (Global)**                                | **Cell-Line Level (Specific)**                                      |
| -------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------- |
| **Statistical Power**      | Higher $n$ (more samples), but higher "noise."              | Lower $n$ (3 replicates), but very "clean" data.                    |
| **Focus**                  | Finds common drug targets.                                  | Finds phenotype-specific mechanisms.                                |
| **Presentation Narrative** | "This drug affects the human splicing machinery generally." | "This drug breaks the machinery specifically in aggressive cancer." |
-> Should pick Cell-Line level.
## Scripts
> [!todo]
> 
> **Prompt:** "Write a Python script that parses a metadata CSV containing sample IDs, treatments, and cell lines. For each cell line (UL1, CNK89, H9, SUPM2), it should identify comparison pairs (e.g., Kromastat_6nM vs. DMSO_Kromastat). For each pair, it must generate two text files: `b1.txt` (Treatment) and `b2.txt` (Control).
> 
> Each text file must contain a single line of comma-separated paths to the BAM files. The BAM paths follow this pattern: `/RNAseq_pipeline/_data/bam/{ID}/{ID}_Aligned.sortedByCoord.out.bam`, where {ID} is the number at the start of the 'File1' column (e.g., '27'). The script should group the 3 biological replicates for each condition into their respective files."
> 

> [!todo]
> **Prompt:** "Write a 2 Bash script to execute an rMATS-turbo pipeline.
> 1. First, call the Python metadata script to generate all `b1.txt` and `b2.txt` files for the 8 comparisons (4 cell lines x 2 drugs).   
> 2. Then, loop through each comparison to run the **rMATS 'Prep' step** using the flag. It must use the correct species-specific GTF file (Human vs. Canine).
> 3. Finally, run the **'Post' step** using the Post script.
>
> Include these specific rMATS parameters:
>- **Statistical Filter**: FDR ≤0.01 and ∣Δψ∣≥0.05.
>- **Read Coverage**: Minimum average of 10 reads per group.
>- **Read Type**: Paired-end (`-t paired`).
>- **Library Type**: Stranded


# Analysis
## **The 5 Splicing Flavors**
Decoding the specific patterns rMATS-turbo detects, such as **Skipped Exons (SE)** or **Retained Introns (RI)**, and what they mean for protein function .
### The 5 Basic Patterns

- **Skipped Exon (SE)** 🦘: Most common flavor where an exon is either included in the final mRNA or completely skipped. 
	- If that exon contained a crucial part of a protein, like a binding site, the protein might still be made but be totally non-functional.
- **Alternative 5' or 3' Splice Sites (A5SS/A3SS)** 📏: Instead of skipping the whole exon, it uses a different "cutting point" at either the start or the end. 
	- This makes the exon slightly longer or shorter, which can subtly shift the protein's shape.
- **Mutually Exclusive Exons (MXE)** ↔️: The cell must choose between two different versions of an exon; it can include Exon A or Exon B, but never both in the same transcript.
- **Retained Intron (RI)** 🛑: This is often the most "destructive" flavor. An intron—which is usually non-coding "junk" mRNA—is kept in the final message.

**HDAC inhibitors** like Romidepsin can cause the splicing machinery to "trip" over certain sequences. For example, if a drug causes a **Retained Intron**, it often introduces a "premature stop codon."

## **Functional Finish Line**
Visualizing the data with **Sashimi plots** and connecting these splicing "crashes" to pathways like cell proliferation

A **Sashimi plot** is the "gold standard" for visualizing alternative splicing because it shows two things at once:
1. **Read Density** 🏔️: The height of the "mountains" tells us how much of an exon was actually expressed.
2. **Junction Curves** 🌉: The curved lines (arcs) between exons tell us how many RNA-seq reads bridged those exons together. The number on the arc represents the count of these "junction reads."