# Upstream Dissection: 05_multiqc.sh

This script is the final "Audit Trail" of your genome mapping. It "vacuum cleans" all the hundreds of log files from STAR, Picard, and featureCounts and turns them into a single, beautiful HTML report.

---

## 1. Output Directory Initialization
```bash
MULTIQC_DIR="$DATA_DIR/multiqc"
...
mkdir -p "$MULTIQC_DIR"
```
- **The Job**: Creates a central folder for the final report.
- **The Reasoning**: MultiQC produces a standalone HTML file AND a data directory. Storing them in a dedicated folder keeps your project root clean and publishable.

---

## 2. MultiQC Aggregation Call
```bash
multiqc "$BAM_DIR" "$QC_DIR" "$QC_RAW_DIR" "$COUNTS_DIR" \
    -o "$MULTIQC_DIR" \
    -f -n "rna_seq_pipeline_summary" \
    -c "$BASE_DIR/../multiqc_config.yaml"
```
- **The Job**: This is the heart of the script. It scans all the subfolders where our previous scripts left their logs.
- **Key Parameters**:
    - `-o "$MULTIQC_DIR"`: Directs the final HTML report to its folder.
    - `-f`: "Force" mode—it will overwrite any previous versions of the report. This is useful if you are re-running a sample to see if a quality issue improved.
    - `-c "...multiqc_config.yaml"`: 
        - **Reasoning**: This configuration file defines the "Look and Feel." It ensures that your Krom vs. Romi samples are grouped in a logical way for your mentor to evaluate easily.

---

## 3. Data Inputs (What the script catches)
- **FastQC**: From Step 01. Shows quality before and after trimming.
- **STAR**: From Step 02. Shows mapping rates.
- **Picard**: From Step 03. Shows EXON/INTRON distribution and RNA health. (Note: RSeQC was originally in this script but is no longer used).
- **featureCounts**: From Step 04. Shows counting summaries and biotype distribution.

---

## 4. Final Output Analysis
- **HTML Report**: `_data/multiqc/rna_seq_pipeline_summary.html`.
- **The Reasoning**: This file is your "Final Exam." Before you look at a single gene in R, you must look at this report. 
- **Success Criteria**:
    - **General Stats**: Uniquely Mapped % > 80%?
    - **Biotype Distribution**: High "Protein Coding" levels?
    - **Exon/Intron**: High "Exon" mapping?

If all these lights are green, you are officially ready for the **Downstream bridge** and R analysis.
