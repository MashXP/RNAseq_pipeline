# TODO

- Refactor the code structure to centralize variables for more future pipelines with more csv. 
    - Reducing hardcodes (e.g. species, group names, etc.)
    - Might consider using a config file (e.g. YAML, JSON) to store these variables
    - Include options to add more complex comparision.
- Allow swappable tool modules (e.g. STAR vs HISAT2, DESeq2 vs edgeR, etc.)
- More thorough clean script.

## COMPLETED ✅
- **Alternative Splicing Pipeline**: Integrated rMATS-turbo for upstream and downstream analysis.
- **Documentation**: Added comprehensive guides for upstream/downstream AS analysis and updated theory notes.
- **Standardization**: Standardized Romidepsin nomenclature and enhanced plot aesthetics across the pipeline.
- **Canine Nomenclature**: Replaced all "DOG" instances with "Canine".

## URGENT - DONE
- Align Dr VuHuynh scripts to current pipeline.
    - Currently, Dr. VuHuynh's upstream scripts have different alignment outputs than the current pipeline: The only difference is the trimming process, which might not be an issue. (FIXED)
    - Need to investigate and make them consistent for paper publification. Then restructure code to merge downstream "bridge" script for less confusion. (done)
    - Re-audit heatmaps based on the specified principles and structured audit format. (no need)
    - Replace All DOG instance to Canine (done)