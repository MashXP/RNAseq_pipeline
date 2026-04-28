#!/usr/bin/env python3
"""
generate_bam_lists.py
---------------------
Parses drPhuong_Sample_Data_Table.csv and generates rMATS-turbo input files
(b1.txt = Treatment, b2.txt = Control) for each cell-line comparison pair.

Output structure:
  altsplicing/bam_lists/{CellLine}/{Treatment}_vs_{Control}/b1.txt
  altsplicing/bam_lists/{CellLine}/{Treatment}_vs_{Control}/b2.txt
"""

import csv
import os
import sys
from pathlib import Path
from collections import defaultdict

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

# Script is in scripts_upstream/utils/, so root is 3 levels up
PIPELINE_ROOT = Path(__file__).resolve().parent.parent.parent

METADATA_CSV  = PIPELINE_ROOT / "drPhuong_Sample_Data_Table.csv"
OUTPUT_BASE   = PIPELINE_ROOT / "_data" / "rMATS-turbo" / "bam_lists"
BAM_TEMPLATE  = str(PIPELINE_ROOT / "_data" / "bam" / "{id}" / "{id}_Aligned.sortedByCoord.out.bam")

# Cell lines of interest and the comparison pairs they support
# Format: (treatment_label, control_label)
CELL_LINE_COMPARISONS = {
    "UL1":   [("Romidepsin_6nM", "DMSO_Romidepsin"),
              ("Kromastat_6nM",  "DMSO_Kromastat"),
              ("Kromastat_6nM",  "Romidepsin_6nM")],
    "CNK89": [("Romidepsin_6nM", "DMSO_Romidepsin"),
              ("Kromastat_6nM",  "DMSO_Kromastat"),
              ("Kromastat_6nM",  "Romidepsin_6nM")],
    "H9":    [("Romidepsin_6nM", "DMSO_Romidepsin"),
              ("Kromastat_6nM",  "DMSO_Kromastat"),
              ("Kromastat_6nM",  "Romidepsin_6nM")],
    "SUPM2": [("Romidepsin_6nM", "DMSO_Romidepsin"),
              ("Kromastat_6nM",  "DMSO_Kromastat"),
              ("Kromastat_6nM",  "Romidepsin_6nM")],
}

# ──────────────────────────────────────────────────────────────────────────────
# PARSE CSV
# ──────────────────────────────────────────────────────────────────────────────

def parse_metadata(csv_path: Path) -> dict:
    """
    Returns a nested dict:  data[cell_line][treatment] = [bam_path, ...]
    Skips the title row ("Sample details...") and any blank rows.
    The real header is the second row: File1, File2, Treatment, Cell line, Group
    """
    data = defaultdict(lambda: defaultdict(list))

    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh)
        header_found = False
        for row in reader:
            # Skip completely blank rows
            if not any(c.strip() for c in row):
                continue
            # First real row is the title; second is the column header
            if not header_found:
                if row[0].strip() == "File1":
                    header_found = True
                continue   # skip title row

            # Stop when we hit the genome-file section
            if row[0].strip() in ("File", "Genome files (Human):", "Genome files (Canine):"):
                break
            if not row[0].strip().endswith(".fastq.gz"):
                continue

            file1     = row[0].strip()   # e.g. "27_R1_001.fastq.gz"
            treatment = row[2].strip()   # e.g. "Romidepsin_6nM"
            cell_line = row[4].strip()   # e.g. "UL1"  (the Group column)

            # Extract the numeric ID from the filename
            sample_id = file1.split("_")[0]   # "27"

            bam_path = BAM_TEMPLATE.format(id=sample_id)
            
            # CRITICAL: Verify BAM existence
            if not os.path.exists(bam_path):
                print(f"[ERROR] BAM file not found: {bam_path}", file=sys.stderr)
                print(f"[ERROR] Sample ID: {sample_id}, Treatment: {treatment}, Cell line: {cell_line}", file=sys.stderr)
                sys.exit(1)

            data[cell_line][treatment].append(bam_path)

    return data

# ──────────────────────────────────────────────────────────────────────────────
# WRITE FILES
# ──────────────────────────────────────────────────────────────────────────────

def write_bam_list(path: Path, bam_paths: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(",".join(bam_paths) + "\n")

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print(f"[INFO] Reading metadata: {METADATA_CSV}")
    if not METADATA_CSV.exists():
        print(f"[ERROR] Metadata CSV not found: {METADATA_CSV}", file=sys.stderr)
        sys.exit(1)

    data = parse_metadata(METADATA_CSV)

    generated = []
    errors     = []

    for cell_line, comparisons in CELL_LINE_COMPARISONS.items():
        if cell_line not in data:
            errors.append(f"  [WARN] Cell line '{cell_line}' not found in CSV.")
            continue

        for treatment_label, control_label in comparisons:
            comp_name = f"{treatment_label}_vs_{control_label}"
            out_dir   = OUTPUT_BASE / cell_line / comp_name

            b1_paths = data[cell_line].get(treatment_label, [])
            b2_paths = data[cell_line].get(control_label,   [])

            # Validation
            ok = True
            if not b1_paths:
                errors.append(f"  [WARN] No samples found for {cell_line} / {treatment_label}")
                ok = False
            if not b2_paths:
                errors.append(f"  [WARN] No samples found for {cell_line} / {control_label}")
                ok = False
            if len(b1_paths) != 3 or len(b2_paths) != 3:
                errors.append(
                    f"  [WARN] Expected 3 replicates for {cell_line}/{comp_name}, "
                    f"got b1={len(b1_paths)}, b2={len(b2_paths)}"
                )

            if not ok:
                continue

            b1_file = out_dir / "b1.txt"
            b2_file = out_dir / "b2.txt"
            write_bam_list(b1_file, b1_paths)
            write_bam_list(b2_file, b2_paths)
            generated.append((cell_line, comp_name, b1_paths, b2_paths))

    # ── Summary report ──────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  rMATS BAM list generation — summary")
    print(f"{'='*70}")
    print(f"  Output directory : {OUTPUT_BASE}\n")

    for cell_line, comp_name, b1, b2 in generated:
        print(f"  [{cell_line}]  {comp_name}")
        print(f"    b1.txt  (Treatment): {', '.join(os.path.basename(p) for p in b1)}")
        print(f"    b2.txt  (Control)  : {', '.join(os.path.basename(p) for p in b2)}")
        print()

    if errors:
        print("Warnings / Errors:")
        for e in errors:
            print(e)

    print(f"  Total comparisons generated: {len(generated)}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
