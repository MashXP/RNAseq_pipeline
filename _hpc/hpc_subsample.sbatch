#!/bin/bash
#SBATCH --job-name=hpc_subsample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=subsample_%j.out
#SBATCH --error=subsample_%j.err

# hpc_subsample.sh
# Run this on the HPC inside the directory containing 'data/drPhuong'
# It creates a flat folder 'test_subset' with randomly sampled reads per sample.

# Set the base directory (where data/ resides)
DATA_ROOT="/mnt/22T/vuphong/cancer_rna_drPhuong_data"
OUTPUT_DIR="/mnt/22T/vuphong/cancer_rna_drPhuong_data/test_subset"
NUM_READS=25000 # Lowered read count to 25k to keep test suite fast

if [ ! -d "$DATA_ROOT" ]; then
    echo "Error: Directory $DATA_ROOT not found."
    exit 1
fi

command -v seqtk >/dev/null 2>&1 || { echo >&2 "Error: seqtk is required but not installed. Please notify the administrator or install it manually."; exit 1; }

mkdir -p "$OUTPUT_DIR"

echo "=== Starting Subsampling for Verification ==="
echo "Output directory: $OUTPUT_DIR"
echo "Target reads per file: $NUM_READS"

# Specifically handling paired reads using seqtk with a fixed seed (-s100)
# We find all R1 files and automatically pair them with R2
find "$DATA_ROOT" -name "*_R1_001.fastq.gz" | while read -r r1_file; do
    r2_file="${r1_file/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    r1_name=$(basename "$r1_file")
    r2_name=$(basename "$r2_file")
    
    echo "Processing pair $r1_name / $r2_name..."
    
    if [ ! -f "$r2_file" ]; then
        echo "Warning: No matching R2 found for $r1_file"
        continue
    fi
    
    # Subsample exactly matching lines due to fixed random seed
    seqtk sample -s100 "$r1_file" $NUM_READS | gzip > "$OUTPUT_DIR/$r1_name"
    seqtk sample -s100 "$r2_file" $NUM_READS | gzip > "$OUTPUT_DIR/$r2_name"
done

echo "=== Subsampling Complete ==="
echo "You can now download the '$OUTPUT_DIR' folder to your local machine."
echo "Suggested download command (from local):"
echo "  rsync -avz [USER]@[HOST]:$OUTPUT_DIR/ /path/to/RNAseq_pipeline/_data/fastq/"
