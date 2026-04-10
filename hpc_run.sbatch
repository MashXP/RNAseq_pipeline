#!/bin/bash
#SBATCH --job-name=MiaPAca_Prod_Run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=hpc_run_log_%j.txt
#SBATCH --error=hpc_run_error_%j.txt

# 1. Initialize Micromamba (Safe method for batch shells)
eval "$(micromamba shell hook --shell bash)"

# 2. Activate the production environment
micromamba activate cancer_rnaseq

# 3. RUN THE PIPELINE (Everything)
# This will run Upstream (Trimming, Aligning, Quantifying) 
# and Downstream (DESeq2, Enrichment, Visualization)
./run all

echo "========================================================="
echo "   MiaPAca-2 Production Run COMPLETED SUCCESSFULLY"
echo "========================================================="
