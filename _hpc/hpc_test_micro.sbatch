#!/bin/bash
#SBATCH --job-name=RNAseq_Test_Micro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=hpc_test_micro_log.txt
#SBATCH --error=hpc_test_micro_error.txt

# 1. Initialize Micromamba logic inside the subshell
# (This ensures it works even if your ~/.bashrc isn't sourced)
eval "$(micromamba shell hook --shell bash)"

# 2. Activate the environment on the 22T drive
micromamba activate cancer_rnaseq

# 3. Run the complete verification suite (Chr21)
# Using 'all' to cover both Upstream (BASH) and Downstream (R)
./test all
