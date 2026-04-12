#!/bin/bash
#SBATCH --job-name=Download_RNASeq_Data
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=hpc_download_log_%j.txt
#SBATCH --error=hpc_download_error_%j.txt

echo "========================================================="
echo "   Starting Google Drive Download job"
echo "========================================================="

# Initialize Micromamba logic inside the subshell
eval "$(micromamba shell hook --shell bash)"

# Activate the environment on the 22T drive
micromamba activate

# Run the download script
bash /mnt/22T/vuphong/cancer_rna_drPhuong_data/download_data.sh

echo "========================================================="
echo "   Download Job COMPLETED"
echo "========================================================="
