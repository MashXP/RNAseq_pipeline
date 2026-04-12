#!/bin/bash

# setup_env.sh
# Automates the creation of the RNA-seq pipeline environment.

set -e

ENV_NAME="cancer_rnaseq"

echo "========================================================="
echo "   RNA-seq Pipeline - Environment Setup"
echo "========================================================="

# 1. Detect Micromamba, Mamba or Conda
if command -v micromamba &> /dev/null; then
    CONDA_CMD="micromamba"
elif command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "Error: Neither micromamba, mamba, nor conda was found."
    echo "Please install one of them first:"
    echo "  - Micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
    echo "  - Mambaforge: https://mamba.readthedocs.io/en/latest/installation/mamba-install.html"
    exit 1
fi

echo "Using $CONDA_CMD to create the environment..."

# 2. Create the environment
# Note: 'env create' works for conda/mamba; micromamba supports it for compatibility.
$CONDA_CMD env create -f environment.yml --name "$ENV_NAME"

echo ""
echo "========================================================="
echo "   Environment Setup SUCCESSFUL"
echo "========================================================="
echo ""
echo "To activate the environment, run:"
echo "    $CONDA_CMD activate $ENV_NAME"
echo ""
echo "After activation, you can run the pipeline with:"
echo "    ./run all"
echo "========================================================="
