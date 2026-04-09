#!/bin/bash

# setup_env.sh
# Automates the creation of the MiaPAca-2 cancer RNA-seq environment.

set -e

ENV_NAME="cancer_rnaseq"

echo "========================================================="
echo "   MiaPAca-2 Cancer RNA-seq - Environment Setup"
echo "========================================================="

# 1. Detect Mamba or Conda
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "Error: Neither mamba nor conda was found."
    echo "Please install Mambaforge or Miniconda first: https://mamba.readthedocs.io/en/latest/installation/mamba-install.html"
    exit 1
fi

echo "Using $CONDA_CMD to create the environment..."

# 2. Create the environment
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
