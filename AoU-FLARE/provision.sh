#!/bin/bash

set -euo pipefail

echo "===== Starting provision.sh ====="

# -------------------------------
# Step 1: Download FLARE
# -------------------------------
echo "Downloading FLARE..."
wget -O flare.jar https://faculty.washington.edu/browning/flare.jar
echo "FLARE downloaded as flare.jar"

# -------------------------------
# Step 2: Download and setup GNU Parallel
# -------------------------------
echo "Downloading GNU Parallel..."
wget -qO parallel.tar.bz2 https://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
tar -xf parallel.tar.bz2

# Detect directory name
PARALLEL_DIR=$(find . -maxdepth 1 -type d -name "parallel-*")
PARALLEL_PATH="$PWD/$PARALLEL_DIR/src"

if [[ -x "$PARALLEL_PATH/parallel" ]]; then
    echo "GNU Parallel successfully extracted."
    export PATH="$PARALLEL_PATH:$PATH"
    echo "export PATH=$PARALLEL_PATH:\$PATH" >> ~/.bashrc
    parallel --version
else
    echo "GNU Parallel installation failed."
    exit 1
fi

# -------------------------------
# Step 3: Install Java 21 via Conda
# -------------------------------
echo "Setting up Java 21 environment using Conda..."

# Ensure LD_LIBRARY_PATH includes libarchive (AoU environment)
export LD_LIBRARY_PATH=/home/jupyter/miniconda/pkgs/libarchive-3.7.4-hfab0078_0/lib:$LD_LIBRARY_PATH

# Set up Conda package cache directory
export CONDA_PKGS_DIRS=$HOME/.conda_pkgs/
mkdir -p "$CONDA_PKGS_DIRS"

# Create and activate Java 21 environment
conda create -n java21_env openjdk=21 -y
conda init bash
source ~/.bashrc
conda activate java21_env

# Verify Java
echo "Java version:"
java -version

# -------------------------------
echo "===== provision.sh completed successfully ====="
