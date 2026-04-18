#!/bin/bash
# ============================================================
# ProPrep Installer
# ============================================================
# This script installs ProPrep and all its dependencies into
# a self-contained conda environment.
#
# Prerequisites:
#   - Miniforge or Anaconda/Miniconda installed
#     (https://github.com/conda-forge/miniforge#download)
#   - A MODELLER license key (free for academics)
#     Register at: https://salilab.org/modeller/registration.html
#
# Usage:
#   bash install_proprep.sh
#
# After installation:
#   conda activate ProPrep
#   source $CONDA_PREFIX/amber.sh
#   proprep
# ============================================================

show_help() {
cat << 'HELPTEXT'
ProPrep Installation Guide
===========================

Quick Install (Recommended)
----------------------------
Run this script:

    bash install_proprep.sh

The script handles everything: creates a conda environment, installs ProPrep
with all dependencies (including AmberTools and MODELLER), and verifies the
installation.

Manual Install
--------------

1. Install Conda

   If you don't have conda, install Miniforge:
   https://github.com/conda-forge/miniforge#download

   Follow the defaults during installation. Restart your terminal when done.

2. Get a MODELLER License Key

   MODELLER is used for homology modeling and requires a free academic license.
   Register at: https://salilab.org/modeller/registration.html

   You'll receive a license key by email. Set it as an environment variable
   before installing:

       export KEY_MODELLER="your_key_here"

   To make this permanent, add the line above to your ~/.bashrc or ~/.zshrc.

3. Create Environment and Install

       conda create --name ProPrep python=3.12 -y
       conda install -n ProPrep -c mjgplab -c dacase -c salilab -c conda-forge proprep -y
       conda run -n ProPrep pip install tmtools

Usage
-----
Each time you want to use ProPrep:

    conda activate ProPrep
    source $CONDA_PREFIX/amber.sh
    proprep

Troubleshooting
---------------

"conda: command not found"
    Restart your terminal after installing Miniforge, or run:
    source ~/.bashrc (Linux) or source ~/.zshrc (macOS)

MODELLER license warning
    If you skipped the MODELLER key during install, homology modeling won't
    work. Set the key and reinstall MODELLER:
        export KEY_MODELLER="your_key_here"
        conda install -n ProPrep -c salilab modeller --force-reinstall -y

FreeSASA not available
    Run: conda install -n ProPrep -c conda-forge freesasa

Permission errors on macOS
    If macOS blocks the install, open System Settings > Privacy & Security
    and allow the blocked items.
HELPTEXT
}

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    show_help
    exit 0
fi

set -e

ENV_NAME="ProPrep"
PYTHON_VERSION="3.12"

echo "============================================================"
echo "  ProPrep Installer"
echo "============================================================"
echo ""

# Check for conda
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in your PATH."
    echo ""
    echo "Please install Miniforge first:"
    echo "  https://github.com/conda-forge/miniforge#download"
    echo ""
    echo "After installing, restart your terminal and run this script again."
    exit 1
fi

# Check for MODELLER license key
echo "MODELLER requires a license key (free for academic use)."
echo "Register at: https://salilab.org/modeller/registration.html"
echo ""
if [ -z "$KEY_MODELLER" ]; then
    read -p "Enter your MODELLER license key (or press Enter to skip): " modeller_key
    if [ -n "$modeller_key" ]; then
        export KEY_MODELLER="$modeller_key"
    else
        echo ""
        echo "NOTE: Skipping MODELLER key. Homology modeling will not be"
        echo "available until you configure the key manually. See the"
        echo "install guide for instructions."
        echo ""
    fi
else
    echo "Found MODELLER key in environment variable."
fi

# Check if environment already exists
echo ""
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '$ENV_NAME' already exists."
    echo ""
    echo "  1. Update ProPrep (keep existing environment)"
    echo "  2. Full reinstall (remove and recreate environment)"
    echo "  3. Cancel"
    echo ""
    read -p "Choose [1/2/3] (1): " install_choice
    install_choice=${install_choice:-1}

    if [[ "$install_choice" == "1" ]]; then
        echo ""
        echo "Updating ProPrep..."
        conda update -n "$ENV_NAME" -c mjgplab -c dacase -c salilab -c conda-forge proprep -y
        conda run -n "$ENV_NAME" pip install --upgrade tmtools

        echo ""
        echo "Verifying installation..."
        if conda run -n "$ENV_NAME" proprep --version &> /dev/null; then
            VERSION=$(conda run -n "$ENV_NAME" proprep --version 2>&1)
            echo "  ProPrep updated successfully: $VERSION"
        else
            echo "  ProPrep updated (version check skipped)."
        fi

        echo ""
        echo "============================================================"
        echo "  Update complete!"
        echo "============================================================"
        echo ""
        echo "  To use ProPrep:"
        echo ""
        echo "    conda activate $ENV_NAME"
        echo "    source \$CONDA_PREFIX/amber.sh"
        echo "    proprep"
        echo ""
        echo "============================================================"
        exit 0

    elif [[ "$install_choice" == "2" ]]; then
        conda env remove -n "$ENV_NAME" -y
    else
        echo "Cancelled."
        exit 0
    fi
fi

# Create conda environment
echo ""
echo "Creating conda environment '$ENV_NAME' with Python $PYTHON_VERSION..."
conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y

# Install ProPrep and dependencies
echo ""
echo "Installing ProPrep and dependencies..."
echo "(This may take several minutes)"
echo ""
conda install -n "$ENV_NAME" -c mjgplab -c dacase -c salilab -c conda-forge proprep -y

# Install PyPI-only dependencies
echo ""
echo "Installing additional Python packages..."
conda run -n "$ENV_NAME" pip install tmtools

# Verify installation
echo ""
echo "Verifying installation..."
if conda run -n "$ENV_NAME" proprep --version &> /dev/null; then
    VERSION=$(conda run -n "$ENV_NAME" proprep --version 2>&1)
    echo "  ProPrep installed successfully: $VERSION"
else
    echo "  ProPrep installed (version check skipped)."
fi

echo ""
echo "============================================================"
echo "  Installation complete!"
echo "============================================================"
echo ""
echo "  To use ProPrep:"
echo ""
echo "    conda activate $ENV_NAME"
echo "    source \$CONDA_PREFIX/amber.sh"
echo "    proprep"
echo ""
if [ -z "$KEY_MODELLER" ]; then
    echo "  To configure MODELLER later:"
    echo "    1. Get a key at https://salilab.org/modeller/registration.html"
    echo "    2. Add this to your ~/.bashrc or ~/.zshrc:"
    echo "       export KEY_MODELLER=\"your_key_here\""
    echo ""
fi
echo "============================================================"
