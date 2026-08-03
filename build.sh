#!/bin/bash
# ProPrep build script
# Builds a standalone executable using PyInstaller

set -e

echo "=== ProPrep Build ==="

# Check if we're in a conda environment
if [ -n "$CONDA_DEFAULT_ENV" ]; then
    echo "Using conda environment: $CONDA_DEFAULT_ENV"
else
    echo "No conda environment detected. Using system Python."
fi

# Install package in development mode (ensures all imports work)
echo "Installing proprep in development mode..."
python3 -m pip install -e . --quiet

# Run the build
echo "Running build script..."
python3 build.py
