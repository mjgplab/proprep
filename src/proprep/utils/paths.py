"""Centralized package directory resolution for ProPrep.

Provides a single function to locate the proprep package root,
handling both source/pip/conda installs and PyInstaller bundles.
"""

import sys
from pathlib import Path


def get_package_dir() -> Path:
    """Return the proprep package root directory.

    In a PyInstaller bundle, returns the extracted bundle path.
    In a source/pip/conda install, returns the actual package directory.
    """
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        return Path(sys._MEIPASS) / "proprep"
    else:
        # This file is at proprep/utils/paths.py
        return Path(__file__).resolve().parent.parent
