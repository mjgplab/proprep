"""
Unified QM/MM Preparation Package

Provides a single entry point for Gaussian (ONIOM) and ORCA QM/MM input
generation with shared topology loading and trajectory frame extraction.
"""

from . import qmmm_preparator  # noqa: F401 — triggers @register_module
