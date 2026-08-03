"""
Redox Site Preparation Module

This package provides functionality for detecting, analyzing, and preparing
redox sites in protein structures for molecular dynamics simulations.
"""

__version__ = "2.0.0"

# Import modern redox site preparation components
from .redoxsite_integration import RedoxSitePreparationModule

__all__ = [
    "RedoxSitePreparationModule",
]
