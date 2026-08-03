#!/usr/bin/env python3
"""
Transformer Creator Package

A simplified, modular system for creating RedoxSite transformers.

This is the new implementation (formerly v2), now integrated into ProPrep.

NOTE: Do NOT import TransformerCreator from this __init__.py
The compatibility wrapper at ../transformer_creator.py should be used instead.
This allows the import path to work correctly:
    from .transformation.transformer_creator import TransformerCreator
"""

# Only export data models, not the TransformerCreator class
# The TransformerCreator should be imported from ../transformer_creator.py wrapper
from .data_models import (
    TransformerSpec,
    SiteRequirements,
    ComponentMatching,
    TransformationSequence,
    Parameters,
    create_empty_spec
)

__all__ = [
    # 'TransformerCreator',  # REMOVED - use wrapper instead
    'TransformerSpec',
    'SiteRequirements',
    'ComponentMatching',
    'TransformationSequence',
    'Parameters',
    'create_empty_spec'
]
