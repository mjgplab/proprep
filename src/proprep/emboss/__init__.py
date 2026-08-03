# src/proprep/emboss/__init__.py
"""
EMBOSS Analysis Package

Comprehensive EMBOSS sequence analysis tools for ProPrep.
"""

from .emboss_module import EMBOSSModule
from .emboss_core import EMBOSSCore
from .sequence_input import SequenceInputManager

__all__ = ['EMBOSSModule', 'EMBOSSCore', 'SequenceInputManager']
