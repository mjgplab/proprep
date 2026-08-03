# Phase 1: Update imports and base class references

# Update import paths
s|from proprep\.metallo_prep\.transformation\.transformers\.base import MetalSiteTransformer|from proprep.redoxsite_prep.transformation.redox_transformer_framework import RedoxSiteTransformerBase, register_redox_transformer|g
s|from proprep\.metallo_prep\.transformation|from proprep.redoxsite_prep.transformation|g
s|proprep\.metallo_prep\.transformation|proprep.redoxsite_prep.transformation|g

# Update base class
s|class \([^ ]*\)(MetalSiteTransformer):|class \1(RedoxSiteTransformerBase):|g
s|(MetalSiteTransformer)|(RedoxSiteTransformerBase)|g

# Update registration method
s|MetalSiteTransformer\.register_transformer|# Registration handled by @register_redox_transformer decorator|g

# Update terminology in comments and docstrings (but preserve when specifically talking about metals)
s|MetalSite object|RedoxSite object|g
s|metal site transformer|redox site transformer|g
s|Metal Site|Redox Site|g

# Add import for CenterType
s|^import logging$|import logging\nfrom proprep.structure_prep.comprehensive_redox_detector import CenterType|
