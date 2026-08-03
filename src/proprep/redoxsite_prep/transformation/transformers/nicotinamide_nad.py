#!/usr/bin/env python3
"""Nicotinamide NAD family transformer (NAD, NDH).

Detects NAD (oxidized) or NAI (PDB-CCD canonical NADH) as an ORGANIC_COFACTOR
and renames to the chosen redox state. NDH is the library's reduced-state code
(distinct from PDB-CCD's NAI); the transformer handles the NAI → NDH rename.
"""

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    register_redox_transformer,
)
from proprep.redoxsite_prep.transformation.transformers._cofactor_base import (
    OrganicCofactorTransformerBase,
)


@register_redox_transformer
class NicotinamideNADTransformer(OrganicCofactorTransformerBase):
    TRANSFORMER_NAME = "nicotinamide_nad"
    DESCRIPTION = "NAD⁺/NADH (NAD, NDH) — nicotinamide adenine dinucleotide (no 2'-phosphate)"
    FORCEFIELD_PATH = "nicotinamide/nad"

    # Accept NAD (PDB-CCD canonical NAD⁺) and NAI (PDB-CCD canonical NADH)
    # as input. The FF-internal code NDH (for the reduced state's output) is
    # NOT included here: PDB-CCD's NDH refers to an unrelated chemistry
    # (1,2-dihydroxy-1,2-dihydronaphthalene), so accepting it as input would
    # surface false-positive compatibility for any structure containing that
    # ligand. ProPrep is run once per structure; reloading a previously-
    # transformed PDB to re-detect the redox site is not a supported
    # workflow.
    CANONICAL_INPUT_RESNAMES = {"NAD", "NAI"}

    REDOX_STATE_OPTIONS = [
        ("oxidized", "NAD",
         "NAD⁺ — aromatic nicotinamide cation, diphosphate dianion (net −1)",
         True),
        ("reduced", "NDH",
         "NADH — 1,4-dihydronicotinamide (4-pro-R) (net −2). Library code is NDH; "
         "PDB-CCD canonical NADH is NAI. Transformer renames as needed.",
         False),
    ]

    # No atom-name renames — the NAD and NDH (or NAI input) atom sets follow
    # PDB CCD canonical naming.
    ATOM_NAME_RENAMES = {}
