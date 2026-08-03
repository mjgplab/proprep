#!/usr/bin/env python3
"""Flavin FMN family transformer (FMN, FMS, FMR, FMH, FMQ).

Detects FMN as an ORGANIC_COFACTOR and renames the residue to match the
chosen redox/protonation state. All 5 states share the same heavy-atom
name set; only HN1 and HN5 protons toggle (handled by the Hydrogen Addition
step using the new residue's .lib as template).
"""

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    register_redox_transformer,
)
from proprep.redoxsite_prep.transformation.transformers._cofactor_base import (
    OrganicCofactorTransformerBase,
)


@register_redox_transformer
class FlavinFMNTransformer(OrganicCofactorTransformerBase):
    TRANSFORMER_NAME = "flavin_fmn"
    DESCRIPTION = "Flavin mononucleotide (FMN) and its 5 redox/protonation variants"
    FORCEFIELD_PATH = "flavin/fmn"

    CANONICAL_INPUT_RESNAMES = {"FMN"}

    REDOX_STATE_OPTIONS = [
        ("oxidized", "FMN", "Quinone (oxidized; PDB-standard FMN)", True),
        ("semiquinone_neutral", "FMS",
         "Neutral blue semiquinone FMNH• (N5 protonated, N1 deprotonated)", False),
        ("semiquinone_anionic", "FMR",
         "Anionic red semiquinone FMN•⁻ (N1 + N5 both deprotonated)", False),
        ("hydroquinone_anionic", "FMH",
         "Anionic hydroquinone FMNH⁻ (N5 protonated, N1 deprotonated)", False),
        ("hydroquinone_neutral", "FMQ",
         "Neutral hydroquinone FMNH₂ (N1 + N5 both protonated)", False),
    ]

    # Atom names in the upstream .lib are PDB-CCD canonical for FMN; no
    # renames needed today. Reserved for future PDB-CCD revisions or for
    # Reduce 4.10-specific quirks if encountered.
    ATOM_NAME_RENAMES = {}
