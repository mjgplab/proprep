#!/usr/bin/env python3
"""Flavin FAD family transformer (FAD, FAS, FAR, FAH, FAQ).

Detects FAD as an ORGANIC_COFACTOR and renames the residue to match the
chosen redox/protonation state. All 5 states share the same heavy-atom
name set; only HN1 and HN5 protons toggle.

FAD differs from FMN by carrying an adenine arm via pyrophosphate. The
upstream .lib uses C4X/C5X for the isoalloxazine ring junction carbons to
disambiguate from the adenine ring's C4A/C5A. Older PDBs may use C4A/C5A
on the isoalloxazine; ATOM_NAME_RENAMES handles that fallback case.

FAR (anionic red semiquinone) Set 2 falls back to Phase-0a charges because
Phase 0c optimization did not converge; documented in the metadata.
"""

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    register_redox_transformer,
)
from proprep.redoxsite_prep.transformation.transformers._cofactor_base import (
    OrganicCofactorTransformerBase,
)


@register_redox_transformer
class FlavinFADTransformer(OrganicCofactorTransformerBase):
    TRANSFORMER_NAME = "flavin_fad"
    DESCRIPTION = "Flavin adenine dinucleotide (FAD) and its 5 redox/protonation variants"
    FORCEFIELD_PATH = "flavin/fad"

    CANONICAL_INPUT_RESNAMES = {"FAD"}

    REDOX_STATE_OPTIONS = [
        ("oxidized", "FAD", "Quinone (oxidized; PDB-standard FAD)", True),
        ("semiquinone_neutral", "FAS",
         "Neutral blue semiquinone FADH• (N5 protonated, N1 deprotonated)", False),
        ("semiquinone_anionic", "FAR",
         "Anionic red semiquinone FAD•⁻ (N1 + N5 both deprotonated). "
         "Set 2 falls back to Phase-0a charges (Phase 0c did not converge).",
         False),
        ("hydroquinone_anionic", "FAH",
         "Anionic hydroquinone FADH⁻ (N5 protonated, N1 deprotonated)", False),
        ("hydroquinone_neutral", "FAQ",
         "Neutral hydroquinone FADH₂ (N1 + N5 both protonated)", False),
    ]

    # PDB-CCD canonical FAD already uses C4X/C5X for the isoalloxazine
    # ring-junction carbons (C4A/C5A are the ADENINE ring's, present
    # simultaneously in the same residue). A blanket rename here would
    # corrupt the adenine ring by colliding atom names. If a user's PDB has
    # legacy non-CCD naming on isoalloxazine (C4A/C5A on the iso side),
    # they must normalize externally; the detector should surface the
    # conflict via evaluate_redox_site. Empty by default.
    ATOM_NAME_RENAMES = {}
