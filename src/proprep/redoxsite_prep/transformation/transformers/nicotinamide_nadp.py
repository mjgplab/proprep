#!/usr/bin/env python3
"""Nicotinamide NADP family transformer (NAP, NDP).

Detects NAP or NDP as an ORGANIC_COFACTOR and renames the residue to match
the chosen redox/protonation state. NAP is oxidized (aromatic nicotinamide
cation); NDP is reduced (1,4-dihydronicotinamide, 4-pro-R hydride).

NP2 (NADPH with mono-protonated 2'-phosphate) was deferred from the v1.0
fragment-typed library — the upstream parameterization seed was missing
HN4 and used an inconsistent charge target. NP2 will rejoin in a follow-up
release; for now use NDP (dianion 2'-phosphate) for physiological pH.
"""

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    register_redox_transformer,
)
from proprep.redoxsite_prep.transformation.transformers._cofactor_base import (
    OrganicCofactorTransformerBase,
)


@register_redox_transformer
class NicotinamideNADPTransformer(OrganicCofactorTransformerBase):
    TRANSFORMER_NAME = "nicotinamide_nadp"
    DESCRIPTION = "NADP⁺/NADPH (NAP, NDP) — nicotinamide adenine dinucleotide phosphate"
    FORCEFIELD_PATH = "nicotinamide/nadp"

    # Accept either NAP (oxidized canonical) or NDP (reduced canonical) as
    # input. Detector reports whichever appears in the PDB; user picks the
    # target state.
    CANONICAL_INPUT_RESNAMES = {"NAP", "NDP"}

    REDOX_STATE_OPTIONS = [
        ("oxidized", "NAP",
         "NADP⁺ — aromatic nicotinamide cation, 2'-phosphate dianion (net −3)",
         True),
        ("reduced", "NDP",
         "NADPH — 1,4-dihydronicotinamide (4-pro-R), 2'-phosphate dianion (net −4)",
         False),
    ]

    # No atom-name renames required — upstream library follows PDB CCD
    # canonical naming for both NAP and NDP atom sets.
    ATOM_NAME_RENAMES = {}
