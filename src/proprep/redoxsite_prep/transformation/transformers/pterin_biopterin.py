#!/usr/bin/env python3
"""Pterin biopterin family transformer (H4B, H4C, H3R, H2Q, H2B, BIO).

Detects H4B as an ORGANIC_COFACTOR and renames to one of six redox/protonation
states spanning the tetrahydro (resting), 1e-oxidized (cation radical and
N5-deprotonated radical), 2e-oxidized (quinonoid and 7,8-dihydro), and
4e-oxidized (fully aromatic biopterin) forms. State protonation differs by
HN5 and HN8; handled by the Hydrogen Addition step using the new lib as
template.
"""

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    register_redox_transformer,
)
from proprep.redoxsite_prep.transformation.transformers._cofactor_base import (
    OrganicCofactorTransformerBase,
)


@register_redox_transformer
class PterinBiopterinTransformer(OrganicCofactorTransformerBase):
    TRANSFORMER_NAME = "pterin_biopterin"
    DESCRIPTION = "Biopterin family (H4B, H4C, H3R, H2Q, H2B, BIO) — 6 redox/protonation variants"
    FORCEFIELD_PATH = "pterin/biopterin"

    # Detector matches input resname against PDB-CCD canonical pterin codes
    # only — H4B (tetrahydrobiopterin), H2B (7,8-dihydrobiopterin), BIO
    # (fully aromatic biopterin). The other four output codes (H4C/H3R/H2Q
    # for cation-radical / N5-radical / quinonoid-dihydro) are FF-internal
    # labels and are NOT included: those PDB-CCD codes refer to unrelated
    # chemistry (a diazepinium, a methylphenyl-pyrazole, and an
    # oxobutanamide respectively), so accepting them as input would surface
    # false-positive transformer compatibility for any structure that
    # happens to contain those unrelated ligands. ProPrep is run once per
    # structure — reloading a previously-transformed PDB to re-detect the
    # redox site is not a supported workflow.
    CANONICAL_INPUT_RESNAMES = {"H4B", "H2B", "BIO"}

    REDOX_STATE_OPTIONS = [
        ("tetrahydro", "H4B",
         "5,6,7,8-tetrahydrobiopterin (resting reduced; PDB-standard H4B)", True),
        ("cation_radical", "H4C",
         "H4B cation radical H4B•⁺ (1e-oxidized; on-cycle NOS electron donor)", False),
        ("n5_radical", "H3R",
         "Trihydrobiopterin neutral N5 radical H3B• (1e-oxidized; N5-deprotonated)",
         False),
        ("quinonoid_dihydro", "H2Q",
         "Quinonoid 6,7-dihydrobiopterin q-H2B (2e-oxidized at C6=C7)", False),
        ("78_dihydro", "H2B",
         "7,8-dihydrobiopterin (2e-oxidized at C7=C8a; q-H2B tautomer)", False),
        ("fully_oxidized", "BIO",
         "Biopterin (fully aromatic pterin; 4e-oxidized; PDB-standard BIO)", False),
    ]

    # No atom-name renames — H4B family follows PDB CCD canonical naming.
    ATOM_NAME_RENAMES = {}
