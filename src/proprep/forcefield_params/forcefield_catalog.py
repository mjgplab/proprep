"""
Curated AMBER force field catalog — single source of truth shared by the
Topology Generator (system building) and the Force Field Explorer (browsing).

This is intentionally a *curated* subset of the mainstream force fields, with
human-readable names, descriptions, and recommended marks. It is NOT a scan of
AMBERHOME; the Force Field Explorer still scans AMBERHOME for full coverage and
uses this catalog only to ENRICH the scanned list (friendly names/descriptions,
recommended marks, and the protein-vs-modified-AA distinction).

Pure data + helpers; no imports from proprep packages, so any module can import
it without circular-import risk.
"""

from typing import Dict, Any


FORCEFIELD_OPTIONS = {'protein': {'title': 'PROTEIN FORCEFIELD',
             'options': [{'name': 'ff19SB',
                          'leaprc': 'leaprc.protein.ff19SB',
                          'recommended': True,
                          'description': 'Latest AMBER protein FF - use with OPC water'},
                         {'name': 'ff14SB',
                          'leaprc': 'leaprc.protein.ff14SB',
                          'description': 'Widely validated, good for most applications'},
                         {'name': 'ff15ipq',
                          'leaprc': 'leaprc.protein.ff15ipq',
                          'description': 'Implicitly polarized - REQUIRES SPC/Eb water'},
                         {'name': 'fb15',
                          'leaprc': 'leaprc.protein.fb15',
                          'description': 'Force Balance optimized - use with TIP3P-FB or TIP4P-FB'},
                         {'name': 'ff03.r1',
                          'leaprc': 'leaprc.protein.ff03.r1',
                          'description': 'Duan et al. (2003) - older alternative'},
                         {'name': 'ff99SBildn',
                          'leaprc': 'oldff/leaprc.ff99SBildn',
                          'description': 'Legacy, improved side-chain dihedrals'},
                         {'name': 'Constant pH (ff10)',
                          'leaprc': 'leaprc.constph',
                          'description': 'Titratable residues - loads ff10 internally'},
                         {'name': 'Constant Redox (ff10)',
                          'leaprc': 'leaprc.conste',
                          'description': 'Redox-active residues - loads ff10 internally'},
                         {'name': 'Constant pH + Redox (ff10)',
                          'leaprc': ['leaprc.constph', 'leaprc.conste'],
                          'description': 'Both titratable and redox-active residues - loads ff10 '
                                         'internally'}],
             'allow_none': True,
             'none_text': 'No protein in system'},
 'modified_aa': {'title': 'MODIFIED AMINO ACIDS',
                 'multi_select': True,
                 'description': 'Phosphorylated residues, spin labels, non-natural amino acids',
                 'options': [{'name': 'Modified AA (ff19SB)',
                              'leaprc': 'leaprc.protein.ff19SB_modAA',
                              'description': 'Selenomethionine, MTSL spin labels, etc.'},
                             {'name': 'Modified AA (ff14SB)',
                              'leaprc': 'leaprc.protein.ff14SB_modAA',
                              'description': 'Same modifications for ff14SB'},
                             {'name': 'Phosphorylated AA (ff19SB)',
                              'leaprc': 'leaprc.phosaa19SB',
                              'description': 'Phosphorylated Ser/Thr/Tyr/His for ff19SB'},
                             {'name': 'Phosphorylated AA (ff14SB)',
                              'leaprc': 'leaprc.phosaa14SB',
                              'description': 'Phosphorylated amino acids for ff14SB'},
                             {'name': 'Phosphorylated AA (ff99SB)',
                              'leaprc': 'leaprc.phosaa10',
                              'description': 'Phosphorylated amino acids for ff99SB and older'},
                             {'name': 'Mimetic residues (ff15ipq)',
                              'leaprc': 'leaprc.mimetic.ff15ipq',
                              'description': 'D-amino acids, beta-residues'},
                             {'name': 'Fluorinated AA (ff15ipq)',
                              'leaprc': 'leaprc.fluorine.ff15ipq',
                              'description': 'Fluorinated aromatic amino acids'}],
                 'allow_none': True,
                 'none_text': 'No modified amino acids'},
 'dna': {'title': 'DNA FORCEFIELD',
         'description': 'NOTE: For protein+DNA, only OL24 or OL21 are compatible with ff19SB!',
         'options': [{'name': 'OL24',
                      'leaprc': 'leaprc.DNA.OL24',
                      'recommended': True,
                      'description': 'Latest with sugar pucker refinements (2024)'},
                     {'name': 'OL21',
                      'leaprc': 'leaprc.DNA.OL21',
                      'description': 'Previous version (also compatible with ff19SB)'},
                     {'name': 'OL15',
                      'leaprc': 'leaprc.DNA.OL15',
                      'description': 'Older version (2015)'},
                     {'name': 'bsc1',
                      'leaprc': 'leaprc.DNA.bsc1',
                      'description': 'Barcelona forcefield'},
                     {'name': 'tumuc1',
                      'leaprc': 'leaprc.DNA.tumuc1',
                      'description': 'QM-refined electrostatic parameters'}],
         'allow_none': True,
         'none_text': 'No DNA in system'},
 'rna': {'title': 'RNA FORCEFIELD',
         'options': [{'name': 'OL3',
                      'leaprc': 'leaprc.RNA.OL3',
                      'recommended': True,
                      'description': 'Standard RNA forcefield'},
                     {'name': 'LJbb',
                      'leaprc': 'leaprc.RNA.LJbb',
                      'description': 'OL3 + improved phosphate LJ (better NMR agreement)'},
                     {'name': 'ROC',
                      'leaprc': 'leaprc.RNA.ROC',
                      'description': 'Rochester torsions'},
                     {'name': 'YIL',
                      'leaprc': 'leaprc.RNA.YIL',
                      'description': 'Yildirim chi modifications'},
                     {'name': 'Shaw',
                      'leaprc': 'leaprc.RNA.Shaw',
                      'description': 'DE Shaw modifications'},
                     {'name': 'modrna08',
                      'leaprc': 'leaprc.modrna08',
                      'description': 'Modified nucleosides (m6A, pseudouridine, etc.)'}],
         'allow_none': True,
         'none_text': 'No RNA in system'},
 'carbohydrates': {'title': 'CARBOHYDRATE FORCEFIELD',
                   'options': [{'name': 'GLYCAM_06j',
                                'leaprc': 'leaprc.GLYCAM_06j-1',
                                'recommended': True,
                                'description': 'GLYCAM force field for glycans'},
                               {'name': 'GLYCAM_06EPb',
                                'leaprc': 'leaprc.GLYCAM_06EPb',
                                'description': 'GLYCAM with extra points (for TIP5P water)'}],
                   'allow_none': True,
                   'none_text': 'No carbohydrates in system'},
 'lipids': {'title': 'LIPID FORCEFIELD',
            'options': [{'name': 'lipid21',
                         'leaprc': 'leaprc.lipid21',
                         'recommended': True,
                         'description': 'Latest lipid forcefield'},
                        {'name': 'lipid17',
                         'leaprc': 'leaprc.lipid17',
                         'description': 'Previous generation'}],
            'allow_none': True,
            'none_text': 'No lipids in system'},
 'small_molecules': {'title': 'SMALL MOLECULE FORCEFIELD',
                     'options': [{'name': 'GAFF2',
                                  'leaprc': 'leaprc.gaff2',
                                  'recommended': True,
                                  'description': 'General AMBER Force Field 2 (recommended)'},
                                 {'name': 'GAFF',
                                  'leaprc': 'leaprc.gaff',
                                  'description': 'Original GAFF (less accurate)'}],
                     'allow_none': True,
                     'none_text': 'No small molecules/ligands'},
 'water': {'title': 'WATER MODEL',
           'description': 'Must be compatible with protein forcefield!',
           'options': [{'name': 'OPC',
                        'leaprc': 'leaprc.water.opc',
                        'recommended': True,
                        'description': '4-point, optimal for ff19SB',
                        'box': 'OPCBOX'},
                       {'name': 'OPC3',
                        'leaprc': 'leaprc.water.opc3',
                        'description': '3-point OPC (faster, good accuracy)',
                        'box': 'OPC3BOX'},
                       {'name': 'OPC3-pol',
                        'leaprc': 'leaprc.water.opc3pol',
                        'description': 'Polarizable OPC3 (for polarization studies)',
                        'box': 'OPC3BOX'},
                       {'name': 'TIP4P-Ew',
                        'leaprc': 'leaprc.water.tip4pew',
                        'description': '4-point, good alternative',
                        'box': 'TIP4PEWBOX'},
                       {'name': 'TIP3P',
                        'leaprc': 'leaprc.water.tip3p',
                        'description': '3-point, traditional (not recommended)',
                        'box': 'TIP3PBOX'},
                       {'name': 'SPC/E',
                        'leaprc': 'leaprc.water.spce',
                        'description': '3-point model',
                        'box': 'SPCBOX'},
                       {'name': 'SPC/Eb',
                        'leaprc': 'leaprc.water.spceb',
                        'description': 'REQUIRED for ff15ipq',
                        'box': 'SPCBOX'},
                       {'name': 'TIP5P',
                        'leaprc': 'leaprc.water.tip5p',
                        'description': '5-point (use with GLYCAM_06EPb)',
                        'box': 'TIP5PBOX'},
                       {'name': 'TIP3P-FB',
                        'frcmod': 'frcmod.tip3pfb',
                        'description': 'For fb15 protein FF',
                        'box': 'TIP3PBOX'},
                       {'name': 'TIP4P-FB',
                        'frcmod': 'frcmod.tip4pfb',
                        'description': 'For fb15 protein FF',
                        'box': 'TIP4PBOX'}],
           'allow_none': True,
           'none_text': 'Implicit solvent (no water in system)'},
 'ions': {'title': 'DIVALENT+ ION PARAMETERS',
          'description': 'For Mg2+, Ca2+, Zn2+, Fe2+/3+, etc. Monovalent ions auto-loaded with '
                         'water.',
          'options': [{'name': '12-6-4 OPC (most accurate)',
                       'frcmod': 'frcmod.ionslm_1264_opc',
                       'recommended': True,
                       'description': 'Best for divalent+ (requires ParmEd add12_6_4)',
                       'for_water': ['OPC']},
                      {'name': '12-6 OPC',
                       'frcmod': 'frcmod.ionslm_126_opc',
                       'description': 'Standard Li/Merz for OPC',
                       'for_water': ['OPC']},
                      {'name': '12-6 OPC3',
                       'frcmod': 'frcmod.ionslm_126_opc3',
                       'description': 'Standard Li/Merz for OPC3',
                       'for_water': ['OPC3']},
                      {'name': 'IOD OPC (structural)',
                       'frcmod': 'frcmod.ionslm_iod_opc',
                       'description': 'Reproduces ion-oxygen distances',
                       'for_water': ['OPC']},
                      {'name': '12-6-4 TIP3P',
                       'frcmod': ['frcmod.ions1lm_1264_tip3p', 'frcmod.ions234lm_1264_tip3p'],
                       'description': 'Best for divalent+ with TIP3P',
                       'for_water': ['TIP3P']},
                      {'name': '12-6-4 SPC/E',
                       'frcmod': ['frcmod.ions1lm_1264_spce', 'frcmod.ions234lm_1264_spce'],
                       'description': 'Best for divalent+ with SPC/E',
                       'for_water': ['SPC/E']},
                      {'name': '12-6-4 TIP4P-Ew',
                       'frcmod': ['frcmod.ions1lm_1264_tip4pew', 'frcmod.ions234lm_1264_tip4pew'],
                       'description': 'Best for divalent+ with TIP4P-Ew',
                       'for_water': ['TIP4P-Ew']},
                      {'name': '12-6 TIP3P',
                       'frcmod': 'frcmod.ions234lm_126_tip3p',
                       'description': 'Standard Li/Merz for TIP3P',
                       'for_water': ['TIP3P']},
                      {'name': 'IOD TIP3P (structural)',
                       'frcmod': 'frcmod.ions234lm_iod_tip3p',
                       'description': 'Reproduces ion-oxygen distances',
                       'for_water': ['TIP3P']},
                      {'name': '12-6 SPC/E',
                       'frcmod': 'frcmod.ions234lm_126_spce',
                       'description': 'Standard Li/Merz for SPC/E',
                       'for_water': ['SPC/E']},
                      {'name': 'Joung-Cheatham TIP3P',
                       'frcmod': 'frcmod.ionsjc_tip3p',
                       'description': 'Alternative monovalent ions',
                       'for_water': ['TIP3P']},
                      {'name': 'Default only',
                       'frcmod': None,
                       'description': 'Use default ions from water model'}],
          'allow_none': True,
          'none_text': 'No divalent+ ions needed'}}


# Map curated catalog category keys -> the Force Field Explorer's display
# category labels. The Explorer scans AMBERHOME and groups files; this lets a
# cataloged leaprc carry the curated grouping (notably modified_aa, which keeps
# phosaa*/phosfb* OUT of "Protein" so they aren't mistaken for standalone FFs).
CATALOG_CATEGORY_TO_EXPLORER = {
    'protein': 'Protein',
    'modified_aa': 'Modified Amino Acids',
    'dna': 'DNA',
    'rna': 'RNA',
    'carbohydrates': 'Carbohydrate',
    'lipids': 'Lipid',
    'small_molecules': 'GAFF (General)',
    'water': 'Water',
    'ions': 'Ions',
}


# Recommended water model per protein force field, keyed by the protein FF's
# leaprc (the catalog's stable identifier) -> the water option ``name`` to mark
# as recommended in the WATER MODEL menu. Each pairing follows the water model
# the protein FF was parameterized/validated against:
#   - ff19SB was fit with OPC (its CMAP correction assumes OPC);
#   - ff14SB and the older/legacy protein FFs (ff03, ff99SBildn) use TIP3P;
#   - ff15ipq REQUIRES SPC/Eb (implicitly-polarized charges tuned to it);
#   - fb15 pairs with the Force-Balance water (TIP3P-FB);
#   - the constant-pH / constant-redox (ff10) sets are validated with TIP3P
#     (their reference titration/redox free energies were computed in TIP3P).
# A protein FF absent from this map has no specific pairing; callers should keep
# the catalog's own default recommendation (OPC) in that case.
PROTEIN_WATER_RECOMMENDATION = {
    'leaprc.protein.ff19SB': 'OPC',
    'leaprc.protein.ff14SB': 'TIP3P',
    'leaprc.protein.ff15ipq': 'SPC/Eb',
    'leaprc.protein.fb15': 'TIP3P-FB',
    'leaprc.protein.ff03.r1': 'TIP3P',
    'oldff/leaprc.ff99SBildn': 'TIP3P',
    'leaprc.constph': 'TIP3P',
    'leaprc.conste': 'TIP3P',
}


def recommended_water_for_protein(protein_leaprc) -> str:
    """Return the recommended water-model ``name`` for a selected protein FF.

    ``protein_leaprc`` may be a single leaprc string or a list of leaprc strings
    (the combined "Constant pH + Redox" option carries
    ``['leaprc.constph', 'leaprc.conste']``). For a list, the first constituent
    with a known pairing wins; both constant-pH and constant-redox map to TIP3P,
    so the combined set resolves to TIP3P either way.

    Returns ``None`` when no specific recommendation is known (e.g. no protein FF
    selected, or a protein FF not in :data:`PROTEIN_WATER_RECOMMENDATION`), so
    callers can fall back to the catalog's default recommendation.
    """
    if isinstance(protein_leaprc, (list, tuple)):
        for lr in protein_leaprc:
            water = PROTEIN_WATER_RECOMMENDATION.get(lr)
            if water:
                return water
        return None
    return PROTEIN_WATER_RECOMMENDATION.get(protein_leaprc)


def recommended_ions_for_water(water_name) -> str:
    """Return the recommended divalent+ ion-parameter option ``name`` for a
    selected water model.

    Divalent+ ion (Li/Merz) parameters are calibrated per water model, so the
    recommended ion set must follow the chosen water — not the catalog's static
    default. This scans the catalog's ion options for those whose ``for_water``
    list includes ``water_name``, preferring the most accurate ``12-6-4``
    variant, then any matching set.

    Falls back to ``'Default only'`` when the water model has no dedicated
    Li/Merz set (e.g. polarizable, 5-point, or Force-Balance waters) or when no
    water model was selected (implicit solvent). ``'Default only'`` is always a
    real option in the catalog, so callers can mark it recommended directly.
    """
    if not water_name:
        return 'Default only'
    matches = [
        opt for opt in FORCEFIELD_OPTIONS['ions']['options']
        if water_name in opt.get('for_water', [])
    ]
    if not matches:
        return 'Default only'
    for opt in matches:
        if opt['name'].startswith('12-6-4'):
            return opt['name']
    return matches[0]['name']


def build_leaprc_index() -> Dict[str, Dict[str, Any]]:
    """Flatten the catalog into a ``leaprc_name -> metadata`` lookup.

    Returns a dict keyed by leaprc filename (e.g. ``leaprc.protein.ff14SB``)
    with curated ``name``, ``description``, ``recommended``, ``category`` and an
    ``is_addon`` flag (True for modified-AA add-ons that require a base protein
    FF). Options keyed only by ``frcmod`` (no leaprc) are skipped; list-valued
    ``leaprc`` entries register every member. The first occurrence of a leaprc
    wins, so a file's primary category is stable.
    """
    index: Dict[str, Dict[str, Any]] = {}
    for cat_key, cat in FORCEFIELD_OPTIONS.items():
        explorer_cat = CATALOG_CATEGORY_TO_EXPLORER.get(cat_key, 'Other')
        is_addon = (cat_key == 'modified_aa')
        for opt in cat.get('options', []):
            leaprc = opt.get('leaprc')
            if not leaprc:
                continue
            names = leaprc if isinstance(leaprc, list) else [leaprc]
            for lr in names:
                index.setdefault(lr, {
                    'name': opt.get('name', lr),
                    'description': opt.get('description', ''),
                    'recommended': bool(opt.get('recommended', False)),
                    'is_addon': is_addon,
                    'category': explorer_cat,
                })
    return index
