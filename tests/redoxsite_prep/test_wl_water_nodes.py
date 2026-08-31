"""Waters as WL fingerprint nodes: the hydrated fingerprint distinguishes
coordinated waters (so MW1/MW2-style renames auto-resolve), while the
anhydrous fingerprint and everything baked before it behave exactly as
before."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from test_auto_rename_transformer import _atom, _bond, _make_site  # noqa: E402

from proprep.redoxsite_prep.transformation.auto_rename import (
    AutoRenameTransformerBase, connectivity_signature)


M1, M2 = (0.0, 0.0, 0.0), (8.0, 0.0, 0.0)
W1_O, W2_O = (0.0, 2.0, 0.0), (8.0, 2.0, 0.0)
HIS_NE2 = (0.0, -2.0, 0.0)          # only M1 has a His: breaks the metal symmetry


def _site(with_water_bonds=True):
    atoms = [
        _atom("A", 1, "MN", "MN", M1, "MN"),
        _atom("A", 2, "MN", "MN", M2, "MN"),
        _atom("A", 3, "HIS", "NE2", HIS_NE2, "N"),
        _atom("A", 4, "HOH", "O", W1_O, "O"),
        _atom("A", 5, "HOH", "O", W2_O, "O"),
    ]
    bonds = [_bond(M1, HIS_NE2, "MN", "N")]
    if with_water_bonds:
        bonds += [_bond(M1, W1_O, "MN", "O"), _bond(M2, W2_O, "MN", "O")]
    return _make_site(atoms, bonds)


def test_hydrated_fingerprint_includes_and_distinguishes_waters():
    site = _site()
    dry = connectivity_signature(site)
    wet = connectivity_signature(site, include_waters=True)
    assert ("A", 4) not in dry and ("A", 5) not in dry
    assert wet[("A", 4)].startswith("L:WAT#") and wet[("A", 5)].startswith("L:WAT#")
    assert wet[("A", 4)] != wet[("A", 5)]          # bonded to distinguishable metals
    assert wet[("A", 1)] != wet[("A", 2)]


def test_anhydrous_fingerprint_is_untouched_by_water_bonds():
    # the invariance that motivated the original exclusion still holds by default
    assert connectivity_signature(_site(True)) == connectivity_signature(_site(False))


def _transformer(table):
    class T(AutoRenameTransformerBase):
        NAME = "t"; DESCRIPTION = "t"; RENAMES = table
    return T


def test_water_roles_auto_resolve_via_hydrated_signatures():
    site = _site()
    wet = connectivity_signature(site, include_waters=True)
    dry = connectivity_signature(site)
    table = [
        {"resname": "MN", "target": "MN1", "signature": dry[("A", 1)], "signature_hydrated": wet[("A", 1)]},
        {"resname": "MN", "target": "MN2", "signature": dry[("A", 2)], "signature_hydrated": wet[("A", 2)]},
        {"resname": "HOH", "target": "MW1", "signature": None, "signature_hydrated": wet[("A", 4)]},
        {"resname": "HOH", "target": "MW2", "signature": None, "signature_hydrated": wet[("A", 5)]},
    ]
    matched, missing = _transformer(table).match_components(site)
    assert not missing and "_ambiguous" not in matched
    assert (matched["MW1_chain"], matched["MW1_id"]) == ("A", 4)
    assert (matched["MW2_chain"], matched["MW2_id"]) == ("A", 5)
    assert (matched["MN1_chain"], matched["MN1_id"]) == ("A", 1)


def test_missing_water_bonds_fall_back_to_anhydrous_for_metals():
    # emitted with water bonds; reused on a site where the user left them out
    baked = _site(True)
    wet = connectivity_signature(baked, include_waters=True)
    dry = connectivity_signature(baked)
    table = [
        {"resname": "MN", "target": "MN1", "signature": dry[("A", 1)], "signature_hydrated": wet[("A", 1)]},
        {"resname": "MN", "target": "MN2", "signature": dry[("A", 2)], "signature_hydrated": wet[("A", 2)]},
        {"resname": "HOH", "target": "MW1", "signature": None, "signature_hydrated": wet[("A", 4)]},
        {"resname": "HOH", "target": "MW2", "signature": None, "signature_hydrated": wet[("A", 5)]},
    ]
    reuse = _site(with_water_bonds=False)
    matched, missing = _transformer(table).match_components(reuse)
    # metals still resolve (anhydrous fallback); the waters go to the user
    assert (matched["MN1_chain"], matched["MN1_id"]) == ("A", 1)
    assert (matched["MN2_chain"], matched["MN2_id"]) == ("A", 2)
    ambiguous = matched.get("_ambiguous") or []
    assert any("HOH" in b["label"] for b in ambiguous)


def test_old_transformers_without_hydrated_signatures_match_as_before():
    site = _site()
    dry = connectivity_signature(site)
    table = [
        {"resname": "MN", "target": "MN1", "signature": dry[("A", 1)]},
        {"resname": "MN", "target": "MN2", "signature": dry[("A", 2)]},
    ]
    matched, missing = _transformer(table).match_components(site)
    assert not missing and "_ambiguous" not in matched
    assert (matched["MN1_chain"], matched["MN1_id"]) == ("A", 1)
