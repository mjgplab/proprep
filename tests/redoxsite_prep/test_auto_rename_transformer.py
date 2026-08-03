"""Tests for auto-emitted residue-rename transformers.

Covers the tiered matcher (unique name, same-name/same-target, same-name/
different-target with signature auto-resolve vs. ambiguous fallback), the
transformation sequence it produces, and the emit -> load -> register
round-trip through ~/.proprep-style user dirs.
"""

from types import SimpleNamespace

import pytest

from proprep.redoxsite_prep.transformation.auto_rename import (
    AutoRenameTransformerBase,
    connectivity_signature,
    coordination_signature,
    emit_rename_transformer,
    load_user_transformers,
)


# ---------------------------------------------------------------------------
# Lightweight RedoxSite stand-in
# ---------------------------------------------------------------------------

def _atom(chain, resid, resname, atom_name, coords, element):
    return SimpleNamespace(chain=chain, resid=resid, resname=resname,
                           atom_name=atom_name, coords=coords, element=element)


def _bond(c1, c2, e1, e2, chemical_type="coordinate"):
    return SimpleNamespace(atom1_coords=c1, atom2_coords=c2,
                           atom1_element=e1, atom2_element=e2,
                           chemical_type=chemical_type)


def _make_site(atoms, bonds):
    coord_to_pdb = {
        a.coords: {"chain": a.chain, "resid": a.resid, "resname": a.resname,
                   "atom_name": a.atom_name, "element": a.element}
        for a in atoms
    }
    return SimpleNamespace(site_id="site_1", atoms=atoms, bonds=bonds,
                           coord_to_pdb=coord_to_pdb, residue_groups={})


def _wl_table(site, mapping):
    """Build a RENAMES table with WL connectivity signatures, mirroring the
    MCPB-4 hook. ``mapping`` is {(chain, resid): target}."""
    labels = connectivity_signature(site)
    resname_by_key = {(a.chain, a.resid): a.resname for a in site.atoms}
    table = []
    for (chain, resid), target in mapping.items():
        entry = {"resname": resname_by_key[(chain, resid)], "target": target}
        if (chain, resid) in labels:
            entry["signature"] = labels[(chain, resid)]
        table.append(entry)
    return table


# Coordinates
M1 = (0.0, 0.0, 0.0)
M2 = (4.0, 0.0, 0.0)
HID_NE2 = (0.0, 2.0, 0.0)
ASP_OD1 = (-1.0, 0.0, 0.0)
GLU65_OE1 = (1.0, 1.0, 0.0)   # bridging Glu: bonds to both metals
GLU65_OE2 = (3.0, 1.0, 0.0)
GLU93_OE1 = (1.0, -1.0, 0.0)  # terminal Glu: bonds to one metal


def _dimn_atoms():
    return [
        _atom("A", 500, "MN", "MN", M1, "MN"),
        _atom("A", 501, "MN", "MN", M2, "MN"),
        _atom("A", 45, "HID", "NE2", HID_NE2, "N"),
        _atom("A", 108, "ASP", "OD1", ASP_OD1, "O"),
        _atom("A", 65, "GLU", "OE1", GLU65_OE1, "O"),
        _atom("A", 65, "GLU", "OE2", GLU65_OE2, "O"),
        _atom("A", 93, "GLU", "OE1", GLU93_OE1, "O"),
    ]


def _dimn_bonds():
    return [
        _bond(HID_NE2, M1, "N", "MN"),
        _bond(ASP_OD1, M1, "O", "MN"),
        _bond(GLU65_OE1, M1, "O", "MN"),   # Glu65 bridges M1 and M2
        _bond(GLU65_OE2, M2, "O", "MN"),
        _bond(GLU93_OE1, M1, "O", "MN"),   # Glu93 terminal on M1
        _bond(M1, M2, "MN", "MN", chemical_type="metal-metal"),
    ]


# ---------------------------------------------------------------------------
# coordination_signature
# ---------------------------------------------------------------------------

def test_signature_distinguishes_bridging_from_terminal():
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    bridging = coordination_signature(site, "A", 65)
    terminal = coordination_signature(site, "A", 93)
    assert bridging == {"n_metals": 2, "coord_atoms": ["OE1", "OE2"]}
    assert terminal == {"n_metals": 1, "coord_atoms": ["OE1"]}


# ---------------------------------------------------------------------------
# Tier 1 — unique residue names
# ---------------------------------------------------------------------------

def _mk(renames, **attrs):
    """Build a throwaway AutoRenameTransformerBase subclass with a table."""
    return type("T", (AutoRenameTransformerBase,), {"RENAMES": renames,
                                                     "DESCRIPTION": "t", **attrs})


def test_tier1_unique_names_resolve_and_rename():
    site = _make_site(
        [_atom("A", 45, "HID", "NE2", HID_NE2, "N"),
         _atom("A", 108, "ASP", "OD1", ASP_OD1, "O")],
        [],
    )
    T = _mk([{"resname": "HID", "target": "HD1"},
             {"resname": "ASP", "target": "AS1"}])
    matched, missing = T.match_components(site)
    assert missing == []
    assert "_ambiguous" not in matched
    seq = T.get_transformation_sequence(matched, {})
    got = {(s["selector"]["residue_id"], s["action"]["change_residue_name"]) for s in seq}
    assert got == {(45, "HD1"), (108, "AS1")}
    # selectors key on chain+resid only (protonation-agnostic); no residue_name
    assert all("residue_name" not in s["selector"] for s in seq)


def test_wl_signature_invariant_to_protonation_and_water_names():
    # Same site, named as at MCPB time (HID, WAT) vs freshly loaded (HIS, HOH).
    # The connectivity fingerprints must be IDENTICAL, else a reuse run can't
    # match the baked signatures and needlessly falls to the prompt.
    def site(his, wat):
        atoms = [
            _atom("A", 1, "MN", "MN", M1, "MN"),
            _atom("A", 2, his, "NE2", HID_NE2, "N"),
            _atom("A", 3, wat, "O", GLU65_OE1, "O"),
        ]
        bonds = [_bond(HID_NE2, M1, "N", "MN"), _bond(GLU65_OE1, M1, "O", "MN")]
        return _make_site(atoms, bonds)

    mcpb = connectivity_signature(site("HID", "WAT"))
    fresh = connectivity_signature(site("HIS", "HOH"))
    # Compare by (chain, resid) — the residues are positionally identical.
    assert mcpb == fresh


def test_wl_signature_invariant_to_water_and_metal_metal_bonds():
    # MCPB-time graph coordinates the metals with restrained waters and a Mn-Mn
    # bond; the reuse site omits both (waters -> MD restraints, Mn-Mn optional).
    # The fingerprint must ignore those edges so the two graphs still match,
    # otherwise every label shifts (WL propagates) and reuse falls to the prompt.
    WAT_O = (2.0, 2.0, 0.0)

    def site(with_incidental):
        atoms = [
            _atom("A", 500, "MN", "MN", M1, "MN"),
            _atom("A", 501, "MN", "MN", M2, "MN"),
            _atom("A", 45, "HIS", "NE2", HID_NE2, "N"),
            _atom("A", 65, "GLU", "OE1", GLU65_OE1, "O"),  # bridges both metals
            _atom("A", 65, "GLU", "OE2", GLU65_OE2, "O"),
            _atom("A", 93, "GLU", "OE1", GLU93_OE1, "O"),  # terminal on M1
        ]
        bonds = [
            _bond(HID_NE2, M1, "N", "MN"),
            _bond(GLU65_OE1, M1, "O", "MN"),
            _bond(GLU65_OE2, M2, "O", "MN"),
            _bond(GLU93_OE1, M1, "O", "MN"),
        ]
        if with_incidental:
            atoms.append(_atom("A", 700, "HOH", "O", WAT_O, "O"))
            bonds.append(_bond(WAT_O, M2, "O", "MN"))                    # restrained water
            bonds.append(_bond(M1, M2, "MN", "MN", chemical_type="metal-metal"))
        return _make_site(atoms, bonds)

    param = connectivity_signature(site(True))    # water + Mn-Mn present
    reuse = connectivity_signature(site(False))   # both omitted
    # The water residue is absent from the reuse signature; every shared residue
    # must carry an identical label.
    shared = set(param) & set(reuse)
    assert ("A", 700) not in reuse                # water excluded from the graph
    assert {k: param[k] for k in shared} == {k: reuse[k] for k in shared}
    # And the two Glu remain distinguishable without the incidental edges.
    assert reuse[("A", 65)] != reuse[("A", 93)]


def test_wl_signature_invariant_to_carboxylate_oxygen_naming():
    # A Glu/Asp carboxylate's two oxygens are chemically equivalent; which one is
    # nearest the metal (OE1 vs OE2 / OD1 vs OD2) flips between structures. The
    # fingerprint must key on the donor ELEMENT, not the atom name, or the two
    # graphs diverge and every reachable label shifts (real 7K0W di-Mn bug).
    def site(glu_oxygen):
        atoms = [
            _atom("A", 500, "MN", "MN", M1, "MN"),
            _atom("A", 45, "HIS", "NE2", HID_NE2, "N"),
            _atom("A", 65, "GLU", glu_oxygen, GLU65_OE1, "O"),
        ]
        bonds = [_bond(HID_NE2, M1, "N", "MN"),
                 _bond(GLU65_OE1, M1, "O", "MN")]
        return _make_site(atoms, bonds)

    assert connectivity_signature(site("OE1")) == connectivity_signature(site("OE2"))


def test_ambiguous_block_carries_debug_diagnostic():
    # When a Tier-3 group can't auto-resolve because the baked signatures don't
    # match the site, the ambiguous block must include a human-readable
    # diagnostic showing expected vs computed fingerprints.
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    # Bake signatures that deliberately do NOT occur in this site.
    T = _mk([{"resname": "GLU", "target": "GL1", "signature": "L:GLU#deadbeefdeadbeef"},
             {"resname": "GLU", "target": "GL2", "signature": "L:GLU#feedfacefeedface"}])
    matched, _ = T.match_components(site)
    blocks = matched.get("_ambiguous")
    assert blocks and "diagnostic" in blocks[0]
    diag = blocks[0]["diagnostic"]
    assert "GL1" in diag and "deadbeef" in diag       # expected role + baked sig
    assert "A:65" in diag or "A:93" in diag           # candidate residues listed
    assert "no baked match" in diag                   # explains the failure


def test_protonation_variant_matches_standard_residue():
    # Transformer baked "HID" (post-protonation MCPB name); the fresh structure
    # still has "HIS" (protonation analyzer runs later). Must still match + rename.
    site = _make_site([_atom("A", 45, "HIS", "NE2", HID_NE2, "N")], [])
    T = _mk([{"resname": "HID", "target": "HD1"}])
    assert T.evaluate_redox_site(site).is_valid is True
    matched, missing = T.match_components(site)
    assert missing == []
    seq = T.get_transformation_sequence(matched, {})
    assert seq[0]["selector"] == {"chain_id": "A", "residue_id": 45}
    assert seq[0]["action"]["change_residue_name"] == "HD1"


def test_atom_renames_emitted_in_action():
    # A small-molecule ligand entry carries an atom_renames map (PDB -> mol2/lib
    # names). The transformation step must rename the residue AND its atoms in one
    # action so a reuse structure (PDB atom names) matches the antechamber-named lib.
    site = _make_site([_atom("A", 182, "E4Z", "O1", (0.0, 0.0, 0.0), "O")], [])
    T = _mk([{"resname": "E4Z", "target": "EZ1",
              "atom_renames": {"C4": "C3", "S3": "S1"}}])
    seq = T.get_transformation_sequence({"EZ1_chain": "A", "EZ1_id": 182}, {})
    assert len(seq) == 1
    action = seq[0]["action"]
    assert action["change_residue_name"] == "EZ1"
    assert action["rename_atoms"] == {"C4": "C3", "S3": "S1"}
    # A residue without atom_renames must NOT get a rename_atoms key.
    T2 = _mk([{"resname": "HID", "target": "HD1"}])
    seq2 = T2.get_transformation_sequence({"HD1_chain": "A", "HD1_id": 45}, {})
    assert "rename_atoms" not in seq2[0]["action"]


def test_missing_residue_reported():
    site = _make_site([_atom("A", 45, "HID", "NE2", HID_NE2, "N")], [])
    T = _mk([{"resname": "HID", "target": "HD1"},
             {"resname": "ASP", "target": "AS1"}])
    matched, missing = T.match_components(site)
    assert "ASP" in missing


# ---------------------------------------------------------------------------
# Tier 2 — several residues of one name, same target
# ---------------------------------------------------------------------------

def test_tier2_same_target_renames_all():
    atoms = [_atom("A", n, "CYS", "SG", (float(n), 0.0, 0.0), "S") for n in (10, 20, 30, 40)]
    site = _make_site(atoms, [])
    T = _mk([{"resname": "CYS", "target": "CYZ"} for _ in range(4)])
    matched, missing = T.match_components(site)
    assert missing == []
    assert "_ambiguous" not in matched
    seq = T.get_transformation_sequence(matched, {})
    assert len(seq) == 4
    assert {s["action"]["change_residue_name"] for s in seq} == {"CYZ"}
    assert {s["selector"]["residue_id"] for s in seq} == {10, 20, 30, 40}


# ---------------------------------------------------------------------------
# Tier 3 — same name, different targets
# ---------------------------------------------------------------------------

def test_tier3_signature_autoresolves():
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    # Two Glu -> GL1 (bridging) / GL2 (terminal); distinct WL neighbourhoods.
    T = _mk(_wl_table(site, {("A", 65): "GL1", ("A", 93): "GL2"}))
    matched, missing = T.match_components(_make_site(_dimn_atoms(), _dimn_bonds()))
    assert missing == []
    assert "_ambiguous" not in matched
    # Bridging Glu 65 -> GL1, terminal Glu 93 -> GL2, resolved automatically.
    assert matched["GL1_id"] == 65
    assert matched["GL2_id"] == 93
    seq = T.get_transformation_sequence(matched, {})
    got = {(s["selector"]["residue_id"], s["action"]["change_residue_name"]) for s in seq}
    assert got == {(65, "GL1"), (93, "GL2")}


def test_two_metals_resolved_by_ligand_set():
    # The WL win: two Mn with DIFFERENT ligand sets are told apart by their
    # coordination neighbourhood — no prompt. (Mn 500 has 4 ligands, Mn 501 has
    # one bridging Glu carboxylate + the Mn-Mn bond.)
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    T = _mk(_wl_table(site, {("A", 500): "MN1", ("A", 501): "MN2"}))
    matched, missing = T.match_components(_make_site(_dimn_atoms(), _dimn_bonds()))
    assert missing == []
    assert "_ambiguous" not in matched
    assert matched["MN1_id"] == 500
    assert matched["MN2_id"] == 501


def test_tier3_symmetric_is_ambiguous():
    # Both Glu terminal on the same metal -> identical WL neighbourhoods -> can't
    # auto-resolve -> declare an _ambiguous block for the manager to prompt.
    atoms = [
        _atom("A", 500, "MN", "MN", M1, "MN"),
        _atom("A", 65, "GLU", "OE1", GLU65_OE1, "O"),
        _atom("A", 93, "GLU", "OE1", GLU93_OE1, "O"),
    ]
    bonds = [_bond(GLU65_OE1, M1, "O", "MN"), _bond(GLU93_OE1, M1, "O", "MN")]
    site = _make_site(atoms, bonds)
    T = _mk(_wl_table(site, {("A", 65): "GL1", ("A", 93): "GL2"}))
    matched, missing = T.match_components(site)
    assert "_ambiguous" in matched
    block = matched["_ambiguous"][0]
    assert block["roles"] == {"GL1": 1, "GL2": 1}
    assert len(block["candidates"]) == 2
    assert {c["resid"] for c in block["candidates"]} == {65, 93}
    # total slots must equal candidate count (resolver contract)
    assert sum(block["roles"].values()) == len(block["candidates"])


def test_tier3_ambiguous_then_resolved_produces_sequence():
    # Simulate the manager writing back resolved <role>_id/<role>_chain.
    T = _mk([
        {"resname": "GLU", "target": "GL1", "signature": {"n_metals": 1, "coord_atoms": ["OE1"]}},
        {"resname": "GLU", "target": "GL2", "signature": {"n_metals": 1, "coord_atoms": ["OE1"]}},
    ])
    resolved = {"GL1_chain": "A", "GL1_id": 93, "GL2_chain": "A", "GL2_id": 65}
    seq = T.get_transformation_sequence(resolved, {})
    got = {(s["selector"]["residue_id"], s["action"]["change_residue_name"]) for s in seq}
    assert got == {(93, "GL1"), (65, "GL2")}


# ---------------------------------------------------------------------------
# ID-mapping remap (microstate reallocation)
# ---------------------------------------------------------------------------

def test_update_components_with_id_mapping_remaps_roles():
    T = _mk([{"resname": "GLU", "target": "GL1"}])
    components = {"GL1_chain": "A", "GL1_id": 65}
    updated = T.update_components_with_id_mapping(components, {("A", 65): 200})
    assert updated["GL1_id"] == 200


# ---------------------------------------------------------------------------
# Emit -> load -> register round-trip
# ---------------------------------------------------------------------------

def test_emit_and_load_roundtrip(tmp_path):
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    table = _wl_table(site, {("A", 45): "HD1", ("A", 65): "GL1", ("A", 93): "GL2"})
    path = emit_rename_transformer(
        table, name="test_dimn_glu_his", description="di-Mn Glu/His reuse",
        redox_state="reduced", spin_state="high_spin",
        forcefield_path="metal_sites/DiMn_Test",
        provenance={"library_path": "/x/y"}, target_dir=tmp_path,
    )
    assert path.exists()
    assert path.name == "test_dimn_glu_his.py"

    newly = load_user_transformers(tmp_path)
    assert "test_dimn_glu_his" in newly

    from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
        redox_transformer_registry,
    )
    cls = redox_transformer_registry.get_transformer("test_dimn_glu_his")
    assert cls is not None
    assert cls.REDOX_STATE == "reduced"
    assert cls.SPIN_STATE == "high_spin"
    # FORCEFIELD_PATH drives the Topology Generator's discover_forcefield_files.
    assert cls.FORCEFIELD_PATH == "metal_sites/DiMn_Test"

    # And it functions: auto-resolves the two Glu on a real site.
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    matched, missing = cls.match_components(site)
    assert missing == []
    assert matched["GL1_id"] == 65 and matched["GL2_id"] == 93


def test_load_is_idempotent(tmp_path):
    emit_rename_transformer([{"resname": "HID", "target": "HD1"}],
                            name="idem_check", description="d", target_dir=tmp_path)
    first = load_user_transformers(tmp_path)
    second = load_user_transformers(tmp_path)
    assert "idem_check" in first
    assert second == []  # already loaded, nothing new


def test_emit_empty_table_rejected(tmp_path):
    with pytest.raises(ValueError):
        emit_rename_transformer([], name="empty", description="d", target_dir=tmp_path)


def test_hook_composition_site_to_reusable_transformer(tmp_path):
    """Mirror what the MCPB-4 hook does: derive the rename table (with
    signatures) from a real site + rename map, emit, load, and confirm it
    re-resolves the two Glu on a fresh copy of the site."""
    site = _make_site(_dimn_atoms(), _dimn_bonds())
    # (chain, resid, old_resname) -> new_name, as the preprocessor builds it.
    pdb_rename_map = {
        ("A", 500, "MN"): "MN1",
        ("A", 501, "MN"): "MN2",
        ("A", 45, "HID"): "HD1",
        ("A", 108, "ASP"): "AS1",
        ("A", 65, "GLU"): "GL1",
        ("A", 93, "GLU"): "GL2",
    }
    mapping = {(chain, resid): new for (chain, resid, _old), new in pdb_rename_map.items()}
    table = _wl_table(site, mapping)

    emit_rename_transformer(table, name="hook_compose", description="d",
                            target_dir=tmp_path)
    load_user_transformers(tmp_path)

    from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
        redox_transformer_registry,
    )
    cls = redox_transformer_registry.get_transformer("hook_compose")
    matched, _ = cls.match_components(_make_site(_dimn_atoms(), _dimn_bonds()))
    # WL neighbourhoods resolve everything on this asymmetric site — including
    # the two Mn (different ligand sets) — with no ambiguity left to prompt.
    assert matched["GL1_id"] == 65 and matched["GL2_id"] == 93
    assert matched["MN1_id"] == 500 and matched["MN2_id"] == 501
    assert "_ambiguous" not in matched
