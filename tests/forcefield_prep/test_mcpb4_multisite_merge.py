"""
Regression test for the mcpb-4 integration multi-site merge.

The original _checklist_mcpb_4_integration parsed only all_fingerprints[0], so a
protein with two independent metal sites (e.g. 7K0W's two Mn) had every residue
of site 2 silently dropped from naming, the FF library, and the PDB rename. The
fix parses EVERY site_*/models/standard.fingerprint and merges residues + M*/Y*
atom-type entries across all of them before generating names.

This test exercises the same parse -> merge -> generate_unique_residue_names
pipeline the checklist runs inline, using real fingerprint files.
"""

from proprep.forcefield_prep.mcpb.integration_utils import (
    parse_fingerprint,
    generate_unique_residue_names,
)


def _write_fingerprint(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _merge(fp_paths):
    """Replicates the inline merge in _checklist_mcpb_4_integration section D."""
    combined_residue_keys = []
    combined_atom_type_entries = []
    seen_keys, seen_entries = set(), set()
    for fp in fp_paths:
        fp_data = parse_fingerprint(str(fp), None)
        for key in fp_data["residues"].keys():
            if key not in seen_keys:
                seen_keys.add(key)
                combined_residue_keys.append(key)
        for entry in fp_data["atom_type_entries"]:
            if entry not in seen_entries:
                seen_entries.add(entry)
                combined_atom_type_entries.append(entry)
    return combined_residue_keys, combined_atom_type_entries


def test_two_sites_both_named_and_types_merged(tmp_path):
    # Site 1: a Mn center (renamed M1) + a coordinating His (HD1 stays orig)
    site1 = tmp_path / "site_1" / "models" / "standard.fingerprint"
    _write_fingerprint(site1, [
        "202-MN-MN     1  MN -> M1",
        "108-ASP-OD1   2  o  -> Y1",
        "108-ASP-CG    3  c  -> c",
    ])
    # Site 2: a second, independent Mn center + a different ligand atom
    site2 = tmp_path / "site_2" / "models" / "standard.fingerprint"
    _write_fingerprint(site2, [
        "203-MN-MN     1  MN -> M2",
        "245-GLU-OE1   2  o  -> Y2",
        "245-GLU-CD    3  c  -> c",
    ])

    keys, entries = _merge([site1, site2])

    # Both sites' residues survive the merge (the [0]-only bug dropped site 2).
    assert (202, "MN") in keys
    assert (108, "ASP") in keys
    assert (203, "MN") in keys
    assert (245, "GLU") in keys

    # Both sites' distinct M*/Y* types reach the merged addAtomTypes set.
    joined = " ".join(entries)
    for t in ("M1", "M2", "Y1", "Y2"):
        assert f'"{t}"' in joined, f"{t} missing from merged atom types: {entries}"

    # Names generated over the union are globally unique (no cross-site clash).
    name_map = generate_unique_residue_names(keys)
    names = list(name_map.values())
    assert len(names) == len(set(names)), f"collision in {names}"
    # Every merged residue got a name.
    assert set(name_map.keys()) == set(keys)


# --------------------------------------------------------------------------- #
# ligand GAFF frcmod collection (mcpb-4 section E)
# --------------------------------------------------------------------------- #
#
# A metal site's organic ligand (e.g. E4Z) needs its own parmchk2 GAFF frcmod
# deposited alongside the MCPB bonded frcmod; otherwise the ligand's atoms type
# at reuse but resolve no vdW/torsion parameters. _collect_ligand_frcmods finds
# that file (small_molecule_params_<RES>/<res>.frcmod) for the site's ligands
# and NOT for protein sidechains or metals. It once crashed with a NameError on
# os.getcwd() (no local `import os`); this locks the import + the discovery in.

from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor


class _WS:
    def __init__(self, d):
        self.d = d

    def get(self, k, default=None):
        return self.d.get(k, default)


def _bare_preprocessor(workspace):
    sp = StructurePreprocessor.__new__(StructurePreprocessor)
    sp.workspace = workspace
    sp._final_pdb = None
    return sp


def test_collect_ligand_frcmods_finds_only_organic_ligand(tmp_path):
    # Lay out the small-molecule parameterizer's on-disk product for one ligand.
    lig_dir = tmp_path / "small_molecule_params_E4Z"
    lig_dir.mkdir()
    (lig_dir / "e4z.frcmod").write_text("MASS\nc 12.01\n\nNONBON\n  c 1.9 0.086\n")
    # A protein sidechain has no small_molecule_params_* dir -> must not resolve.

    sp = _bare_preprocessor(_WS({
        "prepared_pdb": str(tmp_path / "structure.pdb"),
        "preprocessing_organic_ff": {},
        "preprocessing_organometallic_ff": {},
    }))

    found = sp._collect_ligand_frcmods({"HID", "GLU", "ASP", "ILE", "E4Z", "MN"})

    assert set(found) == {"E4Z"}
    assert found["E4Z"] == str(lig_dir / "e4z.frcmod")


def test_collect_ligand_frcmods_prefers_workspace_over_disk(tmp_path):
    # When the live organic-FF result carries a frcmod_file, use it directly.
    ws_frcmod = tmp_path / "ws_e4z.frcmod"
    ws_frcmod.write_text("MASS\n")
    sp = _bare_preprocessor(_WS({
        "prepared_pdb": str(tmp_path / "structure.pdb"),
        "preprocessing_organic_ff": {"E4Z": {"frcmod_file": str(ws_frcmod)}},
        "preprocessing_organometallic_ff": {},
    }))

    found = sp._collect_ligand_frcmods({"E4Z"})
    assert found == {"E4Z": str(ws_frcmod)}


def test_collect_ligand_frcmods_empty_when_no_ligand(tmp_path):
    # All-protein site: nothing to collect, and no crash (the os.getcwd path).
    sp = _bare_preprocessor(_WS({
        "prepared_pdb": str(tmp_path / "structure.pdb"),
        "preprocessing_organic_ff": {},
        "preprocessing_organometallic_ff": {},
    }))
    assert sp._collect_ligand_frcmods({"HID", "GLU", "ASP"}) == {}
