"""Both cysteines of a disulfide become CYX, also when residue IDs are remapped.

Seen on 6R2Q (two disulfides, twenty c-type hemes). The hemes split into three
residues each, which turns the residue ID mapping on for the whole structure.
The mapping then renumbered the first center of every site, the first cysteine
of each disulfide included (C:112 -> 814), but the disulfide's components kept
``cys1_id = 112``: its first rename matched no atom and said nothing. The
structure came out with CYS bonded to CYX, and tLEaP stopped later on "Could
not find bond parameter for atom types: SH - S".

The site below is the real C:112-C:118 disulfide of 6R2Q, loaded through the
same JSON import the Redox Site Preparer uses, and run through the real
analyzer and executor.
"""

import io
import json

import pytest
from rich.console import Console

from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
    RedoxSiteSpaceAnalyzer,
    TransformationExecutor,
)
from proprep.redoxsite_prep.transformation.transformers.disulfide import DisulfideTransformer
from proprep.structure_prep.comprehensive_redox_detector import _import_from_json

ATOMS = [
    (112, "N", "N", (66.74, 9.069, 4.563)),
    (112, "CA", "C", (68.175, 9.27, 4.711)),
    (112, "C", "C", (68.685, 8.063, 5.494)),
    (112, "O", "O", (68.202, 7.795, 6.604)),
    (112, "CB", "C", (68.447, 10.56, 5.488)),
    (112, "SG", "S", (70.176, 10.988, 5.768)),
    (118, "N", "N", (71.755, 9.256, 0.938)),
    (118, "CA", "C", (70.793, 10.165, 1.533)),
    (118, "C", "C", (70.762, 11.503, 0.812)),
    (118, "O", "O", (71.766, 12.22, 0.716)),
    (118, "CB", "C", (71.073, 10.26, 3.037)),
    (118, "SG", "S", (70.511, 11.709, 3.888)),
]


def _site(tmp_path, resnames=("CYS", "CYS")):
    name = dict(zip((112, 118), resnames))
    atom = lambda r, n, e, xyz: {"chain": "C", "resname": name[r], "resid": r, "atom_name": n,
                                 "element": e, "coordinates": list(xyz)}
    sg = {r: xyz for r, n, _, xyz in ATOMS if n == "SG"}
    site = {
        "site_id": "site_1", "structure_id": "temp", "site_type": "disulfide",
        "centers": [dict(atom(r, "SG", "S", sg[r]), center_type="redox_amino_acid") for r in (112, 118)],
        "atoms": [atom(*a) for a in ATOMS],
        "bonds": [{"atom1": {"chain": "C", "resname": name[112], "resid": 112, "atom_name": "SG"},
                   "atom2": {"chain": "C", "resname": name[118], "resid": 118, "atom_name": "SG"},
                   "bond_type": "interresidue", "chemical_type": "disulfide", "distance": 2.04,
                   "atom1_element": "S", "atom2_element": "S",
                   "atom1_coordinates": list(sg[112]), "atom2_coordinates": list(sg[118]),
                   "treatment": "bonded"}],
    }
    path = tmp_path / "sites.json"
    path.write_text(json.dumps({"source_pdb": "x.pdb",
                                "transformer_mappings": {"disulfide": "disulfide"},
                                "sites": [site]}))
    sites, _ = _import_from_json(str(path))
    return sites[0]


def _pdb_lines(resnames=("CYS", "CYS")):
    name = dict(zip((112, 118), resnames))
    return [f"ATOM  {i:5d} {n:<4s} {name[r]} C{r:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {e:>2s}\n"
            for i, (r, n, e, (x, y, z)) in enumerate(ATOMS, 1)]


def _resnames(lines):
    return {int(l[22:26]): l[17:20] for l in lines}


def _run(site, components, lines):
    out = io.StringIO()
    executor = TransformationExecutor(Console(file=out, width=300), None, verbose=False)
    sequence = DisulfideTransformer.get_transformation_sequence(components, {})
    lines, _ = executor.apply_transformation_sequence(site, sequence, lines)
    return lines, sequence, out.getvalue()


def test_both_cysteines_are_renamed(tmp_path):
    site = _site(tmp_path)
    components, missing = DisulfideTransformer.match_components(site)
    assert not missing
    lines, sequence, said = _run(site, components, _pdb_lines())
    assert _resnames(lines) == {112: "CYX", 118: "CYX"}
    assert [t["lines_modified"] for t in sequence] == [6, 6]
    assert said == ""


def test_a_disulfide_asks_for_no_residue_id_space(tmp_path):
    """It renames in place: nothing to renumber, whatever the other sites need."""
    site = _site(tmp_path)
    analysis = RedoxSiteSpaceAnalyzer(Console(file=io.StringIO())).analyze_transformation_space(
        [site], {"site_1": "disulfide"})
    assert analysis["space_requirements"] == {}
    assert analysis["needs_id_mapping"] is False


def test_components_follow_a_remapped_cysteine():
    """The safety net for any transformer: a role not called center/site/main_residue."""
    components = {"cys1_chain": "C", "cys1_id": 112, "cys2_chain": "C", "cys2_id": 118}
    updated = DisulfideTransformer.update_components_with_id_mapping(components, {("C", 112): 814})
    assert (updated["cys1_id"], updated["cys2_id"]) == (814, 118)
    # Same number on another chain is another residue.
    assert DisulfideTransformer.update_components_with_id_mapping(
        components, {("A", 112): 814}) == components


def test_a_required_step_that_matches_nothing_is_reported(tmp_path):
    """What happened on 6R2Q, now said out loud instead of passed over."""
    site = _site(tmp_path)
    components, _ = DisulfideTransformer.match_components(site)
    components["cys1_id"] = 814          # stale: the site still says 112
    lines, _, said = _run(site, components, _pdb_lines())
    assert _resnames(lines) == {112: "CYS", 118: "CYX"}
    assert "rename_cys1_to_cyx" in said and "matched no atoms" in said
    assert "rename_cys2_to_cyx" not in said


def test_a_cysteine_already_named_cyx_is_left_alone_without_a_report(tmp_path):
    site = _site(tmp_path, resnames=("CYX", "CYS"))
    components, _ = DisulfideTransformer.match_components(site)
    lines, sequence, said = _run(site, components, _pdb_lines(("CYX", "CYS")))
    assert [t["id"] for t in sequence] == ["rename_cys2_to_cyx"]
    assert _resnames(lines) == {112: "CYX", 118: "CYX"}
    assert said == ""
