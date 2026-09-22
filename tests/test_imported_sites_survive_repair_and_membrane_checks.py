"""One character, three modules later.

6R2Q, 2026-09-20. Redox sites imported from a site file, structure repaired
(MODELLER renumbers everything), sites prepared, membrane built, tLEaP:
"bond: Argument #2 is of type String". Traced back:

1. The fixer's sync moved each heme site's ATOMS to the new numbers (A:901 to
   A:274) but left its CENTER at A:901. A heme's center is the residue centroid,
   and that branch compared insertion codes as text: BioPython's ' ' against
   the '' an imported center carries (the site file did not store it). Nothing
   was said: a green tick for the atoms.
2. The Redox Site Preparer addresses the site by its center, so every step on
   the heme named a residue the site did not contain, did nothing, and said
   nothing. The structure came out half transformed: 20 residues still HEC.
3. tLEaP failed on them in the Membrane Builder's hydrogen pass. The pass
   tested only that its topology files EXIST; the morning's were still there,
   so the membrane was built around the morning's protein, numbered
   differently from the bond directives.
"""

import io
import types

import pytest
from Bio.PDB import PDBParser
from rich.console import Console

from proprep.membrane_prep.membrane_builder import MembraneBuilderModule
from proprep.redoxsite_prep.transformation.redox_transformer_framework import TransformationExecutor
from proprep.structure_prep import comprehensive_redox_detector as crd
from proprep.structure_prep import structure_completeness as sc

HEME = [("FE", "FE", 23.639, 107.015, 3.850), ("NA", "N", 25.100, 106.000, 4.900), ("NB", "N", 22.400, 105.300, 3.100)]
CYS = [("CA", "C", 68.175, 9.270, 4.711), ("SG", "S", 70.176, 10.988, 5.768)]


def _pdb(path, heme_number, cys_number, heme_name="HEC"):
    lines, i = [], 0
    for name, element, x, y, z in CYS:
        i += 1
        lines.append(f"ATOM  {i:5d} {name:<4s} CYS A{cys_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n")
    for name, element, x, y, z in HEME:
        i += 1
        lines.append(f"HETATM{i:5d} {name:<4s} {heme_name} A{heme_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n")
    path.write_text("".join(lines) + "END\n")
    return str(path)


def _site(insertion_code):
    site = crd.RedoxSite("site_3", "temp")
    site.site_type = "heme"
    centroid = tuple(round(sum(a[k] for a in HEME) / len(HEME), 3) for k in (2, 3, 4))
    site.add_center(crd.RedoxCenter(chain="A", resname="HEC", resid=901, atom_name=None, insertion_code=insertion_code,
                                    altloc="", coords=centroid, center_type=crd.CenterType.ORGANOMETALLIC_COFACTOR,
                                    element=None))
    for name, element, x, y, z in HEME:
        site.add_atom(crd.RedoxSiteAtom(chain="A", resname="HEC", resid=901, atom_name=name, coords=(x, y, z), element=element))
    return site


def _mapper():
    mapper = sc.ResidueMapper(Console(file=io.StringIO()))
    mapper.chain_mapping = {"A": "A"}
    mapper.residue_mappings = {"A": {67: 7, 901: 274}}
    return mapper


def _sync(site, tmp_path, heme_name="HEC"):
    repaired = PDBParser(QUIET=True).get_structure("r", _pdb(tmp_path / "repaired.pdb", 274, 7, heme_name))
    out = io.StringIO()
    sc.RedoxSiteSync(Console(file=out, width=300)).synchronize_sites([site], repaired, _mapper())
    return out.getvalue()


# ── 1. the fixer's sync ─────────────────────────────────────────────────

@pytest.mark.parametrize("insertion_code", ["", " "], ids=["imported from a site file", "detected on a structure"])
def test_a_centroid_center_follows_its_residue_whatever_blank_it_carries(tmp_path, insertion_code):
    site = _site(insertion_code)
    said = _sync(site, tmp_path)
    assert (site.centers[0].chain, site.centers[0].resid) == ("A", 274)
    assert {a.resid for a in site.atoms} == {274} and "could not be placed" not in said


def test_a_center_that_cannot_be_placed_is_not_hidden_under_a_green_tick(tmp_path):
    site = _site(" ")
    said = _sync(site, tmp_path, heme_name="XYZ")            # the repaired structure has no HEC at A:274
    assert site.centers[0].resid == 901
    assert "could not be placed" in said and "KEEP THEIR OLD NUMBERING" in said and "HEC A:901" in said


def test_imported_centers_carry_a_blank_insertion_code_as_detected_ones_do(tmp_path):
    import json
    path = tmp_path / "sites.json"
    path.write_text(json.dumps({"source_pdb": "x", "transformer_mappings": {}, "sites": [{
        "site_id": "s", "structure_id": "t", "site_type": "heme", "atoms": [], "bonds": [],
        "centers": [{"chain": "A", "resname": "HEC", "resid": 901, "atom_name": None, "element": None,
                     "center_type": "organometallic_cofactor", "coordinates": [1.0, 2.0, 3.0]}]}]}))
    sites, _ = crd._import_from_json(str(path))
    assert sites[0].centers[0].insertion_code == " "


# ── 2. the Preparer ─────────────────────────────────────────────────────

def _run_step(site, selector):
    out = io.StringIO()
    executor = TransformationExecutor(Console(file=out, width=400), None, verbose=False)
    step = {"id": "rename", "description": "Rename", "selector": selector, "action": {"change_residue_name": "HCO"}}
    lines = [f"HETATM{i:5d} {n:<4s} HEC A 274    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n" for i, (n, _, x, y, z) in enumerate(HEME, 1)]
    executor.apply_transformation_sequence(site, [step], lines)
    return out.getvalue()


def test_a_step_naming_a_residue_the_site_does_not_contain_is_reported(tmp_path):
    site = _site(" ")
    said = _run_step(site, {"chain_id": "A", "residue_id": 287, "residue_name": "HEC"})
    assert "refers to residue A:287" in said and "does not contain" in said and "A:901" in said


def test_a_step_whose_residue_is_there_but_filtered_out_stays_quiet(tmp_path):
    """rename_hem_to_hec on a heme already named HEC: nothing to do, and nothing wrong."""
    assert _run_step(_site(" "), {"chain_id": "A", "residue_id": 901, "residue_name": "HEM"}) == ""


# ── 3. the Membrane Builder's hydrogen pass ─────────────────────────────

def _atoms_file(path, residues, hydrogens=False):
    lines, i = [], 0
    for number, name, atoms in residues:
        for atom, x in atoms:
            i += 1
            lines.append(f"ATOM  {i:5d} {atom:<4s} {name:>3s}  {number:4d}    {x:8.3f}{1.0:8.3f}{2.0:8.3f}  1.00  0.00\n")
            if hydrogens:
                i += 1
                lines.append(f"ATOM  {i:5d} {'H' + atom:<4s} {name:>3s}  {number:4d}    {x + 0.5:8.3f}{1.0:8.3f}{2.0:8.3f}  1.00  0.00\n")
    path.write_text("".join(lines))
    return str(path)


GIVEN = [(1, "ACE", [("C", 1.0)]), (2, "SER", [("N", 2.0), ("CA", 3.0)]), (3, "HOH", [("O", 9.0)])]
check = MembraneBuilderModule._hydrogen_pass_changed_the_structure


def test_hydrogens_added_and_waters_renamed_is_the_same_structure(tmp_path):
    written = [(1, "ACE", GIVEN[0][2]), (2, "SER", GIVEN[1][2]), (3, "WAT", GIVEN[2][2])]
    assert check(_atoms_file(tmp_path / "a.pdb", GIVEN), _atoms_file(tmp_path / "b.pdb", written, hydrogens=True)) is None


def test_a_stale_output_from_an_earlier_structure_is_refused(tmp_path):
    """The morning's protein: no cap, so one residue fewer."""
    stale = [(1, "SER", GIVEN[1][2]), (2, "WAT", GIVEN[2][2])]
    problem = check(_atoms_file(tmp_path / "a.pdb", GIVEN), _atoms_file(tmp_path / "b.pdb", stale, hydrogens=True))
    assert "1 of its 4 heavy atoms are not in the output" in problem and "ACE" in problem


def test_residues_at_other_positions_are_refused_because_bonds_are_addressed_by_position(tmp_path):
    reordered = [GIVEN[2], GIVEN[0], GIVEN[1]]
    problem = check(_atoms_file(tmp_path / "a.pdb", GIVEN), _atoms_file(tmp_path / "b.pdb", reordered))
    assert "another position in the file" in problem and "bond directives" in problem


def test_unknown_residues_lead_the_tleap_error_summary():
    output = ("Unknown residue: HEC   number: 1532   type: Terminal/beginning\n"
              "Unknown residue: HEC   number: 1533   type: Nonterminal\n"
              "/x/teLeap: Error!\nComparing atoms\n        .R<HCO 1544>.A<NA 2>,\n")
    assert MembraneBuilderModule._tleap_error_messages(output) == [
        "Unknown residue: HEC (2 of them): no loaded library defines it", "Comparing atoms"]


def test_the_hydrogen_pass_removes_an_earlier_runs_outputs_before_it_runs(monkeypatch, tmp_path):
    import subprocess
    from proprep.membrane_prep import membrane_builder as mb
    protein = _atoms_file(tmp_path / "protein.pdb", GIVEN)
    for name in ("pretleap_hydrogen.prmtop", "pretleap_hydrogen.rst7", "protein_with_h.pdb"):
        (tmp_path / name).write_text("from this morning")
    module = mb.MembraneBuilderModule()
    out = io.StringIO()
    module.processor = types.SimpleNamespace(console=Console(file=out, width=300))
    module.config.protein_pdb = protein
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 1, stdout="Unknown residue: HEC   number: 3\n", stderr=""))
    workspace = types.SimpleNamespace(get=lambda key, default=None: default)

    assert module._run_pre_tleap_hydrogen_pass(workspace, str(tmp_path)) is None     # tLEaP failed, and that is what is reported
    assert not (tmp_path / "pretleap_hydrogen.prmtop").exists() and not (tmp_path / "protein_with_h.pdb").exists()
    assert "did not produce topology files" in out.getvalue() and "Unknown residue: HEC" in out.getvalue()
