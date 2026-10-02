#!/usr/bin/env python3
"""
Docking receptor: explicit decisions in, Meeko PDBQT out.

Fixtures (tests/data/docking) are unmodified PDB entries or excerpts:
1HRC (cytochrome c: heme c bonded to Cys14/Cys17, Fe bound by His18 and
Met80, acetylated N-terminus), 2CBA (carbonic anhydrase II: Zn bound by
His94/His96 through NE2 and His119 through ND1, alternate locations,
REMARK 465), and the zinc sites of 2D5M, two of whose LINK records bind a Zn
to a His of a symmetry mate.

What these tests pin:
- nothing Meeko would decide silently is left to it without being reported:
  missing metal charges and alternate locations are refused, unresolved CCD
  problems are refused, residues that took Meeko's default are listed;
- covalent cofactors get templates and padders automatically from Meeko's own
  bond perception (heme c's two thioethers);
- a metal-bound nitrogen never carries the hydrogen (His119 of carbonic
  anhydrase is HIE, whatever pdb2pqr says without the zinc);
- LINK records to symmetry mates are reported, not treated as bonds, and
  follow the chosen alternate location.

Run with: pytest tests/test_docking_receptor.py
"""

import os
import shutil

import pytest

pytest.importorskip("rdkit", reason="RDKit is required for docking receptor preparation")
pytest.importorskip("meeko", reason="Meeko is required for docking receptor preparation")

from proprep.docking_prep import receptor_decisions as rd  # noqa: E402
from proprep.docking_prep.ccd_chemistry import component_from_block  # noqa: E402
from proprep.docking_prep.chemistry_edits import ChemistryEdit  # noqa: E402
from proprep.docking_prep.receptor_prep import prepare_receptor, read_atoms  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data", "docking")
HEME_CARBOXYLATES = [ChemistryEdit("remove_atom", ("H2A",)), ChemistryEdit("remove_atom", ("H2D",))]


@pytest.fixture
def pdb(tmp_path):
    def copy(name):
        target = tmp_path / name.replace(".txt", "")
        shutil.copy(os.path.join(DATA, name), target)
        return str(target)
    return copy


@pytest.fixture(scope="module")
def hec():
    with open(os.path.join(DATA, "ccd_HEC.cif.txt")) as handle:
        return component_from_block(handle.read(), "HEC")


def _cytochrome_c(path, hec, **overrides):
    decisions = dict(protonation={"A:18": "HID", "A:26": "HIE", "A:33": "HID"}, termini={"A:104": "C-terminus"},
                     metal_charges={("A:105", "FE"): 2}, altlocs={},
                     cofactor_edits={"HEC": HEME_CARBOXYLATES}, components={"HEC": hec})
    decisions.update(overrides)
    return prepare_receptor(path, ["A"], ["A:0", "A:105"], **decisions)


# --------------------------------------------------------------- building ---

def test_heme_c_gets_link_templates_from_meekos_own_bond_perception(pdb, hec):
    receptor = _cytochrome_c(pdb("1HRC.pdb.txt"), hec)
    rows = {r.res_id: r for r in receptor.residues}
    assert len(rows) == 107
    assert rows["A:14"].template == rows["A:17"].template == "CYS~SG"
    assert rows["A:105"].template == "HEC"
    assert rows["A:105M"].template == "FE2" and rows["A:105M"].source == "ion, charge chosen"
    assert rows["A:18"].template == "HID"
    assert rows["A:104"].template == "CGLU"
    assert sorted((l["a_atom"], l["b_atom"]) for l in receptor.links) == [("SG", "CAB"), ("SG", "CAC")]
    assert receptor.metals == [{"residue": "HEC A:105", "atom": "FE", "ion_residue": "A:105M",
                                "ion_name": "FE2", "charge": 2}]
    iron = [l for l in receptor.rigid_pdbqt.splitlines() if l[17:20] == "FE2"]
    assert len(iron) == 1 and iron[0][70:76].strip() == "+2.000"


def test_every_residue_meeko_decided_is_listed(pdb):
    path = pdb("2CBA.pdb.txt")
    altlocs = {r: "A" for r in {a.res_id for a in read_atoms(path) if a.altloc}}
    his = {a.res_id: "HIE" for a in read_atoms(path) if a.resname == "HIS"}
    receptor = prepare_receptor(path, ["A"], ["A:262"], protonation=his, termini={},
                                metal_charges={("A:262", "ZN"): 2}, altlocs=altlocs, cofactor_edits={})
    defaults = [r for r in receptor.residues if r.source == "Meeko default"]
    assert [(r.res_id, r.template) for r in defaults] == [("A:261", "CLYS")]   # no terminus decision given
    assert any("Meeko's own choice" in note and "A:261 LYS->CLYS" in note for note in receptor.notes)


def test_every_histidine_needs_a_decision(pdb, hec):
    with pytest.raises(ValueError, match="A:26, A:33"):
        _cytochrome_c(pdb("1HRC.pdb.txt"), hec, protonation={"A:18": "HID"})


def test_missing_metal_charge_is_refused_by_name(pdb, hec):
    with pytest.raises(ValueError, match="FE of HEC A:105"):
        _cytochrome_c(pdb("1HRC.pdb.txt"), hec, metal_charges={})


def test_unresolved_ccd_problem_is_refused(pdb, hec):
    with pytest.raises(ValueError, match="O2A"):
        _cytochrome_c(pdb("1HRC.pdb.txt"), hec, cofactor_edits={})


def test_alternate_locations_are_refused_until_chosen(pdb):
    path = pdb("2CBA.pdb.txt")
    with pytest.raises(ValueError, match=r"A:4 \(A, B\)"):
        prepare_receptor(path, ["A"], ["A:262"], protonation={}, termini={},
                         metal_charges={("A:262", "ZN"): 2}, altlocs={}, cofactor_edits={})


def test_kept_zinc_becomes_an_ion_of_the_chosen_charge(pdb):
    path = pdb("2CBA.pdb.txt")
    altlocs = {r: "A" for r in {a.res_id for a in read_atoms(path) if a.altloc}}
    his = {a.res_id: "HIE" for a in read_atoms(path) if a.resname == "HIS"}
    his.update({"A:94": "HID", "A:96": "HID", "A:119": "HIE"})
    receptor = prepare_receptor(path, ["A"], ["A:262"], protonation=his,
                                termini={}, metal_charges={("A:262", "ZN"): 2}, altlocs=altlocs, cofactor_edits={})
    rows = {r.res_id: r for r in receptor.residues}
    assert rows["A:262"].template == "ZN2"
    assert (rows["A:94"].template, rows["A:96"].template, rows["A:119"].template) == ("HID", "HID", "HIE")
    assert any(r.startswith("A:263 HOH") for r in receptor.removed)


# -------------------------------------------------------------- proposals ---

def test_metal_binding_decides_the_histidine_tautomer(pdb):
    path = pdb("2CBA.pdb.txt")
    altlocs = {r: "A" for r in {a.res_id for a in read_atoms(path) if a.altloc}}
    contacts = rd.metal_contacts(path, ["A:262"], altlocs)
    assert {(c.res_id, c.atom) for c in contacts} == {("A:94", "NE2"), ("A:96", "NE2"), ("A:119", "ND1")}
    assert all("LINK record" in c.source for c in contacts)
    overrides, problems = rd.metal_overrides(contacts)
    assert problems == []
    assert {k: v.value for k, v in overrides.items()} == {"A:94": "HID", "A:96": "HID", "A:119": "HIE"}


def test_links_to_symmetry_mates_are_reported_not_bonded(pdb):
    path = pdb("2D5M_zinc_excerpt.pdb.txt")
    zincs = ["A:201", "A:202", "A:203"]
    contacts = rd.metal_contacts(path, zincs, {"A:40": "A"})
    assert not any(c.resname == "HIS" for c in contacts)            # both His LINKs are to symmetry mates
    notes = rd.symmetry_mate_contacts(path, zincs)
    assert len(notes) == 2 and any("4565" in n for n in notes) and any("4555" in n for n in notes)


def test_link_records_follow_the_chosen_alternate_location(pdb):
    path = pdb("2D5M_zinc_excerpt.pdb.txt")
    zincs = ["A:201", "A:202", "A:203"]
    bound_a = {(c.metal_res, c.atom) for c in rd.metal_contacts(path, zincs, {"A:40": "A"}) if c.res_id == "A:40"}
    bound_b = {(c.metal_res, c.atom) for c in rd.metal_contacts(path, zincs, {"A:40": "B"}) if c.res_id == "A:40"}
    assert bound_a == {("A:201", "SG"), ("A:203", "SG")}
    assert bound_b == {("A:201", "SG")}


def test_termini_come_from_remark_465(pdb):
    proposals = rd.propose_termini(pdb("2CBA.pdb.txt"), ["A"])
    assert proposals["A:3"].value == "chain break"            # acetyl-Ser1 and Ser2 are unmodelled
    assert proposals["A:261"].value == "C-terminus"


def test_a_capped_end_needs_no_terminus_choice(pdb):
    proposals = rd.propose_termini(pdb("1HRC.pdb.txt"), ["A"])
    assert "A:1" not in proposals                              # ACE A:0 caps it


@pytest.mark.skipif(not (shutil.which("pdb2pqr") or shutil.which("pdb2pqr30")
                         or rd.find_executable("pdb2pqr")), reason="pdb2pqr is not installed")
def test_pdb2pqr_proposals_map_to_meeko_templates(pdb):
    path = pdb("2CBA.pdb.txt")
    altlocs = {r: "A" for r in {a.res_id for a in read_atoms(path) if a.altloc}}
    result = rd.propose_protonation(path, ["A"], 7.0, altlocs)
    assert result.ran, result.message
    assert result.proposals["A:94"].value == "HID" and result.proposals["A:96"].value == "HID"
    assert all(p.value in ("HID", "HIE", "HIP", "ASP", "ASH", "GLU", "GLH", "LYS", "LYN", "CYS", "CYX", "CYX-")
               for p in result.proposals.values())


@pytest.mark.skipif(not (shutil.which("pdb2pqr") or rd.find_executable("pdb2pqr")),
                    reason="pdb2pqr is not installed")
def test_pdb2pqr_failure_is_reported_not_raised(pdb):
    # pdb2pqr gives up on 1HRC's acetyl cap only when the cap is included; the protein-only
    # copy leaves HETATM out, so here a structure with no selected chain is the failure
    result = rd.propose_protonation(pdb("1HRC.pdb.txt"), ["Z"], 7.0, {})
    assert not result.ran and result.proposals == {} and "no ATOM records" in result.message


# ------------------------------------------------- clusters and other metals ---

@pytest.fixture(scope="module")
def ccd():
    def load(code):
        with open(os.path.join(DATA, f"ccd_{code}.cif.txt")) as handle:
            return component_from_block(handle.read(), code)
    return load


def test_iron_sulfur_clusters_become_iron_ions_and_sulfide_pieces(pdb, ccd):
    """2FDN: two [4Fe-4S] clusters; removing the irons leaves four separate sulfides each."""
    path = pdb("2FDN.pdb.txt")
    atoms = read_atoms(path)
    altlocs = {r: "A" for r in {a.res_id for a in atoms if a.altloc}}
    clusters = sorted({a.res_id for a in atoms if a.resname == "SF4"})
    irons = {(a.res_id, a.name): 2 for a in atoms if a.resname == "SF4" and a.element == "Fe"}
    overrides, problems = rd.metal_overrides(rd.metal_contacts(path, clusters, altlocs))
    assert problems == [] and len(overrides) == 8 and set(v.value for v in overrides.values()) == {"CYX-"}
    receptor = prepare_receptor(path, ["A"], clusters, protonation={k: v.value for k, v in overrides.items()},
                                termini={}, metal_charges=irons, altlocs=altlocs, cofactor_edits={},
                                components={"SF4": ccd("SF4")})
    rows = [r for r in receptor.residues if r.input_resname in ("SF4", "FE2")]
    assert sum(r.input_resname == "FE2" for r in rows) == 8
    pieces = [r for r in rows if r.input_resname == "SF4"]
    assert sorted(r.template for r in pieces) == sorted(["SF4~1", "SF4~2", "SF4~3", "SF4~4"] * 2)
    assert all("piece" in r.source and ("SF4 A:61" in r.source or "SF4 A:62" in r.source) for r in pieces)
    sulfides = [l for l in receptor.rigid_pdbqt.splitlines() if l.startswith("ATOM") and l[17:20] == "SF4"]
    assert len(sulfides) == 8 and all(l[70:76].strip() == "-2.000" for l in sulfides)


def _molybdenum_site(path, ccd, **extra):
    atoms = read_atoms(path)
    return prepare_receptor(path, ["A"], ["A:501", "A:505"], protonation={"A:140": "HIE", "A:185": "CYX-"},
                            termini={}, metal_charges={("A:505", "MO"): 6}, altlocs={}, cofactor_edits={},
                            components={"MTE": ccd("MTE")}, **extra)


def test_a_metal_vina_cannot_type_is_refused_by_name(pdb, ccd):
    with pytest.raises(ValueError, match=r"Mo at A:505 .*Vina has no Mo type"):
        _molybdenum_site(pdb("1SOX_mo_site_excerpt.pdb.txt"), ccd)


def test_such_a_metal_can_be_left_out_or_given_a_type_explicitly(pdb, ccd):
    path = pdb("1SOX_mo_site_excerpt.pdb.txt")
    without = _molybdenum_site(path, ccd, metal_types={("A:505", "MO"): "omit"})
    assert any("A:505 MO atom MO (metal left out by choice)" in r for r in without.removed)
    assert not any(l[17:20].startswith("MO") for l in without.rigid_pdbqt.splitlines())
    typed = _molybdenum_site(path, ccd, metal_types={("A:505", "MO"): "W"})
    molybdenum = [l for l in typed.rigid_pdbqt.splitlines() if l[17:20] == "MO6"]
    assert len(molybdenum) == 1 and molybdenum[0][77:].strip() == "W"
    assert any("typed W by choice" in note for note in typed.notes)
    assert not typed.links, "a metal must never be perceived as covalently bonded"


def test_vina_is_asked_which_metal_types_it_reads():
    from proprep.docking_prep.receptor_prep import vina_accepts_type
    assert vina_accepts_type("Fe") and vina_accepts_type("Cu")
    assert not vina_accepts_type("Mo")
