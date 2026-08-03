"""
Regression test for the MODELLER chain-rename + terminal-capping redox-site
mapping bug.

Background
----------
When MODELLER repairs a multi-chain structure it renames chains consecutively
(original B,D -> A,B) and numbers residues GLOBALLY across chains (A:1..N,
B:N+1..). Terminal capping then inserts ACE/NME and `_renumber_structure`
restarts each chain at 1. Two defects made redox-site synchronization fail for
every atom in a renamed chain:

  * Bug A: the final (post-cap) numbering map was predicted from the global
    cross-chain numbering, so it never matched the per-chain 1..N numbering on
    disk -> get_final_identity returned a residue number that did not exist.
  * Bug B: cap placement / map-keying used the ORIGINAL chain id against the
    MODELLER-renamed structure, so a renamed chain's caps were dropped or
    mis-applied to a coincidentally-named chain.

The fix makes CappingHandler._renumber_structure the single source of truth for
final numbering (it records old->new per renamed chain as it writes the file)
and translates segment chain ids through mapper.chain_mapping before touching
the renamed structure. _find_atom_coords was also hardened to scan all residues
carrying a number (cap insertion can leave an NME and an ion sharing a number
before renumbering merges them).

This test reproduces the exact mechanism end-to-end through the real classes,
without requiring a MODELLER run.
"""
import io
import os
import tempfile
import unittest
from types import SimpleNamespace

from Bio.PDB import PDBParser
from rich.console import Console

from proprep.structure_prep.structure_completeness import (
    ResidueMapper, CappingHandler, RedoxSiteSync, RepairPlan, MissingSegment,
)


def _res_block(record, serial0, chain, resseq, resname, atoms):
    lines = []
    for i, (name, elem, (x, y, z)) in enumerate(atoms):
        name4 = name if len(name) == 4 else (" %-3s" % name)
        lines.append(
            "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
            % (record, serial0 + i, name4, "", resname, chain, resseq, "",
               x, y, z, 1.0, 0.0, elem)
        )
    return lines, serial0 + len(atoms)


def _backbone(x, oxt=False):
    atoms = [("N", "N", (x, 0., 0.)), ("CA", "C", (x, 1., 0.)),
             ("C", "C", (x, 2., 0.)), ("O", "O", (x, 3., 0.))]
    if oxt:  # free C-terminus carboxylate (must be stripped when NME caps it)
        atoms.append(("OXT", "O", (x, 4., 0.)))
    return atoms


def _build_pdb(chains):
    lines, serial = [], 1
    for chain_id, residues in chains:
        for record, resseq, resname, atoms in residues:
            blk, serial = _res_block(record, serial, chain_id, resseq, resname, atoms)
            lines.extend(blk)
        lines.append("TER\n")
    lines.append("END\n")
    return "".join(lines)


def _parse(pdb_text, name):
    return PDBParser(QUIET=True).get_structure(name, io.StringIO(pdb_text))


class ChainRenameCapSyncTest(unittest.TestCase):
    def setUp(self):
        # Quiet console (suppress rich output during tests).
        self.console = Console(file=io.StringIO())

    def test_renamed_chain_site_resolves_after_capping(self):
        # ORIGINAL: chain B (ALA 1..10), chain D (GLU 50,51,52 + Ca 1005 + HOH 1009).
        # The last protein residue (GLU 52) carries an OXT (free C-terminus).
        original = _parse(_build_pdb([
            ("B", [("ATOM", n, "ALA", _backbone(float(n))) for n in range(1, 11)]),
            ("D", [("ATOM", n, "GLU", _backbone(float(n), oxt=(n == 52))) for n in (50, 51, 52)]
                  + [("HETATM", 1005, "CA", [("CA", "CA", (60., 0., 0.))]),
                     ("HETATM", 1009, "HOH", [("O", "O", (61., 0., 0.))])]),
        ]), "orig")

        # MODELLER OUTPUT: chains A (1..10), B (11,12,13 + Ca 14 + HOH 15) GLOBAL numbering.
        # GLU 13 (the C-terminus that gets NME) keeps its OXT; Ca 14 collides in
        # number with the inserted NME (which is also positioned at 14).
        modeller_out = _parse(_build_pdb([
            ("A", [("ATOM", n, "ALA", _backbone(float(n))) for n in range(1, 11)]),
            ("B", [("ATOM", n, "GLU", _backbone(float(n), oxt=(n == 13))) for n in (11, 12, 13)]
                  + [("HETATM", 14, "CA", [("CA", "CA", (60., 0., 0.))]),
                     ("HETATM", 15, "HOH", [("O", "O", (61., 0., 0.))])]),
        ]), "modeller")

        mapper = ResidueMapper(self.console)
        mapper.build_modeller_mapping(original, {}, {})
        # Sanity: MODELLER renames B->A, D->B and numbers globally.
        self.assertEqual(mapper.chain_mapping, {"B": "A", "D": "B"})
        self.assertEqual(mapper.residue_mappings["B"][1005], 14)

        # Cap chain D (both termini) and chain B (N) — ORIGINAL chain ids.
        plan = RepairPlan()
        mk = lambda c, t: MissingSegment(chain_id=c, residues=[], start_num=0,
                                         end_num=0, is_terminal=True, terminal_type=t)
        plan.segments_to_cap = [mk("D", "N"), mk("D", "C"), mk("B", "N")]

        with tempfile.TemporaryDirectory() as td:
            final_structure = CappingHandler(self.console).add_caps(
                modeller_out, plan, mapper, os.path.join(td, "final.pdb"))

        # Bug A: final_mappings now restart per-chain at 1 (matches on-disk file).
        # Global Ca number 14 must map to a small per-chain final number, NOT 14.
        self.assertIn("B", mapper.final_mappings)
        self.assertIn(14, mapper.final_mappings["B"])
        self.assertLess(mapper.final_mappings["B"][14], 14)

        # Build a redox site referencing ORIGINAL (chain D) identities.
        mkatom = lambda c, r, rn, nm, el, co: SimpleNamespace(
            chain=c, resid=r, resname=rn, atom_name=nm, element=el,
            coords=co, insertion_code=" ")
        site_atoms = [
            mkatom("D", 1005, "CA", "CA", "CA", (60., 0., 0.)),   # Ca: NME/Ca number collision
            mkatom("D", 50, "GLU", "N", "N", (50., 0., 0.)),
            mkatom("D", 50, "GLU", "CA", "C", (50., 1., 0.)),
            mkatom("D", 52, "GLU", "C", "C", (52., 2., 0.)),
            mkatom("D", 1009, "HOH", "O", "O", (61., 0., 0.)),
        ]
        site = SimpleNamespace(
            site_id="site_test", site_type="no_transformation",
            atoms=site_atoms, bonds=[],
            centers=[mkatom("D", 1005, "CA", "CA", "CA", (60., 0., 0.))],
            coord_to_pdb={}, residue_groups={})

        summary = RedoxSiteSync(self.console).synchronize_sites(
            [site], final_structure, mapper)

        # Every atom + the center must resolve (the bug left all "not found").
        self.assertEqual(summary["atoms_updated"], len(site_atoms))
        self.assertEqual(summary["centers_updated"], 1)
        # Atoms must now carry the renamed chain id.
        self.assertTrue(all(a.chain == "B" for a in site_atoms))
        # Ca (which collides in number with the inserted NME) still resolves.
        ca = site_atoms[0]
        self.assertEqual(ca.chain, "B")

        # --- Structural correctness of the capped chain B (renamed orig D) ---
        chain_b = next(iter(final_structure))["B"]
        by_resnum = {}
        for residue in chain_b:
            by_resnum.setdefault(residue.id[1], []).append(residue.resname.strip())

        # (1) NME must exist, and the residue it bonds to must have lost its OXT.
        nme_res = [r for r in chain_b if r.resname.strip() == "NME"]
        self.assertEqual(len(nme_res), 1)
        glu_res = [r for r in chain_b if r.resname.strip() == "GLU"]
        self.assertTrue(glu_res)
        for g in glu_res:
            self.assertNotIn("OXT", {a.name for a in g},
                             "OXT must be stripped from the NME-capped residue")

        # (2) No residue number is shared by two different residues (the NME vs Ca
        #     / cofactor collision used to merge them into one number).
        dups = {n: names for n, names in by_resnum.items() if len(names) > 1}
        self.assertEqual(dups, {}, f"duplicate residue numbers: {dups}")
        nme_num = nme_res[0].id[1]
        ca_num = [r.id[1] for r in chain_b if r.resname.strip() == "CA"][0]
        self.assertNotEqual(nme_num, ca_num)


class CapsOnlyIdentityMappingTest(unittest.TestCase):
    """Regression test for the caps-only (no-MODELLER) double-renumbering bug.

    When MODELLER does not run (termini just need capping), the existing
    residues keep their original numbers until `_renumber_structure` restarts
    each chain at 1 and records that map in `final_mappings` (Step 3 of
    get_final_identity). `build_identity_mapping` must therefore leave Step 2 as
    the IDENTITY. The bug had it ALSO pre-renumber the standard residues to
    1..N, so the two maps composed and double-offset every protein residue:
    e.g. a structure whose present residues start at 83 mapped CYS 86 -> 4
    (Step 2) and then resolved residue 4's coordinates, landing the redox-site
    CYS on the wrong residue (its SG, which only exists on the real CYS, then
    came back "not found"). HETATM cofactors were spared because they are
    excluded from residue_mappings and so skip Step 2 — exactly the asymmetry
    seen in the field report (cofactors correct, CYS doubled).
    """

    def setUp(self):
        self.console = Console(file=io.StringIO())

    def test_caps_only_site_resolves_without_double_offset(self):
        # Present residues start at 83 (N-terminal 1..82 are missing), so a
        # 1..N pre-renumber differs from the identity and the bug bites.
        # ALA 83,84,85 (backbone) + CYS 86 (backbone + SG) + HEM cofactor 1155.
        def cys_atoms(x):
            return _backbone(x) + [("SG", "S", (x, 5., 0.))]

        structure = _parse(_build_pdb([
            ("A",
             [("ATOM", n, "ALA", _backbone(float(n))) for n in (83, 84, 85)]
             + [("ATOM", 86, "CYS", cys_atoms(86.0)),
                ("HETATM", 1155, "HEM", [("FE", "FE", (90., 0., 0.))])]),
        ]), "caps_only")

        # Caps-only path: identity mapping keyed by the REAL chain id.
        mapper = ResidueMapper(self.console)
        mapper.build_identity_mapping(structure)
        # Step 2 must be the identity (orig -> orig), NOT a 1..N renumber.
        self.assertEqual(mapper.residue_mappings["A"][86], 86)

        # Cap the N-terminus only (ACE before residue 83). No NME -> no OXT games.
        plan = RepairPlan()
        plan.segments_to_cap = [MissingSegment(
            chain_id="A", residues=[], start_num=0, end_num=0,
            is_terminal=True, terminal_type="N")]

        with tempfile.TemporaryDirectory() as td:
            final_structure = CappingHandler(self.console).add_caps(
                structure, plan, mapper, os.path.join(td, "final.pdb"))

        # _renumber_structure restarts at 1: ACE->1, 83->2, 84->3, 85->4,
        # CYS 86->5, HEM 1155->6. final_mappings is keyed by ORIGINAL numbers.
        self.assertEqual(mapper.final_mappings["A"][86], 5)
        self.assertEqual(mapper.final_mappings["A"][1155], 6)

        # Build a redox site referencing the ORIGINAL CYS 86 (incl. its SG) and
        # the HEM cofactor.
        mkatom = lambda r, rn, nm, el, co: SimpleNamespace(
            chain="A", resid=r, resname=rn, atom_name=nm, element=el,
            coords=co, insertion_code=" ")
        site_atoms = [
            mkatom(86, "CYS", "N", "N", (86., 0., 0.)),
            mkatom(86, "CYS", "CA", "C", (86., 1., 0.)),
            mkatom(86, "CYS", "SG", "S", (86., 5., 0.)),   # only on the real CYS
            mkatom(1155, "HEM", "FE", "FE", (90., 0., 0.)),
        ]
        site = SimpleNamespace(
            site_id="site_test", site_type="heme_cys_axial",
            atoms=site_atoms, bonds=[], centers=[],
            coord_to_pdb={}, residue_groups={})

        summary = RedoxSiteSync(self.console).synchronize_sites(
            [site], final_structure, mapper)

        # Every atom resolves — the bug left SG "not found" and mis-placed the
        # backbone atoms on residue 4.
        self.assertEqual(summary["atoms_updated"], len(site_atoms))

        # CYS lands on its TRUE post-cap number (5), NOT the double-offset 4.
        cys = [a for a in site_atoms if a.resname == "CYS"]
        self.assertTrue(all(a.resid == 5 for a in cys),
                        f"CYS double-offset: {[a.resid for a in cys]}")
        # SG resolved to the real CYS's SG coordinates (not a backbone atom of a
        # neighbouring ALA).
        sg = next(a for a in cys if a.atom_name == "SG")
        self.assertEqual(sg.coords, (86., 5., 0.))
        # Cofactor still maps correctly (it always did; guards against regression).
        hem = next(a for a in site_atoms if a.resname == "HEM")
        self.assertEqual(hem.resid, 6)


class CapGeometryTest(unittest.TestCase):
    """ACE/NME placeholder geometry must not be collinear: a 180deg layout makes
    bond vectors antiparallel, and planarity tests (e.g. PROPKA SYBYL typing)
    then divide by a zero-length cross product and crash."""

    def setUp(self):
        self.console = Console(file=io.StringIO())

    @staticmethod
    def _coords(lines):
        out = []
        for ln in lines:
            out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        return out

    def _assert_non_collinear(self, pts):
        import numpy as np
        a, b, c = (np.array(p) for p in pts[:3])
        normal = np.cross(b - a, c - a)
        self.assertGreater(float(np.linalg.norm(normal)), 1e-6,
                           "cap atoms are collinear (degenerate geometry)")

    def test_ace_and_nme_caps_non_collinear(self):
        handler = CappingHandler(self.console)
        ref = ("ATOM      1  N   ALA A   2      10.000  20.000  30.000"
               "  1.00  0.00           N  ")
        ace = handler._create_ace_cap("A", 1, 100, ref)
        nme = handler._create_nme_cap("A", 3, 200, ref)
        self.assertEqual(len(ace), 3)
        # NME emits only N + H; tLEaP builds the methyl carbon (named CH3 or C
        # per the loaded FF) and its hydrogens. The off-axis H must still keep
        # N's neighbours (attachment C, H) non-collinear for PROPKA.
        self.assertEqual(len(nme), 2)
        self._assert_non_collinear(self._coords(ace))
        ref_c = (float(ref[30:38]), float(ref[38:46]), float(ref[46:54]))
        self._assert_non_collinear([ref_c] + self._coords(nme))


if __name__ == "__main__":
    unittest.main()
