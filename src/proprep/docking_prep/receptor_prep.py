"""
Receptor PDBQT with Meeko, from decisions made explicitly.

Meeko rebuilds a receptor's hydrogens from residue templates, and left to
itself it chooses silently: every HIS becomes HIE (which on carbonic anhydrase
puts a hydrogen on the zinc-bound nitrogen of two of the three zinc ligands),
an unknown cofactor is downloaded from the CCD on each run and deprotonated by
rule, and iron is Fe3+. Here every such choice is an input:

* ``protonation``   residue -> Meeko template key (HID, ASH, CYX-, ...)
* ``termini``       residue -> "N-terminus" / "C-terminus" / "chain break"
* ``cofactor_edits`` residue name -> ChemistryEdits applied to its CCD entry
* ``metal_charges`` (residue, atom name) -> formal charge of every kept metal
* ``altlocs``       residue -> alternate location, for every residue that has them
* ``altloc_fills``  residue -> {atom: alternate} for atoms its chosen alternate does not model

``prepare_receptor`` refuses rather than guesses when any of these is missing,
and returns a table recording, for every residue, the template used and where
that decision came from. Proposals for the decisions (pdb2pqr/PROPKA, metal
coordination, SEQRES/REMARK 465) are made elsewhere; this module only builds.

Residues are addressed as Meeko does: "chain:number[insertion code]", e.g. "A:105".

Metals are split out of their residues into ions of their own (Meeko does not
match a template in which the metal is a disconnected fragment). A split
metal keeps its chain and residue number and gets a free insertion code, so
the heme iron of A:105 becomes A:105M, residue name FE2 or FE3 by charge.
"""

from __future__ import annotations

import logging
import os
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from rdkit import Chem

from .ccd_chemistry import CcdComponent, fetch_component, is_metal, is_metal_symbol
from .chemistry_edits import ChemistryEdit, apply_edits
from . import receptor_templates as rt

logger = logging.getLogger(__name__)

TERMINUS_CHOICES = ("N-terminus", "C-terminus", "chain break")
_WATER = {"HOH", "WAT", "DOD"}


# ------------------------------------------------------------ PDB records ---

@dataclass
class AtomLine:
    record: str        # "ATOM" or "HETATM"
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    xyz: Tuple[float, float, float]
    element: str
    line: str
    occupancy: Optional[float] = None

    @property
    def res_id(self) -> str:
        return f"{self.chain}:{self.resseq}{self.icode}"


def read_atoms(pdb_path: str) -> List[AtomLine]:
    """ATOM/HETATM records of the first model."""
    from proprep.utils.pdb_format import element_from_name_field
    atoms = []
    with open(pdb_path) as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                element = element_from_name_field(line[12:16])
            atoms.append(AtomLine(
                record=line[:6].strip(), serial=int(line[6:11]), name=line[12:16].strip(),
                altloc=line[16].strip(), resname=line[17:20].strip(), chain=line[21],
                resseq=int(line[22:26]), icode=line[26].strip(),
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                element=element[0].upper() + element[1:].lower() if element else "",
                line=line.rstrip("\n"), occupancy=_occupancy(line)))
    return atoms


def _occupancy(line: str) -> Optional[float]:
    try:
        return float(line[54:60])
    except ValueError:
        return None


def chosen_alternate(atom: AtomLine, altlocs: Dict[str, str],
                     altloc_fills: Optional[Dict[str, Dict[str, str]]] = None) -> bool:
    """Whether ``atom`` survives the alternate-location choices.

    An unlabelled atom always does; a labelled one when its residue's chosen
    alternate is its letter, or when the chosen alternate does not model that
    atom and the fill plan takes it from this letter (see ``altloc_picker``).
    """
    from proprep.structure_prep.altloc_picker import keeps_atom
    return keeps_atom(atom.name, atom.altloc, altlocs.get(atom.res_id), (altloc_fills or {}).get(atom.res_id))


def _format_atom(atom: AtomLine, serial: int, name: str, resname: str, icode: str) -> str:
    from proprep.utils.pdb_format import atom_name_field
    x, y, z = atom.xyz
    return (f"{atom.record:<6}{serial % 100000:5d} {atom_name_field(name, atom.element)} "
            f"{resname:>3} {atom.chain}{atom.resseq:4d}{icode or ' '}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {atom.element.upper():>2}")


# --------------------------------------------------------------- results ---

@dataclass
class ResidueDecision:
    res_id: str
    input_resname: str
    template: str
    source: str        # where the template choice came from


@dataclass
class PreparedReceptor:
    rigid_pdbqt: str
    flex_pdbqt: str
    residues: List[ResidueDecision]
    metals: List[Dict[str, object]]          # original residue/atom -> ion residue, charge
    links: List[Dict[str, object]]           # covalent bonds between residues that needed templates
    removed: List[str]                       # residues left out, with the reason
    notes: List[str] = field(default_factory=list)


# ------------------------------------------------------------------ build ---

def _meeko_default_keys() -> set:
    from meeko import ResidueChemTemplates
    chem = ResidueChemTemplates.create_from_defaults()
    return set(chem.residue_templates) | set(chem.ambiguous)


def _free_icode(taken: set) -> str:
    for code in "MNOPQRSTUVWXYZ":
        if code not in taken:
            return code
    raise ValueError("No free insertion code for a split metal.")


def prepare_receptor(
    pdb_path: str,
    chains: Sequence[str],
    keep_hetero: Sequence[str],
    *,
    protonation: Dict[str, str],
    termini: Dict[str, str],
    metal_charges: Dict[Tuple[str, str], int],
    altlocs: Dict[str, str],
    altloc_fills: Optional[Dict[str, Dict[str, str]]] = None,
    cofactor_edits: Dict[str, List[ChemistryEdit]],
    flexible: Sequence[str] = (),
    components: Optional[Dict[str, CcdComponent]] = None,
    metal_types: Optional[Dict[Tuple[str, str], str]] = None,
) -> PreparedReceptor:
    """Build the receptor PDBQT from ``pdb_path`` and the decisions given.

    ``chains``: polymer chains to include (their ATOM records). ``keep_hetero``:
    HETATM residues to keep, by id; every other HETATM residue (waters, the
    ligand to dock, additives) is left out and listed in ``removed``.
    ``components``: CCD entries by residue name, to use instead of the cached
    dictionary (tests, or entries edited elsewhere).
    ``metal_types``: (residue, atom) -> AutoDock type for a metal, or "omit" to
    leave that metal out. Meeko types only Mg, Ca, Mn, Fe and Zn; any other
    metal is typed as its element when Vina accepts that type (Vina 1.2.7
    accepts Cu, Co, Ni and W but not Mo or V), and otherwise refused until a
    type or "omit" is given here.
    """
    from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer

    atoms = read_atoms(pdb_path)
    keep = set(keep_hetero)
    removed: "OrderedDict[str, str]" = OrderedDict()
    selected: List[AtomLine] = []
    for atom in atoms:
        if atom.record == "ATOM" and atom.chain in chains:
            selected.append(atom)
        elif atom.record == "HETATM" and atom.res_id in keep:
            selected.append(atom)
        else:
            reason = ("water" if atom.resname in _WATER else
                      "chain not selected" if atom.record == "ATOM" else "HETATM not kept")
            removed.setdefault(f"{atom.res_id} {atom.resname}", reason)
    missing_kept = keep - {a.res_id for a in selected}
    metal_types = dict(metal_types or {})
    omitted = [a for a in selected if metal_types.get((a.res_id, a.name)) == "omit"]
    for atom in omitted:
        removed.setdefault(f"{atom.res_id} {atom.resname} atom {atom.name}", "metal left out by choice")
    selected = [a for a in selected if a not in omitted]
    if missing_kept:
        raise LookupError(f"Residues to keep not found: {', '.join(sorted(missing_kept))}.")

    # alternate locations: every residue that has them needs a choice
    by_residue: "OrderedDict[str, List[AtomLine]]" = OrderedDict()
    for atom in selected:
        by_residue.setdefault(atom.res_id, []).append(atom)
    unresolved = {}
    for res_id, residue_atoms in by_residue.items():
        codes = sorted({a.altloc for a in residue_atoms if a.altloc})
        if codes and res_id not in altlocs:
            unresolved[res_id] = codes
        elif codes and altlocs[res_id] not in codes:
            raise ValueError(f"{res_id} has alternate locations {', '.join(codes)}, not {altlocs[res_id]!r}.")
    if unresolved:
        listing = "; ".join(f"{r} ({', '.join(c)})" for r, c in unresolved.items())
        raise ValueError(f"Choose an alternate location for: {listing}.")
    for res_id in list(by_residue):
        by_residue[res_id] = [a for a in by_residue[res_id] if chosen_alternate(a, altlocs, altloc_fills)]

    undecided_his = [res_id for res_id, residue_atoms in by_residue.items()
                     if residue_atoms[0].resname == "HIS" and res_id not in protonation]
    if undecided_his:
        raise ValueError("Choose HID, HIE or HIP for every histidine (Meeko's own choice is a silent HIE, "
                         "or an error at a chain break): " + ", ".join(undecided_his) + ".")
    defaults = _meeko_default_keys()
    incomplete = _incomplete_standard_residues(by_residue)
    if incomplete:
        raise ValueError("Residues missing heavy atoms (rebuild them with the Structure Fixer, or leave "
                         "their chains out): " + "; ".join(incomplete) + ".")
    components = dict(components or {})
    notes: List[str] = []
    templates = rt.TemplateSet()
    metals_report: List[Dict[str, object]] = []
    residue_names = {res_id: residue_atoms[0].resname for res_id, residue_atoms in by_residue.items()}
    renamed: Dict[str, Tuple[str, str]] = {}     # res_id -> (new resname, new icode) for written lines
    atom_renames: Dict[Tuple[str, str], str] = {}  # (res_id, atom) -> new atom name
    split_out: Dict[Tuple[str, str], str] = {}     # (res_id, atom) -> new res_id

    # metals: every metal atom needs a charge; split out of multi-atom residues
    missing_charges = []
    for res_id, residue_atoms in by_residue.items():
        metal_atoms = [a for a in residue_atoms if a.element and is_metal_symbol(a.element)]
        if not metal_atoms:
            continue
        for metal in metal_atoms:
            if (res_id, metal.name) not in metal_charges:
                missing_charges.append(f"{metal.name} of {residue_names[res_id]} {res_id}")
    if missing_charges:
        raise ValueError("Choose a formal charge for every metal: " + "; ".join(missing_charges) + ".")

    new_residues: "OrderedDict[str, List[AtomLine]]" = OrderedDict()
    for res_id, residue_atoms in by_residue.items():
        metal_atoms = [a for a in residue_atoms if a.element and is_metal_symbol(a.element)]
        if not metal_atoms:
            continue
        taken_icodes = {a.icode for a in atoms if a.chain == residue_atoms[0].chain and a.resseq == residue_atoms[0].resseq}
        for metal in metal_atoms:
            charge = metal_charges[(res_id, metal.name)]
            ion_name = rt.ion_residue_name(metal.element, charge)
            templates.merge(rt.ion_template(metal.element, charge))
            if len(residue_atoms) == 1:
                renamed[res_id] = (ion_name, metal.icode)
                atom_renames[(res_id, metal.name)] = metal.element.upper()
                ion_id = res_id
            else:
                code = _free_icode(taken_icodes)
                taken_icodes.add(code)
                ion_id = f"{metal.chain}:{metal.resseq}{code}"
                split_out[(res_id, metal.name)] = ion_id
                atom_renames[(res_id, metal.name)] = metal.element.upper()
                new_residues[ion_id] = [metal]
                renamed[ion_id] = (ion_name, code)
            metals_report.append({"residue": f"{residue_names[res_id]} {res_id}", "atom": metal.name,
                                  "ion_residue": ion_id, "ion_name": ion_name, "charge": charge})
        by_residue[res_id] = [a for a in residue_atoms if a not in metal_atoms] or residue_atoms

    # cofactors: templates from checked CCD chemistry, one part per fragment left after the metals go
    cofactor_parts: Dict[str, List[Tuple[str, Chem.Mol]]] = {}     # resname -> [(template key, mol)]
    for res_id, residue_atoms in list(by_residue.items()):
        resname = residue_names[res_id]
        if res_id in renamed or resname in defaults or resname in cofactor_parts:
            continue
        component = components.get(resname) or fetch_component(resname)
        mol, problems = apply_edits(component.mol, cofactor_edits.get(resname, []))
        if problems:
            raise ValueError(f"CCD {resname} still has problems after the edits given: "
                             + "; ".join(p.describe() for p in problems))
        mol, _ = rt.split_metals(mol)
        expected = {a.GetProp("name") for a in mol.GetAtoms() if a.GetAtomicNum() > 1}
        for other_id, other_atoms in by_residue.items():
            if residue_names[other_id] != resname or other_id in renamed:
                continue
            present = {a.name for a in other_atoms if a.element not in ("H", "D")}
            if present != expected:
                parts = []
                if expected - present:
                    parts.append(f"missing {', '.join(sorted(expected - present))}")
                if present - expected:
                    parts.append(f"not in CCD {resname}: {', '.join(sorted(present - expected))}")
                raise ValueError(f"{resname} {other_id} does not match its CCD entry ({'; '.join(parts)}).")
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(fragments) == 1:
            cofactor_parts[resname] = [(resname, mol)]
        else:
            cofactor_parts[resname] = [(f"{resname}~{k + 1}", frag) for k, frag in enumerate(fragments)]
            notes.append(f"{resname}: {len(fragments)} separate pieces once its metals are split out; "
                         "each is its own residue.")

    set_template: Dict[str, str] = {}
    sources: Dict[str, str] = {}
    part_residue: Dict[Tuple[str, str], str] = {}       # (original res_id, atom) -> res_id of its part
    for res_id, residue_atoms in list(by_residue.items()):
        resname = residue_names[res_id]
        parts = cofactor_parts.get(resname)
        if res_id in renamed or not parts or len(parts) == 1:
            continue
        taken_icodes = {a.icode for a in atoms if a.chain == residue_atoms[0].chain and a.resseq == residue_atoms[0].resseq}
        taken_icodes |= {rid.split(":")[1][-1] for rid in new_residues if rid.startswith(f"{residue_atoms[0].chain}:{residue_atoms[0].resseq}")
                         and rid.split(":")[1][-1].isalpha()}
        kept_here = []
        for k, (key, frag) in enumerate(parts):
            names = {a.GetProp("name") for a in frag.GetAtoms()}
            members = [a for a in residue_atoms if a.name in names]
            if k == 0:
                target = res_id
                kept_here = members
            else:
                code = _free_icode(taken_icodes)
                taken_icodes.add(code)
                target = f"{residue_atoms[0].chain}:{residue_atoms[0].resseq}{code}"
                new_residues[target] = members
                renamed[target] = (resname, code)
                for atom in members:
                    split_out[(res_id, atom.name)] = target
            set_template[target] = key
            sources[target] = f"CCD chemistry, piece {k + 1} of {len(parts)} of {resname} {res_id}"
            for atom in members:
                part_residue[(res_id, atom.name)] = target
        by_residue[res_id] = kept_here

    # the PDB text Meeko reads
    def write_lines() -> str:
        out, serial = [], 0
        blocks = list(by_residue.items()) + list(new_residues.items())
        for res_id, residue_atoms in blocks:
            for atom in residue_atoms:
                serial += 1
                owner = split_out.get((atom.res_id, atom.name), res_id)
                resname, icode = renamed.get(owner, (atom.resname, atom.icode))
                name = atom_renames.get((atom.res_id, atom.name), atom.name)
                out.append(_format_atom(atom, serial, name, resname, icode))
        return "\n".join(out) + "\nEND\n"

    pdb_text = write_lines()

    def part_of(resname: str, atom: str) -> Tuple[str, Chem.Mol]:
        for key, frag in cofactor_parts[resname]:
            if atom in {a.GetProp("name") for a in frag.GetAtoms()}:
                return key, frag
        raise KeyError(f"{atom} is not an atom of {resname}.")

    # covalent links, as Meeko itself perceives them (metals never bonded)
    with _metals_unbonded():
        links = _find_links(pdb_text, defaults, set(cofactor_parts))
    link_atoms: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for link in links:
        for side in ("a", "b"):
            if link[f"{side}_resname"] in cofactor_parts:
                key, _ = part_of(link[f"{side}_resname"], link[f"{side}_atom"])
                link_atoms[key][link[f"{side}_atom"]] += 1

    variant_mols: Dict[str, Chem.Mol] = {}
    for resname, parts in cofactor_parts.items():
        for key, mol in parts:
            per_atom = dict(link_atoms.get(key, {}))
            stripped, removed_h = rt.remove_link_hydrogens(mol, per_atom) if per_atom else (mol, [])
            labels = {atom: f"{key}_{atom}" for atom in per_atom}
            templates.residue_templates[key] = rt.template_dict_from_mol(stripped, labels)
            variant_mols[key] = stripped
            if removed_h:
                notes.append(f"{key}: hydrogens {', '.join(removed_h)} removed where it bonds to "
                             "a neighbouring residue.")

    for link in links:
        for side, other in (("a", "b"), ("b", "a")):
            resname, atom, res_id = link[f"{side}_resname"], link[f"{side}_atom"], link[f"{side}_res"]
            if resname in cofactor_parts:
                continue
            base = protonation.get(res_id, resname)
            key = f"{base}~{atom}"
            if key not in templates.residue_templates:
                templates.residue_templates[key] = rt.linked_variant(base, atom, f"{base}_{atom}", key)
                variant_mols[key] = _template_mol(templates.residue_templates[key])
            set_template[res_id] = key
            sources[res_id] = f"covalent link {atom}-{link[f'{other}_atom']} to {link[f'{other}_resname']} {link[f'{other}_res']}"

    def template_side(resname: str, atom: str, res_id: str) -> Tuple[str, str]:
        """(template key, link label) of one side of a link."""
        if resname in cofactor_parts:
            key, _ = part_of(resname, atom)
            return key, f"{key}_{atom}"
        key = set_template[res_id]
        return key, f"{key.split('~')[0]}_{atom}"

    for link in links:
        for side, other in (("a", "b"), ("b", "a")):
            own_key, label = template_side(link[f"{side}_resname"], link[f"{side}_atom"], link[f"{side}_res"])
            other_key, _ = template_side(link[f"{other}_resname"], link[f"{other}_atom"], link[f"{other}_res"])
            own_mol = variant_mols.get(own_key) or _template_mol(templates.residue_templates[own_key])
            other_mol = variant_mols.get(other_key) or _template_mol(templates.residue_templates[other_key])
            templates.padders[label] = rt.link_padder(own_mol, link[f"{side}_atom"],
                                                      rt.neighbour_smarts(other_mol, link[f"{other}_atom"]))

    # protonation and termini
    for res_id, key in protonation.items():
        if res_id in set_template:
            continue
        set_template[res_id] = key
        sources[res_id] = "protonation choice"
    for res_id, choice in termini.items():
        if choice not in TERMINUS_CHOICES:
            raise ValueError(f"{res_id}: terminus must be one of {', '.join(TERMINUS_CHOICES)}, got {choice!r}.")
        if choice == "chain break":
            sources.setdefault(res_id, "chain break (neutral end)")
            continue
        prefix = "N" if choice == "N-terminus" else "C"
        base = set_template.get(res_id, residue_names.get(res_id))
        if base is None:
            raise LookupError(f"Terminus given for {res_id}, which is not in the receptor.")
        set_template[res_id] = prefix + base
        sources[res_id] = f"{choice}" + (f" + {sources[res_id]}" if res_id in sources else "")

    chem = templates.to_meeko()
    missing_templates = sorted({k for k in set_template.values() if k not in chem.residue_templates})
    if missing_templates:
        raise ValueError(f"No Meeko template named {', '.join(missing_templates)}.")
    mk = MoleculePreparation()
    with _metals_unbonded():
        polymer = Polymer.from_pdb_string(pdb_text, chem, mk, set_template=set_template)
    ignored = polymer.get_ignored_monomers()
    if ignored:
        raise RuntimeError("Meeko matched no template for: " + ", ".join(
            f"{k} ({m.input_resname})" for k, m in ignored.items()))
    notes.extend(_type_metals(polymer, metals_report, metal_types))
    for res_id in flexible:
        polymer.flexibilize_sidechain(res_id, mk)
    rigid, flex = PDBQTWriterLegacy.write_from_polymer(polymer)

    table = []
    for res_id, monomer in polymer.get_valid_monomers().items():
        source = sources.get(res_id)
        if source is None:
            source = "ion, charge chosen" if res_id in renamed else (
                "CCD chemistry" + (" + edits" if cofactor_edits.get(monomer.input_resname) else "")
                if monomer.input_resname in cofactor_parts else "Meeko template (unambiguous)"
                if monomer.residue_template_key == monomer.input_resname else "Meeko default")
        table.append(ResidueDecision(res_id, monomer.input_resname, monomer.residue_template_key, source))
    defaults_used = [d for d in table if d.source == "Meeko default"]
    if defaults_used:
        notes.append(f"{len(defaults_used)} residues took Meeko's own choice among several templates "
                     "(no decision was given for them): "
                     + ", ".join(f"{d.res_id} {d.input_resname}->{d.template}" for d in defaults_used[:12])
                     + (" ..." if len(defaults_used) > 12 else "") + ".")
    return PreparedReceptor(
        rigid_pdbqt=rigid, flex_pdbqt="".join(flex.values()) if flex else "",
        residues=table, metals=metals_report,
        links=[{k: v for k, v in link.items()} for link in links],
        removed=[f"{k} ({v})" for k, v in removed.items()], notes=notes)


_VINA_TYPE_CACHE: Dict[str, bool] = {}


def vina_accepts_type(atom_type: str) -> bool:
    """Ask Vina itself whether it reads ``atom_type`` in a receptor PDBQT (cached)."""
    if atom_type not in _VINA_TYPE_CACHE:
        import tempfile
        from vina import Vina
        with tempfile.NamedTemporaryFile("w", suffix=".pdbqt", delete=False) as handle:
            handle.write("ATOM      1  X   UNK A   1       0.000   0.000   0.000  1.00  0.00     0.000 "
                         f"{atom_type:<2}\n")
            path = handle.name
        try:
            Vina(sf_name="vina", seed=1, verbosity=0).set_receptor(path)
            _VINA_TYPE_CACHE[atom_type] = True
        except Exception:
            _VINA_TYPE_CACHE[atom_type] = False
        finally:
            os.unlink(path)
    return _VINA_TYPE_CACHE[atom_type]


def _type_metals(polymer, metals_report: List[Dict[str, object]],
                 metal_types: Dict[Tuple[str, str], str]) -> List[str]:
    """Give every metal Meeko left untyped a type Vina reads; refuse, by name, the ones it cannot."""
    origin = {m["ion_residue"]: (m["residue"].split()[-1], m["atom"]) for m in metals_report}
    notes, refused = [], []
    for res_id, monomer in polymer.get_valid_monomers().items():
        for atom in monomer.molsetup.atoms:
            if atom.atom_type is not None or atom.is_ignore or atom.is_dummy:
                continue
            element = Chem.GetPeriodicTable().GetElementSymbol(atom.atomic_num)
            if not is_metal_symbol(element):
                refused.append(f"{res_id} atom {atom.index + 1} ({element}): Meeko assigned no type")
                continue
            key = origin.get(res_id, (res_id, element.upper()))
            chosen = metal_types.get(key)
            if chosen is None and vina_accepts_type(element):
                chosen = element
                notes.append(f"{element} at {res_id} typed {element} (Meeko types only Mg, Ca, Mn, Fe, Zn; "
                             "Vina accepts this one).")
            elif chosen is not None:
                if not vina_accepts_type(chosen):
                    refused.append(f"{element} at {res_id}: Vina does not accept the chosen type {chosen!r}")
                    continue
                notes.append(f"{element} at {res_id} typed {chosen} by choice.")
            else:
                refused.append(f"{element} at {res_id} (from {key[0]} atom {key[1]}): Vina has no {element} type")
                continue
            atom.atom_type = chosen
    if refused:
        raise ValueError("Metals without an AutoDock type Vina accepts: " + "; ".join(refused)
                         + ". Give each a type (metal_types) or leave it out ('omit').")
    return notes


def _incomplete_standard_residues(by_residue: "OrderedDict[str, List[AtomLine]]") -> List[str]:
    """Standard residues missing heavy atoms of their Meeko template, as "A:4 ARG: CG, CD, ...".

    Meeko would refuse them with its template statistics; this names the atoms.
    The comparison uses the residue's internal template (no OXT), so a C-terminal
    OXT is never reported missing.
    """
    from meeko import ResidueChemTemplates
    chem = ResidueChemTemplates.create_from_defaults()
    found = []
    for res_id, residue_atoms in by_residue.items():
        resname = residue_atoms[0].resname
        template = chem.residue_templates.get(resname)
        if template is None or residue_atoms[0].record != "ATOM":
            continue
        heavy = [name for name, atom in zip(template.atom_names, template.mol.GetAtoms()) if atom.GetAtomicNum() > 1]
        present = {a.name for a in residue_atoms}
        missing = [name for name in heavy if name not in present]
        if missing:
            found.append(f"{res_id} {resname}: {', '.join(missing)}")
    return found


@contextmanager
def _metals_unbonded():
    """Every metal's covalent radius 0 in Meeko's table while it perceives bonds.

    Meeko already gives Ca, Cu, Fe, K, Mg, Mn, Na and Zn radius 0, so they are
    never bonded, but not Co, Ni, Mo, W, V and the rest, which it would bond to
    their ligands. Here every metal is an ion of its own, so all are treated
    as Meeko treats iron and zinc. The table is restored afterwards.
    """
    from meeko.utils.covalent_radius_table import covalent_radius
    saved = {}
    for symbol, radius in covalent_radius.items():
        if symbol.strip() and is_metal_symbol(symbol) and radius:
            saved[symbol] = radius
            covalent_radius[symbol] = 0.0
    try:
        yield
    finally:
        covalent_radius.update(saved)


def _template_mol(template: dict) -> Chem.Mol:
    from meeko import ResidueTemplate
    mol = ResidueTemplate(template["smiles"], None, template["atom_name"]).mol
    for atom, name in zip(mol.GetAtoms(), template["atom_name"]):
        atom.SetProp("name", name)
    return mol


_PEPTIDE = {("C", "N"), ("N", "C"), ("O3'", "P"), ("P", "O3'")}


def _find_links(pdb_text: str, defaults: set, cofactors: set) -> List[Dict[str, object]]:
    """Bonds between residues that involve a cofactor, or a standard residue bonded to anything
    other than its backbone neighbours and disulfide partner, as Meeko's own perception finds them."""
    from meeko.polymer import Polymer, find_inter_mols_bonds
    raw = Polymer._pdb_to_residue_mols(pdb_text)
    mols = {k: (v[0], v[1]) for k, v in raw.items() if v[0] is not None}
    bonds = find_inter_mols_bonds(mols)
    links = []
    for (res_a, res_b), pairs in bonds.items():
        mol_a, name_a = mols[res_a]
        mol_b, name_b = mols[res_b]
        for i, j in pairs:
            atom_a = mol_a.GetAtomWithIdx(i).GetPDBResidueInfo().GetName().strip()
            atom_b = mol_b.GetAtomWithIdx(j).GetPDBResidueInfo().GetName().strip()
            if name_a not in cofactors and name_b not in cofactors:
                if (atom_a, atom_b) in _PEPTIDE or (atom_a == "SG" and atom_b == "SG"):
                    continue
            links.append({"a_res": res_a, "a_resname": name_a, "a_atom": atom_a,
                          "b_res": res_b, "b_resname": name_b, "b_atom": atom_b})
    return links
