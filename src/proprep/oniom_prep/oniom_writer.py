"""
ONIOM Input Writer - Gaussian Input File Generation

Writes Gaussian ONIOM input files (.com) from ONIOMSetup data structures.
Handles 2-layer and 3-layer ONIOM calculations with proper formatting including
force field parameter sections for non-standard residues.

© 2024 ProPrep Developer. All rights reserved.
"""

import logging
from typing import Dict, List, Optional, Tuple, TextIO, Set
from datetime import datetime
from pathlib import Path
from rich.console import Console

from proprep.structure_prep.comprehensive_redox_detector import RedoxSite, RedoxSiteAtom
from proprep.forcefield_prep.mcpb.ff_parameter_reader import (
    FFParameterReader,
    BondParam,
    AngleParam,
    NonbondedParam,
    DihedralParam
)
from .data_structures import (
    ONIOMSetup,
    ONIOMLayer,
    FreezeFlag,
    LayerAssignment,
    LinkAtom,
    calculate_layer_statistics
)


logger = logging.getLogger(__name__)


# Link atoms get a single dedicated type "HX". Its bond/angle/VDW parameters
# are written explicitly (see _write_link_atom_parameters) rather than reusing
# a real ff14SB H type, which is brittle across force fields.
#
# This is the standard *scaled* hydrogen link atom: Gaussian places the cap at
# a normal C–H distance (default scale factor ~0.7239) and we zero its charge,
# so the cap's parameters are ordinary C–H / N–H values copied from the FF's
# natural H type for each QM-parent type (PROXY_H_TYPE_FOR_LINK). This is NOT
# t2ONIOM's heavy-atom scheme (which leaves the link atom at the heavy-atom
# position with C–C parameters — only consistent with an *unscaled* cut). The
# scaled-H convention matches Gaussian's default and the Schlegel-group TAO
# toolkit (link H, 0.7239 scale, zero charge).
LINK_ATOM_TYPE: str = "HX"

# For each QM-parent atom type, which natural H type's parameters HX *copies*
# when bonded to it (so the cap inherits realistic stiffness/geometry from the
# FF). The carbonyl C ("C", the formamide cap from a single-bond C-cut) has no
# C–H bond in ff14SB — those terms are supplied from GAFF (see below).
PROXY_H_TYPE_FOR_LINK: Dict[str, str] = {
    # Amide/amine nitrogens (N-cut / primary-amide cap).
    "N": "H", "N3": "H", "NA": "H", "NB": "H", "N2": "H", "N*": "H",
    # Protein alpha carbon (sp3 C with adjacent N).
    "CX": "H1",
    # Aliphatic sp3 carbon (sidechain CB and beyond).
    "CT": "HC",
    # Aromatic carbons.
    "CA": "HA", "CB": "HA", "CC": "HA", "CK": "HA", "CM": "HA",
    "CN": "HA", "CQ": "HA", "CR": "HA", "CV": "HA", "CW": "HA",
}

# Cap MM terms ff14SB cannot supply, sourced from GAFF (the parm99 -> GAFF
# "parmlookup" convention used by t2ONIOM/TAO). The carbonyl C ("C") gets an
# aldehyde/formyl C–H for the formamide cap produced by a single-bond C-cut
# (Cα–C). Values are gaff.dat entries: c-h4 (bond), h4-c-o / h4-c2-n (angles).
_SUPPLIED_CAP_BOND: Dict[str, tuple] = {
    "C": (310.7, 1.1121),           # GAFF  c -h4    carbonyl/formyl C–H
}
_SUPPLIED_CAP_ANGLE: Dict[tuple, tuple] = {
    ("O", "C"): (54.2, 120.70),     # GAFF  h4-c -o   carbonyl O – C – cap H
    ("N", "C"): (50.7, 113.44),     # GAFF  h4-c2-n   amide N    – C – cap H
}


def _proxy_h_type_for_link(qm_type: str) -> str:
    """Which natural H type's FF params should HX inherit when bonded to `qm_type`?"""
    return PROXY_H_TYPE_FOR_LINK.get(qm_type, "HC")


class ONIOMWriter:
    """
    Writes Gaussian ONIOM input files from ONIOMSetup objects.

    Gaussian ONIOM Format:
    1. Route section (ONIOM keywords with AMBER=SoftOnly for custom parameters)
    2. Title
    3. Charge and multiplicity
    4. Atom section (element, freeze, x, y, z, layer, atom_type, charge, [link_info])
    5. [Optional] Connectivity section
    6. Two blank lines
    7. Force field parameter section (VDW, bonds, angles, dihedrals)
    8. Two blank lines at end
    """

    def __init__(
        self,
        oniom_setup: ONIOMSetup,
        ff_parameter_reader: Optional[FFParameterReader] = None,
        logger_instance=None,
        console=None
    ):
        """
        Initialize ONIOM writer.

        Args:
            oniom_setup: Complete ONIOM setup configuration
            ff_parameter_reader: FFParameterReader for force field parameters
            logger_instance: Optional logger
            console: Optional Rich Console
        """
        self.oniom_setup = oniom_setup
        self.ff_reader = ff_parameter_reader
        self.logger = logger_instance or logger
        self.console = console or Console()

        # Atom ordering: atom_idx -> 1-based Gaussian index
        self.idx_to_gaussian: Dict[int, int] = {}
        self._build_atom_ordering()

        # Build remapping for atom types that start with digits (breaks Gaussian parser)
        self.type_remap: Dict[str, str] = {}
        self._build_type_remap()

    def write_diagnostic_file(
        self,
        output_path: str,
        title: Optional[str] = None,
    ) -> bool:
        """
        Write a diagnostic Gaussian ONIOM input file.

        Uses ONIOM=OnlyInputFiles to make Gaussian write separate input files
        for each sub-calculation (real-low, model-high, model-low), so you
        can inspect that layer definitions, charges, and atom counts are correct.

        OnlyInputFiles writes the sub-calculation files and stops without
        running the actual calculation (unlike InputFiles which still runs).

        Args:
            output_path: Path to output .com file
            title: Optional title line

        Returns:
            True if successful, False otherwise
        """
        orig_keywords = self.oniom_setup.additional_keywords
        diag_title = title or "ONIOM Diagnostic - OnlyInputFiles (verify layer assignments)"

        # Add OnlyInputFiles as an ONIOM option (will be parsed alongside ScaleCharge)
        extra = "OnlyInputFiles"
        if orig_keywords:
            self.oniom_setup.additional_keywords = f"{orig_keywords} {extra}"
        else:
            self.oniom_setup.additional_keywords = extra

        try:
            success = self._write_file_impl(
                output_path=output_path,
                title=diag_title,
                include_comments=True,
                include_ff_parameters=True,
            )
            return success
        finally:
            # Restore original keywords
            self.oniom_setup.additional_keywords = orig_keywords

    def write_input_file(
        self,
        output_path: str,
        title: Optional[str] = None,
        include_comments: bool = True,
        include_ff_parameters: bool = True
    ) -> bool:
        """
        Write complete Gaussian ONIOM input file.

        Args:
            output_path: Path to output .com file
            title: Optional title line (default: auto-generated)
            include_comments: Include informative comments in file
            include_ff_parameters: Include force field parameter section

        Returns:
            True if successful, False otherwise
        """
        return self._write_file_impl(output_path, title, include_comments, include_ff_parameters)

    def _write_file_impl(
        self,
        output_path: str,
        title: Optional[str] = None,
        include_comments: bool = True,
        include_ff_parameters: bool = True,
        extra_route_keywords: Optional[List[str]] = None,
    ) -> bool:
        """Internal implementation for writing Gaussian input files."""
        try:
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)

            with open(output_file, 'w') as f:
                # Section 1: Route section
                self._write_route_section(f, include_comments, include_ff_parameters,
                                          extra_route_keywords=extra_route_keywords)

                # Section 2: Title
                self._write_title_section(f, title, include_comments)

                # Section 3: Charge and multiplicity
                self._write_charge_multiplicity_section(f, include_comments)

                # Section 4: Atom section
                self._write_atom_section(f, include_comments)

                # Section 5: Connectivity (optional)
                if self.oniom_setup.use_explicit_connectivity:
                    self._write_connectivity_section(f, include_comments)

                # Section 6: Force field parameters (if requested and available)
                if include_ff_parameters and self.ff_reader is not None:
                    # One blank line before parameter section
                    f.write("\n")
                    self._write_force_field_parameters(f, include_comments)

                # Section 7: Final blank lines
                f.write("\n\n")

            self.logger.info(f"Successfully wrote ONIOM input file: {output_path}")
            return True

        except Exception as e:
            self.logger.error(f"Failed to write ONIOM input file: {e}")
            return False

    def _build_atom_ordering(self):
        """
        Build atom_idx -> 1-based Gaussian index mapping for all atoms.

        Link atoms are NOT separate entries; they appear on the same line
        as their QM parent atom. So we only index real atoms.
        """
        for gaussian_idx, atom_idx in enumerate(
            self.oniom_setup.layer_assignments.keys(), start=1
        ):
            self.idx_to_gaussian[atom_idx] = gaussian_idx

        self.logger.debug(f"Built atom ordering: {len(self.idx_to_gaussian)} total atoms")

    def _build_type_remap(self):
        """
        Build remapping for atom types that need Gaussian-compatible names.

        Handles two issues:
        1. Digit-prefix types: Gaussian's parser interprets tokens starting with
           digits as numbers, so GAFF2 types like '2C' are remapped to 'C2'.
        2. Case collisions: AMBER is case-sensitive (CT ≠ ct), but Gaussian is
           case-insensitive. When both ff14SB (uppercase) and GAFF2 (lowercase)
           types are present, colliding lowercase types are renamed (e.g., hc → hc0).
        """
        all_types = set(self.oniom_setup.atom_types.values())

        # --- Pass 0: Types containing dashes (e.g., 'Cl-' for chloride) ---
        # Dashes are the delimiter in Gaussian's Element-Type-Charge format,
        # so atom types containing dashes must be renamed.
        dash_types = [t for t in all_types if '-' in t]
        for orig_type in dash_types:
            candidate = orig_type.replace('-', '')
            if not candidate:
                candidate = orig_type.replace('-', 'x')
            if candidate in all_types or candidate in self.type_remap.values():
                candidate = candidate + "0"
            self.type_remap[orig_type] = candidate
            all_types.add(candidate)

        # --- Pass 1: Digit-prefix remapping ---
        digit_types = [t for t in all_types if t and t[0].isdigit()]
        for orig_type in digit_types:
            # Strategy 1: move leading digits to end (e.g., '2C' -> 'C2')
            i = 0
            while i < len(orig_type) and orig_type[i].isdigit():
                i += 1
            candidate = orig_type[i:] + orig_type[:i]

            if candidate not in all_types and candidate not in self.type_remap.values():
                self.type_remap[orig_type] = candidate
            else:
                # Strategy 2: swap first letter to 'Z' (e.g., 'C2' collision -> 'Z2')
                candidate = "Z" + orig_type[1:]
                suffix = 0
                while candidate in all_types or candidate in self.type_remap.values():
                    suffix += 1
                    candidate = f"Z{orig_type[1:]}{suffix}"
                self.type_remap[orig_type] = candidate

            all_types.add(self.type_remap[orig_type])

        # --- Pass 2: Case-insensitive collision resolution ---
        # Gaussian treats AMBER atom types case-insensitively, but AMBER uses
        # uppercase for protein types (CT, N, O from parm10/ff14SB) and lowercase
        # for ligand types (ct, n, o from GAFF2/amberdyes). We must rename the
        # lowercase variants to avoid "Multiple matching entries" errors.
        from collections import defaultdict
        lower_groups = defaultdict(list)
        for t in all_types:
            lower_groups[t.lower()].append(t)

        # Track all output names (case-insensitive) to avoid creating new collisions
        output_names_lower = set()
        for t in all_types:
            remapped = self.type_remap.get(t, t)
            output_names_lower.add(remapped.lower())

        case_remaps = {}
        for lower_key, group in sorted(lower_groups.items()):
            if len(group) < 2:
                continue  # no collision

            # Sort: uppercase/mixed-case types first (protein/parm10 convention),
            # all-lowercase types second (GAFF2). Uppercase types keep their names.
            group.sort(key=lambda t: (t == t.lower(), t))

            # First type keeps its name; rename the rest (typically lowercase GAFF2)
            for orig_type in group[1:]:
                if orig_type in self.type_remap:
                    continue  # already remapped by digit-prefix logic

                # Generate unique name: append 0, 1, 2... until no collision
                suffix = 0
                candidate = f"{orig_type}{suffix}"
                while (candidate in all_types
                       or candidate in self.type_remap.values()
                       or candidate.lower() in output_names_lower):
                    suffix += 1
                    candidate = f"{orig_type}{suffix}"

                self.type_remap[orig_type] = candidate
                case_remaps[orig_type] = candidate
                all_types.add(candidate)
                output_names_lower.add(candidate.lower())

        # --- Logging ---
        if self.type_remap:
            digit_count = len(digit_types)
            case_count = len(case_remaps)

            if digit_count > 0:
                self.logger.info(
                    f"Remapped {digit_count} digit-prefix atom types: "
                    + ", ".join(f"{k}->{v}" for k, v in sorted(self.type_remap.items())
                               if k and k[0].isdigit())
                )

            if case_count > 0:
                self.logger.info(
                    f"Resolved {case_count} case-insensitive atom type collision(s): "
                    + ", ".join(f"{k}->{v}" for k, v in sorted(case_remaps.items()))
                )
                self.console.print(
                    f"\n  [yellow]⚠ Resolved {case_count} AMBER/GAFF2 atom type "
                    f"case collision(s) for Gaussian:[/yellow]"
                )
                self.console.print(
                    f"    {', '.join(f'{k}→{v}' for k, v in sorted(case_remaps.items()))}",
                    highlight=False
                )
                self.console.print(
                    f"    [grey50]Gaussian is case-insensitive; GAFF2 types renamed "
                    f"to avoid conflicts with ff14SB/parm10 types.[/grey50]"
                )

    def _remap_type(self, atom_type: str) -> str:
        """Apply atom type remapping for Gaussian compatibility."""
        return self.type_remap.get(atom_type, atom_type)

    def _write_route_section(
        self,
        f: TextIO,
        include_comments: bool,
        include_ff_parameters: bool,
        extra_route_keywords: Optional[List[str]] = None,
    ):
        """
        Write Gaussian route section with method keywords.

        Format:
        # method/basis ONIOM(high:medium:low)=options Amber=SoftFirst job_type additional_keywords
        """
        # Note: Gaussian does not support ! comments in input files
        # Comments have been removed to ensure compatibility

        # Memory and processors
        f.write(f"%mem={self.oniom_setup.memory_gb}GB\n")
        f.write(f"%nprocshared={self.oniom_setup.n_processors}\n")

        # Main route line
        route_line = "#P"  # P for verbose output

        # ONIOM method specification
        # SoftOnly goes on the AMBER keyword inside the ONIOM parentheses,
        # NOT as a separate Amber=SoftOnly keyword (which Gaussian treats as
        # a duplicate and may cause unexpected behavior).
        mm_method = self.oniom_setup.mm_forcefield
        if include_ff_parameters and self.ff_reader is not None:
            mm_method = f"{mm_method}=SoftOnly"

        if self.oniom_setup.n_layers == 2:
            oniom_methods = (
                f"{self.oniom_setup.qm_functional}/{self.oniom_setup.qm_basis_set}"
                f":{mm_method}"
            )
        else:
            oniom_methods = (
                f"{self.oniom_setup.qm_functional}/{self.oniom_setup.qm_basis_set}"
                f":{self.oniom_setup.medium_method}"
                f":{mm_method}"
            )

        # Build ONIOM options list — some keywords must be ONIOM options, not standalone
        # ONIOM options: ScaleCharge, InputFiles, SValue, etc.
        oniom_option_prefixes = ("SCALECHARGE", "ONLYINPUTFILES", "INPUTFILES", "SVALUE")
        oniom_options = ["EmbedCharge"]
        additional = self.oniom_setup.additional_keywords
        remaining_keywords = []
        if additional:
            for kw in additional.split():
                if kw.upper().startswith(oniom_option_prefixes):
                    oniom_options.append(kw)
                else:
                    remaining_keywords.append(kw)

        # Add any extra route keywords (e.g., "Test" for diagnostic runs)
        if extra_route_keywords:
            remaining_keywords.extend(extra_route_keywords)

        # Format: ONIOM(methods)=(option1,option2,...)
        route_line += f" ONIOM({oniom_methods})=({','.join(oniom_options)})"

        # Job type
        route_line += f" {self.oniom_setup.job_type}"

        # Connectivity specification
        if self.oniom_setup.use_explicit_connectivity:
            # Using explicit connectivity table - tell Gaussian to read it
            route_line += " Geom=Connect"
        else:
            # Let Gaussian auto-generate connectivity from coordinates
            route_line += " Geom=Connectivity"

        # Additional keywords (ScaleCharge already moved to ONIOM options above)
        if remaining_keywords:
            route_line += f" {' '.join(remaining_keywords)}"

        f.write(route_line + "\n")
        f.write("\n")

    def _write_title_section(self, f: TextIO, title: Optional[str], include_comments: bool):
        """
        Write title section.

        Args:
            f: File handle
            title: Optional custom title
            include_comments: Include metadata comments
        """
        if title is None:
            # Auto-generate title
            stats = calculate_layer_statistics(self.oniom_setup)
            title = (
                f"ONIOM {self.oniom_setup.n_layers}-layer calculation: "
                f"HIGH={stats.n_high} atoms, "
                f"MEDIUM={stats.n_medium} atoms, "
                f"LOW={stats.n_low} atoms, "
                f"LINK={stats.n_link} atoms"
            )

        f.write(title + "\n")
        f.write("\n")

    def _write_charge_multiplicity_section(self, f: TextIO, include_comments: bool):
        """
        Write charge and multiplicity line(s).

        Gaussian ONIOM format for 2-layer (from gaussian.com/oniom):
          chrgreal-low  spinreal-low  chrgmodel-high  spinmodel-high

        Gaussian ONIOM format for 3-layer:
          cRealL sRealL  cIntM sIntM  cIntL sIntL  cModH sModH

        Where real = entire system, model = QM region, int = intermediate.
        """
        # Total system charge = round(sum of all atomic partial charges)
        total_charge = round(sum(self.oniom_setup.charges.values()))

        # Determine real system multiplicity from electron count parity.
        # Total electrons = sum(atomic_numbers) - total_charge.
        # Odd electrons require even multiplicity (minimum: 2).
        parm = self.oniom_setup.parm
        if parm:
            total_electrons = sum(a.atomic_number for a in parm.atoms) - total_charge
            if total_electrons % 2 == 0:
                real_mult = 1  # even electrons → singlet
            else:
                real_mult = 2  # odd electrons → doublet
                self.logger.info(
                    f"Real system has {total_electrons} electrons (odd): "
                    f"using multiplicity 2"
                )
        else:
            real_mult = 1

        if self.oniom_setup.n_layers == 2:
            charge_mult_line = (
                f"{total_charge} {real_mult} "
                f"{self.oniom_setup.qm_charge} {self.oniom_setup.qm_multiplicity}"
            )
        else:
            charge_mult_line = (
                f"{total_charge} {real_mult} "
                f"{self.oniom_setup.medium_charge} {self.oniom_setup.medium_multiplicity} "
                f"{self.oniom_setup.medium_charge} {self.oniom_setup.medium_multiplicity} "
                f"{self.oniom_setup.qm_charge} {self.oniom_setup.qm_multiplicity}"
            )

        f.write(charge_mult_line + "\n")

    def _write_atom_section(self, f: TextIO, include_comments: bool):
        """
        Write atom coordinate section with ONIOM formatting.

        ONIOM Format:
        Element-AtomType-Charge freeze x y z Layer [LinkAtom BondedToIndex]

        Example without link atom:
        C-CA--0.25  0  -4.748381  0.578498  -0.778569  H

        Example with link atom:
        C-CA--0.25  0  -4.748381  0.578498  -0.778569  H  H-HC-0.1  3

        Where:
        - Element-AtomType-Charge is the atom identifier (e.g., C-CA--0.25)
        - freeze is 0 (active) or -1 (frozen)
        - x y z are coordinates
        - Layer is H (high), M (medium), or L (low)
        - LinkAtom is optional link atom spec (e.g., H-HC-0.1)
        - BondedToIndex is the MM parent atom index (1-based)
        """
        # Atom section (comments removed for Gaussian compatibility)

        # Link atom specs go on the LOW-layer atom's line (the atom being
        # replaced by H in the model system), NOT on the HIGH-layer atom.
        # Keyed by mm_parent_idx (the LOW atom).
        link_atoms_by_mm_parent = {}
        for link_atom in self.oniom_setup.link_atoms:
            link_atoms_by_mm_parent.setdefault(link_atom.mm_parent_idx, []).append(link_atom)

        # Write all atoms in layer_assignments order
        for atom_idx, assignment in self.oniom_setup.layer_assignments.items():
            link_atoms = link_atoms_by_mm_parent.get(atom_idx, [])
            self._write_single_atom_from_assignment(f, atom_idx, assignment, link_atoms)

        f.write("\n")

    def _write_single_atom(self, f: TextIO, atom, is_link: bool = False):
        """
        Write a single real atom line (legacy method, rarely used).
        """
        # This method is kept for backward compatibility but the main path
        # uses _write_single_atom_from_assignment which operates on atom_idx.
        pass

    def _write_single_atom_from_assignment(
        self,
        f: TextIO,
        atom_idx: int,
        assignment: LayerAssignment,
        link_atoms: List[LinkAtom] = None,
    ):
        """
        Write a single real atom line from layer assignment.

        Args:
            f: File handle
            atom_idx: Parmed atom index (0-based)
            assignment: LayerAssignment for this atom
            link_atoms: Optional list of LinkAtom objects attached to this atom
        """
        if link_atoms is None:
            link_atoms = []

        atom_type = self._remap_type(
            self.oniom_setup.atom_types.get(atom_idx, assignment.element)
        )
        charge = self.oniom_setup.charges.get(atom_idx, 0.0)

        atom_identifier = f"{assignment.element}-{atom_type}-{charge:.6f}"
        freeze_value = assignment.freeze.value
        x, y, z = assignment.coords
        layer_str = assignment.layer.value

        atom_line = (
            f"{atom_identifier:20s} {freeze_value:2d} "
            f"{x:12.6f} {y:12.6f} {z:12.6f} "
            f"{layer_str:2s}"
        )

        if link_atoms:
            # Gaussian allows one link atom per line. The link atom spec
            # goes on the LOW atom's line; BondedToIndex references the
            # HIGH atom the link H is bonded to in the model system.
            link_atom = link_atoms[0]
            qm_parent_gaussian = self.idx_to_gaussian.get(link_atom.qm_parent_idx)
            if qm_parent_gaussian is not None:
                # Every cap H is type HX. Bond/angle/VDW params for HX paired
                # with each QM-parent type are emitted explicitly in
                # _write_link_atom_parameters().
                link_spec = f"H-{LINK_ATOM_TYPE}-0.0"
                atom_line += f"  {link_spec}  {qm_parent_gaussian}"

            if len(link_atoms) > 1:
                self.logger.warning(
                    f"Atom idx {atom_idx} has {len(link_atoms)} link atoms; "
                    f"only the first is written."
                )

        f.write(atom_line + "\n")

    def _write_connectivity_section(self, f: TextIO, include_comments: bool):
        """
        Write connectivity table (OPTIONAL - only if use_explicit_connectivity=True).

        Format:
        atom_index connected_index1 bond_order1 connected_index2 bond_order2 ...

        Example:
        1 2 1.0 5 1.0
        2 1 1.0 3 2.0

        Note: This is deferred to Phase 4. For now, we rely on Gaussian's
        automatic connectivity generation via Geom=Connectivity.
        """
        # Connectivity section (comments removed for Gaussian compatibility)

        if not self.oniom_setup.connectivity:
            self.logger.warning("Connectivity table requested but not provided")
            return

        for entry in self.oniom_setup.connectivity:
            conn_line = str(entry.atom_index)

            for bonded_index, bond_order in zip(entry.connected_indices, entry.bond_orders):
                conn_line += f" {bonded_index} {bond_order:.1f}"

            f.write(conn_line + "\n")

    def _write_force_field_parameters(self, f: TextIO, include_comments: bool):
        """
        Write force field parameter section for non-standard residues.

        This section is required when using Amber=SoftFirst or Amber=SoftOnly.
        Format follows Gaussian's parameter input specification:
        - VDW atom_type radius well_depth
        - HrmStr1 type1 type2 force_constant eq_length
        - HrmBnd1 type1 type2 type3 force_constant eq_angle
        - Torsion type1 type2 type3 type4 terms...

        Only includes parameters for atom types actually present in the structure.
        """
        if self.ff_reader is None:
            self.logger.warning("Force field parameters requested but no FFParameterReader provided")
            return

        # Force field parameter section (comments removed for Gaussian compatibility)

        # Collect all atom types used in the structure
        atom_types_used = self._get_atom_types_used()

        if not atom_types_used:
            self.logger.warning("No atom types found in structure")
            return

        self.logger.info(f"Writing force field parameters for {len(atom_types_used)} atom types")

        # NonBon header: defines nonbonded interaction model (required before VDW entries)
        # Format: NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.0
        #   3 = geometric combining rule for well depth
        #   1 = arithmetic combining rule for radius
        #   0.5 = 1-4 VDW scale factor
        #   -1.0 = signals Rmin/2 format (not sigma)
        f.write("NonBon 3 1 0 0 0.000 0.000 0.500 0.000 0.000 -1.000\n")

        # Write VDW parameters
        self._write_vdw_parameters(f, atom_types_used, include_comments)

        # Write bond parameters
        self._write_bond_parameters(f, atom_types_used, include_comments)

        # Write angle parameters
        self._write_angle_parameters(f, atom_types_used, include_comments)

        # Write dihedral parameters
        self._write_dihedral_parameters(f, atom_types_used, include_comments)

        # Write link-atom (HX) parameters by copying from the FF's natural H type
        # for each QM-parent atom type the caps actually attach to.
        self._write_link_atom_parameters(f, atom_types_used, include_comments)

    def _write_link_atom_parameters(
        self,
        f: TextIO,
        atom_types_used: Set[str],
        include_comments: bool,
    ):
        """
        Emit MM parameters for the HX link-atom type.

        Strategy: every cap H is type HX (a standard scaled hydrogen link
        atom). The FF has no entries for HX, so we copy them from the FF's
        natural H type for each pairing — `_proxy_h_type_for_link(T)` tells us
        which natural H type to use for bond/angle entries that pair HX with
        QM-parent type T. Where the FF has no such entry (the formamide cap's
        carbonyl C–H, absent from ff14SB), the values are supplied from GAFF
        via _SUPPLIED_CAP_BOND / _SUPPLIED_CAP_ANGLE (parm99 -> GAFF lookup).

        Emits, for each unique QM parent type T appearing in the link atoms:
          - HrmStr1  T  HX  k_b  r_eq       from T-<proxy> (or GAFF fallback)
          - HrmBnd1  X  T  HX  k_θ  θ_eq    for each (X, T) the FF defines
                                            as (X, T, proxy_H), or GAFF-supplied
        Plus a single VDW entry for HX (copied from the FF's plain "H" or
        "H1" — same Lennard-Jones values, just a relabel).

        Dihedrals through HX are not emitted explicitly: ff14SB's wildcard
        torsions (X-CX-N-X, X-CT-CT-X, ...) match any end atoms including HX
        once the central pair is in `atom_types_used`.
        """
        if not self.oniom_setup.link_atoms or self.ff_reader is None:
            return
        if LINK_ATOM_TYPE not in atom_types_used:
            return

        parm = self.oniom_setup.parm

        # Unique QM-parent atom types across all link atoms in this run.
        qm_parent_types: Set[str] = set()
        for link_atom in self.oniom_setup.link_atoms:
            qm_type = parm.atoms[link_atom.qm_parent_idx].type if parm else "CT"
            qm_parent_types.add(qm_type)

        hx_out = self._remap_type(LINK_ATOM_TYPE)

        # ---- VDW ----------------------------------------------------------
        # Borrow Lennard-Jones from a real H type (prefer the most generic).
        vdw_source = None
        for fallback in ("H", "H1", "HC", "HA"):
            vdw_source = self.ff_reader.get_nonbonded_parameter(fallback)
            if vdw_source is not None:
                break
        if vdw_source is not None:
            f.write(
                f"VDW  {hx_out:4s}  "
                f"{vdw_source.radius:8.4f}  {vdw_source.well_depth:10.6f}\n"
            )
        else:
            self.logger.warning(
                f"No FF VDW reference for link atom {LINK_ATOM_TYPE}; "
                "Gaussian may need a NonBon default"
            )

        # ---- Bonds: HX – T ------------------------------------------------
        bond_keys_written: Set[Tuple[str, str]] = set()
        for qm_type in sorted(qm_parent_types):
            proxy = _proxy_h_type_for_link(qm_type)
            bond = self.ff_reader.get_bond_parameter(qm_type, proxy)
            if bond is not None:
                kf, req = bond.force_constant, bond.eq_length
            elif qm_type in _SUPPLIED_CAP_BOND:
                # ff14SB has no qm_type–H bond (e.g. carbonyl C–H for the
                # formamide cap); supply it from GAFF.
                kf, req = _SUPPLIED_CAP_BOND[qm_type]
            else:
                self.logger.warning(
                    f"No FF bond reference for cap on QM parent {qm_type} "
                    f"(tried {qm_type}-{proxy}); cap will be unparameterized"
                )
                continue
            t_out = self._remap_type(qm_type)
            canon = tuple(sorted([t_out, hx_out]))
            if canon in bond_keys_written:
                continue
            f.write(
                f"HrmStr1  {t_out:4s}  {hx_out:4s}  "
                f"{kf:10.4f}  {req:8.4f}\n"
            )
            bond_keys_written.add(canon)

        # ---- Angles: X – T – HX ------------------------------------------
        angle_keys_written: Set[Tuple[str, str, str]] = set()
        for qm_type in sorted(qm_parent_types):
            proxy = _proxy_h_type_for_link(qm_type)
            t_out = self._remap_type(qm_type)
            for partner in sorted(atom_types_used):
                if partner == LINK_ATOM_TYPE:
                    continue
                angle = self.ff_reader.get_angle_parameter(partner, qm_type, proxy)
                if angle is not None:
                    kf, eq = angle.force_constant, angle.eq_angle
                elif (partner, qm_type) in _SUPPLIED_CAP_ANGLE:
                    # ff14SB has no partner–qm_type–H angle (the formamide
                    # cap's O–C–HX / N–C–HX); supply it from GAFF.
                    kf, eq = _SUPPLIED_CAP_ANGLE[(partner, qm_type)]
                else:
                    continue
                p_out = self._remap_type(partner)
                canon = (
                    (p_out, t_out, hx_out) if p_out <= hx_out
                    else (hx_out, t_out, p_out)
                )
                if canon in angle_keys_written:
                    continue
                f.write(
                    f"HrmBnd1  {p_out:4s}  {t_out:4s}  {hx_out:4s}  "
                    f"{kf:10.4f}  {eq:8.4f}\n"
                )
                angle_keys_written.add(canon)

        self.logger.debug(
            f"Link-atom params: VDW={int(vdw_source is not None)}, "
            f"bonds={len(bond_keys_written)}, angles={len(angle_keys_written)} "
            f"(for {len(qm_parent_types)} QM-parent type(s))"
        )

    def _get_atom_types_used(self) -> Set[str]:
        """
        Get set of all atom types used in the structure.

        Returns:
            Set of atom type strings
        """
        atom_types = set()

        # Real atoms
        for coords, atom_type in self.oniom_setup.atom_types.items():
            atom_types.add(atom_type)

        # All link caps are type HX. Standard FF tables don't know HX, so the
        # bond/angle iterators below won't emit anything for it — those entries
        # come from _write_link_atom_parameters() instead.
        if self.oniom_setup.link_atoms:
            atom_types.add(LINK_ATOM_TYPE)

        return atom_types

    def _write_vdw_parameters(
        self,
        f: TextIO,
        atom_types: Set[str],
        include_comments: bool
    ):
        """
        Write VDW (nonbonded) parameters.

        Format: VDW atom_type radius well_depth
        Example: VDW  Zn  1.10  0.0125
        """
        # VDW parameters (comment removed for Gaussian compatibility)
        # Deduplicate by remapped output name (colliding types have been renamed)

        n_written = 0
        seen_output = set()
        for atom_type in sorted(atom_types):
            out_type = self._remap_type(atom_type)
            if out_type in seen_output:
                continue
            param = self.ff_reader.get_nonbonded_parameter(atom_type)
            if param:
                f.write(f"VDW  {out_type:4s}  {param.radius:8.4f}  {param.well_depth:10.6f}\n")
                seen_output.add(out_type)
                n_written += 1

        if n_written == 0:
            self.logger.warning("No VDW parameters found for any atom types")
        else:
            self.logger.debug(f"Wrote {n_written} VDW parameters")

    def _write_bond_parameters(
        self,
        f: TextIO,
        atom_types: Set[str],
        include_comments: bool
    ):
        """
        Write harmonic bond stretch parameters.

        Format: HrmStr1 type1 type2 force_constant eq_length
        Example: HrmStr1  S  O  194.8000  1.8020
        """
        # Bond stretch parameters (comment removed for Gaussian compatibility)

        n_written = 0
        written_pairs = set()  # Deduplicate by remapped output names

        # Generate all pairs of atom types used in structure
        atom_type_list = sorted(atom_types)
        for i, type1 in enumerate(atom_type_list):
            for type2 in atom_type_list[i:]:  # Include self-pairs and avoid duplicates
                out1 = self._remap_type(type1)
                out2 = self._remap_type(type2)
                canon = tuple(sorted([out1, out2]))
                if canon in written_pairs:
                    continue

                param = self.ff_reader.get_bond_parameter(type1, type2)
                if param:
                    f.write(
                        f"HrmStr1  {out1:4s}  {out2:4s}  "
                        f"{param.force_constant:10.4f}  {param.eq_length:8.4f}\n"
                    )
                    written_pairs.add(canon)
                    n_written += 1

        if n_written == 0:
            self.logger.warning("No bond parameters found")
        else:
            self.logger.debug(f"Wrote {n_written} bond parameters")

    def _write_angle_parameters(
        self,
        f: TextIO,
        atom_types: Set[str],
        include_comments: bool
    ):
        """
        Write harmonic angle bend parameters.

        Format: HrmBnd1 type1 type2 type3 force_constant eq_angle
        Example: HrmBnd1  N  CT  HC  0.0000  0.0000
        """
        # Angle bend parameters (comment removed for Gaussian compatibility)

        n_written = 0
        written_triples = set()  # Deduplicate by remapped output names

        # Generate all triples of atom types
        atom_type_list = sorted(atom_types)
        for type1 in atom_type_list:
            for type2 in atom_type_list:
                for type3 in atom_type_list:
                    out1 = self._remap_type(type1)
                    out2 = self._remap_type(type2)
                    out3 = self._remap_type(type3)
                    canon = (out1, out2, out3) if out1 <= out3 else (out3, out2, out1)
                    if canon in written_triples:
                        continue

                    param = self.ff_reader.get_angle_parameter(type1, type2, type3)
                    if param:
                        f.write(
                            f"HrmBnd1  {out1:4s}  {out2:4s}  {out3:4s}  "
                            f"{param.force_constant:10.4f}  {param.eq_angle:8.4f}\n"
                        )
                        written_triples.add(canon)
                        n_written += 1

        # TIP3P water: AMBER uses SHAKE constraints instead of an explicit
        # angle term, so the prmtop has no HW-OW-HW angle. Gaussian needs it.
        if 'HW' in atom_types and 'OW' in atom_types:
            hw_out = self._remap_type('HW')
            ow_out = self._remap_type('OW')
            canon = (hw_out, ow_out, hw_out)
            if canon not in written_triples:
                f.write(
                    f"HrmBnd1  {hw_out:4s}  {ow_out:4s}  {hw_out:4s}  "
                    f"  100.0000  104.5200\n"
                )
                written_triples.add(canon)
                n_written += 1

        if n_written == 0:
            self.logger.warning("No angle parameters found")
        else:
            self.logger.debug(f"Wrote {n_written} angle parameters")

    def _write_dihedral_parameters(
        self,
        f: TextIO,
        atom_types: Set[str],
        include_comments: bool
    ):
        """
        Write dihedral torsion parameters in Gaussian AmbTrs format.

        Gaussian format (from amber.prm):
        AmbTrs t1 t2 t3 t4 phase1 phase2 phase3 phase4 pk1 pk2 pk3 pk4 idivf

        Where phases/pks are indexed by periodicity (1-4), and AMBER wildcard
        'X' is written as '*' for Gaussian.
        """
        n_written = 0

        # Collect unique dihedral keys (avoid writing both A-B-C-D and D-C-B-A)
        seen_keys = set()

        for key, terms in self.ff_reader.dihedral_params.items():
            # Parse the key into 4 atom types
            parts = key.split('-')
            if len(parts) != 4:
                continue

            # Map types to output names (wildcard X → *, others through remap)
            t1 = '*' if parts[0] == 'X' else self._remap_type(parts[0])
            t2 = '*' if parts[1] == 'X' else self._remap_type(parts[1])
            t3 = '*' if parts[2] == 'X' else self._remap_type(parts[2])
            t4 = '*' if parts[3] == 'X' else self._remap_type(parts[3])

            # Canonicalize by remapped output names to avoid duplicates
            canon = (t1, t2, t3, t4) if (t1, t2) <= (t4, t3) else (t4, t3, t2, t1)
            if canon in seen_keys:
                continue

            # Check if relevant: at least one central atom (positions 1,2) must
            # be in our atom type set, or be a wildcard (X)
            central1, central2 = parts[1], parts[2]
            c1_match = central1 == 'X' or central1 in atom_types
            c2_match = central2 == 'X' or central2 in atom_types
            if not (c1_match and c2_match):
                continue

            seen_keys.add(canon)

            # Consolidate multiple terms into periodicity slots 1-4
            phases = [0, 0, 0, 0]      # phases for periodicities 1, 2, 3, 4
            pks = [0.0, 0.0, 0.0, 0.0]  # barrier heights for periodicities 1, 2, 3, 4
            idivf = 1  # use first term's IDIVF (typically all the same)

            for term in terms:
                pn = int(abs(term.pn))
                if 1 <= pn <= 4:
                    phases[pn - 1] = int(term.phase)
                    pks[pn - 1] = term.pk
                    idivf = term.idivf if term.idivf > 0 else 1

            f.write(
                f"AmbTrs {t1:2s} {t2:2s} {t3:2s} {t4:2s}"
                f"  {phases[0]:3d} {phases[1]:3d} {phases[2]:3d} {phases[3]:3d}"
                f"  {pks[0]:6.3f} {pks[1]:6.3f} {pks[2]:6.3f} {pks[3]:6.3f}"
                f" {float(idivf):.1f}\n"
            )
            n_written += 1

        if n_written == 0:
            self.logger.warning("No dihedral parameters found")
        else:
            self.logger.debug(f"Wrote {n_written} dihedral (AmbTrs) parameters")

    def generate_summary_report(self) -> str:
        """
        Generate human-readable summary of ONIOM setup.

        Returns:
            Formatted string report
        """
        stats = calculate_layer_statistics(self.oniom_setup)

        report_lines = [
            "=" * 70,
            "ONIOM Setup Summary",
            "=" * 70,
            "",
            f"Number of layers: {self.oniom_setup.n_layers}",
            f"Real atoms in system: {len(self.idx_to_gaussian)}",
            f"Link H atoms added for QM model: {len(self.oniom_setup.link_atoms)}",
            "",
            "Layer distribution:",
            f"  HIGH:   {stats.n_high:4d} atoms in {len(stats.residues_high):3d} residues (charge: {stats.total_charge_high:+.2f})",
        ]

        if self.oniom_setup.n_layers == 3:
            report_lines.append(
                f"  MEDIUM: {stats.n_medium:4d} atoms in {len(stats.residues_medium):3d} residues (charge: {stats.total_charge_medium:+.2f})"
            )

        report_lines.extend([
            f"  LOW:    {stats.n_low:4d} atoms in {len(stats.residues_low):3d} residues",
            f"  LINK:   {stats.n_link:4d} hydrogen link atoms",
            "",
            f"Frozen atoms: {stats.n_frozen} / {stats.n_active} active",
            "",
            "QM Method:",
            f"  Functional: {self.oniom_setup.qm_functional}",
            f"  Basis set:  {self.oniom_setup.qm_basis_set}",
            f"  Charge:     {self.oniom_setup.qm_charge}",
            f"  Multiplicity: {self.oniom_setup.qm_multiplicity}",
            "",
        ])

        if self.oniom_setup.n_layers == 3:
            report_lines.extend([
                "Medium Method:",
                f"  Method:     {self.oniom_setup.medium_method}",
                f"  Charge:     {self.oniom_setup.medium_charge}",
                f"  Multiplicity: {self.oniom_setup.medium_multiplicity}",
                "",
            ])

        report_lines.extend([
            "MM Method:",
            f"  Forcefield: {self.oniom_setup.mm_forcefield}",
        ])

        if self.ff_reader:
            ff_stats = self.ff_reader.get_statistics()
            n_b = ff_stats.get('n_bonds', ff_stats.get('n_bond', 0))
            n_a = ff_stats.get('n_angles', ff_stats.get('n_angle', 0))
            n_d = ff_stats.get('n_dihedrals', ff_stats.get('n_dihedral', 0))
            report_lines.extend([
                f"  Parameters: {n_b} bonds, {n_a} angles, {n_d} dihedrals",
            ])

        report_lines.extend([
            "",
            "Job Settings:",
            f"  Job type:   {self.oniom_setup.job_type}",
            f"  Processors: {self.oniom_setup.n_processors}",
            f"  Memory:     {self.oniom_setup.memory_gb} GB",
            "",
        ])

        if self.oniom_setup.additional_keywords:
            report_lines.extend([
                f"Additional keywords: {self.oniom_setup.additional_keywords}",
                "",
            ])

        report_lines.append("HIGH layer residues:")
        parm = self.oniom_setup.parm
        if parm:
            for res_idx in sorted(stats.residues_high):
                res = parm.residues[res_idx]
                report_lines.append(f"  {res.name}{res.idx + 1}")
        else:
            report_lines.append(f"  {len(stats.residues_high)} residues")

        if self.oniom_setup.n_layers == 3 and stats.residues_medium:
            report_lines.extend(["", "MEDIUM layer residues:"])
            if parm:
                for res_idx in sorted(stats.residues_medium):
                    res = parm.residues[res_idx]
                    report_lines.append(f"  {res.name}{res.idx + 1}")
            else:
                report_lines.append(f"  {len(stats.residues_medium)} residues")

        report_lines.extend([
            "",
            "Link atoms:",
            f"  Total: {len(self.oniom_setup.link_atoms)}",
        ])

        # Classify link atoms by boundary type
        boundary_counts = {}
        for link_atom in self.oniom_setup.link_atoms:
            btype = link_atom.boundary_type
            boundary_counts[btype] = boundary_counts.get(btype, 0) + 1

        for btype, count in sorted(boundary_counts.items()):
            report_lines.append(f"  {btype}: {count}")

        report_lines.extend([
            "",
            "Validation:",
            f"  Status: {'PASSED' if self.oniom_setup.validation_passed else 'NOT RUN'}",
        ])

        if self.oniom_setup.validation_warnings:
            report_lines.append(f"  Warnings: {len(self.oniom_setup.validation_warnings)}")
            for warning in self.oniom_setup.validation_warnings[:5]:  # Show first 5
                report_lines.append(f"    - {warning}")
            if len(self.oniom_setup.validation_warnings) > 5:
                report_lines.append(f"    ... and {len(self.oniom_setup.validation_warnings) - 5} more")

        report_lines.extend([
            "",
            "=" * 70,
        ])

        return "\n".join(report_lines)
