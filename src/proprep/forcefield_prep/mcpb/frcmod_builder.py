"""
frcmod File Builder
Generates AMBER force field modification files with bonded parameters.

© 2024 ProPrep Developer. All rights reserved.

This module creates AMBER frcmod files containing force field parameters
for metal sites. The frcmod format is used to supplement or override
parameters in standard AMBER force field libraries.

frcmod file format:
    MASS
    <atom_type>  <mass>  <comment>

    BOND
    <type1>-<type2>  <force_constant>  <eq_length>  <comment>

    ANGL
    <type1>-<type2>-<type3>  <force_constant>  <eq_angle>  <comment>

    DIHE
    <type1>-<type2>-<type3>-<type4>  <parameters>  <comment>

    IMPR
    <improper parameters>

    NONB
    <atom_type>  <radius>  <well_depth>  <comment>

PROPRIETARY SOFTWARE: This file contains proprietary code belonging to the ProPrep project.
Unauthorized copying, distribution, or modification is strictly prohibited.
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

from .seminario import BondParameter, AngleParameter


@dataclass
class MassEntry:
    """Mass parameter entry."""
    atom_type: str
    mass: float
    comment: str = ""


@dataclass
class NonbondedEntry:
    """Van der Waals parameter entry."""
    atom_type: str
    radius: float     # Rmin/2 in Angstrom
    well_depth: float # ε in kcal/mol
    comment: str = ""


class FrcmodBuilder:
    """
    Build AMBER frcmod files from computed parameters.

    Assembles force field modification files containing masses, bonded
    parameters (bonds, angles, dihedrals), and non-bonded parameters.
    """

    def __init__(self):
        """Initialize empty frcmod builder."""
        self.masses: List[MassEntry] = []
        self.bonds: List[BondParameter] = []
        self.angles: List[AngleParameter] = []
        self.dihedrals: List[Dict] = []
        self.impropers: List[Dict] = []
        self.nonbonded: List[NonbondedEntry] = []

    def add_mass(self, atom_type: str, mass: float, comment: str = ""):
        """
        Add mass parameter.

        Args:
            atom_type: AMBER atom type (e.g., "M1", "L1")
            mass: Atomic mass in amu
            comment: Description
        """
        self.masses.append(MassEntry(atom_type, mass, comment))

    def add_bond_parameter(self, bond_param: BondParameter):
        """
        Add bond parameter from Seminario calculation.

        Args:
            bond_param: BondParameter object
        """
        self.bonds.append(bond_param)

    def add_angle_parameter(self, angle_param: AngleParameter):
        """
        Add angle parameter from Seminario calculation.

        Args:
            angle_param: AngleParameter object
        """
        self.angles.append(angle_param)

    def add_zero_dihedral(self, atom_types: Tuple[str, str, str, str],
                         comment: str = ""):
        """
        Add zero dihedral parameter for metal-containing torsions.

        MCPB convention: Dihedrals involving metal atoms are set to zero.

        Args:
            atom_types: (type1, type2, type3, type4) tuple
            comment: Description
        """
        self.dihedrals.append({
            "types": atom_types,
            "idivf": 3,
            "pk": 0.00,
            "phase": 0.00,
            "pn": 3,
            "comment": comment if comment else "Metal-containing dihedral (zero)"
        })

    def add_dihedral(self, atom_types: Tuple[str, str, str, str],
                    idivf: int, pk: float, phase: float, pn: float,
                    comment: str = ""):
        """
        Add a dihedral parameter with explicit values.

        Args:
            atom_types: (type1, type2, type3, type4) tuple
            idivf: Divisor for the torsional barrier
            pk: Barrier height in kcal/mol
            phase: Phase offset in degrees
            pn: Periodicity (negative = more terms follow)
            comment: Description
        """
        self.dihedrals.append({
            "types": atom_types,
            "idivf": idivf,
            "pk": pk,
            "phase": phase,
            "pn": pn,
            "comment": comment
        })

    def add_improper(self, atom_types: Tuple[str, str, str, str],
                    pk: float, phase: float, pn: float,
                    comment: str = ""):
        """
        Add an improper dihedral parameter.

        AMBER improper format: 3rd atom is the central atom.
        No IDIVF field in the IMPR section.

        Args:
            atom_types: (type1, type2, type3, type4) — 3rd is central
            pk: Barrier height in kcal/mol
            phase: Phase offset in degrees
            pn: Periodicity
            comment: Description
        """
        self.impropers.append({
            "types": atom_types,
            "pk": pk,
            "phase": phase,
            "pn": pn,
            "comment": comment
        })

    def add_nonbonded_parameter(self, atom_type: str, radius: float,
                               well_depth: float, comment: str = ""):
        """
        Add van der Waals parameter.

        Args:
            atom_type: AMBER atom type
            radius: Rmin/2 in Angstrom
            well_depth: ε in kcal/mol
            comment: Description
        """
        self.nonbonded.append(NonbondedEntry(atom_type, radius, well_depth, comment))

    def write_frcmod(self, output_file: Path, method: str = "Seminario",
                    title: Optional[str] = None):
        """
        Write frcmod file in AMBER format.

        Args:
            output_file: Output frcmod file path
            method: Method name for header comment
            title: Optional custom title line
        """
        with open(output_file, 'w') as f:
            # Title line
            if title:
                f.write(f"{title}\n")
            else:
                f.write(f"Metal center parameters generated by {method} method\n")

            # MASS section
            f.write("MASS\n")
            for entry in self.masses:
                comment_str = f"  {entry.comment}" if entry.comment else ""
                f.write(f"{entry.atom_type:2s}  {entry.mass:8.2f}{comment_str}\n")

            # BOND section
            f.write("\nBOND\n")
            for bond in self.bonds:
                type_str = f"{bond.atom1_type:2s}-{bond.atom2_type:2s}"
                fc = bond.force_constant
                eq = bond.eq_length

                # Generate comment with statistics
                if bond.std_dev > 0.01:
                    comment = f"  {method} (StdDev: {bond.std_dev:.1f})"
                else:
                    comment = f"  {method}"

                f.write(f"{type_str:5s}  {fc:6.1f}   {eq:7.4f}{comment}\n")

            # ANGL section
            f.write("\nANGL\n")
            for angle in self.angles:
                type_str = f"{angle.atom1_type:2s}-{angle.atom2_type:2s}-{angle.atom3_type:2s}"
                fc = angle.force_constant
                eq = angle.eq_angle

                # Generate comment with statistics
                if angle.std_dev > 0.01:
                    comment = f"  {method} (StdDev: {angle.std_dev:.2f})"
                else:
                    comment = f"  {method}"

                f.write(f"{type_str:8s}  {fc:7.2f}  {eq:7.2f}{comment}\n")

            # DIHE section
            f.write("\nDIHE\n")
            for dih in self.dihedrals:
                types = dih['types']
                type_str = f"{types[0]:2s}-{types[1]:2s}-{types[2]:2s}-{types[3]:2s}"
                idivf = dih['idivf']
                pk = dih['pk']
                phase = dih['phase']
                pn = dih['pn']
                comment = dih.get('comment', '')

                # AMBER dihedral format
                f.write(f"{type_str:11s}  {idivf:1d}  {pk:7.2f}  {phase:8.2f}  "
                       f"{pn:4.1f}  {comment}\n")

            # IMPR section
            f.write("\nIMPR\n")
            for imp in self.impropers:
                types = imp['types']
                type_str = f"{types[0]:2s}-{types[1]:2s}-{types[2]:2s}-{types[3]:2s}"
                pk = imp['pk']
                phase = imp['phase']
                pn = imp['pn']
                comment = imp.get('comment', '')
                f.write(f"{type_str:11s}       {pk:7.2f}  {phase:8.2f}  "
                       f"{pn:4.1f}  {comment}\n")

            # NONB section
            f.write("\nNONB\n")
            for nb in self.nonbonded:
                comment_str = f"  {nb.comment}" if nb.comment else ""
                f.write(f"  {nb.atom_type:2s}  {nb.radius:8.4f}  {nb.well_depth:10.6f}{comment_str}\n")

            # End with blank lines
            f.write("\n\n")

    def merge_with_pre_frcmod(self, pre_frcmod: Path, output_file: Path,
                             method: str = "Seminario"):
        """
        Merge bonded parameters with pre-generated frcmod file.

        The pre-frcmod file contains non-metal parameters (masses, VDW, etc.)
        with placeholders marked "NON" for metal-related bonds and angles.
        This method replaces those placeholders with computed values.

        Args:
            pre_frcmod: Path to pre-generated frcmod file
            output_file: Output path for merged file
            method: Method name for comments
        """
        # Build lookup dictionaries for computed parameters
        bond_lookup = {}
        for bond in self.bonds:
            key1 = f"{bond.atom1_type}-{bond.atom2_type}"
            key2 = f"{bond.atom2_type}-{bond.atom1_type}"
            bond_lookup[key1] = bond
            bond_lookup[key2] = bond

        angle_lookup = {}
        for angle in self.angles:
            key1 = f"{angle.atom1_type}-{angle.atom2_type}-{angle.atom3_type}"
            key2 = f"{angle.atom3_type}-{angle.atom2_type}-{angle.atom1_type}"
            angle_lookup[key1] = angle
            angle_lookup[key2] = angle

        # Read and modify pre-frcmod
        with open(pre_frcmod, 'r') as f_in, open(output_file, 'w') as f_out:
            # Update title
            title_line = f_in.readline()
            if "pre.frcmod" in title_line.lower() or "preliminary" in title_line.lower():
                f_out.write(f"Metal center parameters generated by {method} method\n")
            else:
                f_out.write(title_line)

            # Process rest of file
            for line in f_in:
                stripped = line.strip()

                # Drop any stray REMARK/comment line carried over from the
                # pre-frcmod. A frcmod permits exactly ONE free-text title line
                # (already written above); tleap reads any later REMARK line as
                # a section keyword and rejects it ("Unknown keyword"), which
                # desyncs its section parser.
                if stripped.startswith('REMARK'):
                    continue

                # Check for a parameter line still marked "NON" (a placeholder
                # whose QM-computed value must be substituted in).
                if stripped.startswith('NON') and '-' in stripped:
                    # The type string is everything after the "NON" marker and
                    # MAY contain a column-padding space when an atom type is a
                    # single character (e.g. "C -Y6-M1"). Splitting on
                    # whitespace would truncate that to "C" and silently fail
                    # the lookup, leaving the raw NON marker in the output and
                    # crashing tleap. So take the remainder verbatim (keep it
                    # for correct frcmod column alignment on output) and derive
                    # a padding-free key for the lookup.
                    raw_type = stripped[3:].strip()
                    lookup_key = '-'.join(t.strip() for t in raw_type.split('-'))
                    dash_count = lookup_key.count('-')

                    if dash_count == 2:
                        # Angle parameter
                        if lookup_key in angle_lookup:
                            angle = angle_lookup[lookup_key]
                            fc = angle.force_constant
                            eq = angle.eq_angle
                            if angle.std_dev > 0.01:
                                comment = f"  {method} (StdDev: {angle.std_dev:.2f})"
                            else:
                                comment = f"  {method}"
                            f_out.write(f"{raw_type:8s}  {fc:7.2f}  {eq:7.2f}{comment}\n")
                        else:
                            # Genuinely uncomputed — keep marker (surfaces loudly)
                            f_out.write(line)

                    elif dash_count == 1:
                        # Bond parameter
                        if lookup_key in bond_lookup:
                            bond = bond_lookup[lookup_key]
                            fc = bond.force_constant
                            eq = bond.eq_length
                            if bond.std_dev > 0.01:
                                comment = f"  {method} (StdDev: {bond.std_dev:.1f})"
                            else:
                                comment = f"  {method}"
                            f_out.write(f"{raw_type:5s}  {fc:6.1f}   {eq:7.4f}{comment}\n")
                        else:
                            # Genuinely uncomputed — keep marker (surfaces loudly)
                            f_out.write(line)
                    else:
                        # Unrecognized format, keep as is
                        f_out.write(line)

                # Check for lines marked "YES" - strip the marker and keep the content
                elif stripped.startswith('YES '):
                    # Remove "YES " prefix and write the rest
                    # The line format is: "YES <content>"
                    content = stripped[4:]  # Skip "YES "
                    f_out.write(f"{content}\n")

                else:
                    # Regular line (section headers, blank lines, etc.), copy as is
                    f_out.write(line)

                    # After writing a section header, append builder entries for that section
                    if stripped == "DIHE":
                        for dih in self.dihedrals:
                            # Skip zero-metal dihedrals (already in pre-frcmod)
                            if dih.get('comment', '').startswith('Treat as zero'):
                                continue
                            types = dih['types']
                            type_str = f"{types[0]:2s}-{types[1]:2s}-{types[2]:2s}-{types[3]:2s}"
                            f_out.write(
                                f"{type_str:11s}  {dih['idivf']:1d}  {dih['pk']:7.2f}  "
                                f"{dih['phase']:8.2f}  {dih['pn']:4.1f}    {dih.get('comment', '')}\n"
                            )
                    elif stripped == "IMPR":
                        for imp in self.impropers:
                            types = imp['types']
                            type_str = f"{types[0]:2s}-{types[1]:2s}-{types[2]:2s}-{types[3]:2s}"
                            f_out.write(
                                f"{type_str:11s}       {imp['pk']:7.2f}  {imp['phase']:8.2f}  "
                                f"{imp['pn']:4.1f}    {imp.get('comment', '')}\n"
                            )

    def generate_summary(self) -> Dict[str, int]:
        """
        Generate summary statistics of frcmod contents.

        Returns:
            Dict with counts of each parameter type
        """
        return {
            "n_masses": len(self.masses),
            "n_bonds": len(self.bonds),
            "n_angles": len(self.angles),
            "n_dihedrals": len(self.dihedrals),
            "n_impropers": len(self.impropers),
            "n_nonbonded": len(self.nonbonded)
        }

    def validate_parameters(self) -> Tuple[bool, List[str]]:
        """
        Validate parameters for common issues.

        Checks:
        - Non-negative force constants
        - Reasonable bond lengths (0.5-3.0 Å)
        - Reasonable angles (10-170 degrees)

        Returns:
            (is_valid, warnings) tuple
        """
        warnings = []

        # Check bonds
        for bond in self.bonds:
            if bond.force_constant < 0:
                warnings.append(
                    f"Negative bond force constant: {bond.atom1_type}-{bond.atom2_type} "
                    f"= {bond.force_constant:.1f}"
                )
            if bond.eq_length < 0.5 or bond.eq_length > 3.0:
                warnings.append(
                    f"Unusual bond length: {bond.atom1_type}-{bond.atom2_type} "
                    f"= {bond.eq_length:.3f} Å"
                )

        # Check angles
        for angle in self.angles:
            if angle.force_constant < 0:
                warnings.append(
                    f"Negative angle force constant: "
                    f"{angle.atom1_type}-{angle.atom2_type}-{angle.atom3_type} "
                    f"= {angle.force_constant:.2f}"
                )
            atrip = f"{angle.atom1_type}-{angle.atom2_type}-{angle.atom3_type}"
            if angle.eq_angle < 10.0:
                # A very acute equilibrium angle almost always means bad
                # connectivity (e.g. a spurious bond) or a distorted input
                # geometry, not a real chemical angle. Worth investigating.
                warnings.append(
                    f"Very acute equilibrium angle ({atrip} = {angle.eq_angle:.1f}°): "
                    f"likely a connectivity or geometry error, please verify."
                )
            elif angle.eq_angle > 175.0:
                # Near-linear angles are EXPECTED, not errors, for trans
                # ligand-metal-ligand arrangements in octahedral / square-planar /
                # linear coordination (the two donors sit ~180° apart). We flag
                # them only because a harmonic term k*(theta-theta0)^2 is
                # non-differentiable at exactly 180°: as the angle passes through
                # linearity the restoring force can invert and destabilize MD.
                # If equilibration blows up at this angle, that is the cause --
                # options are a larger-model reparameterization, or a restrained /
                # softened treatment of the near-linear angle. The parameter value
                # itself (straight from the QM geometry) is correct.
                warnings.append(
                    f"Near-linear equilibrium angle ({atrip} = {angle.eq_angle:.1f}°): "
                    f"expected for trans metal coordination; the value is fine, but "
                    f"harmonic angle terms can be numerically unstable near 180°."
                )

        is_valid = len(warnings) == 0
        return is_valid, warnings
