"""
Disulfide Bond Detector Worker

Core functionality for detecting and processing disulfide bonds in protein structures.
"""

import logging
import math
import os
from typing import Any, Dict, List, Optional, Set, Tuple

from Bio.PDB import Chain, Model, Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser


logger = logging.getLogger(__name__)


class DisulfideBondDetector:
    """Detects and processes disulfide bonds in protein structures."""

    def __init__(self):
        """Initialize the disulfide bond detector."""
        self.structure = None
        self.input_file = None
        self.selected_chains = []
        self.disulfide_bonds = []
        self.parser = PDBParser(QUIET=True)
        self.distance_threshold = 2.5

    def setup(self, structure=None, input_file=None, selected_chains=None):
        """
        Set up the detector with structure data.

        Args:
            structure: BioPython Structure object
            input_file: Path to PDB file to process
            selected_chains: List of selected chain IDs
        """
        self.structure = structure
        self.input_file = input_file
        self.selected_chains = selected_chains or []
        self.disulfide_bonds = []

        if self.structure is None and self.input_file is None:
            raise ValueError("Either structure or input_file must be provided")

    def detect_disulfide_bonds(self, interactive=True, use_ssbond_records=True):
        """
        Detect disulfide bonds in the structure.

        Args:
            interactive: Whether to prompt for user input
            use_ssbond_records: Whether to try parsing SSBOND records first

        Returns:
            List of disulfide bond tuples (chain_id1, res_id1, chain_id2, res_id2)
        """
        logger.info("Detecting disulfide bonds")

        potential_bonds = []

        if use_ssbond_records:
            bonds_from_records = self._parse_ssbond_records()
            if bonds_from_records:
                logger.info(
                    f"Found {len(bonds_from_records)} disulfide bond(s) from SSBOND records"
                )
                potential_bonds = bonds_from_records
            else:
                logger.info(
                    "No SSBOND records found, attempting distance-based detection"
                )
                potential_bonds = self._detect_bonds_by_distance()
                logger.info(
                    f"Found {len(potential_bonds)} disulfide bond(s) from distance analysis"
                )
        else:
            potential_bonds = self._detect_bonds_by_distance()
            logger.info(
                f"Found {len(potential_bonds)} disulfide bond(s) from distance analysis"
            )

        self.disulfide_bonds = potential_bonds if not interactive else potential_bonds

        return self.disulfide_bonds

    def _parse_ssbond_records(self):
        """
        Parse SSBOND records from the PDB header.

        Returns:
            List of disulfide bond tuples (chain_id1, res_id1, chain_id2, res_id2)
        """
        if not self.input_file:
            return []

        ssbonds = []

        try:
            with open(self.input_file, "r") as f:
                for line in f:
                    if line.startswith("SSBOND"):
                        chain1 = line[15:16].strip()
                        res_num1 = int(line[17:21].strip())
                        chain2 = line[29:30].strip()
                        res_num2 = int(line[31:35].strip())

                        if self.selected_chains and (
                            chain1 not in self.selected_chains
                            or chain2 not in self.selected_chains
                        ):
                            continue

                        ssbonds.append((chain1, res_num1, chain2, res_num2))
                        logger.debug(
                            f"Found SSBOND: Chain {chain1} Res {res_num1} - Chain {chain2} Res {res_num2}"
                        )

                    if line.startswith("ATOM") or line.startswith("MODEL"):
                        break
        except Exception as e:
            logger.error(f"Error parsing SSBOND records: {str(e)}")

        return ssbonds

    def _detect_bonds_by_distance(self):
        """
        Detect disulfide bonds by measuring distances between SG atoms of CYS residues.

        Returns:
            List of disulfide bond tuples (chain_id1, res_id1, chain_id2, res_id2)
        """
        if not self.structure:
            if not self.input_file:
                return []
            self.structure = self.parser.get_structure("structure", self.input_file)

        bonds = []
        cys_residues = []

        for chain in self.structure.get_chains():
            if self.selected_chains and chain.id not in self.selected_chains:
                continue

            for residue in chain:
                if residue.get_resname() in ["CYS", "CYX"]:
                    if "SG" in residue:
                        cys_residues.append((chain.id, residue.id[1], residue))

        for i in range(len(cys_residues)):
            chain_id1, res_num1, res1 = cys_residues[i]

            for j in range(i + 1, len(cys_residues)):
                chain_id2, res_num2, res2 = cys_residues[j]

                sg1 = res1["SG"]
                sg2 = res2["SG"]
                distance = sg1 - sg2

                if distance <= self.distance_threshold:
                    bonds.append((chain_id1, res_num1, chain_id2, res_num2))
                    logger.debug(
                        f"Found disulfide by distance: Chain {chain_id1} Res {res_num1} - "
                        f"Chain {chain_id2} Res {res_num2} (Distance: {distance:.2f} Å)"
                    )

        return bonds

    def detect_bonds_from_file(self, file_path, selected_chains=None):
        """
        Detect disulfide bonds by measuring distances between SG atoms directly from PDB file.

        Args:
            file_path: Path to the PDB file
            selected_chains: List of selected chain IDs

        Returns:
            List of disulfide bond tuples (chain_id1, res_id1, chain_id2, res_id2)
        """
        bonds = []
        cys_residues = []
        sg_coords = {}

        try:
            with open(file_path, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]

                        if selected_chains and chain_id not in selected_chains:
                            continue

                        residue_name = line[17:20].strip()

                        if residue_name in ["CYS", "CYX"]:
                            atom_name = line[12:16].strip()

                            if atom_name == "SG":
                                try:
                                    residue_id = int(line[22:26])
                                    x = float(line[30:38])
                                    y = float(line[38:46])
                                    z = float(line[46:54])

                                    cys_residues.append((chain_id, residue_id))
                                    sg_coords[(chain_id, residue_id)] = (x, y, z)

                                except (ValueError, IndexError):
                                    continue

            for i, (chain1, resid1) in enumerate(cys_residues):
                for j in range(i + 1, len(cys_residues)):
                    chain2, resid2 = cys_residues[j]

                    sg1 = sg_coords[(chain1, resid1)]
                    sg2 = sg_coords[(chain2, resid2)]

                    distance = math.sqrt(
                        (sg1[0] - sg2[0]) ** 2
                        + (sg1[1] - sg2[1]) ** 2
                        + (sg1[2] - sg2[2]) ** 2
                    )

                    if distance <= self.distance_threshold:
                        bonds.append((chain1, resid1, chain2, resid2))

        except Exception as e:
            logger.error(f"Error detecting disulfide bonds from file: {str(e)}")

        return bonds

    def verify_ssbonds_in_file(self, ssbond_records, file_path):
        """
        Verify SSBOND records by checking if both CYS residues exist in the file.

        Args:
            ssbond_records: List of (chain1, res_num1, chain2, res_num2) tuples
            file_path: Path to the PDB file

        Returns:
            List of verified SSBOND records
        """
        verified_bonds = []

        try:
            cys_residues = set()

            with open(file_path, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]
                        residue_name = line[17:20].strip()
                        try:
                            residue_id = int(line[22:26])
                        except ValueError:
                            continue

                        if residue_name in ["CYS", "CYX"]:
                            cys_residues.add((chain_id, residue_id))

            for chain1, res_num1, chain2, res_num2 in ssbond_records:
                if (chain1, res_num1) in cys_residues and (
                    chain2,
                    res_num2,
                ) in cys_residues:
                    verified_bonds.append((chain1, res_num1, chain2, res_num2))

        except Exception as e:
            logger.error(f"Error verifying SSBOND records: {str(e)}")

        return verified_bonds

    def update_cys_residues(self, output_file=None):
        """
        Update CYS residues involved in disulfide bonds to CYX.

        Args:
            output_file: Path to save the updated structure

        Returns:
            Path to the updated structure file
        """
        if not self.disulfide_bonds:
            logger.info("No disulfide bonds to process")
            return None

        if not self.structure:
            if not self.input_file:
                logger.error("No structure to update")
                return None
            self.structure = self.parser.get_structure("structure", self.input_file)

        if not output_file:
            if self.input_file:
                base = os.path.splitext(os.path.basename(self.input_file))[0]
                output_file = f"{base}_disulfide.pdb"
            else:
                output_file = "structure_disulfide.pdb"

        residues_to_update = set()
        for chain1, res_num1, chain2, res_num2 in self.disulfide_bonds:
            residues_to_update.add((chain1, res_num1))
            residues_to_update.add((chain2, res_num2))

        modified = False

        for chain in self.structure.get_chains():
            for residue in chain:
                if (chain.id, residue.id[1]) in residues_to_update:
                    if residue.get_resname() == "CYS" and "SG" in residue:
                        residue.resname = "CYX"
                        modified = True
                        logger.debug(
                            f"Updated Chain {chain.id} Residue {residue.id[1]} from CYS to CYX"
                        )

        if modified:
            io = PDBIO()
            io.set_structure(self.structure)
            output_path_str = str(output_file)
            io.save(output_path_str)
            logger.info(
                f"Saved structure with updated disulfide residues to {output_path_str}"
            )
            return output_path_str
        else:
            logger.warning("No residues were updated")
            return None

    def update_cys_in_file(self, input_file, output_file):
        """
        Update CYS residues to CYX directly in the PDB file for disulfide bonds.

        Args:
            input_file: Path to input PDB file
            output_file: Path to save the updated structure

        Returns:
            Path to saved file or None if no changes made
        """
        residues_to_update = set()
        for chain1, res_num1, chain2, res_num2 in self.disulfide_bonds:
            residues_to_update.add((chain1, res_num1))
            residues_to_update.add((chain2, res_num2))

        if not residues_to_update:
            return None

        modified = False
        updated_lines = []

        try:
            with open(input_file, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]
                        try:
                            residue_id = int(line[22:26])
                            residue_name = line[17:20].strip()

                            if (
                                chain_id,
                                residue_id,
                            ) in residues_to_update and residue_name == "CYS":
                                new_line = line[:17] + "CYX" + line[20:]
                                updated_lines.append(new_line)
                                modified = True
                                continue
                        except ValueError:
                            pass

                    updated_lines.append(line)

            if modified:
                with open(output_file, "w") as f:
                    f.writelines(updated_lines)
                return output_file
            else:
                return None

        except Exception as e:
            logger.error(f"Error updating CYS residues: {str(e)}")
            return None

    def get_tleap_commands(self):
        """
        Generate TLEaP commands for specifying disulfide bonds.

        Returns:
            Dictionary of categorized tLEaP command strings
        """
        if not self.disulfide_bonds:
            return {"disulfide_bonds": []}

        commands = {"disulfide_bonds": []}

        for chain1, res_num1, chain2, res_num2 in self.disulfide_bonds:
            command = f"bond mol.{res_num1}.SG mol.{res_num2}.SG"
            comment = f"Disulfide bond between Chain {chain1} Res {res_num1} and Chain {chain2} Res {res_num2}"
            command += f" # {comment}"
            commands["disulfide_bonds"].append(command)

        return commands

    def add_disulfide_bond(self, chain1, res_num1, chain2, res_num2):
        """Add a disulfide bond to the list."""
        bond = (chain1, res_num1, chain2, res_num2)
        if bond not in self.disulfide_bonds:
            self.disulfide_bonds.append(bond)

    def remove_disulfide_bond(self, chain1, res_num1, chain2, res_num2):
        """Remove a disulfide bond from the list."""
        bond = (chain1, res_num1, chain2, res_num2)
        if bond in self.disulfide_bonds:
            self.disulfide_bonds.remove(bond)

    def clear_disulfide_bonds(self):
        """Clear all disulfide bonds."""
        self.disulfide_bonds = []

    def get_cys_residues_from_file(self, file_path, selected_chains=None):
        """
        Get all CYS/CYX residues from a PDB file.

        Args:
            file_path: Path to the PDB file
            selected_chains: List of selected chain IDs

        Returns:
            List of (chain_id, residue_id, residue_name) tuples
        """
        cys_residues = []

        try:
            with open(file_path, "r") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]

                        if selected_chains and chain_id not in selected_chains:
                            continue

                        residue_name = line[17:20].strip()

                        if residue_name in ["CYS", "CYX"]:
                            try:
                                residue_id = int(line[22:26])
                                residue_key = (chain_id, residue_id)
                                if residue_key not in [
                                    (r[0], r[1]) for r in cys_residues
                                ]:
                                    cys_residues.append(
                                        (chain_id, residue_id, residue_name)
                                    )
                            except ValueError:
                                continue
        except Exception as e:
            logger.error(f"Error reading PDB file: {str(e)}")

        return cys_residues

    def get_cys_residues_from_structure(self, structure, selected_chains=None):
        """
        Get all CYS/CYX residues from a BioPython structure.

        Args:
            structure: BioPython Structure object
            selected_chains: List of selected chain IDs

        Returns:
            List of (chain_id, residue_id, residue_name) tuples
        """
        cys_residues = []

        for chain in structure.get_chains():
            if selected_chains and chain.id not in selected_chains:
                continue

            for residue in chain:
                if residue.get_resname() in ["CYS", "CYX"]:
                    if "SG" in residue:
                        cys_residues.append(
                            (chain.id, residue.id[1], residue.get_resname())
                        )

        return cys_residues
