"""
Protonation State Worker

Core functionality for analyzing protonation states of titratable residues.
"""

import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
from Bio.PDB import PDBIO

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    float_prompt_with_context,
    int_prompt_with_context,
)

logger = logging.getLogger(__name__)


# Standard amino acids for proper terminal classification
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE"  # Selenomethionine
}

PKA_VALUES = {
    "ASP": 3.65,  # Lehninger Principles of Biochemistry
    "GLU": 4.25,  # Lehninger Principles of Biochemistry
    "HIS": 6.00,  # Lehninger Principles of Biochemistry
    "CYS": 8.18,  # Lehninger Principles of Biochemistry
    "TYR": 10.07, # Lehninger Principles of Biochemistry
    "LYS": 10.53, # Lehninger Principles of Biochemistry
    "ARG": 12.48, # Lehninger Principles of Biochemistry
    "N_TERM": 8.0,
    "C_TERM": 3.1,
}

CHARGE_STATES = {
    "ASP": {"protonated": 0, "deprotonated": -1},
    "GLU": {"protonated": 0, "deprotonated": -1},
    "HIS": {"protonated": 1, "deprotonated": 0},
    "CYS": {"protonated": 0, "deprotonated": -1},
    "TYR": {"protonated": 0, "deprotonated": -1},
    "LYS": {"protonated": 1, "deprotonated": 0},
    "ARG": {"protonated": 1, "deprotonated": 0},
    "N_TERM": {"protonated": 1, "deprotonated": 0},
    "C_TERM": {"protonated": 0, "deprotonated": -1},
}


class ProtonationStateAnalyzer:
    """Core analyzer for protonation states of titratable residues."""

    def __init__(self, processor=None):
        """Initialize the protonation state analyzer."""
        self.titratable_residues = {}
        self.excluded_residues = set()
        self.problematic_residues = set()
        self.custom_pkas = {}
        # Per-residue PROPKA determinant breakdown (what drives each predicted
        # pKa: desolvation + side-chain/backbone/Coulombic determinants, with the
        # partner group and — crucially — whether it is a ligand/cofactor).
        # Keyed identically to custom_pkas: f"{chain}_{resnum}_{type}".
        self.determinants = {}
        self.results = {}
        self.net_charge = 0.0
        self.pH = 7.0
        self.threshold = 0.5
        self.structure = None
        self.input_file = None
        self.md_residue_names = {}
        self.propka_run_count = 0
        self.propka_max_attempts = 5
        self.processor = processor

    def setup(
        self,
        structure=None,
        input_file=None,
        excluded_residues=None,
        selected_chains=None,
        selected_residues=None,
        pH=7.0,
        threshold=0.5,
    ):
        """Initialize the analyzer with structure and parameters."""
        self.structure = structure
        self.input_file = input_file
        self.selected_chains = selected_chains or []
        self.selected_residues = selected_residues or {}

        if excluded_residues is not None:
            self.excluded_residues = excluded_residues
            if excluded_residues:
                logger.debug(
                    f"Excluding {len(excluded_residues)} residues from protonation analysis"
                )

        self.pH = pH
        self.threshold = threshold

        if self.structure:
            self.titratable_residues = self._parse_structure()
            logger.debug(
                f"Identified {len(self.titratable_residues)} titratable residues from structure"
            )
        elif self.input_file and os.path.exists(self.input_file):
            self.titratable_residues = self._parse_pdb_file()
            logger.debug(
                f"Identified {len(self.titratable_residues)} titratable residues"
            )
        else:
            logger.warning(
                "No input file or structure provided for protonation analysis"
            )

    def _parse_pdb_file(self):
        """Parse a PDB file to extract titratable residues."""
        titratable_residues = {}
        terminals = defaultdict(lambda: {"N": None, "C": None})
        titratable_types = set(PKA_VALUES.keys()) - {"N_TERM", "C_TERM"}

        with open(self.input_file, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    chain_id = line[21]
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26])
                    icode = line[26].strip()
                    atom_name = line[12:16].strip()

                    # Check exclusion with both empty string and space for icode compatibility
                    residue_id = (chain_id, res_num, icode)
                    residue_id_with_space = (chain_id, res_num, ' ' if icode == '' else icode)
                    if residue_id in self.excluded_residues or residue_id_with_space in self.excluded_residues:
                        continue

                    is_in_selected_chain = chain_id in self.selected_chains
                    is_selected_residue = False

                    if not is_in_selected_chain and chain_id in self.selected_residues:
                        is_selected_residue = any(
                            res_num == r[0] for r in self.selected_residues[chain_id]
                        )

                    if not is_in_selected_chain and not is_selected_residue:
                        continue

                    if res_name in titratable_types:
                        key = f"{chain_id}_{res_num}_{res_name}"
                        if key not in titratable_residues:
                            titratable_residues[key] = {
                                "chain": chain_id,
                                "number": res_num,
                                "type": res_name,
                                "icode": icode,
                                "pKa": PKA_VALUES[res_name],
                            }

                    # Only consider standard amino acids for terminal detection
                    if is_in_selected_chain and res_name in STANDARD_AMINO_ACIDS:
                        if atom_name == "N":
                            if (
                                terminals[chain_id]["N"] is None
                                or res_num < terminals[chain_id]["N"][0]
                            ):
                                terminals[chain_id]["N"] = (res_num, res_name, icode)

                        if atom_name == "C":
                            if (
                                terminals[chain_id]["C"] is None
                                or res_num > terminals[chain_id]["C"][0]
                            ):
                                terminals[chain_id]["C"] = (res_num, res_name, icode)

        for chain in self.selected_chains:
            terms = terminals[chain]
            if terms["N"] is not None:
                n_res_id = (chain, terms["N"][0], terms["N"][2])
                if n_res_id not in self.excluded_residues:
                    n_term = f"{chain}_{terms['N'][0]}_N_TERM"
                    titratable_residues[n_term] = {
                        "chain": chain,
                        "number": terms["N"][0],
                        "type": "N_TERM",
                        "icode": terms["N"][2],
                        "pKa": PKA_VALUES["N_TERM"],
                    }

            if terms["C"] is not None:
                c_res_id = (chain, terms["C"][0], terms["C"][2])
                if c_res_id not in self.excluded_residues:
                    c_term = f"{chain}_{terms['C'][0]}_C_TERM"
                    titratable_residues[c_term] = {
                        "chain": chain,
                        "number": terms["C"][0],
                        "type": "C_TERM",
                        "icode": terms["C"][2],
                        "pKa": PKA_VALUES["C_TERM"],
                    }

        return titratable_residues

    def _parse_structure(self):
        """Parse a BioPython structure object to extract titratable residues."""
        titratable_residues = {}
        titratable_types = set(PKA_VALUES.keys()) - {"N_TERM", "C_TERM"}

        if len(self.structure) == 0:
            logger.error("Structure contains no models")
            return titratable_residues

        model = list(self.structure.get_models())[0]

        for chain in model:
            chain_id = chain.id
            is_selected_chain = chain_id in self.selected_chains
            has_selected_residues = chain_id in self.selected_residues

            if not is_selected_chain and not has_selected_residues:
                continue

            residues = list(chain.get_residues())
            if not residues:
                continue

            residues.sort(key=lambda r: r.id[1])

            for residue in residues:
                if residue.id[0] != " ":
                    continue

                res_name = residue.resname
                res_num = residue.id[1]
                icode = residue.id[2]

                excluded_residue_ids = {
                    (chain, resnum) for chain, resnum, _ in self.excluded_residues
                }

                if (chain_id, res_num) in excluded_residue_ids:
                    logger.debug(
                        f"Excluding metal-coordinating residue: Chain {chain_id}, Residue {res_num} {res_name}"
                    )
                    continue

                if not is_selected_chain:
                    is_selected_residue = False
                    if has_selected_residues:
                        is_selected_residue = any(
                            res_num == r[0] for r in self.selected_residues[chain_id]
                        )

                    if not is_selected_residue:
                        continue

                if res_name in titratable_types:
                    key = f"{chain_id}_{res_num}_{res_name}"
                    if key not in titratable_residues:
                        titratable_residues[key] = {
                            "chain": chain_id,
                            "number": res_num,
                            "type": res_name,
                            "icode": icode,
                            "pKa": PKA_VALUES[res_name],
                        }

            # Find protein termini (only standard amino acids, not metals/waters/ligands)
            if is_selected_chain and residues:
                protein_residues = [res for res in residues if res.resname in STANDARD_AMINO_ACIDS]
                
                if protein_residues:
                    # N-terminal: first protein residue in chain
                    n_term = protein_residues[0]
                    if (chain_id, n_term.id[1]) not in excluded_residue_ids:
                        n_term_key = f"{chain_id}_{n_term.id[1]}_N_TERM"
                        titratable_residues[n_term_key] = {
                            "chain": chain_id,
                            "number": n_term.id[1],
                            "type": "N_TERM",
                            "icode": n_term.id[2],
                            "pKa": PKA_VALUES["N_TERM"],
                        }

                    # C-terminal: last protein residue in chain
                    c_term = protein_residues[-1]
                    if (chain_id, c_term.id[1]) not in excluded_residue_ids:
                        c_term_key = f"{chain_id}_{c_term.id[1]}_C_TERM"
                        titratable_residues[c_term_key] = {
                            "chain": chain_id,
                            "number": c_term.id[1],
                            "type": "C_TERM",
                            "icode": c_term.id[2],
                            "pKa": PKA_VALUES["C_TERM"],
                        }

        return titratable_residues

    def run_propka(self, retry_without=None):
        """Run PROPKA to predict pKa values with error handling and retry logic."""
        if not self.input_file:
            logger.warning("No input file provided for PROPKA analysis")
            return False

        self.propka_run_count += 1

        if self.propka_run_count > self.propka_max_attempts:
            logger.warning(
                f"Maximum PROPKA retry attempts ({self.propka_max_attempts}) reached"
            )
            return False

        filtered_pdb = self._create_filtered_pdb_for_propka()
        if not filtered_pdb:
            logger.warning("Failed to create filtered PDB for PROPKA analysis")
            return False

        input_file = filtered_pdb
        temp_file = None

        if retry_without:
            temp_file = self._create_filtered_pdb(retry_without, base_file=filtered_pdb)
            if temp_file:
                input_file = temp_file
                logger.info(
                    f"Created filtered PDB excluding {len(retry_without)} problematic residues"
                )
            else:
                logger.warning(
                    "Failed to create filtered PDB for retry, using original filtered file"
                )

        try:
            orig_base_name = os.path.splitext(os.path.basename(self.input_file))[0]
            chain_str = (
                "_".join(sorted(self.selected_chains))
                if self.selected_chains
                else "all"
            )
            output_base = f"{orig_base_name}_chains_{chain_str}_pH{self.pH}"

            if retry_without:
                output_base += f"_retry{self.propka_run_count}"

            propka_dir = os.path.join(
                os.path.dirname(self.input_file), "propka_results"
            )
            os.makedirs(propka_dir, exist_ok=True)

            output_pdb = os.path.join(propka_dir, f"{output_base}.pdb")
            output_pka = os.path.join(propka_dir, f"{output_base}.pka")

            shutil.copy2(input_file, output_pdb)

            logger.debug("Trying PROPKA as command-line tool")

            with tempfile.TemporaryDirectory() as temp_dir:
                cmd = ["propka3", output_pdb, "--quiet"]
                logger.debug(f"Running PROPKA on filtered PDB with selected chains")
                result = subprocess.run(cmd, capture_output=True, text=True)

                if result.returncode != 0:
                    error_msg = result.stderr
                    logger.warning(f"PROPKA failed with error: {error_msg}")

                    # Check for PROPKA internal errors that can't be resolved by excluding residues
                    if "ZeroDivisionError" in error_msg or "float division by zero" in error_msg:
                        logger.warning(
                            "PROPKA failed due to internal geometry calculation error. "
                            "This typically occurs with overlapping atoms or degenerate molecular geometry. "
                            "Using default pKa values instead."
                        )
                        
                        # Clean up and return success (using defaults)
                        if temp_file and os.path.exists(temp_file):
                            os.remove(temp_file)
                        if (
                            filtered_pdb
                            and filtered_pdb != self.input_file
                            and os.path.exists(filtered_pdb)
                        ):
                            os.remove(filtered_pdb)
                        
                        return True

                    problematic = self._extract_problematic_residues(error_msg)

                    if problematic:
                        self.problematic_residues.update(problematic)
                        exclude_set = set(self.problematic_residues)
                        logger.info(
                            f"Retrying PROPKA without {len(exclude_set)} problematic residues"
                        )

                        if temp_file and os.path.exists(temp_file):
                            os.remove(temp_file)

                        return self.run_propka(retry_without=exclude_set)

                    logger.warning(
                        "PROPKA failed without identifiable problematic residues, using default pKa values"
                    )

                    if temp_file and os.path.exists(temp_file):
                        os.remove(temp_file)

                    if (
                        filtered_pdb
                        and filtered_pdb != self.input_file
                        and os.path.exists(filtered_pdb)
                    ):
                        os.remove(filtered_pdb)

                    return True

                if os.path.exists(output_pka):
                    propka_pkas = self._parse_propka_output(output_pka)
                    if propka_pkas:
                        self.custom_pkas.update(propka_pkas)
                        self.determinants.update(
                            self._parse_propka_determinants(output_pka)
                        )
                        logger.debug(
                            f"Updated pKa values for {len(propka_pkas)} residues from PROPKA"
                        )
                        logger.debug(f"PROPKA output saved to: {output_pka}")

                        if temp_file and os.path.exists(temp_file):
                            os.remove(temp_file)

                        if (
                            filtered_pdb
                            and filtered_pdb != self.input_file
                            and os.path.exists(filtered_pdb)
                        ):
                            os.remove(filtered_pdb)

                        return True
                    else:
                        logger.warning(
                            "Failed to parse PROPKA output file, will use default pKa values"
                        )

                        if temp_file and os.path.exists(temp_file):
                            os.remove(temp_file)

                        if (
                            filtered_pdb
                            and filtered_pdb != self.input_file
                            and os.path.exists(filtered_pdb)
                        ):
                            os.remove(filtered_pdb)

                        return True
                    
                else:
                    # PROPKA generates output in current working directory with pattern *_pH*.pka
                    current_dir_pka = f"{output_base}.pka"
                    
                    logger.debug(f"Looking for PROPKA output file: {current_dir_pka}")
                    
                    if os.path.exists(current_dir_pka):
                        logger.debug(f"Found PROPKA output file: {current_dir_pka}")
                        # Move the file to propka_results directory
                        shutil.move(current_dir_pka, output_pka)
                        propka_pkas = self._parse_propka_output(output_pka)
                        if propka_pkas:
                            self.custom_pkas.update(propka_pkas)
                            self.determinants.update(
                                self._parse_propka_determinants(output_pka)
                            )
                            logger.debug(
                                f"Updated pKa values for {len(propka_pkas)} residues from PROPKA"
                            )
                            logger.debug(f"PROPKA output saved to: {output_pka}")
                    else:
                        logger.warning(
                            f"Could not find PROPKA output file: {current_dir_pka}"
                        )
                        logger.warning("Will use default pKa values")

                    if temp_file and os.path.exists(temp_file):
                        os.remove(temp_file)

                    if (
                        filtered_pdb
                        and filtered_pdb != self.input_file
                        and os.path.exists(filtered_pdb)
                    ):
                        os.remove(filtered_pdb)

                    return True

        except Exception as e:
            logger.error(f"Error running PROPKA: {str(e)}")

            if temp_file and os.path.exists(temp_file):
                os.remove(temp_file)

            if (
                filtered_pdb
                and filtered_pdb != self.input_file
                and os.path.exists(filtered_pdb)
            ):
                os.remove(filtered_pdb)

            logger.warning("Will use default pKa values")
            return True

    def _extract_problematic_residues(self, error_message):
        """Extract problematic residues from PROPKA error messages."""
        problematic = set()

        patterns = [
            r"Missing atoms or failed protonation for (\w+)\s+(\d+)\s+([A-Za-z])",
            r"Group \((\w+)\) for\s+\d+-\s+\w+\s+(\d+)-(\w+) \(([A-Za-z])\)",
            r"Expected \d+ interaction atoms for ([A-Za-z0-9]+), found:",
            r"([A-Z]{3})\s*([0-9]+)\s*([A-Z])",
        ]

        for pattern in patterns:
            matches = re.findall(pattern, error_message)

            for match in matches:
                try:
                    if len(match) == 3:
                        restype, resnum, chain = match
                        problematic.add((chain, int(resnum), restype))
                    elif len(match) == 4:
                        group, resnum, restype, chain = match
                        problematic.add((chain, int(resnum), restype))
                except (ValueError, IndexError):
                    continue

        if problematic:
            logger.info(
                f"Identified {len(problematic)} problematic residues from PROPKA error message"
            )

        return problematic

    def _create_filtered_pdb_for_propka(self):
        """Create a filtered PDB file containing only selected chains and residues."""
        if self.structure:
            try:
                fd, temp_path = tempfile.mkstemp(suffix=".pdb")
                os.close(fd)

                class ChainSelect:
                    def __init__(self, selected_chains, selected_residues):
                        self.selected_chains = selected_chains
                        self.selected_residues = selected_residues

                    def accept_model(self, model):
                        return True

                    def accept_chain(self, chain):
                        return chain.id in self.selected_chains

                    def accept_residue(self, residue):
                        chain_id = residue.get_parent().id
                        res_num = residue.id[1]

                        if chain_id in self.selected_chains:
                            return True

                        if chain_id in self.selected_residues:
                            return any(
                                res_num == r[0]
                                for r in self.selected_residues[chain_id]
                            )

                        return False

                    def accept_atom(self, atom):
                        return True

                io = PDBIO()
                io.set_structure(self.structure)
                io.save(
                    temp_path, ChainSelect(self.selected_chains, self.selected_residues)
                )

                logger.info(
                    f"Created filtered PDB file from structure for PROPKA: {temp_path}"
                )
                return temp_path

            except Exception as e:
                logger.error(f"Error creating filtered PDB from structure: {str(e)}")
                logger.info("Falling back to creating filtered PDB from input file")

        if not self.input_file or not os.path.exists(self.input_file):
            logger.warning("Cannot create filtered PDB: input file not found")
            return None

        if not self.selected_chains and not self.selected_residues:
            return self.input_file

        try:
            fd, temp_path = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)

            with open(self.input_file, "r") as infile, open(temp_path, "w") as outfile:
                for line in infile:
                    if not line.startswith(("ATOM", "HETATM", "TER", "END")):
                        outfile.write(line)
                    else:
                        break

            with open(self.input_file, "r") as infile, open(temp_path, "a") as outfile:
                for line in infile:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]
                        res_num = int(line[22:26])

                        if chain_id in self.selected_chains:
                            outfile.write(line)
                            continue

                        if chain_id in self.selected_residues:
                            if any(
                                res_num == r[0]
                                for r in self.selected_residues[chain_id]
                            ):
                                outfile.write(line)
                                continue

                    elif line.startswith("TER"):
                        chain_id = line[21] if len(line) >= 22 else ""
                        if chain_id in self.selected_chains:
                            outfile.write(line)
                    elif line.startswith("END"):
                        outfile.write(line)

            logger.debug(
                f"Created filtered PDB file with selected chains for PROPKA: {temp_path}"
            )
            return temp_path

        except Exception as e:
            logger.error(f"Error creating filtered PDB for PROPKA: {str(e)}")
            if "temp_path" in locals() and os.path.exists(temp_path):
                os.remove(temp_path)
            return None

    def _create_filtered_pdb(self, exclude_residues, base_file=None):
        """Create a filtered PDB file excluding specific residues."""
        input_file = base_file if base_file else self.input_file

        if not input_file or not os.path.exists(input_file):
            logger.warning(
                f"Cannot create filtered PDB: base file {input_file} not found"
            )
            return None

        try:
            fd, temp_path = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)

            exclude_ids = {(chain, resnum) for chain, resnum, _ in exclude_residues}

            with open(input_file, "r") as infile, open(temp_path, "w") as outfile:
                for line in infile:
                    if line.startswith(("ATOM", "HETATM")):
                        chain_id = line[21]
                        try:
                            res_num = int(line[22:26])
                            if (chain_id, res_num) not in exclude_ids:
                                outfile.write(line)
                        except ValueError:
                            outfile.write(line)
                    else:
                        outfile.write(line)

            logger.info(
                f"Created filtered PDB file excluding problematic residues: {temp_path}"
            )
            return temp_path

        except Exception as e:
            logger.error(f"Error creating filtered PDB: {str(e)}")
            if "temp_path" in locals() and os.path.exists(temp_path):
                os.remove(temp_path)
            return None

    def _parse_propka_output(self, pka_file):
        """Parse PROPKA output file to extract pKa values.

        PROPKA uses fixed-width label formatting:
          ``{res_type:<3s}{res_num:>4d}{chain:>2s}``
        For residue numbers < 1000 the fields are space-separated, but at
        >= 1000 the residue type and number merge (e.g. ``ASP1012 F``).
        We use a regex to handle both cases reliably.
        """
        import re

        custom_pkas = {}
        parsing_pka_section = False

        res_type_map = {
            "ASP": "ASP",
            "GLU": "GLU",
            "HIS": "HIS",
            "CYS": "CYS",
            "TYR": "TYR",
            "LYS": "LYS",
            "ARG": "ARG",
            "N+": "N_TERM",
            "C-": "C_TERM",
        }

        # Matches: optional whitespace, residue type (2-3 alpha or N+/C-),
        # optional whitespace, residue number, whitespace, chain ID,
        # whitespace, pKa value
        pka_line_re = re.compile(
            r"^\s*([A-Z]{2,3}|[NC][+-])\s*(\d+)\s+([A-Z])\s+([\d.]+)"
        )

        try:
            with open(pka_file, "r") as f:
                for line in f:
                    stripped = line.strip()

                    if "propKa" in stripped or "SUMMARY OF THIS PREDICTION" in stripped:
                        parsing_pka_section = True
                        continue

                    if parsing_pka_section:
                        if stripped.startswith("-") or stripped == "":
                            continue

                        if "Group" in stripped and "pKa" in stripped and "model-pKa" in stripped:
                            continue

                        m = pka_line_re.match(stripped)
                        if m:
                            res_type = m.group(1)
                            res_num = m.group(2)
                            chain = m.group(3)
                            pka_value = float(m.group(4))

                            if res_type in res_type_map:
                                mapped_type = res_type_map[res_type]
                                key = f"{chain}_{res_num}_{mapped_type}"
                                custom_pkas[key] = pka_value

        except Exception as e:
            logger.error(f"Error parsing PROPKA output: {str(e)}")
            return {}

        return custom_pkas

    def _parse_propka_determinants(self, pka_file):
        """Parse the per-residue determinant breakdown from a PROPKA .pka file.

        PROPKA's detailed table explains every predicted pKa as a sum of a
        desolvation penalty plus named side-chain / backbone / Coulombic
        determinants. Each determinant names the partner group; for ligands
        (cofactors, caps, ions) the partner is shown as ``RESNAME ATOM CHAIN``
        rather than ``RESNAME RESNUM CHAIN``, which is how we tell a cofactor
        contribution apart from a protein one.

        Returns {f"{chain}_{resnum}_{type}": {pKa, buried, desolvation_regular,
        desolvation_re, determinants: [ {kind, value, partner_resname,
        partner_id, partner_chain, is_ligand} ... ]}}.
        """
        res_type_map = {
            "ASP": "ASP", "GLU": "GLU", "HIS": "HIS", "CYS": "CYS",
            "TYR": "TYR", "LYS": "LYS", "ARG": "ARG",
            "N+": "N_TERM", "C-": "C_TERM",
        }
        determinants = {}
        try:
            lines = open(pka_file, "r").read().splitlines()
        except Exception as e:
            logger.error(f"Error reading PROPKA file for determinants: {e}")
            return determinants

        # Find the detailed table header ("RESIDUE  pKa  BURIED ...") and use the
        # dashed separator line beneath it to locate column spans robustly.
        start = None
        for i, ln in enumerate(lines):
            if ln.strip().startswith("RESIDUE") and "pKa" in ln and "BURIED" in ln:
                start = i
                break
        if start is None or start + 1 >= len(lines):
            return determinants

        spans = [(m.start(), m.end()) for m in re.finditer(r"-+", lines[start + 1])]
        if len(spans) < 8:
            return determinants
        det_cols = spans[5:8]  # sidechain, backbone, Coulombic
        kinds = ("sidechain", "backbone", "coulombic")

        cur = None
        for ln in lines[start + 2:]:
            if ln.startswith("----") or "SUMMARY OF THIS PREDICTION" in ln:
                break
            if not ln.strip():
                continue
            label = ln[spans[0][0]:spans[0][1] + 1].strip()
            m = re.match(r"([A-Z]{2,3}|[NC][+-])\s+(\S+)\s+([A-Za-z0-9])$", label)
            if not m:
                continue
            resname, resid, chain = m.group(1), m.group(2), m.group(3)
            mapped = res_type_map.get(resname, resname)
            key = f"{chain}_{resid}_{mapped}"

            pka_field = ln[spans[1][0]:spans[1][1]].strip()
            if pka_field:  # header row → new residue block
                desolv = ln[spans[3][0]:spans[4][1]].split()
                cur = determinants.setdefault(key, {
                    "pKa": float(pka_field),
                    "buried": ln[spans[2][0]:spans[2][1]].strip(),
                    "desolvation_regular": float(desolv[0]) if desolv else 0.0,
                    "desolvation_re": float(desolv[2]) if len(desolv) > 2 else 0.0,
                    "determinants": [],
                })
            if cur is None:
                continue

            for ci, (cs, ce) in enumerate(det_cols):
                cell = ln[cs - 1:ce + 4].strip()
                cm = re.match(r"(-?\d+\.\d+)\s+(\S+)\s+(\S+)\s+(\S+)", cell)
                if not cm:
                    continue
                val, p_resname, p_id, p_chain = cm.groups()
                if p_resname == "XXX" or abs(float(val)) < 1e-9:
                    continue
                cur["determinants"].append({
                    "kind": kinds[ci],
                    "value": float(val),
                    "partner_resname": p_resname,
                    "partner_id": p_id,
                    "partner_chain": p_chain,
                    "is_ligand": not p_id.isdigit(),
                })
        return determinants

    def analyze_protonation_states(self):
        """Analyze protonation states of all titratable residues at the specified pH."""
        results = {}
        net_charge = 0.0

        for key, residue in self.titratable_residues.items():
            pka = residue["pKa"]
            if self.custom_pkas and key in self.custom_pkas:
                pka = self.custom_pkas[key]

            is_protonated, probability = self._determine_protonation_state(
                self.pH, pka, self.threshold
            )
            
            # Calculate expected fractional charge based on probability
            expected_charge = (
                probability * CHARGE_STATES[residue["type"]]["protonated"]
                + (1 - probability) * CHARGE_STATES[residue["type"]]["deprotonated"]
            )

            results[key] = {
                "chain": residue["chain"],
                "number": residue["number"],
                "type": residue["type"],
                "pKa": pka,
                "protonated": is_protonated,
                "probability": probability,
                "charge": expected_charge,  # Use fractional charge instead of discrete
                # PROPKA determinant breakdown ({} when PROPKA was not run /
                # fell back to defaults); consumed by the determinant-inspection
                # view in the analyzer.
                "determinants": self.determinants.get(key, {}),
            }

            net_charge += expected_charge

        self.results = results
        self.net_charge = net_charge

        # Calculate ideal charge capacitance

        logger.debug(
            f"Analyzed protonation states at pH {self.pH}, net charge: {net_charge:.2f}"
        )
        return results, net_charge

    def calculate_ideal_capacitance(self, pH=None):
        """
        Calculate the ideal charge capacitance of the system.
        
        Formula:
        C_ideal = Σ_i (10^(pH-pKi)) / (1 + 10^(pH-pKi))^2
        
        Where i runs over all individual titratable residues.
        Results are then grouped by residue type for analysis.
        """
        if pH is None:
            pH = self.pH
        
        total_capacitance = 0.0
        capacitance_per_group = {}
        individual_capacitances = {}  # Store individual residue capacitances
        
        # Calculate capacitance for each individual residue
        for key, residue in self.titratable_residues.items():
            res_type = residue["type"]
            
            # Get pKa for this specific residue
            if key in self.custom_pkas:
                pKa = self.custom_pkas[key]
            elif res_type in PKA_VALUES:
                pKa = PKA_VALUES[res_type]
            else:
                continue
                
            # Calculate the Henderson-Hasselbalch term
            hh_term = 10 ** (pH - pKa)
            
            # Calculate capacitance for this individual residue
            residue_capacitance = hh_term / ((1 + hh_term) ** 2)
            
            # Store individual capacitance
            individual_capacitances[key] = {
                "capacitance": residue_capacitance,
                "pKa": pKa,
                "type": res_type
            }
            
            # Add to total
            total_capacitance += residue_capacitance
            
            # Initialize group if not present
            if res_type not in capacitance_per_group:
                capacitance_per_group[res_type] = {
                    "count": 0,
                    "capacitance": 0.0,
                    "pKa_values": [],
                    "individual_capacitances": []
                }
            
            # Add to group statistics
            capacitance_per_group[res_type]["count"] += 1
            capacitance_per_group[res_type]["capacitance"] += residue_capacitance
            capacitance_per_group[res_type]["pKa_values"].append(pKa)
            capacitance_per_group[res_type]["individual_capacitances"].append(residue_capacitance)
        
        # Calculate additional statistics for each group
        for res_type, data in capacitance_per_group.items():
            count = data["count"]
            if count > 0:
                # Average capacitance per residue of this type
                data["capacitance_per_residue"] = data["capacitance"] / count
                
                # Average pKa for this residue type (useful when using PROPKA)
                data["avg_pKa"] = sum(data["pKa_values"]) / count
                
                # Min and max pKa values
                data["min_pKa"] = min(data["pKa_values"])
                data["max_pKa"] = max(data["pKa_values"])
                
                # Standard deviation of pKa values (if more than one residue)
                if count > 1:
                    avg_pKa = data["avg_pKa"]
                    variance = sum((pKa - avg_pKa) ** 2 for pKa in data["pKa_values"]) / count
                    data["pKa_std"] = variance ** 0.5
                else:
                    data["pKa_std"] = 0.0
        
        # Calculate total number of titratable residues
        total_titratable = len(self.titratable_residues)
        
        # Calculate specific capacitance (capacitance per titratable residue)
        specific_capacitance = total_capacitance / total_titratable if total_titratable > 0 else 0
        
        # Store individual capacitances for potential future use
        self.individual_capacitances = individual_capacitances
        
        return total_capacitance, capacitance_per_group, specific_capacitance, total_titratable

    def generate_capacitance_profile(self, pH_range=(0, 14), step=0.1):
        """
        Generate ideal charge capacitance profile over a pH range.
        
        This method now properly handles individual residue pKa values.
        """
        pH_values = []
        capacitances = []
        specific_capacitances = []
        detailed_profiles = {}
        
        # Initialize detailed profiles for each residue type found in the structure
        unique_res_types = set()
        for residue in self.titratable_residues.values():
            res_type = residue["type"]
            unique_res_types.add(res_type)
        
        for res_type in unique_res_types:
            detailed_profiles[res_type] = []
        
        # Calculate capacitance at each pH
        current_pH = pH_range[0]
        while current_pH <= pH_range[1]:
            pH_values.append(current_pH)
            
            # Calculate total capacitance using the revised method
            total_cap, cap_per_group, specific_cap, total_titratable = self.calculate_ideal_capacitance(current_pH)
            capacitances.append(total_cap)
            specific_capacitances.append(specific_cap)
            
            # Store detailed capacitance for each group
            for res_type in detailed_profiles:
                if res_type in cap_per_group:
                    detailed_profiles[res_type].append(cap_per_group[res_type]["capacitance"])
                else:
                    detailed_profiles[res_type].append(0.0)
            
            current_pH += step
        
        # Find pH of maximum capacitance
        max_cap_idx = capacitances.index(max(capacitances)) if capacitances else 0
        pH_max_cap = pH_values[max_cap_idx] if pH_values else 0
        
        # Also find max specific capacitance
        max_spec_cap_idx = specific_capacitances.index(max(specific_capacitances)) if specific_capacitances else 0
        pH_max_spec_cap = pH_values[max_spec_cap_idx] if pH_values else 0
        
        # Collect pKa statistics for the final report
        pka_statistics = {}
        if hasattr(self, 'custom_pkas') and self.custom_pkas:
            # We have PROPKA values - calculate statistics
            _, final_cap_per_group, _, _ = self.calculate_ideal_capacitance(self.pH)
            for res_type, data in final_cap_per_group.items():
                pka_statistics[res_type] = {
                    "avg_pKa": data.get("avg_pKa", PKA_VALUES.get(res_type, 0)),
                    "min_pKa": data.get("min_pKa", PKA_VALUES.get(res_type, 0)),
                    "max_pKa": data.get("max_pKa", PKA_VALUES.get(res_type, 0)),
                    "pKa_std": data.get("pKa_std", 0),
                    "count": data.get("count", 0)
                }
        else:
            # Using default pKa values
            for res_type in unique_res_types:
                count = sum(1 for r in self.titratable_residues.values() if r["type"] == res_type)
                pka_statistics[res_type] = {
                    "avg_pKa": PKA_VALUES.get(res_type, 0),
                    "min_pKa": PKA_VALUES.get(res_type, 0),
                    "max_pKa": PKA_VALUES.get(res_type, 0),
                    "pKa_std": 0.0,
                    "count": count
                }
        
        return {
            "pH_values": pH_values,
            "total_capacitance": capacitances,
            "specific_capacitance": specific_capacitances,
            "group_capacitances": detailed_profiles,
            "group_counts": {res_type: sum(1 for r in self.titratable_residues.values() if r["type"] == res_type) 
                            for res_type in unique_res_types},
            "pH_max_capacitance": pH_max_cap,
            "max_capacitance": capacitances[max_cap_idx] if capacitances else 0,
            "pH_max_specific_capacitance": pH_max_spec_cap,
            "max_specific_capacitance": specific_capacitances[max_spec_cap_idx] if specific_capacitances else 0,
            "pka_statistics": pka_statistics
        }
        
    def _calculate_protonation_probability(self, pH, pKa):
        """Calculate the probability of protonation using the Henderson-Hasselbalch equation."""
        return 1.0 / (1.0 + 10 ** (pH - pKa))

    def _determine_protonation_state(self, pH, pKa, threshold=0.5):
        """Determine if a residue is protonated based on pH, pKa, and threshold."""
        probability = self._calculate_protonation_probability(pH, pKa)
        is_protonated = probability >= threshold
        return is_protonated, probability

    def recommend_md_residue_names(self, constant_ph_simulation=False):
        """Recommend residue names for MD simulations based on protonation states."""
        recommendations = {}

        # Residues that can be titrated in AMBER constant pH MD
        # See: cpinutil.py --describe
        # AS4 (ASP), GL4 (GLU), CYS, TYR, HIP (HIS), LYS, PRN (heme propionate)
        # N_TERM and C_TERM are NOT supported
        amber_cphmd_titratable = {"ASP", "GLU", "HIS", "TYR", "LYS", "CYS", "PRN"}

        amber_residue_map = {
            "ASP": {
                "protonated": "ASH",
                "deprotonated": "ASP",
                "titratable": "AS4",
            },
            "GLU": {
                "protonated": "GLH",
                "deprotonated": "GLU",
                "titratable": "GL4",
            },
            "HIS": {
                "protonated": "HIP",
                "deprotonated_ND1": "HID",
                "deprotonated_NE2": "HIE",
            },
            "CYS": {
                "protonated": "CYS",
                "deprotonated": "CYM",
            },
            "LYS": {
                "protonated": "LYS",
                "deprotonated": "LYN",
            },
            "TYR": {
                "protonated": "TYR",
                "deprotonated": "TYM",
            },
        }

        for key, res in self.results.items():
            residue_type = res["type"]
            is_protonated = res["protonated"]

            # Skip non-titratable residues in constant pH simulations
            # (ARG, N_TERM, C_TERM, and other residues not supported by AMBER CpHMD)
            if constant_ph_simulation and residue_type not in amber_cphmd_titratable:
                continue

            if residue_type in amber_residue_map:
                if constant_ph_simulation and residue_type in ["ASP", "GLU"]:
                    recommendations[key] = amber_residue_map[residue_type]["titratable"]
                elif residue_type == "HIS" and not is_protonated:
                    recommendations[key] = amber_residue_map[residue_type][
                        "deprotonated_ND1"
                    ]
                else:
                    state = "protonated" if is_protonated else "deprotonated"
                    recommendations[key] = amber_residue_map[residue_type][state]
            else:
                recommendations[key] = residue_type

        return recommendations

    def get_default_recommendations(self, constant_ph_simulation=False,
                                    standard_states=False):
        """Compute default MD residue name recommendations without prompts.

        Uses the same logic as set_md_residue_names but returns the defaults
        without any interactive prompts or console output.

        Args:
            constant_ph_simulation: If True, use constant pH naming conventions
            standard_states: If True, ignore the PROPKA per-residue analysis and
                assign canonical physiological-pH states for every residue:
                ASP/GLU deprotonated, CYS/LYS/ARG/TYR protonated, HIS->HIE, and
                heme propionate PRN->PRD (deprotonated). This is mutually
                exclusive with constant_ph_simulation. PRN->PRD in particular is
                something PROPKA cannot recommend on its own because it does not
                model non-standard groups.

        Returns:
            Dict mapping residue keys (e.g. 'A_94_HIS') to AMBER names (e.g. 'HIE')
        """
        recommendations = {}

        if not self.results:
            return recommendations

        # AMBER residue naming conventions
        amber_residue_map = {
            "ASP": {"protonated": "ASH", "deprotonated": "ASP", "titratable": "AS4"},
            "GLU": {"protonated": "GLH", "deprotonated": "GLU", "titratable": "GL4"},
            "HIS": {"protonated": "HIP", "deprotonated_ND1": "HID",
                     "deprotonated_NE2": "HIE", "titratable": "HIP"},
            "CYS": {"protonated": "CYS", "deprotonated": "CYM", "titratable": "CYS"},
            "TYR": {"protonated": "TYR", "deprotonated": "TYM", "titratable": "TYR"},
            "LYS": {"protonated": "LYS", "deprotonated": "LYN", "titratable": "LYS"},
            "PRN": {"protonated": "PRN", "deprotonated": "PRN", "titratable": "PRN"},
        }

        # Canonical fixed protonation states used when standard_states=True.
        # Acids deprotonated, bases/thiol/phenol protonated, HIS as the more
        # common epsilon tautomer, heme propionate as the deprotonated acid.
        standard_states_map = {
            "ASP": "ASP", "GLU": "GLU", "HIS": "HIE",
            "CYS": "CYS", "TYR": "TYR", "LYS": "LYS",
            "PRN": "PRD",
        }

        amber_cphmd_titratable = {"ASP", "GLU", "HIS", "TYR", "LYS", "CYS", "PRN"}

        # Filter for CpHMD if needed
        working_results = self.results.copy()
        if constant_ph_simulation:
            working_results = {
                key: res for key, res in self.results.items()
                if res["type"] in amber_cphmd_titratable
            }

        # Filter out terminal residues for standard MD
        if not constant_ph_simulation:
            terminal_positions = set()
            for key, res in working_results.items():
                if res["type"] in ("N_TERM", "C_TERM"):
                    terminal_positions.add((res["chain"], res["number"]))

            working_results = {
                key: res for key, res in working_results.items()
                if res["type"] not in ("N_TERM", "C_TERM")
                and not (res["type"] in amber_residue_map
                         and (res["chain"], res["number"]) in terminal_positions)
            }

        # Compute default recommendation for each residue
        for key, res in working_results.items():
            res_type = res["type"]
            is_protonated = res["protonated"]
            probability = res.get("probability", 0.5)

            if res_type not in amber_residue_map:
                recommendations[key] = res_type
                continue

            if standard_states:
                # Force canonical fixed states regardless of the analysis.
                recommendations[key] = standard_states_map.get(res_type, res_type)
                continue

            if constant_ph_simulation and res_type in ["ASP", "GLU"]:
                # Near pKa → titratable; far → fixed state
                if 0.1 <= probability <= 0.9:
                    recommendations[key] = amber_residue_map[res_type]["titratable"]
                elif probability > 0.9:
                    recommendations[key] = amber_residue_map[res_type]["protonated"]
                else:
                    recommendations[key] = amber_residue_map[res_type]["deprotonated"]
            elif constant_ph_simulation and res_type == "HIS":
                recommendations[key] = "HIP"
            elif res_type == "HIS":
                recommendations[key] = "HIP" if is_protonated else "HIE"
            else:
                state = "protonated" if is_protonated else "deprotonated"
                recommendations[key] = amber_residue_map[res_type][state]

        return recommendations

    def set_md_residue_names(self, interactive=True, constant_ph_simulation=False, console=None):
        """
        Interactively set MD residue names for titratable residues.

        Args:
            interactive: If True, prompt user for each residue
            constant_ph_simulation: If True, use constant pH naming conventions
            console: Rich console (optional, will create if not provided)

        Returns:
            Dictionary mapping residue keys to MD residue names
        """
        from rich.prompt import Prompt
        from rich.table import Table
        from rich.console import Console

        if console is None:
            console = Console()
        recommendations = {}

        # Define AMBER residue naming conventions
        amber_residue_map = {
            "ASP": {
                "protonated": "ASH",
                "deprotonated": "ASP",
                "titratable": "AS4",  # For constant pH
            },
            "GLU": {
                "protonated": "GLH",
                "deprotonated": "GLU",
                "titratable": "GL4",  # For constant pH
            },
            "HIS": {
                "protonated": "HIP",
                "deprotonated_ND1": "HID",  # Delta nitrogen protonated
                "deprotonated_NE2": "HIE",  # Epsilon nitrogen protonated
                "titratable": "HIP",  # For constant pH - will be titrated
            },
            "CYS": {
                "protonated": "CYS",
                "deprotonated": "CYM",
                "titratable": "CYS",  # For constant pH - same name
            },
            "TYR": {
                "protonated": "TYR",
                "deprotonated": "TYM",
                "titratable": "TYR",  # For constant pH - same name
            },
            "LYS": {
                "protonated": "LYS",
                "deprotonated": "LYN",
                "titratable": "LYS",  # For constant pH - same name
            },
            "PRN": {
                "protonated": "PRN",
                "deprotonated": "PRN",  # Assuming PRN has only one form
                "titratable": "PRN",  # For constant pH - same name
            },
        }
        
        if not self.results:
            console.print("[yellow]No protonation analysis results available.[/yellow]")
            return recommendations
        
        # Filter out non-titratable residues for constant pH simulations
        # AMBER CpHMD supports: AS4 (ASP), GL4 (GLU), CYS, TYR, HIP (HIS), LYS, PRN (heme propionate)
        # N_TERM and C_TERM are NOT supported
        # See: cpinutil.py --describe
        amber_cphmd_titratable = {"ASP", "GLU", "HIS", "TYR", "LYS", "CYS", "PRN"}

        working_results = self.results.copy()
        if constant_ph_simulation:
            # Exclude residues that cannot be titrated in AMBER CpHMD
            # (N_TERM, C_TERM, ARG, CYS, and other non-titratable residues)
            working_results = {
                key: res for key, res in self.results.items()
                if res["type"] in amber_cphmd_titratable
            }

            excluded_count = len(self.results) - len(working_results)
            if excluded_count > 0:
                excluded_types = set(res["type"] for res in self.results.values() if res["type"] not in amber_cphmd_titratable)
                console.print(f"[yellow]Note: Excluding {excluded_count} non-titratable residues from constant pH simulation[/yellow]")
                console.print(f"[grey50]  Excluded residue types: {', '.join(sorted(excluded_types))}[/grey50]")
        
        # Route to appropriate workflow
        if constant_ph_simulation:
            console.print("\n[bold cyan]Constant pH MD Simulation Mode[/bold cyan]")
            console.print("Note: Terminal residues are excluded from constant pH simulations.")
            console.print()

            # Use new constant pH workflow
            return self._handle_constant_ph_workflow(working_results, amber_residue_map, console)
        else:
            # Use existing standard MD workflow
            console.print("\n[bold cyan]Standard MD Simulation Mode[/bold cyan]")
            console.print("Residues will have fixed protonation states throughout the simulation.")
            console.print()

        # For standard MD: Exclude titratable residues at terminal positions
        # AMBER does not have terminal variants for alternate protonation states
        # (e.g., no NLYN, CLYN, NASH, CASH, NGLH, CGLH, NCYM, CCYM)
        terminal_positions = set()
        for key, res in working_results.items():
            if res["type"] in ("N_TERM", "C_TERM"):
                terminal_positions.add((res["chain"], res["number"]))

        # Filter out titratable residues that are at terminal positions
        # Also filter out N_TERM and C_TERM entries (these are not actual residue names)
        excluded_terminal_titratables = []
        filtered_results = {}
        for key, res in working_results.items():
            res_pos = (res["chain"], res["number"])
            res_type = res["type"]

            # Exclude N_TERM and C_TERM (not actual residue names - just markers)
            if res_type in ("N_TERM", "C_TERM"):
                continue

            # Check if this is a titratable residue at a terminal position
            if res_type in amber_residue_map and res_pos in terminal_positions:
                # This titratable residue is at a terminal - cannot use alternate protonation names
                excluded_terminal_titratables.append((res["chain"], res["number"], res_type))
            else:
                filtered_results[key] = res

        # Update working results
        working_results = filtered_results

        # Warn about excluded terminal titratable residues
        if excluded_terminal_titratables:
            console.print(
                f"[yellow]Note: Excluding {len(excluded_terminal_titratables)} titratable residue(s) "
                f"at terminal positions[/yellow]"
            )
            console.print(
                "[grey50]  AMBER does not support alternate protonation state names (LYN, ASH, GLH, CYM) "
                "for terminal residues.[/grey50]"
            )
            console.print("[grey50]  These residues will retain their standard terminal names (e.g., NLYS, CLYS).[/grey50]")
            if len(excluded_terminal_titratables) <= 5:
                for chain, resnum, restype in excluded_terminal_titratables:
                    console.print(f"[grey50]    - {restype} {chain}:{resnum}[/grey50]")
            console.print()

        # Group residues by type for organized display
        residues_by_type = {}
        for key, res in working_results.items():
            res_type = res["type"]
            if res_type not in residues_by_type:
                residues_by_type[res_type] = []
            residues_by_type[res_type].append((key, res))
        
        # Process each residue type
        for res_type in sorted(residues_by_type.keys()):
            if res_type not in amber_residue_map:
                # Non-titratable residue, keep original name
                for key, res in residues_by_type[res_type]:
                    recommendations[key] = res_type
                continue
            
            residues = residues_by_type[res_type]
            
            # Create a table for this residue type
            table = Table(title=f"\n{res_type} Residues")
            table.add_column("Chain", style="cyan")
            table.add_column("ResID", style="green")
            table.add_column("pKa", style="yellow")
            table.add_column("State at pH {:.1f}".format(self.pH), style="magenta")
            table.add_column("Probability", style="blue")
            table.add_column("Recommended", style="bold green")
            
            # Determine recommendations for each residue
            default_recommendations = {}
            for key, res in residues:
                chain = res["chain"]
                resid = res["number"]
                pka = res["pKa"]
                is_protonated = res["protonated"]
                probability = res["probability"]
                
                # Determine default recommendation based on protonation analysis
                if constant_ph_simulation and res_type in ["ASP", "GLU", "HIS"]:
                    # For constant pH, recommend based on proximity to pKa
                    # Residues near their pKa (probability 0.1-0.9) are good candidates for titration
                    if 0.1 <= probability <= 0.9:
                        # Near pKa - recommend titratable
                        default_name = amber_residue_map[res_type]["titratable"]
                    else:
                        # Far from pKa - consider fixed state
                        # But still default to titratable with a note
                        if res_type in ["ASP", "GLU"]:
                            if probability > 0.9:
                                # Strongly protonated - could use ASH/GLH
                                default_name = amber_residue_map[res_type]["protonated"]
                            else:
                                # Strongly deprotonated - could use ASP/GLU
                                default_name = amber_residue_map[res_type]["deprotonated"]
                        else:  # HIS
                            # For HIS, always recommend HIP for CpHMD regardless of state
                            default_name = "HIP"
                elif res_type == "HIS":
                    # For HIS in standard MD, use protonation state
                    if is_protonated:
                        default_name = "HIP"
                    else:
                        # Default to HIE for neutral HIS (epsilon protonated is more common)
                        default_name = "HIE"
                else:
                    # For other residues, use protonation state
                    state = "protonated" if is_protonated else "deprotonated"
                    default_name = amber_residue_map[res_type][state]
                
                default_recommendations[key] = default_name
                
                state_str = "Protonated" if is_protonated else "Deprotonated"
                
                table.add_row(
                    chain,
                    str(resid),
                    f"{pka:.2f}",
                    state_str,
                    f"{probability:.3f}",
                    default_name
                )
            
            console.print(table)

            # Per-row reps in the viewer, one ball+stick per residue in
            # this table palette-coloured by row index. The user can
            # toggle individual residues in the rep manager — useful
            # for inspecting a single residue's environment without
            # the others crowding the view. Replaces the prior
            # residue-type's reps each iteration.
            self._show_md_name_residue_reps(res_type, residues)

            if interactive:
                # Ask user about this group of residues
                if constant_ph_simulation and res_type in ["ASP", "GLU", "HIS"]:
                    # Special handling for titratable residues in constant pH
                    choices_str = self._get_cphmd_choices(res_type, amber_residue_map)
                    console.print(f"\n[bold]Options for {res_type} residues:[/bold]")
                    console.print(choices_str)

                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n", "individual"],
                        default="y",
                        module="Protonation Worker",
                        description=f"Use defaults for {res_type} (const-pH)",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize", "individual": "Individual selection"}
                    )
                elif res_type == "ASP" and not constant_ph_simulation:
                    # Special handling for ASP protonation states in standard MD
                    console.print(f"\n[bold]Options for ASP residues:[/bold]")
                    console.print("• [green]ASP[/green] - Deprotonated carboxylate (COO⁻, -1 charge)")
                    console.print("• [green]ASH[/green] - Protonated carboxylic acid (COOH, neutral)")
                    console.print("[grey50]Note: ASP is typically deprotonated at physiological pH[/grey50]")

                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n", "individual"],
                        default="y",
                        module="Protonation Worker",
                        description="Use defaults for ASP residues",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize", "individual": "Individual selection"}
                    )
                elif res_type == "GLU" and not constant_ph_simulation:
                    # Special handling for GLU protonation states in standard MD
                    console.print(f"\n[bold]Options for GLU residues:[/bold]")
                    console.print("• [green]GLU[/green] - Deprotonated carboxylate (COO⁻, -1 charge)")
                    console.print("• [green]GLH[/green] - Protonated carboxylic acid (COOH, neutral)")
                    console.print("[grey50]Note: GLU is typically deprotonated at physiological pH[/grey50]")

                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n", "individual"],
                        default="y",
                        module="Protonation Worker",
                        description="Use defaults for GLU residues",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize", "individual": "Individual selection"}
                    )
                elif res_type == "HIS" and not constant_ph_simulation:
                    # Special handling for HIS tautomers in standard MD
                    console.print(f"\n[bold]Options for HIS residues:[/bold]")
                    console.print("• [green]HIP[/green] - Fully protonated (both ND1 and NE2)")
                    console.print("• [green]HIE[/green] - Epsilon nitrogen protonated (NE2)")
                    console.print("• [green]HID[/green] - Delta nitrogen protonated (ND1)")
                    console.print("[grey50]Note: HIE is generally more common than HID[/grey50]")

                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n", "individual"],
                        default="y",
                        module="Protonation Worker",
                        description="Use defaults for HIS residues",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize", "individual": "Individual selection"}
                    )
                elif res_type == "CYS" and not constant_ph_simulation:
                    # Special handling for CYS protonation states in standard MD
                    console.print(f"\n[bold]Options for CYS residues:[/bold]")
                    console.print("• [green]CYS[/green] - Protonated thiol (SH, neutral)")
                    console.print("• [green]CYM[/green] - Deprotonated thiolate (S⁻, -1 charge)")
                    console.print("[grey50]Note: CYS is typically protonated unless in special environments[/grey50]")

                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n", "individual"],
                        default="y",
                        module="Protonation Worker",
                        description="Use defaults for CYS residues",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize", "individual": "Individual selection"}
                    )
                elif res_type == "LYS" and not constant_ph_simulation:
                    # Special handling for LYS protonation states in standard MD
                    console.print(f"\n[bold]Options for LYS residues:[/bold]")
                    console.print("• [green]LYS[/green] - Protonated amine (NH₃⁺, +1 charge)")
                    console.print("• [green]LYN[/green] - Deprotonated amine (NH₂, neutral)")
                    console.print("[grey50]Note: LYS is typically protonated at physiological pH[/grey50]")

                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n", "individual"],
                        default="y",
                        module="Protonation Worker",
                        description="Use defaults for LYS residues",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize", "individual": "Individual selection"}
                    )
                else:
                    use_defaults = prompt_with_context(
                        processor=self.processor,
                        prompt=f"\nUse recommended names for all {res_type} residues?",
                        choices=["y", "n"],
                        default="y",
                        module="Protonation Worker",
                        description=f"Use defaults for {res_type} residues",
                        options_map={"y": "Yes - use all defaults", "n": "No - customize"}
                    )
                
                if use_defaults == "y":
                    # Use all default recommendations
                    for key in default_recommendations:
                        recommendations[key] = default_recommendations[key]
                elif use_defaults == "individual":
                    # Ask for each residue individually
                    for key, res in residues:
                        chain = res["chain"]
                        resid = res["number"]
                        default = default_recommendations[key]
                        
                        # Get available choices for this residue type
                        choices = self._get_residue_choices(
                            res_type, 
                            amber_residue_map, 
                            constant_ph_simulation
                        )
                        
                        console.print(f"\n[bold]Chain {chain}, Residue {resid} ({res_type})[/bold]")
                        console.print(f"pKa: {res['pKa']:.2f}, State at pH {self.pH:.1f}: {'Protonated' if res['protonated'] else 'Deprotonated'}")

                        # Focus the viewer on this residue + its 8Å
                        # neighborhood so the user can decide tautomer
                        # (e.g. HIE vs HID) by inspecting visible
                        # H-bond partners. Single shared label that
                        # replaces on each iteration.
                        self._focus_md_name_residue(chain, resid)

                        # Create options_map from choices
                        options_map = {choice: choice for choice in choices}

                        chosen = prompt_with_context(
                            processor=self.processor,
                            prompt=f"Choose name",
                            choices=choices,
                            default=default,
                            module="Protonation Worker",
                            description=f"Choose name for {res_type} {chain}:{resid}",
                            options_map=options_map
                        )
                        recommendations[key] = chosen
                else:
                    # Allow flexible per-residue customization
                    choices = self._get_residue_choices(
                        res_type,
                        amber_residue_map,
                        constant_ph_simulation
                    )

                    # Build resid -> key mapping for this residue type
                    resid_to_key = {}
                    for key, res in residues:
                        resid = res["number"]
                        chain = res["chain"]
                        # Store with chain prefix for disambiguation
                        resid_to_key[(chain, resid)] = key
                        # Also collect bare-resid keys so a single number
                        # like `110:CYM` applies to every chain that has
                        # that residue (single chain → str, multi-chain → list).
                        existing = resid_to_key.get(resid)
                        if existing is None:
                            resid_to_key[resid] = key
                        elif isinstance(existing, str):
                            resid_to_key[resid] = [existing, key]
                        else:
                            existing.append(key)

                    console.print(f"\n[grey50]Syntax: Enter a name to apply to all, or specify overrides:[/grey50]")
                    console.print(f"[grey50]  • Single name: [green]{choices[0]}[/green] (applies to all)[/grey50]")
                    console.print(f"[grey50]  • Overrides: [green]94,119:{choices[-1]} 122:{choices[0]}[/green] (applies to listed resids in every chain)[/grey50]")
                    console.print(f"[grey50]  • Chain-specific: [green]A:94:{choices[-1]}[/green] (restrict to one chain)[/grey50]")
                    console.print(f"[grey50]  • Press Enter to use all recommendations[/grey50]")

                    user_input = prompt_with_context(
                        processor=self.processor,
                        prompt=f"Customize {res_type} names",
                        default="",
                        module="Protonation Worker",
                        description=f"Customize {res_type} residue names"
                    ).strip().upper()

                    # Start with default recommendations
                    for key in default_recommendations:
                        recommendations[key] = default_recommendations[key]

                    if user_input:
                        # Parse the input
                        parsed = self._parse_residue_name_overrides(
                            user_input, choices, resid_to_key, console
                        )
                        # Apply overrides
                        for key, name in parsed.items():
                            recommendations[key] = name
            else:
                # Non-interactive mode - use all defaults
                for key in default_recommendations:
                    recommendations[key] = default_recommendations[key]
        
        # Summary
        console.print(f"\n[green]Set MD names for {len(recommendations)} residues[/green]")
        
        # Display summary of assignments if interactive
        if interactive and recommendations:
            summary = {}
            for key, name in recommendations.items():
                if name not in summary:
                    summary[name] = 0
                summary[name] += 1
            
            console.print("\n[bold]Summary of assignments:[/bold]")
            for name, count in sorted(summary.items()):
                console.print(f"  {name}: {count} residues")

        # Final cleanup: drop any per-residue reps that survived the
        # last residue type's table (and any focused per-individual
        # picker hook). Subsequent stages start with a clean viewer.
        self._clear_md_name_residue_reps()
        self._clear_md_name_focus()

        return recommendations

    # ----- viewer hooks for set_md_residue_names ---------------------

    def _show_md_name_residue_reps(self, res_type, residues):
        """Add one ball+stick rep per residue in the just-displayed table.

        Replaces the prior residue-type's reps so each table swap
        leaves only the current type on screen. Palette-coloured by
        row index for easy mapping back to table rows. Best-effort —
        a viewer failure here must not block the picker.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
        except Exception as exc:
            logger.debug("md-name reps viewer hook silenced: %s", exc)
            return

        # Make sure the viewer points at the structure that produced
        # this analysis — earlier modules may have swapped it.
        if getattr(self, "input_file", None):
            try:
                _viewer.show_structure(self.input_file)
            except Exception as exc:
                logger.debug("md-name show_structure silenced: %s", exc)

        self._clear_md_name_residue_reps()
        applied: List[str] = []
        for idx, (key, res) in enumerate(residues, 1):
            chain = res.get("chain")
            resid = res.get("number")
            if chain is None or resid is None:
                continue
            label = f"prot_md_{res_type}_{idx}_{chain}_{resid}"
            try:
                _viewer.highlight(
                    f":{chain} and {resid}",
                    style="ball+stick",
                    color=f"palette:{idx}",
                    label=label,
                )
                applied.append(label)
            except Exception as exc:
                logger.debug("md-name highlight silenced: %s", exc)
        self._md_name_residue_labels = applied

    def _clear_md_name_residue_reps(self) -> None:
        """Remove the per-residue reps from the most recent table."""
        labels = getattr(self, "_md_name_residue_labels", None) or []
        if not labels:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            for lbl in labels:
                _viewer.unhighlight(lbl)
        except Exception as exc:
            logger.debug("md-name reps cleanup silenced: %s", exc)
        self._md_name_residue_labels = []

    def _clear_md_name_focus(self) -> None:
        """Remove the focused per-individual-residue rep, if any."""
        label = getattr(self, "_md_name_focus_label", None)
        if not label:
            return
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
            _viewer.unhighlight(label)
            _viewer.unhighlight(f"{label}_neighbors")
        except Exception as exc:
            logger.debug("md-name focus cleanup silenced: %s", exc)
        self._md_name_focus_label = None

    def _focus_md_name_residue(self, chain: str, resid: int, radius: float = 8.0) -> None:
        """Focus the viewer on a single titratable residue + its neighborhood.

        Highlights the residue itself (magenta) and its surrounding
        residues within ``radius`` Å (element colours) and switches
        the viewer to focused-mode so the camera centres on the site.
        Single shared label per call — re-firing replaces the prior
        focus rather than accumulating.
        """
        try:
            from proprep.structure_prep.viewer_coordinator import (
                viewer as _viewer,
            )
        except Exception as exc:
            logger.debug("md-name focus viewer hook silenced: %s", exc)
            return

        neighbors = self._compute_neighbor_residues(chain, resid, radius)
        # Drop the prior focus before drawing the new one so the
        # camera switch is clean.
        self._clear_md_name_focus()

        center_label = f"prot_md_focus_{chain}_{resid}"
        try:
            _viewer.highlight(
                f":{chain} and {resid}",
                style="ball+stick",
                # Hex, not "magenta": NGL treats a bare color word as a
                # (nonexistent) scheme id and silently drops the rep, so the
                # centered residue would never render. #ff00ff == magenta.
                color="#ff00ff",
                label=center_label,
                focused=True,
            )
        except Exception as exc:
            logger.debug("md-name focus center silenced: %s", exc)
            return
        self._md_name_focus_label = center_label

        if neighbors:
            clauses = [f"(:{c} and {r})" for c, r in sorted(neighbors)]
            try:
                _viewer.highlight(
                    " or ".join(clauses),
                    style="ball+stick",
                    color="element",
                    label=f"{center_label}_neighbors",
                )
            except Exception as exc:
                logger.debug("md-name focus neighbors silenced: %s", exc)

    def _compute_neighbor_residues(self, chain: str, resid: int, radius: float):
        """Return the set of (chain, resid) pairs within ``radius`` Å of
        the named residue. Returns an empty set on any failure — the
        caller is responsible for graceful degradation.
        """
        try:
            from Bio.PDB import PDBParser
            from Bio.PDB.NeighborSearch import NeighborSearch
        except Exception as exc:
            logger.debug("NeighborSearch import silenced: %s", exc)
            return set()

        structure = self.structure
        if structure is None and getattr(self, "input_file", None):
            try:
                structure = PDBParser(QUIET=True).get_structure(
                    "protein", self.input_file
                )
            except Exception as exc:
                logger.debug("NeighborSearch parse silenced: %s", exc)
                return set()
        if structure is None:
            return set()

        try:
            model = next(iter(structure.get_models()))
        except StopIteration:
            return set()

        if chain not in model:
            return set()
        target = None
        for residue in model[chain]:
            if residue.id[1] == resid:
                target = residue
                break
        if target is None:
            return set()

        ns = NeighborSearch(list(model.get_atoms()))
        neighbors = set()
        for atom in target.get_atoms():
            for hit in ns.search(atom.coord, radius):
                hit_res = hit.get_parent()
                hit_chain = hit_res.get_parent().id
                hit_num = hit_res.id[1]
                if (hit_chain, hit_num) == (chain, resid):
                    continue
                neighbors.add((hit_chain, hit_num))
        return neighbors

    def _get_residue_choices(self, res_type, amber_map, constant_ph):
        """Get valid residue name choices for a given residue type."""
        if res_type not in amber_map:
            return [res_type]
        
        choices = []
        mapping = amber_map[res_type]
        
        if constant_ph and res_type in ["ASP", "GLU", "HIS"]:
            # For constant pH, include titratable option
            if "titratable" in mapping:
                choices.append(mapping["titratable"])
        
        # Add standard options
        if "protonated" in mapping:
            choices.append(mapping["protonated"])
        if "deprotonated" in mapping:
            choices.append(mapping["deprotonated"])
        if "deprotonated_ND1" in mapping:  # For HIS
            choices.append(mapping["deprotonated_ND1"])
        if "deprotonated_NE2" in mapping:  # For HIS
            choices.append(mapping["deprotonated_NE2"])
        
        # Remove duplicates while preserving order
        seen = set()
        unique_choices = []
        for choice in choices:
            if choice not in seen:
                seen.add(choice)
                unique_choices.append(choice)
        
        return unique_choices
    
    def _get_cphmd_choices(self, res_type, amber_map):
        """Get descriptive text for constant pH choices."""
        if res_type == "ASP":
            return ("• [green]AS4[/green] - Titratable (recommended for CpHMD)\n"
                   "• [yellow]ASP[/yellow] - Fixed deprotonated\n"
                   "• [yellow]ASH[/yellow] - Fixed protonated")
        elif res_type == "GLU":
            return ("• [green]GL4[/green] - Titratable (recommended for CpHMD)\n"
                   "• [yellow]GLU[/yellow] - Fixed deprotonated\n"
                   "• [yellow]GLH[/yellow] - Fixed protonated")
        elif res_type == "HIS":
            return ("• [green]HIP[/green] - Titratable (recommended for CpHMD)\n"
                   "• [yellow]HIE[/yellow] - Fixed epsilon-protonated (not titrated)\n"
                   "• [yellow]HID[/yellow] - Fixed delta-protonated (not titrated)")
        return ""

    def generate_ph_profile(self, pH_range=(0, 14), step=0.5):
        """Generate pH titration profile for the analyzed residues."""
        if not self.titratable_residues:
            logger.warning("No titratable residues found for pH profile generation")
            return {"pH": [], "net_charge": []}

        current_pH = self.pH
        min_pH, max_pH = pH_range
        pH_values = np.arange(min_pH, max_pH + step, step)

        results = {"pH": [], "net_charge": []}

        for pH in pH_values:
            self.pH = pH
            net_charge = 0.0

            for key, residue in self.titratable_residues.items():
                pka = residue["pKa"]
                if self.custom_pkas and key in self.custom_pkas:
                    pka = self.custom_pkas[key]

                probability = self._calculate_protonation_probability(pH, pka)

                expected_charge = (
                    probability * CHARGE_STATES[residue["type"]]["protonated"]
                    + (1 - probability) * CHARGE_STATES[residue["type"]]["deprotonated"]
                )

                net_charge += expected_charge

            results["pH"].append(pH)
            results["net_charge"].append(net_charge)

        self.pH = current_pH

        logger.debug(f"Generated pH profile with {len(pH_values)} points")
        return results

    def save_ph_profile(self, ph_profile, output_file):
        """Save pH profile to CSV file."""
        try:
            import csv

            with open(output_file, "w", newline="") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["pH", "Net_Charge"])

                for pH, charge in zip(ph_profile["pH"], ph_profile["net_charge"]):
                    writer.writerow([f"{pH:.2f}", f"{charge:.4f}"])

            logger.debug(f"pH profile saved to {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"Error saving pH profile: {str(e)}")
            return None

    def save_results(self, output_file, output_format="text"):
        """Save analysis results to file."""
        try:
            if output_format == "text":
                with open(output_file, "w") as f:
                    f.write(f"Protonation State Analysis Results\n")
                    f.write(f"pH: {self.pH}\n")
                    f.write(f"Net Charge: {self.net_charge:.2f}\n")
                    f.write(f"Threshold: {self.threshold}\n\n")

                    f.write(
                        "Residue\tChain\tNumber\tType\tpKa\tProtonated\tProbability\tCharge\n"
                    )
                    for key, res in self.results.items():
                        f.write(
                            f"{key}\t{res['chain']}\t{res['number']}\t{res['type']}\t"
                            f"{res['pKa']:.2f}\t{res['protonated']}\t{res['probability']:.4f}\t"
                            f"{res['charge']:.2f}\n"
                        )

            elif output_format == "csv":
                import csv

                with open(output_file, "w", newline="") as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(
                        [
                            "Residue_Key",
                            "Chain",
                            "Number",
                            "Type",
                            "pKa",
                            "Protonated",
                            "Probability",
                            "Charge",
                        ]
                    )
                    for key, res in self.results.items():
                        writer.writerow(
                            [
                                key,
                                res["chain"],
                                res["number"],
                                res["type"],
                                res["pKa"],
                                res["protonated"],
                                res["probability"],
                                res["charge"],
                            ]
                        )

            logger.info(f"Results saved to {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            return None

    def _handle_constant_ph_workflow(self, working_results, amber_residue_map, console):
        """
        Handle the redesigned constant pH workflow with titration selection.
        
        Returns:
            dict: MD residue names with titration data stored in self.constant_ph_data
        """
        from rich.table import Table
        from rich.prompt import Prompt
        
        recommendations = {}
        
        # Initialize constant pH data storage
        self.constant_ph_data = {
            'enabled': True,
            'titratable_selections': {},
            'resnames': [],
            'resnums': {}
        }
        
        # Group residues by type
        residues_by_type = {}
        for key, res in working_results.items():
            res_type = res["type"]
            if res_type not in residues_by_type:
                residues_by_type[res_type] = []
            residues_by_type[res_type].append((key, res))
        
        # Define titratable types for AMBER constant pH MD (7 types supported)
        # AS4 (ASP), GL4 (GLU), HIS (HIP), CYS, TYR, LYS, PRN (heme propionate)
        # See: cpinutil.py --describe
        titratable_types = ["ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "PRN"]
        
        # Process titratable types with new interface
        for res_type in titratable_types:
            if res_type in residues_by_type:
                residue_list = residues_by_type[res_type]
                titration_selections, type_recommendations = self._handle_titratable_type_selection(
                    res_type, residue_list, amber_residue_map, console
                )
                
                # Store titration data
                if titration_selections:
                    titratable_name = amber_residue_map[res_type]["titratable"]
                    self.constant_ph_data['titratable_selections'][titratable_name] = titration_selections
                    self.constant_ph_data['resnames'].append(titratable_name)
                    self.constant_ph_data['resnums'][titratable_name] = titration_selections
                
                # Merge recommendations
                recommendations.update(type_recommendations)
        
        # Process non-titratable types with standard group interface
        non_titratable_types = [res_type for res_type in residues_by_type.keys() 
                               if res_type not in titratable_types]
        
        for res_type in sorted(non_titratable_types):
            residue_list = residues_by_type[res_type]
            type_recommendations = self._handle_non_titratable_type(
                res_type, residue_list, amber_residue_map, console
            )
            recommendations.update(type_recommendations)
        
        # Display final summary
        self._display_constant_ph_summary(recommendations, console)
        
        return recommendations

    def _handle_titratable_type_selection(self, res_type, residue_list, amber_residue_map, console):
        """
        Handle titration selection for a titratable residue type (ASP, GLU, HIS).

        Returns:
            tuple: (titration_selections_list, recommendations_dict)
        """
        from rich.table import Table
        from rich.prompt import Prompt

        console.print(f"\n[bold]{res_type} Residues:[/bold]")
        
        # Create indexed table
        table = Table(show_header=True)
        table.add_column("Index", style="cyan", width=8)
        table.add_column("Chain", style="blue", width=8)
        table.add_column("Residue", style="white", width=10)
        table.add_column("pKa", style="yellow", width=8)
        table.add_column("State at pH 7.0", style="green", width=18)
        
        indexed_residues = []
        for i, (key, res) in enumerate(residue_list, 1):
            chain = res["chain"]
            resid = res["number"]
            pka = res["pKa"]
            state = "Protonated" if res["protonated"] else "Deprotonated"
            
            table.add_row(str(i), chain, str(resid), f"{pka:.1f}", state)
            indexed_residues.append((key, res))
        
        console.print(table)
        
        # Get titration selections
        titratable_name = amber_residue_map[res_type]["titratable"]
        console.print(f"\n[bold]Select {res_type} residues to TITRATE[/bold] (will use {titratable_name} name):")
        console.print("Enter indices (e.g., 1,3 or 1-3), 'all', or 'none':")

        selection = prompt_with_context(
            processor=self.processor,
            prompt="Indices to titrate",
            default="",
            module="Protonation Worker",
            description=f"Select {res_type} indices to titrate (const-pH)"
        ).strip().lower()

        # Parse selection
        selected_indices = []
        if selection == "all":
            selected_indices = list(range(1, len(indexed_residues) + 1))
        elif selection == "none" or selection == "":
            selected_indices = []
        elif selection:
            try:
                selected_indices = self._parse_index_selection(selection, len(indexed_residues))
            except ValueError as e:
                console.print(f"[red]Error: {e}[/red]")
                console.print("[yellow]No residues selected for titration[/yellow]")
                selected_indices = []
        
        # Build results
        titration_selections = []
        recommendations = {}
        
        # First, handle selected (titratable) residues
        for i, (key, res) in enumerate(indexed_residues, 1):
            resid = res["number"]
            
            if i in selected_indices:
                # Titratable - automatic name
                recommendations[key] = titratable_name
                titration_selections.append(resid)
                console.print(f"[green]  {res['chain']}:{resid} → {titratable_name} (titratable)[/green]")
        
        # Then, handle non-selected (non-titratable) residues with explanation
        non_selected_residues = [
            (i, key, res) for i, (key, res) in enumerate(indexed_residues, 1)
            if i not in selected_indices
        ]
        
        if non_selected_residues:
            self._display_non_titratable_explanation(res_type, amber_residue_map, console)

            # Compute defaults for the group
            default_names = {}
            for i, key, res in non_selected_residues:
                mapping = amber_residue_map[res_type]
                is_protonated = res["protonated"]
                if res_type in ["ASP", "GLU"]:
                    default_names[key] = mapping["protonated"] if is_protonated else mapping["deprotonated"]
                elif res_type == "HIS":
                    default_names[key] = "HIP" if is_protonated else "HIE"
                elif res_type in ["CYS", "TYR", "LYS"]:
                    default_names[key] = mapping["protonated"] if is_protonated else mapping["deprotonated"]
                else:
                    default_names[key] = mapping.get("protonated", res_type)

            use_defaults = prompt_with_context(
                processor=self.processor,
                prompt=f"Use recommended names for all non-titratable {res_type} residues?",
                choices=["y", "n", "individual"],
                default="y",
                module="Protonation Worker",
                description=f"Use defaults for non-titratable {res_type} (const-pH)",
                options_map={"y": "Yes - use all defaults", "n": "No - set one name for all", "individual": "Individual selection"}
            )

            if use_defaults == "y":
                recommendations.update(default_names)
            elif use_defaults == "n":
                # Let user pick a single name for the whole group
                mapping = amber_residue_map[res_type]
                if res_type in ["ASP", "GLU"]:
                    choices = [mapping["deprotonated"], mapping["protonated"]]
                elif res_type == "HIS":
                    choices = [mapping["deprotonated_NE2"], mapping["deprotonated_ND1"], mapping["protonated"]]
                elif res_type in ["CYS", "TYR", "LYS"]:
                    choices = [mapping["protonated"], mapping["deprotonated"]]
                else:
                    choices = [mapping.get("protonated", res_type)]

                group_name = prompt_with_context(
                    processor=self.processor,
                    prompt=f"Name for all non-titratable {res_type} residues",
                    choices=choices,
                    default=choices[0],
                    module="Protonation Worker",
                    description=f"Group name for non-titratable {res_type}"
                )
                for key in default_names:
                    recommendations[key] = group_name
            else:
                # Individual selection
                for i, key, res in non_selected_residues:
                    non_titratable_name = self._get_non_titratable_name(
                        res_type, res, amber_residue_map, console
                    )
                    recommendations[key] = non_titratable_name
        
        if titration_selections:
            console.print(f"[cyan]Selected {len(titration_selections)} {res_type} residues for titration[/cyan]")
        
        return titration_selections, recommendations

    def _display_non_titratable_explanation(self, res_type, amber_residue_map, console):
        """Display explanation for non-titratable residue naming."""
        console.print(f"\n[bold yellow]Setting names for non-titratable {res_type} residues:[/bold yellow]")
        
        if res_type == "ASP":
            console.print("• [green]ASP[/green] = Deprotonated (fixed -1 charge)")
            console.print("• [green]ASH[/green] = Protonated (fixed 0 charge)")
        elif res_type == "GLU":
            console.print("• [green]GLU[/green] = Deprotonated (fixed -1 charge)")
            console.print("• [green]GLH[/green] = Protonated (fixed 0 charge)")
        elif res_type == "HIS":
            console.print("• [green]HIE[/green] = Epsilon nitrogen protonated (NE2-H)")
            console.print("• [green]HID[/green] = Delta nitrogen protonated (ND1-H)")
            console.print("• [green]HIP[/green] = Fully protonated (both ND1-H and NE2-H)")
        elif res_type == "CYS":
            console.print("• [green]CYS[/green] = Protonated thiol (fixed 0 charge)")
            console.print("• [green]CYM[/green] = Deprotonated thiolate (fixed -1 charge)")
        elif res_type == "TYR":
            console.print("• [green]TYR[/green] = Protonated phenol (fixed 0 charge)")
            console.print("• [green]TYM[/green] = Deprotonated phenolate (fixed -1 charge)")
        elif res_type == "LYS":
            console.print("• [green]LYS[/green] = Protonated amine (fixed +1 charge)")
            console.print("• [green]LYN[/green] = Deprotonated amine (fixed 0 charge)")
        elif res_type == "PRN":
            console.print("• [green]PRN[/green] = Heme propionic acid (fixed protonation state)")
        
        console.print("[grey50]These residues will maintain fixed protonation states during simulation.[/grey50]")
        console.print()

    def _get_non_titratable_name(self, res_type, res, amber_residue_map, console):
        """Get name choice for a non-titratable residue."""
        from rich.prompt import Prompt
        
        chain = res["chain"]
        resid = res["number"]
        is_protonated = res["protonated"]
        
        # Get available choices
        mapping = amber_residue_map[res_type]
        choices = []
        
        if res_type in ["ASP", "GLU"]:
            choices = [mapping["deprotonated"], mapping["protonated"]]
            default = mapping["deprotonated"] if not is_protonated else mapping["protonated"]
        elif res_type == "HIS":
            choices = [mapping["deprotonated_NE2"], mapping["deprotonated_ND1"], mapping["protonated"]]
            default = mapping["deprotonated_NE2"]  # HIE is more common
        elif res_type in ["CYS", "TYR"]:
            choices = [mapping["protonated"], mapping["deprotonated"]]
            default = mapping["protonated"] if is_protonated else mapping["deprotonated"]
        elif res_type == "LYS":
            choices = [mapping["protonated"], mapping["deprotonated"]]
            default = mapping["protonated"] if is_protonated else mapping["deprotonated"]
        elif res_type == "PRN":
            # PRN typically only has one form - just return it
            return mapping.get("protonated", res_type)
        else:
            # Other types - just return the protonated state
            return mapping.get("protonated", res_type)
        
        # Create options_map from choices
        options_map = {choice: choice for choice in choices}

        choice = prompt_with_context(
            processor=self.processor,
            prompt=f"Choose name for {chain}:{resid} ({res_type})",
            choices=choices,
            default=default,
            module="Protonation Worker",
            description=f"Choose name for non-titratable {res_type} {chain}:{resid}",
            options_map=options_map
        )

        return choice

    def _handle_non_titratable_type(self, res_type, residue_list, amber_residue_map, console):
        """Handle non-titratable residue types with group selection."""
        from rich.table import Table
        from rich.prompt import Prompt
        
        recommendations = {}
        
        # Create table
        table = Table(show_header=True)
        table.add_column("Chain", style="blue", width=8)
        table.add_column("Residue", style="white", width=10)
        table.add_column("State at pH 7.0", style="green", width=18)
        table.add_column("Recommended", style="yellow", width=15)
        
        for key, res in residue_list:
            chain = res["chain"]
            resid = res["number"]
            state = "Protonated" if res["protonated"] else "Deprotonated"
            recommended = res_type  # Non-titratable types keep their name
            
            table.add_row(chain, str(resid), state, recommended)
        
        console.print(f"\n[bold]{res_type} Residues (Cannot be titrated):[/bold]")
        console.print(table)

        # Group recommendation
        use_defaults = prompt_with_context(
            processor=self.processor,
            prompt=f"Use recommended names for all {res_type} residues?",
            choices=["y", "n"],
            default="y",
            module="Protonation Worker",
            description=f"Use defaults for non-titratable {res_type}",
            options_map={"y": "Yes - use defaults", "n": "No - customize"}
        )
        
        if use_defaults == "y":
            for key, res in residue_list:
                recommendations[key] = res_type
        else:
            # Individual choices (simplified for non-titratable)
            for key, res in residue_list:
                recommendations[key] = res_type  # They don't have alternatives
        
        return recommendations

    def _parse_index_selection(self, selection, max_index):
        """Parse comma-separated index selection with range support."""
        indices = set()
        
        for part in selection.split(','):
            part = part.strip()
            
            if '-' in part:
                # Handle range like "5-7"
                try:
                    start, end = part.split('-', 1)
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())
                    
                    if start_idx < 1 or end_idx > max_index or start_idx > end_idx:
                        raise ValueError(f"Invalid range {part} (valid range: 1-{max_index})")
                    
                    for i in range(start_idx, end_idx + 1):
                        indices.add(i)
                        
                except ValueError as ve:
                    raise ValueError(f"Invalid range format: {part}")
            else:
                # Handle single index
                try:
                    idx = int(part)
                    if idx < 1 or idx > max_index:
                        raise ValueError(f"Index {idx} out of range (valid range: 1-{max_index})")
                    indices.add(idx)
                except ValueError:
                    raise ValueError(f"Invalid index: {part}")
        
        return sorted(list(indices))

    def _parse_residue_name_overrides(self, user_input, valid_choices, resid_to_key, console):
        """
        Parse flexible residue name override syntax.

        Supports:
        - Single name: "HID" → applies to all residues
        - Overrides: "94,119:HID 122:HIE" → specific residues get specific names
        - Chain-specific: "A:94:HID" → for multi-chain disambiguation

        Args:
            user_input: User's input string (already uppercased)
            valid_choices: List of valid residue names
            resid_to_key: Dict mapping (chain, resid) or resid to key
            console: Console for output

        Returns:
            Dict of key -> residue name overrides
        """
        overrides = {}

        # Check if it's a single name (applies to all)
        if user_input in valid_choices:
            for lookup, key in resid_to_key.items():
                if isinstance(key, str):  # Only actual keys, not disambiguation markers
                    overrides[key] = user_input
            console.print(f"[green]Applying {user_input} to all residues[/green]")
            return overrides

        # Parse override syntax: "94,119:HID 122:HIE" or "A:94:HID"
        # Split by whitespace to get individual override groups
        groups = user_input.split()

        for group in groups:
            if ':' not in group:
                console.print(f"[yellow]Skipping invalid syntax: {group} (expected resid:NAME)[/yellow]")
                continue

            resid_part, name = group.rsplit(':', 1)
            name = name.strip()

            if name not in valid_choices:
                console.print(f"[yellow]Invalid name '{name}' - valid choices: {', '.join(valid_choices)}[/yellow]")
                continue

            # Parse resid_part: "94,119" or "A:94,A:119"
            for resid_spec in resid_part.split(','):
                resid_spec = resid_spec.strip()
                if not resid_spec:
                    continue

                # Check for chain prefix (e.g., "A:94")
                if ':' in resid_spec:
                    chain, resid_str = resid_spec.split(':', 1)
                    try:
                        resid = int(resid_str)
                        lookup_key = (chain, resid)
                    except ValueError:
                        console.print(f"[yellow]Invalid resid: {resid_spec}[/yellow]")
                        continue
                else:
                    try:
                        resid = int(resid_spec)
                        lookup_key = resid
                    except ValueError:
                        console.print(f"[yellow]Invalid resid: {resid_spec}[/yellow]")
                        continue

                # Look up the key
                if lookup_key in resid_to_key:
                    matched = resid_to_key[lookup_key]
                    if isinstance(matched, str):
                        overrides[matched] = name
                        console.print(f"[grey50]  {resid_spec} → {name}[/grey50]")
                    elif isinstance(matched, list):
                        for k in matched:
                            overrides[k] = name
                        chains = sorted({k.split("_", 1)[0] for k in matched})
                        console.print(
                            f"[grey50]  {resid_spec} → {name} (chains {', '.join(chains)})[/grey50]"
                        )
                    else:
                        console.print(f"[yellow]Ambiguous resid {resid_spec} - use chain prefix (e.g., A:{resid})[/yellow]")
                else:
                    console.print(f"[yellow]Resid {resid_spec} not found in this residue type[/yellow]")

        return overrides

    def _display_constant_ph_summary(self, recommendations, console):
        """Display final summary of constant pH selections."""
        from collections import defaultdict
        
        console.print("\n[bold]Summary of assignments:[/bold]")
        
        # Count by residue name
        name_counts = defaultdict(int)
        for name in recommendations.values():
            name_counts[name] += 1
        
        for name, count in sorted(name_counts.items()):
            console.print(f"  {name}: {count} residues")
        
        total = len(recommendations)
        console.print(f"Set names for {total} residues")
        
        # Show titration summary if any
        if hasattr(self, 'constant_ph_data') and self.constant_ph_data['titratable_selections']:
            console.print("\n[bold cyan]Titration Summary:[/bold cyan]")
            for titratable_name, residue_list in self.constant_ph_data['titratable_selections'].items():
                if residue_list:
                    console.print(f"  {titratable_name}: {len(residue_list)} residues will be titrated")
                    console.print(f"    Residues: {', '.join(map(str, residue_list))}")
        
        console.print()
