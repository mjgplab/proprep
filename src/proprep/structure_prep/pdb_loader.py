"""
PDB Metadata Extractor Module

Extracts and displays header records and metadata from PDB files.
This module can run as a standalone program or be integrated into other applications.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from proprep.utils.debug_utils import debug_workspace
from proprep.utils.module_registry import ProcessingModule, register_module, registry
from proprep.utils.prompts import prompt_with_context

# Configure logging
logger = logging.getLogger(__name__)


class PDBMetadataExtractor:
    """Class for extracting and displaying PDB file metadata."""

    def __init__(self, pdb_file=None):
        """Initialize the PDB metadata extractor."""
        self.pdb_file = pdb_file
        self.header_info = {}
        self.title_lines = []
        self.author_lines = []
        self.compnd_lines = []
        self.source_lines = []
        self.keywords = []
        self.experimental_method = ""
        self.resolution = None
        self.revisions = []
        self.remark_records = {}
        self.ssbond_records = []
        self.link_records = []
        self.seqres_records = {}
        self.modres_records = {}
        self.missing_res_records = {}
        self.helix_records = []
        self.sheet_records = []
        self.het_records = {}
        self.hetnam_records = {}
        self.formul_records = {}
        self.site_records = {}
        self.cryst_info = None
        self.origx_matrices = []
        self.scale_matrices = []

        # Safety/Status records
        self.caveat_records = []
        self.obslte_info = None
        self.sprsde_info = None
        self.split_info = None

        # Database cross-references and sequence info
        self.dbref_records = []
        self.seqadv_records = []
        self.cispep_records = []

        # Citation/Journal info
        self.jrnl_records = []

        # H++ structure detection flag
        self.is_hplusplus_structure = False

        if pdb_file:
            self.parse_pdb_file(pdb_file)

    def parse_pdb_file(self, pdb_file):
        """Parse a PDB file and extract metadata."""
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")

        logger.info(f"Loading PDB file: {pdb_file}")
        self.pdb_file = pdb_file

        # Parse the PDB file
        self._extract_header_records(pdb_file)
        self._extract_title_records(pdb_file)
        self._extract_author_records(pdb_file)
        self._extract_compnd_records(pdb_file)
        self._extract_source_records(pdb_file)
        self._extract_keyword_records(pdb_file)
        self._extract_expdta_records(pdb_file)
        self._extract_remark_records(pdb_file)
        self._extract_revdat_records(pdb_file)
        self._extract_ssbond_records(pdb_file)
        self._extract_link_records(pdb_file)
        self._extract_seqres_records(pdb_file)
        self._extract_modres_records(pdb_file)
        self._extract_missing_residue_records(pdb_file)
        self._extract_helix_sheet_records(pdb_file)
        self._extract_het_records(pdb_file)
        self._extract_hetnam_records(pdb_file)
        self._extract_formul_records(pdb_file)
        self._extract_site_records(pdb_file)
        self._extract_crystal_records(pdb_file)
        self._extract_origx_records(pdb_file)
        self._extract_scale_records(pdb_file)

        # Safety/Status records
        self._extract_caveat_records(pdb_file)
        self._extract_obslte_records(pdb_file)
        self._extract_sprsde_records(pdb_file)
        self._extract_split_records(pdb_file)

        # Database and sequence records
        self._extract_dbref_records(pdb_file)
        self._extract_seqadv_records(pdb_file)
        self._extract_cispep_records(pdb_file)

        # Citation records
        self._extract_jrnl_records(pdb_file)

        # Detect H++ structures
        self._detect_hplusplus_structure()

        # Log what we found
        logger.info(f"Extracted metadata from PDB file: {pdb_file}")
        if self.is_hplusplus_structure:
            logger.info("H++ structure detected (AMBER forcefield naming)")

        return self

    def _extract_header_records(self, pdb_file):
        """Extract HEADER record from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("HEADER"):
                        classification = line[10:50].strip()
                        deposition_date = line[50:59].strip()
                        pdb_id = line[62:66].strip()

                        self.header_info = {
                            "classification": classification,
                            "deposition_date": deposition_date,
                            "pdb_id": pdb_id,
                        }
                        break
        except Exception as e:
            logger.warning(f"Error extracting HEADER record: {str(e)}")

    def _extract_title_records(self, pdb_file):
        """Extract TITLE records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("TITLE"):
                        # Skip the "TITLE" part and any continuation number
                        title_text = line[10:].strip()
                        self.title_lines.append(title_text)
        except Exception as e:
            logger.warning(f"Error extracting TITLE records: {str(e)}")

    def _extract_author_records(self, pdb_file):
        """Extract AUTHOR records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("AUTHOR"):
                        # Skip the "AUTHOR" part and any continuation number
                        author_text = line[10:].strip()
                        self.author_lines.append(author_text)
        except Exception as e:
            logger.warning(f"Error extracting AUTHOR records: {str(e)}")

    def _extract_compnd_records(self, pdb_file):
        """Extract COMPND records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("COMPND"):
                        # Skip the "COMPND" part and any continuation number
                        compnd_text = line[10:].strip()
                        self.compnd_lines.append(compnd_text)
        except Exception as e:
            logger.warning(f"Error extracting COMPND records: {str(e)}")

    def _extract_source_records(self, pdb_file):
        """Extract SOURCE records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SOURCE"):
                        # Skip the "SOURCE" part and any continuation number
                        source_text = line[10:].strip()
                        self.source_lines.append(source_text)
        except Exception as e:
            logger.warning(f"Error extracting SOURCE records: {str(e)}")

    def _extract_keyword_records(self, pdb_file):
        """Extract KEYWDS records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                keyword_text = ""
                for line in f:
                    if line.startswith("KEYWDS"):
                        # Skip the "KEYWDS" part and any continuation number
                        keyword_text += line[10:].strip() + " "

                if keyword_text:
                    # Split by commas and strip whitespace
                    self.keywords = [k.strip() for k in keyword_text.split(",")]
        except Exception as e:
            logger.warning(f"Error extracting KEYWDS records: {str(e)}")

    def _extract_expdta_records(self, pdb_file):
        """Extract EXPDTA records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("EXPDTA"):
                        # Skip the "EXPDTA" part and any continuation number
                        expdta_text = line[10:].strip()
                        self.experimental_method = expdta_text
                        break
        except Exception as e:
            logger.warning(f"Error extracting EXPDTA record: {str(e)}")

    def _extract_remark_records(self, pdb_file):
        """Extract important REMARK records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                current_remark = None
                current_text = []

                for line in f:
                    if line.startswith("REMARK"):
                        try:
                            remark_num = int(line[7:10].strip())

                            # Check if this is a new remark number
                            if (
                                current_remark is not None
                                and current_remark != remark_num
                            ):
                                # Store the previous remark
                                self.remark_records[current_remark] = current_text
                                current_text = []

                            current_remark = remark_num

                            # Extract text part of the remark
                            remark_text = line[11:].strip()
                            if remark_text:  # Only add non-empty lines
                                current_text.append(remark_text)

                            # Special handling for resolution information in REMARK 2
                            if remark_num == 2 and "RESOLUTION" in line:
                                try:
                                    resolution_parts = line.split()
                                    for i, part in enumerate(resolution_parts):
                                        if part == "RESOLUTION.":
                                            self.resolution = float(
                                                resolution_parts[i + 1]
                                            )
                                            break
                                except (ValueError, IndexError) as e:
                                    pass

                        except ValueError:
                            # Not a numeric remark, just continue
                            pass

                # Store the last remark
                if current_remark is not None:
                    self.remark_records[current_remark] = current_text
        except Exception as e:
            logger.warning(f"Error extracting REMARK records: {str(e)}")

    def _extract_revdat_records(self, pdb_file):
        """Extract REVDAT records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("REVDAT"):
                        try:
                            mod_num = int(line[7:10].strip())
                            mod_date = line[13:22].strip()
                            mod_id = line[23:28].strip()
                            mod_type = line[31:35].strip()
                            mod_details = line[40:].strip()

                            self.revisions.append(
                                {
                                    "mod_num": mod_num,
                                    "date": mod_date,
                                    "id": mod_id,
                                    "type": mod_type,
                                    "details": mod_details,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing REVDAT line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting REVDAT records: {str(e)}")

    def _extract_ssbond_records(self, pdb_file):
        """Extract SSBOND records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SSBOND"):
                        try:
                            serial = int(line[7:10].strip())
                            res1_name = line[11:14].strip()
                            chain1 = line[15:16].strip()
                            # Handle insertion codes in residue numbers
                            res1_num_str = line[17:21].strip()
                            res1_icode = line[21:22].strip() if len(line) > 21 else ''
                            try:
                                res1_num = int(res1_num_str)
                                if res1_icode:
                                    res1_num = f"{res1_num}{res1_icode}"
                            except ValueError:
                                res1_num = res1_num_str

                            res2_name = line[25:28].strip()
                            chain2 = line[29:30].strip()
                            res2_num_str = line[31:35].strip()
                            res2_icode = line[35:36].strip() if len(line) > 35 else ''
                            try:
                                res2_num = int(res2_num_str)
                                if res2_icode:
                                    res2_num = f"{res2_num}{res2_icode}"
                            except ValueError:
                                res2_num = res2_num_str

                            # Extract distance if available
                            distance = None
                            if len(line) >= 73:
                                try:
                                    distance = float(line[73:78].strip())
                                except ValueError:
                                    pass

                            self.ssbond_records.append(
                                {
                                    "serial": serial,
                                    "res1_name": res1_name,
                                    "chain1": chain1,
                                    "res1_num": res1_num,
                                    "res2_name": res2_name,
                                    "chain2": chain2,
                                    "res2_num": res2_num,
                                    "distance": distance,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing SSBOND line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting SSBOND records: {str(e)}")

    def _extract_link_records(self, pdb_file):
        """Extract LINK records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("LINK"):
                        # Extract information from LINK record
                        # First atom
                        try:
                            name1 = line[12:16].strip()
                            resname1 = line[17:20].strip()
                            chain1 = line[21:22].strip()
                            try:
                                resid1 = int(line[22:26].strip())
                            except ValueError:
                                resid1 = line[22:26].strip()

                            # Second atom
                            name2 = line[42:46].strip()
                            resname2 = line[47:50].strip()
                            chain2 = line[51:52].strip()
                            try:
                                resid2 = int(line[52:56].strip())
                            except ValueError:
                                resid2 = line[52:56].strip()

                            # Distance might be included in some PDB files in later columns
                            distance = None
                            if len(line) >= 73:
                                try:
                                    distance = float(line[73:78].strip())
                                except (ValueError, IndexError):
                                    pass

                            self.link_records.append(
                                {
                                    "name1": name1,
                                    "resname1": resname1,
                                    "chain1": chain1,
                                    "resid1": resid1,
                                    "name2": name2,
                                    "resname2": resname2,
                                    "chain2": chain2,
                                    "resid2": resid2,
                                    "distance": distance,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing LINK line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting LINK records: {str(e)}")

    def _extract_seqres_records(self, pdb_file):
        """Extract SEQRES records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SEQRES"):
                        # Format of SEQRES:
                        # SEQRES  serial chainID numRes residue names (up to 13)
                        try:
                            chain_id = line[11:12].strip()
                            residues = line[19:].strip().split()

                            if chain_id not in self.seqres_records:
                                self.seqres_records[chain_id] = []

                            self.seqres_records[chain_id].extend(residues)
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing SEQRES line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting SEQRES records: {str(e)}")

    def _extract_modres_records(self, pdb_file):
        """Extract MODRES records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("MODRES"):
                        # Format of MODRES:
                        # MODRES recordName resName chainID resSeq iCode stdRes
                        try:
                            chain_id = line[16:17].strip()
                            # Handle insertion codes
                            res_num_str = line[18:22].strip()
                            res_icode = line[22:23].strip() if len(line) > 22 else ''
                            try:
                                res_num = int(res_num_str)
                                if res_icode:
                                    res_num = f"{res_num}{res_icode}"
                            except ValueError:
                                res_num = res_num_str
                            res_name = line[12:15].strip()
                            std_name = line[24:27].strip()

                            if chain_id not in self.modres_records:
                                self.modres_records[chain_id] = []

                            self.modres_records[chain_id].append(
                                {
                                    "residue_number": res_num,
                                    "residue_name": res_name,
                                    "standard_name": std_name,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing MODRES line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting MODRES records: {str(e)}")

    def _extract_missing_residue_records(self, pdb_file):
        """Extract REMARK 465 (missing residues) records from a PDB file."""
        reading_missing = False
        header_skipped = False

        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("REMARK 465"):
                        if "MISSING RESIDUES" in line:
                            reading_missing = True
                            continue

                        if reading_missing:
                            if not header_skipped and (
                                "RES C SSSEQI" in line or "M RES C SSSEQI" in line
                            ):
                                header_skipped = True
                                continue

                            # Skip lines until we get past the header
                            if not header_skipped:
                                continue

                            # Check if we have a valid line with residue info
                            if len(line) > 20:
                                try:
                                    # Format may vary between PDB files, but common format is:
                                    # REMARK 465   residue_name chain_id residue_number insertion_code
                                    chain_id = line[19:20].strip()
                                    if not chain_id:  # Skip blank lines
                                        continue

                                    res_name = line[15:18].strip()
                                    # Handle insertion codes
                                    res_num_str = line[21:26].strip()
                                    try:
                                        res_num = int(res_num_str)
                                    except ValueError:
                                        res_num = res_num_str

                                    if chain_id not in self.missing_res_records:
                                        self.missing_res_records[chain_id] = []

                                    self.missing_res_records[chain_id].append(
                                        {
                                            "residue_number": res_num,
                                            "residue_name": res_name,
                                        }
                                    )
                                except (ValueError, IndexError) as e:
                                    logger.debug(
                                        f"Error parsing REMARK 465 line: {line.strip()}, {str(e)}"
                                    )
                                    continue
                    elif reading_missing:
                        # A non-REMARK 465 line signals the end of the missing residues section
                        reading_missing = False
                        break
        except Exception as e:
            logger.warning(f"Error extracting missing residue records: {str(e)}")

    def _extract_helix_sheet_records(self, pdb_file):
        """Extract HELIX and SHEET records from a PDB file."""
        helix_types = {
            1: "Right-handed alpha",
            2: "Right-handed omega",
            3: "Right-handed pi",
            4: "Right-handed gamma",
            5: "Right-handed 3-10",
            6: "Left-handed alpha",
            7: "Left-handed omega",
            8: "Left-handed gamma",
            9: "2-7 ribbon/helix",
            10: "Polyproline",
        }

        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    # Parse HELIX records
                    if line.startswith("HELIX"):
                        try:
                            serial = int(line[7:10].strip())
                            helix_id = line[11:14].strip()
                            init_res_name = line[15:18].strip()
                            init_chain = line[19:20].strip()
                            # Handle insertion codes
                            init_seq_str = line[21:25].strip()
                            init_icode = line[25:26].strip() if len(line) > 25 else ''
                            try:
                                init_seq_num = int(init_seq_str)
                                if init_icode:
                                    init_seq_num = f"{init_seq_num}{init_icode}"
                            except ValueError:
                                init_seq_num = init_seq_str

                            term_res_name = line[27:30].strip()
                            term_chain = line[31:32].strip()
                            term_seq_str = line[33:37].strip()
                            term_icode = line[37:38].strip() if len(line) > 37 else ''
                            try:
                                term_seq_num = int(term_seq_str)
                                if term_icode:
                                    term_seq_num = f"{term_seq_num}{term_icode}"
                            except ValueError:
                                term_seq_num = term_seq_str
                            try:
                                helix_class = int(line[38:40].strip())
                            except ValueError:
                                helix_class = 1  # Default to alpha helix

                            length = term_seq_num - init_seq_num + 1

                            self.helix_records.append(
                                {
                                    "serial": serial,
                                    "id": helix_id,
                                    "start_res_name": init_res_name,
                                    "start_res": init_seq_num,
                                    "start_chain": init_chain,
                                    "end_res_name": term_res_name,
                                    "end_res": term_seq_num,
                                    "end_chain": term_chain,
                                    "class": helix_class,
                                    "type": helix_types.get(helix_class, "Unknown"),
                                    "length": length,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing HELIX record: {line.strip()}, {str(e)}"
                            )

                    # Parse SHEET records
                    elif line.startswith("SHEET"):
                        try:
                            strand = line[7:10].strip()
                            sheet_id = line[11:14].strip()
                            try:
                                num_strands = int(line[14:16].strip())
                            except ValueError:
                                num_strands = 1

                            init_res_name = line[17:20].strip()
                            init_chain = line[21:22].strip()
                            # Handle insertion codes
                            init_seq_str = line[22:26].strip()
                            init_icode = line[26:27].strip() if len(line) > 26 else ''
                            try:
                                init_seq_num = int(init_seq_str)
                                if init_icode:
                                    init_seq_num = f"{init_seq_num}{init_icode}"
                            except ValueError:
                                init_seq_num = init_seq_str

                            term_res_name = line[28:31].strip()
                            term_chain = line[32:33].strip()
                            term_seq_str = line[33:37].strip()
                            term_icode = line[37:38].strip() if len(line) > 37 else ''
                            try:
                                term_seq_num = int(term_seq_str)
                                if term_icode:
                                    term_seq_num = f"{term_seq_num}{term_icode}"
                            except ValueError:
                                term_seq_num = term_seq_str

                            # Sense of strand with respect to previous strand in the sheet
                            sense = 0
                            if len(line) >= 40:
                                try:
                                    sense = int(line[38:40].strip())
                                except ValueError:
                                    pass

                            length = term_seq_num - init_seq_num + 1

                            self.sheet_records.append(
                                {
                                    "strand": strand,
                                    "sheet_id": sheet_id,
                                    "num_strands": num_strands,
                                    "start_res_name": init_res_name,
                                    "start_chain": init_chain,
                                    "start_res": init_seq_num,
                                    "end_res_name": term_res_name,
                                    "end_chain": term_chain,
                                    "end_res": term_seq_num,
                                    "sense": sense,
                                    "length": length,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing SHEET record: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting HELIX/SHEET records: {str(e)}")

    def get_formatted_title(self):
        """Get the complete title formatted as a single string."""
        return " ".join(self.title_lines)

    def get_formatted_author(self):
        """Get the complete author list formatted as a single string."""
        return " ".join(self.author_lines)

    def get_formatted_compound(self):
        """Get the complete compound information formatted as a single string."""
        return " ".join(self.compnd_lines)

    def get_formatted_source(self):
        """Get the complete source information formatted as a single string."""
        return " ".join(self.source_lines)

    def get_structure_summary(self) -> Dict[str, Any]:
        """Get a summary of the structure for integration with other modules."""
        summary = {
            "pdb_id": self.header_info.get("pdb_id"),
            "title": self.get_formatted_title(),
            "experimental_method": self.experimental_method,
            "resolution": self.resolution,
            "chains": list(self.seqres_records.keys()),
            "has_ligands": len(self.het_records) > 0,
            "has_missing_residues": len(self.missing_res_records) > 0,
            "helix_count": len(self.helix_records),
            "sheet_count": len(self.sheet_records),
            "ssbond_count": len(self.ssbond_records),
        }
        return summary

    def get_chain_sequences(self) -> Dict[str, str]:
        """Get sequences for each chain."""
        sequences = {}
        for chain_id, residues in self.seqres_records.items():
            sequences[chain_id] = "".join(
                self._convert_three_to_one(res) for res in residues
            )
        return sequences

    def _convert_three_to_one(self, three_letter_code: str) -> str:
        """Convert three letter amino acid code to one letter code."""
        conversion = {
            "ALA": "A",
            "ARG": "R",
            "ASN": "N",
            "ASP": "D",
            "CYS": "C",
            "GLN": "Q",
            "GLU": "E",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LEU": "L",
            "LYS": "K",
            "MET": "M",
            "PHE": "F",
            "PRO": "P",
            "SER": "S",
            "THR": "T",
            "TRP": "W",
            "TYR": "Y",
            "VAL": "V",
            "MSE": "M",  # Selenomethionine as methionine
        }
        return conversion.get(three_letter_code, "X")

    def export_metadata_json(self, output_file: str = None) -> str:
        """
        Export metadata to a JSON file

        :param output_file: Optional output file path
        :return: Path to the saved file or JSON string if no file specified
        """
        import json

        # Create a serializable version of the metadata
        metadata = {
            "pdb_id": self.header_info.get("pdb_id"),
            "classification": self.header_info.get("classification"),
            "deposition_date": self.header_info.get("deposition_date"),
            "title": self.get_formatted_title(),
            "authors": self.get_formatted_author(),
            "experimental_method": self.experimental_method,
            "resolution": self.resolution,
            "compound": self.get_formatted_compound(),
            "source": self.get_formatted_source(),
            "keywords": self.keywords,
            "chains": list(self.seqres_records.keys()),
            "secondary_structure": {
                "helix_count": len(self.helix_records),
                "sheet_count": len(self.sheet_records),
                "ssbond_count": len(self.ssbond_records),
            },
            "missing_residues": {
                chain_id: len(residues)
                for chain_id, residues in self.missing_res_records.items()
            },
            "hetero_compounds": [het_id for het_id in self.hetnam_records.keys()],
        }

        # Add safety/status information
        if self.caveat_records or self.obslte_info or self.sprsde_info or self.split_info:
            metadata["safety_status"] = {}
            if self.caveat_records:
                metadata["safety_status"]["caveats"] = self.caveat_records
            if self.obslte_info:
                metadata["safety_status"]["obsolete"] = self.obslte_info
            if self.sprsde_info:
                metadata["safety_status"]["supersedes"] = self.sprsde_info
            if self.split_info:
                metadata["safety_status"]["split"] = self.split_info

        # Add database references
        if self.dbref_records:
            metadata["database_references"] = self.dbref_records
        if self.seqadv_records:
            metadata["sequence_differences"] = self.seqadv_records
        if self.cispep_records:
            metadata["cis_peptides"] = self.cispep_records

        # Add citations
        if self.jrnl_records:
            metadata["citations"] = self.jrnl_records

        # Add crystal information if available
        if self.cryst_info:
            metadata["crystal"] = {
                "space_group": self.cryst_info.get("space_group"),
                "cell_dimensions": {
                    "a": self.cryst_info.get("a"),
                    "b": self.cryst_info.get("b"),
                    "c": self.cryst_info.get("c"),
                    "alpha": self.cryst_info.get("alpha"),
                    "beta": self.cryst_info.get("beta"),
                    "gamma": self.cryst_info.get("gamma"),
                },
            }

        # Export to file or return as string
        if output_file:
            with open(output_file, "w") as f:
                json.dump(metadata, f, indent=2)
            return output_file
        else:
            return json.dumps(metadata, indent=2)

    def get_structure_validation(self) -> Dict[str, Any]:
        """
        Extract structure validation information from REMARK records

        :return: Dictionary with validation metrics
        """
        validation = {
            "resolution": self.resolution,
            "r_factor": None,
            "r_free": None,
            "clashscore": None,
            "ramachandran_outliers": None,
            "rotamer_outliers": None,
        }

        # Extract R-factor and R-free from REMARK 3
        if 3 in self.remark_records:
            for line in self.remark_records[3]:
                if "R VALUE" in line and "FREE" not in line:
                    parts = line.split(":")
                    if len(parts) > 1:
                        try:
                            validation["r_factor"] = float(parts[1].strip().split()[0])
                        except (ValueError, IndexError):
                            pass

                if "FREE R VALUE" in line:
                    parts = line.split(":")
                    if len(parts) > 1:
                        try:
                            validation["r_free"] = float(parts[1].strip().split()[0])
                        except (ValueError, IndexError):
                            pass

        # Extract validation metrics from other REMARK records if available
        if 500 in self.remark_records:
            for line in self.remark_records[500]:
                if "CLASHSCORE" in line:
                    parts = line.split(":")
                    if len(parts) > 1:
                        try:
                            validation["clashscore"] = float(
                                parts[1].strip().split()[0]
                            )
                        except (ValueError, IndexError):
                            pass

                if "RAMACHANDRAN OUTLIERS" in line:
                    parts = line.split(":")
                    if len(parts) > 1:
                        try:
                            validation["ramachandran_outliers"] = float(
                                parts[1].strip().split()[0]
                            )
                        except (ValueError, IndexError):
                            pass

                if "ROTAMER OUTLIERS" in line:
                    parts = line.split(":")
                    if len(parts) > 1:
                        try:
                            validation["rotamer_outliers"] = float(
                                parts[1].strip().split()[0]
                            )
                        except (ValueError, IndexError):
                            pass

        return validation

    def display_header_records(self):
        """Display detailed information about header records."""
        print("\n===== PDB HEADER RECORDS =====")

        # Display HEADER
        if self.header_info:
            print("\n----- HEADER -----")
            for key, value in self.header_info.items():
                print(f"{key}: {value}")

        # Display TITLE
        if self.title_lines:
            print("\n----- TITLE -----")
            for i, line in enumerate(self.title_lines):
                print(f"Line {i+1}: {line}")

        # Display AUTHOR
        if self.author_lines:
            print("\n----- AUTHOR -----")
            for i, line in enumerate(self.author_lines):
                print(f"Line {i+1}: {line}")

        # Display COMPND
        if self.compnd_lines:
            print("\n----- COMPOUND -----")
            for i, line in enumerate(self.compnd_lines):
                print(f"Line {i+1}: {line}")

        # Display SOURCE
        if self.source_lines:
            print("\n----- SOURCE -----")
            for i, line in enumerate(self.source_lines):
                print(f"Line {i+1}: {line}")

        # Display KEYWDS
        if self.keywords:
            print("\n----- KEYWORDS -----")
            print(", ".join(self.keywords))

        # Display EXPDTA
        if self.experimental_method:
            print("\n----- EXPERIMENTAL DATA -----")
            print(self.experimental_method)

        # Display REVDAT
        if self.revisions:
            print("\n----- REVISION HISTORY -----")
            for rev in self.revisions:
                print(
                    f"Modification {rev['mod_num']}: {rev['date']} - Type: {rev['type']} - Details: {rev['details']}"
                )

    def display_secondary_structure(self):
        """Display detailed information about secondary structure records."""
        print("\n===== SECONDARY STRUCTURE =====")

        # Display HELIX records
        if self.helix_records:
            print("\n----- HELICES -----")
            for helix in self.helix_records:
                print(f"Helix {helix['serial']} (ID: {helix['id']}):")
                print(f"  Type: {helix['type']} (Class: {helix['class']})")
                print(
                    f"  Location: Chain {helix['start_chain']} {helix['start_res_name']}{helix['start_res']} to "
                    + f"Chain {helix['end_chain']} {helix['end_res_name']}{helix['end_res']}"
                )
                print(f"  Length: {helix['length']} residues")

        # Display SHEET records
        if self.sheet_records:
            # Group sheets by sheet_id
            sheets = {}
            for strand in self.sheet_records:
                sheet_id = strand["sheet_id"]
                if sheet_id not in sheets:
                    sheets[sheet_id] = []
                sheets[sheet_id].append(strand)

            print("\n----- SHEETS -----")
            for sheet_id, strands in sheets.items():
                print(f"Sheet {sheet_id} ({len(strands)} strands):")
                for strand in strands:
                    sense_str = ""
                    if strand["sense"] == 0:
                        sense_str = "first strand"
                    elif strand["sense"] == 1:
                        sense_str = "parallel"
                    elif strand["sense"] == -1:
                        sense_str = "anti-parallel"

                    print(f"  Strand {strand['strand']} ({sense_str}):")
                    print(
                        f"    Location: Chain {strand['start_chain']} {strand['start_res_name']}{strand['start_res']} to "
                        + f"Chain {strand['end_chain']} {strand['end_res_name']}{strand['end_res']}"
                    )
                    print(f"    Length: {strand['length']} residues")

        # Display SSBOND records
        if self.ssbond_records:
            print("\n----- DISULFIDE BONDS -----")
            for bond in self.ssbond_records:
                distance_str = (
                    f", Distance: {bond['distance']:.2f} Å" if bond["distance"] else ""
                )
                print(
                    f"Bond {bond['serial']}: Chain {bond['chain1']} {bond['res1_name']}{bond['res1_num']} to "
                    + f"Chain {bond['chain2']} {bond['res2_name']}{bond['res2_num']}{distance_str}"
                )

        # Display LINK records
        if self.link_records:
            print("\n----- OTHER LINKS -----")
            for i, link in enumerate(self.link_records):
                distance_str = (
                    f", Distance: {link['distance']:.2f} Å" if link["distance"] else ""
                )
                print(
                    f"Link {i+1}: {link['name1']} in {link['resname1']}{link['resid1']} (Chain {link['chain1']}) to "
                    + f"{link['name2']} in {link['resname2']}{link['resid2']} (Chain {link['chain2']}){distance_str}"
                )

    def display_sequence_information(self):
        """Display sequence information from SEQRES, MODRES, and missing residues records."""
        print("\n===== SEQUENCE INFORMATION =====")

        # Display SEQRES records
        if self.seqres_records:
            print("\n----- SEQUENCE -----")
            for chain_id, residues in self.seqres_records.items():
                print(f"Chain {chain_id} ({len(residues)} residues):")

                # Format the sequence in groups of 10 for readability
                seq_lines = []
                for i in range(0, len(residues), 10):
                    chunk = residues[i : i + 10]
                    seq_lines.append(" ".join(chunk))

                for i, line in enumerate(seq_lines):
                    print(f"  {i*10+1:4d}-{min((i+1)*10, len(residues)):4d}: {line}")

        # Display MODRES records
        if self.modres_records:
            print("\n----- MODIFIED RESIDUES -----")
            for chain_id, residues in self.modres_records.items():
                print(f"Chain {chain_id} ({len(residues)} modified residues):")
                for res in residues:
                    print(
                        f"  Position {res['residue_number']}: {res['residue_name']} (Standard: {res['standard_name']})"
                    )

        # Display missing residues
        if self.missing_res_records:
            print("\n----- MISSING RESIDUES -----")
            for chain_id, residues in self.missing_res_records.items():
                print(f"Chain {chain_id} ({len(residues)} missing residues):")

                # Group consecutive missing residues for cleaner output
                groups = []
                current_group = []

                sorted_residues = sorted(residues, key=lambda x: x["residue_number"])

                for i, res in enumerate(sorted_residues):
                    if (
                        i == 0
                        or res["residue_number"]
                        != sorted_residues[i - 1]["residue_number"] + 1
                    ):
                        if current_group:
                            groups.append(current_group)
                        current_group = [res]
                    else:
                        current_group.append(res)

                if current_group:
                    groups.append(current_group)

                for group in groups:
                    if len(group) == 1:
                        res = group[0]
                        print(
                            f"  Position {res['residue_number']}: {res['residue_name']}"
                        )
                    else:
                        first = group[0]
                        last = group[-1]
                        print(
                            f"  Positions {first['residue_number']}-{last['residue_number']}: "
                            + f"{first['residue_name']}-{last['residue_name']} ({len(group)} residues)"
                        )

    def display_remarks(self):
        """Display important REMARK records."""
        print("\n===== REMARK RECORDS =====")

        # Specifically display important REMARKs
        important_remarks = [2, 3, 4, 200, 350, 465, 470, 500]

        for remark_num in important_remarks:
            if remark_num in self.remark_records:
                print(f"\n----- REMARK {remark_num} -----")
                for line in self.remark_records[remark_num]:
                    print(f"  {line}")

        # List other available remarks
        other_remarks = [
            num for num in self.remark_records.keys() if num not in important_remarks
        ]
        if other_remarks:
            print("\n----- OTHER REMARKS -----")
            print(
                f"Other available REMARK records: {', '.join(map(str, sorted(other_remarks)))}"
            )

            # Option to display a specific remark
            print(
                "\nTo display a specific REMARK, run the program with the --remarks flag and use:"
            )
            print("  python pdb_metadata.py <pdb_file> --display-remark <remark_num>")

    def _extract_het_records(self, pdb_file):
        """Extract HET records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("HET "):
                        # Format of HET:
                        # HET    het-id chain-id seq-num num-het-atoms
                        try:
                            het_id = line[7:10].strip()
                            chain_id = line[12:13].strip()
                            seq_num = line[13:17].strip()

                            # Combine chain and seq_num for unique identifier
                            het_key = f"{chain_id}_{seq_num}"

                            # Number of atoms in the hetero group
                            try:
                                num_atoms = int(line[20:25].strip())
                            except ValueError:
                                num_atoms = 0

                            # Description text if available
                            description = line[30:].strip()

                            self.het_records[het_key] = {
                                "het_id": het_id,
                                "chain_id": chain_id,
                                "seq_num": seq_num,
                                "num_atoms": num_atoms,
                                "description": description,
                            }
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing HET line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting HET records: {str(e)}")

    def _extract_hetnam_records(self, pdb_file):
        """Extract HETNAM records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("HETNAM"):
                        # Format of HETNAM:
                        # HETNAM het-id het-name
                        try:
                            het_id = line[11:14].strip()
                            het_name = line[15:].strip()

                            if het_id in self.hetnam_records:
                                # Append to existing name (continuation record)
                                self.hetnam_records[het_id] += " " + het_name
                            else:
                                self.hetnam_records[het_id] = het_name
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing HETNAM line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting HETNAM records: {str(e)}")

    def _extract_formul_records(self, pdb_file):
        """Extract FORMUL records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("FORMUL"):
                        # Format of FORMUL:
                        # FORMUL  comp-num  het-id  chemical-formula
                        try:
                            comp_num = line[8:10].strip()
                            het_id = line[12:15].strip()
                            formula = line[19:].strip()

                            if het_id in self.formul_records:
                                # Append to existing formula (continuation record)
                                self.formul_records[het_id] += " " + formula
                            else:
                                self.formul_records[het_id] = formula
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing FORMUL line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting FORMUL records: {str(e)}")

    def _extract_site_records(self, pdb_file):
        """Extract SITE records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                current_site = None

                for line in f:
                    if line.startswith("SITE"):
                        # Format of SITE:
                        # SITE   seq-num site-id num-residues
                        # followed by residue information
                        try:
                            seq_num = line[7:10].strip()
                            site_id = line[11:14].strip()

                            try:
                                num_residues = int(line[15:17].strip())
                            except ValueError:
                                num_residues = 0

                            # Create new site entry if it doesn't exist
                            if site_id not in self.site_records:
                                self.site_records[site_id] = {
                                    "seq_num": seq_num,
                                    "num_residues": num_residues,
                                    "residues": [],
                                }

                            current_site = site_id

                            # Extract residue information from this line
                            # SITE can have up to 4 residues per line
                            for i in range(4):
                                start_idx = 18 + (i * 11)
                                if start_idx + 10 <= len(line):
                                    res_name = line[start_idx : start_idx + 3].strip()
                                    chain_id = line[
                                        start_idx + 4 : start_idx + 5
                                    ].strip()

                                    try:
                                        res_num = int(
                                            line[start_idx + 5 : start_idx + 9].strip()
                                        )
                                    except ValueError:
                                        continue

                                    if res_name and chain_id:
                                        self.site_records[current_site][
                                            "residues"
                                        ].append(
                                            {
                                                "res_name": res_name,
                                                "chain_id": chain_id,
                                                "res_num": res_num,
                                            }
                                        )
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing SITE line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting SITE records: {str(e)}")

    def _extract_crystal_records(self, pdb_file):
        """Extract CRYST1 record from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("CRYST1"):
                        # Format of CRYST1:
                        # CRYST1   a     b     c     alpha  beta  gamma  space-group  z
                        try:
                            a = float(line[6:15].strip())
                            b = float(line[15:24].strip())
                            c = float(line[24:33].strip())
                            alpha = float(line[33:40].strip())
                            beta = float(line[40:47].strip())
                            gamma = float(line[47:54].strip())
                            space_group = line[55:66].strip()
                            z_value = line[66:70].strip()

                            self.cryst_info = {
                                "a": a,
                                "b": b,
                                "c": c,
                                "alpha": alpha,
                                "beta": beta,
                                "gamma": gamma,
                                "space_group": space_group,
                                "z": z_value,
                            }

                            # We only need one CRYST1 record
                            break
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing CRYST1 line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting CRYST1 record: {str(e)}")

    def _extract_origx_records(self, pdb_file):
        """Extract ORIGXn records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("ORIGX"):
                        # Format of ORIGXn:
                        # ORIGXn   o[n][1]  o[n][2]  o[n][3]  t[n]
                        try:
                            n = int(line[5:6])
                            o1 = float(line[10:20].strip())
                            o2 = float(line[20:30].strip())
                            o3 = float(line[30:40].strip())
                            t = float(line[45:55].strip())

                            # Make sure we have a list long enough
                            while len(self.origx_matrices) < n:
                                self.origx_matrices.append(None)

                            self.origx_matrices[n - 1] = {
                                "o1": o1,
                                "o2": o2,
                                "o3": o3,
                                "t": t,
                            }
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing ORIGX line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting ORIGX records: {str(e)}")

    def _extract_scale_records(self, pdb_file):
        """Extract SCALEn records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SCALE"):
                        # Format of SCALEn:
                        # SCALEn   s[n][1]  s[n][2]  s[n][3]  u[n]
                        try:
                            n = int(line[5:6])
                            s1 = float(line[10:20].strip())
                            s2 = float(line[20:30].strip())
                            s3 = float(line[30:40].strip())
                            u = float(line[45:55].strip())

                            # Make sure we have a list long enough
                            while len(self.scale_matrices) < n:
                                self.scale_matrices.append(None)

                            self.scale_matrices[n - 1] = {
                                "s1": s1,
                                "s2": s2,
                                "s3": s3,
                                "u": u,
                            }
                        except (ValueError, IndexError) as e:
                            logger.debug(
                                f"Error parsing SCALE line: {line.strip()}, {str(e)}"
                            )
        except Exception as e:
            logger.warning(f"Error extracting SCALE records: {str(e)}")

    def _extract_caveat_records(self, pdb_file):
        """Extract CAVEAT records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("CAVEAT"):
                        # CAVEAT records warn about severe errors in the entry
                        continuation = line[8:10].strip()
                        pdb_id = line[11:15].strip()
                        comment = line[19:79].strip()
                        self.caveat_records.append({
                            "continuation": continuation,
                            "pdb_id": pdb_id,
                            "comment": comment,
                        })
        except Exception as e:
            logger.warning(f"Error extracting CAVEAT records: {str(e)}")

    def _extract_obslte_records(self, pdb_file):
        """Extract OBSLTE record from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("OBSLTE"):
                        # OBSLTE indicates this entry is obsolete
                        continuation = line[8:10].strip()
                        replacement_date = line[11:20].strip()
                        old_pdb_id = line[21:25].strip()
                        # Replacement PDB IDs start at position 31 and repeat every 10 chars
                        replacement_ids = []
                        for i in range(31, 71, 10):
                            if i + 4 <= len(line):
                                rep_id = line[i:i+4].strip()
                                if rep_id:
                                    replacement_ids.append(rep_id)

                        self.obslte_info = {
                            "continuation": continuation,
                            "date": replacement_date,
                            "old_id": old_pdb_id,
                            "replacement_ids": replacement_ids,
                        }
                        break  # Only one OBSLTE record per file
        except Exception as e:
            logger.warning(f"Error extracting OBSLTE record: {str(e)}")

    def _extract_sprsde_records(self, pdb_file):
        """Extract SPRSDE record from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SPRSDE"):
                        # SPRSDE indicates this entry supersedes others
                        continuation = line[8:10].strip()
                        supersede_date = line[11:20].strip()
                        this_pdb_id = line[21:25].strip()
                        # Superseded PDB IDs start at position 31
                        superseded_ids = []
                        for i in range(31, 71, 10):
                            if i + 4 <= len(line):
                                old_id = line[i:i+4].strip()
                                if old_id:
                                    superseded_ids.append(old_id)

                        self.sprsde_info = {
                            "continuation": continuation,
                            "date": supersede_date,
                            "this_id": this_pdb_id,
                            "superseded_ids": superseded_ids,
                        }
                        break
        except Exception as e:
            logger.warning(f"Error extracting SPRSDE record: {str(e)}")

    def _extract_split_records(self, pdb_file):
        """Extract SPLIT record from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SPLIT "):
                        # SPLIT indicates this entry was split into multiple entries
                        continuation = line[8:10].strip()
                        # Split PDB IDs start at position 11
                        split_ids = []
                        for i in range(11, 71, 10):
                            if i + 4 <= len(line):
                                split_id = line[i:i+4].strip()
                                if split_id:
                                    split_ids.append(split_id)

                        self.split_info = {
                            "continuation": continuation,
                            "split_ids": split_ids,
                        }
                        break
        except Exception as e:
            logger.warning(f"Error extracting SPLIT record: {str(e)}")

    def _extract_dbref_records(self, pdb_file):
        """Extract DBREF records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("DBREF "):
                        # DBREF provides cross-references to sequence databases
                        pdb_id = line[7:11].strip()
                        chain_id = line[12:13].strip()
                        seq_begin = line[14:18].strip()
                        insert_begin = line[18:19].strip()
                        seq_end = line[20:24].strip()
                        insert_end = line[24:25].strip()
                        database = line[26:32].strip()
                        db_accession = line[33:41].strip()
                        db_id_code = line[42:54].strip()
                        db_seq_begin = line[55:60].strip()
                        db_insert_begin = line[60:61].strip()
                        db_seq_end = line[62:67].strip()
                        db_insert_end = line[67:68].strip()

                        self.dbref_records.append({
                            "pdb_id": pdb_id,
                            "chain_id": chain_id,
                            "seq_begin": seq_begin,
                            "insert_begin": insert_begin,
                            "seq_end": seq_end,
                            "insert_end": insert_end,
                            "database": database,
                            "db_accession": db_accession,
                            "db_id_code": db_id_code,
                            "db_seq_begin": db_seq_begin,
                            "db_insert_begin": db_insert_begin,
                            "db_seq_end": db_seq_end,
                            "db_insert_end": db_insert_end,
                        })
        except Exception as e:
            logger.warning(f"Error extracting DBREF records: {str(e)}")

    def _extract_seqadv_records(self, pdb_file):
        """Extract SEQADV records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("SEQADV"):
                        # SEQADV indicates differences between PDB and database sequences
                        pdb_id = line[7:11].strip()
                        res_name = line[12:15].strip()
                        chain_id = line[16:17].strip()
                        seq_num = line[18:22].strip()
                        insert_code = line[22:23].strip()
                        database = line[24:28].strip()
                        db_accession = line[29:38].strip()
                        db_res = line[39:42].strip()
                        db_seq = line[43:48].strip()
                        conflict = line[49:70].strip()

                        self.seqadv_records.append({
                            "pdb_id": pdb_id,
                            "res_name": res_name,
                            "chain_id": chain_id,
                            "seq_num": seq_num,
                            "insert_code": insert_code,
                            "database": database,
                            "db_accession": db_accession,
                            "db_res": db_res,
                            "db_seq": db_seq,
                            "conflict": conflict,
                        })
        except Exception as e:
            logger.warning(f"Error extracting SEQADV records: {str(e)}")

    def _extract_cispep_records(self, pdb_file):
        """Extract CISPEP records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("CISPEP"):
                        # CISPEP identifies cis peptide bonds
                        ser_num = line[7:10].strip()
                        pep1 = line[11:14].strip()
                        chain1 = line[15:16].strip()
                        seq1 = line[17:21].strip()
                        icode1 = line[21:22].strip()
                        pep2 = line[25:28].strip()
                        chain2 = line[29:30].strip()
                        seq2 = line[31:35].strip()
                        icode2 = line[35:36].strip()
                        modnum = line[43:46].strip()
                        measure = line[53:59].strip()

                        self.cispep_records.append({
                            "ser_num": ser_num,
                            "pep1": pep1,
                            "chain1": chain1,
                            "seq1": seq1,
                            "icode1": icode1,
                            "pep2": pep2,
                            "chain2": chain2,
                            "seq2": seq2,
                            "icode2": icode2,
                            "modnum": modnum,
                            "measure": measure,
                        })
        except Exception as e:
            logger.warning(f"Error extracting CISPEP records: {str(e)}")

    def _extract_jrnl_records(self, pdb_file):
        """Extract JRNL records from a PDB file."""
        try:
            with open(pdb_file, "r") as f:
                current_jrnl = {}
                for line in f:
                    if line.startswith("JRNL"):
                        # JRNL records contain citation information
                        subrecord = line[12:16].strip()
                        content = line[19:79].strip()

                        if subrecord == "AUTH":
                            if "authors" not in current_jrnl:
                                current_jrnl["authors"] = []
                            current_jrnl["authors"].append(content)
                        elif subrecord == "TITL":
                            if "title" not in current_jrnl:
                                current_jrnl["title"] = []
                            current_jrnl["title"].append(content)
                        elif subrecord == "EDIT":
                            if "editors" not in current_jrnl:
                                current_jrnl["editors"] = []
                            current_jrnl["editors"].append(content)
                        elif subrecord == "REF":
                            if "reference" not in current_jrnl:
                                current_jrnl["reference"] = []
                            current_jrnl["reference"].append(content)
                        elif subrecord == "PUBL":
                            if "publisher" not in current_jrnl:
                                current_jrnl["publisher"] = []
                            current_jrnl["publisher"].append(content)
                        elif subrecord == "REFN":
                            current_jrnl["refn"] = content
                        elif subrecord == "PMID":
                            current_jrnl["pmid"] = content
                        elif subrecord == "DOI":
                            current_jrnl["doi"] = content
                    elif current_jrnl:
                        # End of JRNL section
                        self.jrnl_records.append(current_jrnl)
                        current_jrnl = {}

                # Don't forget the last one if file ends with JRNL
                if current_jrnl:
                    self.jrnl_records.append(current_jrnl)

        except Exception as e:
            logger.warning(f"Error extracting JRNL records: {str(e)}")

    def _detect_hplusplus_structure(self):
        """
        Detect if this PDB file was generated by the H++ server.

        H++ structures have:
        - REMARK containing "Created by http://biophysics.cs.vt.edu/H++"
        - Hydrogens already added
        - AMBER forcefield residue naming conventions (HID, HIE, HIP, CYX, etc.)
        """
        try:
            # Check all remark records for H++ signature
            for remark_num, remark_lines in self.remark_records.items():
                for line in remark_lines:
                    if "biophysics.cs.vt.edu/H++" in line or "H++ server" in line.lower():
                        self.is_hplusplus_structure = True
                        return
        except Exception as e:
            logger.warning(f"Error detecting H++ structure: {str(e)}")
            # Default to False if detection fails
            self.is_hplusplus_structure = False

    def display_hetero_compounds(self):
        """Display information about hetero compounds."""
        print("\n===== HETERO COMPOUNDS =====")

        # Check if we have any hetero compound information
        if not self.het_records and not self.hetnam_records and not self.formul_records:
            print("No hetero compound information found.")
            return

        # Create a comprehensive view by combining information from different records
        het_ids = set()
        for record_dict in [self.het_records, self.hetnam_records, self.formul_records]:
            for key in record_dict:
                if isinstance(key, str) and len(key) == 3:  # HET ID is 3 characters
                    het_ids.add(key)

        # Also add HET IDs from het_records which are stored with chain/seq info
        for het_key in self.het_records:
            if "_" in het_key:  # This is a chain_seqnum key
                het_id = self.het_records[het_key]["het_id"]
                het_ids.add(het_id)

        # Display information for each hetero compound
        for het_id in sorted(het_ids):
            print(f"\n----- {het_id} -----")

            # Display name if available
            if het_id in self.hetnam_records:
                print(f"Name: {self.hetnam_records[het_id]}")

            # Display formula if available
            if het_id in self.formul_records:
                print(f"Formula: {self.formul_records[het_id]}")

            # Display occurrences of this hetero compound
            occurrences = []
            for het_key, het_info in self.het_records.items():
                if het_info["het_id"] == het_id:
                    occurrences.append(het_info)

            if occurrences:
                print(f"Occurrences: {len(occurrences)}")
                for occ in occurrences:
                    atoms_str = (
                        f", {occ['num_atoms']} atoms" if occ["num_atoms"] > 0 else ""
                    )
                    desc_str = f" - {occ['description']}" if occ["description"] else ""
                    print(
                        f"  Chain {occ['chain_id']} {occ['seq_num']}{atoms_str}{desc_str}"
                    )

    def display_crystallographic_info(self):
        """Display crystallographic information."""
        print("\n===== CRYSTALLOGRAPHIC INFORMATION =====")

        # Display CRYST1 information
        if self.cryst_info:
            print("\n----- UNIT CELL PARAMETERS -----")
            print(f"a = {self.cryst_info['a']:.3f} Å")
            print(f"b = {self.cryst_info['b']:.3f} Å")
            print(f"c = {self.cryst_info['c']:.3f} Å")
            print(f"α = {self.cryst_info['alpha']:.2f}°")
            print(f"β = {self.cryst_info['beta']:.2f}°")
            print(f"γ = {self.cryst_info['gamma']:.2f}°")
            print(f"Space Group: {self.cryst_info['space_group']}")
            if self.cryst_info["z"]:
                print(f"Z Value: {self.cryst_info['z']}")
        else:
            print("No unit cell parameters found.")

        # Display ORIGX matrices
        if self.origx_matrices:
            print("\n----- ORIGX MATRICES -----")
            print("Original coordinate transformation:")
            for i, matrix in enumerate(self.origx_matrices):
                if matrix:
                    print(
                        f"  ORIGX{i+1}: {matrix['o1']:10.6f} {matrix['o2']:10.6f} {matrix['o3']:10.6f} {matrix['t']:10.6f}"
                    )

        # Display SCALE matrices
        if self.scale_matrices:
            print("\n----- SCALE MATRICES -----")
            print("Fractional to orthogonal transformation:")
            for i, matrix in enumerate(self.scale_matrices):
                if matrix:
                    print(
                        f"  SCALE{i+1}: {matrix['s1']:10.6f} {matrix['s2']:10.6f} {matrix['s3']:10.6f} {matrix['u']:10.6f}"
                    )

    def display_site_records(self):
        """Display SITE records."""
        print("\n===== SITE INFORMATION =====")

        if not self.site_records:
            print("No site information found.")
            return

        for site_id, site_info in self.site_records.items():
            print(f"\n----- SITE {site_id} -----")
            print(f"Number of residues: {site_info['num_residues']}")

            if site_info["residues"]:
                print("Residues:")
                for res in site_info["residues"]:
                    print(f"  {res['res_name']} {res['chain_id']} {res['res_num']}")

    def display_safety_status_records(self):
        """Display safety and status records (CAVEAT, OBSLTE, SPRSDE, SPLIT)."""
        print("\n===== SAFETY & STATUS RECORDS =====")

        has_content = False

        if self.caveat_records:
            has_content = True
            print("\n** CAVEAT - Severe errors/warnings **")
            for caveat in self.caveat_records:
                print(f"  {caveat['comment']}")

        if self.obslte_info:
            has_content = True
            print("\n** OBSOLETE ENTRY **")
            print(f"  Date: {self.obslte_info['date']}")
            print(f"  Superseded by: {', '.join(self.obslte_info['replacement_ids'])}")

        if self.sprsde_info:
            has_content = True
            print("\n** SUPERSEDES **")
            print(f"  Date: {self.sprsde_info['date']}")
            print(f"  Supersedes: {', '.join(self.sprsde_info['superseded_ids'])}")

        if self.split_info:
            has_content = True
            print("\n** SPLIT ENTRY **")
            print(f"  Split into: {', '.join(self.split_info['split_ids'])}")

        if not has_content:
            print("No safety or status warnings found.")

    def display_database_references(self):
        """Display database cross-references and sequence information."""
        print("\n===== DATABASE REFERENCES =====")

        has_content = False

        if self.dbref_records:
            has_content = True
            print("\n** DBREF - Database Cross-References **")
            for dbref in self.dbref_records:
                print(f"  Chain {dbref['chain_id']}: {dbref['database']} {dbref['db_accession']} ({dbref['db_id_code']})")
                print(f"    PDB range: {dbref['seq_begin']}-{dbref['seq_end']}")
                print(f"    DB range: {dbref['db_seq_begin']}-{dbref['db_seq_end']}")

        if self.seqadv_records:
            has_content = True
            print("\n** SEQADV - Sequence Differences **")
            for seqadv in self.seqadv_records:
                print(f"  Chain {seqadv['chain_id']}, {seqadv['res_name']} {seqadv['seq_num']}: {seqadv['conflict']}")
                print(f"    Database: {seqadv['database']} {seqadv['db_accession']}")

        if self.cispep_records:
            has_content = True
            print("\n** CISPEP - Cis Peptide Bonds **")
            for cispep in self.cispep_records:
                print(f"  {cispep['pep1']} {cispep['chain1']}:{cispep['seq1']} - {cispep['pep2']} {cispep['chain2']}:{cispep['seq2']}")
                if cispep['measure']:
                    print(f"    Angle: {cispep['measure']}°")

        if not has_content:
            print("No database reference information found.")

    def display_citations(self):
        """Display citation and journal information."""
        print("\n===== CITATIONS =====")

        if not self.jrnl_records:
            print("No citation information found.")
            return

        for i, jrnl in enumerate(self.jrnl_records, 1):
            if len(self.jrnl_records) > 1:
                print(f"\n** Citation {i} **")

            if "authors" in jrnl:
                authors = " ".join(jrnl["authors"])
                print(f"  Authors: {authors}")

            if "title" in jrnl:
                title = " ".join(jrnl["title"])
                print(f"  Title: {title}")

            if "reference" in jrnl:
                reference = " ".join(jrnl["reference"])
                print(f"  Reference: {reference}")

            if "publisher" in jrnl:
                publisher = " ".join(jrnl["publisher"])
                print(f"  Publisher: {publisher}")

            if "pmid" in jrnl:
                print(f"  PubMed ID: {jrnl['pmid']}")

            if "doi" in jrnl:
                print(f"  DOI: {jrnl['doi']}")

    def display_all(self):
        """Display all metadata information."""
        self.display_metadata_summary()
        self.display_safety_status_records()
        self.display_header_records()
        self.display_secondary_structure()
        self.display_sequence_information()
        self.display_database_references()
        self.display_citations()
        self.display_hetero_compounds()
        self.display_crystallographic_info()
        self.display_site_records()
        self.display_remarks()

    def display_metadata_summary(self, pdb_id=None):
        """Display a summary of the extracted metadata.

        Args:
            pdb_id: Optional PDB ID fallback (used when header contains XXXX)
        """
        print("\n===== PDB METADATA SUMMARY =====")

        # Display basic info
        if self.header_info:
            extracted_id = self.header_info.get('pdb_id', '')
            display_id = extracted_id if extracted_id and extracted_id != 'XXXX' else (pdb_id or 'Unknown')
            print(f"\nID: {display_id}")
            print(
                f"Classification: {self.header_info.get('classification', 'Unknown')}"
            )
            print(
                f"Deposition Date: {self.header_info.get('deposition_date', 'Unknown')}"
            )

        # Display title
        title = self.get_formatted_title()
        if title:
            print(f"\nTitle: {title}")

        # Display authors
        authors = self.get_formatted_author()
        if authors:
            print(f"\nAuthors: {authors}")

        # Display experimental method
        if self.experimental_method:
            print(f"\nExperimental Method: {self.experimental_method}")

        # Display resolution if available
        if self.resolution:
            print(f"Resolution: {self.resolution} Å")

        # Display compound info
        compound = self.get_formatted_compound()
        if compound:
            print(f"\nCompound: {compound}")

        # Display source info
        source = self.get_formatted_source()
        if source:
            print(f"\nSource: {source}")

        # Display keywords
        if self.keywords:
            print(f"\nKeywords: {', '.join(self.keywords)}")

        # Display crystallographic information if available
        if self.cryst_info:
            print(f"\nCrystallographic Information:")
            print(f"  Space Group: {self.cryst_info['space_group']}")
            print(
                f"  Cell Dimensions: a={self.cryst_info['a']:.1f}, b={self.cryst_info['b']:.1f}, c={self.cryst_info['c']:.1f}"
            )
            print(
                f"  Cell Angles: α={self.cryst_info['alpha']:.1f}°, β={self.cryst_info['beta']:.1f}°, γ={self.cryst_info['gamma']:.1f}°"
            )

        # Display revision history
        if self.revisions:
            print("\nRevision History:")
            for rev in self.revisions:
                print(
                    f"  {rev['mod_num']}: {rev['date']} - {rev['type']} {rev['details']}"
                )

        # Display hetero compound summary
        het_count = (
            len(set([info["het_id"] for info in self.het_records.values()]))
            if self.het_records
            else 0
        )
        if het_count > 0:
            print(f"\nHetero Compounds: {het_count} unique types")

        # Display site information
        if self.site_records:
            print(f"\nSites: {len(self.site_records)} defined")

        # Display secondary structure summary
        print("\nSecondary Structure Summary:")
        print(f"  Helices: {len(self.helix_records)}")
        print(f"  Sheets: {len(self.sheet_records)}")
        print(f"  Disulfide Bonds: {len(self.ssbond_records)}")
        print(f"  Other Links: {len(self.link_records)}")

        # Display chain info
        print("\nChain Information:")
        for chain_id, residues in self.seqres_records.items():
            print(f"  Chain {chain_id}: {len(residues)} residues")

            # Check for modified residues
            mod_res = self.modres_records.get(chain_id, [])
            if mod_res:
                print(f"    Modified Residues: {len(mod_res)}")

            # Check for missing residues
            missing_res = self.missing_res_records.get(chain_id, [])
            if missing_res:
                print(f"    Missing Residues: {len(missing_res)}")

        # Summary of remarks
        remark_counts = len(self.remark_records)
        if remark_counts > 0:
            print(f"\nRemarks: {remark_counts} types")
            if 2 in self.remark_records:
                print("  Including resolution information (REMARK 2)")
            if 3 in self.remark_records:
                print("  Including refinement information (REMARK 3)")
            if 350 in self.remark_records:
                print("  Including biological assembly information (REMARK 350)")
            if 465 in self.remark_records:
                print("  Including missing residue information (REMARK 465)")


@register_module
# Updated PDBLoaderModule class for pdb_loader.py


@register_module
class PDBLoaderModule(ProcessingModule):
    """Module for loading and analyzing PDB files"""

    NAME = "PDB Loader"
    DESCRIPTION = "Load PDB structures and extract metadata"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    PRIORITY = 1

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module Display Hook

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for PDB Loader module data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        from rich.console import Console
        from rich.panel import Panel
        from rich.table import Table
        from rich.text import Text

        # Handle PDB metadata display
        if workspace_key == "metadata" and hasattr(value, "header_info"):
            # Create a separate console for rendering tables
            string_console = Console(width=75, highlight=False)

            # Basic information section
            basic_info = []
            if hasattr(value, "header_info") and value.header_info:
                basic_info.append(
                    f"[cyan]PDB ID:[/cyan] {value.header_info.get('pdb_id', 'Unknown')}"
                )
                basic_info.append(
                    f"[cyan]Classification:[/cyan] {value.header_info.get('classification', 'Unknown')}"
                )
                basic_info.append(
                    f"[cyan]Deposition Date:[/cyan] {value.header_info.get('deposition_date', 'Unknown')}"
                )

            # Title
            title = (
                value.get_formatted_title()
                if hasattr(value, "get_formatted_title")
                else "Unknown"
            )
            basic_info.append(f"[cyan]Title:[/cyan] {title}")

            # Experimental information
            if hasattr(value, "experimental_method") and value.experimental_method:
                basic_info.append(
                    f"[cyan]Experimental Method:[/cyan] {value.experimental_method}"
                )

            if hasattr(value, "resolution") and value.resolution:
                basic_info.append(f"[cyan]Resolution:[/cyan] {value.resolution} Å")

            # Authors
            authors = (
                value.get_formatted_author()
                if hasattr(value, "get_formatted_author")
                else None
            )
            if authors:
                basic_info.append(f"[cyan]Authors:[/cyan] {authors}")

            # Structure validation
            if hasattr(value, "get_structure_validation"):
                validation = value.get_structure_validation()
                if validation.get("r_factor") is not None:
                    basic_info.append(
                        f"[cyan]R-factor:[/cyan] {validation['r_factor']}"
                    )
                if validation.get("r_free") is not None:
                    basic_info.append(f"[cyan]R-free:[/cyan] {validation['r_free']}")

            # Create secondary structure section
            ss_section = []
            if (
                hasattr(value, "helix_records")
                or hasattr(value, "sheet_records")
                or hasattr(value, "ssbond_records")
            ):
                ss_section.append("[cyan]Secondary Structure:[/cyan]")
                ss_section.append(
                    f"  [cyan]Helices:[/cyan] {len(value.helix_records) if hasattr(value, 'helix_records') else 0}"
                )
                ss_section.append(
                    f"  [cyan]Sheets:[/cyan] {len(value.sheet_records) if hasattr(value, 'sheet_records') else 0}"
                )
                ss_section.append(
                    f"  [cyan]Disulfide Bonds:[/cyan] {len(value.ssbond_records) if hasattr(value, 'ssbond_records') else 0}"
                )
                ss_section.append(
                    f"  [cyan]Other Links:[/cyan] {len(value.link_records) if hasattr(value, 'link_records') else 0}"
                )

            # Create chain information section
            chain_section = []
            if hasattr(value, "seqres_records") and value.seqres_records:
                chain_section.append("[cyan]Chain Information:[/cyan]")
                for chain_id, residues in value.seqres_records.items():
                    mod_count = (
                        len(value.modres_records.get(chain_id, []))
                        if hasattr(value, "modres_records")
                        else 0
                    )
                    missing_count = (
                        len(value.missing_res_records.get(chain_id, []))
                        if hasattr(value, "missing_res_records")
                        else 0
                    )
                    chain_section.append(
                        f"  [cyan]Chain {chain_id}:[/cyan] {len(residues)} residues, {mod_count} modified, {missing_count} missing"
                    )

            # Hetero compounds section
            het_section = []
            if hasattr(value, "hetnam_records") and value.hetnam_records:
                het_section.append(
                    f"[cyan]Hetero Compounds:[/cyan] {len(value.hetnam_records)} unique types"
                )

                for het_id, het_name in list(value.hetnam_records.items())[
                    :5
                ]:  # Show first 5
                    het_section.append(f"  [cyan]{het_id}:[/cyan] {het_name}")

                if len(value.hetnam_records) > 5:
                    het_section.append(
                        f"  ... and {len(value.hetnam_records) - 5} more"
                    )

            # Crystallographic information section
            cryst_section = []
            if hasattr(value, "cryst_info") and value.cryst_info:
                cryst_section.append("[cyan]Crystallographic Information:[/cyan]")
                cryst_section.append(
                    f"  [cyan]Space Group:[/cyan] {value.cryst_info.get('space_group', 'Unknown')}"
                )
                cryst_section.append(
                    f"  [cyan]Cell Dimensions:[/cyan] a={value.cryst_info.get('a', 0):.1f}, b={value.cryst_info.get('b', 0):.1f}, c={value.cryst_info.get('c', 0):.1f}"
                )
                cryst_section.append(
                    f"  [cyan]Cell Angles:[/cyan] α={value.cryst_info.get('alpha', 0):.1f}°, β={value.cryst_info.get('beta', 0):.1f}°, γ={value.cryst_info.get('gamma', 0):.1f}°"
                )

            # Combine all sections with proper spacing
            all_sections = []
            if basic_info:
                all_sections.append("\n".join(basic_info))
            if ss_section:
                all_sections.append("\n".join(ss_section))
            if chain_section:
                all_sections.append("\n".join(chain_section))
            if het_section:
                all_sections.append("\n".join(het_section))
            if cryst_section:
                all_sections.append("\n".join(cryst_section))

            content = "\n\n".join(all_sections)

            # Create the main panel
            metadata_panel = Panel(
                content,
                title=f"PDB Metadata: {value.header_info.get('pdb_id', 'Unknown') if hasattr(value, 'header_info') else 'Unknown'}",
                border_style="green",
                expand=False,
            )

            console.print(metadata_panel)
            return True

        return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Helper methods to centralize interactions with the workspace

    def get_workspace(self):
        """
        Helper method to get the appropriate workspace object

        Returns:
            The current workspace object
        """
        return self.processor.workspace

    def get_from_workspace(self, key, default=None):
        """
        Helper method to get values from the processor's workspace

        Args:
            key: Key to retrieve
            default: Default value if key not found

        Returns:
            The value for the given key or default
        """
        return self.processor.workspace.get(key, default)

    def get_from_workspace_obj(self, workspace, key, default=None):
        """
        Helper method to get values from a workspace object

        Args:
            workspace: Workspace object
            key: Key to retrieve
            default: Default value if key not found

        Returns:
            The value for the given key or default
        """
        return workspace.get(key, default)

    def update_workspace(self, key, value):
        """
        Helper method to update the processor's workspace

        Args:
            key: Key to update
            value: New value
        """
        old_value = self.processor.workspace.get(key)
        self.processor.workspace.set(key, value)
        debug_enabled = self.processor.workspace.get("debug", False)

        # Debug output if enabled
        if debug_enabled and hasattr(self.processor, "console"):
            self._debug_value(key, value, old_value, "updated")

    def update_workspace_obj(self, workspace, key, value):
        """
        Helper method to update a workspace object

        Args:
            workspace: Workspace object
            key: Key to update
            value: New value

        Returns:
            Updated workspace object
        """
        debug_enabled = workspace.get("debug", False)
        old_value = workspace.get(key)

        workspace.set(key, value)

        # Debug output if enabled
        if debug_enabled and hasattr(self.processor, "console"):
            self._debug_value(key, value, old_value, "updated")

        return workspace

    # Migration helper methods for workspace key transition
    # TODO: Remove these methods and old key support after all modules are updated
    # Key mappings:
    # - pdb_file -> original_pdb_file
    # - metadata -> original_metadata  
    # - structure -> original_structure
    # - residue_numbering_info -> original_residue_numbering_info
    
    def _get_workspace_value(self, workspace, old_key, new_key, default=None):
        """
        Helper method to get values with migration support
        Checks new key first, falls back to old key with warning
        
        Args:
            workspace: Workspace object
            old_key: Old key name (deprecated)
            new_key: New key name (preferred)
            default: Default value if neither key exists
            
        Returns:
            Value from workspace
        """
        if workspace.get(new_key) is not None:
            return workspace.get(new_key)
        elif workspace.get(old_key) is not None:
            if hasattr(self.processor, "console") and workspace.get("debug", False):
                self.processor.console.print(f"[yellow]Warning: Using deprecated workspace key '{old_key}'. Please use '{new_key}' instead.[/yellow]")
            return workspace.get(old_key)
        else:
            return default
    
    def _set_workspace_value(self, workspace, old_key, new_key, value):
        """
        Helper method to set values with migration support
        Sets both old and new keys during transition period
        
        Args:
            workspace: Workspace object
            old_key: Old key name (for backward compatibility)
            new_key: New key name (preferred)
            value: Value to set
            
        Returns:
            Updated workspace object
        """
        # Set new key (preferred)
        workspace.set(new_key, value)
        # Set old key (for backward compatibility)
        workspace.set(old_key, value)
        
        debug_enabled = workspace.get("debug", False)
        if debug_enabled and hasattr(self.processor, "console"):
            self.processor.console.print(f"[grey50]Set workspace keys: '{new_key}' and '{old_key}' (deprecated)[/grey50]")
        
        return workspace

    def _debug_value(self, key, value, old_value, action):
        """
        Debug output for value updates

        Args:
            key: Key being updated
            value: New value
            old_value: Previous value
            action: Action being performed
        """
        if key == "debug":
            return  # Don't debug the debug flag itself

        value_type = type(value).__name__
        if hasattr(value, "__len__") and not isinstance(value, (str, bytes)):
            value_info = f"{value_type}({len(value)})"
        else:
            value_info = value_type

        if old_value is None:
            self.processor.console.print(
                f"[yellow]DEBUG: {self.NAME} {action} new value for '{key}': {value_info}[/yellow]"
            )
        else:
            old_type = type(old_value).__name__
            if hasattr(old_value, "__len__") and not isinstance(
                old_value, (str, bytes)
            ):
                old_info = f"{old_type}({len(old_value)})"

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Input from/Ouput to the workspace

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements"""
        return []  # PDB Loader can work without any pre-existing workspace keys

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "original_metadata", "original_structure",
            "metadata", "structure",  # deprecated aliases
        ]

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Module Menu

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "load": "Load a PDB file",
            "metadata": "View PDB metadata",
        }

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status and dependencies
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus, determine_option_status

        options = []

        # Option 1: Load a PDB file
        # Always available (no prerequisites)
        pdb_loaded = workspace.get("pdb_file") is not None or workspace.get("structure") is not None
        if pdb_loaded:
            load_status = OptionStatus.COMPLETED
        else:
            load_status = OptionStatus.AVAILABLE

        options.append(MenuOption(
            key="1",
            description="Load a PDB file",
            status=load_status
        ))

        # Option 2: View PDB metadata
        # Requires PDB to be loaded first
        if pdb_loaded:
            metadata_status = OptionStatus.READY
        else:
            metadata_status = OptionStatus.BLOCKED

        options.append(MenuOption(
            key="2",
            description="View PDB metadata",
            status=metadata_status,
            dependency_text="[Need to load first] ○" if not pdb_loaded else ""
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.

        Args:
            workspace: Current workspace

        Returns:
            Suggestion text or None
        """
        pdb_loaded = workspace.get("pdb_file") is not None or workspace.get("structure") is not None

        # Check if this is the first time seeing metadata as available
        shown_metadata_info = workspace.get("_metadata_info_shown", False)

        if not pdb_loaded:
            return "Start by loading a PDB file (option 1)"
        elif pdb_loaded and not shown_metadata_info:
            # Mark as shown for this session
            workspace.set("_metadata_info_shown", True)

            return ("PDB files contain detailed information beyond atomic coordinates, "
                    "including experimental methods (X-ray/NMR/Cryo-EM), structure quality metrics, sequences, "
                    "secondary structure, bound ligands, citations, and more. Use option 2 to explore.")
        else:
            return "View metadata with option 2, or press [m] to return to the main menu"

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        if option == "load":
            if self.processor:
                self.processor.load_pdb_interactive()
            return True
        elif option == "metadata":
            if self.processor:
                self.processor.run_metadata_module()
            return True
        return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Get Structure

    def get_structure(self, workspace):
        """
        Get the Biopython structure object from a loaded PDB file

        Args:
            workspace: Current workspace

        Returns:
            Structure: Biopython Structure object or None
        """
        pdb_file = self._get_workspace_value(workspace, "pdb_file", "original_pdb_file")
        if not pdb_file:
            return None

        try:
            from Bio.PDB import PDBParser

            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", pdb_file)
            return structure
        except Exception as e:
            if hasattr(self, "processor") and self.processor:
                self.processor.console.print(
                    f"[red]Error parsing PDB file: {str(e)}[/red]"
                )
            return None

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#
    # Process

    def can_process(self, workspace) -> bool:
        """
        Check if the module can process the current workspace

        Args:
            workspace: Current workspace (dict or Workspace)
        """
        # PDB Loader doesn't have any specific requirements
        # It can process any workspace
        return True

    def process(self, workspace) -> Any:
        """
        Process the workspace with the PDB Loader module

        Args:
            workspace: Current workspace (dict or Workspace object)

        Returns:
            Updated workspace
        """

        # Debug the incoming workspace
        if (
            hasattr(self, "processor")
            and self.processor
            and workspace.get("debug", False)
        ):
            debug_workspace("PDBLoader", "receiving", workspace, self.processor.console)
            self.logger.debug("PDB Loader process method started")

        pdb_file = self._get_workspace_value(workspace, "pdb_file", "original_pdb_file")
        if pdb_file:
            if hasattr(self, "processor") and self.processor:
                self.logger.debug(f"Loading PDB file: {pdb_file}")

            # Load PDB file and extract metadata
            metadata = PDBMetadataExtractor(pdb_file)

            # Update workspace with metadata (both old and new keys)
            workspace = self._set_workspace_value(workspace, "metadata", "original_metadata", metadata)

            if hasattr(self, "processor") and self.processor:
                self.logger.debug("Metadata extracted successfully")

            # Get structure using Biopython
            try:
                if hasattr(self, "processor") and self.processor:
                    self.logger.debug("Loading structure with Biopython PDBParser")

                from Bio.PDB import PDBParser

                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("structure", pdb_file)

                # Update workspace with structure (both old and new keys)
                workspace = self._set_workspace_value(workspace, "structure", "original_structure", structure)

                if hasattr(self, "processor") and self.processor:
                    self.logger.debug("Structure loaded successfully")

                # Check residue numbering
                try:
                    if hasattr(self, "processor") and self.processor:
                        self.processor.check_residue_numbering()

                except Exception as e:
                    if hasattr(self, "processor") and self.processor:
                        self.logger.error(f"Error in residue numbering check: {str(e)}")

            except Exception as e:
                if hasattr(self, "processor") and self.processor:
                    self.logger.error(f"Error loading structure: {str(e)}")

        # Debug the outgoing workspace
        if (
            hasattr(self, "processor")
            and self.processor
            and workspace.get("debug", False)
        ):
            debug_workspace("PDBLoader", "returning", workspace, self.processor.console)
            self.logger.debug("PDB Loader process method completed")

        return workspace


# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#+#+#+#


def main():
    """Main function to run the PDB metadata extractor as a standalone program."""
    parser = argparse.ArgumentParser(
        description="Extract and display metadata from PDB files"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--pdbid", help="PDB ID to download and analyze")
    group.add_argument("--pdbfile", help="Path to local PDB file to analyze")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument("--export", help="Export metadata to JSON file")
    parser.add_argument(
        "--display-remark", type=int, help="Display specific REMARK record"
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger(__name__).setLevel(logging.CRITICAL)

    try:
        # Interactive menu for loading PDB
        pdb_file = None

        if args.pdbid:
            # Download the specified PDB ID
            pdb_file = download_pdb(args.pdbid)
        elif args.pdbfile:
            # Use the specified local file
            if not os.path.exists(args.pdbfile):
                raise FileNotFoundError(f"PDB file not found: {args.pdbfile}")
            pdb_file = args.pdbfile
        else:
            # No arguments provided, display menu
            pdb_file = load_pdb_interactive()

        if not pdb_file:
            print("No PDB file specified. Exiting.")
            sys.exit(0)

        # Load the PDB file
        print(f"\nLoading PDB file: {pdb_file}")
        extractor = PDBMetadataExtractor(pdb_file)

        # Process specific options
        if args.export:
            extractor.export_metadata_json(args.export)
            print(f"Exported metadata to {args.export}")
        elif args.display_remark is not None:
            remark_num = args.display_remark
            if remark_num in extractor.remark_records:
                print(f"\n----- REMARK {remark_num} -----")
                for line in extractor.remark_records[remark_num]:
                    print(f"  {line}")
            else:
                print(f"REMARK {remark_num} not present in this PDB file.")
        else:
            # Interactive menu for viewing metadata
            display_metadata_menu(extractor)

    except FileNotFoundError as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


def download_pdb(pdb_id):
    """Download a PDB file from RCSB PDB database."""
    import urllib.request
    from urllib.error import HTTPError

    pdb_id = pdb_id.upper()
    output_file = f"{pdb_id}.pdb"

    if os.path.exists(output_file):
        print(f"Using existing PDB file: {output_file}")
        return output_file

    print(f"Downloading PDB file {pdb_id}...")
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        urllib.request.urlretrieve(pdb_url, output_file)
        print(f"Downloaded PDB file to {output_file}")
        return output_file
    except HTTPError as e:
        if e.code == 404:
            # Display informative message about large structure limitation
            print(f"\n{'='*80}")
            print(f"LARGE STRUCTURE LIMITATION")
            print(f"{'='*80}")
            print(f"PDB {pdb_id} is not available in PDB format (likely a large structure).")
            print(f"RCSB only provides mmCIF format for structures with >99,999 atoms.")
            print(f"")
            print(f"CURRENT LIMITATION:")
            print(f"  • ProPrep's metal site transformation pipeline requires PDB format")
            print(f"  • This pipeline uses direct line editing to preserve non-standard residues")
            print(f"  • BioPython drops many non-standard groups needed for forcefields")
            print(f"")
            print(f"WORKAROUNDS:")
            print(f"  1. Manual download: https://files.rcsb.org/download/{pdb_id}.cif")
            print(f"  2. Convert CIF→PDB using VMD, PyMOL, or other tools")
            print(f"  3. Use for topology analysis only (structure filtering works)")
            print(f"")
            print(f"TO FULLY SUPPORT CIF FORMAT:")
            print(f"  • Rewrite metal site transformation pipeline (~2000 lines)")
            print(f"  • Replace all direct PDB line editing with CIF parsing")
            print(f"  • Maintain dual format support for backward compatibility")
            print(f"  • Extensive testing of metal site coordinate transformations")
            print(f"")
            print(f"IMPACT: Metal site preparation will not work with CIF files")
            print(f"SAFE: Structure loading, filtering, topology analysis work fine")
            print(f"{'='*80}")
            return None
        else:
            print(f"Error downloading PDB file: {str(e)}")
            return None
    except Exception as e:
        print(f"Error downloading PDB file: {str(e)}")
        return None


def list_pdb_files(directory="."):
    """
    List all PDB files in the specified directory

    Args:
        directory: Directory path to search for PDB files (default: current directory)

    Returns:
        list: List of PDB files
    """
    import glob
    import os

    # Get all files with .pdb extension (case insensitive)
    pdb_files = []

    # Try both .pdb and .PDB extensions
    pdb_files.extend(glob.glob(os.path.join(directory, "*.pdb")))
    pdb_files.extend(glob.glob(os.path.join(directory, "*.PDB")))

    # Sort the files
    pdb_files.sort()

    return pdb_files


def display_pdb_file_menu(directory=".", console=None, processor=None, allow_multiple=False):
    """
    Display a menu of available PDB files with directory navigation.

    A thin wrapper over :func:`proprep.utils.file_browser.file_browser`; the
    actual UX (bare-``N`` selection, ``q`` to cancel, comma/range multi-select,
    filename-based session replay) lives there and is shared with every other
    ProPrep file browser.

    Args:
        directory: Directory path to search for PDB files (default: current directory)
        console: Rich console object for display (optional)
        allow_multiple: When True, the selection accepts several indices at once
            (e.g. ``1,3,5`` or ``1-5,7``) and the function returns a list of file
            paths. When False (default) only a single index is accepted and a
            single path string is returned, so existing single-file callers are
            unaffected.

    Returns:
        - allow_multiple=False (default): str or None — selected file path,
          or None if canceled.
        - allow_multiple=True: list[str] or None — list of selected file paths,
          or None if canceled.
    """
    # Thin wrapper over the shared interactive browser so the pdb picker shares
    # the unified UX (bare-N selection, q to cancel, comma/range multi-select)
    # and filename-based session replay with every other ProPrep file browser.
    from proprep.utils.file_browser import file_browser, default_size_detail

    return file_browser(
        directory=directory,
        extensions=[".pdb"],
        console=console,
        processor=processor,
        multi=allow_multiple,
        label="PDB file",
        entry_detail=default_size_detail,
        allow_path_jump=True,
        module="PDB Loader",
    )


def load_pdb_interactive(processor=None):
    """Interactive menu for loading a PDB file."""
    import os
    import shutil
    from pathlib import Path

    from rich.console import Console
    from proprep.utils.prompts import prompt_with_context

    console = Console()

    while True:
        console.print("\n===== LOAD PDB FILE =====")
        console.print("1. Download PDB by ID", highlight=False)
        console.print("2. Load local PDB file", highlight=False)
        console.print("3. Exit", highlight=False)

        try:
            if processor:
                choice = prompt_with_context(
                    processor=processor,
                    prompt="\nEnter your choice",
                    choices=["1", "2", "3"],
                    default="1",
                    module="PDB Loader",
                    description="Select PDB loading method",
                    options_map={
                        "1": "Download PDB by ID",
                        "2": "Load local PDB file",
                        "3": "Exit"
                    }
                )
            else:
                choice = prompt_with_context(None,
                    "\nEnter your choice", choices=["1", "2", "3"], default="1"
                )

            if choice == "1":
                if processor:
                    pdb_id = prompt_with_context(
                        processor=processor,
                        prompt="Enter PDB ID (e.g., 1HHO)",
                        module="PDB Loader",
                        description="Enter PDB ID to download"
                    ).strip()
                else:
                    pdb_id = prompt_with_context(None, "Enter PDB ID (e.g., 1HHO)").strip()
                if pdb_id:
                    # Import the setup function from main module
                    from proprep.main import setup_results_directory
                    
                    # Setup results directory for PDB ID
                    setup_results_directory(pdb_id, is_pdb_id=True)
                    
                    pdb_file = download_pdb(pdb_id)
                    if pdb_file:
                        return pdb_file
                else:
                    console.print("[yellow]Invalid PDB ID.[/yellow]")
            elif choice == "2":
                # Show the file selection menu
                file_path = display_pdb_file_menu(directory=".", console=console, processor=processor)
                if file_path:
                    # Import the setup function from main module
                    from proprep.main import setup_results_directory
                    
                    # Setup results directory for PDB file
                    original_dir = os.getcwd()
                    setup_results_directory(file_path, is_pdb_id=False)
                    
                    # Copy the PDB file to the results directory
                    pdb_filename = Path(file_path).name
                    full_path = Path(original_dir) / file_path
                    shutil.copy2(full_path, pdb_filename)
                    
                    return pdb_filename
                else:
                    # If user canceled or no files found, allow manual entry
                    if processor:
                        file_path = prompt_with_context(
                            processor=processor,
                            prompt="Enter path to PDB file",
                            module="PDB Loader",
                            description="Enter local PDB file path"
                        ).strip()
                    else:
                        file_path = prompt_with_context(None, "Enter path to PDB file").strip()
                    if file_path and os.path.exists(file_path):
                        # Import the setup function from main module
                        from proprep.main import setup_results_directory
                        
                        # Setup results directory for PDB file
                        setup_results_directory(file_path, is_pdb_id=False)
                        
                        # Copy the PDB file to the results directory
                        pdb_filename = Path(file_path).name
                        shutil.copy2(file_path, pdb_filename)
                        
                        return pdb_filename
                    else:
                        console.print("[yellow]File not found.[/yellow]")
            elif choice == "3":
                return None
            else:
                console.print("[yellow]Invalid choice. Please try again.[/yellow]")
        except KeyboardInterrupt:
            console.print("\n[yellow]Operation cancelled by user.[/yellow]")
            return None
        except Exception as e:
            console.print(f"[red]Error: {str(e)}[/red]")

def display_metadata_menu(extractor, console=None, processor=None, pdb_id=None):
    """
    Interactive menu for viewing metadata.

    :param extractor: PDBMetadataExtractor instance
    :param console: Optional Rich console instance
    :param processor: Optional processor instance for session recording
    :param pdb_id: Optional PDB ID fallback (used when header contains XXXX, e.g. biological assemblies)
    """
    # Use provided console or create a new one
    if console is None:
        from rich.console import Console
        console = Console()

    while True:
        console.print("\n===== PDB METADATA VIEWER =====")
        console.print(f"File: {extractor.pdb_file}")
        if extractor.header_info:
            extracted_id = extractor.header_info.get('pdb_id', '')
            display_id = extracted_id if extracted_id and extracted_id != 'XXXX' else (pdb_id or 'Unknown')
            console.print(f"ID: {display_id}")

        console.print("\nAvailable Sections:")

        # Calculate the width needed for right alignment (15 is the highest number)
        max_option_num = 15
        num_width = len(str(max_option_num))

        # Display options with right-aligned numbers and inline descriptions
        options = [
            ("1", "Summary", "Quick overview of structure and experimental details"),
            ("2", "Safety & Status Records", "Obsolete/superseded status and warnings"),
            ("3", "Header Records", "Classification, deposition date, and keywords"),
            ("4", "Secondary Structure", "Helices, sheets, and turns"),
            ("5", "Sequence Information", "Amino acid sequences and chain details"),
            ("6", "Database References", "Links to UniProt, GenBank, etc."),
            ("7", "Citations", "Publications and authors"),
            ("8", "Hetero Compounds", "Ligands, cofactors, and non-standard residues"),
            ("9", "Crystallographic Information", "Space group, unit cell, resolution"),
            ("10", "Site Information", "Active sites, binding sites, metal centers"),
            ("11", "Remark Records", "Quality metrics, R-factors, refinement details"),
            ("12", "View All", "Display all metadata sections at once"),
            ("13", "Export Metadata", "Save metadata to JSON or text file"),
            ("14", "Help", "Detailed explanations of what each section contains"),
            ("15", "Exit", "Return to previous menu")
        ]

        for option_num, option_text, description in options:
            formatted_num = f"{option_num}.".rjust(num_width + 1)
            console.print(f"{formatted_num} {option_text:<35} [grey50]- {description}[/grey50]")

        try:
            from proprep.utils.prompts import prompt_with_context

            if processor:
                choice = prompt_with_context(
                    processor=processor,
                    prompt="\nEnter your choice",
                    choices=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"],
                    default="1",
                    module="PDB Metadata Viewer",
                    description="Select metadata section to view",
                    options_map={
                        "1": "Summary",
                        "2": "Safety & Status Records",
                        "3": "Header Records",
                        "4": "Secondary Structure",
                        "5": "Sequence Information",
                        "6": "Database References",
                        "7": "Citations",
                        "8": "Hetero Compounds",
                        "9": "Crystallographic Information",
                        "10": "Site Information",
                        "11": "Remark Records",
                        "12": "View All",
                        "13": "Export Metadata",
                        "14": "Help",
                        "15": "Exit"
                    }
                )
            else:
                choice = prompt_with_context(None,
                    "\nEnter your choice",
                    choices=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"],
                    default="1",
                )

            if choice == "1":
                extractor.display_metadata_summary()
            elif choice == "2":
                extractor.display_safety_status_records()
            elif choice == "3":
                extractor.display_header_records()
            elif choice == "4":
                extractor.display_secondary_structure()
            elif choice == "5":
                extractor.display_sequence_information()
            elif choice == "6":
                extractor.display_database_references()
            elif choice == "7":
                extractor.display_citations()
            elif choice == "8":
                extractor.display_hetero_compounds()
            elif choice == "9":
                extractor.display_crystallographic_info()
            elif choice == "10":
                extractor.display_site_records()
            elif choice == "11":
                display_remarks_menu(extractor, console, processor)
            elif choice == "12":
                extractor.display_all()
            elif choice == "13":
                _export_metadata_menu(extractor, console, processor)
            elif choice == "14":
                _display_metadata_help(console)
            elif choice == "15":
                console.print("Exiting metadata viewer.")
                return
            else:
                console.print("Invalid choice. Please try again.")

        except KeyboardInterrupt:
            console.print("\nOperation cancelled by user.")
            return
        except Exception as e:
            console.print(f"Error: {str(e)}")


def _display_metadata_help(console):
    """
    Display detailed help about metadata sections.

    :param console: Rich console instance
    """
    from rich.panel import Panel
    from rich.table import Table

    help_table = Table(show_header=True, header_style="bold cyan", show_lines=True)
    help_table.add_column("Section", style="cyan", width=30)
    help_table.add_column("What It Contains", width=50)
    help_table.add_column("When To Use", width=40)

    help_table.add_row(
        "1. Summary",
        "Quick overview including structure ID, experimental method, resolution, organism, and chain information",
        "First stop for understanding the structure basics"
    )

    help_table.add_row(
        "2. Safety & Status",
        "Obsolete/superseded status, warnings about structure quality or data issues",
        "Check if this is the most current version before using the structure"
    )

    help_table.add_row(
        "3. Header Records",
        "Classification (enzyme, DNA, membrane protein, etc.), deposition/release dates, keywords, compound names",
        "Understanding the biological context and finding related structures"
    )

    help_table.add_row(
        "4. Secondary Structure",
        "Location and details of helices, beta sheets, and turns as assigned by the author or DSSP",
        "Analyzing protein fold, comparing structural elements, or extracting specific motifs"
    )

    help_table.add_row(
        "5. Sequence Information",
        "Complete amino acid sequences for all chains, modified residues, sequence conflicts with databases",
        "Preparing for mutation studies, checking for non-standard residues, or sequence alignment"
    )

    help_table.add_row(
        "6. Database References",
        "Cross-references to UniProt, GenBank, PubMed, and other databases with accession codes",
        "Tracing the structure back to sequence databases or finding related data"
    )

    help_table.add_row(
        "7. Citations",
        "Primary publication, authors, journal, DOI, and other related publications",
        "Citing the structure properly or finding the methods paper"
    )

    help_table.add_row(
        "8. Hetero Compounds",
        "Ligands, cofactors, metal ions, modified residues, and their chemical formulas",
        "Critical for MD setup - identifies non-standard components needing special parameters"
    )

    help_table.add_row(
        "9. Crystallographic Info",
        "Space group, unit cell dimensions, crystal symmetry, resolution, R-factors",
        "Assessing structure quality, understanding crystal packing effects"
    )

    help_table.add_row(
        "10. Site Information",
        "Active sites, binding sites, metal coordination sites annotated by authors",
        "Identifying functionally important regions, especially for enzyme studies"
    )

    help_table.add_row(
        "11. Remark Records",
        "Quality metrics (R-free, R-work, Ramachandran statistics), refinement details, model accuracy",
        "Evaluating structure reliability before using for modeling or docking"
    )

    help_table.add_row(
        "12. View All",
        "Displays all metadata sections in sequence",
        "Comprehensive review or when searching for specific information"
    )

    help_table.add_row(
        "13. Export Metadata",
        "Save metadata to JSON or text file for external analysis",
        "Batch processing, archiving, or programmatic access to metadata"
    )

    console.print("\n")
    console.print(Panel(help_table, title="[bold]PDB Metadata Guide[/bold]", border_style="cyan"))


def _export_metadata_menu(extractor, console, processor=None):
    """
    Menu for exporting metadata

    :param extractor: PDBMetadataExtractor instance
    :param console: Rich console instance
    """
    from proprep.utils.prompts import prompt_with_context

    console.print("\n===== EXPORT METADATA =====")
    console.print("1. Export as JSON", highlight=False)
    console.print("2. Export as plain text summary", highlight=False)
    console.print("3. Return to metadata menu", highlight=False)

    if processor:
        choice = prompt_with_context(
            processor=processor,
            prompt="\nEnter your choice",
            choices=["1", "2", "3"],
            default="1",
            module="PDB Metadata Viewer",
            description="Select export format",
            options_map={
                "1": "Export as JSON",
                "2": "Export as plain text summary",
                "3": "Return to metadata menu"
            }
        )
    else:
        choice = prompt_with_context(None, "\nEnter your choice", choices=["1", "2", "3"], default="1")

    if choice == "3":
        return

    # Get output filename
    default_filename = f"{extractor.header_info.get('pdb_id', 'pdb')}_metadata"
    default_filename += ".json" if choice == "1" else ".txt"

    if processor:
        output_file = prompt_with_context(
            processor=processor,
            prompt="Enter output filename",
            default=default_filename,
            module="PDB Metadata Viewer",
            description="Enter filename for exported metadata"
        )
    else:
        output_file = prompt_with_context(None, "Enter output filename", default=default_filename)

    try:
        if choice == "1":
            # Export as JSON
            extractor.export_metadata_json(output_file)
            console.print(f"[green]Exported metadata to {output_file}[/green]")
        elif choice == "2":
            # Export as plain text summary
            with open(output_file, "w") as f:
                # Redirect stdout temporarily
                import sys

                old_stdout = sys.stdout
                sys.stdout = f

                # Write summary to file
                extractor.display_metadata_summary()

                # Restore stdout
                sys.stdout = old_stdout

            console.print(f"[green]Exported summary to {output_file}[/green]")

    except Exception as e:
        console.print(f"[red]Error exporting metadata: {str(e)}[/red]")


def display_remarks_menu(extractor, console=None, processor=None):
    """
    Menu for viewing specific remark records

    :param extractor: PDBMetadataExtractor instance
    :param console: Optional Rich console instance
    :param processor: Optional processor instance for session recording
    """
    # Use provided console or create a new one
    if console is None:
        from rich.console import Console

        console = Console()

    from proprep.utils.prompts import prompt_with_context

    # Define REMARK descriptions
    remark_descriptions = {
        0: "Re-refinement notice",
        1: "Related publications",
        2: "Resolution (Å)",
        3: "Refinement statistics",
        4: "Version and compliance",
        100: "Processing site details",
        200: "X-ray diffraction experiment",
        205: "Fiber diffraction details",
        210: "NMR experiment details",
        215: "Solution NMR details",
        217: "Solid state NMR details",
        230: "Neutron diffraction details",
        240: "Electron crystallography",
        245: "Electron microscopy (EM)",
        247: "EM coordinate generation",
        250: "Other experimental techniques",
        265: "Solution scattering",
        280: "Crystal properties (Matthews, solvent)",
        285: "CRYST1 unit cell details",
        290: "Crystallographic symmetry operators",
        300: "Biomolecule description",
        350: "Biomolecule generation",
        375: "Special position atoms",
        400: "Macromolecular compound details",
        450: "Biological source details",
        465: "Missing residues",
        470: "Missing atoms",
        475: "Zero occupancy residues",
        480: "Zero occupancy atoms",
        500: "Geometry and stereochemistry issues",
        525: "Solvent atoms distance info",
        600: "Heterogen details",
        610: "Missing atoms in heterogens",
        615: "Zero occupancy in heterogens",
        620: "Metal coordination",
        630: "Inhibitor/compound description",
        650: "Helix details",
        700: "Sheet details",
        800: "Important functional sites",
        900: "Related PDB entries",
        999: "Unusual sequence information",
    }

    while True:
        console.print("\n===== REMARK RECORDS =====")

        # List all available remarks with descriptions
        available_remarks = sorted(extractor.remark_records.keys())

        if not available_remarks:
            console.print("\n[yellow]No REMARK records found in this PDB file.[/yellow]")
            console.print("\nOptions:")
            console.print("0. Back to main menu", highlight=False)

            try:
                if processor:
                    choice = prompt_with_context(
                        processor=processor,
                        prompt="\nEnter your choice",
                        choices=["0"],
                        default="0",
                        module="PDB Metadata Viewer",
                        description="Return to main menu (no REMARKs found)"
                    )
                else:
                    choice = prompt_with_context(None, "\nEnter your choice", choices=["0"], default="0")
                if choice == "0":
                    return
            except KeyboardInterrupt:
                console.print("\nOperation cancelled by user.")
                return
            continue

        console.print("\nAvailable Remarks in this file:")
        console.print()

        # Display in a nice table format
        from rich.table import Table
        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("REMARK", style="green", width=8)
        table.add_column("Description", style="white")

        for remark_num in available_remarks:
            description = remark_descriptions.get(remark_num, "")
            if description:
                table.add_row(str(remark_num), description)
            else:
                table.add_row(str(remark_num), "[grey50]No description available[/grey50]")

        console.print(table)

        console.print("\nEnter a REMARK number to view, or 0 to go back to main menu")

        try:
            if processor:
                choice = prompt_with_context(
                    processor=processor,
                    prompt="\nEnter REMARK number (or 0 to exit)",
                    default="0",
                    module="PDB Metadata Viewer",
                    description="Enter REMARK number to view or 0 to exit"
                )
            else:
                choice = prompt_with_context(None,
                    "\nEnter REMARK number (or 0 to exit)",
                    default="0",
                )

            if choice == "0":
                return

            try:
                remark_num = int(choice)
                if remark_num in extractor.remark_records:
                    description = remark_descriptions.get(remark_num, "No description available")
                    console.print(f"\n----- REMARK {remark_num}: {description} -----")
                    for line in extractor.remark_records[remark_num]:
                        console.print(f"  {line}")
                else:
                    console.print(
                        f"[yellow]REMARK {remark_num} not present in this PDB file.[/yellow]"
                    )
            except ValueError:
                console.print("[red]Invalid input. Please enter a valid REMARK number.[/red]")

            # Wait for user to press enter before continuing
            if processor:
                prompt_with_context(
                    processor=processor,
                    prompt="\nPress Enter to continue",
                    default="",
                    module="PDB Metadata Viewer",
                    description="Press Enter to return to REMARKs menu"
                )
            else:
                prompt_with_context(None, "\nPress Enter to continue", default="")

        except KeyboardInterrupt:
            console.print("\nOperation cancelled by user.")
            return
        except Exception as e:
            console.print(f"Error: {str(e)}")
            if processor:
                prompt_with_context(
                    processor=processor,
                    prompt="\nPress Enter to continue",
                    default="",
                    module="PDB Metadata Viewer",
                    description="Press Enter after error"
                )
            else:
                prompt_with_context(None, "\nPress Enter to continue", default="")


if __name__ == "__main__":
    main()
