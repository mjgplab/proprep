"""
BLAST Worker

Core functionality for performing NCBI BLAST searches with detailed parameter control.
"""

import logging
import math
import os
import re
import sys
import tempfile
import time
import urllib.request
from typing import Any, Dict, List, Optional

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger(__name__)


class BLASTWorker:
    """Handles BLAST searches with detailed parameter control."""

    def __init__(self):
        """Initialize the BLAST parameters."""
        self.program = "blastn"
        self.database = "nt"
        self.sequence = ""
        self.e_value = 10.0
        self.word_size = None
        self.max_hits = 50
        self.megablast = True
        self.output_format = "XML"
        self.results = None
        self.additional_params = {}

    def setup(self, program=None, database=None, sequence=None, **kwargs):
        """
        Set up the BLAST parameters.

        Args:
            program: BLAST program type
            database: Target database
            sequence: Query sequence
            **kwargs: Additional parameters
        """
        if program:
            self.program = program
        if database:
            self.database = database
        if sequence:
            self.sequence = sequence

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.additional_params[key] = value

    def download_fasta(self, pdb_id: str) -> Dict[str, str]:
        """Download FASTA sequence for a PDB ID."""
        pdb_id = pdb_id.lower()
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/download"

        try:
            with urllib.request.urlopen(fasta_url) as response:
                fasta_content = response.read().decode("utf-8")

                sequences = {}
                current_chains = []
                current_seq = ""

                for line in fasta_content.splitlines():
                    if line.startswith(">"):
                        if current_chains and current_seq:
                            for chain in current_chains:
                                sequences[chain] = current_seq

                        header_match = re.search(r"\|Chains\s+([^|]+)\|", line)
                        if header_match:
                            chains_part = header_match.group(1)
                            chain_matches = re.findall(
                                r"([A-Za-z0-9])(?:\[auth [A-Za-z0-9]\])?", chains_part
                            )
                            current_chains = chain_matches if chain_matches else []
                        else:
                            current_chains = []

                        current_seq = ""
                    else:
                        if current_chains:
                            current_seq += line.strip()

                if current_chains and current_seq:
                    for chain in current_chains:
                        sequences[chain] = current_seq

                return sequences

        except Exception as e:
            logger.warning(f"Failed to download FASTA for {pdb_id}: {str(e)}")
            return {}

    def run_blast(self):
        """
        Execute the BLAST search with selected parameters.

        Returns:
            bool: True if BLAST search was successful, False otherwise
        """
        if not self.sequence:
            logger.error("No sequence provided")
            return False

        try:
            blast_params = {
                "expect": self.e_value,
                "hitlist_size": self.max_hits,
            }

            if self.word_size:
                blast_params["word_size"] = self.word_size

            if self.program == "blastn" and self.megablast:
                blast_params["megablast"] = "on"

            if self.additional_params:
                if (
                    "gapopen" in self.additional_params
                    and "gapextend" in self.additional_params
                ):
                    gap_open = self.additional_params["gapopen"]
                    gap_extend = self.additional_params["gapextend"]
                    blast_params["gapcosts"] = f"{gap_open} {gap_extend}"
                    del self.additional_params["gapopen"]
                    del self.additional_params["gapextend"]

            param_mapping = {
                "matrix": "matrix_name",
                "filter": "filter",
                "query_genetic_code": "genetic_code",
                "db_genetic_code": "db_genetic_code",
                "perc_identity": "perc_ident",
                "entrez_query": "entrez_query",
                "comp_based_stats": "composition_based_statistics",
            }

            for param, value in self.additional_params.items():
                api_param = param_mapping.get(param, param)
                if param in ["comp_based_stats", "template_type"]:
                    blast_params[api_param] = str(value)
                else:
                    blast_params[api_param] = value

            start_time = time.time()
            result_handle = NCBIWWW.qblast(
                program=self.program,
                database=self.database,
                sequence=self.sequence,
                **blast_params,
            )

            self.results = result_handle.read()
            result_handle.close()

            self.elapsed_time = time.time() - start_time
            logger.info(f"BLAST search completed in {self.elapsed_time:.2f} seconds")

            return True

        except Exception as e:
            logger.error(f"Error running BLAST: {str(e)}")
            return False

    def save_results(self, filename=None):
        """
        Save BLAST results to file.

        Args:
            filename: Output filename

        Returns:
            str: Path to the saved results file, or None if failed
        """
        if not self.results:
            logger.error("No results to save")
            return None

        if not filename:
            filename = f"blast_results_{int(time.time())}.xml"

        try:
            with open(filename, "w") as f:
                f.write(
                    self.results.decode("utf-8")
                    if isinstance(self.results, bytes)
                    else self.results
                )

            logger.info(f"Results saved to {filename}")
            return filename

        except Exception as e:
            logger.error(f"Error saving results: {str(e)}")
            return None

    def parse_results(self):
        """
        Parse BLAST results.

        Returns:
            list: List of BLAST records or None if failed
        """
        if not self.results:
            logger.error("No results to parse")
            return None

        try:
            with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp:
                temp_filename = temp.name
                temp.write(
                    self.results.decode("utf-8")
                    if isinstance(self.results, bytes)
                    else self.results
                )

            with open(temp_filename) as result_file:
                blast_records = list(NCBIXML.parse(result_file))

            os.unlink(temp_filename)
            return blast_records

        except Exception as e:
            logger.error(f"Error parsing results: {str(e)}")
            return None

    def get_results_summary(self):
        """
        Get a summary of BLAST results.

        Returns:
            dict: Results summary or None if failed
        """
        blast_records = self.parse_results()
        if not blast_records:
            return None

        record = blast_records[0]

        summary = {
            "query": record.query,
            "database": self.database,
            "program": self.program,
            "num_hits": len(record.alignments),
            "alignments": [],
        }

        for i, alignment in enumerate(record.alignments):
            if i >= self.max_hits:
                break

            for hsp in alignment.hsps:
                identity_pct = (hsp.identities / hsp.align_length) * 100

                positives = getattr(hsp, 'positives', None)
                gaps = getattr(hsp, 'gaps', 0)

                alignment_info = {
                    "title": alignment.title,
                    "length": alignment.length,
                    "e_value": hsp.expect,
                    "score": hsp.score,
                    "bits": hsp.bits,
                    "identities": hsp.identities,
                    "align_length": hsp.align_length,
                    "identity_percent": identity_pct,
                    "query_start": hsp.query_start,
                    "query_end": hsp.query_end,
                    "sbjct_start": hsp.sbjct_start,
                    "sbjct_end": hsp.sbjct_end,
                    "query_seq": hsp.query,
                    "match_seq": hsp.match,
                    "sbjct_seq": hsp.sbjct,
                    "positives": positives,
                    "gaps": gaps,
                }

                if positives is not None:
                    alignment_info["positives_percent"] = (positives / hsp.align_length) * 100

                summary["alignments"].append(alignment_info)
                break

        return summary

    def extract_sequences_from_pdb(self, pdb_file):
        """
        Extract amino acid sequences from PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            dict: Dictionary of chain_id -> sequence
        """
        if not os.path.exists(pdb_file):
            logger.error(f"PDB file not found: {pdb_file}")
            return {}

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_file)

            sequences = {}
            ppb = PPBuilder()

            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    pp_list = ppb.build_peptides(chain)

                    if pp_list:
                        sequence = "".join(str(pp.get_sequence()) for pp in pp_list)
                        sequences[chain_id] = sequence

            return sequences
        except Exception as e:
            logger.error(f"Error extracting sequences from PDB: {str(e)}")
            return {}

    def set_program_defaults(self, program):
        """Set appropriate defaults for a BLAST program."""
        self.program = program

        default_word_sizes = {
            "blastn": 11,
            "blastp": 3,
            "blastx": 3,
            "tblastn": 3,
            "tblastx": 3,
        }

        if program in ["blastn", "tblastn", "tblastx"]:
            self.database = "nt"
        else:
            self.database = "nr"

        self.word_size = default_word_sizes.get(program)

    def set_database_for_program(self, database_choice, program):
        """Set database based on choice and program type."""
        if program in ["blastn", "tblastn", "tblastx"]:
            database_map = {
                "1": "nt",
                "2": "refseq_rna",
                "3": "refseq_genomic",
                "4": "pdbnt",
            }
        else:
            database_map = {
                "1": "nr",
                "2": "refseq_protein",
                "3": "swissprot",
                "4": "pdb",
            }

        self.database = database_map.get(
            database_choice,
            "nt" if program in ["blastn", "tblastn", "tblastx"] else "nr",
        )

    def set_matrix_and_gap_costs(self, matrix_choice):
        """Set scoring matrix and compatible gap costs."""
        matrix_map = {
            "1": "BLOSUM62",
            "2": "BLOSUM45",
            "3": "BLOSUM80",
            "4": "BLOSUM90",
            "5": "PAM30",
            "6": "PAM70",
            "7": "PAM250",
        }

        matrix_gap_costs = {
            "BLOSUM62": [(11, 1), (9, 2), (8, 2), (7, 2), (6, 2), (13, 1)],
            "BLOSUM45": [(14, 2), (13, 3), (12, 3), (11, 3), (10, 3), (15, 2)],
            "BLOSUM80": [(10, 1), (9, 2), (8, 2), (7, 2), (6, 2)],
            "BLOSUM90": [(10, 1), (9, 2), (8, 2), (7, 2), (6, 2)],
            "PAM30": [(9, 1), (7, 2), (6, 2), (5, 2), (10, 1)],
            "PAM70": [(10, 1), (8, 2), (7, 2), (6, 2), (11, 1)],
            "PAM250": [(14, 2), (15, 2), (13, 3), (12, 3), (11, 3)],
        }

        selected_matrix = matrix_map.get(matrix_choice, "BLOSUM62")
        self.additional_params["matrix"] = selected_matrix

        compatible_gaps = matrix_gap_costs.get(selected_matrix, [(11, 1)])
        default_gap = compatible_gaps[0]

        self.additional_params["gapopen"] = default_gap[0]
        self.additional_params["gapextend"] = default_gap[1]

        return compatible_gaps

    def set_genetic_code(self, code_choice):
        """Set genetic code for translated searches."""
        if self.program == "blastx":
            self.additional_params["query_genetic_code"] = int(code_choice)
        elif self.program == "tblastn":
            self.additional_params["db_genetic_code"] = int(code_choice)
        elif self.program == "tblastx":
            self.additional_params["query_genetic_code"] = int(code_choice)
            self.additional_params["db_genetic_code"] = int(code_choice)

    def clear_parameters(self):
        """Reset all parameters to defaults."""
        self.program = "blastn"
        self.database = "nt"
        self.sequence = ""
        self.e_value = 10.0
        self.word_size = None
        self.max_hits = 50
        self.megablast = True
        self.additional_params = {}
        self.results = None
