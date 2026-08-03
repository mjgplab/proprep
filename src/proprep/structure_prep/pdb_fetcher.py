"""
RCSB PDB Fetcher

Specialist module for downloading and processing structures from RCSB Protein Data Bank.
Extracted from pdb_loader.py to provide reusable PDB download and metadata extraction.

Author: ProPrep Development Team
Date: 2025-11-08
"""

import logging
import os
import urllib.request
from urllib.error import HTTPError
from typing import Optional, Dict, Any, List

logger = logging.getLogger(__name__)


class PDBFetcher:
    """
    Handles downloading and metadata extraction from RCSB Protein Data Bank.

    This class provides utilities for:
    - Downloading PDB files by ID from RCSB
    - Extracting metadata from PDB headers
    - Validating PDB files
    """

    def __init__(self, output_dir: str = "."):
        """
        Initialize the PDB fetcher.

        Args:
            output_dir: Directory to store downloaded PDB files
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def download_by_id(self, pdb_id: str, overwrite: bool = False) -> Optional[str]:
        """
        Download a PDB file from RCSB PDB database.

        Args:
            pdb_id: PDB identifier (e.g., '1ABC', '7xyz')
            overwrite: If True, download even if file exists

        Returns:
            Path to downloaded PDB file, or None if download failed
        """
        pdb_id = pdb_id.upper().strip()
        output_file = os.path.join(self.output_dir, f"{pdb_id}.pdb")

        # Check if file already exists
        if os.path.exists(output_file) and not overwrite:
            logger.info(f"Using existing PDB file: {output_file}")
            return output_file

        logger.info(f"Downloading PDB file {pdb_id} from RCSB...")
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

        try:
            urllib.request.urlretrieve(pdb_url, output_file)
            logger.info(f"Downloaded PDB file to {output_file}")
            return output_file

        except HTTPError as e:
            if e.code == 404:
                logger.error(f"PDB {pdb_id} not found (404). May be a large structure only available in mmCIF format.")
                logger.warning("RCSB only provides mmCIF format for structures with >99,999 atoms.")
                logger.warning("ProPrep's metal site transformation currently requires PDB format.")
                return None
            else:
                logger.error(f"HTTP error downloading PDB file: {e}")
                return None

        except Exception as e:
            logger.error(f"Error downloading PDB file {pdb_id}: {e}")
            return None

    def download_biological_assembly(
        self, pdb_id: str, assembly_num: int = 1, overwrite: bool = False
    ) -> Optional[str]:
        """
        Download a biological assembly PDB file from RCSB PDB database.

        RCSB provides pre-built biological assemblies that don't require
        applying transformation matrices locally.

        Args:
            pdb_id: PDB identifier (e.g., '1ABC', '7xyz')
            assembly_num: Assembly number (1, 2, etc.). Default is 1.
            overwrite: If True, download even if file exists

        Returns:
            Path to downloaded PDB file, or None if download failed
        """
        import gzip
        import shutil

        pdb_id = pdb_id.upper().strip()
        output_file = os.path.join(self.output_dir, f"{pdb_id}_assembly{assembly_num}.pdb")

        # Check if file already exists
        if os.path.exists(output_file) and not overwrite:
            logger.info(f"Using existing biological assembly file: {output_file}")
            return output_file

        logger.info(f"Downloading biological assembly {assembly_num} for {pdb_id} from RCSB...")

        # RCSB provides biological assemblies as gzipped PDB files
        # Format: https://files.rcsb.org/download/{pdb_id}.pdb{assembly_num}.gz
        assembly_url = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_num}.gz"
        temp_gz_file = output_file + ".gz"

        try:
            # Download the gzipped file
            urllib.request.urlretrieve(assembly_url, temp_gz_file)

            # Decompress
            with gzip.open(temp_gz_file, 'rb') as f_in:
                with open(output_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Remove temp gzip file
            os.remove(temp_gz_file)

            logger.info(f"Downloaded biological assembly to {output_file}")
            return output_file

        except HTTPError as e:
            if e.code == 404:
                logger.error(
                    f"Biological assembly {assembly_num} for {pdb_id} not found (404). "
                    f"This assembly may not exist or may only be available in mmCIF format."
                )
                return None
            else:
                logger.error(f"HTTP error downloading biological assembly: {e}")
                return None

        except Exception as e:
            logger.error(f"Error downloading biological assembly {pdb_id}.{assembly_num}: {e}")
            # Clean up temp file if it exists
            if os.path.exists(temp_gz_file):
                os.remove(temp_gz_file)
            return None

    def get_available_assemblies(self, pdb_id: str) -> List[int]:
        """
        Check which biological assemblies are available for a PDB ID.

        Args:
            pdb_id: PDB identifier

        Returns:
            List of available assembly numbers (e.g., [1, 2, 3])
        """
        pdb_id = pdb_id.upper().strip()
        available = []

        # Try assembly numbers 1-10 (most structures have fewer)
        for assembly_num in range(1, 11):
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_num}.gz"
            try:
                # Use HEAD request to check if file exists without downloading
                request = urllib.request.Request(url, method='HEAD')
                urllib.request.urlopen(request, timeout=5)
                available.append(assembly_num)
            except HTTPError as e:
                if e.code == 404:
                    # No more assemblies
                    break
                # Other errors - continue checking
            except Exception:
                # Timeout or other error - continue checking
                continue

        return available

    def get_metadata(self, pdb_file: str) -> Dict[str, Any]:
        """
        Extract basic metadata from a PDB file.

        This is a lightweight version that extracts essential metadata.
        For full metadata extraction, use PDBMetadataExtractor from pdb_loader.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary with basic metadata
        """
        if not os.path.exists(pdb_file):
            logger.error(f"PDB file not found: {pdb_file}")
            return {}

        metadata = {
            "file": pdb_file,
            "header": {},
            "title": "",
            "experimental_method": "",
            "resolution": None,
        }

        try:
            with open(pdb_file, "r") as f:
                title_lines = []

                for line in f:
                    # HEADER record
                    if line.startswith("HEADER"):
                        classification = line[10:50].strip()
                        deposition_date = line[50:59].strip()
                        pdb_id = line[62:66].strip()

                        metadata["header"] = {
                            "classification": classification,
                            "deposition_date": deposition_date,
                            "pdb_id": pdb_id,
                        }

                    # TITLE records
                    elif line.startswith("TITLE"):
                        title_text = line[10:].strip()
                        title_lines.append(title_text)

                    # EXPDTA record (experimental method)
                    elif line.startswith("EXPDTA"):
                        metadata["experimental_method"] = line[10:].strip()

                    # REMARK   2 RESOLUTION record
                    elif line.startswith("REMARK   2 RESOLUTION"):
                        try:
                            # Extract resolution value
                            res_str = line.split()[-2]  # Get second-to-last word
                            metadata["resolution"] = float(res_str)
                        except (ValueError, IndexError):
                            pass

                # Combine title lines
                if title_lines:
                    metadata["title"] = " ".join(title_lines)

            logger.debug(f"Extracted metadata from {pdb_file}")
            return metadata

        except Exception as e:
            logger.error(f"Error extracting metadata from {pdb_file}: {e}")
            return metadata

    def validate_pdb_file(self, pdb_file: str) -> tuple[bool, str]:
        """
        Validate that a file is a valid PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Tuple of (is_valid, error_message)
        """
        if not os.path.exists(pdb_file):
            return False, f"File not found: {pdb_file}"

        if not pdb_file.endswith(('.pdb', '.PDB', '.ent', '.ENT')):
            return False, f"File does not have PDB extension: {pdb_file}"

        try:
            with open(pdb_file, "r") as f:
                # Read first few lines
                has_header = False
                has_atoms = False

                for i, line in enumerate(f):
                    if i > 1000:  # Check first 1000 lines
                        break

                    if line.startswith(("HEADER", "TITLE", "REMARK")):
                        has_header = True

                    if line.startswith("ATOM"):
                        has_atoms = True
                        break

                if not has_atoms:
                    return False, "No ATOM records found in file"

                return True, ""

        except Exception as e:
            return False, f"Error reading file: {e}"
