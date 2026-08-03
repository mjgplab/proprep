"""
AlphaFold Worker

Core functionality for running AlphaFold predictions via ColabFold.
Handles GPU detection, sequence input, parameter configuration, and output parsing.
"""

import logging
import os
import subprocess
import json
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio.PDB import PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

logger = logging.getLogger(__name__)


class AlphaFoldWorker:
    """Handles AlphaFold predictions with detailed parameter control."""

    def __init__(self):
        """Initialize the AlphaFold worker."""
        self.backend = self._detect_backend()
        self.gpu_info = self._detect_gpu()

        # Default parameters (aligned with ColabFold defaults)
        self.num_models = 5
        self.num_recycle = 3
        self.use_templates = True
        self.amber_relax = True
        self.use_gpu_relax = True
        self.max_msa = "auto"
        self.msa_mode = "mmseqs2"
        self.pair_mode = "unpaired"
        self.random_seed = None
        self.save_all_models = True
        self.save_msa = False
        self.save_pae_json = True
        self.custom_template_path = None

        # Results storage
        self.results = None
        self.output_dir = None

    def _detect_gpu(self) -> Dict[str, Any]:
        """
        Detect GPU availability (NVIDIA or Apple Silicon).

        Returns:
            Dict with keys: has_gpu, gpu_type, gpu_details
        """
        gpu_info = {
            "has_gpu": False,
            "gpu_type": None,
            "gpu_details": None
        }

        # Check NVIDIA GPU
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0 and result.stdout.strip():
                gpu_info["has_gpu"] = True
                gpu_info["gpu_type"] = "NVIDIA"
                gpu_info["gpu_details"] = result.stdout.strip()
                logger.debug(f"Detected NVIDIA GPU: {gpu_info['gpu_details']}")
                return gpu_info
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        # Check Apple Silicon GPU
        try:
            result = subprocess.run(
                ['system_profiler', 'SPDisplaysDataType'],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                output = result.stdout
                # Look for M1/M2/M3/M4 chips or Apple GPU
                if any(chip in output for chip in ['M1', 'M2', 'M3', 'M4', 'Apple']):
                    gpu_info["has_gpu"] = True
                    gpu_info["gpu_type"] = "Apple Silicon"

                    # Extract chip model from output
                    for line in output.split('\n'):
                        if 'Chipset Model' in line or 'Chip Model' in line:
                            gpu_info["gpu_details"] = line.split(':')[1].strip()
                            break

                    if not gpu_info["gpu_details"]:
                        # Fallback: look for any Apple chip mention
                        for line in output.split('\n'):
                            if any(chip in line for chip in ['M1', 'M2', 'M3', 'M4']):
                                gpu_info["gpu_details"] = line.strip()
                                break

                    logger.debug(f"Detected Apple Silicon GPU: {gpu_info['gpu_details']}")
                    return gpu_info
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        logger.debug("No GPU detected, will use CPU")
        return gpu_info

    def _detect_backend(self) -> str:
        """
        Check if ColabFold is installed locally.

        Returns:
            'local' if colabfold_batch is available, 'colab' otherwise
        """
        try:
            result = subprocess.run(
                ['colabfold_batch', '--help'],
                capture_output=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.debug("Local ColabFold installation detected")
                return 'local'
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        logger.debug("ColabFold not found locally, will use Google Colab instructions")
        return 'colab'

    def setup(self, **kwargs):
        """
        Set up prediction parameters.

        Args:
            **kwargs: Parameter overrides (num_models, num_recycle, etc.)
        """
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                logger.warning(f"Unknown parameter: {key}")

    def extract_sequences_from_pdb(self, pdb_file: str) -> Dict[str, str]:
        """
        Extract amino acid sequences from a PDB file.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dictionary mapping chain IDs to sequences
        """
        if not os.path.exists(pdb_file):
            logger.error(f"PDB file not found: {pdb_file}")
            return {}

        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", pdb_file)

            ppb = PPBuilder()
            sequences = {}

            for model in structure:
                for chain in model:
                    chain_id = chain.get_id()
                    # Build polypeptide sequences
                    peptides = ppb.build_peptides(chain)
                    if peptides:
                        seq = ''.join([str(pp.get_sequence()) for pp in peptides])
                        if seq:  # Only add if sequence is not empty
                            sequences[chain_id] = seq

            logger.info(f"Extracted {len(sequences)} sequence(s) from {pdb_file}")
            return sequences

        except Exception as e:
            logger.error(f"Error extracting sequences from PDB: {str(e)}")
            return {}

    def validate_sequence(self, sequence: str) -> Tuple[bool, str]:
        """
        Validate an amino acid sequence.

        Args:
            sequence: Amino acid sequence string

        Returns:
            Tuple of (is_valid, error_message)
        """
        if not sequence:
            return False, "Sequence is empty"

        # Remove whitespace and newlines
        sequence = sequence.replace('\n', '').replace(' ', '').strip()

        # Check if all characters are valid amino acids
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        sequence_set = set(sequence.upper())

        invalid_chars = sequence_set - valid_aa
        if invalid_chars:
            return False, f"Invalid characters in sequence: {invalid_chars}"

        # Check minimum length
        if len(sequence) < 10:
            return False, "Sequence too short (minimum 10 residues)"

        # Check maximum length (ColabFold limit ~3000)
        if len(sequence) > 3000:
            return False, f"Sequence too long ({len(sequence)} residues, maximum 3000)"

        return True, ""

    def predict(self, sequence: str, output_dir: str, job_name: str = "alphafold_prediction",
                **params) -> Dict[str, Any]:
        """
        Run AlphaFold prediction.

        Args:
            sequence: Amino acid sequence
            output_dir: Directory for output files
            job_name: Name for this prediction job
            **params: Additional parameters to override defaults

        Returns:
            Dictionary with results including:
                - best_model_pdb: Path to best ranked model
                - all_models: List of all model PDB files
                - confidence_metrics: Dict with pLDDT, pTM scores
                - pae_matrix: Predicted aligned error data
                - output_dir: Full output directory path
                - structure_object: BioPython Structure object
        """
        # Validate sequence
        is_valid, error_msg = self.validate_sequence(sequence)
        if not is_valid:
            raise ValueError(f"Invalid sequence: {error_msg}")

        # Update parameters with any overrides
        self.setup(**params)

        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        self.output_dir = str(output_path)

        if self.backend == 'local':
            return self._run_local(sequence, output_path, job_name)
        else:
            return self._provide_colab_instructions(sequence, job_name)

    def _run_local(self, sequence: str, output_dir: Path, job_name: str) -> Dict[str, Any]:
        """
        Run ColabFold locally.

        Args:
            sequence: Amino acid sequence
            output_dir: Output directory path
            job_name: Job name

        Returns:
            Results dictionary
        """
        # Write sequence to FASTA file
        fasta_file = output_dir / f"{job_name}.fasta"
        with open(fasta_file, 'w') as f:
            f.write(f">{job_name}\n{sequence}\n")

        logger.info(f"Wrote sequence to {fasta_file}")

        # Build colabfold_batch command
        cmd = [
            "colabfold_batch",
            str(fasta_file),
            str(output_dir),
            "--num-models", str(self.num_models),
            "--num-recycle", str(self.num_recycle),
        ]

        # Add optional parameters
        if self.amber_relax:
            cmd.append("--amber")

        if self.use_gpu_relax and self.gpu_info["has_gpu"]:
            cmd.append("--use-gpu-relax")

        if self.use_templates:
            cmd.append("--templates")

        if self.custom_template_path:
            cmd.extend(["--custom-template-path", str(self.custom_template_path)])

        if self.random_seed is not None:
            cmd.extend(["--random-seed", str(self.random_seed)])

        if self.max_msa != "auto":
            cmd.extend(["--max-msa", str(self.max_msa)])

        if not self.save_all_models:
            # ColabFold doesn't have a direct flag for this,
            # we'll handle cleanup after
            pass

        # Add overwrite flag to avoid prompts
        cmd.append("--overwrite-existing-results")

        logger.info(f"Running ColabFold command: {' '.join(cmd)}")

        # Run the command
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )

            if result.returncode != 0:
                logger.error(f"ColabFold failed with return code {result.returncode}")
                logger.error(f"STDERR: {result.stderr}")
                raise RuntimeError(f"ColabFold prediction failed: {result.stderr}")

            logger.info("ColabFold prediction completed successfully")
            logger.debug(f"STDOUT: {result.stdout}")

        except subprocess.TimeoutExpired:
            raise RuntimeError("ColabFold prediction timed out after 2 hours")

        # Parse results
        return self._parse_results(output_dir, job_name)

    def _parse_results(self, output_dir: Path, job_name: str) -> Dict[str, Any]:
        """
        Parse ColabFold output files.

        Args:
            output_dir: Output directory path
            job_name: Job name

        Returns:
            Results dictionary
        """
        results = {
            "best_model_pdb": None,
            "all_models": [],
            "confidence_metrics": {},
            "pae_matrix": None,
            "output_dir": str(output_dir),
            "structure_object": None
        }

        # Find all output files
        # ColabFold naming: {job_name}_relaxed_rank_001_*.pdb
        pdb_files = list(output_dir.glob(f"{job_name}*.pdb"))

        if not pdb_files:
            raise RuntimeError(f"No PDB files found in {output_dir}")

        # Find best model (rank_001 with relaxed if available)
        best_model = None
        for pdb_file in pdb_files:
            if "rank_001" in pdb_file.name or "rank_1" in pdb_file.name:
                if "relaxed" in pdb_file.name:
                    best_model = pdb_file
                    break
                elif best_model is None:
                    best_model = pdb_file

        if best_model is None:
            # Fallback to first file
            best_model = sorted(pdb_files)[0]

        results["best_model_pdb"] = str(best_model)
        results["all_models"] = [str(f) for f in sorted(pdb_files)]

        logger.info(f"Found {len(pdb_files)} model(s), best: {best_model.name}")

        # Parse confidence scores from JSON files
        json_files = list(output_dir.glob(f"{job_name}*scores*.json"))
        if json_files:
            for json_file in json_files:
                if "rank_001" in json_file.name or "rank_1" in json_file.name:
                    try:
                        with open(json_file, 'r') as f:
                            scores = json.load(f)
                            results["confidence_metrics"] = scores
                            logger.info(f"Loaded confidence scores from {json_file.name}")
                            break
                    except Exception as e:
                        logger.warning(f"Could not parse {json_file}: {e}")

        # Parse PAE matrix
        pae_files = list(output_dir.glob(f"{job_name}*pae*.json"))
        if pae_files:
            for pae_file in pae_files:
                if "rank_001" in pae_file.name or "rank_1" in pae_file.name:
                    try:
                        with open(pae_file, 'r') as f:
                            pae_data = json.load(f)
                            results["pae_matrix"] = pae_data
                            logger.info(f"Loaded PAE matrix from {pae_file.name}")
                            break
                    except Exception as e:
                        logger.warning(f"Could not parse {pae_file}: {e}")

        # Load structure object
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(job_name, str(best_model))
            results["structure_object"] = structure
            logger.info("Loaded BioPython Structure object")
        except Exception as e:
            logger.error(f"Could not load structure: {e}")

        return results

    def _provide_colab_instructions(self, sequence: str, job_name: str) -> Dict[str, Any]:
        """
        Provide instructions for running on Google Colab.

        Args:
            sequence: Amino acid sequence
            job_name: Job name

        Returns:
            Dictionary with instructions
        """
        # Format sequence for display (truncate if very long)
        if len(sequence) > 100:
            seq_display = sequence[:100] + f"... ({len(sequence)} total residues)"
        else:
            seq_display = sequence

        instructions = f"""
╔════════════════════════════════════════════════════════════════╗
║          ColabFold Not Installed Locally                       ║
╚════════════════════════════════════════════════════════════════╝

To run AlphaFold prediction, use Google Colab:

1. Go to: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

2. In the "query_sequence" field, paste:
   >{job_name}
   {seq_display}

3. Configure parameters in the notebook:
   - num_models: {self.num_models}
   - num_recycle: {self.num_recycle}
   - template_mode: "{"pdb100" if self.use_templates else "none"}"
   - use_amber: {self.amber_relax}
   - msa_mode: "{self.msa_mode}"
   - pair_mode: "{self.pair_mode}"
   - max_msa: "{self.max_msa}"

4. Run all cells (Runtime → Run all) and wait for completion
   - This typically takes 10-60 minutes depending on sequence length
   - GPU is provided free by Google Colab

5. Download the results:
   - Look for files ending in "_relaxed_rank_001_*.pdb"
   - Download the best model (rank 001)

6. Load the downloaded PDB file in ProPrep:
   - Use "PDB Loader" → "Load from file"
   - Or drag and drop into ProPrep

Alternative: Install ColabFold locally for faster predictions:
   pip install colabfold[alphafold]>=1.5.5

   Or visit: https://github.com/YoshitakaMo/localcolabfold

Note: Your configured parameters have been applied to these instructions.
      You can reconfigure them using "Configure prediction parameters".
"""

        return {
            "backend": "colab",
            "instructions": instructions,
            "sequence": sequence,
            "job_name": job_name,
            "parameters": {
                "num_models": self.num_models,
                "num_recycle": self.num_recycle,
                "use_templates": self.use_templates,
                "amber_relax": self.amber_relax,
                "msa_mode": self.msa_mode,
                "pair_mode": self.pair_mode,
                "max_msa": self.max_msa
            }
        }

    def get_status(self) -> Dict[str, Any]:
        """
        Get current status of the AlphaFold worker.

        Returns:
            Status dictionary
        """
        return {
            "backend": self.backend,
            "gpu_info": self.gpu_info,
            "has_results": self.results is not None,
            "output_dir": self.output_dir
        }
