"""
EMBOSS Core Execution Engine

Core functionality for executing EMBOSS tools and parsing results.
Handles tool execution, file management, and result parsing.
"""

import logging
import os
import tempfile
import subprocess
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import re

logger = logging.getLogger(__name__)


class EMBOSSCore:
    """Core engine for EMBOSS tool execution and result parsing"""

    def __init__(self):
        """Initialize EMBOSS core with tool discovery and validation"""
        self.emboss_available = self._check_emboss_installation()
        self.temp_dir = None
        self.tool_paths = self._discover_emboss_tools()
        
    def _check_emboss_installation(self) -> bool:
        """Check if EMBOSS is installed and accessible"""
        try:
            result = subprocess.run(
                ["embossversion"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            if result.returncode == 0:
                logger.info(f"EMBOSS found: {result.stdout.strip()}")
                return True
            else:
                logger.warning("EMBOSS installation found but embossversion failed")
                return False
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
            logger.warning("EMBOSS not found in PATH")
            return False

    def _discover_emboss_tools(self) -> Dict[str, str]:
        """Discover available EMBOSS tools"""
        tools = {}
        emboss_tools = [
            "pepstats", "pepinfo", "charge", "iep",
            "needle", "water", "stretcher", "diffseq",
            "patmatmotifs", "fuzzpro", "sigcleave",
            "dotplot", "matcher"
        ]
        
        for tool in emboss_tools:
            try:
                result = subprocess.run(
                    [tool, "-help"], 
                    capture_output=True, 
                    text=True, 
                    timeout=5
                )
                if result.returncode == 0 or "Usage:" in result.stderr:
                    tools[tool] = shutil.which(tool) or tool
                    logger.debug(f"Found EMBOSS tool: {tool}")
            except (subprocess.TimeoutExpired, FileNotFoundError):
                logger.debug(f"EMBOSS tool not found: {tool}")
                
        return tools

    def _setup_temp_directory(self) -> str:
        """Create and return temporary directory for EMBOSS operations"""
        if not self.temp_dir or not os.path.exists(self.temp_dir):
            self.temp_dir = tempfile.mkdtemp(prefix="emboss_")
        return self.temp_dir

    def _cleanup_temp_directory(self):
        """Clean up temporary directory"""
        if self.temp_dir and os.path.exists(self.temp_dir):
            try:
                shutil.rmtree(self.temp_dir)
                self.temp_dir = None
            except Exception as e:
                logger.warning(f"Failed to cleanup temp directory: {e}")

    def _write_sequence_file(self, sequence: str, name: str = "sequence") -> str:
        """Write sequence to temporary FASTA file"""
        temp_dir = self._setup_temp_directory()
        seq_file = os.path.join(temp_dir, f"{name}.fasta")
        
        # Clean sequence name for filename safety
        safe_name = re.sub(r'[^\w\-_\.]', '_', name)
        
        with open(seq_file, 'w') as f:
            f.write(f">{safe_name}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")
                
        return seq_file

    def _run_emboss_tool(self, tool_name: str, args: List[str], input_files: List[str] = None) -> Dict[str, Any]:
        """
        Execute an EMBOSS tool with given arguments
        
        Args:
            tool_name: Name of the EMBOSS tool
            args: Command line arguments
            input_files: List of input files to track
            
        Returns:
            Dictionary with execution results
        """
        if not self.emboss_available:
            return {
                "success": False,
                "error": "EMBOSS not available",
                "output": "",
                "stderr": ""
            }
            
        if tool_name not in self.tool_paths:
            return {
                "success": False,
                "error": f"Tool {tool_name} not found",
                "output": "",
                "stderr": ""
            }

        try:
            # Prepare command
            cmd = [self.tool_paths[tool_name]] + args
            
            logger.debug(f"Running EMBOSS command: {' '.join(cmd)}")
            
            # Execute command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout
                cwd=self.temp_dir
            )
            
            return {
                "success": result.returncode == 0,
                "returncode": result.returncode,
                "output": result.stdout,
                "stderr": result.stderr,
                "command": ' '.join(cmd)
            }
            
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "error": f"Tool {tool_name} timed out",
                "output": "",
                "stderr": "Command timed out after 5 minutes"
            }
        except Exception as e:
            return {
                "success": False,
                "error": f"Failed to execute {tool_name}: {str(e)}",
                "output": "",
                "stderr": str(e)
            }

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Sequence Analysis Tools
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def run_pepstats(self, sequence: str, name: str = "sequence") -> Dict[str, Any]:
        """
        Run pepstats to get protein statistics
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            
        Returns:
            Dictionary with pepstats results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_pepstats.txt")
        
        args = [
            "-sequence", input_file,
            "-outfile", output_file,
            "-auto"  # Accept default parameters
        ]
        
        result = self._run_emboss_tool("pepstats", args, [input_file])
        
        if result["success"]:
            # Parse pepstats output
            parsed_results = self._parse_pepstats_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_pepinfo(self, sequence: str, name: str = "sequence") -> Dict[str, Any]:
        """
        Run pepinfo for detailed residue analysis
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            
        Returns:
            Dictionary with pepinfo results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_pepinfo.txt")
        
        args = [
            "-sequence", input_file,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("pepinfo", args, [input_file])
        
        if result["success"]:
            parsed_results = self._parse_pepinfo_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_charge(self, sequence: str, name: str = "sequence") -> Dict[str, Any]:
        """
        Run charge to calculate net charge at different pH values
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            
        Returns:
            Dictionary with charge results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_charge.txt")
        
        args = [
            "-sequence", input_file,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("charge", args, [input_file])
        
        if result["success"]:
            parsed_results = self._parse_charge_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_iep(self, sequence: str, name: str = "sequence") -> Dict[str, Any]:
        """
        Run iep to calculate isoelectric point
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            
        Returns:
            Dictionary with iep results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_iep.txt")
        
        args = [
            "-sequence", input_file,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("iep", args, [input_file])
        
        if result["success"]:
            parsed_results = self._parse_iep_output(output_file)
            result.update(parsed_results)
        
        return result

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Alignment Tools
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def run_needle(self, seq1: str, seq2: str, name1: str = "seq1", name2: str = "seq2") -> Dict[str, Any]:
        """
        Run needle for global sequence alignment
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            name1: First sequence name
            name2: Second sequence name
            
        Returns:
            Dictionary with alignment results
        """
        temp_dir = self._setup_temp_directory()
        input1 = self._write_sequence_file(seq1, name1)
        input2 = self._write_sequence_file(seq2, name2)
        output_file = os.path.join(temp_dir, f"{name1}_vs_{name2}_needle.txt")
        
        args = [
            "-asequence", input1,
            "-bsequence", input2,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("needle", args, [input1, input2])
        
        if result["success"]:
            parsed_results = self._parse_alignment_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_water(self, seq1: str, seq2: str, name1: str = "seq1", name2: str = "seq2") -> Dict[str, Any]:
        """
        Run water for local sequence alignment
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            name1: First sequence name
            name2: Second sequence name
            
        Returns:
            Dictionary with alignment results
        """
        temp_dir = self._setup_temp_directory()
        input1 = self._write_sequence_file(seq1, name1)
        input2 = self._write_sequence_file(seq2, name2)
        output_file = os.path.join(temp_dir, f"{name1}_vs_{name2}_water.txt")
        
        args = [
            "-asequence", input1,
            "-bsequence", input2,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("water", args, [input1, input2])
        
        if result["success"]:
            parsed_results = self._parse_alignment_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_stretcher(self, seq1: str, seq2: str, name1: str = "seq1", name2: str = "seq2") -> Dict[str, Any]:
        """
        Run stretcher for global alignment optimized for long sequences
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            name1: First sequence name
            name2: Second sequence name
            
        Returns:
            Dictionary with alignment results
        """
        temp_dir = self._setup_temp_directory()
        input1 = self._write_sequence_file(seq1, name1)
        input2 = self._write_sequence_file(seq2, name2)
        output_file = os.path.join(temp_dir, f"{name1}_vs_{name2}_stretcher.txt")
        
        args = [
            "-asequence", input1,
            "-bsequence", input2,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("stretcher", args, [input1, input2])
        
        if result["success"]:
            parsed_results = self._parse_alignment_output(output_file)
            result.update(parsed_results)
        
        return result

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Pattern Search Tools
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def run_patmatmotifs(self, sequence: str, name: str = "sequence") -> Dict[str, Any]:
        """
        Run patmatmotifs to search PROSITE motif database
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            
        Returns:
            Dictionary with motif search results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_motifs.txt")
        
        args = [
            "-sequence", input_file,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("patmatmotifs", args, [input_file])
        
        if result["success"]:
            parsed_results = self._parse_motif_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_fuzzpro(self, sequence: str, name: str = "sequence", pattern: str = None) -> Dict[str, Any]:
        """
        Run fuzzpro to search for patterns with mismatches
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            pattern: Pattern to search for
            
        Returns:
            Dictionary with pattern search results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_fuzzpro.txt")
        
        # Use default pattern if none provided
        if not pattern:
            pattern = "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H"  # Zinc finger pattern
        
        args = [
            "-sequence", input_file,
            "-pattern", pattern,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("fuzzpro", args, [input_file])
        
        if result["success"]:
            parsed_results = self._parse_fuzzpro_output(output_file)
            result.update(parsed_results)
        
        return result

    def run_sigcleave(self, sequence: str, name: str = "sequence") -> Dict[str, Any]:
        """
        Run sigcleave to find signal peptide cleavage sites
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            
        Returns:
            Dictionary with signal peptide results
        """
        temp_dir = self._setup_temp_directory()
        input_file = self._write_sequence_file(sequence, name)
        output_file = os.path.join(temp_dir, f"{name}_sigcleave.txt")
        
        args = [
            "-sequence", input_file,
            "-outfile", output_file,
            "-auto"
        ]
        
        result = self._run_emboss_tool("sigcleave", args, [input_file])
        
        if result["success"]:
            parsed_results = self._parse_sigcleave_output(output_file)
            result.update(parsed_results)
        
        return result

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Output Parsers
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def _parse_pepstats_output(self, output_file: str) -> Dict[str, Any]:
        """Parse pepstats output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            properties = {}
            
            # Extract molecular weight
            mw_match = re.search(r'Molecular weight = ([\d.]+)', content)
            if mw_match:
                properties['molecular_weight'] = float(mw_match.group(1))
            
            # Extract length
            length_match = re.search(r'Residues = (\d+)', content)
            if length_match:
                properties['length'] = int(length_match.group(1))
            
            # Extract average residue weight
            avg_match = re.search(r'Average weight = ([\d.]+)', content)
            if avg_match:
                properties['average_weight'] = float(avg_match.group(1))
            
            # Extract charge (if present)
            charge_match = re.search(r'Charge = ([+-]?[\d.]+)', content)
            if charge_match:
                properties['net_charge'] = float(charge_match.group(1))
            
            # Extract isoelectric point (if present)
            pi_match = re.search(r'Isoelectric Point = ([\d.]+)', content)
            if pi_match:
                properties['isoelectric_point'] = float(pi_match.group(1))
            
            return {
                "properties": properties,
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse pepstats output: {e}")
            return {"properties": {}, "parse_error": str(e)}

    def _parse_charge_output(self, output_file: str) -> Dict[str, Any]:
        """Parse charge output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            charges = {}
            
            # Look for charge vs pH data
            for line in content.split('\n'):
                # Match lines like "pH 7.0: charge = -2.50"
                match = re.search(r'pH\s+([\d.]+).*?charge\s*=\s*([+-]?[\d.]+)', line)
                if match:
                    ph = float(match.group(1))
                    charge = float(match.group(2))
                    charges[ph] = charge
            
            return {
                "charges": charges,
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse charge output: {e}")
            return {"charges": {}, "parse_error": str(e)}

    def _parse_iep_output(self, output_file: str) -> Dict[str, Any]:
        """Parse iep output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract isoelectric point
            iep_match = re.search(r'Isoelectric point = ([\d.]+)', content)
            isoelectric_point = float(iep_match.group(1)) if iep_match else None
            
            return {
                "isoelectric_point": isoelectric_point,
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse iep output: {e}")
            return {"isoelectric_point": None, "parse_error": str(e)}

    def _parse_alignment_output(self, output_file: str) -> Dict[str, Any]:
        """Parse alignment output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Extract alignment statistics
            summary = {}
            
            # Score
            score_match = re.search(r'Score:\s*([\d.]+)', content)
            if score_match:
                summary['score'] = float(score_match.group(1))
            
            # Identity
            identity_match = re.search(r'Identity:\s*(\d+)/(\d+)\s*\(([\d.]+)%\)', content)
            if identity_match:
                summary['identity'] = float(identity_match.group(3))
                summary['identity_count'] = int(identity_match.group(1))
                summary['identity_total'] = int(identity_match.group(2))
            
            # Similarity
            similarity_match = re.search(r'Similarity:\s*(\d+)/(\d+)\s*\(([\d.]+)%\)', content)
            if similarity_match:
                summary['similarity'] = float(similarity_match.group(3))
            
            # Extract alignment sequences
            alignment = self._extract_alignment_sequences(content)
            
            return {
                "summary": summary,
                "alignment": alignment,
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse alignment output: {e}")
            return {"summary": {}, "parse_error": str(e)}

    def _parse_motif_output(self, output_file: str) -> Dict[str, Any]:
        """Parse motif search output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            motifs = []
            
            # Parse motif matches
            for line in content.split('\n'):
                # Look for motif match lines
                match = re.search(r'(\d+)\s+(\d+)\s+([A-Z]+)\s+(.+)', line.strip())
                if match:
                    motifs.append({
                        'start': int(match.group(1)),
                        'end': int(match.group(2)),
                        'sequence': match.group(3),
                        'description': match.group(4).strip()
                    })
            
            return {
                "motifs": motifs,
                "motif_count": len(motifs),
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse motif output: {e}")
            return {"motifs": [], "parse_error": str(e)}

    def _parse_fuzzpro_output(self, output_file: str) -> Dict[str, Any]:
        """Parse fuzzpro output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            matches = []
            
            # Parse pattern matches
            for line in content.split('\n'):
                match = re.search(r'(\d+)\s+([A-Z]+)', line.strip())
                if match:
                    matches.append({
                        'position': int(match.group(1)),
                        'sequence': match.group(2)
                    })
            
            return {
                "matches": matches,
                "match_count": len(matches),
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse fuzzpro output: {e}")
            return {"matches": [], "parse_error": str(e)}

    def _parse_sigcleave_output(self, output_file: str) -> Dict[str, Any]:
        """Parse sigcleave output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            cleavage_sites = []
            
            # Parse cleavage site predictions
            for line in content.split('\n'):
                match = re.search(r'(\d+)\s+([\d.]+)', line.strip())
                if match:
                    cleavage_sites.append({
                        'position': int(match.group(1)),
                        'score': float(match.group(2))
                    })
            
            return {
                "cleavage_sites": cleavage_sites,
                "site_count": len(cleavage_sites),
                "raw_output": content
            }
            
        except Exception as e:
            logger.error(f"Failed to parse sigcleave output: {e}")
            return {"cleavage_sites": [], "parse_error": str(e)}

    def _extract_alignment_sequences(self, content: str) -> Dict[str, str]:
        """Extract aligned sequences from alignment output"""
        alignment = {}
        
        try:
            # Look for alignment blocks
            lines = content.split('\n')
            seq1_lines = []
            seq2_lines = []
            
            for line in lines:
                if line.startswith('seq1') or line.startswith('asequence'):
                    seq1_lines.append(line.split()[-1])
                elif line.startswith('seq2') or line.startswith('bsequence'):
                    seq2_lines.append(line.split()[-1])
            
            alignment['sequence1'] = ''.join(seq1_lines)
            alignment['sequence2'] = ''.join(seq2_lines)
            
        except Exception as e:
            logger.error(f"Failed to extract alignment sequences: {e}")
            
        return alignment

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Utility Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def get_timestamp(self) -> str:
        """Get current timestamp for results"""
        return datetime.now().isoformat()

    def is_available(self) -> bool:
        """Check if EMBOSS is available"""
        return self.emboss_available

    def get_available_tools(self) -> List[str]:
        """Get list of available EMBOSS tools"""
        return list(self.tool_paths.keys())

    def __del__(self):
        """Cleanup on destruction"""
        self._cleanup_temp_directory()
