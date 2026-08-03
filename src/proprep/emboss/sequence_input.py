"""
Sequence Input Manager

Handles flexible sequence input from multiple sources:
- Direct paste/manual entry
- FASTA file upload
- Workspace sequences (PDB, BLAST results, etc.)
- UniProt ID fetching
- Batch processing from multiple sources

Provides intelligent workspace scanning and unified interface.
"""

import logging
import re
import requests
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.table import Table
from rich.panel import Panel
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

logger = logging.getLogger(__name__)


class SequenceInputManager:
    """Manages flexible sequence input from multiple sources"""

    def __init__(self, module):
        """
        Initialize sequence input manager
        
        Args:
            module: Parent EMBOSS module instance
        """
        self.module = module
        self.console = module.processor.console if hasattr(module, 'processor') else None

    def get_sequence_input(self, label: str = "Sequence") -> Optional[Tuple[str, str, str]]:
        """
        Get sequence input from user with maximum flexibility
        
        Args:
            label: Label for the sequence (e.g., "First Sequence")
            
        Returns:
            Tuple of (sequence, name, source) or None if cancelled
        """
        if not self.console:
            return None

        # Scan workspace for available sequences
        available_sequences = self._scan_workspace_sequences()

        # Build options menu
        options = {
            "1": "Paste sequence manually",
            "2": "Load from FASTA file",
            "3": "Fetch from UniProt ID"
        }

        # Add workspace options
        workspace_start = 4
        for i, (name, seq_info) in enumerate(available_sequences.items(), workspace_start):
            seq_len = len(seq_info['sequence'])
            source = seq_info['source']
            options[str(i)] = f"Use {name} ({seq_len} residues) - {source}"

        # Display options
        self.console.print(f"\n[bold]Select source for {label}:[/bold]")
        for key, desc in options.items():
            self.console.print(f"  {key}. {desc}")

        # Get user choice
        choice = prompt_with_context(
            self.module.processor,
            "Select option",
            choices=list(options.keys()),
            default="1",
            module="EMBOSS - Sequence Input",
            description=f"Select source for {label}",
            options_map=options,
        )

        # Handle choice
        return self._handle_sequence_choice(choice, available_sequences, workspace_start)

    def _handle_sequence_choice(self, choice: str, available_sequences: Dict, workspace_start: int) -> Optional[Tuple[str, str, str]]:
        """Handle user's sequence source choice"""
        
        if choice == "1":
            return self._get_manual_sequence()
        elif choice == "2":
            return self._get_fasta_file_sequence()
        elif choice == "3":
            return self._get_uniprot_sequence()
        else:
            # Workspace sequence choice
            try:
                seq_index = int(choice) - workspace_start
                seq_names = list(available_sequences.keys())
                if 0 <= seq_index < len(seq_names):
                    seq_name = seq_names[seq_index]
                    seq_info = available_sequences[seq_name]
                    return seq_info['sequence'], seq_name, seq_info['source']
            except (ValueError, IndexError):
                pass
                
        self.console.print("[red]Invalid choice[/red]")
        return None

    def _get_manual_sequence(self) -> Optional[Tuple[str, str, str]]:
        """Get sequence from manual paste/entry"""
        self.console.print("\n[bold]Manual Sequence Entry[/bold]")
        self.console.print("Paste your sequence (FASTA format or raw sequence):")
        self.console.print("[grey50]Press Enter twice when done[/grey50]")
        
        lines = []
        empty_lines = 0
        
        while empty_lines < 2:
            try:
                line = input()
                if line.strip():
                    lines.append(line.strip())
                    empty_lines = 0
                else:
                    empty_lines += 1
            except EOFError:
                break
        
        if not lines:
            self.console.print("[yellow]No sequence entered[/yellow]")
            return None
        
        # Parse input
        sequence, name = self._parse_sequence_input('\n'.join(lines))
        
        if not sequence:
            self.console.print("[red]No valid sequence found[/red]")
            return None
        
        # Get sequence name if not from FASTA header
        if not name:
            name = prompt_with_context(
                self.module.processor,
                "Enter sequence name",
                default="manual_sequence",
                module="EMBOSS - Sequence Input",
                description="Name for manually entered sequence",
            )
        
        self.console.print(f"[green]✓ Sequence loaded: {name} ({len(sequence)} residues)[/green]")
        return sequence, name, "Manual Entry"

    def _get_fasta_file_sequence(self) -> Optional[Tuple[str, str, str]]:
        """Get sequence from FASTA file"""
        file_path = prompt_with_context(
            self.module.processor,
            "Enter FASTA file path",
            module="EMBOSS - Sequence Input",
            description="FASTA file path",
        )
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            sequences = self._parse_fasta_content(content)
            
            if not sequences:
                self.console.print("[red]No valid sequences found in file[/red]")
                return None
            
            if len(sequences) == 1:
                # Single sequence
                seq, name = sequences[0]
                self.console.print(f"[green]✓ Loaded: {name} ({len(seq)} residues)[/green]")
                return seq, name, f"File: {Path(file_path).name}"
            else:
                # Multiple sequences - let user choose
                self.console.print(f"\n[bold]Found {len(sequences)} sequences:[/bold]")
                for i, (seq, name) in enumerate(sequences, 1):
                    self.console.print(f"  {i}. {name} ({len(seq)} residues)")
                
                choices = [str(i) for i in range(1, len(sequences) + 1)]
                fasta_options_map = {
                    str(i): f"{name} ({len(seq)} residues)"
                    for i, (seq, name) in enumerate(sequences, 1)
                }
                choice = prompt_with_context(
                    self.module.processor,
                    "Select sequence",
                    choices=choices,
                    default="1",
                    module="EMBOSS - Sequence Input",
                    description=f"Select sequence from FASTA file ({Path(file_path).name})",
                    options_map=fasta_options_map,
                )
                
                seq_index = int(choice) - 1
                seq, name = sequences[seq_index]
                self.console.print(f"[green]✓ Selected: {name} ({len(seq)} residues)[/green]")
                return seq, name, f"File: {Path(file_path).name}"
                
        except FileNotFoundError:
            self.console.print(f"[red]File not found: {file_path}[/red]")
        except Exception as e:
            self.console.print(f"[red]Error reading file: {str(e)}[/red]")
            
        return None

    def _get_uniprot_sequence(self) -> Optional[Tuple[str, str, str]]:
        """Get sequence from UniProt ID"""
        uniprot_id = prompt_with_context(
            self.module.processor,
            "Enter UniProt ID",
            module="EMBOSS - Sequence Input",
            description="UniProt ID to fetch sequence",
        )
        
        try:
            # Fetch from UniProt
            url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                sequences = self._parse_fasta_content(response.text)
                if sequences:
                    seq, name = sequences[0]
                    self.console.print(f"[green]✓ Downloaded: {name} ({len(seq)} residues)[/green]")
                    return seq, name, f"UniProt: {uniprot_id}"
            
            self.console.print(f"[red]Failed to fetch UniProt entry: {uniprot_id}[/red]")
            
        except requests.RequestException as e:
            self.console.print(f"[red]Network error: {str(e)}[/red]")
        except Exception as e:
            self.console.print(f"[red]Error fetching UniProt sequence: {str(e)}[/red]")
            
        return None

    def _scan_workspace_sequences(self) -> Dict[str, Dict[str, str]]:
        """
        Scan workspace for available sequences
        
        Returns:
            Dictionary mapping sequence names to {sequence, source} info
        """
        sequences = {}
        
        if not hasattr(self.module, 'get_from_workspace'):
            return sequences

        # Check for PDB structure sequences
        structure = self.module.get_from_workspace("structure")
        if structure:
            pdb_sequences = self._extract_pdb_sequences(structure)
            for chain_id, seq in pdb_sequences.items():
                sequences[f"PDB Chain {chain_id}"] = {
                    'sequence': seq,
                    'source': 'PDB Structure'
                }

        # Check for filtered structure sequences
        filtered_structure = self.module.get_from_workspace("filtered_structure")
        if filtered_structure:
            filtered_sequences = self._extract_pdb_sequences(filtered_structure)
            for chain_id, seq in filtered_sequences.items():
                sequences[f"Filtered Chain {chain_id}"] = {
                    'sequence': seq,
                    'source': 'Filtered PDB'
                }

        # Check for BLAST results
        blast_results = self.module.get_from_workspace("blast_results")
        if blast_results and isinstance(blast_results, list):
            for i, hit in enumerate(blast_results[:5]):  # Top 5 hits
                if 'sequence' in hit and hit['sequence']:
                    hit_id = hit.get('id', f'hit_{i+1}')
                    sequences[f"BLAST Hit {i+1} ({hit_id})"] = {
                        'sequence': hit['sequence'],
                        'source': 'BLAST Search'
                    }

        # Check for amino acid mutator sequences
        mutated_sequence = self.module.get_from_workspace("mutated_sequence")
        if mutated_sequence:
            sequences["Mutated Sequence"] = {
                'sequence': mutated_sequence,
                'source': 'Amino Acid Mutator'
            }

        # Check for any custom sequences stored in workspace
        custom_sequences = self.module.get_from_workspace("custom_sequences", {})
        for name, seq_data in custom_sequences.items():
            if isinstance(seq_data, dict) and 'sequence' in seq_data:
                sequences[name] = {
                    'sequence': seq_data['sequence'],
                    'source': seq_data.get('source', 'Custom')
                }
            elif isinstance(seq_data, str):
                sequences[name] = {
                    'sequence': seq_data,
                    'source': 'Custom'
                }

        return sequences

    def get_batch_sequences(self) -> List[Tuple[str, str, str]]:
        """
        Get multiple sequences for batch processing
        
        Returns:
            List of (name, sequence, source) tuples
        """
        sequences = []
        
        # Get all available sequences from workspace
        available_sequences = self._scan_workspace_sequences()
        
        if available_sequences:
            self.console.print("\n[bold]Available sequences for batch processing:[/bold]")
            seq_list = list(available_sequences.items())
            
            for i, (name, seq_info) in enumerate(seq_list, 1):
                seq_len = len(seq_info['sequence'])
                source = seq_info['source']
                self.console.print(f"  {i}. {name} ({seq_len} residues) - {source}")
            
            # Ask user which sequences to include
            if confirm_with_context(
                self.module.processor,
                "\nInclude all workspace sequences?",
                default=True,
                module="EMBOSS - Batch Input",
                description="Include all workspace sequences in batch",
            ):
                for name, seq_info in available_sequences.items():
                    sequences.append((name, seq_info['sequence'], seq_info['source']))
            else:
                # Let user select specific sequences
                choices = prompt_with_context(
                    self.module.processor,
                    "Enter sequence numbers (comma-separated)",
                    default="all",
                    module="EMBOSS - Batch Input",
                    description="Specific workspace sequence numbers to include",
                )
                
                if choices.lower() != "all":
                    try:
                        indices = [int(x.strip()) - 1 for x in choices.split(',')]
                        for i in indices:
                            if 0 <= i < len(seq_list):
                                name, seq_info = seq_list[i]
                                sequences.append((name, seq_info['sequence'], seq_info['source']))
                    except (ValueError, IndexError):
                        self.console.print("[red]Invalid selection, using all sequences[/red]")
                        for name, seq_info in available_sequences.items():
                            sequences.append((name, seq_info['sequence'], seq_info['source']))
                else:
                    for name, seq_info in available_sequences.items():
                        sequences.append((name, seq_info['sequence'], seq_info['source']))

        # Offer to add additional sequences
        if confirm_with_context(
            self.module.processor,
            "\nAdd additional sequences from files?",
            default=False,
            module="EMBOSS - Batch Input",
            description="Add additional sequences from files to batch",
        ):
            while True:
                additional_seq = self._get_fasta_file_sequence()
                if additional_seq:
                    sequences.append(additional_seq)
                    if not confirm_with_context(
                        self.module.processor,
                        "Add another sequence?",
                        default=False,
                        module="EMBOSS - Batch Input",
                        description="Add another sequence to batch",
                    ):
                        break
                else:
                    break

        return sequences

    def _extract_pdb_sequences(self, structure) -> Dict[str, str]:
        """
        Extract protein sequences from PDB structure
        
        Args:
            structure: Bio.PDB Structure object
            
        Returns:
            Dictionary mapping chain IDs to sequences
        """
        sequences = {}
        
        try:
            # Standard amino acid mapping
            aa_map = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                'MSE': 'M'  # Selenomethionine
            }
            
            for model in structure:
                for chain in model:
                    sequence = []
                    for residue in chain:
                        if residue.id[0] == ' ':  # Standard residue
                            resname = residue.resname
                            if resname in aa_map:
                                sequence.append(aa_map[resname])
                            else:
                                sequence.append('X')  # Unknown residue
                    
                    if sequence:
                        sequences[chain.id] = ''.join(sequence)
                
                break  # Only process first model
                
        except Exception as e:
            logger.error(f"Failed to extract PDB sequences: {e}")
            
        return sequences

    def _parse_sequence_input(self, input_text: str) -> Tuple[str, str]:
        """
        Parse sequence input (FASTA or raw sequence)
        
        Args:
            input_text: Input text from user
            
        Returns:
            Tuple of (sequence, name)
        """
        input_text = input_text.strip()
        
        if input_text.startswith('>'):
            # FASTA format
            lines = input_text.split('\n')
            name = lines[0][1:].strip()
            sequence = ''.join(lines[1:])
        else:
            # Raw sequence
            name = ""
            sequence = input_text
        
        # Clean sequence - remove whitespace and non-amino acid characters
        sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
        
        return sequence, name

    def _parse_fasta_content(self, content: str) -> List[Tuple[str, str]]:
        """
        Parse FASTA file content
        
        Args:
            content: FASTA file content
            
        Returns:
            List of (sequence, name) tuples
        """
        sequences = []
        
        current_name = ""
        current_seq = []
        
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                # New sequence header
                if current_name and current_seq:
                    # Save previous sequence
                    seq = ''.join(current_seq)
                    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())
                    if seq:
                        sequences.append((seq, current_name))
                
                current_name = line[1:].strip()
                current_seq = []
            elif line and not line.startswith(';'):
                # Sequence line (ignore comment lines starting with ;)
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_name and current_seq:
            seq = ''.join(current_seq)
            seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())
            if seq:
                sequences.append((seq, current_name))
        
        return sequences

    def validate_sequence(self, sequence: str) -> bool:
        """
        Validate that sequence contains only valid amino acids
        
        Args:
            sequence: Protein sequence to validate
            
        Returns:
            True if valid, False otherwise
        """
        if not sequence:
            return False
        
        # Check for valid amino acid characters
        valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
        seq_chars = set(sequence.upper())
        
        return seq_chars.issubset(valid_chars)

    def store_sequence_in_workspace(self, sequence: str, name: str, source: str):
        """
        Store a sequence in workspace for future use
        
        Args:
            sequence: Protein sequence
            name: Sequence name
            source: Source description
        """
        if not hasattr(self.module, 'update_workspace'):
            return
        
        # Get existing custom sequences
        custom_sequences = self.module.get_from_workspace("custom_sequences", {})
        
        # Add new sequence
        custom_sequences[name] = {
            'sequence': sequence,
            'source': source,
            'timestamp': self._get_timestamp()
        }
        
        # Update workspace
        self.module.update_workspace("custom_sequences", custom_sequences)

    def _get_timestamp(self) -> str:
        """Get current timestamp"""
        from datetime import datetime
        return datetime.now().isoformat()

    def get_sequence_statistics(self, sequence: str) -> Dict[str, Any]:
        """
        Get basic statistics for a sequence
        
        Args:
            sequence: Protein sequence
            
        Returns:
            Dictionary with sequence statistics
        """
        if not sequence:
            return {}
        
        # Count amino acids
        aa_counts = {}
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            aa_counts[aa] = sequence.count(aa)
        
        # Calculate basic properties
        length = len(sequence)
        composition = {aa: (count / length) * 100 for aa, count in aa_counts.items()}
        
        return {
            'length': length,
            'composition': composition,
            'aa_counts': aa_counts
        }
