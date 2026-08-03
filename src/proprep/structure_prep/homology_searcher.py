"""
BLAST Integration Module

Module for performing BLAST homology searches on PDB structure sequences.
Integrates with the MPSA processor workflow using the command pattern.
"""

import logging
import os
from typing import Any, Dict, List, Optional

from rich.panel import Panel
from rich.status import Status
from rich.table import Table

from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.application.processor_command import ModuleActionCommand

from .blast_worker import BLASTWorker
from .blast_commands import (
    RunBLASTCommand,
    ViewBLASTResultsCommand,
    ExportBLASTResultsCommand,
    BuildHomologyModelCommand,
)
from .homology_modeler import HomologyModeler, HAS_MODELLER


logger = logging.getLogger(__name__)


@register_module
class BLASTIntegrationModule(ProcessingModule):
    """Module for BLAST sequence analysis of PDB structures"""

    NAME = "Homology Searcher"
    DESCRIPTION = "Run BLAST homology searches"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    PRIORITY = 50

    def initialize(self):
        """Initialize module resources"""
        self.blast_worker = BLASTWorker()
        self.blast_results = None
        self.blast_result_file = None
        self.raw_results = None
        self._local_command_history = []

    @property
    def console(self):
        """Get console from processor or create new one"""
        if self.processor and hasattr(self.processor, 'console'):
            return self.processor.console
        else:
            from rich.console import Console
            return Console()

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for BLAST module data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        if workspace_key == "blast_results" and value is not None:
            panel_content = []

            panel_content.append(
                f"[cyan]Program:[/cyan] {value.get('program', 'Unknown')}"
            )
            panel_content.append(
                f"[cyan]Database:[/cyan] {value.get('database', 'Unknown')}"
            )
            panel_content.append(
                f"[cyan]Sequence Length:[/cyan] {value.get('sequence_length', 0)} residues"
            )
            panel_content.append(
                f"[cyan]Number of Hits:[/cyan] {value.get('hits_count', 0)}"
            )

            params = value.get("parameters", {})
            if params:
                panel_content.append("\n[cyan]Parameters:[/cyan]")
                panel_content.append(
                    f"  E-value threshold: {params.get('e_value', 'N/A')}"
                )
                panel_content.append(
                    f"  Word size: {params.get('word_size', 'Default')}"
                )
                panel_content.append(f"  Maximum hits: {params.get('max_hits', 'N/A')}")

                if "megablast" in params and params["megablast"] is not None:
                    panel_content.append(
                        f"  Megablast: {'Enabled' if params['megablast'] else 'Disabled'}"
                    )

            result_file = value.get("result_file")
            if result_file:
                panel_content.append(f"\n[cyan]Result File:[/cyan] {result_file}")

            console.print(
                Panel(
                    "\n".join(panel_content),
                    title="BLAST Results",
                    border_style="green",
                )
            )
            return True

        return False

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "run_blast": "Run BLAST homology search",
            "view_results": "View BLAST results",
            "retrieve_alphafold": "Retrieve AlphaFold structures for hits",
            "export_alignments": "Export alignments and annotations",
            "build_homology_model": "Build homology model (MODELLER)",
        }

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.

        Args:
            workspace: Current workspace

        Returns:
            List of MenuOption objects with status
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        # Check if BLAST has been run
        blast_results = workspace.get("blast_results")
        blast_run = blast_results is not None

        # Option 1: Run BLAST - needs a loaded structure (the query
        # sequence is read from it). Mirror option 5's can_process gate so
        # the indicator can't say ✓ while the action would fail.
        if blast_run:
            blast_status, blast_dep = OptionStatus.COMPLETED, ""
        elif self.can_process(workspace):
            blast_status, blast_dep = OptionStatus.AVAILABLE, ""
        else:
            blast_status = OptionStatus.BLOCKED
            blast_dep = self.availability_note(workspace) or "Load a structure first"
        options.append(MenuOption(
            key="1",
            description="Run BLAST homology search",
            status=blast_status,
            dependency_text=blast_dep,
        ))

        # Option 2: View results - requires BLAST to be run
        options.append(MenuOption(
            key="2",
            description="View BLAST results",
            status=OptionStatus.READY if blast_run else OptionStatus.BLOCKED,
            dependency_text="[Need to run BLAST first] ○" if not blast_run else ""
        ))

        # Option 3: Retrieve AlphaFold structures - requires BLAST to be run
        alphafold_homologs = workspace.get("alphafold_homologs")
        alphafold_retrieved = alphafold_homologs is not None and len(alphafold_homologs) > 0

        if alphafold_retrieved:
            af_status = OptionStatus.COMPLETED
            af_desc = f"Retrieve AlphaFold structures for hits [{len(alphafold_homologs)} retrieved] ✓"
        elif blast_run:
            af_status = OptionStatus.READY
            af_desc = "Retrieve AlphaFold structures for hits"
        else:
            af_status = OptionStatus.BLOCKED
            af_desc = "Retrieve AlphaFold structures for hits"

        options.append(MenuOption(
            key="3",
            description=af_desc,
            status=af_status,
            dependency_text="[Need to run BLAST first] ○" if not blast_run else ""
        ))

        # Option 4: Export alignments - requires BLAST to be run
        options.append(MenuOption(
            key="4",
            description="Export alignments and annotations",
            status=OptionStatus.READY if blast_run else OptionStatus.BLOCKED,
            dependency_text="[Need to run BLAST first] ○" if not blast_run else ""
        ))

        # Option 5: Build homology model - requires MODELLER and a loaded structure
        has_structure = self.can_process(workspace)
        homology_model_done = workspace.get("homology_model_pdb_file") is not None

        if not HAS_MODELLER:
            hm_status = OptionStatus.BLOCKED
            hm_dep = "[MODELLER not installed] ○"
        elif not has_structure:
            hm_status = OptionStatus.BLOCKED
            hm_dep = "[Need to load a structure first] ○"
        elif homology_model_done:
            hm_status = OptionStatus.COMPLETED
            hm_dep = ""
        else:
            hm_status = OptionStatus.READY
            hm_dep = ""

        options.append(MenuOption(
            key="5",
            description="Build homology model (MODELLER)",
            status=hm_status,
            dependency_text=hm_dep,
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
        blast_results = workspace.get("blast_results")

        if not blast_results:
            if not self.can_process(workspace):
                return f"{self.availability_note(workspace) or 'A structure is required'}. Load one via the Structure Loader."
            return "Start by running a BLAST search (option 1) to find homologous structures, or build a homology model (option 5)"
        else:
            return "View results with option 2, export with option 4, build a homology model with option 5, or press [m] to return to the main menu"

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "run_blast":
                command = RunBLASTCommand(self.processor, interactive=True)
                return command.execute_with_error_handling()

            elif option == "view_results":
                command = ViewBLASTResultsCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "retrieve_alphafold":
                self.retrieve_alphafold_for_hits()
                return True

            elif option == "export_alignments":
                command = ExportBLASTResultsCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "build_homology_model":
                command = BuildHomologyModelCommand(self.processor)
                return command.execute_with_error_handling()

            else:
                logger.warning(f"Unknown menu option: {option}")
                return False

        except Exception as e:
            logger.error(f"Error handling menu option '{option}': {str(e)}")
            return False

    def run_blast_search(
        self, workspace: Dict[str, Any] = None, interactive=True
    ) -> Dict[str, Any]:
        """Run BLAST homology search on PDB structure."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        # Use structure selector to get the structure file
        from proprep.utils.structure_selector import get_interactive_pdb_file

        pdb_file = get_interactive_pdb_file(
            workspace,
            self.console,
            processor=self.processor
        )

        if not pdb_file:
            self.processor.console.print(
                "[yellow]No PDB file loaded. Please load a PDB file first.[/yellow]"
            )
            return workspace

        structure_id = self._get_structure_id(workspace, pdb_file)
        chain_sequences = self.blast_worker.extract_sequences_from_pdb(pdb_file)

        if not chain_sequences:
            self.processor.console.print(
                "[yellow]No sequences found in PDB file.[/yellow]"
            )
            return workspace

        if interactive:
            success = self._run_interactive_blast(structure_id, chain_sequences)
        else:
            success = self._run_automated_blast(chain_sequences)

        if success:
            results_info = self._save_blast_results()
            if results_info:
                workspace = self.update_workspace(
                    workspace, "blast_results", results_info
                )
                workspace = self.update_workspace(
                    workspace, "blast_raw_results", self.blast_worker.results
                )
                elapsed = getattr(self.blast_worker, 'elapsed_time', None)
                time_str = f" in {elapsed:.0f}s" if elapsed else ""
                self.processor.console.print(
                    f"[green]BLAST search completed{time_str}. "
                    f"Found {results_info['hits_count']} hits.[/green]"
                )

                # Automatically display results
                self._display_blast_results()

        return workspace

    def view_blast_results(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """View BLAST results."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        raw_results = workspace.get("blast_raw_results")
        if not raw_results:
            self.processor.console.print(
                "[yellow]No BLAST results available. Run a BLAST search first.[/yellow]"
            )
            return workspace

        self.blast_worker.results = raw_results
        self._display_blast_results()
        return workspace

    def export_blast_results(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Export BLAST alignments and annotations."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        raw_results = workspace.get("blast_raw_results")
        if not raw_results:
            self.processor.console.print(
                "[yellow]No BLAST results available. Run a BLAST search first.[/yellow]"
            )
            return workspace

        output_file = prompt_with_context(
            processor=self.processor,
            prompt="Enter output filename for BLAST results",
            default="blast_export.xml",
            module="Homology Searcher",
            description="Enter output filename for BLAST results"
        )

        try:
            with open(output_file, "w") as f:
                f.write(
                    raw_results.decode("utf-8")
                    if isinstance(raw_results, bytes)
                    else raw_results
                )

            self.processor.console.print(
                f"[green]BLAST results exported to: {output_file}[/green]"
            )

            if confirm_with_context(
                processor=self.processor,
                prompt="Export summary table as well?",
                default=True,
                module="Homology Searcher",
                description="Export summary table"
            ):
                self._export_summary_table(output_file.replace(".xml", "_summary.tsv"))

        except Exception as e:
            self.processor.console.print(
                f"[red]Error exporting results: {str(e)}[/red]"
            )

        return workspace

    def _get_structure_id(self, workspace, pdb_file):
        """Extract structure identifier from metadata or filename.

        For RCSB/Local PDB: Returns PDB ID (e.g., '4LM8')
        For AlphaFold: Returns UniProt ID (e.g., 'P00520')
        """
        structure_id = None

        # Check for metadata in workspace (could be from any source)
        metadata = workspace.get("metadata")
        if (
            metadata
            and hasattr(metadata, "header_info")
            and "pdb_id" in metadata.header_info
        ):
            structure_id = metadata.header_info["pdb_id"]

        # Try to extract from filename
        if not structure_id and pdb_file:
            basename = os.path.basename(pdb_file)

            # Check if it's an AlphaFold structure (AF-UNIPROT-F1 pattern)
            if basename.startswith("AF-") and "-F" in basename:
                # Extract UniProt ID from AF-P00520-F1.cif
                parts = basename.split("-")
                if len(parts) >= 2:
                    structure_id = parts[1]  # The UniProt ID
                    self.processor.console.print(
                        f"[yellow]Using UniProt ID {structure_id} from AlphaFold structure[/yellow]"
                    )
            # Check if it's a PDB file (4-character alphanumeric)
            elif len(basename) >= 4 and basename[:4].isalnum():
                structure_id = basename[:4].upper()
                self.processor.console.print(
                    f"[yellow]Using PDB ID {structure_id} extracted from filename[/yellow]"
                )

        return structure_id

    def _run_interactive_blast(self, structure_id, chain_sequences):
        """Run interactive BLAST with user parameter selection."""
        self.processor.console.print(
            "\n[bold]===== BLAST Interactive Search Tool =====[/bold]"
        )
        self.processor.console.print("Enter parameters for your BLAST search")

        sequence = self._get_sequence_input(structure_id, chain_sequences)
        if not sequence:
            return False

        self.blast_worker.sequence = sequence

        self._get_blast_program()
        self._get_blast_database()
        self._get_blast_parameters()

        self._display_blast_parameters()

        if not confirm_with_context(
            processor=self.processor,
            prompt="Run BLAST with these parameters?",
            default=True,
            module="Homology Searcher",
            description="Confirm BLAST search parameters"
        ):
            return False

        with Status(
            "[bold green]Submitting BLAST search to NCBI server...[/bold green]"
        ):
            self.processor.console.print(
                "This may take several minutes depending on sequence length and server load."
            )
            return self.blast_worker.run_blast()

    def _run_automated_blast(self, chain_sequences):
        """Run automated BLAST with default parameters."""
        if not chain_sequences:
            return False

        first_chain = list(chain_sequences.keys())[0]
        self.blast_worker.sequence = chain_sequences[first_chain]
        self.blast_worker.set_program_defaults("blastp")

        return self.blast_worker.run_blast()

    def _get_sequence_input(self, structure_id, chain_sequences):
        """Get sequence input from user."""
        self.processor.console.print("\n[bold]Sequence input methods:[/bold]")
        self.processor.console.print("1. Enter sequence directly", highlight=False)
        self.processor.console.print("2. Load from FASTA file", highlight=False)

        options = ["1", "2"]
        if structure_id and chain_sequences:
            self.processor.console.print(f"3. Use sequence from loaded structure ({structure_id})")
            options.append("3")
        elif structure_id:
            # Only offer download option for PDB IDs (4-character codes)
            if len(structure_id) == 4 and structure_id.isalnum():
                self.processor.console.print(f"3. Download sequences for PDB {structure_id}")
                options.append("3")

        # Build options map dynamically
        options_map = {
            "1": "Enter sequence directly",
            "2": "Load from FASTA file"
        }
        if "3" in options:
            if chain_sequences:
                options_map["3"] = f"Use sequence from loaded structure ({structure_id})"
            else:
                options_map["3"] = f"Download sequences for PDB {structure_id}"

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Choose input method",
            choices=options,
            default="1",
            module="Homology Searcher",
            description="Select sequence input method",
            options_map=options_map
        )

        if choice == "1":
            return prompt_with_context(
                processor=self.processor,
                prompt="Enter your sequence",
                module="Homology Searcher",
                description="Enter protein/nucleotide sequence"
            )
        elif choice == "2":
            return self._load_sequence_from_file()
        elif choice == "3":
            return self._get_structure_sequence(structure_id, chain_sequences)

        return None

    def _load_sequence_from_file(self):
        """Load sequence from FASTA file."""
        while True:
            filename = prompt_with_context(
                processor=self.processor,
                prompt="Enter FASTA filename",
                module="Homology Searcher",
                description="Enter FASTA filename to load"
            )
            if os.path.isfile(filename):
                try:
                    with open(filename, "r") as f:
                        lines = f.readlines()
                        if lines[0].startswith(">"):
                            sequence = "".join(line.strip() for line in lines[1:])
                        else:
                            sequence = "".join(line.strip() for line in lines)

                    self.processor.console.print(f"Loaded {len(sequence)} characters")
                    return sequence
                except Exception as e:
                    self.processor.console.print(
                        f"[red]Error reading file: {str(e)}[/red]"
                    )
            else:
                self.processor.console.print("[red]File not found. Try again.[/red]")

    def _get_structure_sequence(self, structure_id, chain_sequences):
        """Get sequence from structure chains (PDB or AlphaFold)."""
        if structure_id and chain_sequences:
            # Sequences already extracted from loaded structure
            self.processor.console.print(
                f"\n[bold]Available chains in structure {structure_id}:[/bold]"
            )
            for chain_id, sequence in chain_sequences.items():
                self.processor.console.print(
                    f"Chain {chain_id}: {len(sequence)} residues"
                )

            chain_choices = list(chain_sequences.keys())
            options_map = {chain_id: f"Chain {chain_id}" for chain_id in chain_choices}
            selected_chain = prompt_with_context(
                processor=self.processor,
                prompt="Select chain",
                choices=chain_choices,
                default=chain_choices[0],
                module="Homology Searcher",
                description="Select chain from structure",
                options_map=options_map
            )

            if selected_chain in chain_sequences:
                sequence = chain_sequences[selected_chain]
                self.processor.console.print(
                    f"Using sequence from chain {selected_chain} ({len(sequence)} residues)"
                )
                # Halo the selected chain in the viewer so the user sees
                # exactly which chain is about to be BLAST'd. Stable label
                # so a re-run replaces rather than accumulates.
                try:
                    from proprep.structure_prep.viewer_coordinator import viewer as _viewer
                    _viewer.highlight(
                        f":{selected_chain}",
                        style="halo",
                        color="#1f78b4",
                        label="homology_blast_chain",
                    )
                except Exception:
                    pass
                return sequence
        elif structure_id and len(structure_id) == 4 and structure_id.isalnum():
            # Only download if it's a PDB ID (4-character code)
            with Status(
                f"[bold green]Downloading sequences for PDB {structure_id}...[/bold green]"
            ):
                downloaded_sequences = self.blast_worker.download_fasta(structure_id)

            if downloaded_sequences:
                return self._select_downloaded_sequence(structure_id, downloaded_sequences)

        return None

    def _select_downloaded_sequence(self, pdb_id, sequences):
        """Select from downloaded sequences."""
        self.processor.console.print(
            f"\n[bold]Available chains in PDB {pdb_id}:[/bold]"
        )
        for chain_id, sequence in sequences.items():
            self.processor.console.print(f"Chain {chain_id}: {len(sequence)} residues")

        chain_choices = list(sequences.keys())
        options_map = {chain_id: f"Chain {chain_id}" for chain_id in chain_choices}
        selected_chain = prompt_with_context(
            processor=self.processor,
            prompt="Select chain",
            choices=chain_choices,
            default=chain_choices[0],
            module="Homology Searcher",
            description="Select chain from downloaded sequences",
            options_map=options_map
        )

        if selected_chain in sequences:
            sequence = sequences[selected_chain]
            self.processor.console.print(
                f"Using sequence from chain {selected_chain} ({len(sequence)} residues)"
            )
            return sequence

        return None

    def _get_blast_program(self):
        """Get BLAST program selection."""
        self.processor.console.print("\n[bold]BLAST Program Options:[/bold]")
        self.processor.console.print("1. blastn - Nucleotide to nucleotide", highlight=False)
        self.processor.console.print("2. blastp - Protein to protein", highlight=False)
        self.processor.console.print("3. blastx - Nucleotide to protein", highlight=False)
        self.processor.console.print("4. tblastn - Protein to nucleotide", highlight=False)
        self.processor.console.print(
            "5. tblastx - Translated nucleotide to translated nucleotide"
        )

        program_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select BLAST program",
            choices=["1", "2", "3", "4", "5"],
            default="2",
            module="Homology Searcher",
            description="Select BLAST program type",
            options_map={
                "1": "blastn - Nucleotide to nucleotide",
                "2": "blastp - Protein to protein",
                "3": "blastx - Nucleotide to protein",
                "4": "tblastn - Protein to nucleotide",
                "5": "tblastx - Translated nucleotide to translated nucleotide"
            }
        )

        program_map = {
            "1": "blastn",
            "2": "blastp",
            "3": "blastx",
            "4": "tblastn",
            "5": "tblastx",
        }

        self.blast_worker.set_program_defaults(
            program_map.get(program_choice, "blastp")
        )

    def _get_blast_database(self):
        """Get database selection."""
        self.processor.console.print("\n[bold]Database Options:[/bold]")

        if self.blast_worker.program in ["blastn", "tblastn", "tblastx"]:
            self.processor.console.print("1. nt - Nucleotide collection", highlight=False)
            self.processor.console.print("2. refseq_rna - NCBI RNA reference sequences", highlight=False)
            self.processor.console.print(
                "3. refseq_genomic - NCBI genomic reference sequences"
            )
            self.processor.console.print("4. pdbnt - PDB nucleotide database", highlight=False)
            database_choices = ["1", "2", "3", "4"]
            options_map = {
                "1": "nt - Nucleotide collection",
                "2": "refseq_rna - NCBI RNA reference sequences",
                "3": "refseq_genomic - NCBI genomic reference sequences",
                "4": "pdbnt - PDB nucleotide database"
            }
        else:
            self.processor.console.print("1. nr - Non-redundant protein sequences [partial AlphaFold support]", highlight=False)
            self.processor.console.print(
                "2. refseq_protein - NCBI protein reference sequences"
            )
            self.processor.console.print("3. swissprot - UniProtKB/Swiss-Prot [best for AlphaFold retrieval]", highlight=False)
            self.processor.console.print("4. pdb - Protein Data Bank proteins", highlight=False)

            self.processor.console.print(
                "\n[grey50]Note: Options 1 and 3 allow retrieval of AlphaFold structures for homologs.[/grey50]"
            )

            database_choices = ["1", "2", "3", "4"]
            options_map = {
                "1": "nr - Non-redundant protein sequences [partial AlphaFold support]",
                "2": "refseq_protein - NCBI protein reference sequences",
                "3": "swissprot - UniProtKB/Swiss-Prot [best for AlphaFold retrieval]",
                "4": "pdb - Protein Data Bank proteins"
            }

        database_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select database",
            choices=database_choices,
            default="1",
            module="Homology Searcher",
            description="Select BLAST database",
            options_map=options_map
        )
        self.blast_worker.set_database_for_program(
            database_choice, self.blast_worker.program
        )

    def _get_blast_parameters(self):
        """Get additional BLAST parameters."""
        e_value_input = prompt_with_context(
            processor=self.processor,
            prompt="Enter E-value threshold",
            default=str(self.blast_worker.e_value),
            module="Homology Searcher",
            description="Enter E-value threshold for BLAST"
        )
        try:
            self.blast_worker.e_value = float(e_value_input)
        except ValueError:
            self.processor.console.print(
                "[yellow]Invalid E-value. Using default.[/yellow]"
            )

        max_hits_input = prompt_with_context(
            processor=self.processor,
            prompt="Enter maximum number of hits to return",
            default=str(self.blast_worker.max_hits),
            module="Homology Searcher",
            description="Enter maximum number of BLAST hits"
        )
        try:
            self.blast_worker.max_hits = int(max_hits_input)
        except ValueError:
            self.processor.console.print(
                "[yellow]Invalid number. Using default.[/yellow]"
            )

        if self.blast_worker.program in ["blastp", "blastx", "tblastn"]:
            self._get_protein_parameters()

        if self.blast_worker.program == "blastn":
            self.blast_worker.megablast = confirm_with_context(
                processor=self.processor,
                prompt="Use MegaBLAST (optimized for highly similar sequences)?",
                default=True,
                module="Homology Searcher",
                description="Use MegaBLAST for nucleotide search"
            )

    def _get_protein_parameters(self):
        """Get protein-specific parameters."""
        self.processor.console.print("\n[bold]Scoring Matrix Options:[/bold]")
        self.processor.console.print(
            "[grey50]Lower numbers = more divergent sequences, Higher numbers = more similar sequences[/grey50]\n"
        )
        self.processor.console.print(
            "1. BLOSUM62 - General purpose (~62% identity)"
        )
        self.processor.console.print(
            "   [grey50]Best for: Most searches, remote homologs, default choice[/grey50]"
        )
        self.processor.console.print(
            "2. BLOSUM45 - Divergent sequences (~45% identity)"
        )
        self.processor.console.print(
            "   [grey50]Best for: Distant evolutionary relationships, deep homology[/grey50]"
        )
        self.processor.console.print(
            "3. BLOSUM80 - Closely related (~80% identity)"
        )
        self.processor.console.print(
            "   [grey50]Best for: Recent duplications, orthologs within same genus[/grey50]"
        )
        self.processor.console.print(
            "4. BLOSUM90 - Very similar (~90% identity)"
        )
        self.processor.console.print(
            "   [grey50]Best for: Variants, isoforms, species-specific differences[/grey50]"
        )

        matrix_choice = prompt_with_context(
            processor=self.processor,
            prompt="Select scoring matrix",
            choices=["1", "2", "3", "4"],
            default="1",
            module="Homology Searcher",
            description="Select scoring matrix for protein BLAST",
            options_map={
                "1": "BLOSUM62 - General purpose (~62% identity)",
                "2": "BLOSUM45 - Divergent sequences (~45% identity)",
                "3": "BLOSUM80 - Closely related (~80% identity)",
                "4": "BLOSUM90 - Very similar (~90% identity)"
            }
        )
        self.blast_worker.set_matrix_and_gap_costs(matrix_choice)

    def _display_blast_parameters(self):
        """Display current BLAST parameters."""
        self.processor.console.print(
            "\n[bold]Ready to run BLAST with the following parameters:[/bold]"
        )
        self.processor.console.print(
            f"  [cyan]Program:[/cyan] {self.blast_worker.program}"
        )
        self.processor.console.print(
            f"  [cyan]Database:[/cyan] {self.blast_worker.database}"
        )
        self.processor.console.print(
            f"  [cyan]Sequence length:[/cyan] {len(self.blast_worker.sequence)} characters"
        )
        self.processor.console.print(
            f"  [cyan]E-value threshold:[/cyan] {self.blast_worker.e_value}"
        )
        self.processor.console.print(
            f"  [cyan]Maximum hits:[/cyan] {self.blast_worker.max_hits}"
        )

        if self.blast_worker.program == "blastn":
            self.processor.console.print(
                f"  [cyan]Megablast:[/cyan] {'enabled' if self.blast_worker.megablast else 'disabled'}"
            )

    def _save_blast_results(self):
        """Save BLAST results with auto-generated filename and return summary."""
        if not self.blast_worker.results:
            return None

        import time as _time
        filename = f"blast_results_{self.blast_worker.program}.xml"
        if os.path.exists(filename):
            filename = f"blast_results_{self.blast_worker.program}_{int(_time.time())}.xml"

        result_file = self.blast_worker.save_results(filename)
        if not result_file:
            return None

        self.processor.console.print(
            f"[grey50]Results saved to: {result_file}[/grey50]"
        )

        summary = self.blast_worker.get_results_summary()
        if summary:
            results_info = {
                "program": self.blast_worker.program,
                "database": self.blast_worker.database,
                "sequence_length": len(self.blast_worker.sequence),
                "hits_count": summary["num_hits"],
                "result_file": result_file,
                "parameters": {
                    "e_value": self.blast_worker.e_value,
                    "word_size": self.blast_worker.word_size,
                    "max_hits": self.blast_worker.max_hits,
                    "megablast": (
                        self.blast_worker.megablast
                        if self.blast_worker.program == "blastn"
                        else None
                    ),
                },
            }
            return results_info

        return None

    def _display_blast_results(self):
        """Display BLAST results in formatted tables."""
        summary = self.blast_worker.get_results_summary()
        if not summary:
            self.processor.console.print("[yellow]No results to display.[/yellow]")
            return

        self.processor.console.print("\n[bold]===== BLAST Results Summary =====[/bold]")
        self.processor.console.print(
            f"[bold]Query:[/bold] {summary['query'][:50]}{'...' if len(summary['query']) > 50 else ''}"
        )
        self.processor.console.print(f"[bold]Database:[/bold] {summary['database']}")
        self.processor.console.print(
            f"[bold]Number of hits:[/bold] {summary['num_hits']}"
        )

        if not summary["alignments"]:
            self.processor.console.print(
                "[yellow]No significant alignments found.[/yellow]"
            )
            return

        self.processor.console.print("\n[bold]Top Hits:[/bold]")
        table = Table(show_header=True, header_style="bold", expand=True)
        table.add_column("Hit", style="cyan", width=4)
        table.add_column("Description", style="green", no_wrap=False, min_width=30, max_width=80)
        table.add_column("Score", style="white", width=8)
        table.add_column("E-value", style="yellow", width=10)
        table.add_column("Identity", style="magenta", width=10)
        table.add_column("Length", style="blue", width=8)

        for i, alignment in enumerate(summary["alignments"], 1):
            if i > self.blast_worker.max_hits:
                break

            title = alignment["title"]
            description = title[:77] + "..." if len(title) > 80 else title

            table.add_row(
                f"{i}",
                description,
                f"{alignment['bits']:.0f}",
                f"{alignment['e_value']:.2e}",
                f"{alignment['identity_percent']:.1f}%",
                f"{alignment['align_length']}",
            )

        self.processor.console.print(table)

        if confirm_with_context(
            processor=self.processor,
            prompt="\nView detailed alignment for specific hits?",
            default=True,
            module="Homology Searcher",
            description="View detailed alignment for hits"
        ):
            hit_nums_input = prompt_with_context(
                processor=self.processor,
                prompt="Enter hit numbers (comma-separated)",
                default="1",
                module="Homology Searcher",
                description="Enter hit numbers to view details"
            )
            try:
                hit_nums = [int(num.strip()) for num in hit_nums_input.split(",")]
                self._display_detailed_alignments(summary["alignments"], hit_nums)
            except ValueError:
                self.processor.console.print("[yellow]Invalid input format.[/yellow]")

    def _display_detailed_alignments(self, alignments, hit_nums):
        """Display detailed alignment information for selected hits."""
        for hit_num in hit_nums:
            if 1 <= hit_num <= len(alignments):
                alignment = alignments[hit_num - 1]

                self.processor.console.print(
                    f"\n[bold]Detailed view of hit #{hit_num}:[/bold]"
                )
                self.processor.console.print(
                    f"[bold]Description:[/bold] {alignment['title']}"
                )
                self.processor.console.print(
                    f"[bold]Score:[/bold] {alignment['score']} bits ({alignment['bits']})"
                )
                self.processor.console.print(
                    f"[bold]E-value:[/bold] {alignment['e_value']}"
                )
                self.processor.console.print(
                    f"[bold]Identities:[/bold] {alignment['identities']}/{alignment['align_length']} ({alignment['identity_percent']:.1f}%)"
                )

                if alignment.get("positives") is not None:
                    self.processor.console.print(
                        f"[bold]Positives:[/bold] {alignment['positives']}/{alignment['align_length']} "
                        f"({alignment.get('positives_percent', 0):.1f}%)"
                    )

                gaps = alignment.get("gaps", 0)
                if gaps:
                    gap_pct = (gaps / alignment['align_length']) * 100
                    self.processor.console.print(
                        f"[bold]Gaps:[/bold] {gaps}/{alignment['align_length']} ({gap_pct:.1f}%)"
                    )

                self.processor.console.print(
                    f"[bold]Query range:[/bold] {alignment['query_start']}-{alignment['query_end']}"
                )
                self.processor.console.print(
                    f"[bold]Subject range:[/bold] {alignment['sbjct_start']}-{alignment['sbjct_end']}"
                )

                self._format_alignment_block(alignment)
            else:
                self.processor.console.print(
                    f"[yellow]Hit #{hit_num} not found.[/yellow]"
                )

    def _format_alignment_block(self, alignment, wrap_width=60):
        """Format and display alignment as classic BLAST-style wrapped blocks."""
        query_seq = alignment.get("query_seq", "")
        match_seq = alignment.get("match_seq", "")
        sbjct_seq = alignment.get("sbjct_seq", "")

        if not query_seq:
            return

        query_pos = alignment["query_start"]
        sbjct_pos = alignment["sbjct_start"]

        self.processor.console.print("\n[bold]Alignment:[/bold]")

        for i in range(0, len(query_seq), wrap_width):
            q_chunk = query_seq[i:i + wrap_width]
            m_chunk = match_seq[i:i + wrap_width]
            s_chunk = sbjct_seq[i:i + wrap_width]

            q_end = query_pos + len(q_chunk.replace("-", "")) - 1
            s_end = sbjct_pos + len(s_chunk.replace("-", "")) - 1

            q_label = f"Query  {query_pos:>4d}  "
            s_label = f"Sbjct  {sbjct_pos:>4d}  "
            m_label = " " * len(q_label)

            self.processor.console.print(
                f"[cyan]{q_label}[/cyan]{q_chunk}  [cyan]{q_end}[/cyan]"
            )
            self.processor.console.print(f"{m_label}{m_chunk}")
            self.processor.console.print(
                f"[green]{s_label}[/green]{s_chunk}  [green]{s_end}[/green]"
            )
            self.processor.console.print("")

            query_pos = q_end + 1
            sbjct_pos = s_end + 1

    def _export_summary_table(self, filename):
        """Export summary table to TSV format."""
        summary = self.blast_worker.get_results_summary()
        if not summary:
            return

        try:
            with open(filename, "w") as f:
                f.write(
                    "Hit\tDescription\tScore\tE-value\tIdentity_%\tPositives_%\tGaps\tLength\t"
                    "Query_Start\tQuery_End\tSubject_Start\tSubject_End\n"
                )

                for i, alignment in enumerate(summary["alignments"], 1):
                    positives_pct = alignment.get('positives_percent', '')
                    if positives_pct != '':
                        positives_pct = f"{positives_pct:.1f}"
                    gaps = alignment.get('gaps', 0)

                    f.write(
                        f"{i}\t{alignment['title']}\t{alignment['bits']:.0f}\t"
                        f"{alignment['e_value']}\t{alignment['identity_percent']:.1f}\t"
                        f"{positives_pct}\t{gaps}\t"
                        f"{alignment['align_length']}\t{alignment['query_start']}\t{alignment['query_end']}\t"
                        f"{alignment['sbjct_start']}\t{alignment['sbjct_end']}\n"
                    )

            self.processor.console.print(
                f"[green]Summary table exported to: {filename}[/green]"
            )
        except Exception as e:
            self.processor.console.print(
                f"[red]Error exporting summary table: {str(e)}[/red]"
            )

    def retrieve_alphafold_for_hits(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Retrieve AlphaFold structures for BLAST hits with UniProt IDs."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        # Extract hits with UniProt IDs
        uniprot_hits = self._extract_uniprot_ids_from_blast(workspace)

        if not uniprot_hits:
            blast_results = workspace.get("blast_results", {})
            database = blast_results.get("database", "unknown")

            self.processor.console.print(
                "[yellow]No UniProt IDs found in BLAST results.[/yellow]"
            )
            self.processor.console.print(f"Current database: {database}")
            self.processor.console.print(
                "\n[grey50]AlphaFold retrieval requires UniProt accessions, "
                "which are available when searching:[/grey50]"
            )
            self.processor.console.print("  • SwissProt (swissprot)")
            self.processor.console.print("  • Non-redundant (nr) - partial coverage")
            self.processor.console.print(
                "\nConsider re-running BLAST against SwissProt for best AlphaFold coverage."
            )
            return workspace

        # Show hits with UniProt IDs
        self.processor.console.print(
            f"\n[green]Found {len(uniprot_hits)} BLAST hits with UniProt IDs[/green]\n"
        )

        # Let user select which hits to retrieve
        selected_hits = self._select_hits_for_retrieval(uniprot_hits)

        if not selected_hits:
            self.processor.console.print("[yellow]No hits selected for retrieval.[/yellow]")
            return workspace

        # Create homologs subdirectory
        import os
        from pathlib import Path
        homologs_dir = Path("homologs")
        homologs_dir.mkdir(exist_ok=True)
        self.processor.console.print(f"[grey50]Saving structures to: {homologs_dir.absolute()}[/grey50]\n")

        # Retrieve AlphaFold structures
        alphafold_homologs = {}
        failed_retrievals = []

        from proprep.structure_prep.alphafold_fetcher import AlphaFoldFetcher
        from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn

        # Suppress AlphaFold fetcher logging
        fetcher_logger = logging.getLogger('proprep.structure_prep.alphafold_fetcher')
        original_level = fetcher_logger.level
        fetcher_logger.setLevel(logging.WARNING)

        try:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                console=self.processor.console
            ) as progress:
                task = progress.add_task(
                    "[cyan]Retrieving AlphaFold structures...",
                    total=len(selected_hits)
                )

                for hit_info in selected_hits:
                    uniprot_id = hit_info["uniprot_id"]
                    progress.update(
                        task,
                        description=f"[cyan]Retrieving {uniprot_id}..."
                    )

                    try:
                        # Use AlphaFoldFetcher to download to homologs directory
                        fetcher = AlphaFoldFetcher(output_dir=str(homologs_dir))
                        result = fetcher.retrieve_structure(uniprot_id)

                        if result:
                            # Store with absolute path
                            alphafold_homologs[f"hit_{hit_info['rank']}"] = {
                                "uniprot_id": uniprot_id,
                                "blast_rank": hit_info["rank"],
                                "e_value": hit_info["e_value"],
                                "identity_percent": hit_info["identity_percent"],
                                "description": hit_info["description"],
                                "pdb_file": str(Path(result["pdb_file"]).absolute()),
                                "structure": result["structure"],
                                "confidence": result.get("confidence", {}),
                                "alphafold_id": result.get("alphafold_id", f"AF-{uniprot_id}-F1")
                            }
                        else:
                            failed_retrievals.append(uniprot_id)

                    except Exception as e:
                        logger.warning(f"Failed to retrieve {uniprot_id}: {e}")
                        failed_retrievals.append(uniprot_id)

                    progress.advance(task)

        finally:
            fetcher_logger.setLevel(original_level)

        # Report results
        self.processor.console.print(
            f"\n[bold green]✓ Successfully retrieved {len(alphafold_homologs)} AlphaFold structures[/bold green]"
        )

        if failed_retrievals:
            self.processor.console.print(
                f"[yellow]Failed to retrieve {len(failed_retrievals)} structures:[/yellow]"
            )
            for uniprot_id in failed_retrievals[:5]:  # Show first 5
                self.processor.console.print(f"  • {uniprot_id}")
            if len(failed_retrievals) > 5:
                self.processor.console.print(f"  ... and {len(failed_retrievals) - 5} more")

        # Prompt user to select a homolog for downstream processing
        if alphafold_homologs:
            selected_homolog = self._select_homolog_for_processing(alphafold_homologs)
            if selected_homolog:
                # Save selected homolog to dedicated workspace keys for downstream processing
                workspace = self.update_workspace(
                    workspace,
                    "alphafold_homolog_pdb_file",
                    selected_homolog["pdb_file"]
                )
                workspace = self.update_workspace(
                    workspace,
                    "alphafold_homolog_structure",
                    selected_homolog["structure"]
                )
                self.processor.console.print(
                    f"[green]✓ Selected {selected_homolog['uniprot_id']} for downstream processing[/green]"
                )

        # Save to workspace
        workspace = self.update_workspace(workspace, "alphafold_homologs", alphafold_homologs)

        return workspace

    def _extract_uniprot_ids_from_blast(self, workspace: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Extract UniProt IDs from BLAST results.

        Returns:
            List of dicts with keys: rank, uniprot_id, description, e_value, identity_percent
        """
        raw_results = workspace.get("blast_raw_results")
        if not raw_results:
            return []

        self.blast_worker.results = raw_results
        summary = self.blast_worker.get_results_summary()

        if not summary or not summary.get("alignments"):
            return []

        uniprot_hits = []

        for rank, alignment in enumerate(summary["alignments"], 1):
            title = alignment.get("title", "")
            uniprot_id = self._extract_uniprot_from_title(title)

            if uniprot_id:
                uniprot_hits.append({
                    "rank": rank,
                    "uniprot_id": uniprot_id,
                    "description": title,
                    "e_value": alignment.get("e_value", 0),
                    "identity_percent": alignment.get("identity_percent", 0)
                })

        return uniprot_hits

    def _extract_uniprot_from_title(self, title: str) -> Optional[str]:
        """
        Extract UniProt ID from BLAST hit title.

        Handles formats:
        - sp|P12345|GENE_ORG ... (SwissProt)
        - sp|P12345.1|GENE_ORG ... (SwissProt with version)
        - tr|A0A123|GENE_ORG ... (TrEMBL)

        Returns:
            UniProt accession (without version number) or None
        """
        if '|' not in title:
            return None

        parts = title.split('|')
        if len(parts) >= 2:
            prefix = parts[0].strip()
            if prefix in ['sp', 'tr']:  # SwissProt or TrEMBL
                uniprot_id = parts[1].strip()
                # Remove version number if present (e.g., P12345.1 -> P12345)
                if '.' in uniprot_id:
                    uniprot_id = uniprot_id.split('.')[0]
                return uniprot_id

        return None

    def _select_hits_for_retrieval(self, uniprot_hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Let user select which hits to retrieve.

        Returns:
            List of selected hit dicts
        """
        # Display table of available hits
        self.processor.console.print("[bold]BLAST Hits with UniProt IDs:[/bold]\n")

        table = Table(show_header=True, header_style="bold")
        table.add_column("Rank", style="cyan", width=6)
        table.add_column("UniProt ID", style="green", width=12)
        table.add_column("E-value", style="yellow", width=12)
        table.add_column("Identity", style="magenta", width=10)
        table.add_column("Description", style="white", no_wrap=False, min_width=30, max_width=80)

        for hit in uniprot_hits[:20]:  # Show top 20
            desc = hit["description"][:77] + "..." if len(hit["description"]) > 80 else hit["description"]
            table.add_row(
                str(hit["rank"]),
                hit["uniprot_id"],
                f"{hit['e_value']:.2e}",
                f"{hit['identity_percent']:.1f}%",
                desc
            )

        self.processor.console.print(table)

        if len(uniprot_hits) > 20:
            self.processor.console.print(
                f"\n[grey50]Showing top 20 of {len(uniprot_hits)} hits with UniProt IDs[/grey50]"
            )

        # Selection method
        self.processor.console.print("\n[bold]Selection Methods:[/bold]")
        self.processor.console.print("1. Top N hits", highlight=False)
        self.processor.console.print("2. By E-value threshold", highlight=False)
        self.processor.console.print("3. By identity % threshold", highlight=False)
        self.processor.console.print("4. Manual selection (enter ranks)", highlight=False)
        self.processor.console.print("5. All hits", highlight=False)

        choice = prompt_with_context(
            processor=self.processor,
            prompt="Select method",
            choices=["1", "2", "3", "4", "5"],
            default="1",
            module="Homology Searcher - AlphaFold",
            description="Select hit selection method",
            options_map={
                "1": "Top N hits",
                "2": "By E-value threshold",
                "3": "By identity % threshold",
                "4": "Manual selection",
                "5": "All hits"
            }
        )

        if choice == "1":
            # Top N
            max_n = min(len(uniprot_hits), 20)
            n_str = prompt_with_context(
                processor=self.processor,
                prompt=f"Enter number of top hits to retrieve (max {max_n})",
                default="5",
                module="Homology Searcher - AlphaFold",
                description="Enter number of top hits"
            )
            try:
                n = int(n_str)
                return uniprot_hits[:min(n, max_n)]
            except ValueError:
                self.processor.console.print("[yellow]Invalid number, using top 5[/yellow]")
                return uniprot_hits[:5]

        elif choice == "2":
            # E-value threshold
            threshold_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter E-value threshold (e.g., 1e-10)",
                default="1e-10",
                module="Homology Searcher - AlphaFold",
                description="Enter E-value threshold"
            )
            try:
                threshold = float(threshold_str)
                selected = [hit for hit in uniprot_hits if hit["e_value"] <= threshold]
                self.processor.console.print(f"[green]{len(selected)} hits meet threshold[/green]")
                return selected
            except ValueError:
                self.processor.console.print("[yellow]Invalid threshold[/yellow]")
                return []

        elif choice == "3":
            # Identity % threshold
            threshold_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter identity % threshold (e.g., 50)",
                default="50",
                module="Homology Searcher - AlphaFold",
                description="Enter identity % threshold"
            )
            try:
                threshold = float(threshold_str)
                selected = [hit for hit in uniprot_hits if hit["identity_percent"] >= threshold]
                self.processor.console.print(f"[green]{len(selected)} hits meet threshold[/green]")
                return selected
            except ValueError:
                self.processor.console.print("[yellow]Invalid threshold[/yellow]")
                return []

        elif choice == "4":
            # Manual selection
            ranks_str = prompt_with_context(
                processor=self.processor,
                prompt="Enter rank numbers (comma-separated, e.g., 1,3,5)",
                default="1",
                module="Homology Searcher - AlphaFold",
                description="Enter rank numbers to retrieve"
            )
            try:
                ranks = [int(r.strip()) for r in ranks_str.split(",")]
                selected = [hit for hit in uniprot_hits if hit["rank"] in ranks]
                self.processor.console.print(f"[green]Selected {len(selected)} hits[/green]")
                return selected
            except ValueError:
                self.processor.console.print("[yellow]Invalid input[/yellow]")
                return []

        elif choice == "5":
            # All hits
            if len(uniprot_hits) > 50:
                if not confirm_with_context(
                    processor=self.processor,
                    prompt=f"Retrieve all {len(uniprot_hits)} hits? This may take a while.",
                    default=False,
                    module="Homology Searcher - AlphaFold",
                    description="Confirm retrieving all hits"
                ):
                    return []
            return uniprot_hits

        return []

    def _select_homolog_for_processing(self, alphafold_homologs: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Prompt user to select one homolog for downstream processing.

        Args:
            alphafold_homologs: Dict of retrieved homologs from AlphaFold

        Returns:
            Selected homolog dict or None if user skips
        """
        # Ask if user wants to select a homolog
        self.processor.console.print("\n[bold cyan]Select Homolog for Downstream Processing[/bold cyan]")
        self.processor.console.print(
            "[grey50]You can select one homolog to use in downstream modules (structure preparation, MD setup, etc.)[/grey50]\n"
        )

        if not confirm_with_context(
            processor=self.processor,
            prompt="Would you like to select a homolog for downstream processing?",
            default=True,
            module="Homology Searcher - Select Homolog",
            description="Choose whether to select homolog for processing"
        ):
            return None

        # Build sorted list of homologs by BLAST rank
        homolog_list = []
        for key, homolog in alphafold_homologs.items():
            homolog_list.append(homolog)

        homolog_list.sort(key=lambda x: x["blast_rank"])

        # Display table of available homologs
        self.processor.console.print("\n[bold]Retrieved AlphaFold Homologs:[/bold]\n")

        table = Table(show_header=True, header_style="bold")
        table.add_column("Option", style="cyan", width=8)
        table.add_column("BLAST Rank", style="yellow", width=12)
        table.add_column("UniProt ID", style="green", width=12)
        table.add_column("E-value", style="yellow", width=12)
        table.add_column("Identity", style="magenta", width=10)
        table.add_column("Description", style="white", no_wrap=False, min_width=30, max_width=80)

        options_map = {}
        for idx, homolog in enumerate(homolog_list, 1):
            desc = homolog["description"][:77] + "..." if len(homolog["description"]) > 80 else homolog["description"]
            table.add_row(
                str(idx),
                str(homolog["blast_rank"]),
                homolog["uniprot_id"],
                f"{homolog['e_value']:.2e}",
                f"{homolog['identity_percent']:.1f}%",
                desc
            )
            options_map[str(idx)] = f"{homolog['uniprot_id']} (rank {homolog['blast_rank']})"

        self.processor.console.print(table)

        # Prompt for selection
        choice = prompt_with_context(
            processor=self.processor,
            prompt="\nSelect homolog (enter number)",
            choices=[str(i) for i in range(1, len(homolog_list) + 1)],
            default="1",
            module="Homology Searcher - Select Homolog",
            description="Select homolog for downstream processing",
            options_map=options_map
        )

        selected_idx = int(choice) - 1
        return homolog_list[selected_idx]

    def build_homology_model(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Build a homology model using MODELLER.

        Orchestrates: template selection -> chain extraction -> target sequence input
        -> alignment -> PIR file creation -> MODELLER execution -> quality assessment.
        """
        if workspace is None:
            workspace = self.processor._get_workspace()

        if not HAS_MODELLER:
            self.console.print(
                "[red]MODELLER is not installed. Cannot build homology model.[/red]\n"
                "[grey50]Install with: conda install -c salilab modeller[/grey50]"
            )
            return workspace

        # 1. Get template PDB
        from proprep.utils.structure_selector import get_interactive_pdb_file

        self.console.print("\n[bold]===== Build Homology Model (MODELLER) =====[/bold]")
        self.console.print("[grey50]Select the template structure to model from.[/grey50]\n")

        pdb_file = get_interactive_pdb_file(
            workspace, self.console, processor=self.processor
        )
        if not pdb_file:
            self.console.print("[yellow]No PDB file selected. Aborting.[/yellow]")
            return workspace

        # 2. Extract chain sequences and HETATM counts
        modeler = HomologyModeler()
        chain_sequences = modeler.extract_template_sequences(pdb_file)
        if not chain_sequences:
            self.console.print("[yellow]No protein sequences found in template.[/yellow]")
            return workspace

        hetatm_counts = modeler.count_hetatm_per_chain(pdb_file)

        # Display available chains
        self.console.print("\n[bold]Template chains:[/bold]")
        chain_table = Table(show_header=True, header_style="bold")
        chain_table.add_column("Chain", style="cyan", width=8)
        chain_table.add_column("Residues", style="yellow", width=10)
        chain_table.add_column("HETATMs", style="green", width=10)
        chain_table.add_column("Waters", style="blue", width=10)

        for cid in sorted(chain_sequences.keys()):
            counts = hetatm_counts.get(cid, {"hetatm": 0, "water": 0})
            chain_table.add_row(
                cid,
                str(len(chain_sequences[cid])),
                str(counts["hetatm"]),
                str(counts["water"]),
            )
        self.console.print(chain_table)

        # 3. Let user select chains
        all_chain_ids = sorted(chain_sequences.keys())
        if len(all_chain_ids) == 1:
            selected_chains = all_chain_ids
            self.console.print(f"\nUsing chain {selected_chains[0]}.")
        else:
            chains_input = prompt_with_context(
                processor=self.processor,
                prompt="Select chains to model (comma-separated, or 'all')",
                default="all",
                module="Homology Searcher",
                description="Select template chains to include in the model",
            )
            if chains_input.strip().lower() == "all":
                selected_chains = all_chain_ids
            else:
                selected_chains = [
                    c.strip().upper() for c in chains_input.split(",")
                    if c.strip().upper() in chain_sequences
                ]
            if not selected_chains:
                self.console.print("[yellow]No valid chains selected. Aborting.[/yellow]")
                return workspace

        # 4. Get target sequence for each chain
        from proprep.emboss.sequence_input import SequenceInputManager

        seq_manager = SequenceInputManager(self)
        chain_alignments = {}

        for chain_id in selected_chains:
            template_seq = chain_sequences[chain_id]
            self.console.print(
                f"\n[bold cyan]Target sequence for chain {chain_id}[/bold cyan] "
                f"(template: {len(template_seq)} residues)"
            )

            result = seq_manager.get_sequence_input(
                label=f"Target sequence for chain {chain_id}"
            )
            if result is None:
                self.console.print(f"[yellow]No target sequence provided for chain {chain_id}. Aborting.[/yellow]")
                return workspace

            target_seq, seq_name, seq_source = result
            self.console.print(
                f"[grey50]Target: {seq_name} ({len(target_seq)} residues, source: {seq_source})[/grey50]"
            )

            # 5. Align
            aln_result = modeler.align_sequences(template_seq, target_seq)
            modeler.display_alignment(aln_result, self.console, chain_id=chain_id)

            # Warn on low identity
            if aln_result["identity_pct"] < 25:
                self.console.print(
                    f"[red]Warning: Identity {aln_result['identity_pct']:.1f}% is below 25%. "
                    f"Model quality may be unreliable.[/red]"
                )
            elif aln_result["identity_pct"] < 30:
                self.console.print(
                    f"[yellow]Note: Identity {aln_result['identity_pct']:.1f}% is below 30%. "
                    f"Model should be used with caution.[/yellow]"
                )

            chain_alignments[chain_id] = {
                "template_aligned": aln_result["template_aligned"],
                "target_aligned": aln_result["target_aligned"],
                "identity_pct": aln_result["identity_pct"],
                "gaps_pct": aln_result["gaps_pct"],
                "score": aln_result["score"],
            }

        # 6. Confirm before running MODELLER
        self.console.print("\n[bold]Summary:[/bold]")
        for chain_id, aln in chain_alignments.items():
            self.console.print(
                f"  Chain {chain_id}: identity={aln['identity_pct']:.1f}%, "
                f"gaps={aln['gaps_pct']:.1f}%"
            )

        if not confirm_with_context(
            processor=self.processor,
            prompt="Proceed with MODELLER to build homology model?",
            default=True,
            module="Homology Searcher",
            description="Confirm homology model building",
        ):
            self.console.print("[yellow]Cancelled.[/yellow]")
            return workspace

        # 7. Create temp directory, copy template PDB, write PIR, run MODELLER
        import tempfile
        import shutil

        work_dir = tempfile.mkdtemp(prefix="proprep_homology_")
        template_basename = "input.pdb"
        shutil.copy2(pdb_file, os.path.join(work_dir, template_basename))

        alignment_file = os.path.join(work_dir, "alignment.ali")
        modeler.create_pir_file(
            output_file=alignment_file,
            template_pdb=template_basename,
            chain_ids=selected_chains,
            chain_alignments=chain_alignments,
            hetatm_counts=hetatm_counts,
        )

        self.console.print("\n[bold green]Running MODELLER...[/bold green]")
        success, result_path, model_structure = modeler.run_modeller(
            work_dir=work_dir,
            template_pdb=os.path.join(work_dir, template_basename),
            alignment_file=alignment_file,
            console=self.console,
        )

        if not success:
            self.console.print(f"[red]MODELLER failed: {result_path}[/red]")
            return workspace

        # 8. Copy output PDB to user-chosen location
        output_filename = prompt_with_context(
            processor=self.processor,
            prompt="Save homology model as",
            default="homology_model.pdb",
            module="Homology Searcher",
            description="Enter filename for the homology model PDB",
        )

        output_path = os.path.abspath(output_filename)
        shutil.copy2(result_path, output_path)
        self.console.print(f"[green]Homology model saved to: {output_path}[/green]")

        # 9. Store in workspace
        alignment_info = {
            "chains": {},
            "template_pdb": pdb_file,
            "selected_chains": selected_chains,
        }
        for chain_id, aln in chain_alignments.items():
            alignment_info["chains"][chain_id] = {
                "identity_pct": aln["identity_pct"],
                "gaps_pct": aln["gaps_pct"],
                "score": aln["score"],
            }

        workspace = self.update_workspace(workspace, "homology_model_pdb_file", output_path)
        workspace = self.update_workspace(workspace, "homology_model_structure", model_structure)
        workspace = self.update_workspace(workspace, "homology_model_alignment", alignment_info)

        self.console.print(
            "[green]Homology model results saved to workspace.[/green]"
        )

        # Clean up temp directory
        try:
            shutil.rmtree(work_dir)
        except OSError:
            pass

        return workspace

    def _record_local_command(self, command):
        """Record command in local history for module-specific tracking"""
        self._local_command_history.append(command)

    def get_breadcrumb_for_command(self, command):
        """Return breadcrumb string for a command"""
        action_to_run = getattr(command, "_action_to_run", None)
        if action_to_run:
            action_map = {
                "run_blast_search": "Run BLAST",
                "view_blast_results": "View Results",
                "export_blast_results": "Export",
                "build_homology_model": "Build Homology Model",
            }
            action_name = action_map.get(action_to_run, action_to_run)
            return f"{self.NAME} > {action_name}"
        return self.NAME

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements - need at least one structure loaded"""
        return ["rcsb_pdb_file | local_pdb_file | alphafold_pdb_file | alphafill_pdb_file"]

    def get_workspace_outputs(self) -> List[str]:
        """Get workspace outputs"""
        return [
            "blast_results",
            "blast_raw_results",
            "alphafold_homologs",
            "alphafold_homolog_pdb_file",
            "alphafold_homolog_structure",
            "homology_model_pdb_file",
            "homology_model_structure",
            "homology_model_alignment",
        ]

    def can_process(self, workspace: Dict[str, Any]) -> bool:
        """Check if the module can process the current workspace.

        Uses StructureSelector to check for any available structure.
        """
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def process(self, workspace):
        """Process the workspace"""
        if self.can_process(workspace):
            blast_results = workspace.get("blast_results")
            if blast_results:
                self.blast_results = blast_results
                self.blast_result_file = blast_results.get("result_file")

            raw_results = workspace.get("blast_raw_results")
            if raw_results:
                self.raw_results = raw_results
                self.blast_worker.results = raw_results

        return workspace

    def cleanup(self):
        """Clean up module resources"""
        self.blast_worker = None
        self.blast_results = None
        self.blast_result_file = None
        self.raw_results = None
        self._local_command_history.clear()
