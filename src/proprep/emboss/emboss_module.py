"""
EMBOSS Analysis Module

A comprehensive EMBOSS toolkit for sequence analysis with intelligent workspace integration.
Provides maximum flexibility - works standalone or integrates with existing ProPrep workflows.

Features:
- Universal sequence input (paste, file, workspace, UniProt)
- Sequence analysis (pepstats, pepinfo, charge, iep)
- Pairwise alignment (needle, water, stretcher, diffseq)
- Pattern/motif detection (patmatmotifs, fuzzpro, sigcleave)
- Batch processing capabilities
- Smart workspace awareness and contribution
"""

import logging
import os
import tempfile
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from rich.console import Console
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.status import Status

from proprep.utils.module_registry import ProcessingModule, register_module
from .emboss_commands import (
    SequenceAnalysisCommand,
    PairwiseAlignmentCommand,
    PatternSearchCommand,
    BatchAnalysisCommand,
    ViewEMBOSSResultsCommand,
    ExportEMBOSSResultsCommand,
)
from .emboss_core import EMBOSSCore
from .sequence_input import SequenceInputManager

logger = logging.getLogger(__name__)


@register_module
class EMBOSSModule(ProcessingModule):
    """Comprehensive EMBOSS sequence analysis toolkit"""

    NAME = "EMBOSS Analysis"
    DESCRIPTION = "Sequence analysis tools from the EMBOSS suite"
    VERSION = "1.0.0"
    CATEGORY = "analysis"
    REQUIRES = []  # No hard requirements - fully flexible!
    PRIORITY = 45  # Between homology and other analysis tools

    def initialize(self):
        """Initialize module resources"""
        self.emboss_core = EMBOSSCore()
        self.sequence_manager = SequenceInputManager(self)
        self.results = {}
        self.sequence_cache = {}
        self._local_command_history = []

    def get_workspace_display(self, workspace_key, value, console):
        """
        Custom display method for EMBOSS module data

        Args:
            workspace_key: The key in the workspace dictionary
            value: The value stored in the workspace
            console: The rich console object for output

        Returns:
            bool: True if the module handled the display, False otherwise
        """
        if workspace_key.startswith("emboss_"):
            self._display_emboss_results(workspace_key, value, console)
            return True
        return False

    def _display_emboss_results(self, key, value, console):
        """Display EMBOSS results in formatted tables"""
        if key == "emboss_pepstats" and isinstance(value, dict):
            self._display_pepstats_results(value, console)
        elif key == "emboss_alignments" and isinstance(value, list):
            self._display_alignment_results(value, console)
        elif key == "emboss_motifs" and isinstance(value, list):
            self._display_motif_results(value, console)
        else:
            console.print(f"[cyan]{key}:[/cyan] {value}")

    def get_menu_options(self) -> Dict[str, str]:
        """Get available menu options for EMBOSS analysis"""
        return {
            "sequence_analysis": "Protein Statistics & Properties (pepstats, pepinfo, charge)",
            "pairwise_alignment": "Sequence Alignment (needle, water, stretcher)",
            "pattern_search": "Motif & Pattern Detection (patmatmotifs, fuzzpro)",
            "batch_analysis": "Batch Processing of Multiple Sequences",
            "view_results": "View Previous EMBOSS Results",
            "export_results": "Export Results to Files"
        }

    def get_workspace_requirements(self) -> List[str]:
        """Get workspace requirements - none required, but we'll use what's available"""
        return []  # Fully flexible - no requirements

    def can_process(self, workspace) -> bool:
        """Check if module can process - always true for maximum flexibility"""
        return True

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection"""
        try:
            if option == "sequence_analysis":
                from .emboss_commands import SequenceAnalysisCommand
                command = SequenceAnalysisCommand(self.processor)
                return self.processor.execute_command(command)

            elif option == "pairwise_alignment":
                from .emboss_commands import PairwiseAlignmentCommand
                command = PairwiseAlignmentCommand(self.processor)
                return self.processor.execute_command(command)

            elif option == "pattern_search":
                from .emboss_commands import PatternSearchCommand
                command = PatternSearchCommand(self.processor)
                return self.processor.execute_command(command)

            elif option == "batch_analysis":
                from .emboss_commands import BatchAnalysisCommand
                command = BatchAnalysisCommand(self.processor)
                return self.processor.execute_command(command)

            elif option == "view_results":
                from .emboss_commands import ViewEMBOSSResultsCommand
                command = ViewEMBOSSResultsCommand(self.processor)
                return self.processor.execute_command(command)

            elif option == "export_results":
                from .emboss_commands import ExportEMBOSSResultsCommand
                command = ExportEMBOSSResultsCommand(self.processor)
                return self.processor.execute_command(command)

            else:
                self.processor.console.print(f"[red]Unknown option: {option}[/red]")
                return False

        except Exception as e:
            self.processor.console.print(f"[red]Error executing command: {str(e)}[/red]")
            logger.error(f"Command execution error in EMBOSS module: {str(e)}")
            return False

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Core EMBOSS Actions
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def _sequence_analysis(self, interactive=True):
        """Perform sequence analysis using EMBOSS tools"""
        self.processor.console.print(Panel.fit(
            "[bold]EMBOSS Sequence Analysis[/bold]\n"
            "Analyze protein sequences with pepstats, pepinfo, charge, and iep tools.",
            border_style="blue"
        ))

        # Get sequence input
        sequence_data = self.sequence_manager.get_sequence_input("Analysis Sequence")
        if not sequence_data:
            return False

        sequence, seq_name, seq_source = sequence_data

        # Show analysis options
        analysis_options = {
            "1": "pepstats - Basic protein statistics",
            "2": "pepinfo - Detailed residue-by-residue analysis", 
            "3": "charge - Net charge at different pH values",
            "4": "iep - Isoelectric point calculation",
            "5": "all - Run all analyses"
        }

        self.processor.console.print("\n[bold]Available Analysis Tools:[/bold]")
        for key, desc in analysis_options.items():
            self.processor.console.print(f"  {key}. {desc}")

        if interactive:
            choice = prompt_with_context(
                self.processor,
                "Select analysis",
                choices=list(analysis_options.keys()),
                default="1",
                module="EMBOSS - Sequence Analysis",
                description="Select sequence analysis tool",
                options_map=analysis_options,
            )
        else:
            choice = "5"  # Run all for non-interactive

        # Execute selected analyses
        results = {}

        if choice in ["1", "5"]:
            results["pepstats"] = self.emboss_core.run_pepstats(sequence, seq_name)
            
        if choice in ["2", "5"]:
            results["pepinfo"] = self.emboss_core.run_pepinfo(sequence, seq_name)
            
        if choice in ["3", "5"]:
            results["charge"] = self.emboss_core.run_charge(sequence, seq_name)
            
        if choice in ["4", "5"]:
            results["iep"] = self.emboss_core.run_iep(sequence, seq_name)

        # Display results
        self._display_analysis_results(results, seq_name)

        # Store results in workspace
        self.update_workspace("emboss_sequence_analysis", {
            "sequence_name": seq_name,
            "sequence_source": seq_source,
            "results": results,
            "timestamp": self.emboss_core.get_timestamp()
        })

        return True

    def _pairwise_alignment(self, interactive=True):
        """Perform pairwise sequence alignment"""
        self.processor.console.print(Panel.fit(
            "[bold]EMBOSS Pairwise Alignment[/bold]\n"
            "Align two sequences using needle, water, or stretcher algorithms.",
            border_style="blue"
        ))

        # Get first sequence
        seq1_data = self.sequence_manager.get_sequence_input("First Sequence")
        if not seq1_data:
            return False
        seq1, name1, source1 = seq1_data

        # Get second sequence
        seq2_data = self.sequence_manager.get_sequence_input("Second Sequence")
        if not seq2_data:
            return False
        seq2, name2, source2 = seq2_data

        # Show alignment options
        alignment_options = {
            "1": "needle - Global alignment (Needleman-Wunsch)",
            "2": "water - Local alignment (Smith-Waterman)",
            "3": "stretcher - Global alignment optimized for long sequences",
            "4": "all - Run all alignment methods"
        }

        self.processor.console.print("\n[bold]Available Alignment Methods:[/bold]")
        for key, desc in alignment_options.items():
            self.processor.console.print(f"  {key}. {desc}")

        if interactive:
            choice = prompt_with_context(
                self.processor,
                "Select alignment method",
                choices=list(alignment_options.keys()),
                default="1",
                module="EMBOSS - Pairwise Alignment",
                description="Select pairwise alignment method",
                options_map=alignment_options,
            )
        else:
            choice = "1"  # Default to needle

        # Execute alignment
        results = {}
        
        if choice in ["1", "4"]:
            results["needle"] = self.emboss_core.run_needle(seq1, seq2, name1, name2)
            
        if choice in ["2", "4"]:
            results["water"] = self.emboss_core.run_water(seq1, seq2, name1, name2)
            
        if choice in ["3", "4"]:
            results["stretcher"] = self.emboss_core.run_stretcher(seq1, seq2, name1, name2)

        # Display results
        self._display_alignment_results_detailed(results, name1, name2)

        # Store results in workspace
        alignment_result = {
            "sequence1": {"name": name1, "source": source1},
            "sequence2": {"name": name2, "source": source2},
            "results": results,
            "timestamp": self.emboss_core.get_timestamp()
        }
        
        # Add to list of alignments
        existing_alignments = self.get_from_workspace("emboss_alignments", [])
        existing_alignments.append(alignment_result)
        self.update_workspace("emboss_alignments", existing_alignments)

        return True

    def _pattern_search(self, interactive=True):
        """Search for patterns and motifs in sequences"""
        self.processor.console.print(Panel.fit(
            "[bold]EMBOSS Pattern & Motif Search[/bold]\n"
            "Search for functional motifs and patterns using PROSITE and custom patterns.",
            border_style="blue"
        ))

        # Get sequence input
        sequence_data = self.sequence_manager.get_sequence_input("Search Sequence")
        if not sequence_data:
            return False

        sequence, seq_name, seq_source = sequence_data

        # Show search options
        search_options = {
            "1": "patmatmotifs - Search PROSITE motif database",
            "2": "fuzzpro - Search for patterns with mismatches",
            "3": "sigcleave - Signal peptide cleavage sites",
            "4": "all - Run all searches"
        }

        self.processor.console.print("\n[bold]Available Search Tools:[/bold]")
        for key, desc in search_options.items():
            self.processor.console.print(f"  {key}. {desc}")

        if interactive:
            choice = prompt_with_context(
                self.processor,
                "Select search method",
                choices=list(search_options.keys()),
                default="1",
                module="EMBOSS - Pattern Search",
                description="Select pattern search method",
                options_map=search_options,
            )
        else:
            choice = "4"  # Run all for non-interactive

        # Execute searches
        results = {}

        if choice in ["1", "4"]:
            results["patmatmotifs"] = self.emboss_core.run_patmatmotifs(sequence, seq_name)

        if choice in ["2", "4"]:
            # Get custom pattern if interactive
            pattern = None
            if interactive and choice == "2":
                pattern = prompt_with_context(
                    self.processor,
                    "Enter pattern to search (or press Enter for default)",
                    module="EMBOSS - Pattern Search",
                    description="Enter custom pattern for fuzzpro",
                )
            results["fuzzpro"] = self.emboss_core.run_fuzzpro(sequence, seq_name, pattern)
            
        if choice in ["3", "4"]:
            results["sigcleave"] = self.emboss_core.run_sigcleave(sequence, seq_name)

        # Display results
        self._display_pattern_results(results, seq_name)

        # Store results in workspace
        pattern_result = {
            "sequence_name": seq_name,
            "sequence_source": seq_source,
            "results": results,
            "timestamp": self.emboss_core.get_timestamp()
        }
        
        existing_motifs = self.get_from_workspace("emboss_motifs", [])
        existing_motifs.append(pattern_result)
        self.update_workspace("emboss_motifs", existing_motifs)

        return True

    def _batch_analysis(self, interactive=True):
        """Perform batch analysis on multiple sequences"""
        self.processor.console.print(Panel.fit(
            "[bold]EMBOSS Batch Analysis[/bold]\n"
            "Process multiple sequences from workspace or files simultaneously.",
            border_style="blue"
        ))

        # Get sequences for batch processing
        sequences = self.sequence_manager.get_batch_sequences()
        if not sequences:
            self.processor.console.print("[yellow]No sequences available for batch processing.[/yellow]")
            return False

        self.processor.console.print(f"\n[green]Found {len(sequences)} sequences for batch processing:[/green]")
        for i, (name, seq, source) in enumerate(sequences, 1):
            self.processor.console.print(f"  {i}. {name} ({len(seq)} residues) - {source}")

        if interactive:
            proceed = confirm_with_context(
                self.processor,
                "\nProceed with batch analysis?",
                default=True,
                module="EMBOSS - Batch Analysis",
                description="Confirm proceed with batch sequence analysis",
            )
            if not proceed:
                return False

        # Run batch analysis
        batch_results = []
        
        with Status("[bold green]Running batch analysis...") as status:
            for i, (name, seq, source) in enumerate(sequences):
                status.update(f"Processing {name} ({i+1}/{len(sequences)})...")
                
                seq_results = {
                    "name": name,
                    "source": source,
                    "length": len(seq),
                    "pepstats": self.emboss_core.run_pepstats(seq, name),
                    "charge": self.emboss_core.run_charge(seq, name),
                    "patmatmotifs": self.emboss_core.run_patmatmotifs(seq, name)
                }
                batch_results.append(seq_results)

        # Display summary
        self._display_batch_results(batch_results)

        # Store results in workspace
        self.update_workspace("emboss_batch_analysis", {
            "results": batch_results,
            "sequence_count": len(sequences),
            "timestamp": self.emboss_core.get_timestamp()
        })

        return True

    def _view_results(self):
        """View previous EMBOSS results from workspace"""
        # Implementation for viewing results
        pass

    def _export_results(self):
        """Export EMBOSS results to files"""
        # Implementation for exporting results
        pass

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Display Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def _display_analysis_results(self, results, seq_name):
        """Display sequence analysis results"""
        self.processor.console.print(f"\n[bold]Analysis Results for {seq_name}:[/bold]")
        
        for tool, result in results.items():
            if result and result.get("success"):
                self.processor.console.print(f"\n[green]✓ {tool.upper()} Results:[/green]")
                if tool == "pepstats":
                    self._display_pepstats_table(result)
                elif tool == "charge":
                    self._display_charge_table(result)
                # Add more display methods as needed
            else:
                self.processor.console.print(f"[red]✗ {tool.upper()} failed[/red]")

    def _display_pepstats_table(self, result):
        """Display pepstats results in a formatted table"""
        if not result.get("properties"):
            return
            
        table = Table(title="Protein Statistics")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")
        
        props = result["properties"]
        table.add_row("Molecular Weight", f"{props.get('molecular_weight', 'N/A')} Da")
        table.add_row("Isoelectric Point", f"{props.get('isoelectric_point', 'N/A')}")
        table.add_row("Length", f"{props.get('length', 'N/A')} residues")
        table.add_row("Net Charge", f"{props.get('net_charge', 'N/A')}")
        
        self.processor.console.print(table)

    def _display_charge_table(self, result):
        """Display charge analysis results"""
        if not result.get("charges"):
            return
            
        table = Table(title="Charge vs pH")
        table.add_column("pH", style="cyan")
        table.add_column("Net Charge", style="green")
        
        for ph, charge in result["charges"].items():
            table.add_row(f"{ph:.1f}", f"{charge:.2f}")
            
        self.processor.console.print(table)

    def _display_alignment_results_detailed(self, results, name1, name2):
        """Display alignment results"""
        self.processor.console.print(f"\n[bold]Alignment Results: {name1} vs {name2}[/bold]")
        
        for method, result in results.items():
            if result and result.get("success"):
                self.processor.console.print(f"\n[green]✓ {method.upper()} Alignment:[/green]")
                if result.get("summary"):
                    summary = result["summary"]
                    self.processor.console.print(f"Score: {summary.get('score', 'N/A')}")
                    self.processor.console.print(f"Identity: {summary.get('identity', 'N/A')}%")
                    self.processor.console.print(f"Similarity: {summary.get('similarity', 'N/A')}%")
            else:
                self.processor.console.print(f"[red]✗ {method.upper()} failed[/red]")

    def _display_pattern_results(self, results, seq_name):
        """Display pattern search results"""
        self.processor.console.print(f"\n[bold]Pattern Search Results for {seq_name}:[/bold]")
        
        for tool, result in results.items():
            if result and result.get("success"):
                motifs = result.get("motifs", [])
                self.processor.console.print(f"\n[green]✓ {tool.upper()}: Found {len(motifs)} matches[/green]")
                
                if motifs:
                    table = Table()
                    table.add_column("Position", style="cyan")
                    table.add_column("Motif", style="green")
                    table.add_column("Description", style="yellow")
                    
                    for motif in motifs[:5]:  # Show first 5
                        table.add_row(
                            f"{motif.get('start', '')}-{motif.get('end', '')}",
                            motif.get('sequence', ''),
                            motif.get('description', '')
                        )
                    
                    self.processor.console.print(table)
                    
                    if len(motifs) > 5:
                        self.processor.console.print(f"[grey50]... and {len(motifs) - 5} more matches[/grey50]")
            else:
                self.processor.console.print(f"[red]✗ {tool.upper()} failed[/red]")

    def _display_batch_results(self, batch_results):
        """Display batch analysis results summary"""
        table = Table(title="Batch Analysis Summary")
        table.add_column("Sequence", style="cyan")
        table.add_column("Length", style="green")
        table.add_column("MW (Da)", style="yellow")
        table.add_column("pI", style="magenta")
        table.add_column("Motifs", style="blue")
        
        for result in batch_results:
            pepstats = result.get("pepstats", {})
            props = pepstats.get("properties", {})
            motifs = result.get("patmatmotifs", {}).get("motifs", [])
            
            table.add_row(
                result["name"],
                str(result["length"]),
                str(props.get("molecular_weight", "N/A")),
                str(props.get("isoelectric_point", "N/A")),
                str(len(motifs))
            )
        
        self.processor.console.print(table)

    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
    # Utility Methods
    # =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

    def get_workspace_outputs(self) -> List[str]:
        """List what this module contributes to workspace"""
        return [
            "emboss_sequence_analysis",
            "emboss_alignments",
            "emboss_motifs",
            "emboss_batch_analysis",
        ]

    def get_breadcrumb_for_command(self, command):
        """Return breadcrumb string for a command"""
        action_to_run = getattr(command, "_action_to_run", None)
        if action_to_run:
            action_map = {
                "_sequence_analysis": "Sequence Analysis",
                "_pairwise_alignment": "Pairwise Alignment", 
                "_pattern_search": "Pattern Search",
                "_batch_analysis": "Batch Analysis",
                "_view_results": "View Results",
                "_export_results": "Export Results"
            }
            action_name = action_map.get(action_to_run, action_to_run)
            return f"{self.NAME} > {action_name}"
        return self.NAME

    def _record_local_command(self, command):
        """Record command in local history for module-specific tracking"""
        self._local_command_history.append(command)
