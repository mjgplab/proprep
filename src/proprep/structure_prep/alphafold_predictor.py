"""
AlphaFold Predictor Module

Predicts protein structures from amino acid sequences using AlphaFold/ColabFold.
Integrates with ProPrep workspace for sequence extraction and template selection.
"""

import logging
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.debug_utils import debug_workspace
from proprep.utils.module_registry import ProcessingModule, register_module
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context
from proprep.utils.file_browser import remap_recorded_index_by_key, annotate_recorded_key
from proprep.structure_prep.alphafold_worker import AlphaFoldWorker
from proprep.structure_prep.alphafold_commands import (
    PredictStructureCommand,
    ConfigureParametersCommand,
    ViewResultsCommand,
    ExportBestModelCommand,
    CheckInstallationCommand,
)

logger = logging.getLogger(__name__)


@register_module
class AlphaFoldPredictor(ProcessingModule):
    """Module for predicting protein structures using AlphaFold"""

    NAME = "AlphaFold Predictor"
    DESCRIPTION = "Predict protein structures from sequence using AlphaFold"
    VERSION = "1.0.0"
    CATEGORY = "preparation"
    PRIORITY = 1.5  # Between PDB Loader (1) and PDB Filter (2)

    def initialize(self):
        """Initialize module resources"""
        self.worker = None  # Lazy initialization
        self.prediction_results = None
        self._local_command_history = []
        self._structure_type_registered = False

    def _ensure_worker_initialized(self):
        """Ensure worker is initialized (lazy initialization)."""
        if self.worker is None:
            self.worker = AlphaFoldWorker()
            self._register_structure_type()

    def _register_structure_type(self):
        """Register AlphaFold structure type with the structure selector."""
        if self._structure_type_registered:
            return

        try:
            from proprep.utils.structure_selector import register_structure_type

            register_structure_type(
                workspace_key="alphafold_pdb_file",
                display_name="AlphaFold-Predicted",
                priority=5,  # Same as original_pdb_file
                description="Structure predicted by AlphaFold from sequence"
            )
            self._structure_type_registered = True
            logger.debug("Registered AlphaFold structure type")
        except Exception as e:
            logger.warning(f"Could not register structure type: {e}")

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "check_install": "Check AlphaFold installation",
            "configure": "View/Edit configuration",
            "predict": "Predict structure from sequence",
            "view_results": "View prediction results",
            "export": "Export best model",
        }

    def get_enhanced_menu_options(self, workspace):
        """
        Get menu options with enhanced status information.
        """
        from proprep.utils.enhanced_menu import MenuOption, OptionStatus

        options = []

        # Check if AlphaFold has been run
        alphafold_results = workspace.get("alphafold_pdb_file")
        has_results = alphafold_results is not None

        # Check if ColabFold is installed (lazy check)
        self._ensure_worker_initialized()
        is_installed = self.worker.backend == 'local'

        # Option 1: Check installation
        options.append(MenuOption(
            key="1",
            description="Check AlphaFold installation",
            status=OptionStatus.AVAILABLE
        ))

        # Option 2: View/Edit configuration - only available after first prediction
        config_exists = workspace.get("alphafold_config") is not None
        options.append(MenuOption(
            key="2",
            description="View/Edit configuration",
            status=OptionStatus.COMPLETED if config_exists else OptionStatus.BLOCKED,
            dependency_text="[Run prediction first to set configuration] ○" if not config_exists else ""
        ))

        # Option 3: Predict structure
        if not is_installed:
            predict_status = OptionStatus.AVAILABLE
            predict_desc = "Predict structure from sequence (will use Google Colab)"
        else:
            predict_status = OptionStatus.COMPLETED if has_results else OptionStatus.AVAILABLE
            predict_desc = "Predict structure from sequence"

        options.append(MenuOption(
            key="3",
            description=predict_desc,
            status=predict_status
        ))

        # Option 4: View results - requires prediction to be run
        options.append(MenuOption(
            key="4",
            description="View prediction results",
            status=OptionStatus.READY if has_results else OptionStatus.BLOCKED,
            dependency_text="[Need to run prediction first] ○" if not has_results else ""
        ))

        # Option 5: Export best model - requires prediction to be run
        options.append(MenuOption(
            key="5",
            description="Export best model",
            status=OptionStatus.READY if has_results else OptionStatus.BLOCKED,
            dependency_text="[Need to run prediction first] ○" if not has_results else ""
        ))

        return options

    def get_menu_suggestion(self, workspace):
        """
        Get a suggestion for the next recommended action.
        """
        alphafold_results = workspace.get("alphafold_pdb_file")

        # Ensure worker is initialized to check backend
        self._ensure_worker_initialized()
        is_installed = self.worker.backend == 'local'

        # Check if we have BLAST results (contextual suggestion)
        blast_results = workspace.get("blast_results")

        # Check if parameters are configured
        has_config = workspace.get("alphafold_config") is not None

        if not alphafold_results:
            # No prediction yet - guide user to run first prediction
            if not is_installed:
                return "Check installation status (option 1), then run first prediction (option 3)"
            elif blast_results:
                return "You have BLAST results - run prediction on homologous sequences (option 3)"
            else:
                return "Run first prediction (option 3) - parameters will be set interactively"
        else:
            # Have prediction - now can tweak config and run again
            return "Review/tweak configuration (option 2) or view results (option 4)"

    def handle_menu_option(self, option: str) -> bool:
        """Handle a menu option selection"""
        try:
            if option == "check_install":
                command = CheckInstallationCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "configure":
                command = ConfigureParametersCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "predict":
                command = PredictStructureCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "view_results":
                command = ViewResultsCommand(self.processor)
                return command.execute_with_error_handling()

            elif option == "export":
                command = ExportBestModelCommand(self.processor)
                return command.execute_with_error_handling()

            else:
                logger.warning(f"Unknown menu option: {option}")
                return False

        except Exception as e:
            logger.error(f"Error handling menu option '{option}': {str(e)}")
            return False

    def get_workspace_requirements(self) -> List[str]:
        """No hard requirements - can work standalone"""
        return []

    def get_workspace_outputs(self) -> List[str]:
        """Structures/data this module produces"""
        return [
            "alphafold_pdb_file",
            "alphafold_structure",
            "alphafold_confidence",
            "alphafold_output_dir",
            "alphafold_sequence_source",
            "alphafold_used_template",
            "alphafold_all_models",
            "alphafold_config",
            "alphafold_pae_matrix",
        ]

    # ========================================================================
    # MAIN PREDICTION WORKFLOW
    # ========================================================================

    def predict_structure_interactive(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Interactive AlphaFold prediction workflow.

        This method:
        1. Offers sequence sources (PDB, BLAST, custom)
        2. Offers template structures from workspace
        3. Configures all prediction parameters (transparent)
        4. Runs prediction
        5. Stores results in workspace
        """
        if workspace is None:
            workspace = self.processor._get_workspace()

        # Initialize worker and show status
        self._ensure_worker_initialized()

        self.processor.console.print("\n[bold cyan]═══ AlphaFold Structure Prediction ═══[/bold cyan]\n")

        # Show installation and GPU status
        if self.worker.backend == 'local':
            self.processor.console.print("[green]✓ ColabFold installed locally[/green]")
        else:
            self.processor.console.print("[yellow]ℹ ColabFold not installed locally[/yellow]")
            self.processor.console.print("  Will provide Google Colab instructions\n")

        if self.worker.gpu_info["has_gpu"]:
            gpu_type = self.worker.gpu_info["gpu_type"]
            gpu_details = self.worker.gpu_info["gpu_details"] or "Unknown"
            self.processor.console.print(f"[green]✓ GPU detected: {gpu_type} ({gpu_details})[/green]\n")
        else:
            self.processor.console.print("[grey50]ℹ No GPU detected - predictions will run on CPU[/grey50]\n")

        # Step 1: Get sequence
        sequence, sequence_source = self._get_sequence_interactive(workspace)
        if not sequence:
            self.processor.console.print("[red]No sequence provided. Aborting.[/red]")
            return workspace

        # Step 2: Get template (optional)
        template_file = self._get_template_interactive(workspace)

        # Step 3: Configure parameters
        params = self._configure_all_parameters_interactive(
            workspace,
            has_template=template_file is not None
        )

        # Save configured parameters to workspace for future use
        workspace = self.update_workspace(workspace, "alphafold_config", params)

        # Step 4: Determine output directory
        output_dir = self._get_output_directory()

        # Step 5: Run prediction
        try:
            self.processor.console.print(f"\n[bold]Running AlphaFold prediction...[/bold]")
            self.processor.console.print(f"Sequence length: {len(sequence)} residues")
            if template_file:
                self.processor.console.print(f"Using template: {template_file}")

            results = self.worker.predict(
                sequence=sequence,
                output_dir=output_dir,
                job_name="alphafold_prediction",
                custom_template_path=os.path.dirname(template_file) if template_file else None,
                **params
            )

            # Check if we got Colab instructions instead of results
            if results.get("backend") == "colab":
                self.processor.console.print(Panel(
                    results["instructions"],
                    title="Google Colab Instructions",
                    border_style="yellow"
                ))
                return workspace

            # Step 6: Store results in workspace
            workspace = self._store_results_in_workspace(
                workspace,
                results,
                sequence_source,
                template_file
            )

            self.processor.console.print("\n[green]✓ Prediction completed successfully![/green]")
            self._display_results_summary(results)

            self.prediction_results = results

        except Exception as e:
            self.processor.console.print(f"\n[red]Error during prediction: {str(e)}[/red]")
            logger.error(f"Prediction error: {str(e)}", exc_info=True)

        return workspace

    def _get_sequence_interactive(self, workspace) -> Tuple[Optional[str], str]:
        """
        Interactively get sequence from user with multiple source options.

        Returns:
            Tuple of (sequence, source_description)
        """
        seq_options = self._get_sequence_options(workspace)

        self.processor.console.print("[bold]Sequence Source Options:[/bold]\n")

        choices = []
        choice_map = {}
        choice_num = 1

        # Option: From loaded PDB
        if "pdb" in seq_options:
            pdb_info = seq_options["pdb"]
            self.processor.console.print(f"  {choice_num}. Extract from loaded PDB")
            self.processor.console.print(f"     File: {os.path.basename(pdb_info['file'])}")
            self.processor.console.print(f"     Chains: {', '.join(pdb_info['sequences'].keys())}")
            choices.append(str(choice_num))
            choice_map[str(choice_num)] = "pdb"
            choice_num += 1

        # Option: From BLAST results
        if "blast" in seq_options:
            blast_hits = seq_options["blast"]
            self.processor.console.print(f"  {choice_num}. Select from BLAST homology results")
            self.processor.console.print(f"     Found {len(blast_hits)} hits")
            choices.append(str(choice_num))
            choice_map[str(choice_num)] = "blast"
            choice_num += 1

        # Option: Custom sequence
        self.processor.console.print(f"  {choice_num}. Enter custom sequence (FASTA or plain text)")
        choices.append(str(choice_num))
        choice_map[str(choice_num)] = "custom"

        choice = prompt_with_context(
            self.processor,
            "\nSelect sequence source",
            choices=choices,
            default=choices[0],
            module="AlphaFold Predictor",
            description="Sequence source selection",
            options_map=choice_map
        )

        source_type = choice_map[choice]

        # Get actual sequence based on choice
        if source_type == "pdb":
            return self._select_pdb_sequence(seq_options["pdb"])
        elif source_type == "blast":
            return self._select_blast_sequence(seq_options["blast"])
        else:
            return self._get_custom_sequence()

    def _get_sequence_options(self, workspace) -> Dict[str, Any]:
        """Get available sequence sources from workspace."""
        options = {}

        # Option 1: From loaded PDB structure
        from proprep.utils.structure_selector import get_priority_pdb_file
        pdb_file = get_priority_pdb_file(workspace, self.processor.console, silent=True)

        if pdb_file and os.path.exists(pdb_file):
            sequences = self.worker.extract_sequences_from_pdb(pdb_file)
            if sequences:
                options["pdb"] = {
                    "file": pdb_file,
                    "sequences": sequences
                }

        # Option 2: From BLAST results
        blast_results = workspace.get("blast_results")
        if blast_results and blast_results.get("hits"):
            hits = []
            for hit in blast_results["hits"][:20]:  # Top 20 hits
                if hit.get("sequence"):
                    hits.append({
                        "id": hit.get("accession", hit.get("hit_id", "Unknown")),
                        "description": hit.get("description", ""),
                        "sequence": hit.get("sequence"),
                        "e_value": hit.get("e_value", "N/A"),
                        "identity": hit.get("identity", "N/A")
                    })
            if hits:
                options["blast"] = hits

        # Option 3: Custom input always available
        options["custom"] = True

        return options

    def _select_pdb_sequence(self, pdb_info: Dict) -> Tuple[str, str]:
        """Select a sequence from PDB chains."""
        sequences = pdb_info["sequences"]
        pdb_file = pdb_info["file"]

        if len(sequences) == 1:
            # Only one chain, use it directly
            chain_id = list(sequences.keys())[0]
            sequence = sequences[chain_id]
            source = f"pdb_{os.path.basename(pdb_file)}_chain{chain_id}"

            self.processor.console.print(f"\n[green]Using chain {chain_id} ({len(sequence)} residues)[/green]")
            return sequence, source

        # Multiple chains - let user select
        self.processor.console.print("\n[bold]Available chains:[/bold]")
        chain_list = list(sequences.keys())

        for i, chain_id in enumerate(chain_list, 1):
            seq = sequences[chain_id]
            self.processor.console.print(f"  {i}. Chain {chain_id}: {len(seq)} residues")
            if len(seq) <= 60:
                self.processor.console.print(f"     {seq}")
            else:
                self.processor.console.print(f"     {seq[:60]}...")

        choices = [str(i) for i in range(1, len(chain_list) + 1)]
        choice = prompt_with_context(
            self.processor,
            "Select chain",
            choices=choices,
            default="1",
            module="AlphaFold Predictor",
            description="Chain selection from PDB"
        )

        choice = remap_recorded_index_by_key(self.processor, chain_list, lambda c: c, str(choice))
        chain_id = chain_list[int(choice) - 1]
        annotate_recorded_key(self.processor, chain_id)
        sequence = sequences[chain_id]
        source = f"pdb_{os.path.basename(pdb_file)}_chain{chain_id}"

        return sequence, source

    def _select_blast_sequence(self, blast_hits: List[Dict]) -> Tuple[str, str]:
        """Select a sequence from BLAST hits."""
        self.processor.console.print("\n[bold]BLAST Homology Hits:[/bold]")

        table = Table()
        table.add_column("#", style="cyan")
        table.add_column("ID", style="green")
        table.add_column("Description", style="white")
        table.add_column("E-value", style="yellow")
        table.add_column("Length", style="magenta")

        for i, hit in enumerate(blast_hits[:15], 1):  # Show top 15
            table.add_row(
                str(i),
                hit["id"][:20],
                hit["description"][:40],
                str(hit["e_value"]),
                str(len(hit["sequence"]))
            )

        self.processor.console.print(table)

        choices = [str(i) for i in range(1, min(len(blast_hits), 15) + 1)]
        choice = prompt_with_context(
            self.processor,
            "\nSelect hit to predict",
            choices=choices,
            default="1",
            module="AlphaFold Predictor",
            description="BLAST hit selection"
        )

        choice = remap_recorded_index_by_key(self.processor, blast_hits, lambda h: h["id"], str(choice))
        hit = blast_hits[int(choice) - 1]
        annotate_recorded_key(self.processor, hit["id"])
        sequence = hit["sequence"]
        source = f"blast_{hit['id']}"

        return sequence, source

    def _get_custom_sequence(self) -> Tuple[str, str]:
        """Get custom sequence from user input."""
        self.processor.console.print("\n[bold]Enter sequence (FASTA format or plain text):[/bold]")
        self.processor.console.print("[grey50]Press Enter twice when done[/grey50]\n")

        lines = []
        while True:
            line = prompt_with_context(
                self.processor, "",
                default="",
                module="AlphaFold Predictor",
                description="Custom sequence line (blank to finish)",
            )
            if not line:
                break
            lines.append(line)

        input_text = "\n".join(lines)

        # Parse FASTA or plain sequence
        sequence = ""
        if input_text.startswith(">"):
            # FASTA format
            sequence_lines = []
            for line in lines[1:]:  # Skip header
                if not line.startswith(">"):
                    sequence_lines.append(line.strip())
            sequence = "".join(sequence_lines)
        else:
            # Plain sequence
            sequence = input_text.replace("\n", "").replace(" ", "").strip()

        # Validate
        is_valid, error_msg = self.worker.validate_sequence(sequence)
        if not is_valid:
            self.processor.console.print(f"[red]Invalid sequence: {error_msg}[/red]")
            return None, ""

        return sequence, "user_input"

    def _get_template_interactive(self, workspace) -> Optional[str]:
        """Interactively select template structure from workspace."""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.processor.console)
        available = selector.get_available_structures()

        if not available:
            self.processor.console.print("\n[grey50]No structures available in workspace for templates[/grey50]")
            return None

        self.processor.console.print("\n[bold]Template Structure Options:[/bold]")
        self.processor.console.print("Available structures in workspace:\n")

        table = Table()
        table.add_column("#", style="cyan")
        table.add_column("Type", style="green")
        table.add_column("File", style="white")

        for i, struct_info in enumerate(available, 1):
            table.add_row(
                str(i),
                struct_info.structure_type.display_name,
                os.path.basename(struct_info.file_path)
            )

        self.processor.console.print(table)

        if not confirm_with_context(
            self.processor,
            "\nUse one of these structures as a template?",
            default=True,
            module="AlphaFold Predictor",
            description="Template structure selection"
        ):
            return None

        choices = [str(i) for i in range(1, len(available) + 1)]
        choice = prompt_with_context(
            self.processor,
            "Select structure",
            choices=choices,
            default="1",
            module="AlphaFold Predictor",
            description="Template selection"
        )

        choice = remap_recorded_index_by_key(self.processor, available, lambda s: str(s.file_path), str(choice))
        selected = available[int(choice) - 1]
        annotate_recorded_key(self.processor, str(selected.file_path))
        self.processor.console.print(f"[green]✓ Will use {selected.structure_type.display_name} as template[/green]")

        return selected.file_path

    def _configure_all_parameters_interactive(self, workspace, has_template=False) -> Dict[str, Any]:
        """
        Configure all AlphaFold parameters with full transparency.
        No hidden "Advanced" options - everything is shown.
        """
        self.processor.console.print("\n[bold cyan]═══ Prediction Parameters ═══[/bold cyan]")
        self.processor.console.print("[grey50]All parameters shown for full transparency[/grey50]\n")

        # Check if there's a saved configuration
        saved_config = workspace.get("alphafold_config")
        if saved_config:
            self.processor.console.print("[green]✓ Found saved parameter configuration[/green]")
            if confirm_with_context(
                self.processor,
                "Use saved configuration?",
                default=True,
                module="AlphaFold Predictor",
                description="Use saved AlphaFold configuration"
            ):
                self.processor.console.print("[green]Using saved configuration[/green]\n")
                return saved_config.copy()
            else:
                self.processor.console.print("[yellow]Configuring new parameters...[/yellow]\n")

        params = {}

        # 1. Number of Models
        self.processor.console.print("[bold]1. Number of Models[/bold]")
        self.processor.console.print("   AlphaFold generates multiple predictions, ranked by confidence")
        self.processor.console.print("   More models = better sampling, longer runtime")
        self.processor.console.print("   • 5 models: Best quality (default, ~30-60min on GPU)")
        self.processor.console.print("   • 1-2 models: Quick preview (~10-15min on GPU)")

        params['num_models'] = int(prompt_with_context(
            self.processor,
            "   Number of models",
            choices=["1", "2", "3", "4", "5"],
            default="5",
            module="AlphaFold Predictor",
            description="Number of models to generate"
        ))

        # 2. Recycling Iterations
        self.processor.console.print("\n[bold]2. Recycling Iterations[/bold]")
        self.processor.console.print("   Number of refinement passes through the neural network")
        self.processor.console.print("   More recycles = better geometry (diminishing returns after 5)")
        self.processor.console.print("   • 3: Default, good balance")
        self.processor.console.print("   • 10-20: Maximum quality, much slower")

        params['num_recycle'] = int_prompt_with_context(
            self.processor,
            "   Number of recycles",
            default=3,
            module="AlphaFold Predictor",
            description="Recycling iterations"
        )

        # 3. Template Search
        if not has_template:
            self.processor.console.print("\n[bold]3. Template Search[/bold]")
            self.processor.console.print("   Search PDB for homologous structures to use as templates")
            self.processor.console.print("   Helpful if similar structures exist in PDB")
            self.processor.console.print("   • Yes: Include templates (adds ~5-10min)")
            self.processor.console.print("   • No: Pure ab initio prediction")

            params['use_templates'] = confirm_with_context(
                self.processor,
                "   Use template search",
                default=True,
                module="AlphaFold Predictor",
                description="Template search toggle"
            )
        else:
            # Template provided, force templates on
            params['use_templates'] = True
            self.processor.console.print("\n[bold]3. Template Search[/bold]")
            self.processor.console.print("   [green]Enabled (custom template provided)[/green]")

        # 4. AMBER Relaxation
        self.processor.console.print("\n[bold]4. AMBER Energy Minimization[/bold]")
        self.processor.console.print("   Relax structure using AMBER force field")
        self.processor.console.print("   Fixes clashes, improves geometry")
        self.processor.console.print("   [yellow]HIGHLY RECOMMENDED for MD preparation[/yellow]")
        self.processor.console.print("   • Yes: Better geometry (adds ~30s-2min per model)")
        self.processor.console.print("   • No: Faster, use raw prediction")

        params['amber_relax'] = confirm_with_context(
            self.processor,
            "   Run AMBER relaxation",
            default=True,
            module="AlphaFold Predictor",
            description="AMBER relaxation toggle"
        )

        # 5. GPU Relaxation (if GPU available)
        if self.worker.gpu_info["has_gpu"] and params['amber_relax']:
            self.processor.console.print("\n[bold]5. GPU-Accelerated Relaxation[/bold]")
            self.processor.console.print(f"   GPU detected: {self.worker.gpu_info['gpu_details']}")
            self.processor.console.print("   Use GPU for AMBER relaxation (much faster)")

            params['use_gpu_relax'] = confirm_with_context(
                self.processor,
                "   Use GPU for relaxation",
                default=True,
                module="AlphaFold Predictor",
                description="GPU relaxation toggle"
            )

        # 6. MSA Depth
        self.processor.console.print("\n[bold]6. MSA (Multiple Sequence Alignment) Depth[/bold]")
        self.processor.console.print("   Maximum number of sequences to use from MSA")
        self.processor.console.print("   • auto: Let AlphaFold decide (recommended)")
        self.processor.console.print("   • Custom number: e.g., 512 (faster), 5120 (slower, more accurate)")

        max_msa = prompt_with_context(
            self.processor,
            "   Max MSA sequences",
            default="auto",
            module="AlphaFold Predictor",
            description="MSA depth"
        )
        params['max_msa'] = max_msa

        # 7. Random Seed
        self.processor.console.print("\n[bold]7. Random Seed (Reproducibility)[/bold]")
        self.processor.console.print("   Seed for random number generator")
        self.processor.console.print("   • None: Different results each run")
        self.processor.console.print("   • Fixed number: Reproducible results")

        if confirm_with_context(
            self.processor,
            "   Use fixed random seed for reproducibility",
            default=False,
            module="AlphaFold Predictor",
            description="Random seed toggle"
        ):
            params['random_seed'] = int_prompt_with_context(
                self.processor,
                "   Enter random seed",
                default=42,
                module="AlphaFold Predictor",
                description="Random seed value"
            )

        # 8. Output Options
        self.processor.console.print("\n[bold]8. Output Options[/bold]")

        params['save_all_models'] = confirm_with_context(
            self.processor,
            "   Save all models (not just top-ranked)",
            default=True,
            module="AlphaFold Predictor",
            description="Save all models toggle"
        )

        params['save_pae_json'] = confirm_with_context(
            self.processor,
            "   Save PAE (Predicted Aligned Error) matrix",
            default=True,
            module="AlphaFold Predictor",
            description="Save PAE toggle"
        )

        self.processor.console.print("\n[green]✓ Parameters configured[/green]")

        return params

    def _get_output_directory(self) -> str:
        """Get output directory for prediction results."""
        # Use current working directory + alphafold_results
        output_dir = Path.cwd() / "alphafold_results"
        output_dir.mkdir(parents=True, exist_ok=True)
        return str(output_dir)

    def _store_results_in_workspace(self, workspace, results, sequence_source, template_file):
        """Store AlphaFold results in workspace."""
        workspace = self.update_workspace(workspace, "alphafold_pdb_file", results["best_model_pdb"])
        workspace = self.update_workspace(workspace, "alphafold_structure", results["structure_object"])
        workspace = self.update_workspace(workspace, "alphafold_confidence", results["confidence_metrics"])
        workspace = self.update_workspace(workspace, "alphafold_output_dir", results["output_dir"])
        workspace = self.update_workspace(workspace, "alphafold_sequence_source", sequence_source)
        workspace = self.update_workspace(workspace, "alphafold_all_models", results["all_models"])

        if template_file:
            workspace = self.update_workspace(workspace, "alphafold_used_template", template_file)

        if results.get("pae_matrix"):
            workspace = self.update_workspace(workspace, "alphafold_pae_matrix", results["pae_matrix"])

        return workspace

    def _display_results_summary(self, results):
        """Display summary of prediction results."""
        self.processor.console.print("\n[bold cyan]═══ Prediction Results ═══[/bold cyan]\n")

        # Model information
        self.processor.console.print(f"[bold]Best Model:[/bold] {os.path.basename(results['best_model_pdb'])}")
        self.processor.console.print(f"[bold]Total Models:[/bold] {len(results['all_models'])}")
        self.processor.console.print(f"[bold]Output Directory:[/bold] {results['output_dir']}")

        # Confidence metrics
        if results.get("confidence_metrics"):
            metrics = results["confidence_metrics"]
            self.processor.console.print("\n[bold]Confidence Metrics:[/bold]")

            if "plddt" in metrics:
                plddt = metrics["plddt"]
                if isinstance(plddt, list):
                    avg_plddt = sum(plddt) / len(plddt)
                    self.processor.console.print(f"  Average pLDDT: {avg_plddt:.2f}")
                else:
                    self.processor.console.print(f"  pLDDT: {plddt}")

            if "ptm" in metrics:
                self.processor.console.print(f"  pTM Score: {metrics['ptm']:.3f}")

            if "iptm" in metrics:
                self.processor.console.print(f"  iPTM Score: {metrics['iptm']:.3f}")

        self.processor.console.print("\n[grey50]Use 'View prediction results' to see detailed analysis[/grey50]")

    def _display_config_summary(self, config: Dict[str, Any]):
        """Display a summary of configuration parameters."""
        from rich.table import Table

        table = Table(title="AlphaFold Configuration", show_header=True)
        table.add_column("Parameter", style="cyan")
        table.add_column("Value", style="green")

        # Display parameters in a user-friendly way
        if "num_models" in config:
            table.add_row("Number of Models", str(config["num_models"]))

        if "num_recycle" in config:
            table.add_row("Recycling Iterations", str(config["num_recycle"]))

        if "use_templates" in config:
            table.add_row("Template Search", "Yes" if config["use_templates"] else "No")

        if "amber_relax" in config:
            table.add_row("AMBER Relaxation", "Yes" if config["amber_relax"] else "No")

        if "use_gpu_relax" in config:
            table.add_row("GPU Relaxation", "Yes" if config["use_gpu_relax"] else "No")

        if "max_msa" in config:
            table.add_row("MSA Depth", str(config["max_msa"]))

        if "random_seed" in config:
            table.add_row("Random Seed", str(config["random_seed"]))

        if "save_all_models" in config:
            table.add_row("Save All Models", "Yes" if config["save_all_models"] else "No")

        if "save_pae_json" in config:
            table.add_row("Save PAE Matrix", "Yes" if config["save_pae_json"] else "No")

        self.processor.console.print(table)

    # ========================================================================
    # ADDITIONAL MODULE METHODS
    # ========================================================================

    def configure_parameters(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """View and edit prediction parameters (only available after first prediction)."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        # Ensure worker is initialized
        self._ensure_worker_initialized()

        # Check if configuration exists - if not, guide user
        saved_config = workspace.get("alphafold_config")
        if not saved_config:
            self.processor.console.print("\n[yellow]No configuration found.[/yellow]")
            self.processor.console.print("Please run a prediction first (option 3).")
            self.processor.console.print("Parameters will be configured interactively during your first prediction.\n")
            return workspace

        # Show current backend status
        self.processor.console.print("\n[bold cyan]═══ View/Edit AlphaFold Configuration ═══[/bold cyan]\n")

        if self.worker.backend == 'local':
            self.processor.console.print("[green]✓ Parameters for local ColabFold execution[/green]\n")
        else:
            self.processor.console.print("[yellow]ℹ Parameters for Google Colab instructions[/yellow]\n")

        # Show current configuration
        self.processor.console.print("[green]✓ Current configuration:[/green]\n")
        self._display_config_summary(saved_config)

        if not confirm_with_context(
            self.processor,
            "\nEdit configuration?",
            default=False,
            module="AlphaFold Predictor",
            description="Edit current configuration"
        ):
            self.processor.console.print("[grey50]Configuration unchanged[/grey50]")
            return workspace

        self.processor.console.print("\n[yellow]Updating configuration...[/yellow]\n")

        params = self._configure_all_parameters_interactive(workspace, has_template=False)

        self.processor.console.print("\n[green]✓ Configuration updated[/green]")

        # Update worker defaults
        self.worker.setup(**params)

        # Store in workspace for persistence
        workspace = self.update_workspace(workspace, "alphafold_config", params)

        return workspace

    def view_results(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """View detailed prediction results."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        alphafold_pdb = workspace.get("alphafold_pdb_file")

        if not alphafold_pdb:
            self.processor.console.print("[yellow]No AlphaFold prediction results available[/yellow]")
            return workspace

        self.processor.console.print("\n[bold cyan]═══ AlphaFold Prediction Results ═══[/bold cyan]\n")

        # Display stored information
        self.processor.console.print(f"[bold]Best Model:[/bold] {alphafold_pdb}")
        self.processor.console.print(f"[bold]Sequence Source:[/bold] {workspace.get('alphafold_sequence_source', 'Unknown')}")

        template = workspace.get("alphafold_used_template")
        if template:
            self.processor.console.print(f"[bold]Template Used:[/bold] {template}")

        output_dir = workspace.get("alphafold_output_dir")
        if output_dir:
            self.processor.console.print(f"[bold]Output Directory:[/bold] {output_dir}")

        # Confidence metrics
        confidence = workspace.get("alphafold_confidence")
        if confidence:
            self.processor.console.print("\n[bold]Confidence Metrics:[/bold]")

            if isinstance(confidence, dict):
                for key, value in confidence.items():
                    if isinstance(value, (int, float)):
                        self.processor.console.print(f"  {key}: {value:.3f}")
                    elif isinstance(value, list):
                        avg = sum(value) / len(value)
                        self.processor.console.print(f"  {key} (average): {avg:.3f}")

        # All models
        all_models = workspace.get("alphafold_all_models")
        if all_models:
            self.processor.console.print(f"\n[bold]All Models ({len(all_models)}):[/bold]")
            for i, model in enumerate(all_models, 1):
                self.processor.console.print(f"  {i}. {os.path.basename(model)}")

        return workspace

    def export_best_model(self, workspace: Dict[str, Any] = None, output_file: str = None) -> Dict[str, Any]:
        """Export the best AlphaFold model to a specified location."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        alphafold_pdb = workspace.get("alphafold_pdb_file")

        if not alphafold_pdb:
            self.processor.console.print("[yellow]No AlphaFold prediction results available[/yellow]")
            return workspace

        if not output_file:
            output_file = prompt_with_context(
                self.processor,
                "Enter output file path",
                default="alphafold_best_model.pdb",
                module="AlphaFold Predictor",
                description="Output file path for best AlphaFold model",
            )

        try:
            shutil.copy(alphafold_pdb, output_file)
            self.processor.console.print(f"[green]✓ Exported best model to: {output_file}[/green]")
        except Exception as e:
            self.processor.console.print(f"[red]Error exporting model: {str(e)}[/red]")

        return workspace

    def check_installation(self, workspace: Dict[str, Any] = None) -> Dict[str, Any]:
        """Check AlphaFold/ColabFold installation status."""
        if workspace is None:
            workspace = self.processor._get_workspace()

        # Initialize worker to check status
        self._ensure_worker_initialized()

        self.processor.console.print("\n[bold cyan]═══ AlphaFold Installation Status ═══[/bold cyan]\n")

        # Backend
        if self.worker.backend == 'local':
            self.processor.console.print("[green]✓ ColabFold is installed locally[/green]")
        else:
            self.processor.console.print("[yellow]✗ ColabFold not found locally[/yellow]")
            self.processor.console.print("\n[bold yellow]Important:[/bold yellow] ColabFold cannot be installed in the same")
            self.processor.console.print("environment as ProPrep due to incompatible dependencies")
            self.processor.console.print("(BioPython, TensorFlow, pandas, etc.)")
            self.processor.console.print("\n[bold]Option 1: Use Google Colab (Recommended)[/bold]")
            self.processor.console.print("  • Free GPU access")
            self.processor.console.print("  • No installation required")
            self.processor.console.print("  • Works on any platform")
            self.processor.console.print("  → ProPrep will provide instructions when you run predictions")
            self.processor.console.print("\n[bold]Option 2: Install ColabFold separately (Advanced)[/bold]")
            self.processor.console.print("  Create a separate conda environment:")
            self.processor.console.print("  1. conda create -n colabfold_env python=3.11")
            self.processor.console.print("  2. conda activate colabfold_env")
            self.processor.console.print("  3. conda install -c conda-forge -c bioconda colabfold=1.5.5")
            self.processor.console.print("\n  Or use localcolabfold: https://github.com/YoshitakaMo/localcolabfold")

        # GPU
        gpu_info = self.worker.gpu_info
        if gpu_info["has_gpu"]:
            self.processor.console.print(f"\n[green]✓ GPU detected: {gpu_info['gpu_type']}[/green]")
            if gpu_info["gpu_details"]:
                self.processor.console.print(f"  Details: {gpu_info['gpu_details']}")
        else:
            self.processor.console.print("\n[yellow]✗ No GPU detected[/yellow]")
            self.processor.console.print("  Predictions will run on CPU (much slower)")

        return workspace

    def process(self, workspace: Dict[str, Any]) -> Dict[str, Any]:
        """Process the workspace (for automated workflows)."""
        debug_workspace("AlphaFold", "receiving", workspace, self.processor.console)

        # This module is primarily interactive, so process() is minimal
        # Could add automated prediction from workspace.get("target_sequence") if needed

        debug_workspace("AlphaFold", "returning", workspace, self.processor.console)
        return workspace

    def cleanup(self):
        """Clean up module resources."""
        self.worker = None
        self.prediction_results = None
        self._local_command_history.clear()
