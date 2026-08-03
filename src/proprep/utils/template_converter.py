"""
Template Converter for ProPrep Session Files

Converts recorded session files into reusable templates with variable substitution.
Enables batch processing of multiple proteins using the same workflow.
"""

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

from rich.console import Console
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.panel import Panel
from rich.table import Table


class TemplateConverter:
    """Convert recorded sessions into reusable templates."""

    def __init__(self, console: Optional[Console] = None):
        """Initialize the template converter."""
        self.console = console or Console()

    def convert_to_template(
        self,
        session_file: str,
        output_file: Optional[str] = None,
        auto_detect: bool = True,
        interactive: bool = True
    ) -> str:
        """
        Convert a recorded session to a template.

        Args:
            session_file: Path to the recorded session file
            output_file: Path for the output template file (optional)
            auto_detect: Automatically detect variables to templatize
            interactive: Prompt user for confirmation

        Returns:
            Path to the generated template file
        """
        # Load the session
        session_data = self._load_session(session_file)

        # Detect potential variables
        detected_vars = self._detect_variables(session_data) if auto_detect else {}

        if not detected_vars:
            self.console.print("[yellow]No variables detected for templatization[/yellow]")
            self.console.print("[grey50]Make sure the session includes protein loading[/grey50]")
            return None

        # Interactive confirmation if requested
        if interactive and detected_vars:
            detected_vars = self._confirm_variables(detected_vars, session_data)

        if not detected_vars:
            self.console.print("[yellow]No variables confirmed. Template creation cancelled.[/yellow]")
            return None

        # Apply variable substitution
        template_data = self._create_template(session_data, detected_vars)

        # Generate output filename if not provided
        if not output_file:
            session_path = Path(session_file)
            output_file = str(session_path.parent / f"{session_path.stem}_template.json")

        # Save template
        self._save_template(template_data, output_file)

        self.console.print(f"\n[bold green]✓ Template created: {output_file}[/bold green]")
        # Surface --web when the source session was recorded under proprep-web
        # so the batch replays in the same mode (mode-gated prompts must match).
        web_recorded = bool(session_data.get("metadata", {}).get("web_shell_mode", False))
        web_flag = " --web" if web_recorded else ""
        self.console.print(f"[grey50]Use with: proprep --batch-replay {output_file} "
                           f"--batch-list proteins.txt{web_flag}[/grey50]")

        return output_file

    def _load_session(self, session_file: str) -> Dict[str, Any]:
        """Load session data from file."""
        try:
            with open(session_file, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            raise ValueError(f"Session file not found: {session_file}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid session file format: {e}")

    def _detect_variables(self, session_data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """
        Automatically detect variables that should be templatized.

        Looks for protein input in the session based on the Structure Loader flow:
        1. Find "Select source" interaction
        2. Determine which source was chosen (RCSB/AlphaFold/Local)
        3. Find the corresponding protein input

        Returns:
            Dictionary of variable definitions with their metadata
        """
        variables = {}
        interactions = session_data.get("interactions", [])

        # Find the structure loading sequence
        for i, interaction in enumerate(interactions):
            context = interaction.get("context", {})

            # Look for "Select structure source" prompt
            if (context.get("description") == "Select structure source" and
                context.get("module") == "Structure Loader"):

                source_choice = interaction.get("response")

                # Now find the protein input based on the source choice
                protein_input = self._find_protein_input(interactions, i, source_choice)

                if protein_input:
                    variables["input_protein"] = {
                        "description": "PDB ID, AlphaFold ID, or local file path",
                        "type": "protein_input",
                        "example": protein_input["value"],
                        "required": True,
                        "interaction_index": protein_input["index"],
                        "prompt": protein_input["prompt"],
                        "original_value": protein_input["value"],
                        "source_type": self._classify_source(source_choice, protein_input["value"]),
                        "source_choice": source_choice
                    }
                    # Found the protein input, we're done
                    break

        return variables

    def _find_protein_input(
        self,
        interactions: List[Dict[str, Any]],
        source_index: int,
        source_choice: str
    ) -> Optional[Dict[str, Any]]:
        """
        Find the protein input interaction following the source selection.

        Args:
            interactions: List of all interactions
            source_index: Index of the "Select source" interaction
            source_choice: The choice made ("1" for RCSB, "2" for AlphaFold, "3" for local)

        Returns:
            Dictionary with protein input info or None
        """
        # Search forward from the source selection
        for i in range(source_index + 1, len(interactions)):
            interaction = interactions[i]
            context = interaction.get("context", {})
            description = context.get("description", "")
            module = context.get("module", "")

            # RCSB PDB path (choice "1")
            if source_choice == "1":
                # Look for "Enter PDB ID to download"
                if "Enter PDB ID to download" in description:
                    return {
                        "index": i,
                        "value": interaction.get("response"),
                        "prompt": interaction.get("prompt"),
                        "description": description
                    }
                # Also check for search results selection (if user searched first)
                # This would be more complex and we'll handle in a future version

            # AlphaFold path (choice "2")
            elif source_choice == "2":
                # Look for UniProt ID or gene name prompts
                if ("UniProt" in description or
                    "gene name" in description or
                    "AlphaFold" in module):
                    # AlphaFold prompts vary - look for the actual input
                    if interaction.get("type") == "prompt" and not interaction.get("choices"):
                        return {
                            "index": i,
                            "value": interaction.get("response"),
                            "prompt": interaction.get("prompt"),
                            "description": description
                        }

            # Local file path (choice "3")
            elif source_choice == "3":
                # Look for file path entry or file selection
                if ("file path" in description.lower() or
                    "select file" in description.lower() or
                    "local" in description.lower()):
                    # File path input
                    if interaction.get("type") == "prompt" and not interaction.get("choices"):
                        value = interaction.get("response")
                        # Check if it looks like a file path or file selection
                        if value and ("/" in value or "\\" in value or value.endswith(".pdb")):
                            return {
                                "index": i,
                                "value": value,
                                "prompt": interaction.get("prompt"),
                                "description": description
                            }

            # Stop if we hit another major menu (went too far)
            if context.get("description") == "Select structure source":
                break

        return None

    def _classify_source(self, source_choice: str, value: str) -> str:
        """
        Classify what type of protein input this is.

        Args:
            source_choice: The source menu choice ("1", "2", "3")
            value: The input value

        Returns:
            Type classification: "rcsb_pdb", "alphafold", "local_file"
        """
        if source_choice == "1":
            return "rcsb_pdb"
        elif source_choice == "2":
            return "alphafold"
        elif source_choice == "3":
            return "local_file"

        # Fallback: try to guess from the value
        if re.match(r'^[A-Za-z0-9]{4}$', value.strip()):
            return "rcsb_pdb"
        elif "/" in value or "\\" in value or value.endswith(".pdb"):
            return "local_file"
        else:
            return "unknown"

    def _confirm_variables(
        self,
        detected_vars: Dict[str, Dict[str, Any]],
        session_data: Dict[str, Any]
    ) -> Dict[str, Dict[str, Any]]:
        """
        Interactively confirm detected variables with user.

        Returns:
            Confirmed variables (may be modified or removed)
        """
        self.console.print("\n[bold cyan]Template Conversion: Variable Detection[/bold cyan]\n")

        confirmed_vars = {}

        for var_name, var_info in detected_vars.items():
            # Display what was detected
            source_labels = {
                "rcsb_pdb": "RCSB Protein Data Bank",
                "alphafold": "AlphaFold Database",
                "local_file": "Local PDB File"
            }

            self.console.print(Panel(
                f"[bold]Variable:[/bold] {var_name}\n"
                f"[bold]Source Type:[/bold] {source_labels.get(var_info['source_type'], 'Unknown')}\n"
                f"[bold]Detected Value:[/bold] {var_info['original_value']}\n"
                f"[bold]Prompt:[/bold] {var_info['prompt']}\n"
                f"[bold]Description:[/bold] {var_info.get('description', 'N/A')}",
                title="Detected Protein Input",
                border_style="cyan"
            ))

            # Confirm with user
            if confirm_with_context(None,
                f"\n[yellow]Make this a template variable?[/yellow]",
                default=True
            ):
                # Allow customization of variable name
                custom_name = prompt_with_context(None,
                    "Variable name",
                    default=var_name
                )

                confirmed_vars[custom_name] = var_info
                self.console.print(f"[green]✓ Added variable: {custom_name}[/green]")
            else:
                self.console.print(f"[yellow]Skipped variable: {var_name}[/yellow]")

        return confirmed_vars

    def _create_template(
        self,
        session_data: Dict[str, Any],
        variables: Dict[str, Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Create template data structure with variable substitution.

        Args:
            session_data: Original session data
            variables: Variable definitions

        Returns:
            Template data structure
        """
        # Get original metadata
        original_metadata = session_data.get("metadata", {})

        # Create template metadata
        template_data = {
            "version": "1.2",
            "template": True,
            "template_metadata": {
                "created_at": datetime.now().isoformat(),
                "created_from": Path(original_metadata.get("working_directory", "")).name,
                "description": "ProPrep batch processing template",
                "proprep_version": original_metadata.get("proprep_version", "unknown")
            },
            "template_variables": {},
            "metadata": original_metadata.copy(),
            "interactions": session_data.get("interactions", []).copy()
        }

        # Add variable definitions (cleaned for template)
        for var_name, var_info in variables.items():
            template_data["template_variables"][var_name] = {
                "description": var_info["description"],
                "type": var_info["type"],
                "example": var_info["example"],
                "required": var_info["required"],
                "source_type": var_info.get("source_type", "unknown")
            }

            # Substitute in interactions
            interaction_idx = var_info.get("interaction_index")
            if interaction_idx is not None and interaction_idx < len(template_data["interactions"]):
                # Replace the response with template variable
                template_data["interactions"][interaction_idx]["response"] = f"{{{{ {var_name} }}}}"

                # Add template marker to context
                if "context" not in template_data["interactions"][interaction_idx]:
                    template_data["interactions"][interaction_idx]["context"] = {}
                template_data["interactions"][interaction_idx]["context"]["template_variable"] = var_name

        return template_data

    def _save_template(self, template_data: Dict[str, Any], output_file: str):
        """Save template to file."""
        with open(output_file, 'w') as f:
            json.dump(template_data, f, indent=2)

    def validate_template(self, template_file: str) -> Tuple[bool, List[str]]:
        """
        Validate a template file.

        Args:
            template_file: Path to template file

        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []

        try:
            with open(template_file, 'r') as f:
                template_data = json.load(f)
        except FileNotFoundError:
            return False, [f"Template file not found: {template_file}"]
        except json.JSONDecodeError as e:
            return False, [f"Invalid JSON format: {e}"]

        # Check required fields
        if not template_data.get("template"):
            errors.append("Missing 'template' flag (should be true)")

        if "template_variables" not in template_data:
            errors.append("Missing 'template_variables' section")

        if "interactions" not in template_data:
            errors.append("Missing 'interactions' section")

        # Validate template variables
        variables = template_data.get("template_variables", {})
        for var_name, var_def in variables.items():
            if "type" not in var_def:
                errors.append(f"Variable '{var_name}' missing 'type' field")
            if "description" not in var_def:
                errors.append(f"Variable '{var_name}' missing 'description' field")

        # Check that variables are actually used in interactions
        interactions = template_data.get("interactions", [])
        for var_name in variables.keys():
            var_pattern = f"{{{{ {var_name} }}}}"
            found = any(var_pattern in str(interaction.get("response", ""))
                       for interaction in interactions)
            if not found:
                errors.append(f"Variable '{var_name}' defined but not used in interactions")

        return (len(errors) == 0, errors)

    def display_template_info(self, template_file: str):
        """Display information about a template file."""
        try:
            with open(template_file, 'r') as f:
                template_data = json.load(f)
        except Exception as e:
            self.console.print(f"[red]Error loading template: {e}[/red]")
            return

        # Display basic info
        metadata = template_data.get("template_metadata", {})
        self.console.print(f"\n[bold cyan]Template Information[/bold cyan]")
        self.console.print(f"[bold]File:[/bold] {template_file}")
        self.console.print(f"[bold]Version:[/bold] {template_data.get('version', 'unknown')}")
        self.console.print(f"[bold]Created:[/bold] {metadata.get('created_at', 'unknown')}")
        self.console.print(f"[bold]Description:[/bold] {metadata.get('description', 'No description')}")

        # Display variables
        variables = template_data.get("template_variables", {})
        if variables:
            self.console.print(f"\n[bold]Template Variables:[/bold]")
            table = Table(show_header=True, header_style="bold magenta")
            table.add_column("Variable", style="cyan")
            table.add_column("Type", style="green")
            table.add_column("Source", style="yellow")
            table.add_column("Example", style="blue")

            source_labels = {
                "rcsb_pdb": "RCSB PDB",
                "alphafold": "AlphaFold",
                "local_file": "Local File"
            }

            for var_name, var_def in variables.items():
                table.add_row(
                    var_name,
                    var_def.get("type", "unknown"),
                    source_labels.get(var_def.get("source_type"), "Unknown"),
                    str(var_def.get("example", ""))
                )

            self.console.print(table)

        # Display interaction count
        interactions = template_data.get("interactions", [])
        self.console.print(f"\n[bold]Total Interactions:[/bold] {len(interactions)}")


PLACEHOLDER_PATTERN = r'\{\{\s*(\w+)\s*\}\}'


def find_response_placeholders(template_data: Dict[str, Any]) -> Dict[str, List[int]]:
    """
    Map each {{ variable }} placeholder to the interactions that use it.

    Substitution only ever rewrites interaction responses (see
    SessionReplayer._substitute_variables), so a response is the only place a
    placeholder can do any work.

    Args:
        template_data: Parsed template file

    Returns:
        Dict of variable name -> list of interaction indices using it
    """
    found: Dict[str, List[int]] = {}
    for idx, interaction in enumerate(template_data.get("interactions", [])):
        response = interaction.get("response", "")
        if not isinstance(response, str):
            continue
        for var_name in re.findall(PLACEHOLDER_PATTERN, response):
            found.setdefault(var_name, []).append(idx)
    return found


def validate_variable_coverage(
    template_data: Dict[str, Any],
    rows: List[Dict[str, str]]
) -> Tuple[List[str], List[str]]:
    """
    Check a template against the runs that will drive it.

    validate_template checks the template's shape in isolation, which cannot
    catch a template that is well-formed but unsatisfiable by a given input
    list. Both failures below otherwise surface only once a batch is underway,
    as a prompt with no recorded answer:

      - a variable no row supplies. If declared, SessionReplayer raises during
        construction, before the batch loop's restoring finally block is in
        scope; if undeclared, the placeholder survives substitution and replay
        reaches a prompt it cannot answer.
      - a variable used as a placeholder but never declared, which the
        required-variable check never sees because it iterates declarations
        rather than uses.

    Args:
        template_data: Parsed template file
        rows: Normalized runs, each a dict of variable name -> value

    Returns:
        Tuple of (errors, warnings). Errors mean the batch cannot run to
        completion; warnings are likely mistakes that will not stall it.
    """
    errors: List[str] = []
    warnings: List[str] = []

    declared = template_data.get("template_variables", {})
    placeholders = find_response_placeholders(template_data)

    for var_name, indices in sorted(placeholders.items()):
        if var_name not in declared:
            warnings.append(
                f"Placeholder '{{{{ {var_name} }}}}' (interaction {indices[0]}) "
                f"is not declared in template_variables"
            )

    # Every placeholder needs a value whether or not it was declared, and every
    # variable declared required needs one whether or not it is used -- the
    # latter because SessionReplayer enforces declarations, not uses.
    needed = set(placeholders)
    needed.update(name for name, d in declared.items() if d.get("required", True))

    for var_name in sorted(needed):
        missing = [i for i, row in enumerate(rows, 1)
                   if not str(row.get(var_name, "")).strip()]
        if not missing:
            continue

        used_at = placeholders.get(var_name)
        role = (f"substituted into interaction {used_at[0]}" if used_at
                else "declared required but never used")
        scope = ("no run supplies a value" if len(missing) == len(rows)
                 else "no value in run(s) " + ", ".join(str(i) for i in missing))
        errors.append(f"Variable '{var_name}' ({role}): {scope}")

    known = set(declared) | set(placeholders)
    for var_name in sorted({key for row in rows for key in row} - known):
        warnings.append(
            f"Input column '{var_name}' matches no template variable "
            f"and will be ignored (check the header for a typo)"
        )

    return (errors, warnings)


def create_template_from_session(
    session_file: str,
    output_file: Optional[str] = None,
    interactive: bool = True
) -> str:
    """
    Convenience function to create a template from a session file.

    Args:
        session_file: Path to recorded session file
        output_file: Output template file path (optional)
        interactive: Whether to prompt for confirmation

    Returns:
        Path to created template file
    """
    converter = TemplateConverter()
    return converter.convert_to_template(
        session_file=session_file,
        output_file=output_file,
        interactive=interactive
    )
