"""
Interactive Session Editor for ProPrep

Provides a terminal-based editor for modifying recorded session files.
Supports two modes:
  - "edit": Edit/delete interactions, save as modified session
  - "template": Mark interactions as variables, save as template with batch file generation

Inspired by git interactive rebase — mark actions on interactions, then apply.
"""

import csv
import json
import os
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.utils.session_recorder import (
    format_interaction_display,
    safe_load_session_file,
)
from proprep.utils.template_converter import TemplateConverter


class InteractiveSessionEditor:
    """Interactive terminal editor for ProPrep session files."""

    PAGE_SIZE = 20

    def __init__(self, session_file: str, mode: str = "edit"):
        """
        Initialize the session editor.

        Args:
            session_file: Path to session file to edit
            mode: "edit" for session editing, "template" for template creation
        """
        self.session_file = session_file
        self.mode = mode
        self.console = Console()

        # Load session
        self.session_data = safe_load_session_file(session_file, auto_recover=True)
        if self.session_data is None:
            raise ValueError(f"Could not load session file: {session_file}")

        interactions = self.session_data.get("interactions", [])
        self.metadata = self.session_data.get("metadata", {})

        # Build display list (filter out standalone INPUT entries)
        self.display_items: List[Tuple[int, int, dict]] = []
        for actual_idx, interaction in enumerate(interactions):
            if interaction.get("type", "").lower() != "input":
                disp_idx = len(self.display_items)
                self.display_items.append((disp_idx, actual_idx, interaction))

        # Edit state
        self.edits: Dict[int, str] = {}  # display_idx -> new_value
        self.deletes: Set[int] = set()  # display indices to skip
        self.variables: Dict[int, dict] = {}  # display_idx -> variable info

        # Pagination
        self.current_page = 0
        self.total_pages = max(1, (len(self.display_items) + self.PAGE_SIZE - 1) // self.PAGE_SIZE)

        # Auto-detect template variables in template mode
        # Disabled: don't auto-mark the protein input (e.g. PDB ID) as a
        # variable. Users mark variables manually with 'v #'.
        # if self.mode == "template":
        #     self._auto_detect_variables()

    # ── Core loop ──────────────────────────────────────────────

    def run(self) -> Optional[str]:
        """
        Run the interactive editor.

        Returns:
            Path to saved file, or None if user quit without saving.
        """
        self._show_header()
        self.console.print(
            "[yellow]Note:[/yellow] Editing a menu choice that changes the workflow path (e.g., picking\n"
            "a different tool) may invalidate later interactions. During replay, ProPrep\n"
            "will fall through to live input when recorded interactions no longer match.\n"
        )
        self._show_page()

        while True:
            cmd = self._get_command()
            if cmd is None:
                continue

            action = cmd[0]
            args = cmd[1:]

            if action == "edit" and args:
                self._edit_interaction(args[0])
                self._show_page()
            elif action == "delete" and args:
                if len(args) == 2:
                    self._delete_range(args[0], args[1])
                else:
                    self._delete_interaction(args[0])
                self._show_page()
            elif action == "undo" and args:
                if len(args) == 2:
                    self._undo_range(args[0], args[1])
                else:
                    self._undo_change(args[0])
                self._show_page()
            elif action == "variable" and args:
                if self.mode != "template":
                    self.console.print("[yellow]Variable marking is only available in template mode[/yellow]")
                else:
                    self._mark_variable(args[0])
                    self._show_page()
            elif action == "next":
                self._next_page()
            elif action == "back":
                self._prev_page()
            elif action == "jump" and args:
                self._jump_to(args[0])
            elif action == "preview":
                self._show_preview()
            elif action == "save":
                return self._save()
            elif action == "quit":
                if self._has_changes():
                    if not confirm_with_context(None, "Discard unsaved changes?", default=False):
                        continue
                return None
            elif action == "help":
                self._show_help()

    # ── Display ────────────────────────────────────────────────

    def _show_header(self):
        """Display editor header with session info and change summary."""
        filename = os.path.basename(self.session_file)
        info_parts = [f"File: {filename}"]

        if self.metadata.get("pdb_id"):
            info_parts.append(f"PDB: {self.metadata['pdb_id']}")
        elif self.metadata.get("pdb_file"):
            info_parts.append(f"File: {os.path.basename(self.metadata['pdb_file'])}")

        info_parts.append(f"{len(self.display_items)} interactions")

        # Change summary
        change_parts = []
        if self.edits:
            change_parts.append(f"{len(self.edits)} edit{'s' if len(self.edits) != 1 else ''}")
        if self.deletes:
            change_parts.append(f"{len(self.deletes)} delete{'s' if len(self.deletes) != 1 else ''}")
        if self.variables:
            change_parts.append(f"{len(self.variables)} variable{'s' if len(self.variables) != 1 else ''}")

        if change_parts:
            info_parts.append(", ".join(change_parts))
        else:
            info_parts.append("no changes")

        mode_label = "Template Creator" if self.mode == "template" else "Session Editor"
        self.console.print(Panel(
            "  |  ".join(info_parts),
            title=f"[bold]{mode_label}[/bold]",
            border_style="bright_blue",
            expand=False,
        ))

    def _show_page(self):
        """Display current page of interactions."""
        start = self.current_page * self.PAGE_SIZE
        end = min(start + self.PAGE_SIZE, len(self.display_items))

        table = Table(show_header=True, header_style="bold magenta", expand=True)
        table.add_column("#", style="bold bright_blue", width=5, justify="right")
        table.add_column("St", width=2, justify="center")
        table.add_column("Module", style="#0f7f99", width=22, overflow="ellipsis")
        table.add_column("Question", overflow="fold", no_wrap=False)
        table.add_column("Answer", width=30, overflow="ellipsis")

        for i in range(start, end):
            disp_idx, actual_idx, interaction = self.display_items[i]
            context = interaction.get("context", {})
            response = interaction.get("response") or ""

            # Status marker and styling
            if disp_idx in self.variables:
                status = Text("V", style="bold #0f7f99")
                var_name = self.variables[disp_idx]["name"]
                answer = Text(f"{{{{ {var_name} }}}}", style="#0f7f99")
            elif disp_idx in self.deletes:
                status = Text("D", style="bold red")
                answer = Text("[deleted]", style="red")
            elif disp_idx in self.edits:
                status = Text("E", style="bold yellow")
                new_val = self.edits[disp_idx]
                # Show new value with old value reference
                option_label = self._resolve_option_label(interaction, new_val)
                old_label = self._resolve_option_label(interaction, response)
                if option_label != new_val:
                    answer_str = f"{option_label} (was: {old_label})"
                else:
                    answer_str = f"{new_val} (was: {response})"
                answer = Text(answer_str, style="yellow")
            else:
                status = Text("·", style="grey50")
                option_label = context.get("option_label") or self._resolve_option_label(interaction, response)
                if interaction.get("type") == "confirm":
                    option_label = "Yes" if response.lower() in ("yes", "y", "true", "1") else "No"
                answer = Text(option_label)

            # Module and question
            module = context.get("module", "")
            question = context.get("description") or interaction.get("prompt", "")

            table.add_row(str(disp_idx + 1), status, module, question, answer)

        self.console.print(table)
        self.console.print(f"[grey50]Page {self.current_page + 1}/{self.total_pages}[/grey50]")

        # Command hint. highlight=False disables Rich's automatic repr
        # highlighter, which otherwise bolds the [ ] brackets (repr.brace
        # style) — that bold renders invisibly on a white background, hiding
        # the shortcut brackets. Plain text keeps them in the default
        # foreground, legible on both light and dark backgrounds.
        if self.mode == "template":
            self.console.print(
                "\\[e]dit #  \\[d]elete #|#-#|#-end  \\[u]ndo #  \\[v]ariable #  "
                "\\[n]ext  \\[b]ack  \\[j]ump #  \\[p]review  save  quit  help",
                highlight=False,
            )
        else:
            self.console.print(
                "\\[e]dit #  \\[d]elete #|#-#|#-end  \\[u]ndo #|#-#  "
                "\\[n]ext  \\[b]ack  \\[j]ump #  \\[p]review  save  quit  help",
                highlight=False,
            )

    def _get_command(self) -> Optional[Tuple]:
        """Parse user command input. Returns (action, *args) or None."""
        try:
            raw = prompt_with_context(None, "> ")
        except (EOFError, KeyboardInterrupt):
            return ("quit",)

        raw = raw.strip().lower()
        if not raw:
            return None

        parts = raw.split()
        cmd = parts[0]

        # Single-letter shortcuts
        shortcuts = {
            "e": "edit", "d": "delete", "u": "undo", "v": "variable",
            "n": "next", "b": "back", "p": "preview", "j": "jump",
        }
        cmd = shortcuts.get(cmd, cmd)

        # Commands that take a number argument (or range for delete/undo).
        # Indices are shown to the user 1-based but stored 0-based, so this is
        # the single place we translate input: subtract 1.
        if cmd in ("edit", "delete", "undo", "variable", "jump"):
            if len(parts) < 2:
                self.console.print(f"[yellow]Usage: {cmd} <number>[/yellow]")
                return None

            n = len(self.display_items)

            # Range syntax for delete/undo: "delete 175-554". The end may be
            # "end"/"last" (or omitted, "175-") to mean "through the last".
            if cmd in ("delete", "undo") and "-" in parts[1]:
                range_parts = parts[1].split("-", 1)
                end_token = range_parts[1].strip()
                try:
                    start = int(range_parts[0]) - 1
                    if end_token in ("", "end", "last"):
                        end = n - 1
                    else:
                        end = int(end_token) - 1
                except ValueError:
                    self.console.print("[yellow]Please enter a valid range (e.g. 175-554 or 175-end)[/yellow]")
                    return None
                if start > end:
                    start, end = end, start
                if start < 0 or end > n - 1:
                    self.console.print(f"[yellow]Invalid range. Valid indices: 1-{n}[/yellow]")
                    return None
                return (cmd, start, end)

            try:
                idx = int(parts[1]) - 1
            except ValueError:
                self.console.print("[yellow]Please enter a valid number[/yellow]")
                return None
            if cmd != "jump" and not self._valid_display_idx(idx):
                self.console.print(f"[yellow]Invalid index. Valid range: 1-{n}[/yellow]")
                return None
            return (cmd, idx)

        if cmd in ("next", "back", "preview", "save", "quit", "help"):
            return (cmd,)

        self.console.print(f"[yellow]Unknown command: {raw}. Type 'help' for commands.[/yellow]")
        return None

    # ── Actions ────────────────────────────────────────────────

    def _edit_interaction(self, display_idx: int):
        """Prompt user for a new response value at the given interaction."""
        disp_idx, actual_idx, interaction = self.display_items[display_idx]
        context = interaction.get("context", {})
        response = interaction.get("response") or ""
        description = context.get("description") or interaction.get("prompt", "")

        self.console.print(f"\n[bold]Editing #{display_idx + 1}:[/bold] {description}")

        # Show available options if this was a choice prompt
        options_map = context.get("options_map", {})
        choices = interaction.get("choices")
        if options_map:
            self.console.print("[grey50]Available options:[/grey50]")
            for key, label in options_map.items():
                marker = " [yellow]<-- current[/yellow]" if key == response else ""
                self.console.print(f"  [bright_blue]{key}[/bright_blue]: {label}{marker}")

        current_label = self._resolve_option_label(interaction, response)
        self.console.print(f"[grey50]Current value: {response} ({current_label})[/grey50]")

        new_value = prompt_with_context(None, "New value", default=response)

        if new_value == response:
            self.console.print("[grey50]No change[/grey50]")
            return

        new_label = self._resolve_option_label(interaction, new_value)
        self.edits[display_idx] = new_value

        # Remove from deletes/variables if it was there
        self.deletes.discard(display_idx)
        self.variables.pop(display_idx, None)

        self.console.print(f"[green]Marked #{display_idx + 1} for edit: {current_label} -> {new_label}[/green]")

    def _delete_interaction(self, display_idx: int):
        """Mark an interaction for deletion."""
        disp_idx, actual_idx, interaction = self.display_items[display_idx]
        context = interaction.get("context", {})
        description = context.get("description") or interaction.get("prompt", "")

        self.deletes.add(display_idx)
        # Remove from edits/variables if it was there
        self.edits.pop(display_idx, None)
        self.variables.pop(display_idx, None)

        self.console.print(f"[red]Deleted #{display_idx + 1}: {description}[/red] [grey50](use 'u {display_idx + 1}' to restore)[/grey50]")

    def _delete_range(self, start: int, end: int):
        """Mark a range of interactions for deletion (inclusive)."""
        count = 0
        for idx in range(start, end + 1):
            self.deletes.add(idx)
            self.edits.pop(idx, None)
            self.variables.pop(idx, None)
            count += 1
        self.console.print(
            f"[red]Deleted #{start + 1}-#{end + 1} ({count} interactions)[/red] "
            f"[grey50](use 'u {start + 1}-{end + 1}' to restore)[/grey50]"
        )

    def _undo_change(self, display_idx: int):
        """Remove any pending change at the given interaction."""
        changed = False
        if display_idx in self.edits:
            del self.edits[display_idx]
            changed = True
        if display_idx in self.deletes:
            self.deletes.discard(display_idx)
            changed = True
        if display_idx in self.variables:
            del self.variables[display_idx]
            changed = True

        if changed:
            self.console.print(f"[green]Restored #{display_idx + 1} to original[/green]")
        else:
            self.console.print(f"[grey50]#{display_idx + 1} has no pending changes[/grey50]")

    def _undo_range(self, start: int, end: int):
        """Undo pending changes for a range of interactions (inclusive)."""
        count = 0
        for idx in range(start, end + 1):
            if idx in self.edits or idx in self.deletes or idx in self.variables:
                self.edits.pop(idx, None)
                self.deletes.discard(idx)
                self.variables.pop(idx, None)
                count += 1
        if count:
            self.console.print(f"[green]Restored #{start + 1}-#{end + 1} ({count} interactions) to original[/green]")
        else:
            self.console.print(f"[grey50]No pending changes in range #{start + 1}-#{end + 1}[/grey50]")

    def _mark_variable(self, display_idx: int):
        """Mark an interaction as a template variable."""
        disp_idx, actual_idx, interaction = self.display_items[display_idx]
        context = interaction.get("context", {})
        response = interaction.get("response") or ""
        description = context.get("description") or interaction.get("prompt", "")

        self.console.print(f"\n[bold cyan]Marking #{display_idx + 1} as template variable[/bold cyan]")
        self.console.print(f"[grey50]Interaction: {description}[/grey50]")

        # Show available options for context
        options_map = context.get("options_map", {})
        choices = interaction.get("choices")
        if options_map:
            self.console.print("[grey50]This interaction had these options:[/grey50]")
            for key, label in options_map.items():
                marker = " [yellow]<-- original[/yellow]" if key == response else ""
                self.console.print(f"  [bright_blue]{key}[/bright_blue]: {label}{marker}")

        current_label = self._resolve_option_label(interaction, response)
        self.console.print(f"[grey50]Original value: {response} ({current_label})[/grey50]")

        # Suggest a variable name from description
        suggested_name = self._suggest_variable_name(description, display_idx)
        var_name = prompt_with_context(None, "Variable name", default=suggested_name)

        # Check for duplicate variable names
        for other_idx, other_var in self.variables.items():
            if other_var["name"] == var_name and other_idx != display_idx:
                self.console.print(f"[red]Variable name '{var_name}' already used at #{other_idx + 1}[/red]")
                return

        var_description = prompt_with_context(None, "Description", default=description)
        required = confirm_with_context(None, "Required?", default=True)

        self.variables[display_idx] = {
            "name": var_name,
            "description": var_description,
            "required": required,
            "choices": choices,
            "options_map": options_map,
            "original_value": response,
            "original_label": current_label,
        }

        # Remove from edits/deletes if it was there
        self.edits.pop(display_idx, None)
        self.deletes.discard(display_idx)

        self.console.print(f"[cyan]Marked as variable: {{{{ {var_name} }}}}[/cyan]")

    # ── Navigation ─────────────────────────────────────────────

    def _next_page(self):
        if self.current_page < self.total_pages - 1:
            self.current_page += 1
            self._show_page()
        else:
            self.console.print("[grey50]Already on the last page[/grey50]")

    def _prev_page(self):
        if self.current_page > 0:
            self.current_page -= 1
            self._show_page()
        else:
            self.console.print("[grey50]Already on the first page[/grey50]")

    def _jump_to(self, display_idx: int):
        if not self._valid_display_idx(display_idx):
            self.console.print(f"[yellow]Invalid index. Valid range: 1-{len(self.display_items)}[/yellow]")
            return
        self.current_page = display_idx // self.PAGE_SIZE
        self._show_page()

    # ── Preview & Save ─────────────────────────────────────────

    def _show_preview(self):
        """Show summary of all pending changes."""
        if not self._has_changes():
            self.console.print("[grey50]No changes to preview[/grey50]")
            return

        table = Table(title="Pending Changes", show_header=True, header_style="bold magenta")
        table.add_column("#", style="bold bright_blue", width=5, justify="right")
        table.add_column("Action", width=10)
        table.add_column("Description", width=40)
        table.add_column("Detail", width=30)

        for disp_idx in sorted(set(list(self.edits.keys()) + list(self.deletes) + list(self.variables.keys()))):
            _, _, interaction = self.display_items[disp_idx]
            context = interaction.get("context", {})
            description = context.get("description") or interaction.get("prompt", "")
            response = interaction.get("response") or ""

            if disp_idx in self.variables:
                var = self.variables[disp_idx]
                table.add_row(
                    str(disp_idx + 1),
                    Text("VARIABLE", style="#0f7f99"),
                    description,
                    Text(f"{{{{ {var['name']} }}}}", style="#0f7f99"),
                )
            elif disp_idx in self.deletes:
                old_label = self._resolve_option_label(interaction, response)
                table.add_row(
                    str(disp_idx + 1),
                    Text("DELETE", style="red"),
                    description,
                    Text(old_label, style="red"),
                )
            elif disp_idx in self.edits:
                old_label = self._resolve_option_label(interaction, response)
                new_label = self._resolve_option_label(interaction, self.edits[disp_idx])
                table.add_row(
                    str(disp_idx + 1),
                    Text("EDIT", style="yellow"),
                    description,
                    Text(f"{old_label} -> {new_label}", style="yellow"),
                )

        self.console.print(table)

    def _save(self) -> Optional[str]:
        """Apply changes and save. Returns path to saved file."""
        if not self._has_changes() and self.mode == "edit":
            self.console.print("[yellow]No changes to save[/yellow]")
            if not confirm_with_context(None, "Save unchanged copy anyway?", default=False):
                return None

        # Warn about edits that may cause downstream divergence
        if self.edits or self.deletes:
            self.console.print(
                "\n[yellow]Reminder: If any edited or deleted interaction changes the "
                "workflow path (e.g., a different menu selection), later interactions "
                "may not match during replay. ProPrep will switch to live input at "
                "that point.[/yellow]"
            )
            if not confirm_with_context(None, "Continue with save?", default=True):
                return None

        if self.mode == "template":
            if not self.variables:
                self.console.print("[yellow]No variables marked. A template needs at least one variable.[/yellow]")
                if not confirm_with_context(None, "Save as edited session instead?", default=True):
                    return None
                return self._save_session()
            return self._save_template()
        else:
            return self._save_session()

    def _save_session(self) -> str:
        """Apply edits and save as modified session file."""
        modified_data = deepcopy(self.session_data)
        all_interactions = modified_data.get("interactions", [])

        # Build the set of actual indices to delete
        actual_deletes = set()
        for disp_idx in self.deletes:
            _, actual_idx, _ = self.display_items[disp_idx]
            actual_deletes.add(actual_idx)

        # Apply edits
        for disp_idx, new_value in self.edits.items():
            _, actual_idx, _ = self.display_items[disp_idx]
            if actual_idx < len(all_interactions):
                old_value = all_interactions[actual_idx].get("response", "")
                all_interactions[actual_idx]["response"] = new_value
                # Update option_label in context if options_map exists
                context = all_interactions[actual_idx].get("context", {})
                options_map = context.get("options_map", {})
                if new_value in options_map:
                    context["option_label"] = options_map[new_value]

        # Remove deleted interactions and re-index
        filtered = []
        for i, interaction in enumerate(all_interactions):
            if i not in actual_deletes:
                interaction["index"] = len(filtered)
                filtered.append(interaction)

        modified_data["interactions"] = filtered

        # Add edit metadata
        modified_data.setdefault("metadata", {})
        modified_data["metadata"]["edited_from"] = os.path.basename(self.session_file)
        modified_data["metadata"]["edit_time"] = datetime.now().isoformat()
        modified_data["metadata"]["changes_applied"] = {
            "edits": len(self.edits),
            "deletes": len(self.deletes),
        }

        # Generate output filename
        stem = Path(self.session_file).stem
        parent = Path(self.session_file).parent
        output_file = str(parent / f"{stem}_edited.json")

        with open(output_file, "w") as f:
            json.dump(modified_data, f, indent=2)

        self.console.print(f"\n[bold green]Session saved: {output_file}[/bold green]")
        self.console.print(f"[grey50]{len(self.edits)} edit(s), {len(self.deletes)} deletion(s) applied[/grey50]")
        return output_file

    def _save_template(self) -> str:
        """Apply edits + variables and save as template file."""
        modified_data = deepcopy(self.session_data)
        all_interactions = modified_data.get("interactions", [])

        # Build the set of actual indices to delete
        actual_deletes = set()
        for disp_idx in self.deletes:
            _, actual_idx, _ = self.display_items[disp_idx]
            actual_deletes.add(actual_idx)

        # Apply edits
        for disp_idx, new_value in self.edits.items():
            _, actual_idx, _ = self.display_items[disp_idx]
            if actual_idx < len(all_interactions):
                all_interactions[actual_idx]["response"] = new_value
                context = all_interactions[actual_idx].get("context", {})
                options_map = context.get("options_map", {})
                if new_value in options_map:
                    context["option_label"] = options_map[new_value]

        # Apply variable substitutions
        template_variables = {}
        for disp_idx, var_info in self.variables.items():
            _, actual_idx, interaction = self.display_items[disp_idx]
            if actual_idx < len(all_interactions):
                var_name = var_info["name"]
                all_interactions[actual_idx]["response"] = f"{{{{ {var_name} }}}}"

                # Add template marker to context
                ctx = all_interactions[actual_idx].setdefault("context", {})
                ctx["template_variable"] = var_name

                # Build template variable definition
                var_def = {
                    "description": var_info["description"],
                    "type": "choice" if var_info.get("choices") else "text",
                    "example": var_info["original_value"],
                    "required": var_info["required"],
                }
                # Include choices/options_map for context-aware prompts
                if var_info.get("choices"):
                    var_def["choices"] = var_info["choices"]
                if var_info.get("options_map"):
                    var_def["options_map"] = var_info["options_map"]
                if var_info.get("original_label"):
                    var_def["original_label"] = var_info["original_label"]

                template_variables[var_name] = var_def

        # Remove deleted interactions and re-index
        filtered = []
        for i, interaction in enumerate(all_interactions):
            if i not in actual_deletes:
                interaction["index"] = len(filtered)
                filtered.append(interaction)

        modified_data["interactions"] = filtered

        # Set template fields
        modified_data["version"] = "1.2"
        modified_data["template"] = True
        modified_data["template_variables"] = template_variables
        modified_data["template_metadata"] = {
            "created_at": datetime.now().isoformat(),
            "created_from": os.path.basename(self.session_file),
            "description": "ProPrep batch processing template",
            "proprep_version": self.metadata.get("proprep_version", "unknown"),
        }

        if self.edits:
            modified_data.setdefault("metadata", {})
            modified_data["metadata"]["edits_applied"] = len(self.edits)

        # Generate output filename
        stem = Path(self.session_file).stem
        parent = Path(self.session_file).parent
        output_file = str(parent / f"{stem}_template.json")

        with open(output_file, "w") as f:
            json.dump(modified_data, f, indent=2)

        self.console.print(f"\n[bold green]Template saved: {output_file}[/bold green]")
        self.console.print(f"[grey50]{len(template_variables)} variable(s), "
                           f"{len(self.edits)} edit(s), {len(self.deletes)} deletion(s)[/grey50]")

        # Show variables summary
        self.console.print("\n[bold]Template variables:[/bold]")
        for var_name, var_def in template_variables.items():
            req = "[red]*[/red]" if var_def["required"] else ""
            self.console.print(f"  {req}[cyan]{var_name}[/cyan]: {var_def['description']} "
                               f"[grey50](example: {var_def['example']})[/grey50]")

        # Offer to generate batch input file
        self._generate_batch_file(output_file, template_variables)

        return output_file

    # ── Batch file generation ──────────────────────────────────

    def _generate_batch_file(self, template_file: str, template_variables: Dict[str, dict]):
        """Prompt user to create a CSV batch input file after template creation."""
        self.console.print()
        if not confirm_with_context(None, "Create a batch input file now?", default=True):
            self.console.print(
                f"[grey50]You can create one later. Format: CSV with headers matching variable names.[/grey50]"
            )
            return

        var_names = list(template_variables.keys())
        rows = []
        run_num = 1

        self.console.print(
            "\n[bold]Enter variable values for each run.[/bold] "
            "[grey50]Press Enter with no input to finish.[/grey50]\n"
        )

        while True:
            self.console.print(f"[bold bright_blue]Run {run_num}:[/bold bright_blue]")
            row = {}
            cancelled = False

            for var_name in var_names:
                var_def = template_variables[var_name]

                # Show available options if this was a choice variable
                options_map = var_def.get("options_map", {})
                if options_map:
                    self.console.print(f"  [grey50]{var_name} options:[/grey50]")
                    for key, label in options_map.items():
                        self.console.print(f"    [bright_blue]{key}[/bright_blue]: {label}")

                example = var_def.get("example", "")
                value = prompt_with_context(
                    None,
                    f"  {var_name}",
                    default="" if run_num > 1 else example,
                )

                if not value and run_num > 1:
                    # Empty input on non-first run = done
                    cancelled = True
                    break

                if not value and var_def.get("required", True):
                    self.console.print(f"  [yellow]{var_name} is required[/yellow]")
                    cancelled = True
                    break

                row[var_name] = value

            if cancelled or not row:
                break

            rows.append(row)
            run_num += 1
            self.console.print()

        if not rows:
            self.console.print("[grey50]No runs entered. Batch file not created.[/grey50]")
            return

        # Write CSV
        stem = Path(template_file).stem.replace("_template", "")
        parent = Path(template_file).parent
        csv_file = str(parent / f"{stem}_batch.csv")

        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=var_names)
            writer.writeheader()
            writer.writerows(rows)

        self.console.print(f"\n[bold green]Batch file saved: {csv_file}[/bold green] ({len(rows)} run{'s' if len(rows) != 1 else ''})")
        # The batch runner reads the template's recorded mode automatically, but
        # surface it in the command too so it's explicit (and works even if the
        # template metadata is ever stripped). A web-recorded template needs the
        # matching replay mode or mode-gated prompts (e.g. the PDB Filter "view
        # structure?" Y/N) desync the replay.
        web_recorded = bool(self.metadata.get("web_shell_mode", False))
        web_flag = " --web" if web_recorded else ""
        self.console.print(f"[grey50]Run with: proprep --batch-replay {os.path.basename(template_file)} "
                           f"--batch-list {os.path.basename(csv_file)}{web_flag}[/grey50]")
        if web_recorded:
            self.console.print("[grey50]  (--web: this template was recorded under proprep-web; "
                               "the batch replays in browser mode to stay in sync.)[/grey50]")

    # ── Auto-detection ─────────────────────────────────────────

    def _auto_detect_variables(self):
        """Auto-detect template variables using TemplateConverter logic."""
        try:
            converter = TemplateConverter(self.console)
            detected = converter._detect_variables(self.session_data)
        except Exception:
            return

        if not detected:
            return

        # Map detected variables to display indices
        for var_name, var_info in detected.items():
            interaction_idx = var_info.get("interaction_index")
            if interaction_idx is None:
                continue

            # Find the display index for this actual index
            for disp_idx, actual_idx, interaction in self.display_items:
                if actual_idx == interaction_idx:
                    context = interaction.get("context", {})
                    self.variables[disp_idx] = {
                        "name": var_name,
                        "description": var_info.get("description", ""),
                        "required": var_info.get("required", True),
                        "choices": interaction.get("choices"),
                        "options_map": context.get("options_map", {}),
                        "original_value": var_info.get("original_value", ""),
                        "original_label": var_info.get("original_value", ""),
                    }
                    break

        if self.variables:
            self.console.print(f"[cyan]Auto-detected {len(self.variables)} variable(s). "
                               f"Use 'u #' to unmark or 'v #' to mark more.[/cyan]")

    # ── Helpers ────────────────────────────────────────────────

    def _valid_display_idx(self, idx: int) -> bool:
        return 0 <= idx < len(self.display_items)

    def _has_changes(self) -> bool:
        return bool(self.edits or self.deletes or self.variables)

    def _resolve_option_label(self, interaction: dict, value: str) -> str:
        """Resolve a response value to its display label using options_map."""
        context = interaction.get("context", {})
        options_map = context.get("options_map", {})
        if value in options_map:
            return options_map[value]
        option_label = context.get("option_label", "")
        if option_label and value == interaction.get("response"):
            return option_label
        return value or ""

    def _suggest_variable_name(self, description: str, display_idx: int) -> str:
        """Suggest a variable name from the interaction description."""
        # Common patterns
        desc_lower = description.lower()
        if "pdb" in desc_lower and ("id" in desc_lower or "download" in desc_lower):
            return "input_protein"
        if "uniprot" in desc_lower or "alphafold" in desc_lower:
            return "input_protein"
        if "file" in desc_lower and "path" in desc_lower:
            return "input_file"
        if "force field" in desc_lower or "forcefield" in desc_lower:
            return "forcefield"
        if "water" in desc_lower and "model" in desc_lower:
            return "water_model"

        # Fallback: clean up description
        clean = description.lower().replace(" ", "_")
        clean = "".join(c for c in clean if c.isalnum() or c == "_")
        if clean and len(clean) <= 30:
            return clean
        return f"var_{display_idx + 1}"

    def _show_help(self):
        """Display help text for editor commands."""
        help_text = [
            "[bold]Commands:[/bold]",
            "",
            "  [bright_blue]e <#>[/bright_blue]      Edit interaction response",
            "  [bright_blue]d <#>[/bright_blue]      Delete interaction (skip during replay)",
            "  [bright_blue]d <#>-<#>[/bright_blue]  Delete range (e.g. d 175-554)",
            "  [bright_blue]u <#>[/bright_blue]      Undo change (restore to original)",
            "  [bright_blue]u <#>-<#>[/bright_blue]  Undo range",
        ]
        if self.mode == "template":
            help_text.append("  [bright_blue]v <#>[/bright_blue]  Mark as template variable")
        help_text.extend([
            "",
            "  [bright_blue]n[/bright_blue]      Next page",
            "  [bright_blue]b[/bright_blue]      Back (previous page)",
            "  [bright_blue]j <#>[/bright_blue]  Jump to page containing interaction #",
            "",
            "  [bright_blue]p[/bright_blue]      Preview all pending changes",
            "  [bright_blue]save[/bright_blue]     Apply changes and save",
            "  [bright_blue]quit[/bright_blue]     Discard and exit",
            "  [bright_blue]help[/bright_blue]     Show this help",
        ])
        self.console.print(Panel("\n".join(help_text), title="Help", border_style="bright_blue", expand=False))
