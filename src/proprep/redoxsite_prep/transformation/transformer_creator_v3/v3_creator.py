"""
Transformer Creator v3 — Main Orchestrator

Guided transformer creation through:
  1. Automatic site requirements (from RedoxSite)
  2. Role labeling (user names each residue)
  3. Relationship phrases (constraint-based resolution)
  4. PDB editing commands (sequential passes)
  5. Parameter definition (state-dependent name grid)
  6. Metadata + code generation
"""

import logging
import readline
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .command_parser import CommandParser, CommandParseError, CommandValidationError
from .constraint_resolver import ConstraintResolver, ResolutionResult
from .data_models import (
    EditCommand,
    EditingPass,
    ParameterOption,
    ParameterSpec,
    RoleDefinition,
    TransformerSpecV3,
)
from .phrase_parser import PhraseParser, ParseError, ValidationError, BOND_TYPES

logger = logging.getLogger(__name__)


class _PhraseCompleter:
    """Readline completer for relationship phrase input."""

    def __init__(self, roles: Dict[str, RoleDefinition]):
        self.roles = roles
        self._matches: List[str] = []

    def complete(self, text: str, state: int) -> Optional[str]:
        if state == 0:
            self._matches = self._get_matches(text)
        if state < len(self._matches):
            return self._matches[state]
        return None

    def _get_matches(self, text: str) -> List[str]:
        line = readline.get_line_buffer()

        # After "with " → suggest bond types
        if " with " in line and line.rstrip().endswith(text):
            return [bt + " " for bt in sorted(BOND_TYPES) if bt.startswith(text)]

        # After "role:" or "of " context → suggest atom names
        # Find if we're in an "atom on role" or "role:atom" context
        for label, role in self.roles.items():
            # After "[role]:" in command context
            if f"{label}:" in line and line.rstrip().endswith(text):
                prefix = text.upper()
                return [a for a in sorted(role.atom_names) if a.startswith(prefix)]
            # After "[atom] of " or "[atom] on " — suggest role names
            if line.rstrip().endswith("of " + text) or line.rstrip().endswith("on " + text):
                return [r + " " for r in sorted(self.roles.keys()) if r.startswith(text)]

        # Default: suggest role names and keywords
        matches = []
        for label in sorted(self.roles.keys()):
            if label.startswith(text):
                matches.append(label + " ")

        keywords = ["bonds", "through", "to", "on", "with",
                     "is", "closest", "farthest", "within", "from",
                     "has", "a", "lower", "higher", "resid", "than",
                     "the", "same", "different", "chain", "as"]
        for kw in keywords:
            if kw.startswith(text) and kw + " " not in matches:
                matches.append(kw + " ")

        return matches


class _CommandCompleter:
    """Readline completer for PDB editing command input."""

    def __init__(self, cmd_parser: CommandParser, roles: Dict[str, RoleDefinition]):
        self.cmd_parser = cmd_parser
        self.roles = roles
        self._matches: List[str] = []

    def complete(self, text: str, state: int) -> Optional[str]:
        if state == 0:
            self._matches = self._get_matches(text)
        if state < len(self._matches):
            return self._matches[state]
        return None

    def _get_matches(self, text: str) -> List[str]:
        line = readline.get_line_buffer()

        # At start of line → command verbs
        if not line.strip() or line.strip() == text:
            verbs = ["rename ", "move ", "hetatm ", "atom "]
            return [v for v in verbs if v.startswith(text)]

        tokens = line.split()

        # After verb → suggest role names (possibly with :atoms)
        if len(tokens) == 1 and line.endswith(" " + text if text else " "):
            return [r + " " for r in sorted(self.cmd_parser._role_atoms.keys())
                    if r.startswith(text)]

        # After "role:" → suggest atom names
        for label in self.cmd_parser._role_atoms:
            if f"{label}:" in line and line.rstrip().endswith(text):
                atoms = sorted(self.cmd_parser.get_role_atoms(label))
                prefix = text.upper()
                return [a for a in atoms if a.startswith(prefix)]

        # After "to " → suggest target roles or "new"
        if " to " in line and line.rstrip().endswith(text):
            targets = sorted(self.cmd_parser._role_atoms.keys()) + ["new "]
            return [t + " " if not t.endswith(" ") else t
                    for t in targets if t.startswith(text)]

        # After "as " → nothing specific, user types new label
        # After "new RESNAME " → suggest "as "
        if " to new " in line and "as" not in line.lower() and text == "":
            return ["as "]

        # Default: role names
        return [r + " " for r in sorted(self.cmd_parser._role_atoms.keys())
                if r.startswith(text)]


class TransformerCreatorV3:
    """
    Main orchestrator for creating transformers via the v3 guided flow.

    Args:
        redox_site: a detected RedoxSite object
        output_dir: directory to write generated transformer (default: transformers/)
        console: Rich Console for output
    """

    def __init__(self, redox_site, output_dir: Optional[Path] = None,
                 console: Optional[Console] = None):
        self.site = redox_site
        self.console = console or Console()
        self.spec = TransformerSpecV3()

        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            this_file = Path(__file__).resolve()
            self.output_dir = this_file.parent.parent / "transformers"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Working state
        self._roles: Dict[str, RoleDefinition] = {}
        self._resolved: Dict[str, Tuple[str, int]] = {}  # label -> (chain, resid)
        self._old_completer = None
        self._old_delims = None

    def _install_completer(self, completer):
        """Install a readline completer, saving the previous one."""
        self._old_completer = readline.get_completer()
        self._old_delims = readline.get_completer_delims()
        readline.set_completer(completer.complete)
        readline.set_completer_delims(' \t\n')
        readline.parse_and_bind('tab: complete')

    def _restore_completer(self):
        """Restore the previous readline completer."""
        readline.set_completer(self._old_completer)
        if self._old_delims is not None:
            readline.set_completer_delims(self._old_delims)

    def create_transformer(self) -> Optional[Path]:
        """
        Run the full guided creation flow.

        Returns path to generated transformer file, or None if cancelled.
        """
        self.console.print(Panel(
            "[bold]Transformer Creator v3[/bold]\n"
            "Create a transformer by labeling residues, defining relationships,\n"
            "and demonstrating PDB edits.",
            title="Welcome",
        ))

        # Step 1: Site requirements (automatic)
        self._step_site_requirements()

        # Step 2: Role labeling
        if not self._step_role_labeling():
            return None

        # Step 3: Relationship definition
        if not self._step_relationships():
            return None

        # Step 4: PDB editing
        if not self._step_pdb_editing():
            return None

        # Step 5: Parameters
        self._step_parameters()

        # Step 6: Metadata + generate
        return self._step_metadata_and_generate()

    # ===== STEP 1: SITE REQUIREMENTS =====

    def _step_site_requirements(self):
        """Extract site requirements automatically from the RedoxSite."""
        self.console.print("\n[bold]Step 1: Site Requirements[/bold] (automatic)")

        # Count residue types
        resname_counts: Dict[str, int] = {}
        for atom in self.site.atoms:
            key = (atom.chain, atom.resid, atom.resname)
            resname_counts[atom.resname] = resname_counts.get(atom.resname, 0)
        # Deduplicate by residue, not atom
        residues_seen = set()
        resname_counts = {}
        for atom in self.site.atoms:
            rkey = (atom.chain, atom.resid)
            if rkey not in residues_seen:
                residues_seen.add(rkey)
                resname_counts[atom.resname] = resname_counts.get(atom.resname, 0) + 1

        # Display
        self.console.print(f"  Centers: {len(self.site.centers)}")
        self.console.print(f"  Residues: {len(residues_seen)}")
        for resname, count in sorted(resname_counts.items()):
            self.console.print(f"    {resname}: {count}")
        self.console.print(f"  Bonds: {len(self.site.bonds)}")
        self.console.print("  [grey50]Requirements will be derived automatically.[/grey50]")

    # ===== STEP 2: ROLE LABELING =====

    def _step_role_labeling(self) -> bool:
        """User assigns role labels to each residue in the site."""
        self.console.print("\n[bold]Step 2: Role Labeling[/bold]")
        self.console.print("Assign a label to each residue in the site.\n")

        # Build residue list
        residues = self._get_residue_list()

        # Display residues
        table = Table(title="Residues in site")
        table.add_column("#", style="grey50")
        table.add_column("Chain")
        table.add_column("ResID")
        table.add_column("ResName")
        table.add_column("Atoms")
        table.add_column("Center?")

        for i, res in enumerate(residues, 1):
            is_center = "yes" if res["is_center"] else ""
            atoms_str = ", ".join(res["atom_names"][:8])
            if len(res["atom_names"]) > 8:
                atoms_str += f", ... ({len(res['atom_names'])} total)"
            table.add_row(
                str(i), res["chain"], str(res["resid"]),
                res["resname"], atoms_str, is_center,
            )

        self.console.print(table)
        self.console.print()

        # Collect labels
        for i, res in enumerate(residues, 1):
            default_label = ""
            if res["is_center"]:
                default_label = "center"

            prompt = f"  Label for {res['resname']} {res['chain']}:{res['resid']}"
            if default_label:
                prompt += f" [{default_label}]"
            prompt += ": "

            label = input(prompt).strip()
            if not label and default_label:
                label = default_label
            if not label:
                self.console.print("  [red]Label required. Cancelled.[/red]")
                return False

            # Validate snake_case
            if not label.replace("_", "").isalnum():
                self.console.print(f"  [red]Label must be snake_case: '{label}'[/red]")
                return False

            if label in self._roles:
                self.console.print(f"  [red]Duplicate label: '{label}'[/red]")
                return False

            role = RoleDefinition(
                label=label,
                resname=res["resname"],
                chain=res["chain"],
                resid=res["resid"],
                is_center=res["is_center"],
                atom_names=res["atom_names"],
            )

            # Ask about alternative resnames
            alt = input(f"    Could this be a different resname in other structures? "
                       f"[N/comma-separated names]: ").strip()
            if alt and alt.lower() != "n":
                role.alt_resnames = [a.strip().upper() for a in alt.split(",")]

            self._roles[label] = role

        self.console.print(f"\n  [green]Labeled {len(self._roles)} roles.[/green]")
        return True

    # ===== STEP 3: RELATIONSHIPS =====

    def _step_relationships(self) -> bool:
        """User defines relationship phrases, system resolves roles."""
        self.console.print("\n[bold]Step 3: Relationship Definition[/bold]")
        self.console.print("Define relationships between roles using phrases.")
        self.console.print("Categories:")
        self.console.print("  [cyan]Connectivity:[/cyan] X bonds [through ATOM] to [ATOM on] Y [with TYPE]")
        self.console.print("  [cyan]Spatial:[/cyan]       [ATOM of] X is closest to/farthest from/within N Å of [ATOM on] Y")
        self.console.print("  [cyan]Sequence:[/cyan]      X is +N/-N from Y")
        self.console.print("                   X has a lower/higher resid than Y")
        self.console.print("                   X is on the same/different chain as Y")
        self.console.print("\n  Type phrases (blank line when done):\n")

        parser = PhraseParser(self._roles, self.site.bonds)

        # Install tab completion for role names, atoms, bond types
        completer = _PhraseCompleter(self._roles)
        self._install_completer(completer)
        self.console.print("  [grey50]Tab completion available for role names, atoms, and bond types.[/grey50]\n")

        try:
            while True:
                text = input("  > ").strip()
                if not text:
                    break

                try:
                    phrase = parser.parse(text)
                    # Add phrase to the subject role
                    role = self._roles[phrase.subject_role]
                    role.phrases.append(phrase)

                    # Show immediate evaluation
                    resolver = ConstraintResolver(self._roles, self.site)
                    result = resolver.resolve()
                    self._display_resolution_status(result)

                except ParseError as e:
                    self.console.print(f"  [red]Parse error: {e}[/red]")
                except ValidationError as e:
                    self.console.print(f"  [red]Validation error: {e}[/red]")
        finally:
            self._restore_completer()

        # Final resolution
        resolver = ConstraintResolver(self._roles, self.site)
        result = resolver.resolve()

        if not result.is_complete:
            self.console.print("\n  [yellow]Not all roles resolved:[/yellow]")
            self._display_resolution_status(result)
            proceed = input("  Continue anyway? [y/N]: ").strip().lower()
            if proceed != "y":
                return False

        self._resolved = result.resolved
        self.console.print(f"\n  [green]All {len(result.resolved)} roles resolved.[/green]")
        return True

    def _display_resolution_status(self, result: ResolutionResult):
        """Display current resolution status."""
        for label, rkey in sorted(result.resolved.items()):
            role = self._roles[label]
            self.console.print(
                f"    [green]✓[/green] {label}: {role.resname} {rkey[0]}:{rkey[1]}"
            )
        for label, candidates in sorted(result.ambiguous.items()):
            role = self._roles[label]
            cand_str = ", ".join(f"{c[0]}:{c[1]}" for c in candidates)
            self.console.print(
                f"    [yellow]?[/yellow] {label}: {len(candidates)} candidates — {cand_str}"
            )
        for label in result.conflicting:
            self.console.print(
                f"    [red]✗[/red] {label}: no candidates (conflicting constraints)"
            )

    # ===== STEP 4: PDB EDITING =====

    def _step_pdb_editing(self) -> bool:
        """User enters editing commands in sequential passes."""
        self.console.print("\n[bold]Step 4: PDB Editing[/bold]")
        self.console.print("Enter editing commands in sequential passes.")
        self.console.print("Commands:")
        self.console.print("  [cyan]rename[/cyan] [role] [new_resname]")
        self.console.print("  [cyan]rename[/cyan] [role]:[old_atoms] [new_atoms]")
        self.console.print("  [cyan]move[/cyan]   [role]:[atoms] to [role]")
        self.console.print("  [cyan]move[/cyan]   [role]:[atoms] to new [resname] as [new_role]")
        self.console.print("  [cyan]hetatm[/cyan] [role]")
        self.console.print("  [cyan]atom[/cyan]   [role]")
        self.console.print("\n  Type commands. Blank line ends a pass. 'done' ends editing.\n")

        cmd_parser = CommandParser(self._roles)
        pass_num = 1

        # Install tab completion for command verbs, roles, atoms
        completer = _CommandCompleter(cmd_parser, self._roles)
        self._install_completer(completer)
        self.console.print("  [grey50]Tab completion available for commands, roles, and atoms.[/grey50]\n")

        try:
            while True:
                desc = input(f"  Pass {pass_num} description (or 'done'): ").strip()
                if desc.lower() == "done":
                    break

                current_pass = EditingPass(description=desc)

                # Show current state
                self._display_site_state(cmd_parser)

                while True:
                    text = input("    > ").strip()
                    if not text:
                        break

                    try:
                        cmd = cmd_parser.parse(text)
                        cmd_parser.apply(cmd)
                        current_pass.commands.append(cmd)
                        self.console.print(f"    [green]✓[/green] {text}")

                        # If move to new, register the new role
                        if cmd.command_type.value == "move_atoms_new":
                            # Calculate new resid: source role resid + offset
                            source_role = self._roles[cmd.source_role]
                            new_count = sum(
                                1 for p in self.spec.editing_passes
                                for c in p.commands
                                if c.command_type.value == "move_atoms_new"
                            ) + sum(
                                1 for c in current_pass.commands
                                if c.command_type.value == "move_atoms_new"
                            )
                            new_resid = source_role.resid + new_count

                            new_role = RoleDefinition(
                                label=cmd.new_role_label,
                                resname=cmd.new_residue_name,
                                chain=source_role.chain,
                                resid=new_resid,
                                atom_names=list(cmd.atom_names),
                            )
                            self._roles[cmd.new_role_label] = new_role

                            # Update completer with new role
                            completer.roles = self._roles

                        # Show updated state
                        self._display_site_state(cmd_parser)

                    except (CommandParseError, CommandValidationError) as e:
                        self.console.print(f"    [red]{e}[/red]")

                if current_pass.commands:
                    self.spec.editing_passes.append(current_pass)
                    self.console.print(
                        f"  [green]Pass {pass_num}: {len(current_pass.commands)} commands recorded.[/green]"
                    )
                    pass_num += 1
        finally:
            self._restore_completer()

        total_cmds = sum(len(p.commands) for p in self.spec.editing_passes)
        self.console.print(
            f"\n  [green]{len(self.spec.editing_passes)} passes, "
            f"{total_cmds} commands total.[/green]"
        )
        return total_cmds > 0

    def _display_site_state(self, cmd_parser: CommandParser):
        """Display current state of all roles as a Rich table."""
        table = Table(title="Current Site State", show_lines=True)
        table.add_column("Role", style="cyan")
        table.add_column("Atoms")

        for label in sorted(cmd_parser._role_atoms.keys()):
            atoms = sorted(cmd_parser.get_role_atoms(label))
            atoms_str = ", ".join(atoms) if atoms else "[grey50](empty)[/grey50]"
            table.add_row(label, atoms_str)

        self.console.print(table)

    # ===== STEP 5: PARAMETERS =====

    def _step_parameters(self):
        """Collect parameter definitions if any names are state-dependent."""
        self.console.print("\n[bold]Step 5: Parameters[/bold]")

        # Collect all new resnames from editing passes
        new_names: Dict[str, str] = {}  # role_label -> new resname
        for editing_pass in self.spec.editing_passes:
            for cmd in editing_pass.commands:
                if cmd.command_type.value == "rename_residue" and cmd.new_resname:
                    new_names[cmd.source_role] = cmd.new_resname
                elif cmd.command_type.value == "move_atoms_new" and cmd.new_residue_name:
                    new_names[cmd.new_role_label] = cmd.new_residue_name

        if not new_names:
            self.console.print("  No residue renames — no parameters needed.")
            return

        self.console.print("  New residue names from editing:")
        for role, name in sorted(new_names.items()):
            self.console.print(f"    {role}: {name}")

        vary = input("\n  Do any of these names vary by electronic state? [y/N]: ").strip()
        if vary.lower() != "y":
            return

        # Define parameter(s)
        param_spec = ParameterSpec()

        while True:
            pname = input("  Parameter name (e.g., redox_state) [blank to stop]: ").strip()
            if not pname:
                break

            pdesc = input(f"  Description for {pname}: ").strip()
            ptype = input(f"  Type — (c)hoice or (f)ixed? [c]: ").strip().lower()

            if ptype == "f":
                value = input(f"  Fixed value: ").strip()
                note = input(f"  Note (why fixed): ").strip()
                param_spec.parameters.append(ParameterOption(
                    name=pname, description=pdesc, param_type="fixed",
                    fixed_value=value, note=note or None,
                ))
            else:
                options_str = input(f"  Options (comma-separated): ").strip()
                options = [o.strip() for o in options_str.split(",")]
                default = input(f"  Default [{options[0]}]: ").strip() or options[0]
                param_spec.parameters.append(ParameterOption(
                    name=pname, description=pdesc, param_type="choice",
                    options=options, default=default,
                ))

        # Identify which roles are state-dependent
        self.console.print("\n  Which roles have state-dependent names?")
        for role, name in sorted(new_names.items()):
            dep = input(f"    {role} ({name}) — state-dependent? [y/N]: ").strip()
            if dep.lower() == "y":
                param_spec.state_dependent_roles.append(role)

        # Build the name grid
        if param_spec.state_dependent_roles:
            choice_params = [p for p in param_spec.parameters if p.param_type == "choice"]
            if choice_params:
                # Generate combinations
                from itertools import product
                all_options = [p.options for p in choice_params]
                combos = list(product(*all_options))

                # The demonstrated case is the first combo (default values)
                demonstrated_key = "_".join(p.default for p in choice_params)
                param_spec.name_mappings[demonstrated_key] = {
                    role: new_names[role] for role in param_spec.state_dependent_roles
                }

                self.console.print(f"\n  You demonstrated the '{demonstrated_key}' case.")
                self.console.print("  Fill in names for other combinations:\n")

                for combo in combos:
                    combo_key = "_".join(combo)
                    if combo_key == demonstrated_key:
                        continue

                    self.console.print(f"  Combination: {combo_key}")
                    mapping = {}
                    for role in param_spec.state_dependent_roles:
                        alt = input(f"    {role} name: ").strip().upper()
                        mapping[role] = alt
                    param_spec.name_mappings[combo_key] = mapping

        self.spec.parameters = param_spec

    # ===== STEP 6: METADATA + GENERATE =====

    def _step_metadata_and_generate(self) -> Optional[Path]:
        """Collect metadata and generate the transformer."""
        self.console.print("\n[bold]Step 6: Metadata[/bold]")

        self.spec.name = input("  Transformer name (snake_case): ").strip()
        if not self.spec.name:
            self.console.print("  [red]Name required. Cancelled.[/red]")
            return None
        self.spec.description = input("  Description: ").strip()
        ff_path = input("  Forcefield path (e.g., heme/bis_his_c_type) [none]: ").strip()
        if ff_path and ff_path.lower() != "none":
            self.spec.forcefield_path = ff_path

        self.spec.roles = list(self._roles.values())

        # Calculate residue space from new roles created during editing
        new_roles = [r for r in self._roles.values() if not r.is_center
                     and r.label not in {rl.label for rl in self.spec.roles
                                         if any(a.chain == r.chain for a in self.site.atoms)}]
        # Count new residues created via move...to new
        new_residue_count = sum(
            1 for p in self.spec.editing_passes
            for c in p.commands
            if c.command_type.value == "move_atoms_new"
        )
        self.spec.required_residue_count = 1 + new_residue_count

        # Summary
        self.console.print("\n[bold]Summary:[/bold]")
        self.console.print(f"  Name: {self.spec.name}")
        self.console.print(f"  Description: {self.spec.description}")
        self.console.print(f"  Roles: {len(self.spec.roles)}")
        self.console.print(f"  Editing passes: {len(self.spec.editing_passes)}")
        self.console.print(f"  Residue IDs needed: {self.spec.required_residue_count}")
        if self.spec.parameters:
            self.console.print(f"  Parameters: {len(self.spec.parameters.parameters)}")

        confirm = input("\n  Generate transformer? [Y/n]: ").strip().lower()
        if confirm == "n":
            self.console.print("  Cancelled.")
            return None

        # Serialize to a JSON transformer spec — reusable DATA, not generated
        # code — written to the user transformer dir where SpecTransformer
        # discovers it. Falls back to the legacy code generator only if
        # serialization fails, so an authoring session is never lost.
        import json as _json
        from proprep.redoxsite_prep.transformation.auto_rename import (
            DEFAULT_USER_TRANSFORMER_DIR,
        )
        try:
            from proprep.redoxsite_prep.transformation.spec_serializer import (
                build_json_spec,
            )
            json_spec = build_json_spec(self.spec, self._roles, self._resolved, self.site)
            DEFAULT_USER_TRANSFORMER_DIR.mkdir(parents=True, exist_ok=True)
            output_file = DEFAULT_USER_TRANSFORMER_DIR / f"{self.spec.name}.json"
            with open(output_file, "w") as f:
                _json.dump(json_spec, f, indent=2)
            self.console.print(f"\n  [green]✓ Transformer spec saved: {output_file}[/green]")
            self.console.print("  [dim]Reusable on any protein with this site type; "
                               "discovered automatically by the Redox Site Preparer.[/dim]")
            return output_file
        except Exception as exc:
            logger.warning("JSON spec serialization failed (%s); falling back to "
                           "legacy code generation.", exc)
            from .code_generator import CodeGeneratorV3
            code = CodeGeneratorV3(self.spec, self._roles, self._resolved).generate()
            output_file = self.output_dir / f"{self.spec.name}.py"
            with open(output_file, "w") as f:
                f.write(code)
            self.console.print(f"\n  [yellow]✓ Transformer generated (legacy .py): "
                               f"{output_file}[/yellow]")
            return output_file

    # ===== HELPERS =====

    def _get_residue_list(self) -> List[Dict[str, Any]]:
        """Build a list of residues in the site with metadata."""
        residues = []
        seen = set()

        # Identify center residues
        center_keys = set()
        for center in self.site.centers:
            center_keys.add((center.chain, center.resid))

        for atom in self.site.atoms:
            rkey = (atom.chain, atom.resid)
            if rkey in seen:
                # Add atom name to existing entry
                for res in residues:
                    if (res["chain"], res["resid"]) == rkey:
                        if atom.atom_name not in res["atom_names"]:
                            res["atom_names"].append(atom.atom_name)
                        break
            else:
                seen.add(rkey)
                residues.append({
                    "chain": atom.chain,
                    "resid": atom.resid,
                    "resname": atom.resname,
                    "atom_names": [atom.atom_name],
                    "is_center": rkey in center_keys,
                })

        return residues
