"""General-purpose hydrogen addition + interactive editing.

This module centralizes the hydrogen workflow that was originally embedded in
``SmallMolWorkflowRunner``: analyze H content, add hydrogens with AmberTools
``reduce``, inspect in the structure viewer, and interactively add/remove
individual H atoms until the protonation state is correct.

It is deliberately chemistry-agnostic — the user curates the protonation state
visually — so it works for any residue or covalent adduct, not just small
molecules. The small-molecule parameterizer and the modified-amino-acid
"from structure" route both drive the same :class:`HydrogenEditor`.

The ``module`` argument on the reduce/config helpers and on
:class:`HydrogenEditor` names the calling module in interactive prompts; it
feeds ProPrep's session-recording keys, so each caller passes its own module
name to keep replay matching stable.
"""

import os
import re
import subprocess

from rich.panel import Panel

from proprep.utils.prompts import prompt_with_context, confirm_with_context

# Session-recording module label used by the original small-molecule caller.
# Kept as the default so existing callers (small molecule parameterizer,
# structure preprocessor) preserve their exact prompt/recording keys.
_DEFAULT_MODULE = "Small Molecule Parameterizer"


# ── reduce availability / configuration / execution ────────────────────────

def check_reduce_availability():
    """Check if reduce program is available."""
    try:
        # Try different common flags to check reduce availability
        for flag in ['-version', '-Version', '-h']:
            try:
                result = subprocess.run(['reduce', flag],
                                      capture_output=True, text=True, timeout=10)

                # Check if reduce executed (even with non-zero exit code)
                # and produced expected output in either stdout or stderr
                output = result.stdout + result.stderr
                if 'reduce' in output.lower() and ('version' in output.lower() or 'usage' in output.lower()):
                    # Extract version info from either stdout or stderr
                    version_info = result.stderr.strip() if result.stderr else result.stdout.strip()
                    if version_info:
                        version_line = version_info.split('\n')[0]
                        return True, version_line
                    else:
                        return True, f"reduce available (tested with {flag})"

            except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
                continue

        # If we get here, none of the flags worked
        return False, "reduce command failed with all test flags"

    except FileNotFoundError:
        return False, "reduce not found in PATH"
    except Exception as e:
        return False, f"error testing reduce: {str(e)}"


def configure_reduce_options_aligned(interactive=True, console=None, processor=None,
                                     module=_DEFAULT_MODULE):
    """Configure reduce options with aligned text formatting."""
    if not interactive:
        return ["-NOFLIP", "-OH"]

    # Display comprehensive options explanation
    options_help = (
        "[bold]Reduce Hydrogen Addition Options:[/bold]\n\n"
        "[cyan]NOFLIP vs FLIP mode:[/cyan]\n"
        "  [green]NOFLIP[/green] (Recommended for small molecules)\n"
        "    Adds hydrogens in their initial orientations without optimization.\n"
        "    Safer for ligands/cofactors where orientation matters.\n\n"
        "  [yellow]FLIP[/yellow] (Use with caution)\n"
        "    Optimizes orientations of Asn, Gln, His sidechains by flipping.\n"
        "    Designed for proteins; may incorrectly alter small molecule geometry.\n\n"
        "[cyan]Additional options:[/cyan]\n"
        "  [grey50]-OH[/grey50]      Add hydrogens to hydroxyl (OH) and thiol (SH) groups\n"
        "            Important for alcohols, phenols, and thiols.\n\n"
        "  [grey50]-NUCLEAR[/grey50] Use nuclear X-H distances (~1.08 Å for C-H)\n"
        "            Default uses electron cloud distances (~1.0 Å for C-H).\n"
        "            Nuclear distances better match neutron/spectroscopy/QM data.\n\n"
        "  [grey50]-Keep[/grey50]    Preserve existing bond lengths from input structure\n"
        "            Useful if your input has optimized geometry."
    )

    console.print(Panel(options_help, title="Hydrogen Addition Options",
                        border_style="cyan", expand=False))

    options = []

    # Main mode with clearer prompt
    console.print("\n[bold]Select hydrogen addition mode:[/bold]")
    console.print("  [green]NOFLIP[/green] - Conservative, preserves input orientations (recommended)")
    console.print("  [yellow]FLIP[/yellow]   - Aggressive, optimizes orientations (protein-focused)")
    flip_mode = prompt_with_context(
        processor,
        "Mode",
        choices=["NOFLIP", "FLIP"],
        default="NOFLIP",
        module=module,
        description="Select reduce hydrogen addition mode"
    )
    options.append(f"-{flip_mode}")

    # OH/SH hydrogens - important for most small molecules
    console.print("\n[bold]Hydroxyl and thiol groups:[/bold]")
    add_oh = confirm_with_context(
        processor,
        "Add hydrogens to OH and SH groups? (recommended for most molecules)",
        default=True,
        module=module,
        description="Add hydrogens to OH/SH groups"
    )
    if add_oh:
        options.append("-OH")

    # Nuclear distances - explain the difference
    console.print("\n[bold]Bond length convention:[/bold]")
    console.print("[grey50]  Electron cloud: ~1.0 Å for C-H (default, X-ray crystallography convention)[/grey50]")
    console.print("[grey50]  Nuclear:        ~1.08 Å for C-H (neutron diffraction/spectroscopy/QM)[/grey50]")
    use_nuclear = confirm_with_context(
        processor,
        "Use nuclear X-H distances? (recommended for QM calculations)",
        default=False,
        module=module,
        description="Use nuclear X-H distances"
    )
    if use_nuclear:
        options.append("-NUCLEAR")

    # Keep bond lengths
    keep_bonds = confirm_with_context(
        processor,
        "Keep original bond lengths from input structure?",
        default=False,
        module=module,
        description="Keep original bond lengths"
    )
    if keep_bonds:
        options.append("-Keep")

    return options


def run_reduce_aligned(input_pdb, output_pdb, options, console):
    """Run reduce with aligned command explanations."""
    cmd = ["reduce"] + options + [input_pdb]

    console.print(f"\n[cyan]Running reduce:[/cyan]")
    console.print(f"[grey50]  {' '.join(cmd)} > {output_pdb}[/grey50]")

    # Aligned explanations for selected options
    explanations = {
        "-NOFLIP": "Conservative mode - preserves input orientations",
        "-FLIP": "Aggressive mode - optimizes Asn/Gln/His orientations",
        "-OH": "Adding hydrogens to hydroxyl and thiol groups",
        "-NUCLEAR": "Using nuclear X-H distances (~1.08 Å)",
        "-Keep": "Preserving original bond lengths"
    }

    console.print(f"\n[cyan]Selected options:[/cyan]")
    for option in options:
        if option in explanations:
            console.print(f"  {option}: {explanations[option]}")

    try:
        with open(output_pdb, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE,
                                  text=True, timeout=60)

        if result.returncode == 0:
            if result.stderr:
                # Parse and summarize reduce output
                stderr_lines = result.stderr.strip().split('\n')
                warnings = []
                info_lines = []
                skipped = []
                added_count = 0

                for line in stderr_lines:
                    if 'WARNING:' in line:
                        # Deduplicate common warnings
                        if 'will be treated as hydrogen' in line:
                            continue  # Skip these noisy warnings
                        warnings.append(line)
                    elif 'SKIPPED' in line:
                        skipped.append(line)
                    elif 'Added' in line and 'hydrogens' in line:
                        # Extract hydrogen count
                        match = re.search(r'Added (\d+) hydrogens', line)
                        if match:
                            added_count = int(match.group(1))
                        info_lines.append(line)
                    elif any(kw in line for kw in ['orientation', 'Processing', 'Building']):
                        continue  # Skip verbose processing info
                    elif line.strip() and not line.startswith('If you publish'):
                        info_lines.append(line)

                # Display summarized output
                if added_count > 0:
                    console.print(f"\n[green]✓ Added {added_count} hydrogen atoms[/green]")

                if skipped:
                    console.print(f"\n[yellow]Skipped atoms ({len(skipped)}):[/yellow]")
                    for skip in skipped[:3]:  # Show first 3
                        console.print(f"[grey50]  {skip}[/grey50]")
                    if len(skipped) > 3:
                        console.print(f"[grey50]  ... and {len(skipped)-3} more[/grey50]")

                if warnings:
                    unique_warnings = list(set(warnings))[:5]  # Dedupe and limit
                    console.print(f"\n[yellow]Warnings ({len(unique_warnings)}):[/yellow]")
                    for warn in unique_warnings:
                        console.print(f"[grey50]  {warn}[/grey50]")

            return True, "Hydrogens added successfully"
        else:
            error_msg = result.stderr if result.stderr else "Unknown error"
            return False, f"Reduce failed: {error_msg}"

    except subprocess.TimeoutExpired:
        return False, "Reduce timed out"
    except Exception as e:
        return False, f"Error running reduce: {str(e)}"


# ── hydrogen analysis + viewer ─────────────────────────────────────────────

def analyze_hydrogens(pdb_file):
    """Count hydrogens vs total atoms in a PDB file.

    Returns:
        (hydrogen_count, total_atoms, hydrogen_atom_names)
    """
    hydrogen_count = 0
    total_atoms = 0
    hydrogen_atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                total_atoms += 1
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                if element == 'H' or (not element and atom_name.startswith('H')):
                    hydrogen_count += 1
                    hydrogen_atoms.append(atom_name)
    return hydrogen_count, total_atoms, hydrogen_atoms


def launch_viewer(pdb_file: str, residue_name: str, console) -> bool:
    """Launch the interactive structure viewer for a single-residue PDB.

    Routed through the ``ViewerCoordinator`` singleton so the launch shares
    state with all other coordinator-driven hooks across the session.

    Args:
        pdb_file: Path to the PDB file to display
        residue_name: Name of the residue being viewed
        console: Rich console for output

    Returns:
        True if viewer launched successfully
    """
    try:
        from proprep.structure_prep.viewer_coordinator import viewer as _viewer
    except ImportError:
        console.print("[yellow]Structure viewer not available[/yellow]")
        console.print("[grey50]Install with: pip install nglview[/grey50]")
        return False

    if not os.path.exists(pdb_file):
        console.print(f"[yellow]Structure file not found: {pdb_file}[/yellow]")
        return False

    # Convert to absolute path - important because the working directory
    # may change after the viewer is launched (in a background thread)
    abs_pdb_file = os.path.abspath(pdb_file)

    # force=True is correct here: the caller gates this function behind a
    # "Launch interactive structure viewer?" Y/N prompt, so by the time we
    # arrive the user has explicitly opted in. The ball+stick overlay across
    # 'all' atoms gives the same single-residue view the original did.
    try:
        _viewer.show_structure(abs_pdb_file, force=True)
        _viewer.highlight(
            "all",
            style="ball+stick",
            color="element",
            label="small_molecule",
            force=True,
        )
        return True
    except Exception as e:
        console.print(f"[yellow]Could not launch viewer: {e}[/yellow]")
        return False


# ── interactive hydrogen editor ────────────────────────────────────────────

class HydrogenEditor:
    """Add hydrogens to a PDB and let the user curate the protonation state.

    The editor operates in place on ``pdb_file`` (which may be reassigned to a
    reduce-produced ``*_H.pdb``). Call :meth:`run` for the full add-then-curate
    workflow; :meth:`add_interactive` / :meth:`remove_interactive` are the
    individual editing loops.

    Args:
        pdb_file: Path to the heavy-atom (or partially protonated) PDB.
        name: Output basename stem (reduce writes ``{name}_H.pdb``).
        console: Rich console.
        processor: ProPrep processor (for session recording); optional.
        interactive: Whether to prompt the user.
        residue_name: Display name for the viewer; defaults to ``name``.
        module: Session-recording module label for prompts.
    """

    # Standard nuclear X-H bond lengths (Å) — match reduce -NUCLEAR
    BOND_LENGTHS = {'C': 1.08, 'N': 1.01, 'O': 0.96, 'S': 1.34}
    DEFAULT_BOND_LENGTH = 1.08

    def __init__(self, pdb_file, name, *, console, processor=None,
                 interactive=True, residue_name=None, module=_DEFAULT_MODULE):
        self.pdb_file = pdb_file
        self.name = name
        self.residue_name = residue_name or name
        self.console = console
        self.processor = processor
        self.interactive = interactive
        self.module = module

    def run(self):
        """Analyze H content, add via reduce, inspect, and edit interactively.

        Returns:
            dict with keys ``pdb_file`` (possibly updated), ``h_count``,
            ``total``, and ``summary``.
        """
        h_count, total_count, h_atoms = analyze_hydrogens(self.pdb_file)
        h_percentage = (h_count / total_count * 100) if total_count > 0 else 0

        self.console.print(f"[cyan]Hydrogen analysis:[/cyan]")
        self.console.print(f"  Total atoms: {total_count}")
        self.console.print(f"  Hydrogen atoms: {h_count} ({h_percentage:.1f}%)")

        if h_count > 0:
            self.console.print(f"  Hydrogen atoms found: {', '.join(h_atoms[:10])}")
            if len(h_atoms) > 10:
                self.console.print(f"    ... and {len(h_atoms)-10} more")

        likely_missing = h_percentage < 30 and total_count > 3

        if likely_missing:
            self.console.print(f"[yellow]Low hydrogen content detected ({h_percentage:.1f}%)[/yellow]")
            self.console.print(f"[yellow]   Small organic molecules typically have ~40-60% hydrogen atoms[/yellow]")

            if self.interactive:
                add_hydrogens = confirm_with_context(
                    self.processor, "Would you like to add hydrogens using reduce?",
                    default=True, module=self.module,
                    description="Add hydrogens using reduce",
                )
            else:
                add_hydrogens = False
                self.console.print("[yellow]Non-interactive mode: skipping hydrogen addition[/yellow]")

            if add_hydrogens:
                reduce_available, reduce_info = check_reduce_availability()

                if reduce_available:
                    self.console.print(f"[green]✓ Found reduce: {reduce_info}[/green]")
                    reduce_options = configure_reduce_options_aligned(
                        self.interactive, self.console, self.processor, module=self.module)

                    h_added_file = f"{self.name}_H.pdb"
                    success, message = run_reduce_aligned(self.pdb_file, h_added_file, reduce_options, self.console)

                    if success:
                        self.pdb_file = h_added_file

                        # Offer structure visualization
                        self.console.print(f"[cyan]Verify hydrogen placement before proceeding.[/cyan]")
                        self.console.print(f"[grey50]You can also inspect {h_added_file} with GaussView, PyMOL, or similar tools.[/grey50]")

                        if self.interactive:
                            view_structure = confirm_with_context(
                                self.processor, "Launch interactive structure viewer?",
                                default=True, module=self.module,
                                description="Launch structure viewer",
                            )
                            if view_structure:
                                launch_viewer(h_added_file, self.residue_name, self.console)

                        # Offer atom editing (loop until done)
                        if self.interactive:
                            edited = False
                            while True:
                                edit_choice = prompt_with_context(
                                    self.processor,
                                    "Edit H atoms?",
                                    choices=["add", "remove", "done"],
                                    default="done", module=self.module,
                                    description="Edit H atoms",
                                ).strip().lower()
                                if edit_choice == 'done':
                                    break
                                elif edit_choice == 'add':
                                    if self.add_interactive():
                                        edited = True
                                elif edit_choice == 'remove':
                                    if self.remove_interactive():
                                        edited = True

                            if edited:
                                relaunch = confirm_with_context(
                                    self.processor, "Relaunch structure viewer to verify?",
                                    default=True, module=self.module,
                                    description="Relaunch structure viewer",
                                )
                                if relaunch:
                                    launch_viewer(self.pdb_file, self.residue_name, self.console)

                        # Recount after any interactive edits
                        final_h_count, final_total, _ = analyze_hydrogens(self.pdb_file)
                        total_added = final_h_count - h_count
                        return {'pdb_file': self.pdb_file, 'h_count': final_h_count,
                                'total': final_total,
                                'summary': f'Added {total_added} hydrogens ({final_total} atoms total)'}
                    else:
                        self.console.print(f"[red]✗ {message}[/red]")
                        self.console.print(f"[yellow]Proceeding with original structure[/yellow]")
                else:
                    self.console.print(f"[red]✗ Reduce not available: {reduce_info}[/red]")
                    self.console.print(
                        f"[yellow]Install reduce via conda ([bold]conda install dacase::ambertools-dac=25[/bold]) "
                        f"or add hydrogens manually in GaussView before running Gaussian[/yellow]"
                    )
            else:
                self.console.print(f"[yellow]Skipping hydrogen addition[/yellow]")
        else:
            self.console.print(f"[green]✓ Structure appears to have adequate hydrogens ({h_percentage:.1f}%)[/green]")
            if h_percentage > 70:
                self.console.print(f"[cyan]Note: High hydrogen percentage - verify this is expected[/cyan]")

        return {'pdb_file': self.pdb_file, 'h_count': h_count, 'total': total_count,
                'summary': f'Hydrogen analysis complete ({h_count}/{total_count} atoms are H)'}

    def add_interactive(self) -> bool:
        """Add a hydrogen atom bonded to a specified heavy atom.

        The user provides the heavy atom name and its bonded neighbors. The H
        is placed along the direction opposite the existing bonds, at a
        standard bond length for the element type.

        Returns:
            True if any atoms were added.
        """
        import numpy as np

        BOND_LENGTHS = self.BOND_LENGTHS
        DEFAULT_BOND_LENGTH = self.DEFAULT_BOND_LENGTH

        # Parse atoms from PDB. A capped conjugate is a MULTI-residue model
        # (ACE / amino acid / NME / ligand), so bare atom names can repeat — e.g.
        # backbone "N" in both the AA and NME, "CH3" in ACE and NME, "C"/"O" in
        # ACE and the AA. Keying by name alone would silently collapse these to
        # one pickable atom (the bug where two peptide "N" atoms are
        # indistinguishable). Label each atom by its bare name when that name is
        # unique in the file, else by "name@resSeq" so the user can target the
        # right one. A single-residue small molecule has all-unique names, so its
        # labels are just the plain names (unchanged UX).
        raw = []  # (name, resseq, x, y, z, element, line)
        name_counts = {}
        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    name = line[12:16].strip()
                    resseq = line[22:26].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip()
                    raw.append((name, resseq, x, y, z, element, line))
                    name_counts[name] = name_counts.get(name, 0) + 1

        atoms = {}  # label -> (x, y, z, element, line)
        for name, resseq, x, y, z, element, line in raw:
            label = name if name_counts[name] == 1 else f"{name}@{resseq}"
            atoms[label] = (x, y, z, element, line)

        ambiguous = any(c > 1 for c in name_counts.values())
        non_h_names = sorted(n for n, (_, _, _, e, _) in atoms.items()
                             if not (e == 'H' or (not e and n.split('@')[0].startswith('H'))))
        self.console.print(f"\n[cyan]Heavy atoms:[/cyan] {', '.join(non_h_names)}")
        if ambiguous:
            self.console.print(
                "[grey50]Repeated atom names are shown as name@residue "
                "(e.g. N@44 vs N@45); use that exact label to pick one.[/grey50]")

        added_any = False
        while True:
            heavy_name = prompt_with_context(
                self.processor,
                "Heavy atom to bond H to (e.g. 'N1'), or 'done'",
                default="done", module=self.module,
                description="Heavy atom for H addition",
            ).strip()
            if heavy_name.lower() == 'done':
                break

            if heavy_name not in atoms:
                self.console.print(f"[red]Atom '{heavy_name}' not found.[/red]")
                continue

            neighbors_input = prompt_with_context(
                self.processor,
                f"Atoms bonded to {heavy_name} (space-separated, e.g. 'C8A C2')",
                default="", module=self.module,
                description="Bonded neighbor atoms",
            ).strip()
            if not neighbors_input:
                self.console.print("[yellow]No neighbors specified, skipping.[/yellow]")
                continue

            neighbor_names = neighbors_input.split()
            invalid = [n for n in neighbor_names if n not in atoms]
            if invalid:
                self.console.print(f"[red]Atom names not found: {', '.join(invalid)}[/red]")
                continue

            # Compute H position
            hx, hy, hz, element, _ = atoms[heavy_name]
            heavy_pos = np.array([hx, hy, hz])

            bond_vectors = []
            for nb in neighbor_names:
                nx, ny, nz, _, _ = atoms[nb]
                v = np.array([nx, ny, nz]) - heavy_pos
                norm = np.linalg.norm(v)
                if norm > 0:
                    bond_vectors.append(v / norm)

            if not bond_vectors:
                self.console.print("[red]Could not compute bond direction.[/red]")
                continue

            # H direction is opposite the sum of normalized bond vectors
            direction = -np.sum(bond_vectors, axis=0)
            norm = np.linalg.norm(direction)
            if norm < 1e-6:
                self.console.print("[red]Bond vectors cancel out — cannot determine H placement.[/red]")
                continue
            direction = direction / norm

            bond_length = BOND_LENGTHS.get(element, DEFAULT_BOND_LENGTH)
            h_pos = heavy_pos + direction * bond_length

            # Generate unique H name
            h_num = 1
            while f"H{h_num}" in atoms:
                h_num += 1
            h_name = f"H{h_num}"

            # Build HETATM line using the heavy atom's line as template
            # for residue name, chain, etc.
            template = atoms[heavy_name][4]
            res_name = template[17:20]
            chain = template[21:22]
            res_seq = template[22:26]
            serial = len(atoms) + 1
            h_line = (
                f"HETATM{serial:5d} {h_name:<4s} {res_name}{chain}{res_seq}"
                f"    {h_pos[0]:8.3f}{h_pos[1]:8.3f}{h_pos[2]:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}          {'H':>2s}  \n"
            )

            # Insert H line after the last HETATM/ATOM line
            lines = []
            last_atom_idx = -1
            with open(self.pdb_file, 'r') as f:
                lines = f.readlines()
            for i, line in enumerate(lines):
                if line.startswith(('ATOM', 'HETATM')):
                    last_atom_idx = i
            lines.insert(last_atom_idx + 1, h_line)

            with open(self.pdb_file, 'w') as f:
                f.writelines(lines)

            atoms[h_name] = (h_pos[0], h_pos[1], h_pos[2], 'H', h_line)
            added_any = True
            self.console.print(
                f"[green]✓ Added {h_name} bonded to {heavy_name} "
                f"(bond length {bond_length:.2f} Å, {len(atoms)} atoms total)[/green]"
            )

        return added_any

    def remove_interactive(self) -> bool:
        """Interactive loop to remove misplaced atoms from the PDB.

        Returns:
            True if any atoms were removed.
        """
        atom_names = []
        with open(self.pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    name = line[12:16].strip()
                    element = line[76:78].strip()
                    atom_names.append((name, element))

        h_names = [n for n, e in atom_names if e == 'H' or (not e and n.startswith('H'))]
        self.console.print(f"\n[cyan]Hydrogens ({len(h_names)}):[/cyan] {', '.join(h_names)}")

        removed_any = False
        while True:
            atoms_input = prompt_with_context(
                self.processor,
                "Atom names to remove (space-separated, e.g. 'H7 H14'), or 'done'",
                default="done", module=self.module,
                description="Atom names to remove",
            )
            if atoms_input.strip().lower() == 'done':
                break

            names_to_remove = atoms_input.strip().split()
            all_atom_names = [n for n, _ in atom_names]

            invalid = [n for n in names_to_remove if n not in all_atom_names]
            if invalid:
                self.console.print(f"[red]Atom names not found: {', '.join(invalid)}[/red]")
                self.console.print(f"[yellow]Available: {', '.join(all_atom_names)}[/yellow]")
                continue

            kept_lines = []
            removed_count = 0
            removed_set = set(names_to_remove)
            with open(self.pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        atom_name = line[12:16].strip()
                        if atom_name in removed_set:
                            removed_count += 1
                            removed_set.discard(atom_name)
                            continue
                    kept_lines.append(line)

            with open(self.pdb_file, 'w') as f:
                f.writelines(kept_lines)

            removed_any = True
            self.console.print(f"[green]Removed {removed_count} atom(s): {', '.join(names_to_remove)}[/green]")

            atom_names = [(n, e) for n, e in atom_names if n not in set(names_to_remove)]
            h_names = [n for n, e in atom_names if e == 'H' or (not e and n.startswith('H'))]
            self.console.print(f"  Remaining atoms: {len(atom_names)} ({len(h_names)} H)")

        return removed_any
