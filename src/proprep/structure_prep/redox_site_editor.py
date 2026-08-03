"""
Redox Site Editor

Comprehensive editor for RedoxSite objects, allowing detailed modification of:
- Centers (redox-active centers)
- Residues (add/remove entire residues)
- Atoms (fine-grained atom operations)
- Bonds (add/remove/reclassify bonds)
- Metadata (site type, detection parameters)
"""

import copy
import logging
from typing import Dict, List, Optional, Tuple, Any

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.prompts import prompt_with_context, confirm_with_context

logger = logging.getLogger(__name__)


class RedoxSiteEditor:
    """Comprehensive editor for RedoxSite objects with user-driven operations"""

    def __init__(self, sites: List, console: Console, processor=None, workspace=None):
        """
        Initialize the RedoxSite editor.

        Args:
            sites: List of RedoxSite objects to edit
            console: Rich console for output
            processor: ProPrep processor (optional, for session recording)
            workspace: ProPrep workspace (for loading structure)
        """
        self.original_sites = sites
        self.working_sites = None  # Deep copy for editing
        self.console = console
        self.processor = processor
        self.workspace = workspace
        self.current_site = None
        self.current_site_index = None
        self._pdb_file = None  # Cached PDB file path
        self._structure = None  # Cached BioPython structure

    def _get_structure(self):
        """Get BioPython structure using StructureSelector priority system"""
        if self._structure is not None:
            return self._structure

        if not self.workspace:
            return None

        # Use StructureSelector to get PDB file
        from proprep.utils.structure_selector import StructureSelector
        selector = StructureSelector(self.workspace, self.console, self.processor)

        self._pdb_file = selector.get_structure(interactive=False, silent=True)

        if not self._pdb_file:
            return None

        # Load structure
        try:
            from Bio.PDB import PDBParser
            parser = PDBParser(QUIET=True)
            self._structure = parser.get_structure("loaded", self._pdb_file)
            return self._structure
        except Exception as e:
            self.console.print(f"[red]Error loading structure: {e}[/red]")
            return None

    def run(self) -> Optional[List]:
        """
        Main editor loop.

        Returns:
            Modified sites list if saved, None if cancelled
        """
        # Create deep copy for editing
        self.working_sites = copy.deepcopy(self.original_sites)

        while True:
            choice = self._show_site_selection_menu()

            if choice == 'save':
                return self.working_sites
            elif choice == 'cancel':
                if self._confirm_cancel():
                    return None
                else:
                    continue
            elif choice.startswith('remove '):
                # Parse inline selection: "remove 3,27,47" or "remove 1-5" or "remove all"
                selection = choice[7:].strip()  # Remove "remove " prefix
                if selection:
                    self._remove_sites(selection)
                else:
                    self.console.print("[red]Please specify sites to remove (e.g., 'remove 3' or 'remove 1-5')[/red]")
            elif choice == 'remove':
                self.console.print("[yellow]Please specify sites to remove (e.g., 'remove 3', 'remove 1,3,5', 'remove 1-5', 'remove all')[/yellow]")
            elif choice.startswith('edit '):
                # Parse inline selection: "edit 3" or "edit 1,3,5"
                selection = choice[5:].strip()  # Remove "edit " prefix
                if selection:
                    self._edit_sites_by_selection(selection)
                else:
                    self.console.print("[red]Please specify site(s) to edit (e.g., 'edit 3' or 'edit 1,3,5')[/red]")
            elif choice == 'edit':
                self.console.print("[yellow]Please specify site(s) to edit (e.g., 'edit 3', 'edit 1,3,5')[/yellow]")
            elif choice == 'validate' or choice == 'validate_all':
                self._validate_all_sites()
            else:
                self.console.print("[red]Invalid choice. Use 'edit <N>', 'remove <N>', 'validate', 'save', or 'cancel'[/red]")

    def _show_site_selection_menu(self) -> str:
        """Show site selection menu and get user choice"""
        self.console.print("\n" + "=" * 60)
        self.console.print("[bold cyan]RedoxSite Editor - Site Selection[/bold cyan]")
        self.console.print("=" * 60 + "\n")

        # Show sites table
        table = Table(title="Detected Redox Sites")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Site ID", style="green")
        table.add_column("Centers", style="yellow", justify="right")
        table.add_column("Residues", style="yellow", justify="right")
        table.add_column("Atoms", style="yellow", justify="right")
        table.add_column("Bonds", style="yellow", justify="right")
        table.add_column("Site Type", style="magenta")

        for i, site in enumerate(self.working_sites, 1):
            num_centers = len(site.centers) if hasattr(site, 'centers') else 0
            num_residues = len(site.residue_groups) if hasattr(site, 'residue_groups') else 0
            num_atoms = len(site.atoms) if hasattr(site, 'atoms') else 0
            num_bonds = len(site.bonds) if hasattr(site, 'bonds') else 0
            site_type = site.site_type if hasattr(site, 'site_type') and site.site_type else "N/A"

            table.add_row(
                str(i),
                site.site_id,
                str(num_centers),
                str(num_residues),
                str(num_atoms),
                str(num_bonds),
                site_type
            )

        self.console.print(table)

        # Menu options
        self.console.print("\n[bold]Options:[/bold]")
        self.console.print("  [cyan]edit <N>[/cyan]     Edit site(s) by number (e.g., 'edit 3' or 'edit 1,3,5')")
        self.console.print("  [cyan]remove <N>[/cyan]   Remove site(s) (e.g., 'remove 3', 'remove 1-5', 'remove all')")
        self.console.print("  [cyan]validate[/cyan]     Validate all sites")
        self.console.print("  [cyan]save[/cyan]         Save changes and exit")
        self.console.print("  [cyan]cancel[/cyan]       Cancel and discard changes")

        choice = prompt_with_context(
            self.processor,
            "\n[bold]Choice[/bold]",
            default="cancel",
            module="RedoxSite Editor",
            description="Select site to edit or perform operation"
        )
        return choice.strip().lower()

    def _edit_sites_by_selection(self, selection: str):
        """Edit one or more sites by selection string"""
        if not self.working_sites:
            self.console.print("[yellow]No sites to edit[/yellow]")
            return

        try:
            indices = self._parse_selection(selection, len(self.working_sites))
        except ValueError as e:
            self.console.print(f"[red]Invalid selection: {e}[/red]")
            return

        if not indices:
            self.console.print("[yellow]No sites selected for editing[/yellow]")
            return

        # Edit each site in sequence
        for idx in sorted(indices):
            site_index = idx - 1  # Convert to 0-based
            if 0 <= site_index < len(self.working_sites):
                self._edit_single_site(site_index)

    def _remove_sites(self, selection: str = None):
        """Remove one or more sites from the list"""
        if not self.working_sites:
            self.console.print("[yellow]No sites to remove[/yellow]")
            return

        # If no selection provided, prompt for it (fallback for backward compatibility)
        if selection is None:
            self.console.print("\n[bold cyan]Remove Site(s)[/bold cyan]\n")
            self.console.print("[grey50]Enter site index/indices to remove (e.g., '6', '6,18,30', '1-3', or 'all')[/grey50]")
            self.console.print("[grey50]Enter 'cancel' to go back[/grey50]")

            selection = prompt_with_context(
                self.processor,
                "\n[bold]Site(s) to remove[/bold]",
                default="cancel",
                module="RedoxSite Editor",
                description="Select sites to remove from workspace"
            )

            if selection.lower() == 'cancel':
                return

        # Parse the selection (reuse the parsing logic from redox_detector_module)
        try:
            indices_to_remove = self._parse_selection(selection, len(self.working_sites))
        except ValueError as e:
            self.console.print(f"[red]Invalid selection: {e}[/red]")
            return

        if not indices_to_remove:
            self.console.print("[yellow]No sites selected for removal[/yellow]")
            return

        # Show what will be removed
        self.console.print("\n[bold]Sites to be removed:[/bold]")
        for idx in sorted(indices_to_remove):
            site = self.working_sites[idx - 1]
            self.console.print(f"  {idx}. {site.site_id} ({len(site.centers)} center(s), {len(site.atoms)} atoms)")

        # Confirm removal
        if confirm_with_context(
            self.processor,
            "\nConfirm removal of these site(s)?",
            default=False,
            module="RedoxSite Editor",
            description="Confirm site removal"
        ):
            # Remove the sites (in reverse order to maintain indices)
            for idx in sorted(indices_to_remove, reverse=True):
                removed_site = self.working_sites.pop(idx - 1)
                self.console.print(f"[green]✓ Removed {removed_site.site_id}[/green]")

            remaining = len(self.working_sites)
            self.console.print(f"\n[green]✓ Successfully removed {len(indices_to_remove)} site(s)[/green]")
            self.console.print(f"[cyan]Remaining sites: {remaining}[/cyan]")
        else:
            self.console.print("[grey50]Removal cancelled[/grey50]")

    def _parse_selection(self, selection: str, max_index: int) -> set:
        """
        Parse user selection string into a set of indices.

        Supports:
        - Single index: '1'
        - Multiple indices: '1,3,5'
        - Ranges: '1-3' (inclusive)
        - Combined: '1,3-5,7'
        - All: 'all'

        Returns:
            Set of 1-based indices

        Raises:
            ValueError: If selection is invalid
        """
        selection = selection.strip().lower()

        if selection == 'all':
            return set(range(1, max_index + 1))

        indices = set()
        parts = selection.split(',')

        for part in parts:
            part = part.strip()

            if '-' in part:
                # Range
                try:
                    start, end = part.split('-')
                    start_idx = int(start.strip())
                    end_idx = int(end.strip())

                    if start_idx < 1 or end_idx > max_index:
                        raise ValueError(f"Range {start_idx}-{end_idx} is out of bounds (1-{max_index})")
                    if start_idx > end_idx:
                        raise ValueError(f"Invalid range: {start_idx}-{end_idx}")

                    indices.update(range(start_idx, end_idx + 1))
                except ValueError as e:
                    if "invalid literal" in str(e):
                        raise ValueError(f"Invalid range format: '{part}'")
                    else:
                        raise
            else:
                # Single index
                try:
                    idx = int(part)
                    if idx < 1 or idx > max_index:
                        raise ValueError(f"Index {idx} is out of bounds (1-{max_index})")
                    indices.add(idx)
                except ValueError:
                    raise ValueError(f"Invalid index: '{part}'")

        return indices

    def _edit_single_site(self, site_index: int):
        """Edit a single RedoxSite with hierarchical menu"""
        self.current_site = self.working_sites[site_index]
        self.current_site_index = site_index

        while True:
            choice = self._show_site_editor_menu()

            if choice == 'back':
                break
            elif choice in ['1', 'centers']:
                self._edit_centers_menu()
            elif choice in ['2', 'residues']:
                self._edit_residues_menu()
            elif choice in ['3', 'atoms']:
                self._edit_atoms_menu()
            elif choice in ['4', 'bonds']:
                self._edit_bonds_menu()
            elif choice in ['5', 'metadata']:
                self._edit_metadata_menu()
            elif choice == 'validate':
                self._validate_site(self.current_site)
            elif choice == 'summary':
                self._show_site_summary(self.current_site)

    def _show_site_editor_menu(self) -> str:
        """Show editor menu for current site"""
        self.console.print("\n" + "=" * 60)
        self.console.print(f"[bold cyan]Editing Site: {self.current_site.site_id}[/bold cyan]")
        self.console.print("=" * 60 + "\n")

        # Quick summary
        num_centers = len(self.current_site.centers)
        num_residues = len(self.current_site.residue_groups)
        num_atoms = len(self.current_site.atoms)
        num_bonds = len(self.current_site.bonds)

        self.console.print(f"  Centers: {num_centers} | Residues: {num_residues} | Atoms: {num_atoms} | Bonds: {num_bonds}\n")

        # Menu options
        self.console.print("[bold]Edit Operations:[/bold]")
        self.console.print("  [cyan]1. centers[/cyan]   Edit redox centers")
        self.console.print("  [cyan]2. residues[/cyan]  Edit residues (add/remove/replace)")
        self.console.print("  [cyan]3. atoms[/cyan]     Edit individual atoms")
        self.console.print("  [cyan]4. bonds[/cyan]     Edit bonds")
        self.console.print("  [cyan]5. metadata[/cyan]  Edit site metadata")
        self.console.print("\n[bold]Other Operations:[/bold]")
        self.console.print("  [cyan]summary[/cyan]   View detailed site summary")
        self.console.print("  [cyan]validate[/cyan]  Validate site consistency")
        self.console.print("  [cyan]back[/cyan]      Return to site selection")

        choice = prompt_with_context(
            self.processor,
            "\n[bold]Choice[/bold]",
            default="back",
            module="RedoxSite Editor",
            description=f"Edit site {self.current_site.site_id}"
        )
        return choice.strip().lower()

    # ========================================================================
    # CENTERS EDITING
    # ========================================================================

    def _edit_centers_menu(self):
        """Menu for editing redox centers"""
        while True:
            choice = self._show_centers_menu()

            if choice == 'back':
                break
            elif choice == 'list':
                self._list_centers()
            elif choice == 'remove':
                self._remove_center()
            elif choice == 'modify':
                self._modify_center()

    def _show_centers_menu(self) -> str:
        """Show centers editing menu"""
        self.console.print("\n" + "-" * 60)
        self.console.print("[bold yellow]Edit Centers[/bold yellow]")
        self.console.print("-" * 60 + "\n")

        self.console.print("[bold]Operations:[/bold]")
        self.console.print("  [cyan]list[/cyan]    List all centers")
        self.console.print("  [cyan]remove[/cyan]  Remove a center")
        self.console.print("  [cyan]modify[/cyan]  Modify center properties")
        self.console.print("  [cyan]back[/cyan]    Return to site editor")

        choice = prompt_with_context(
            self.processor,
            "\n[bold]Choice[/bold]",
            default="back",
            module="RedoxSite Editor - Centers",
            description="Edit redox centers"
        )
        return choice.strip().lower()

    def _list_centers(self):
        """List all centers in current site"""
        if not self.current_site.centers:
            self.console.print("[yellow]No centers in this site[/yellow]")
            return

        table = Table(title=f"Centers in {self.current_site.site_id}")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Chain", style="green")
        table.add_column("Resname", style="green")
        table.add_column("ResID", style="yellow")
        table.add_column("Type", style="magenta")
        table.add_column("Atom", style="blue")
        table.add_column("Element", style="blue")

        for i, center in enumerate(self.current_site.centers, 1):
            center_type = center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type)
            table.add_row(
                str(i),
                center.chain,
                center.resname,
                str(center.resid),
                center_type,
                center.atom_name if center.atom_name else "N/A",
                center.element if center.element else "N/A"
            )

        self.console.print(table)

    def _remove_center(self):
        """Remove a center from the site"""
        if not self.current_site.centers:
            self.console.print("[yellow]No centers to remove[/yellow]")
            return

        self._list_centers()

        self.console.print("\n[grey50]Enter center index to remove, or 'cancel'[/grey50]")
        choice = prompt_with_context(
            self.processor,
            "[bold]Index[/bold]",
            default="cancel",
            module="RedoxSite Editor - Centers",
            description="Select center to remove"
        )

        if choice.lower() == 'cancel':
            return

        try:
            index = int(choice) - 1
            if 0 <= index < len(self.current_site.centers):
                center = self.current_site.centers[index]

                # Show what will be removed
                self.console.print(f"\n[yellow]Remove center:[/yellow] {center.chain}:{center.resname}:{center.resid}")

                if confirm_with_context(
                    self.processor,
                    "Confirm removal?",
                    default=False,
                    module="RedoxSite Editor - Centers",
                    description="Confirm center removal"
                ):
                    # Remove center
                    coords = center.coords
                    removed_center = self.current_site.centers.pop(index)

                    # Remove from coord_to_pdb mapping
                    if coords in self.current_site.coord_to_pdb:
                        del self.current_site.coord_to_pdb[coords]

                    self.console.print(f"[green]✓ Removed center {removed_center.chain}:{removed_center.resname}:{removed_center.resid}[/green]")
            else:
                self.console.print("[red]Invalid index[/red]")
        except ValueError:
            self.console.print("[red]Invalid input[/red]")

    def _modify_center(self):
        """Modify center properties"""
        if not self.current_site.centers:
            self.console.print("[yellow]No centers to modify[/yellow]")
            return

        self._list_centers()

        self.console.print("\n[grey50]Enter center index to modify, or 'cancel'[/grey50]")
        choice = prompt_with_context(
            self.processor,
            "[bold]Index[/bold]",
            default="cancel",
            module="RedoxSite Editor - Centers",
            description="Select center to modify"
        )

        if choice.lower() == 'cancel':
            return

        try:
            index = int(choice) - 1
            if 0 <= index < len(self.current_site.centers):
                center = self.current_site.centers[index]
                self._modify_center_properties(center)
            else:
                self.console.print("[red]Invalid index[/red]")
        except ValueError:
            self.console.print("[red]Invalid input[/red]")

    def _modify_center_properties(self, center):
        """Modify properties of a specific center"""
        self.console.print(f"\n[bold]Modifying center: {center.chain}:{center.resname}:{center.resid}[/bold]")
        self.console.print("\n[grey50]Press Enter to keep current value[/grey50]\n")

        # Chain
        new_chain = prompt_with_context(
            self.processor,
            f"Chain [{center.chain}]",
            default=center.chain,
            module="RedoxSite Editor - Centers",
            description="Modify center chain"
        )
        center.chain = new_chain

        # Resname
        new_resname = prompt_with_context(
            self.processor,
            f"Resname [{center.resname}]",
            default=center.resname,
            module="RedoxSite Editor - Centers",
            description="Modify center resname"
        )
        center.resname = new_resname

        # ResID
        new_resid = prompt_with_context(
            self.processor,
            f"ResID [{center.resid}]",
            default=str(center.resid),
            module="RedoxSite Editor - Centers",
            description="Modify center resid"
        )
        center.resid = int(new_resid)

        # Center type
        from .comprehensive_redox_detector import CenterType
        self.console.print("\nCenter types:")
        self.console.print("  1. metal_ion (isolated metal like Zn, Cu)")
        self.console.print("  2. organometallic_cofactor (metal in multi-atom residue like heme)")
        self.console.print("  3. organic_cofactor (non-standard residue without metal)")
        self.console.print("  4. redox_amino_acid")

        current_type = center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type)
        type_choice = prompt_with_context(
            self.processor,
            f"Center type [{current_type}]",
            default=current_type,
            module="RedoxSite Editor - Centers",
            description="Modify center type",
            options_map={'1': 'metal_ion', '2': 'organometallic_cofactor', '3': 'organic_cofactor', '4': 'redox_amino_acid'}
        )

        # Parse center type
        if type_choice in ['1', 'metal_ion']:
            center.center_type = CenterType.METAL_ION
        elif type_choice in ['2', 'organometallic_cofactor']:
            center.center_type = CenterType.ORGANOMETALLIC_COFACTOR
        elif type_choice in ['3', 'organic_cofactor']:
            center.center_type = CenterType.ORGANIC_COFACTOR
        elif type_choice in ['4', 'redox_amino_acid']:
            center.center_type = CenterType.REDOX_AMINO_ACID

        # Update coord_to_pdb mapping
        if center.coords in self.current_site.coord_to_pdb:
            self.current_site.coord_to_pdb[center.coords].update({
                'chain': center.chain,
                'resname': center.resname,
                'resid': center.resid,
                'center_type': center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type)
            })

        self.console.print("[green]✓ Center updated[/green]")

    # ========================================================================
    # RESIDUES EDITING
    # ========================================================================

    def _edit_residues_menu(self):
        """Menu for editing residues"""
        while True:
            choice = self._show_residues_menu()

            if choice == 'back':
                break
            elif choice == 'list':
                self._list_residues()
            elif choice == 'add':
                self._add_residue()
            elif choice == 'remove':
                self._remove_residue()
            elif choice == 'replace':
                self._replace_residue()

    def _show_residues_menu(self) -> str:
        """Show residues editing menu"""
        self.console.print("\n" + "-" * 60)
        self.console.print("[bold yellow]Edit Residues[/bold yellow]")
        self.console.print("-" * 60 + "\n")

        self.console.print("[bold]Operations:[/bold]")
        self.console.print("  [cyan]list[/cyan]     List all residues")
        self.console.print("  [cyan]add[/cyan]      Add residue from structure")
        self.console.print("  [cyan]remove[/cyan]   Remove a residue")
        self.console.print("  [cyan]replace[/cyan]  Replace a residue with another")
        self.console.print("  [cyan]back[/cyan]     Return to site editor")

        choice = prompt_with_context(
            self.processor,
            "\n[bold]Choice[/bold]",
            default="back",
            module="RedoxSite Editor - Residues",
            description="Edit residues in site"
        )
        return choice.strip().lower()

    def _list_residues(self):
        """List all residues in current site"""
        if not self.current_site.residue_groups:
            self.console.print("[yellow]No residues in this site[/yellow]")
            return

        table = Table(title=f"Residues in {self.current_site.site_id}")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Chain", style="green")
        table.add_column("Resname", style="green")
        table.add_column("ResID", style="yellow")
        table.add_column("Atoms", style="blue", justify="right")
        table.add_column("Is Center", style="magenta")

        # Get list of residues
        residues = []
        for (chain, resid, insertion), coords_list in self.current_site.residue_groups.items():
            # Get resname from first atom
            if coords_list:
                first_atom_info = self.current_site.coord_to_pdb.get(coords_list[0], {})
                resname = first_atom_info.get('resname', 'UNK')
            else:
                resname = 'UNK'

            # Check if this residue is a center
            is_center = any(
                c.chain == chain and c.resid == resid
                for c in self.current_site.centers
            )

            residues.append({
                'chain': chain,
                'resname': resname,
                'resid': resid,
                'insertion': insertion,
                'num_atoms': len(coords_list),
                'is_center': is_center
            })

        # Sort by chain, then resid
        residues.sort(key=lambda r: (r['chain'], r['resid']))

        for i, res in enumerate(residues, 1):
            table.add_row(
                str(i),
                res['chain'],
                res['resname'],
                str(res['resid']),
                str(res['num_atoms']),
                "Yes" if res['is_center'] else "No"
            )

        self.console.print(table)

    def _add_residue(self):
        """Add a residue from the structure file"""
        structure = self._get_structure()

        if not structure:
            self.console.print("[yellow]No structure available[/yellow]")
            self.console.print("[grey50]Workspace does not contain a valid PDB file[/grey50]")
            return

        self.console.print("\n[bold]Add Residue from Structure[/bold]")
        self.console.print("Enter residue details:\n")

        chain = prompt_with_context(
            self.processor,
            "Chain",
            module="RedoxSite Editor - Residues",
            description="Enter chain ID for residue to add"
        )
        resid = int(prompt_with_context(
            self.processor,
            "ResID",
            module="RedoxSite Editor - Residues",
            description="Enter residue ID to add"
        ))
        insertion_code = prompt_with_context(
            self.processor,
            "Insertion code",
            default="",
            module="RedoxSite Editor - Residues",
            description="Enter insertion code (optional)"
        )

        # Find residue in structure
        try:

            # Find the residue
            target_residue = None
            for model in structure:
                for chain_obj in model:
                    if chain_obj.id == chain:
                        for residue in chain_obj:
                            res_id = residue.id
                            if res_id[1] == resid and res_id[2].strip() == insertion_code.strip():
                                target_residue = residue
                                break
                    if target_residue:
                        break
                if target_residue:
                    break

            if not target_residue:
                self.console.print(f"[red]Residue not found in structure: {chain}:{resid}[/red]")
                return

            resname = target_residue.resname.strip()
            atoms = list(target_residue.get_atoms())

            self.console.print(f"\n[cyan]Found residue:[/cyan] {chain}:{resname}:{resid}")
            self.console.print(f"  Atoms: {len(atoms)}")

            if confirm_with_context(
                self.processor,
                "Add this residue to the site?",
                default=True,
                module="RedoxSite Editor - Residues",
                description="Confirm adding residue to site"
            ):
                from .comprehensive_redox_detector import RedoxSiteAtom

                # Add each atom
                added_count = 0
                for atom in atoms:
                    coords = tuple(atom.coord)

                    # Check if atom already exists
                    if coords in self.current_site.coord_to_pdb:
                        self.console.print(f"  [yellow]Skipping {atom.name} (already in site)[/yellow]")
                        continue

                    # Create RedoxSiteAtom
                    site_atom = RedoxSiteAtom(
                        chain=chain,
                        resname=resname,
                        resid=resid,
                        atom_name=atom.name.strip(),
                        coords=coords,
                        element=atom.element.strip(),
                        altloc=atom.altloc.strip(),
                        insertion_code=insertion_code,
                        occupancy=atom.occupancy,
                        bfactor=atom.bfactor
                    )

                    # Add to site
                    self.current_site.add_atom(site_atom)
                    added_count += 1

                self.console.print(f"[green]✓ Added {added_count} atom(s) from {chain}:{resname}:{resid}[/green]")
                self.console.print("[yellow]Note: No bonds added - use 'add bond' to create bonds[/yellow]")

        except ImportError:
            self.console.print("[red]BioPython not available - cannot load from PDB file[/red]")
        except Exception as e:
            self.console.print(f"[red]Error loading residue: {e}[/red]")

    def _remove_residue(self):
        """Remove a residue and associated bonds"""
        if not self.current_site.residue_groups:
            self.console.print("[yellow]No residues to remove[/yellow]")
            return

        self._list_residues()

        self.console.print("\n[grey50]Enter residue index to remove, or 'cancel'[/grey50]")
        choice = prompt_with_context(
            self.processor,
            "[bold]Index[/bold]",
            default="cancel",
            module="RedoxSite Editor - Residues",
            description="Select residue to remove"
        )

        if choice.lower() == 'cancel':
            return

        try:
            index = int(choice) - 1
            residues_list = sorted(self.current_site.residue_groups.keys(), key=lambda r: (r[0], r[1]))

            if 0 <= index < len(residues_list):
                res_key = residues_list[index]
                chain, resid, insertion = res_key
                coords_list = self.current_site.residue_groups[res_key]

                # Get resname
                first_atom_info = self.current_site.coord_to_pdb.get(coords_list[0], {})
                resname = first_atom_info.get('resname', 'UNK')

                # Show what will be removed
                self.console.print(f"\n[yellow]Remove residue:[/yellow] {chain}:{resname}:{resid}")
                self.console.print(f"  Atoms: {len(coords_list)}")

                # Find bonds involving this residue
                bonds_to_remove = []
                for bond in self.current_site.bonds:
                    if bond.atom1_coords in coords_list or bond.atom2_coords in coords_list:
                        bonds_to_remove.append(bond)

                if bonds_to_remove:
                    self.console.print(f"  [yellow]Bonds to remove: {len(bonds_to_remove)}[/yellow]")

                if confirm_with_context(
                    self.processor,
                    "Confirm removal?",
                    default=False,
                    module="RedoxSite Editor - Residues",
                    description="Confirm residue removal"
                ):
                    # Remove atoms
                    atoms_removed = 0
                    for coords in coords_list:
                        # Remove from atoms list
                        self.current_site.atoms = [a for a in self.current_site.atoms if a.coords != coords]
                        atoms_removed += 1

                        # Remove from coord_to_pdb
                        if coords in self.current_site.coord_to_pdb:
                            del self.current_site.coord_to_pdb[coords]

                    # Remove from residue_groups
                    del self.current_site.residue_groups[res_key]

                    # Remove bonds
                    for bond in bonds_to_remove:
                        self.current_site.bonds.remove(bond)

                    self.console.print(f"[green]✓ Removed residue {chain}:{resname}:{resid}[/green]")
                    self.console.print(f"  Atoms removed: {atoms_removed}")
                    self.console.print(f"  Bonds removed: {len(bonds_to_remove)}")
            else:
                self.console.print("[red]Invalid index[/red]")
        except ValueError:
            self.console.print("[red]Invalid input[/red]")

    def _replace_residue(self):
        """Replace a residue with another from the structure"""
        if not self.current_site.residue_groups:
            self.console.print("[yellow]No residues to replace[/yellow]")
            return

        self._list_residues()

        self.console.print("\n[grey50]Enter residue index to replace, or 'cancel'[/grey50]")
        choice = prompt_with_context(
            self.processor,
            "[bold]Index[/bold]",
            default="cancel",
            module="RedoxSite Editor - Residues",
            description="Select residue to replace"
        )

        if choice.lower() == 'cancel':
            return

        try:
            index = int(choice) - 1
            residues_list = sorted(self.current_site.residue_groups.keys(), key=lambda r: (r[0], r[1]))

            if 0 <= index < len(residues_list):
                old_res_key = residues_list[index]
                old_chain, old_resid, old_insertion = old_res_key
                old_coords_list = self.current_site.residue_groups[old_res_key]

                # Get old resname
                first_atom_info = self.current_site.coord_to_pdb.get(old_coords_list[0], {})
                old_resname = first_atom_info.get('resname', 'UNK')

                self.console.print(f"\n[yellow]Replacing:[/yellow] {old_chain}:{old_resname}:{old_resid}")

                # Get new residue info
                self.console.print("\n[bold]New residue:[/bold]")
                new_chain = prompt_with_context(
                    self.processor,
                    "Chain",
                    default=old_chain,
                    module="RedoxSite Editor - Residues",
                    description="Enter chain for replacement residue"
                )
                new_resname = prompt_with_context(
                    self.processor,
                    "Resname",
                    default=old_resname,
                    module="RedoxSite Editor - Residues",
                    description="Enter resname for replacement residue"
                )
                new_resid = prompt_with_context(
                    self.processor,
                    "ResID",
                    default=str(old_resid),
                    module="RedoxSite Editor - Residues",
                    description="Enter resid for replacement residue"
                )
                new_resid = int(new_resid)

                # Find bonds involving old residue
                bonds_involving_old = []
                for bond in self.current_site.bonds:
                    if bond.atom1_coords in old_coords_list or bond.atom2_coords in old_coords_list:
                        bonds_involving_old.append(bond)

                if bonds_involving_old:
                    self.console.print(f"\n[yellow]Found {len(bonds_involving_old)} bond(s) involving this residue[/yellow]")
                    self.console.print("[grey50]You will need to edit bonds after replacement[/grey50]")

                if confirm_with_context(
                    self.processor,
                    "\nProceed with replacement?",
                    default=False,
                    module="RedoxSite Editor - Residues",
                    description="Confirm residue replacement"
                ):
                    # Load new residue from structure
                    structure = self._get_structure()

                    if not structure:
                        self.console.print("\n[red]Cannot load structure from workspace[/red]")
                        self.console.print("[yellow]Falling back to metadata update only[/yellow]")

                        # Update metadata of existing atoms
                        for coords in old_coords_list:
                            # Update coord_to_pdb
                            if coords in self.current_site.coord_to_pdb:
                                self.current_site.coord_to_pdb[coords].update({
                                    'chain': new_chain,
                                    'resname': new_resname,
                                    'resid': new_resid
                                })

                            # Update atom objects
                            for atom in self.current_site.atoms:
                                if atom.coords == coords:
                                    atom.chain = new_chain
                                    atom.resname = new_resname
                                    atom.resid = new_resid
                                    break

                        # Update residue_groups
                        del self.current_site.residue_groups[old_res_key]
                        new_res_key = (new_chain, new_resid, old_insertion)
                        self.current_site.residue_groups[new_res_key] = old_coords_list

                        self.console.print(f"[green]✓ Updated residue metadata to {new_chain}:{new_resname}:{new_resid}[/green]")
                    else:
                        # Find the new residue in the structure
                        target_residue = None
                        for model in structure:
                            for chain_obj in model:
                                if chain_obj.id == new_chain:
                                    for residue in chain_obj:
                                        res_id = residue.id
                                        if res_id[1] == new_resid and res_id[2].strip() == old_insertion.strip():
                                            target_residue = residue
                                            break
                                if target_residue:
                                    break
                            if target_residue:
                                break

                        if not target_residue:
                            self.console.print(f"\n[red]Residue {new_chain}:{new_resname}:{new_resid} not found in structure[/red]")
                            return

                        # Remove old residue completely
                        from .comprehensive_redox_detector import RedoxSiteAtom

                        # Remove old atoms
                        for coords in old_coords_list:
                            self.current_site.atoms = [a for a in self.current_site.atoms if a.coords != coords]
                            if coords in self.current_site.coord_to_pdb:
                                del self.current_site.coord_to_pdb[coords]

                        # Remove old residue group
                        del self.current_site.residue_groups[old_res_key]

                        # Remove old bonds
                        for bond in bonds_involving_old:
                            self.current_site.bonds.remove(bond)

                        # Add new residue atoms
                        added_count = 0
                        for atom in target_residue.get_atoms():
                            coords = tuple(atom.coord)

                            # Create RedoxSiteAtom
                            site_atom = RedoxSiteAtom(
                                chain=new_chain,
                                resname=new_resname,
                                resid=new_resid,
                                atom_name=atom.name.strip(),
                                coords=coords,
                                element=atom.element.strip(),
                                altloc=atom.altloc.strip(),
                                insertion_code=old_insertion,
                                occupancy=atom.occupancy,
                                bfactor=atom.bfactor
                            )

                            # Add to site
                            self.current_site.add_atom(site_atom)
                            added_count += 1

                        self.console.print(f"[green]✓ Replaced {old_chain}:{old_resname}:{old_resid} with {new_chain}:{new_resname}:{new_resid}[/green]")
                        self.console.print(f"  Added {added_count} atoms from structure")
                        self.console.print(f"  Removed {len(bonds_involving_old)} bond(s)")

                    if bonds_involving_old:
                        self.console.print("\n[yellow]⚠ Remember to add bonds for the new residue[/yellow]")
            else:
                self.console.print("[red]Invalid index[/red]")
        except ValueError:
            self.console.print("[red]Invalid input[/red]")

    # ========================================================================
    # ATOMS EDITING
    # ========================================================================

    def _edit_atoms_menu(self):
        """Menu for editing individual atoms"""
        while True:
            choice = self._show_atoms_menu()

            if choice == 'back':
                break
            elif choice == 'list':
                self._list_atoms()
            elif choice == 'remove':
                self._remove_atom()
            elif choice == 'modify':
                self._modify_atom()

    def _show_atoms_menu(self) -> str:
        """Show atoms editing menu"""
        self.console.print("\n" + "-" * 60)
        self.console.print("[bold yellow]Edit Atoms[/bold yellow]")
        self.console.print("-" * 60 + "\n")

        self.console.print("[bold]Operations:[/bold]")
        self.console.print("  [cyan]list[/cyan]    List atoms (with filtering)")
        self.console.print("  [cyan]remove[/cyan]  Remove specific atom(s)")
        self.console.print("  [cyan]modify[/cyan]  Modify atom metadata")
        self.console.print("  [cyan]back[/cyan]    Return to site editor")

        choice = prompt_with_context(
            self.processor,
            "\n[bold]Choice[/bold]",
            default="back",
            module="RedoxSite Editor - Atoms",
            description="Edit individual atoms"
        )
        return choice.strip().lower()

    def _list_atoms(self):
        """List atoms with optional filtering"""
        if not self.current_site.atoms:
            self.console.print("[yellow]No atoms in this site[/yellow]")
            return

        self.console.print("\n[bold]Filter options:[/bold]")
        self.console.print("  [cyan]all[/cyan]       Show all atoms")
        self.console.print("  [cyan]residue[/cyan]   Filter by residue")
        self.console.print("  [cyan]element[/cyan]   Filter by element")

        filter_choice = prompt_with_context(
            self.processor,
            "[bold]Filter[/bold]",
            default="all",
            module="RedoxSite Editor - Atoms",
            description="Choose atom filter",
            options_map={'all': 'Show all atoms', 'residue': 'Filter by residue', 'element': 'Filter by element'}
        )

        atoms_to_show = self.current_site.atoms

        if filter_choice == 'residue':
            chain = prompt_with_context(
                self.processor,
                "Chain",
                module="RedoxSite Editor - Atoms",
                description="Filter atoms by chain"
            )
            resid = int(prompt_with_context(
                self.processor,
                "ResID",
                module="RedoxSite Editor - Atoms",
                description="Filter atoms by residue ID"
            ))
            atoms_to_show = [a for a in self.current_site.atoms if a.chain == chain and a.resid == resid]
        elif filter_choice == 'element':
            element = prompt_with_context(
                self.processor,
                "Element",
                module="RedoxSite Editor - Atoms",
                description="Filter atoms by element"
            ).upper()
            atoms_to_show = [a for a in self.current_site.atoms if a.element == element]

        if not atoms_to_show:
            self.console.print("[yellow]No atoms match filter[/yellow]")
            return

        # Show atoms table
        table = Table(title=f"Atoms in {self.current_site.site_id}")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Chain", style="green")
        table.add_column("Resname", style="green")
        table.add_column("ResID", style="yellow")
        table.add_column("Atom", style="blue")
        table.add_column("Element", style="magenta")

        for i, atom in enumerate(atoms_to_show, 1):
            table.add_row(
                str(i),
                atom.chain,
                atom.resname,
                str(atom.resid),
                atom.atom_name,
                atom.element
            )

        self.console.print(table)
        self.console.print(f"\n[grey50]Showing {len(atoms_to_show)} of {len(self.current_site.atoms)} total atoms[/grey50]")

    def _remove_atom(self):
        """Remove specific atom(s)"""
        if not self.current_site.atoms:
            self.console.print("[yellow]No atoms to remove[/yellow]")
            return

        # First, let user filter/select atoms
        self.console.print("\n[bold]Select atom to remove:[/bold]")
        self.console.print("Enter chain, resname, resid, and atom name")

        chain = prompt_with_context(
            self.processor,
            "Chain",
            module="RedoxSite Editor - Atoms",
            description="Enter chain for atom to remove"
        )
        resname = prompt_with_context(
            self.processor,
            "Resname",
            module="RedoxSite Editor - Atoms",
            description="Enter resname for atom to remove"
        )
        resid = int(prompt_with_context(
            self.processor,
            "ResID",
            module="RedoxSite Editor - Atoms",
            description="Enter resid for atom to remove"
        ))
        atom_name = prompt_with_context(
            self.processor,
            "Atom name",
            module="RedoxSite Editor - Atoms",
            description="Enter atom name to remove"
        )

        # Find matching atom
        matching_atom = None
        for atom in self.current_site.atoms:
            if (atom.chain == chain and atom.resname == resname and
                atom.resid == resid and atom.atom_name == atom_name):
                matching_atom = atom
                break

        if not matching_atom:
            self.console.print(f"[red]Atom not found: {chain}:{resname}:{resid} {atom_name}[/red]")
            return

        # Show what will be removed
        self.console.print(f"\n[yellow]Remove atom:[/yellow] {chain}:{resname}:{resid} {atom_name}")

        # Find bonds involving this atom
        bonds_involving_atom = []
        for bond in self.current_site.bonds:
            if bond.atom1_coords == matching_atom.coords or bond.atom2_coords == matching_atom.coords:
                bonds_involving_atom.append(bond)

        if bonds_involving_atom:
            self.console.print(f"  [yellow]Bonds to remove: {len(bonds_involving_atom)}[/yellow]")
            for bond in bonds_involving_atom:
                atom1_info = self.current_site.coord_to_pdb.get(bond.atom1_coords, {})
                atom2_info = self.current_site.coord_to_pdb.get(bond.atom2_coords, {})
                self.console.print(f"    • {atom1_info.get('resname','?')} {atom1_info.get('atom_name','?')} ↔ {atom2_info.get('resname','?')} {atom2_info.get('atom_name','?')}")

        if confirm_with_context(
            self.processor,
            "Confirm removal?",
            default=False,
            module="RedoxSite Editor - Atoms",
            description="Confirm atom removal"
        ):
            coords = matching_atom.coords

            # Remove from atoms list
            self.current_site.atoms.remove(matching_atom)

            # Remove from coord_to_pdb
            if coords in self.current_site.coord_to_pdb:
                del self.current_site.coord_to_pdb[coords]

            # Remove from residue_groups
            res_key = (chain, resid, matching_atom.insertion_code)
            if res_key in self.current_site.residue_groups:
                self.current_site.residue_groups[res_key].remove(coords)
                # If residue now empty, remove the key
                if not self.current_site.residue_groups[res_key]:
                    del self.current_site.residue_groups[res_key]

            # Remove bonds
            for bond in bonds_involving_atom:
                self.current_site.bonds.remove(bond)

            self.console.print(f"[green]✓ Removed atom {chain}:{resname}:{resid} {atom_name}[/green]")
            if bonds_involving_atom:
                self.console.print(f"  Bonds removed: {len(bonds_involving_atom)}")

    def _modify_atom(self):
        """Modify atom metadata"""
        if not self.current_site.atoms:
            self.console.print("[yellow]No atoms to modify[/yellow]")
            return

        # Select atom
        self.console.print("\n[bold]Select atom to modify:[/bold]")
        self.console.print("Enter chain, resname, resid, and atom name")

        chain = prompt_with_context(
            self.processor,
            "Chain",
            module="RedoxSite Editor - Atoms",
            description="Enter chain for atom to modify"
        )
        resname = prompt_with_context(
            self.processor,
            "Resname",
            module="RedoxSite Editor - Atoms",
            description="Enter resname for atom to modify"
        )
        resid = int(prompt_with_context(
            self.processor,
            "ResID",
            module="RedoxSite Editor - Atoms",
            description="Enter resid for atom to modify"
        ))
        atom_name = prompt_with_context(
            self.processor,
            "Atom name",
            module="RedoxSite Editor - Atoms",
            description="Enter atom name to modify"
        )

        # Find matching atom
        matching_atom = None
        for atom in self.current_site.atoms:
            if (atom.chain == chain and atom.resname == resname and
                atom.resid == resid and atom.atom_name == atom_name):
                matching_atom = atom
                break

        if not matching_atom:
            self.console.print(f"[red]Atom not found: {chain}:{resname}:{resid} {atom_name}[/red]")
            return

        # Modify properties
        self.console.print(f"\n[bold]Modifying atom: {chain}:{resname}:{resid} {atom_name}[/bold]")
        self.console.print("[grey50]Press Enter to keep current value[/grey50]\n")
        self.console.print("[yellow]Note: Coordinates remain unchanged (permanent identifier)[/yellow]\n")

        # New chain
        new_chain = prompt_with_context(
            self.processor,
            f"Chain [{matching_atom.chain}]",
            default=matching_atom.chain,
            module="RedoxSite Editor - Atoms",
            description="Modify atom chain"
        )

        # New resname
        new_resname = prompt_with_context(
            self.processor,
            f"Resname [{matching_atom.resname}]",
            default=matching_atom.resname,
            module="RedoxSite Editor - Atoms",
            description="Modify atom resname"
        )

        # New resid
        new_resid = prompt_with_context(
            self.processor,
            f"ResID [{matching_atom.resid}]",
            default=str(matching_atom.resid),
            module="RedoxSite Editor - Atoms",
            description="Modify atom resid"
        )
        new_resid = int(new_resid)

        # New atom name
        new_atom_name = prompt_with_context(
            self.processor,
            f"Atom name [{matching_atom.atom_name}]",
            default=matching_atom.atom_name,
            module="RedoxSite Editor - Atoms",
            description="Modify atom name"
        )

        # New element
        new_element = prompt_with_context(
            self.processor,
            f"Element [{matching_atom.element}]",
            default=matching_atom.element,
            module="RedoxSite Editor - Atoms",
            description="Modify atom element"
        )

        if confirm_with_context(
            self.processor,
            "\nApply changes?",
            default=True,
            module="RedoxSite Editor - Atoms",
            description="Confirm atom modifications"
        ):
            old_res_key = (matching_atom.chain, matching_atom.resid, matching_atom.insertion_code)

            # Update atom object
            matching_atom.chain = new_chain
            matching_atom.resname = new_resname
            matching_atom.resid = new_resid
            matching_atom.atom_name = new_atom_name
            matching_atom.element = new_element

            # Update coord_to_pdb
            if matching_atom.coords in self.current_site.coord_to_pdb:
                self.current_site.coord_to_pdb[matching_atom.coords].update({
                    'chain': new_chain,
                    'resname': new_resname,
                    'resid': new_resid,
                    'atom_name': new_atom_name,
                    'element': new_element
                })

            # Update residue_groups if residue changed
            new_res_key = (new_chain, new_resid, matching_atom.insertion_code)
            if old_res_key != new_res_key:
                # Remove from old group
                if old_res_key in self.current_site.residue_groups:
                    self.current_site.residue_groups[old_res_key].remove(matching_atom.coords)
                    if not self.current_site.residue_groups[old_res_key]:
                        del self.current_site.residue_groups[old_res_key]

                # Add to new group
                if new_res_key not in self.current_site.residue_groups:
                    self.current_site.residue_groups[new_res_key] = []
                self.current_site.residue_groups[new_res_key].append(matching_atom.coords)

            self.console.print(f"[green]✓ Atom updated to {new_chain}:{new_resname}:{new_resid} {new_atom_name}[/green]")

    # ========================================================================
    # BONDS EDITING
    # ========================================================================

    def _edit_bonds_menu(self):
        """Menu for editing bonds"""
        while True:
            choice = self._show_bonds_menu()

            if choice == 'back':
                break
            elif choice == 'list':
                self._list_bonds()
            elif choice == 'add':
                self._add_bond()
            elif choice == 'remove':
                self._remove_bond()
            elif choice == 'reclassify':
                self._reclassify_bond()

    def _show_bonds_menu(self) -> str:
        """Show bonds editing menu"""
        self.console.print("\n" + "-" * 60)
        self.console.print("[bold yellow]Edit Bonds[/bold yellow]")
        self.console.print("-" * 60 + "\n")

        self.console.print("[bold]Operations:[/bold]")
        self.console.print("  [cyan]list[/cyan]        List all bonds")
        self.console.print("  [cyan]add[/cyan]         Add a new bond")
        self.console.print("  [cyan]remove[/cyan]      Remove a bond")
        self.console.print("  [cyan]reclassify[/cyan]  Change bond type")
        self.console.print("  [cyan]back[/cyan]        Return to site editor")

        choice = prompt_with_context(
            self.processor,
            "\n[bold]Choice[/bold]",
            default="back",
            module="RedoxSite Editor - Bonds",
            description="Edit bonds in site"
        )
        return choice.strip().lower()

    def _list_bonds(self):
        """List all bonds in current site"""
        if not self.current_site.bonds:
            self.console.print("[yellow]No bonds in this site[/yellow]")
            return

        # Group bonds by chemical type
        from collections import defaultdict
        bonds_by_type = defaultdict(list)
        for bond in self.current_site.bonds:
            bonds_by_type[bond.chemical_type].append(bond)

        for chem_type, bonds in bonds_by_type.items():
            table = Table(title=f"{chem_type.capitalize()} Bonds ({len(bonds)})")
            table.add_column("Index", style="cyan", width=6)
            table.add_column("Atom 1", style="green")
            table.add_column("Atom 2", style="yellow")
            table.add_column("Distance", style="blue", justify="right")
            table.add_column("Bond Type", style="magenta")

            for i, bond in enumerate(bonds, 1):
                # Get atom info
                atom1_info = self.current_site.coord_to_pdb.get(bond.atom1_coords, {})
                atom2_info = self.current_site.coord_to_pdb.get(bond.atom2_coords, {})

                atom1_str = f"{atom1_info.get('chain', '?')}:{atom1_info.get('resname', '?')}:{atom1_info.get('resid', '?')} {atom1_info.get('atom_name', '?')}"
                atom2_str = f"{atom2_info.get('chain', '?')}:{atom2_info.get('resname', '?')}:{atom2_info.get('resid', '?')} {atom2_info.get('atom_name', '?')}"

                table.add_row(
                    str(i),
                    atom1_str,
                    atom2_str,
                    f"{bond.distance:.2f}Å",
                    bond.bond_type
                )

            self.console.print(table)
            self.console.print()

    def _remove_bond(self):
        """Remove one or more bonds by index."""
        if not self.current_site.bonds:
            self.console.print("[yellow]No bonds to remove[/yellow]")
            return

        # Display flat numbered list of all bonds
        all_bonds = list(self.current_site.bonds)
        table = Table(title=f"All Bonds ({len(all_bonds)})")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Atom 1", style="green")
        table.add_column("Atom 2", style="yellow")
        table.add_column("Distance", style="blue", justify="right")
        table.add_column("Type", style="magenta")

        for i, bond in enumerate(all_bonds, 1):
            atom1_info = self.current_site.coord_to_pdb.get(bond.atom1_coords, {})
            atom2_info = self.current_site.coord_to_pdb.get(bond.atom2_coords, {})
            atom1_str = f"{atom1_info.get('chain', '?')}:{atom1_info.get('resname', '?')}:{atom1_info.get('resid', '?')} {atom1_info.get('atom_name', '?')}"
            atom2_str = f"{atom2_info.get('chain', '?')}:{atom2_info.get('resname', '?')}:{atom2_info.get('resid', '?')} {atom2_info.get('atom_name', '?')}"
            table.add_row(str(i), atom1_str, atom2_str, f"{bond.distance:.2f}Å", bond.chemical_type)

        self.console.print(table)

        self.console.print("\n[grey50]Enter bond number(s) to remove (e.g., '2' or '1,3,4'), or 'cancel'[/grey50]")
        choice = prompt_with_context(
            self.processor,
            "[bold]Bond[/bold]",
            default="cancel",
            module="RedoxSite Editor - Bonds",
            description="Select bond(s) to remove"
        )

        if choice.lower() == 'cancel':
            return

        try:
            indices = [int(x.strip()) - 1 for x in choice.split(',')]
        except ValueError:
            self.console.print("[red]Invalid input. Enter number(s) separated by commas.[/red]")
            return

        # Validate all indices first
        invalid = [i + 1 for i in indices if i < 0 or i >= len(all_bonds)]
        if invalid:
            self.console.print(f"[red]Invalid index(es): {', '.join(str(x) for x in invalid)}. Valid range: 1-{len(all_bonds)}[/red]")
            return

        # Show bonds to remove
        bonds_to_remove = [all_bonds[i] for i in indices]
        for bond in bonds_to_remove:
            atom1_info = self.current_site.coord_to_pdb.get(bond.atom1_coords, {})
            atom2_info = self.current_site.coord_to_pdb.get(bond.atom2_coords, {})
            atom1_str = f"{atom1_info.get('chain', '?')}:{atom1_info.get('resname', '?')}:{atom1_info.get('resid', '?')} {atom1_info.get('atom_name', '?')}"
            atom2_str = f"{atom2_info.get('chain', '?')}:{atom2_info.get('resname', '?')}:{atom2_info.get('resid', '?')} {atom2_info.get('atom_name', '?')}"
            self.console.print(f"  [yellow]Remove:[/yellow] {atom1_str} ↔ {atom2_str}")

        if confirm_with_context(
            self.processor,
            f"Remove {len(bonds_to_remove)} bond(s)?",
            default=False,
            module="RedoxSite Editor - Bonds",
            description="Confirm bond removal"
        ):
            for bond in bonds_to_remove:
                self.current_site.bonds.remove(bond)
            self.console.print(f"[green]✓ {len(bonds_to_remove)} bond(s) removed[/green]")

    def _add_bond(self):
        """Add a new bond between two atoms"""
        if not self.current_site.atoms:
            self.console.print("[yellow]No atoms in site to bond[/yellow]")
            return

        self.console.print("\n[bold]Add Bond[/bold]")
        self.console.print("Enter details for the two atoms to bond:\n")

        # Get first atom
        self.console.print("[cyan]Atom 1:[/cyan]")
        chain1 = prompt_with_context(
            self.processor,
            "  Chain",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 1 chain"
        )
        resname1 = prompt_with_context(
            self.processor,
            "  Resname",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 1 resname"
        )
        resid1 = int(prompt_with_context(
            self.processor,
            "  ResID",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 1 resid"
        ))
        atom_name1 = prompt_with_context(
            self.processor,
            "  Atom name",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 1 atom name"
        )

        # Find first atom
        atom1 = None
        for atom in self.current_site.atoms:
            if (atom.chain == chain1 and atom.resname == resname1 and
                atom.resid == resid1 and atom.atom_name == atom_name1):
                atom1 = atom
                break

        if not atom1:
            self.console.print(f"[red]Atom 1 not found: {chain1}:{resname1}:{resid1} {atom_name1}[/red]")
            return

        # Get second atom
        self.console.print("\n[cyan]Atom 2:[/cyan]")
        chain2 = prompt_with_context(
            self.processor,
            "  Chain",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 2 chain"
        )
        resname2 = prompt_with_context(
            self.processor,
            "  Resname",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 2 resname"
        )
        resid2 = int(prompt_with_context(
            self.processor,
            "  ResID",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 2 resid"
        ))
        atom_name2 = prompt_with_context(
            self.processor,
            "  Atom name",
            module="RedoxSite Editor - Bonds",
            description="Add bond - specify atom 2 atom name"
        )

        # Find second atom
        atom2 = None
        for atom in self.current_site.atoms:
            if (atom.chain == chain2 and atom.resname == resname2 and
                atom.resid == resid2 and atom.atom_name == atom_name2):
                atom2 = atom
                break

        if not atom2:
            self.console.print(f"[red]Atom 2 not found: {chain2}:{resname2}:{resid2} {atom_name2}[/red]")
            return

        # Calculate distance
        import math
        coords1 = atom1.coords
        coords2 = atom2.coords
        distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(coords1, coords2)))

        self.console.print(f"\n[yellow]Distance between atoms: {distance:.2f}Å[/yellow]")

        # Get bond type
        self.console.print("\n[bold]Bond type:[/bold]")
        self.console.print("  1. covalent")
        self.console.print("  2. coordinate")
        self.console.print("  3. disulfide")
        self.console.print("  4. metal-metal")

        bond_type_choice = prompt_with_context(
            self.processor,
            "Choice",
            choices=['1', '2', '3', '4'],
            default='1',
            module="RedoxSite Editor - Bonds",
            description="Add bond - select bond type"
        )

        bond_type_map = {
            '1': 'covalent',
            '2': 'coordinate',
            '3': 'disulfide',
            '4': 'metal-metal'
        }
        chemical_type = bond_type_map[bond_type_choice]

        # Determine bond_type (intra vs inter-residue)
        if chain1 == chain2 and resid1 == resid2 and resname1 == resname2:
            bond_type = "intraresidue"
        else:
            bond_type = "interresidue"

        if confirm_with_context(
            self.processor,
            f"\nAdd {chemical_type} bond ({distance:.2f}Å)?",
            default=True,
            module="RedoxSite Editor - Bonds",
            description="Confirm adding bond"
        ):
            # Create bond using RedoxSiteBond dataclass
            from .comprehensive_redox_detector import RedoxSiteBond

            bond = RedoxSiteBond(
                atom1_coords=coords1,
                atom2_coords=coords2,
                bond_type=bond_type,
                chemical_type=chemical_type,
                distance=distance,
                atom1_element=atom1.element,
                atom2_element=atom2.element,
                atom1_residue_info=self.current_site.coord_to_pdb.get(coords1, {}),
                atom2_residue_info=self.current_site.coord_to_pdb.get(coords2, {})
            )

            self.current_site.bonds.append(bond)
            self.console.print(f"[green]✓ Bond added: {chain1}:{resname1}:{resid1} {atom_name1} ↔ {chain2}:{resname2}:{resid2} {atom_name2}[/green]")

    def _reclassify_bond(self):
        """Change bond chemical type"""
        if not self.current_site.bonds:
            self.console.print("[yellow]No bonds to reclassify[/yellow]")
            return

        self._list_bonds()

        self.console.print("\n[grey50]Enter bond index to reclassify (use chemical type + index, e.g., 'coordinate 1'), or 'cancel'[/grey50]")
        choice = prompt_with_context(
            self.processor,
            "[bold]Bond[/bold]",
            default="cancel",
            module="RedoxSite Editor - Bonds",
            description="Reclassify bond - select bond to reclassify"
        )

        if choice.lower() == 'cancel':
            return

        # Parse input
        parts = choice.split()
        if len(parts) != 2:
            self.console.print("[red]Invalid format. Use: <type> <index>[/red]")
            return

        chem_type, index_str = parts

        try:
            index = int(index_str) - 1

            # Find bonds of this type
            bonds_of_type = [b for b in self.current_site.bonds if b.chemical_type == chem_type]

            if 0 <= index < len(bonds_of_type):
                bond = bonds_of_type[index]

                # Get atom info for display
                atom1_info = self.current_site.coord_to_pdb.get(bond.atom1_coords, {})
                atom2_info = self.current_site.coord_to_pdb.get(bond.atom2_coords, {})

                atom1_str = f"{atom1_info.get('chain', '?')}:{atom1_info.get('resname', '?')}:{atom1_info.get('resid', '?')} {atom1_info.get('atom_name', '?')}"
                atom2_str = f"{atom2_info.get('chain', '?')}:{atom2_info.get('resname', '?')}:{atom2_info.get('resid', '?')} {atom2_info.get('atom_name', '?')}"

                self.console.print(f"\n[yellow]Reclassify bond:[/yellow] {atom1_str} ↔ {atom2_str}")
                self.console.print(f"Current type: [cyan]{bond.chemical_type}[/cyan]")
                self.console.print(f"Distance: [cyan]{bond.distance:.2f}Å[/cyan]")

                # Get new chemical type
                self.console.print("\n[bold]New bond type:[/bold]")
                self.console.print("  1. covalent")
                self.console.print("  2. coordinate")
                self.console.print("  3. disulfide")
                self.console.print("  4. metal-metal")

                type_choice = prompt_with_context(
                    self.processor,
                    "Choice",
                    choices=['1', '2', '3', '4'],
                    module="RedoxSite Editor - Bonds",
                    description="Reclassify bond - select new bond type"
                )

                type_map = {
                    '1': 'covalent',
                    '2': 'coordinate',
                    '3': 'disulfide',
                    '4': 'metal-metal'
                }
                new_chemical_type = type_map[type_choice]

                if confirm_with_context(
                    self.processor,
                    f"\nChange bond type to {new_chemical_type}?",
                    default=True,
                    module="RedoxSite Editor - Bonds",
                    description="Confirm bond reclassification"
                ):
                    bond.chemical_type = new_chemical_type
                    self.console.print(f"[green]✓ Bond reclassified to {new_chemical_type}[/green]")
            else:
                self.console.print("[red]Invalid index[/red]")
        except ValueError:
            self.console.print("[red]Invalid input[/red]")

    # ========================================================================
    # METADATA EDITING
    # ========================================================================

    def _edit_metadata_menu(self):
        """Menu for editing site metadata"""
        self.console.print("\n" + "-" * 60)
        self.console.print("[bold yellow]Edit Metadata[/bold yellow]")
        self.console.print("-" * 60 + "\n")

        self.console.print(f"[bold]Current metadata for {self.current_site.site_id}:[/bold]\n")
        self.console.print(f"  Site Type: [cyan]{self.current_site.site_type if self.current_site.site_type else 'N/A'}[/cyan]")
        self.console.print(f"  Detection Method: [cyan]{self.current_site.detection_method if self.current_site.detection_method else 'N/A'}[/cyan]")
        self.console.print(f"  Boundary Definition: [cyan]{self.current_site.boundary_definition if self.current_site.boundary_definition else 'N/A'}[/cyan]")

        self.console.print("\n[grey50]Press Enter to keep current value[/grey50]\n")

        # Edit site type
        new_site_type = prompt_with_context(
            self.processor,
            f"Site Type [{self.current_site.site_type if self.current_site.site_type else 'N/A'}]",
            default=self.current_site.site_type if self.current_site.site_type else "",
            module="RedoxSite Editor - Metadata",
            description="Edit site type"
        )

        if new_site_type:
            self.current_site.site_type = new_site_type
            self.console.print("[green]✓ Site type updated[/green]")

    # ========================================================================
    # VALIDATION
    # ========================================================================

    def _validate_site(self, site):
        """Validate site consistency"""
        self.console.print("\n" + "-" * 60)
        self.console.print(f"[bold yellow]Validating {site.site_id}[/bold yellow]")
        self.console.print("-" * 60 + "\n")

        issues = []

        # Check 1: All atoms in atoms list have coords in coord_to_pdb
        for atom in site.atoms:
            if atom.coords not in site.coord_to_pdb:
                issues.append(f"Atom {atom.chain}:{atom.resname}:{atom.resid} {atom.atom_name} missing from coord_to_pdb")

        # Check 2: All centers have coords in coord_to_pdb
        for center in site.centers:
            if center.coords not in site.coord_to_pdb:
                issues.append(f"Center {center.chain}:{center.resname}:{center.resid} missing from coord_to_pdb")

        # Check 3: All bonds reference atoms in coord_to_pdb
        for i, bond in enumerate(site.bonds):
            if bond.atom1_coords not in site.coord_to_pdb:
                issues.append(f"Bond {i+1} atom1 coords not in site")
            if bond.atom2_coords not in site.coord_to_pdb:
                issues.append(f"Bond {i+1} atom2 coords not in site")

        # Check 4: residue_groups consistency
        for res_key, coords_list in site.residue_groups.items():
            for coords in coords_list:
                if coords not in site.coord_to_pdb:
                    chain, resid, insertion = res_key
                    issues.append(f"Residue group {chain}:{resid} contains coords not in coord_to_pdb")
                    break

        # Report results
        if issues:
            self.console.print(f"[red]✗ Found {len(issues)} issue(s):[/red]\n")
            for issue in issues:
                self.console.print(f"  • {issue}")
        else:
            self.console.print("[green]✓ Site validation passed - no issues found[/green]")

        return len(issues) == 0

    def _validate_all_sites(self):
        """Validate all sites"""
        self.console.print("\n[bold cyan]Validating All Sites[/bold cyan]\n")

        total_issues = 0
        for site in self.working_sites:
            is_valid = self._validate_site(site)
            if not is_valid:
                total_issues += 1
            self.console.print()

        if total_issues == 0:
            self.console.print("[green]✓ All sites passed validation[/green]")
        else:
            self.console.print(f"[red]✗ {total_issues} site(s) have validation issues[/red]")

    # ========================================================================
    # SUMMARY AND UTILITIES
    # ========================================================================

    def _show_site_summary(self, site):
        """Show detailed summary of site"""
        panel_content = []

        panel_content.append(f"[bold]Site ID:[/bold] {site.site_id}")
        panel_content.append(f"[bold]Structure:[/bold] {site.structure_id}")
        panel_content.append(f"[bold]Site Type:[/bold] {site.site_type if site.site_type else 'N/A'}")
        panel_content.append("")

        panel_content.append(f"[bold cyan]Centers ({len(site.centers)}):[/bold cyan]")
        for center in site.centers:
            center_type = center.center_type.value if hasattr(center.center_type, 'value') else str(center.center_type)
            panel_content.append(f"  • {center.chain}:{center.resname}:{center.resid} ({center_type})")

        panel_content.append("")
        panel_content.append(f"[bold cyan]Residues ({len(site.residue_groups)}):[/bold cyan]")
        for (chain, resid, insertion), coords_list in sorted(site.residue_groups.items(), key=lambda x: (x[0][0], x[0][1])):
            first_atom_info = site.coord_to_pdb.get(coords_list[0], {})
            resname = first_atom_info.get('resname', 'UNK')
            panel_content.append(f"  • {chain}:{resname}:{resid} ({len(coords_list)} atoms)")

        panel_content.append("")
        panel_content.append(f"[bold cyan]Statistics:[/bold cyan]")
        panel_content.append(f"  Total Atoms: {len(site.atoms)}")
        panel_content.append(f"  Total Bonds: {len(site.bonds)}")

        # Bond breakdown
        from collections import Counter
        bond_types = Counter(b.chemical_type for b in site.bonds)
        for bond_type, count in bond_types.items():
            panel_content.append(f"    {bond_type}: {count}")

        panel = Panel(
            "\n".join(panel_content),
            title=f"[bold]Site Summary: {site.site_id}[/bold]",
            border_style="cyan"
        )

        self.console.print(panel)

    def _confirm_cancel(self) -> bool:
        """Confirm cancellation and discarding changes"""
        self.console.print("\n[yellow]⚠ You have unsaved changes[/yellow]")
        return confirm_with_context(
            self.processor,
            "Discard all changes and exit?",
            default=False,
            module="RedoxSite Editor",
            description="Confirm discarding changes"
        )
