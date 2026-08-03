"""
Shared Bond Management Base Class

Provides common functionality for tLEaP bond command generation and management
across different ProPrep modules (disulfide, redox sites, etc.).
"""

import logging
from typing import Dict, List, Tuple, Any, Optional
from rich.console import Console
from rich.table import Table
from pathlib import Path

from proprep.utils.prompts import prompt_with_context, confirm_with_context
import re

logger = logging.getLogger(__name__)


class BondManagementBase:
    """Base class for tLEaP bond command generation and management."""
    
    def __init__(self, processor, module_name: str):
        self.processor = processor
        self.module_name = module_name
        self.console = processor.console if processor else Console()
        self.bonds_by_category = {}
        self.original_bonds = {}
        
    def initialize_bond_categories(self, categories: Dict[str, str]):
        """Initialize bond categories with descriptions.
        
        Args:
            categories: Dict mapping category keys to descriptions
                       e.g., {"disulfide_bonds": "Inter-residue Disulfide Bonds"}
        """
        self.bond_categories = categories
        self.bonds_by_category = {cat: [] for cat in categories.keys()}
        
    def add_bonds_to_category(self, category: str, bonds: List[str]):
        """Add bonds to a specific category."""
        if category not in self.bonds_by_category:
            self.bonds_by_category[category] = []
        self.bonds_by_category[category].extend(bonds)
        
    def backup_original_bonds(self):
        """Create a backup of the original bonds for reset functionality."""
        self.original_bonds = {
            cat: bonds.copy() for cat, bonds in self.bonds_by_category.items()
        }
        
    def get_total_bond_count(self) -> int:
        """Get total number of bonds across all categories."""
        return sum(len(bonds) for bonds in self.bonds_by_category.values())
        
    def display_bond_summary(self):
        """Display summary of bonds by category."""
        total_bonds = self.get_total_bond_count()
        
        if total_bonds == 0:
            self.console.print("[yellow]No bond commands generated.[/yellow]")
            return
            
        self.console.print(f"\n[bold]Generated {total_bonds} bond commands across {len([c for c in self.bonds_by_category.values() if c])} categories:[/bold]")
        
        for category, bonds in self.bonds_by_category.items():
            if bonds:
                category_name = self.bond_categories.get(category, category)
                self.console.print(f"• {category_name}: {len(bonds)} commands")
                
                
    def _review_and_edit_bonds(self):
        """Review and edit all bond commands."""
        total_bonds = self.get_total_bond_count()
        if total_bonds == 0:
            self.console.print("[yellow]No bonds to review.[/yellow]")
            return
            
        # Create numbered list of all bonds with categories
        all_bonds = []
        for category, bonds in self.bonds_by_category.items():
            for bond in bonds:
                all_bonds.append((category, bond))
                
        # Display all bonds in a table
        table = Table(title="All Bond Commands")
        table.add_column("Index", style="cyan", width=6)
        table.add_column("Category", style="magenta", width=20)
        table.add_column("Bond Command", style="green")
        table.add_column("Comment", style="yellow")
        
        for i, (category, bond) in enumerate(all_bonds, 1):
            category_name = self.bond_categories.get(category, category)
            parts = bond.split("#", 1)
            command = parts[0].strip()
            comment = parts[1].strip() if len(parts) > 1 else ""
            table.add_row(str(i), category_name, command, comment)
            
        self.console.print(table)
        
        while True:
            action = prompt_with_context(
                self.processor,
                "Action",
                choices=["e", "d", "b"],
                default="b",
                module=f"{self.module_name} Bond Management",
                description="Bond review action",
                options_map={"e": "Edit specific bond", "d": "Delete specific bond", "b": "Back"},
            ).lower()

            if action == "e":
                self._edit_specific_bond(all_bonds)
            elif action == "d":
                self._delete_specific_bond(all_bonds)
            elif action == "b":
                break
                
    def _edit_specific_bond(self, all_bonds: List[Tuple[str, str]]):
        """Edit a specific bond command."""
        if not all_bonds:
            return
            
        try:
            index = int(prompt_with_context(
                self.processor,
                f"Enter bond index to edit [1-{len(all_bonds)}]",
                module=f"{self.module_name} Bond Management",
                description="Bond index to edit",
            )) - 1
            if 0 <= index < len(all_bonds):
                category, old_bond = all_bonds[index]

                # Show current command
                parts = old_bond.split("#", 1)
                current_command = parts[0].strip()
                current_comment = parts[1].strip() if len(parts) > 1 else ""

                self.console.print(f"Current command: [green]{current_command}[/green]")
                self.console.print(f"Current comment: [yellow]{current_comment}[/yellow]")

                # Get new command
                new_command = prompt_with_context(
                    self.processor,
                    "Enter new command",
                    default=current_command,
                    module=f"{self.module_name} Bond Management",
                    description="New bond command text",
                )
                new_comment = prompt_with_context(
                    self.processor,
                    "Enter new comment",
                    default=current_comment,
                    module=f"{self.module_name} Bond Management",
                    description="New bond comment text",
                )
                
                # Validate command
                if self._validate_bond_command(new_command):
                    new_bond = f"{new_command} # {new_comment}" if new_comment else new_command
                    
                    # Update the bond in the appropriate category
                    old_index = self.bonds_by_category[category].index(old_bond)
                    self.bonds_by_category[category][old_index] = new_bond
                    
                    self.console.print("[green]Bond updated successfully.[/green]")
                else:
                    self.console.print("[red]Invalid bond command format.[/red]")
            else:
                self.console.print("[red]Invalid index.[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number.[/red]")
            
    def _delete_specific_bond(self, all_bonds: List[Tuple[str, str]]):
        """Delete a specific bond command."""
        if not all_bonds:
            return
            
        try:
            index = int(prompt_with_context(
                self.processor,
                f"Enter bond index to delete [1-{len(all_bonds)}]",
                module=f"{self.module_name} Bond Management",
                description="Bond index to delete",
            )) - 1
            if 0 <= index < len(all_bonds):
                category, bond_to_delete = all_bonds[index]

                if confirm_with_context(
                    self.processor,
                    f"Delete bond: {bond_to_delete}",
                    module=f"{self.module_name} Bond Management",
                    description="Confirm delete bond",
                ):
                    self.bonds_by_category[category].remove(bond_to_delete)
                    self.console.print("[green]Bond deleted successfully.[/green]")
            else:
                self.console.print("[red]Invalid index.[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number.[/red]")
            
    def _add_new_bond(self):
        """Add a new bond command."""
        # Show available categories
        categories = [(k, v) for k, v in self.bond_categories.items()]
        
        self.console.print("\n[bold]Available Categories:[/bold]")
        for i, (key, desc) in enumerate(categories, 1):
            self.console.print(f"{i}. {desc}")
            
        try:
            cat_options_map = {str(i + 1): desc for i, (_, desc) in enumerate(categories)}
            choice = int(prompt_with_context(
                self.processor,
                f"Select category [1-{len(categories)}]",
                module=f"{self.module_name} Bond Management",
                description="Select bond category",
                options_map=cat_options_map,
            )) - 1
            if 0 <= choice < len(categories):
                category_key = categories[choice][0]

                # Get bond command from user
                command = prompt_with_context(
                    self.processor,
                    "Enter bond command (e.g., 'bond mol.42.SG mol.156.SG')",
                    module=f"{self.module_name} Bond Management",
                    description="Enter bond command",
                )

                if self._validate_bond_command(command):
                    comment = prompt_with_context(
                        self.processor,
                        "Enter comment (optional)",
                        default="",
                        module=f"{self.module_name} Bond Management",
                        description="Bond comment",
                    )
                    new_bond = f"{command} # {comment}" if comment else command
                    
                    # Check for duplicates
                    if new_bond not in self.bonds_by_category[category_key]:
                        self.bonds_by_category[category_key].append(new_bond)
                        self.console.print("[green]Bond added successfully.[/green]")
                    else:
                        self.console.print("[yellow]Bond already exists.[/yellow]")
                else:
                    self.console.print("[red]Invalid bond command format.[/red]")
            else:
                self.console.print("[red]Invalid category selection.[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number.[/red]")
            
    def _remove_bonds(self):
        """Remove bonds by category or individually."""
        # Show categories with bond counts
        categories_with_bonds = [(k, v, len(self.bonds_by_category[k])) 
                                for k, v in self.bond_categories.items() 
                                if self.bonds_by_category[k]]
        
        if not categories_with_bonds:
            self.console.print("[yellow]No bonds to remove.[/yellow]")
            return
            
        self.console.print("\n[bold]Categories with bonds:[/bold]")
        for i, (key, desc, count) in enumerate(categories_with_bonds, 1):
            self.console.print(f"{i}. {desc} ({count} bonds)")
            
        try:
            rm_options_map = {
                str(i + 1): f"{desc} ({count} bonds)"
                for i, (_, desc, count) in enumerate(categories_with_bonds)
            }
            choice = int(prompt_with_context(
                self.processor,
                f"Select category to remove from [1-{len(categories_with_bonds)}]",
                module=f"{self.module_name} Bond Management",
                description="Select bond category to clear",
                options_map=rm_options_map,
            )) - 1
            if 0 <= choice < len(categories_with_bonds):
                category_key = categories_with_bonds[choice][0]

                if confirm_with_context(
                    self.processor,
                    f"Remove all bonds from {categories_with_bonds[choice][1]}?",
                    module=f"{self.module_name} Bond Management",
                    description=f"Confirm remove all bonds from {categories_with_bonds[choice][1]}",
                ):
                    self.bonds_by_category[category_key].clear()
                    self.console.print("[green]Bonds removed successfully.[/green]")
            else:
                self.console.print("[red]Invalid category selection.[/red]")
        except ValueError:
            self.console.print("[red]Please enter a valid number.[/red]")
            
    def _view_bonds_by_category(self):
        """View bonds organized by category."""
        for category, bonds in self.bonds_by_category.items():
            if not bonds:
                continue
                
            category_name = self.bond_categories.get(category, category)
            table = Table(title=f"{category_name} ({len(bonds)} commands)")
            table.add_column("Index", style="cyan", width=6)
            table.add_column("Bond Command", style="green")
            table.add_column("Comment", style="yellow")
            
            for i, bond in enumerate(bonds, 1):
                parts = bond.split("#", 1)
                command = parts[0].strip()
                comment = parts[1].strip() if len(parts) > 1 else ""
                table.add_row(str(i), command, comment)
                
            self.console.print(table)
            
        prompt_with_context(
            self.processor,
            "Press Enter to continue",
            default="",
            module=f"{self.module_name} Bond Management",
            description="Pause after bond view",
        )

    def _reset_bonds(self):
        """Reset bonds to original generated state."""
        if confirm_with_context(
            self.processor,
            "Reset to original generated bonds?",
            module=f"{self.module_name} Bond Management",
            description="Confirm reset bonds to original",
        ):
            self.bonds_by_category = {
                cat: bonds.copy() for cat, bonds in self.original_bonds.items()
            }
            self.console.print("[green]Bonds reset to original state.[/green]")
            
    def _validate_bond_command(self, command: str) -> bool:
        """Validate tLEaP bond command format."""
        # Basic validation - should start with 'bond' and have proper format
        pattern = r'^bond\s+mol\.\w+\.\w+\s+mol\.\w+\.\w+.*$'
        return bool(re.match(pattern, command.strip(), re.IGNORECASE))
        
    def save_bonds_to_file(self, filename: str, pdb_name: str = "unknown") -> bool:
        """Save bonds to file organized by category."""
        try:
            filepath = Path(filename)
            
            with open(filepath, 'w') as f:
                f.write(f"# {self.module_name} Bond Commands Generated by ProPrep\n")
                f.write(f"# PDB: {pdb_name}\n")
                f.write(f"# Generated: {self._get_timestamp()}\n\n")
                
                total_bonds = 0
                for category, bonds in self.bonds_by_category.items():
                    if not bonds:
                        continue
                        
                    category_name = self.bond_categories.get(category, category)
                    f.write(f"#=== {category_name} ===\n")
                    
                    for bond in bonds:
                        f.write(f"{bond}\n")
                        total_bonds += 1
                    f.write("\n")
                    
                f.write(f"# Total: {total_bonds} bond commands generated\n")
                
            self.console.print(f"[green]Bond commands saved to {filepath}[/green]")
            return True
            
        except Exception as e:
            self.console.print(f"[red]Error saving file: {e}[/red]")
            return False
            
    def _get_timestamp(self) -> str:
        """Get current timestamp for file headers."""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
    # Note: Removed get_intra_residue_bonds_from_tleap_generator method
    # The correct flow is: bond generators → workspace → Topology Generator gathers bonds
        
    def _is_intra_residue_bond(self, bond_command: str) -> bool:
        """Check if a bond command represents an intra-residue bond."""
        # Parse bond command to extract residue identifiers
        # Format: bond mol.resid1.atom1 mol.resid2.atom2
        pattern = r'bond\s+mol\.(\w+)\.(\w+)\s+mol\.(\w+)\.(\w+)'
        match = re.match(pattern, bond_command.strip(), re.IGNORECASE)
        
        if match:
            resid1, atom1, resid2, atom2 = match.groups()
            # Same residue ID indicates intra-residue bond
            return resid1 == resid2
            
        return False