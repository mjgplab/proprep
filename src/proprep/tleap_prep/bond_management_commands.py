"""
Bond Management Menu Commands

Provides MenuCommand-based bond management interface that follows 
standard ProPrep navigation patterns.
"""

from typing import Dict, List, Tuple, Any
from proprep.utils.prompts import prompt_with_context, confirm_with_context
from rich.table import Table

from proprep.application.processor_command import MenuCommand, ProcessorCommand


class BondManagementMenuCommand(MenuCommand):
    """MenuCommand for bond management that follows standard ProPrep patterns."""
    
    def __init__(self, processor, bond_manager, menu_category: str = None):
        super().__init__(
            processor, 
            f"{bond_manager.module_name} Bond Management",
            menu_category
        )
        self.bond_manager = bond_manager
        
    def execute(self) -> None:
        """Execute the bond management menu using standard MenuCommand pattern."""
        self.enter_menu()
        
        try:
            while True:
                self.console.print()
                self.bond_manager.display_bond_summary()
                
                self.show_breadcrumbs()
                self.console.print(f"\n[bold]===== {self.bond_manager.module_name} Bond Management =====[/bold]")
                
                # Show menu options in standard format
                self.console.print("1. Review and edit bond commands", highlight=False)
                self.console.print("2. Add new bond command", highlight=False)
                self.console.print("3. Remove bond commands", highlight=False)
                self.console.print("4. View bonds by category", highlight=False)
                self.console.print("5. Reset to original generated bonds", highlight=False)
                self.console.print("6. Save to file and continue", highlight=False)
                
                # Create options map for navigation
                options_map = {
                    "1": "review",
                    "2": "add",
                    "3": "remove",
                    "4": "view_by_category",
                    "5": "reset",
                    "6": "save"
                }

                # Add standard navigation options
                self.add_navigation_options(options_map)

                # Create choices list
                choices = ["1", "2", "3", "4", "5", "6", "b", "s", "m", "x"]

                display_options_map = {
                    "1": "Review and edit bond commands",
                    "2": "Add new bond command",
                    "3": "Remove bond commands",
                    "4": "View bonds by category",
                    "5": "Reset to original generated bonds",
                    "6": "Save to file and continue",
                    "b": "Back",
                    "s": "Section Menu",
                    "m": "Main Menu",
                    "x": "Exit",
                }

                choice = prompt_with_context(
                    self.processor,
                    "Enter your choice",
                    choices=choices,
                    module=f"{self.bond_manager.module_name} Bond Management",
                    description="Bond management menu choice",
                    options_map=display_options_map,
                )
                
                # Handle navigation choices first
                nav_action = self.handle_navigation_choice(choice, options_map)
                if nav_action:
                    if nav_action == "back":
                        self.exit_menu()
                        return "back"
                    else:
                        # Execute navigation (section/main/exit)
                        self.execute_navigation(nav_action)
                        self.exit_menu()
                        return nav_action
                
                # Handle menu options
                if choice == "1":
                    self.bond_manager._review_and_edit_bonds()
                elif choice == "2":
                    self.bond_manager._add_new_bond()
                elif choice == "3":
                    self.bond_manager._remove_bonds()
                elif choice == "4":
                    self.bond_manager._view_bonds_by_category()
                elif choice == "5":
                    self.bond_manager._reset_bonds()
                elif choice == "6":
                    self.exit_menu()
                    return True  # Save and continue
                    
        except Exception as e:
            self.console.print(f"[red]Error in bond management: {str(e)}[/red]")
            self.exit_menu()
            return False
        finally:
            self.exit_menu()