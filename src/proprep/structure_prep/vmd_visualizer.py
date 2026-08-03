"""
Enhanced VMD Visualization Module - Complete Workspace-Aware Version

This is a complete replacement for vmd_visualizer.py that includes all existing components
plus the new workspace-aware functionality.
"""

import datetime
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
import urllib.request
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

# Rich UI library
try:
    from rich.console import Console
    from rich.panel import Panel
    from proprep.utils.prompts import prompt_with_context, confirm_with_context
    from rich.status import Status
    from rich.table import Table
    has_rich = True
except ImportError:
    has_rich = False
    print("Installing rich for better UI...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rich"])
    from rich.console import Console
    from rich.panel import Panel
    from proprep.utils.prompts import prompt_with_context, confirm_with_context
    from rich.status import Status
    from rich.table import Table

# Biopython for structure handling
try:
    from Bio.PDB import PDBIO, PDBParser, Selection
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Structure import Structure
    has_biopython = True
except ImportError:
    has_biopython = False
    print("Installing biopython for structure analysis...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
    from Bio.PDB import PDBIO, PDBParser, Selection
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Model import Model
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Structure import Structure

# Import module registry functionality
from proprep.structure_prep.chem_comp_dict_fetcher import CCDParser
from proprep.utils.module_registry import ProcessingModule, register_module, registry


# ============================================================================
# EXISTING VMD COMPONENTS (Keep all existing functionality)
# ============================================================================

class EnhancedSelectionBuilder:
    """Enhanced Builder for VMD atom selections with intuitive menu interface"""

    def __init__(self, component_classifier, console, processor=None):
        """Initialize the selection builder"""
        self.component_classifier = component_classifier
        self.console = console
        self.processor = processor
        self.current_selection = "all"
        self.selection_history = []
        self.selection_text = "All atoms in the structure"

    def build_selection_interactive(self, structure=None):
        """Interactive selection builder with comprehensive options"""
        from rich.panel import Panel
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        # Start with all atoms selected
        self.current_selection = "all"
        self.selection_text = "All atoms in the structure"
        self.selection_history = [("all", "All atoms in the structure")]

        while True:
            # Show current selection
            self._display_current_selection()

            # Main selection categories
            self.console.print("\n[bold]Selection Categories:[/bold]")
            self.console.print("1. Basic atom/residue selection", highlight=False)
            self.console.print("2. Protein structure elements", highlight=False)
            self.console.print("3. Chemical properties", highlight=False)
            self.console.print("4. Spatial relationships", highlight=False)
            self.console.print("5. Combine with current selection", highlight=False)
            self.console.print("6. Advanced/Custom selection", highlight=False)
            self.console.print("7. Use current selection", highlight=False)
            self.console.print("8. Reset to 'all'", highlight=False)

            choice = prompt_with_context(
                self.processor,
                "Choose category",
                choices=["1", "2", "3", "4", "5", "6", "7", "8"],
                default="7",
                module="Structure Viewer - Selection Builder",
                description="Selection category",
                options_map={
                    "1": "Basic atom/residue selection",
                    "2": "Protein structure elements",
                    "3": "Chemical properties",
                    "4": "Spatial relationships",
                    "5": "Combine with current selection",
                    "6": "Advanced/Custom selection",
                    "7": "Use current selection",
                    "8": "Reset to 'all'",
                },
            )

            if choice == "1":
                self._basic_selection_menu()
            elif choice == "2":
                self._protein_structure_menu()
            elif choice == "3":
                self._chemical_properties_menu()
            elif choice == "4":
                self._spatial_relationships_menu()
            elif choice == "5":
                self._combine_selections_menu()
            elif choice == "6":
                self._custom_selection_menu()
            elif choice == "7":
                break
            elif choice == "8":
                self._update_selection("all", "All atoms in the structure")

        return self.current_selection

    def _display_current_selection(self):
        """Display current selection status"""
        panel_content = f"[bold]Current Selection:[/bold]\n{self.selection_text}\n\n[bold]VMD Syntax:[/bold]\n{self.current_selection}"
        self.console.print(Panel(panel_content, title="Selection Status"))

    def _basic_selection_menu(self):
        """Basic atom/residue selection menu"""
        self.console.print("\n[bold]Basic Selection Options:[/bold]")
        self.console.print("1. All atoms", highlight=False)
        self.console.print("2. Protein only", highlight=False)
        self.console.print("3. Water molecules", highlight=False)
        self.console.print("4. Specific chain(s)", highlight=False)
        self.console.print("5. Specific residue(s)", highlight=False)
        self.console.print("6. Specific atom name(s)", highlight=False)
        self.console.print("7. Back to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4", "5", "6", "7"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Basic selection option",
            options_map={
                "1": "All atoms",
                "2": "Protein only",
                "3": "Water molecules",
                "4": "Specific chain(s)",
                "5": "Specific residue(s)",
                "6": "Specific atom name(s)",
                "7": "Back to main menu",
            },
        )

        if choice == "1":
            self._update_selection("all", "All atoms in the structure")
        elif choice == "2":
            self._update_selection("protein", "Protein atoms only")
        elif choice == "3":
            self._update_selection("water", "Water molecules")
        elif choice == "4":
            self._chain_selection_menu()
        elif choice == "5":
            self._residue_selection_menu()
        elif choice == "6":
            self._atom_name_selection_menu()

    def _protein_structure_menu(self):
        """Protein structure elements menu"""
        self.console.print("\n[bold]Protein Structure Elements:[/bold]")
        self.console.print("1. Backbone atoms (CA, N, C, O)", highlight=False)
        self.console.print("2. Side chains", highlight=False)
        self.console.print("3. Alpha carbons only", highlight=False)
        self.console.print("4. Secondary structure elements", highlight=False)
        self.console.print("5. N-terminus", highlight=False)
        self.console.print("6. C-terminus", highlight=False)
        self.console.print("7. Back to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4", "5", "6", "7"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Protein structure element",
            options_map={
                "1": "Backbone atoms (CA, N, C, O)",
                "2": "Side chains",
                "3": "Alpha carbons only",
                "4": "Secondary structure elements",
                "5": "N-terminus",
                "6": "C-terminus",
                "7": "Back to main menu",
            },
        )

        if choice == "1":
            self._update_selection("backbone", "Protein backbone atoms")
        elif choice == "2":
            self._update_selection("sidechain", "Protein side chain atoms")
        elif choice == "3":
            self._update_selection("name CA", "Alpha carbon atoms")
        elif choice == "4":
            self._secondary_structure_menu()
        elif choice == "5":
            self._update_selection("protein and resid 1", "N-terminus residue")
        elif choice == "6":
            # This would need to be dynamic based on actual structure
            self._update_selection("protein and resid last", "C-terminus residue")

    def _chemical_properties_menu(self):
        """Chemical properties selection menu"""
        self.console.print("\n[bold]Chemical Properties:[/bold]")
        self.console.print("1. Charged residues", highlight=False)
        self.console.print("2. Hydrophobic residues", highlight=False)
        self.console.print("3. Polar residues", highlight=False)
        self.console.print("4. Aromatic residues", highlight=False)
        self.console.print("5. Metal ions", highlight=False)
        self.console.print("6. Hetero compounds", highlight=False)
        self.console.print("7. Back to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4", "5", "6", "7"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Chemical property",
            options_map={
                "1": "Charged residues",
                "2": "Hydrophobic residues",
                "3": "Polar residues",
                "4": "Aromatic residues",
                "5": "Metal ions",
                "6": "Hetero compounds",
                "7": "Back to main menu",
            },
        )

        if choice == "1":
            charged_residues = "resname ARG LYS ASP GLU HIS"
            self._update_selection(charged_residues, "Charged residues")
        elif choice == "2":
            hydrophobic_residues = "resname ALA VAL LEU ILE MET PHE TRP PRO"
            self._update_selection(hydrophobic_residues, "Hydrophobic residues")
        elif choice == "3":
            polar_residues = "resname SER THR ASN GLN TYR CYS"
            self._update_selection(polar_residues, "Polar residues")
        elif choice == "4":
            aromatic_residues = "resname PHE TYR TRP HIS"
            self._update_selection(aromatic_residues, "Aromatic residues")
        elif choice == "5":
            metal_selection = "name ZN MG CA FE MN CU NI CO"
            self._update_selection(metal_selection, "Metal ions")
        elif choice == "6":
            self._update_selection("not protein and not water", "Hetero compounds")

    def _spatial_relationships_menu(self):
        """Spatial relationships menu"""
        self.console.print("\n[bold]Spatial Relationships:[/bold]")
        self.console.print("1. Within distance of selection", highlight=False)
        self.console.print("2. Same residue as selection", highlight=False)
        self.console.print("3. Same chain as selection", highlight=False)
        self.console.print("4. Contact with selection", highlight=False)
        self.console.print("5. Back to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4", "5"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Spatial relationship",
            options_map={
                "1": "Within distance of selection",
                "2": "Same residue as selection",
                "3": "Same chain as selection",
                "4": "Contact with selection",
                "5": "Back to main menu",
            },
        )

        if choice == "1":
            self._within_distance_menu()
        elif choice == "2":
            # Would need to be based on current selection
            self.console.print("[yellow]This requires a current selection with residues[/yellow]")
        elif choice == "3":
            self.console.print("[yellow]This requires a current selection with chains[/yellow]")
        elif choice == "4":
            self._contact_selection_menu()

    def _chain_selection_menu(self):
        """Chain selection submenu"""
        chain_input = prompt_with_context(
            self.processor,
            "Enter chain ID(s) (e.g., A or A B C)",
            default="A",
            module="Structure Viewer - Selection Builder",
            description="Chain IDs to select",
        )
        chains = chain_input.strip().split()

        if len(chains) == 1:
            selection = f"chain {chains[0]}"
            description = f"Chain {chains[0]}"
        else:
            chain_list = " ".join(chains)
            selection = f"chain {chain_list}"
            description = f"Chains {', '.join(chains)}"

        self._update_selection(selection, description)

    def _residue_selection_menu(self):
        """Residue selection submenu"""
        self.console.print("\n[bold]Residue Selection:[/bold]")
        self.console.print("1. Specific residue numbers", highlight=False)
        self.console.print("2. Residue range", highlight=False)
        self.console.print("3. Specific residue types", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Residue selection method",
            options_map={
                "1": "Specific residue numbers",
                "2": "Residue range",
                "3": "Specific residue types",
            },
        )

        if choice == "1":
            residue_input = prompt_with_context(
                self.processor,
                "Enter residue number(s) (e.g., 123 or 123 456 789)",
                default="1",
                module="Structure Viewer - Selection Builder",
                description="Residue numbers",
            )
            residues = residue_input.strip().split()
            if len(residues) == 1:
                selection = f"resid {residues[0]}"
                description = f"Residue {residues[0]}"
            else:
                residue_list = " ".join(residues)
                selection = f"resid {residue_list}"
                description = f"Residues {', '.join(residues)}"

        elif choice == "2":
            start = prompt_with_context(
                self.processor,
                "Start residue number",
                default="1",
                module="Structure Viewer - Selection Builder",
                description="Residue range start",
            )
            end = prompt_with_context(
                self.processor,
                "End residue number",
                default="10",
                module="Structure Viewer - Selection Builder",
                description="Residue range end",
            )
            selection = f"resid {start} to {end}"
            description = f"Residues {start} to {end}"

        elif choice == "3":
            restype_input = prompt_with_context(
                self.processor,
                "Enter residue type(s) (e.g., ALA or ALA GLY PRO)",
                default="ALA",
                module="Structure Viewer - Selection Builder",
                description="Residue types (three-letter codes)",
            )
            restypes = restype_input.strip().split()
            if len(restypes) == 1:
                selection = f"resname {restypes[0]}"
                description = f"Residue type {restypes[0]}"
            else:
                restype_list = " ".join(restypes)
                selection = f"resname {restype_list}"
                description = f"Residue types {', '.join(restypes)}"

        self._update_selection(selection, description)

    def _atom_name_selection_menu(self):
        """Atom name selection submenu"""
        atom_input = prompt_with_context(
            self.processor,
            "Enter atom name(s) (e.g., CA or CA CB CG)",
            default="CA",
            module="Structure Viewer - Selection Builder",
            description="Atom name(s) to select",
        )
        atoms = atom_input.strip().split()

        if len(atoms) == 1:
            selection = f"name {atoms[0]}"
            description = f"Atom {atoms[0]}"
        else:
            atom_list = " ".join(atoms)
            selection = f"name {atom_list}"
            description = f"Atoms {', '.join(atoms)}"

        self._update_selection(selection, description)

    def _secondary_structure_menu(self):
        """Secondary structure selection menu"""
        self.console.print("\n[bold]Secondary Structure:[/bold]")
        self.console.print("1. Alpha helices", highlight=False)
        self.console.print("2. Beta sheets", highlight=False)
        self.console.print("3. Loops/coils", highlight=False)
        self.console.print("4. Turns", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Secondary structure element",
            options_map={
                "1": "Alpha helices",
                "2": "Beta sheets",
                "3": "Loops/coils",
                "4": "Turns",
            },
        )

        # Note: These selections require VMD to compute secondary structure
        if choice == "1":
            self._update_selection("structure H", "Alpha helices")
        elif choice == "2":
            self._update_selection("structure E", "Beta sheets")
        elif choice == "3":
            self._update_selection("structure C", "Loops and coils")
        elif choice == "4":
            self._update_selection("structure T", "Turns")

    def _within_distance_menu(self):
        """Within distance selection menu"""
        distance = prompt_with_context(
            self.processor,
            "Enter distance in Angstroms",
            default="5.0",
            module="Structure Viewer - Selection Builder",
            description="Within-distance cutoff (Å)",
        )
        target = prompt_with_context(
            self.processor,
            "Enter target selection (e.g., 'resname ZN' for zinc)",
            default="resname ZN",
            module="Structure Viewer - Selection Builder",
            description="Target selection for within-distance",
        )

        selection = f"within {distance} of ({target})"
        description = f"Within {distance} Å of {target}"
        self._update_selection(selection, description)

    def _contact_selection_menu(self):
        """Contact selection menu"""
        distance = prompt_with_context(
            self.processor,
            "Enter contact distance in Angstroms",
            default="3.5",
            module="Structure Viewer - Selection Builder",
            description="Contact distance (Å)",
        )
        target = prompt_with_context(
            self.processor,
            "Enter target selection for contacts",
            default="protein",
            module="Structure Viewer - Selection Builder",
            description="Target selection for contacts",
        )

        selection = f"within {distance} of ({target}) and not ({target})"
        description = f"In contact with {target} (≤{distance} Å)"
        self._update_selection(selection, description)

    def _custom_selection_menu(self):
        """Custom selection input menu"""
        self.console.print("\n[bold]Custom Selection:[/bold]")
        self.console.print("Enter a custom VMD selection string")
        self.console.print("Examples:")
        self.console.print("  protein and resid 1 to 50")
        self.console.print("  name CA and chain A")
        self.console.print("  resname ALA and sidechain")

        custom_selection = prompt_with_context(
            self.processor,
            "Enter VMD selection",
            default=self.current_selection,
            module="Structure Viewer - Selection Builder",
            description="Custom VMD selection string",
        )
        description = prompt_with_context(
            self.processor,
            "Enter description",
            default="Custom selection",
            module="Structure Viewer - Selection Builder",
            description="Description for custom selection",
        )

        self._update_selection(custom_selection, description)

    def _combine_selections_menu(self):
        """Combine current selection with another"""
        self.console.print("\n[bold]Combine Selections:[/bold]")
        self.console.print("1. AND with another selection", highlight=False)
        self.console.print("2. OR with another selection", highlight=False)
        self.console.print("3. NOT (exclude) another selection", highlight=False)
        self.console.print("4. Back to main menu", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose operation",
            choices=["1", "2", "3", "4"],
            default="1",
            module="Structure Viewer - Selection Builder",
            description="Selection combination operation",
            options_map={
                "1": "AND with another selection",
                "2": "OR with another selection",
                "3": "NOT (exclude) another selection",
                "4": "Back to main menu",
            },
        )

        if choice == "4":
            return

        operations = {"1": ("and", "AND"), "2": ("or", "OR"), "3": ("and not", "NOT")}

        # First save current selection
        current_sel = self.current_selection
        current_desc = self.selection_text

        # Now build a new selection to combine with
        self.console.print(
            f"\n[bold]Build selection to {operations[choice][1]} with current selection:[/bold]"
        )

        # Reset selection for building second part
        tmp_sel = self.current_selection
        tmp_text = self.selection_text
        tmp_history = self.selection_history.copy()

        self.current_selection = "all"
        self.selection_text = "All atoms in the structure"
        self.selection_history = [("all", "All atoms in the structure")]

        # Build the second selection
        self.build_selection_interactive()
        second_sel = self.current_selection
        second_desc = self.selection_text

        # Restore and create the combined selection
        self.current_selection = tmp_sel
        self.selection_text = tmp_text
        self.selection_history = tmp_history

        combined_sel = f"({current_sel}) {operations[choice][0]} ({second_sel})"
        combined_desc = f"{current_desc} {operations[choice][1]} {second_desc}"

        self._update_selection(combined_sel, combined_desc)

    def _update_selection(self, selection, description, operation="and"):
        """Update the current selection"""
        # For first selection, just set it
        if self.current_selection == "all" and selection != "all":
            self.current_selection = selection
            self.selection_text = description
        else:
            # Combine with existing selection
            self.current_selection = (
                f"({self.current_selection}) {operation} ({selection})"
            )
            self.selection_text = f"{self.selection_text} {operation} {description}"

        # Add to history
        self.selection_history.append((selection, description))


class RepresentationManager:
    """Manages VMD representations with comprehensive drawing methods and parameters"""

    # Drawing Method mappings with default parameter values
    DRAWING_METHODS = {
        "Lines": {
            "description": "Simple line representation",
            "parameters": {"thickness": 1.0},
        },
        "Bonds": {
            "description": "Cylinder bonds with spherical atoms",
            "parameters": {"radius": 0.3, "resolution": 10},
        },
        "DynamicBonds": {
            "description": "Bonds computed dynamically",
            "parameters": {"radius": 0.3, "resolution": 10, "distance": 1.6},
        },
        "HBonds": {
            "description": "Hydrogen bonds",
            "parameters": {
                "angle_cutoff": 20.0,
                "distance_cutoff": 3.0,
                "angle_mode": "max",
            },
        },
        "Points": {
            "description": "Point representation",
            "parameters": {"size": 1.0},
        },
        "VDW": {
            "description": "Van der Waals spheres",
            "parameters": {"sphere_scale": 1.0, "sphere_resolution": 10},
        },
        "CPK": {
            "description": "Space-filling with bonds",
            "parameters": {
                "sphere_scale": 1.0,
                "cylinder_radius": 0.3,
                "sphere_resolution": 10,
            },
        },
        "Licorice": {
            "description": "Stick representation",
            "parameters": {"bond_radius": 0.3, "sphere_resolution": 10},
        },
        "Polyhedra": {
            "description": "Coordination polyhedra",
            "parameters": {"distance_cutoff": 3.0},
        },
        "Trace": {
            "description": "CA trace",
            "parameters": {"bond_radius": 0.3},
        },
        "Tube": {
            "description": "Smooth tube through CA atoms",
            "parameters": {"radius": 0.3, "resolution": 10},
        },
        "Ribbons": {
            "description": "Ribbon through backbone",
            "parameters": {"width": 3.0, "thickness": 0.3},
        },
        "Cartoon": {
            "description": "Cartoon representation",
            "parameters": {
                "thickness": 0.3,
                "resolution": 10,
                "spline": "Catmull-Rom",
                "aspect_ratio": 4.1,
            },
        },
        "NewCartoon": {
            "description": "Enhanced cartoon representation",
            "parameters": {
                "thickness": 0.3,
                "resolution": 10,
                "spline": "Catmull-Rom",
                "aspect_ratio": 4.1,
            },
        },
        "PaperChain": {
            "description": "Paper chain representation",
            "parameters": {"thickness": 0.3, "resolution": 10},
        },
        "Twister": {
            "description": "Twisted ribbon representation",
            "parameters": {"thickness": 0.3, "resolution": 10},
        },
        "QuickSurf": {
            "description": "Quick molecular surface",
            "parameters": {
                "resolution": 1.0,
                "radius_scale": 1.0,
                "density_isovalue": 0.5,
                "grid_spacing": 1.0,
            },
        },
        "MSMS": {
            "description": "MSMS molecular surface",
            "parameters": {"sample_density": 1.5, "probe_radius": 1.5},
        },
        "Surf": {
            "description": "Solvent accessible surface",
            "parameters": {"probe_radius": 1.4},
        },
        "VolumeSlice": {
            "description": "Volume slice representation",
            "parameters": {"slice_axis": 2, "slice_position": 0.5},
        },
        "Isosurface": {
            "description": "Isosurface representation",
            "parameters": {"isovalue": 0.5, "draw_box": False, "box_color": "black"},
        },
        "FieldLines": {
            "description": "Field lines representation",
            "parameters": {"seed_resolution": 1, "line_thickness": 1},
        },
        "Orbital": {
            "description": "Molecular orbital representation",
            "parameters": {"isovalue": 0.05, "grid_spacing": 0.1},
        },
    }

    # Coloring method mappings
    COLORING_METHODS = {
        "Name": "Name",
        "Type": "Type",
        "Element": "Element",
        "ResName": "ResName",
        "ResType": "ResType",
        "ResID": "ResID",
        "Chain": "Chain",
        "SegName": "SegName",
        "Conformation": "Conformation",
        "Molecule": "Molecule",
        "Structure": "Structure",
        "ColorID": "ColorID",
        "Beta": "Beta",
        "Occupancy": "Occupancy",
        "Mass": "Mass",
        "Charge": "Charge",
        "Pos": "Pos",
        "PosX": "PosX",
        "PosY": "PosY",
        "PosZ": "PosZ",
        "User": "User",
        "User2": "User2",
        "User3": "User3",
        "User4": "User4",
        "Fragment": "Fragment",
        "Index": "Index",
        "Backbone": "Backbone",
        "Throb": "Throb",
        "PhysicalTime": "PhysicalTime",
        "Timestep": "Timestep",
        "Velocity": "Velocity",
        "Volume": "Volume",
    }

    # Material type mappings
    MATERIAL_TYPES = {
        "Opaque": "Opaque",
        "Transparent": "Transparent",
        "BrushedMetal": "BrushedMetal",
        "Diffuse": "Diffuse",
        "Ghost": "Ghost",
        "Glass1": "Glass1",
        "Glass2": "Glass2",
        "Glass3": "Glass3",
        "Glossy": "Glossy",
        "HardPlastic": "HardPlastic",
        "MetallicPastel": "MetallicPastel",
        "Steel": "Steel",
        "Translucent": "Translucent",
        "Edgy": "Edgy",
        "EdgyShiny": "EdgyShiny",
        "EdgyGlass": "EdgyGlass",
        "Goodsell": "Goodsell",
        "AOShiny": "AOShiny",
        "AOChalky": "AOChalky",
        "AOEdgy": "AOEdgy",
        "BlownGlass": "BlownGlass",
        "GlassBubble": "GlassBubble",
        "RTChrome": "RTChrome",
    }

    # Color ID mappings
    COLOR_IDS = {
        "Blue": 0,
        "Red": 1,
        "Gray": 2,
        "Orange": 3,
        "Yellow": 4,
        "Tan": 5,
        "Silver": 6,
        "Green": 7,
        "White": 8,
        "Pink": 9,
        "Cyan": 10,
        "Purple": 11,
        "Lime": 12,
        "Mauve": 13,
        "Ochre": 14,
        "Iceblue": 15,
        "Black": 16,
        "Yellow2": 17,
        "Yellow3": 18,
        "Green2": 19,
        "Green3": 20,
        "Cyan2": 21,
        "Cyan3": 22,
        "Blue2": 23,
        "Blue3": 24,
        "Violet": 25,
        "Violet2": 26,
        "Magenta": 27,
        "Magenta2": 28,
        "Red2": 29,
        "Red3": 30,
        "Orange2": 31,
        "Orange3": 32,
    }

    def __init__(self):
        """Initialize representation manager"""
        self.representations = []

    def create_representation(self, rep_params):
        """Create a representation from parameters"""
        # Verify required parameters
        required_params = ["drawing_method", "coloring_method", "material"]
        for param in required_params:
            if param not in rep_params:
                raise ValueError(f"Missing required parameter: {param}")

        # Translate to VMD-specific values
        drawing_method = rep_params["drawing_method"]
        vmd_drawing = self.DRAWING_METHODS.get(drawing_method, {}).get("parameters", {})
        vmd_coloring = self.COLORING_METHODS.get(
            rep_params["coloring_method"], rep_params["coloring_method"]
        )
        vmd_material = self.MATERIAL_TYPES.get(
            rep_params["material"], rep_params["material"]
        )

        # Create representation
        rep = {
            "drawing_method": drawing_method,
            "drawing_parameters": vmd_drawing.copy(),  # Copy default parameters
            "coloring_method": vmd_coloring,
            "material": vmd_material,
            "selection": rep_params.get("selection", "all"),
            "selection_name": rep_params.get("selection_name", "Selection"),
        }

        # Override default drawing parameters with any provided
        if "drawing_parameters" in rep_params:
            for param, value in rep_params["drawing_parameters"].items():
                rep["drawing_parameters"][param] = value

        # Add optional parameters
        if "opacity" in rep_params:
            rep["opacity"] = rep_params["opacity"]

        if "color_value" in rep_params:
            rep["color_value"] = rep_params["color_value"]

        # Add representation
        self.representations.append(rep)

        return rep

    def generate_vmd_commands(self, rep, index):
        """Generate VMD Tcl commands for a representation"""
        commands = []
        drawing_method = rep["drawing_method"]
        drawing_params = rep.get("drawing_parameters", {})

        # Generate the appropriate representation command based on drawing method
        if drawing_method == "Lines":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("thickness", 1.0)}'
            )

        elif drawing_method in ["Bonds", "DynamicBonds"]:
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("radius", 0.3)} {drawing_params.get("resolution", 10)}'
            )

        elif drawing_method == "HBonds":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("angle_cutoff", 20.0)} {drawing_params.get("distance_cutoff", 3.0)}'
            )

        elif drawing_method == "Points":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("size", 1.0)}'
            )

        elif drawing_method == "VDW":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("sphere_scale", 1.0)} {drawing_params.get("sphere_resolution", 10)}'
            )

        elif drawing_method == "CPK":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("sphere_scale", 1.0)} {drawing_params.get("cylinder_radius", 0.3)} {drawing_params.get("sphere_resolution", 10)}'
            )

        elif drawing_method == "Licorice":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("bond_radius", 0.3)} {drawing_params.get("sphere_resolution", 10)}'
            )

        elif drawing_method == "Polyhedra":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("distance_cutoff", 3.0)}'
            )

        elif drawing_method == "Trace":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("bond_radius", 0.3)}'
            )

        elif drawing_method == "Tube":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("radius", 0.3)} {drawing_params.get("resolution", 10)}'
            )

        elif drawing_method == "Ribbons":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("width", 3.0)} {drawing_params.get("thickness", 0.3)}'
            )

        elif drawing_method in ["Cartoon", "NewCartoon", "PaperChain", "Twister"]:
            # Correct Cartoon syntax: just the drawing method name
            commands.append(f'mol representation {drawing_method}')

        elif drawing_method == "MSMS":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("sample_density", 1.5)} {drawing_params.get("probe_radius", 1.5)}'
            )

        elif drawing_method == "QuickSurf":
            commands.append(
                f'mol representation {drawing_method} {drawing_params.get("resolution", 1.0)} {drawing_params.get("radius_scale", 1.0)} {drawing_params.get("density_isovalue", 0.5)} {drawing_params.get("grid_spacing", 1.0)}'
            )

        else:
            # Default for other types
            commands.append(f"mol representation {drawing_method}")

        # Set coloring method
        coloring_method = rep["coloring_method"]
        if coloring_method == "ColorID" and "color_value" in rep:
            # For ColorID, we need to use the numeric color value
            color_value = rep["color_value"]
            if isinstance(color_value, str):
                # Convert color name to number
                color_value = self.COLOR_IDS.get(color_value, 1)
            commands.append(f'mol color ColorID {color_value}')
        else:
            commands.append(f'mol color {coloring_method}')

        # Set selection
        commands.append(f'mol selection {{{rep["selection"]}}}')

        # Set material
        commands.append(f'mol material {rep["material"]}')

        # Add the representation
        commands.append(f"mol addrep top")

        # Handle opacity after adding the representation
        if "opacity" in rep and rep["opacity"] < 1.0:
            opacity = rep["opacity"]
            material = rep["material"]
            # Correct VMD syntax for changing material opacity
            commands.append(f'material change opacity {material} {opacity}')

        return "\n".join(commands)

    def get_drawing_method_options(self):
        """Return UI-friendly drawing method options"""
        return list(self.DRAWING_METHODS.keys())

    def get_drawing_method_parameters(self, method):
        """Return parameters for a specific drawing method"""
        return self.DRAWING_METHODS.get(method, {}).get("parameters", {}).copy()

    def get_coloring_method_options(self):
        """Return UI-friendly coloring method options"""
        return list(self.COLORING_METHODS.keys())

    def get_material_options(self):
        """Return UI-friendly material options"""
        return list(self.MATERIAL_TYPES.keys())

    def get_color_id_options(self):
        """Return UI-friendly color ID options"""
        return list(self.COLOR_IDS.keys())

    def get_color_id_value(self, color_name):
        """Return the numeric value for a color ID name"""
        return self.COLOR_IDS.get(color_name, 0)


class VMDStateGenerator:
    """Generates VMD state files from selections and representations"""

    def __init__(self):
        """Initialize state generator"""
        self.rep_manager = RepresentationManager()

    def generate_state_file(
        self, filepath, structure, representations, view_settings=None
    ):
        """Generate a complete VMD state file"""
        # Create base PDB file
        pdb_path = f"{filepath}.pdb"
        self._write_pdb(pdb_path, structure)

        # Begin state file
        state = self._generate_header()

        # Add molecule loading command
        state += f"mol new {os.path.basename(pdb_path)} type pdb first 0 last -1 step 1 waitfor 1\n"

        # Remove the default representation that VMD automatically creates
        state += "mol delrep 0 top\n\n"

        # Add representations
        for i, rep in enumerate(representations):
            # Make sure selection exists
            if "selection" in rep and rep["selection"]:
                state += f'# Representation {i+1}: {rep.get("selection_name", "Selection")}\n'
                state += self.rep_manager.generate_vmd_commands(rep, i) + "\n\n"

        # Add view settings
        state += self._generate_view_settings(view_settings)

        # Write state file
        state_file = f"{filepath}.vmd"
        with open(state_file, "w") as f:
            f.write(state)

        # Generate a Tcl script to load everything (optional)
        self._generate_loader_script(filepath)

        return pdb_path, state_file

    def _write_pdb(self, filepath, structure):
        """Write structure to a PDB file"""
        io = PDBIO()
        io.set_structure(structure)
        io.save(filepath)

    def _generate_header(self):
        """Generate VMD state file header"""
        header = f"""# VMD visualization state file
# Generated by MPSA VMD Visualization Module
# Date: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
#
# To use this file:
# 1. Start VMD
# 2. From the VMD Main menu, choose "File" -> "Load Visualization State"
# 3. Select this file
#
# Or run the companion loader script if it was generated

# Turn off autocentering
axes location off

# Load structure and remove default representation
"""
        return header

    def _generate_view_settings(self, view_settings=None):
        """Generate view setting commands"""
        settings = """
# Configure default view
display projection Orthographic
display depthcue off
display cuemode Linear
display cuestart 0.5
display cueend 10.0
display cuedensity 0.4
display cuemode Exp2
"""

        # Add background color
        if view_settings and "background" in view_settings:
            bg = view_settings["background"].lower()
            if bg == "white":
                settings += "color Display Background white\n"
            elif bg == "gray" or bg == "grey":
                settings += "color Display Background gray\n"
            else:
                # Default to black
                settings += "color Display Background black\n"
        else:
            settings += "color Display Background black\n"

        # Add axes display
        if view_settings and "axes" in view_settings and view_settings["axes"]:
            settings += "axes location LowerRight\n"
        else:
            settings += "axes location Off\n"

        # Add other view settings as needed

        return settings

    def _generate_loader_script(self, filepath):
        """Generate a Tcl script to load the visualization"""
        script = f"""#!/usr/bin/tclsh
# VMD loader script for {filepath}.vmd
# Generated by MPSA VMD Visualization Module

# Check if VMD is in the path
set vmd_cmd "vmd"
if {{[catch {{exec which $vmd_cmd}} result]}} {{
    puts "Error: VMD not found in path"
    puts "Please make sure VMD is installed and in your PATH"
    exit 1
}}

# Load visualization state
exec $vmd_cmd -e {os.path.basename(filepath)}.vmd
"""

        # Write script
        script_file = f"{filepath}_load.tcl"
        with open(script_file, "w") as f:
            f.write(script)

        # Make executable
        try:
            os.chmod(script_file, 0o755)
        except:
            pass

        return script_file


class RepresentationsManager:
    """
    Interactive manager for creating and manipulating multiple representations
    for VMD visualization
    """

    def __init__(self, console, component_classifier, structure=None, processor=None):
        """Initialize the representations manager"""
        self.console = console
        self.component_classifier = component_classifier
        self.structure = structure
        self.processor = processor
        self.representations = []
        self.rep_manager = RepresentationManager()

    def run_representations_manager(self):
        """Run the interactive representations manager"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        while True:
            self._display_representations_summary()

            # Options
            self.console.print("\n[bold]Options:[/bold]")
            self.console.print("1. Add new representation", highlight=False)
            self.console.print("2. Edit existing representation", highlight=False)
            self.console.print("3. Delete representation", highlight=False)
            self.console.print("4. Preview representations", highlight=False)
            self.console.print("5. Finish and generate visualization", highlight=False)
            self.console.print("6. Cancel", highlight=False)

            choice = prompt_with_context(
                self.processor,
                "Choose option",
                choices=["1", "2", "3", "4", "5", "6"],
                default="1",
                module="Structure Viewer - Representations",
                description="Representations manager option",
                options_map={
                    "1": "Add new representation",
                    "2": "Edit existing representation",
                    "3": "Delete representation",
                    "4": "Preview representations",
                    "5": "Finish and generate visualization",
                    "6": "Cancel",
                },
            )

            if choice == "1":
                self._add_representation()
            elif choice == "2":
                self._edit_representation()
            elif choice == "3":
                self._delete_representation()
            elif choice == "4":
                self._preview_representations()
            elif choice == "5":
                return self.representations
            elif choice == "6":
                return []

    def _display_representations_summary(self):
        """Display summary of current representations"""
        if not self.representations:
            self.console.print("\n[yellow]No representations defined yet[/yellow]")
            return

        table = Table(title="Current Representations")
        table.add_column("#", width=3)
        table.add_column("Name", width=20)
        table.add_column("Drawing", width=15)
        table.add_column("Coloring", width=15)
        table.add_column("Material", width=15)
        table.add_column("Selection", width=25)

        for i, rep in enumerate(self.representations, 1):
            # Handle material with opacity
            material = rep.get("material", "Opaque")
            if "opacity" in rep:
                material = f"{material} ({rep['opacity']:.1f})"

            # Truncate selection for display
            selection = rep.get("selection", "all")
            if len(selection) > 30:
                selection = selection[:27] + "..."

            table.add_row(
                str(i),
                rep.get("selection_name", f"Rep {i}"),
                rep.get("drawing_method", "Unknown"),
                rep.get("coloring_method", "Default"),
                material,
                selection,
            )

        self.console.print(table)

    def _add_representation(self):
        """Add a new representation"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        self.console.print("\n[bold]===== Create New Representation =====[/bold]")

        # Step 1: Define what to show (atom selection)
        self.console.print("\n[bold]Step 1: Define what to show (Selection)[/bold]")

        selection_builder = EnhancedSelectionBuilder(
            self.component_classifier, self.console, processor=self.processor
        )
        selection = selection_builder.build_selection_interactive(self.structure)
        selection_text = selection_builder.selection_text

        # Step 2: Define visual representation (drawing method, coloring)
        self.console.print("\n[bold]Step 2: Define how to show it[/bold]")

        # Drawing method
        drawing_methods = self.rep_manager.get_drawing_method_options()
        self.console.print("\n[bold]Available Drawing Methods:[/bold]")

        # Group methods by category for better organization
        protein_methods = [
            "Cartoon",
            "NewCartoon",
            "Ribbons",
            "NewRibbons",
            "Tube",
            "Trace",
        ]
        surface_methods = ["QuickSurf", "MSMS", "Surf"]
        atomic_methods = ["Lines", "Bonds", "Licorice", "VDW", "CPK"]

        # Create a mapping from display index to actual method index
        current_index = 1
        method_indices = {}

        self.console.print("\n[bold cyan]Protein Backbone Methods:[/bold cyan]")
        for method in protein_methods:
            if method in drawing_methods:
                self.console.print(f"{current_index}. {method}")
                method_indices[current_index] = drawing_methods.index(method)
                current_index += 1

        self.console.print("\n[bold cyan]Surface Methods:[/bold cyan]")
        for method in surface_methods:
            if method in drawing_methods:
                self.console.print(f"{current_index}. {method}")
                method_indices[current_index] = drawing_methods.index(method)
                current_index += 1

        self.console.print("\n[bold cyan]Atomic Detail Methods:[/bold cyan]")
        for method in atomic_methods:
            if method in drawing_methods:
                self.console.print(f"{current_index}. {method}")
                method_indices[current_index] = drawing_methods.index(method)
                current_index += 1

        # Get user choice and map back to actual method
        drawing_options_map = {
            str(k): drawing_methods[v] for k, v in method_indices.items()
        }
        method_choice = prompt_with_context(
            self.processor,
            "Select drawing method",
            choices=[str(i) for i in method_indices.keys()],
            default="1",
            module="Structure Viewer - Representations",
            description="VMD drawing method",
            options_map=drawing_options_map,
        )

        method_idx = method_indices[int(method_choice)]
        drawing_method = drawing_methods[method_idx]

        # Configure drawing parameters
        drawing_params = self.rep_manager.get_drawing_method_parameters(drawing_method)

        if confirm_with_context(
            self.processor,
            "Configure drawing method parameters?",
            default=False,
            module="Structure Viewer - Representations",
            description="Configure drawing method parameters",
        ):
            drawing_params = self._configure_drawing_parameters(
                drawing_method, drawing_params
            )

        # Coloring method
        coloring_methods = self.rep_manager.get_coloring_method_options()
        self.console.print("\n[bold]Coloring Method:[/bold]")

        # Recommend coloring method based on selection and drawing
        recommended_coloring = self._recommend_coloring_method(
            selection, drawing_method
        )
        recommended_idx = (
            coloring_methods.index(recommended_coloring)
            if recommended_coloring in coloring_methods
            else 0
        )

        self.console.print(
            f"\n[yellow]Recommended for this selection: {recommended_coloring}[/yellow]"
        )

        for i, method in enumerate(coloring_methods):
            if method == recommended_coloring:
                self.console.print(f"{i+1}. {method} [green](Recommended)[/green]")
            else:
                self.console.print(f"{i+1}. {method}")

        coloring_options_map = {
            str(i + 1): method for i, method in enumerate(coloring_methods)
        }
        color_idx = prompt_with_context(
            self.processor,
            "Select coloring method",
            choices=[str(i + 1) for i in range(len(coloring_methods))],
            default=str(recommended_idx + 1),
            module="Structure Viewer - Representations",
            description="VMD coloring method",
            options_map=coloring_options_map,
        )

        coloring_method = coloring_methods[int(color_idx) - 1]

        # Handle ColorID selection if needed
        color_value = None
        if coloring_method == "ColorID":
            color_options = self.rep_manager.get_color_id_options()
            self.console.print("\n[bold]Available Colors:[/bold]")

            # Display color options in a grid
            num_columns = 3
            for i in range(0, len(color_options), num_columns):
                row = color_options[i : i + num_columns]
                row_display = "  ".join(
                    [f"{i+j+1}. {color}" for j, color in enumerate(row)]
                )
                self.console.print(row_display)

            color_choice_options_map = {
                str(i + 1): color for i, color in enumerate(color_options)
            }
            color_choice = prompt_with_context(
                self.processor,
                "Select color",
                choices=[str(i + 1) for i in range(len(color_options))],
                default="1",
                module="Structure Viewer - Representations",
                description="ColorID value",
                options_map=color_choice_options_map,
            )

            selected_color = color_options[int(color_choice) - 1]
            color_value = self.rep_manager.get_color_id_value(selected_color)

        # Step 3: Material properties
        self.console.print("\n[bold]Step 3: Material properties[/bold]")

        materials = self.rep_manager.get_material_options()
        self.console.print("\n[bold]Available Materials:[/bold]")

        # Recommend material based on drawing method
        recommended_material = self._recommend_material(drawing_method)
        recommended_idx = (
            materials.index(recommended_material)
            if recommended_material in materials
            else 0
        )

        for i, material in enumerate(materials):
            if material == recommended_material:
                self.console.print(f"{i+1}. {material} [green](Recommended)[/green]")
            else:
                self.console.print(f"{i+1}. {material}")

        material_options_map = {str(i + 1): m for i, m in enumerate(materials)}
        material_idx = prompt_with_context(
            self.processor,
            "Select material",
            choices=[str(i + 1) for i in range(len(materials))],
            default=str(recommended_idx + 1),
            module="Structure Viewer - Representations",
            description="Material for representation",
            options_map=material_options_map,
        )

        material = materials[int(material_idx) - 1]

        # Transparency (opacity)
        opacity = 1.0
        if material == "Transparent" or "Glass" in material:
            opacity = float(prompt_with_context(
                self.processor,
                "Enter opacity (0.0-1.0)",
                default="0.5",
                module="Structure Viewer - Representations",
                description="Opacity for transparent/glass material",
            ))

        # Step 4: Name the representation
        self.console.print("\n[bold]Step 4: Name the representation[/bold]")

        # Generate a default name based on selection and drawing method
        default_name = self._generate_representation_name(
            selection_text, drawing_method
        )

        name = prompt_with_context(
            self.processor,
            "Enter a name for this representation",
            default=default_name,
            module="Structure Viewer - Representations",
            description="Representation name",
        )

        # Create representation
        representation = {
            "selection_name": name,
            "selection": selection,
            "drawing_method": drawing_method,
            "drawing_parameters": drawing_params,
            "coloring_method": coloring_method,
            "material": material,
        }

        # Add opacity if not 1.0
        if opacity < 1.0:
            representation["opacity"] = opacity

        # Add color value for ColorID
        if color_value is not None:
            representation["color_value"] = color_value

        # Add to representations
        self.representations.append(representation)

        self.console.print(f"[green]Added representation: {name}[/green]")
        return True

    def _edit_representation(self):
        """Edit an existing representation"""
        from proprep.utils.prompts import prompt_with_context, confirm_with_context

        if not self.representations:
            self.console.print("[yellow]No representations to edit[/yellow]")
            return

        # Display representations with numbers
        self._display_representations_summary()

        # Get representation to edit
        rep_edit_options_map = {
            str(i + 1): rep.get("selection_name", f"Rep {i + 1}")
            for i, rep in enumerate(self.representations)
        }
        rep_idx = prompt_with_context(
            self.processor,
            "\nEnter the number of the representation to edit",
            choices=[str(i + 1) for i in range(len(self.representations))],
            default="1",
            module="Structure Viewer - Representations",
            description="Select representation to edit",
            options_map=rep_edit_options_map,
        )

        rep_idx = int(rep_idx) - 1
        representation = self.representations[rep_idx]

        self.console.print(
            f"\n[bold]Editing: {representation.get('selection_name', f'Rep {rep_idx+1}')}[/bold]"
        )

        # Show current representation details
        self.console.print("\n[bold]Current Settings:[/bold]")
        self.console.print(
            f"Drawing Method: {representation.get('drawing_method', 'Unknown')}"
        )

        # Display drawing parameters if any
        drawing_params = representation.get("drawing_parameters", {})
        if drawing_params:
            self.console.print("\nDrawing Parameters:")
            for param, value in drawing_params.items():
                display_param = param.replace("_", " ").title()
                self.console.print(f"  {display_param}: {value}")

        # Display coloring method with color id if applicable
        coloring_method = representation.get("coloring_method", "Default")
        self.console.print(f"\nColoring Method: {coloring_method}")
        if coloring_method == "ColorID" and "color_value" in representation:
            # Find the color name for this value
            color_name = "Unknown"
            for name, value in self.rep_manager.COLOR_IDS.items():
                if value == representation["color_value"]:
                    color_name = name
                    break
            self.console.print(f"Color: {color_name}")

        self.console.print(f"Material: {representation.get('material', 'Default')}")
        if "opacity" in representation:
            self.console.print(f"Opacity: {representation['opacity']}")

        self.console.print(f"Selection: {representation.get('selection', 'all')}")

        # What to edit
        self.console.print("\n[bold]What would you like to modify?[/bold]")
        self.console.print("1. Selection (what atoms are shown)", highlight=False)
        self.console.print("2. Drawing Method", highlight=False)
        self.console.print("3. Coloring Method", highlight=False)
        self.console.print("4. Material", highlight=False)
        self.console.print("5. Representation Name", highlight=False)
        self.console.print("6. All of the above", highlight=False)
        self.console.print("7. Cancel", highlight=False)

        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4", "5", "6", "7"],
            default="1",
            module="Structure Viewer - Representations",
            description="Representation field to edit",
            options_map={
                "1": "Selection (what atoms are shown)",
                "2": "Drawing Method",
                "3": "Coloring Method",
                "4": "Material",
                "5": "Representation Name",
                "6": "All of the above",
                "7": "Cancel",
            },
        )

        if choice == "7":
            return

        # Edit selection
        if choice in ["1", "6"]:
            self.console.print("\n[bold]Edit Selection:[/bold]")

            selection_builder = EnhancedSelectionBuilder(
                self.component_classifier, self.console, processor=self.processor
            )
            selection_builder.current_selection = representation.get("selection", "all")
            selection_builder.selection_text = "Current selection"

            selection = selection_builder.build_selection_interactive(self.structure)
            selection_text = selection_builder.selection_text

            representation["selection"] = selection

        # Edit drawing method
        if choice in ["2", "6"]:
            self.console.print("\n[bold]Edit Drawing Method:[/bold]")

            drawing_methods = self.rep_manager.get_drawing_method_options()
            current_method = representation.get("drawing_method")

            self.console.print("\n[bold]Available Drawing Methods:[/bold]")
            for i, method in enumerate(drawing_methods):
                if method == current_method:
                    self.console.print(f"{i+1}. {method} [green](Current)[/green]")
                else:
                    self.console.print(f"{i+1}. {method}")

            edit_drawing_map = {str(i + 1): m for i, m in enumerate(drawing_methods)}
            method_idx = prompt_with_context(
                self.processor,
                "Select drawing method",
                choices=[str(i + 1) for i in range(len(drawing_methods))],
                default=str(drawing_methods.index(current_method) + 1)
                if current_method in drawing_methods
                else "1",
                module="Structure Viewer - Representations",
                description="Drawing method for edit",
                options_map=edit_drawing_map,
            )

            drawing_method = drawing_methods[int(method_idx) - 1]
            representation["drawing_method"] = drawing_method

            # Update parameters for new drawing method
            drawing_params = self.rep_manager.get_drawing_method_parameters(
                drawing_method
            )
            representation["drawing_parameters"] = drawing_params

        # Edit coloring method
        if choice in ["3", "6"]:
            self.console.print("\n[bold]Edit Coloring Method:[/bold]")

            coloring_methods = self.rep_manager.get_coloring_method_options()
            current_coloring = representation.get("coloring_method")

            self.console.print("\n[bold]Available Coloring Methods:[/bold]")
            for i, method in enumerate(coloring_methods):
                if method == current_coloring:
                    self.console.print(f"{i+1}. {method} [green](Current)[/green]")
                else:
                    self.console.print(f"{i+1}. {method}")

            edit_coloring_map = {str(i + 1): m for i, m in enumerate(coloring_methods)}
            color_idx = prompt_with_context(
                self.processor,
                "Select coloring method",
                choices=[str(i + 1) for i in range(len(coloring_methods))],
                default=str(coloring_methods.index(current_coloring) + 1)
                if current_coloring in coloring_methods
                else "1",
                module="Structure Viewer - Representations",
                description="Coloring method for edit",
                options_map=edit_coloring_map,
            )

            coloring_method = coloring_methods[int(color_idx) - 1]
            representation["coloring_method"] = coloring_method

            # Handle ColorID
            if coloring_method == "ColorID":
                color_options = self.rep_manager.get_color_id_options()
                self.console.print("\n[bold]Available Colors:[/bold]")

                for i, color in enumerate(color_options):
                    self.console.print(f"{i+1}. {color}")

                edit_color_map = {str(i + 1): c for i, c in enumerate(color_options)}
                color_choice = prompt_with_context(
                    self.processor,
                    "Select color",
                    choices=[str(i + 1) for i in range(len(color_options))],
                    default="1",
                    module="Structure Viewer - Representations",
                    description="ColorID for edit",
                    options_map=edit_color_map,
                )

                selected_color = color_options[int(color_choice) - 1]
                color_value = self.rep_manager.get_color_id_value(selected_color)
                representation["color_value"] = color_value

        # Edit material
        if choice in ["4", "6"]:
            self.console.print("\n[bold]Edit Material:[/bold]")

            materials = self.rep_manager.get_material_options()
            current_material = representation.get("material")

            self.console.print("\n[bold]Available Materials:[/bold]")
            for i, material in enumerate(materials):
                if material == current_material:
                    self.console.print(f"{i+1}. {material} [green](Current)[/green]")
                else:
                    self.console.print(f"{i+1}. {material}")

            edit_material_map = {str(i + 1): m for i, m in enumerate(materials)}
            material_idx = prompt_with_context(
                self.processor,
                "Select material",
                choices=[str(i + 1) for i in range(len(materials))],
                default=str(materials.index(current_material) + 1)
                if current_material in materials
                else "1",
                module="Structure Viewer - Representations",
                description="Material for edit",
                options_map=edit_material_map,
            )

            material = materials[int(material_idx) - 1]
            representation["material"] = material

            # Handle opacity
            if material == "Transparent" or "Glass" in material:
                current_opacity = representation.get("opacity", 1.0)
                opacity = float(
                    prompt_with_context(
                        self.processor,
                        "Enter opacity (0.0-1.0)",
                        default=str(current_opacity),
                        module="Structure Viewer - Representations",
                        description="Opacity for transparent/glass material (edit)",
                    )
                )
                representation["opacity"] = opacity

        # Edit name
        if choice in ["5", "6"]:
            current_name = representation.get("selection_name", "Representation")
            new_name = prompt_with_context(
                self.processor,
                "Enter new name",
                default=current_name,
                module="Structure Viewer - Representations",
                description="New representation name",
            )
            representation["selection_name"] = new_name

        self.console.print("[green]Representation updated successfully[/green]")

    def _delete_representation(self):
        """Delete a representation"""
        if not self.representations:
            self.console.print("[yellow]No representations to delete[/yellow]")
            return

        # Display representations
        self._display_representations_summary()

        # Get representation to delete
        rep_delete_options_map = {
            str(i + 1): rep.get("selection_name", f"Rep {i + 1}")
            for i, rep in enumerate(self.representations)
        }
        rep_idx = prompt_with_context(
            self.processor,
            "\nEnter the number of the representation to delete",
            choices=[str(i + 1) for i in range(len(self.representations))],
            default="1",
            module="Structure Viewer - Representations",
            description="Select representation to delete",
            options_map=rep_delete_options_map,
        )

        rep_idx = int(rep_idx) - 1
        representation = self.representations[rep_idx]

        if confirm_with_context(
            self.processor,
            f"Delete '{representation.get('selection_name', f'Rep {rep_idx+1}')}' ?",
            default=False,
            module="Structure Viewer - Representations",
            description="Confirm delete representation",
        ):
            deleted_rep = self.representations.pop(rep_idx)
            self.console.print(
                f"[green]Deleted representation: {deleted_rep.get('selection_name', 'Unknown')}[/green]"
            )

    def _preview_representations(self):
        """Preview representations without generating files"""
        if not self.representations:
            self.console.print("[yellow]No representations to preview[/yellow]")
            return

        self.console.print("\n[bold]Representation Preview:[/bold]")

        for i, rep in enumerate(self.representations, 1):
            self.console.print(f"\n[bold cyan]Representation {i}:[/bold cyan]")

            panel_content = f"""[bold]Name:[/bold] {rep.get('selection_name', f'Rep {i}')}
[bold]Selection:[/bold] {rep.get('selection', 'all')}
[bold]Drawing:[/bold] {rep.get('drawing_method', 'Unknown')}
[bold]Coloring:[/bold] {rep.get('coloring_method', 'Default')}
[bold]Material:[/bold] {rep.get('material', 'Default')}"""

            if "opacity" in rep:
                panel_content += f"\n[bold]Opacity:[/bold] {rep['opacity']}"

            if "color_value" in rep:
                # Find color name
                color_name = "Unknown"
                for name, value in self.rep_manager.COLOR_IDS.items():
                    if value == rep["color_value"]:
                        color_name = name
                        break
                panel_content += f"\n[bold]Color:[/bold] {color_name}"

            self.console.print(Panel(panel_content))

    def _configure_drawing_parameters(self, drawing_method, current_params):
        """Configure drawing method parameters"""
        self.console.print(f"\n[bold]Configure {drawing_method} Parameters:[/bold]")

        new_params = current_params.copy()

        for param, current_value in current_params.items():
            display_param = param.replace("_", " ").title()
            new_value = prompt_with_context(
                self.processor,
                f"Enter {display_param}",
                default=str(current_value),
                module="Structure Viewer - Representations",
                description=f"Drawing-method parameter {display_param}",
            )

            # Try to convert to appropriate type
            try:
                if isinstance(current_value, float):
                    new_params[param] = float(new_value)
                elif isinstance(current_value, int):
                    new_params[param] = int(new_value)
                else:
                    new_params[param] = new_value
            except ValueError:
                self.console.print(f"[yellow]Invalid value, keeping default[/yellow]")

        return new_params

    def _recommend_coloring_method(self, selection, drawing_method):
        """Recommend coloring method based on selection and drawing method"""
        selection_lower = selection.lower()

        # Specific recommendations based on selection
        if "protein" in selection_lower:
            if "backbone" in selection_lower or drawing_method in [
                "Cartoon",
                "NewCartoon",
                "Ribbons",
            ]:
                return "Structure"  # Good for secondary structure
            else:
                return "ResType"  # Good for residue types

        elif "metal" in selection_lower or "name" in selection_lower:
            return "Element"  # Good for individual atoms

        elif "water" in selection_lower:
            return "ColorID"  # Single color for water

        elif "chain" in selection_lower:
            return "Chain"  # Different colors per chain

        else:
            # Default recommendations
            if drawing_method in ["VDW", "CPK", "Licorice"]:
                return "Element"
            elif drawing_method in ["Cartoon", "NewCartoon", "Ribbons"]:
                return "Structure"
            else:
                return "ResType"

    def _recommend_material(self, drawing_method):
        """Recommend material based on drawing method"""
        if drawing_method in ["QuickSurf", "MSMS", "Surf"]:
            return "Transparent"
        elif drawing_method in ["VDW", "CPK"]:
            return "Glossy"
        elif drawing_method in ["Cartoon", "NewCartoon"]:
            return "Opaque"
        else:
            return "Opaque"

    def _generate_representation_name(self, selection_text, drawing_method):
        """Generate a default name for the representation"""
        # Simplify selection text for name
        if "All atoms" in selection_text:
            base_name = "All"
        elif "Protein" in selection_text:
            base_name = "Protein"
        elif "Chain" in selection_text:
            base_name = selection_text.split()[1] if len(selection_text.split()) > 1 else "Chain"
        elif "Metal" in selection_text:
            base_name = "Metal"
        elif "Water" in selection_text:
            base_name = "Water"
        else:
            # Take first word or use "Selection"
            words = selection_text.split()
            base_name = words[0] if words else "Selection"

        return f"{base_name} {drawing_method}"


# ============================================================================
# ENHANCED WORKSPACE-AWARE VMD VISUALIZATION MODULE
# ============================================================================
# DEPRECATED: This module is kept in codebase for reference but is no longer
# accessible from menus. It has been replaced by InteractiveStructureViewer.
# The @register_module decorator has been removed to prevent conflicts.
# ============================================================================

# @register_module  # DEPRECATED - Commented out to disable menu access
class EnhancedVMDVisualizationModule(ProcessingModule):
    """
    Enhanced workspace-aware VMD visualization module.

    DEPRECATED: This module has been replaced by InteractiveStructureViewer
    which provides a modern web-based viewer without requiring VMD installation.
    This code is kept for reference only.
    """

    NAME = "Structure Viewer (VMD - DEPRECATED)"
    DESCRIPTION = "[DEPRECATED] Create intelligent VMD visualizations based on workspace analysis"
    VERSION = "2.0.0"
    CATEGORY = "visualization"
    REQUIRES = ["pdb_file"]
    PRIORITY = 10  # Run after analysis modules

    def initialize(self):
        """Initialize module resources"""
        self.console = Console()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Initialize with existing robust components
        self.rep_manager = RepresentationManager()
        self.state_generator = VMDStateGenerator()
        
        # Component classifier setup
        self.component_classifier = self._get_component_classifier()
        
        # Workspace analysis results
        self.workspace_features = {}
        self.available_structures = {}
        self.suggested_representations = []
        
        # User preferences
        self.visualization_mode = None
        self.selected_features = []

    def _get_component_classifier(self):
        """Get component classifier from existing infrastructure"""
        try:
            # Try to get the component classifier from the processor
            if hasattr(self, 'processor') and hasattr(self.processor, 'component_classifier'):
                return self.processor.component_classifier
            
            # Try to get from registry or create simple one
            from proprep.structure_prep.chem_comp_dict_fetcher import CCDParser
            
            class SimpleComponentClassifier:
                def __init__(self, ccd_parser=None):
                    self.ccd_parser = ccd_parser
                
                def classify_residue(self, residue, ccd_parser=None):
                    """Simple residue classification"""
                    res_name = residue.get_resname()
                    
                    # Standard amino acids
                    standard_aa = {
                        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                        "THR", "TRP", "TYR", "VAL"
                    }
                    
                    if res_name in standard_aa:
                        return "amino_acid"
                    elif res_name in {"HOH", "WAT", "TIP3", "TIP4"}:
                        return "water"
                    else:
                        return "hetero"
            
            return SimpleComponentClassifier()
            
        except Exception as e:
            self.logger.warning(f"Could not get component classifier: {e}")
            return None

    def get_menu_options(self) -> Dict[str, str]:
        """Get module menu options"""
        return {
            "create": "Create intelligent VMD visualization",
            "analyze": "Analyze workspace features",
            "presets": "Use visualization presets",
        }

    def handle_menu_option(self, option: str) -> bool:
        """Handle menu option selection"""
        if option == "create":
            self._run_intelligent_visualization_creator()
            return True
        elif option == "analyze":
            self._display_workspace_analysis()
            return True
        elif option == "presets":
            self._use_visualization_presets()
            return True
        return False

    def can_process(self, workspace) -> bool:
        """Check if module can process the current workspace"""
        from proprep.utils.structure_selector import StructureSelector

        selector = StructureSelector(workspace, self.console)
        status = selector.get_structure_status()
        return status.get("has_any", False)

    def _run_intelligent_visualization_creator(self):
        """Run the enhanced workspace-aware visualization creator"""
        self.console.print("\n[bold cyan]===== Intelligent VMD Visualization Creator =====[/bold cyan]")
        
        # Step 1: Analyze workspace content
        self.console.print("\n[bold]Step 1: Analyzing workspace content...[/bold]")
        self.workspace_features = self._analyze_workspace_features()
        self.available_structures = self.workspace_features['structures']
        
        if not any(self.available_structures.values()):
            self.console.print("[red]No structures available for visualization.[/red]")
            return
        
        # Display analysis summary
        self._display_feature_summary(self.workspace_features)
        
        # Step 2: Select visualization mode
        self.console.print("\n[bold]Step 2: Select visualization approach[/bold]")
        self.visualization_mode = self._select_visualization_mode()
        if not self.visualization_mode:
            return
        
        # Step 3: Choose features to highlight
        self.console.print("\n[bold]Step 3: Select features to highlight[/bold]")
        self.selected_features = self._select_features_to_highlight()
        
        # Step 4: Generate suggested representations
        self.console.print("\n[bold]Step 4: Generate representation suggestions[/bold]")
        self.suggested_representations = self._generate_feature_representations()
        
        # Step 5: Review and customize representations
        self.console.print("\n[bold]Step 5: Review and customize representations[/bold]")
        final_representations = self._review_and_customize_representations()
        
        if not final_representations:
            self.console.print("[yellow]No representations selected. Canceling visualization.[/yellow]")
            return
        
        # Step 6: Generate VMD files
        self.console.print("\n[bold]Step 6: Generate VMD visualization[/bold]")
        self._generate_workspace_aware_visualization(final_representations)

    def _analyze_workspace_features(self) -> Dict[str, Any]:
        """Analyze workspace to determine available visualization features"""
        features = {
            'structures': {},
            'mutations': {},
            'structural_issues': {},
            'chemical_features': {},
            'analysis_results': {}
        }
        
        # Structure analysis
        features['structures'] = {
            'original': self._get_workspace_structure('structure') is not None,
            'filtered': self._get_workspace_structure('filtered_structure') is not None,
            'repaired': self._get_workspace_structure('repaired_structure') is not None,
            'processed': self._get_workspace_structure('processed_structure') is not None
        }
        
        # Mutation analysis
        pending_mutations = self.get_from_workspace('pending_mutations', [])
        applied_mutations = self.get_from_workspace('mutations_applied', [])
        features['mutations'] = {
            'pending': len(pending_mutations) > 0,
            'applied': len(applied_mutations) > 0,
            'pending_count': len(pending_mutations),
            'applied_count': len(applied_mutations),
            'pending_data': pending_mutations,
            'applied_data': applied_mutations
        }
        
        # Structural completeness issues
        completeness_results = self.get_from_workspace('completeness_results', {})
        features['structural_issues'] = {
            'missing_residues': self._has_structural_issues(completeness_results, 'missing_residues'),
            'missing_atoms': self._has_structural_issues(completeness_results, 'missing_atoms'),
            'alternate_locations': self._has_structural_issues(completeness_results, 'alternate_locations'),
            'repaired_regions': self.get_from_workspace('repaired_structure') is not None
        }
        
        # Chemical features
        metal_sites = self.get_from_workspace('metal_sites', [])
        disulfide_bonds = self.get_from_workspace('disulfide_bonds', [])
        features['chemical_features'] = {
            'metal_sites': len(metal_sites) > 0,
            'disulfide_bonds': len(disulfide_bonds) > 0,
            'metal_sites_count': len(metal_sites) if isinstance(metal_sites, list) else 0,
            'disulfide_count': len(disulfide_bonds),
            'metal_sites_data': metal_sites,
            'disulfide_data': disulfide_bonds
        }
        
        # Analysis results
        protonation_results = self.get_from_workspace('protonation_results', {})
        filter_selections = self.get_from_workspace('filter_selections', {})
        features['analysis_results'] = {
            'protonation_states': len(protonation_results) > 0,
            'filtered_regions': len(filter_selections) > 0,
            'titratable_count': len(protonation_results),
            'protonation_data': protonation_results,
            'filter_data': filter_selections
        }
        
        return features

    def _display_feature_summary(self, features: Dict[str, Any]):
        """Display a summary of detected workspace features"""
        table = Table(title="Workspace Analysis Summary", show_header=True)
        table.add_column("Category", style="cyan", width=20)
        table.add_column("Feature", style="green", width=25)
        table.add_column("Status", style="yellow", width=15)
        table.add_column("Details", style="white", width=30)
        
        # Structures
        for struct_type, available in features['structures'].items():
            status = "✓ Available" if available else "✗ Not available"
            table.add_row("Structures", struct_type.replace('_', ' ').title(), status, "")
        
        # Mutations
        if features['mutations']['pending'] or features['mutations']['applied']:
            if features['mutations']['pending']:
                table.add_row("Mutations", "Pending", "✓ Detected", 
                            f"{features['mutations']['pending_count']} mutations")
            if features['mutations']['applied']:
                table.add_row("Mutations", "Applied", "✓ Detected",
                            f"{features['mutations']['applied_count']} mutations")
        else:
            table.add_row("Mutations", "None", "✗ Not detected", "")
        
        # Chemical features
        if features['chemical_features']['metal_sites']:
            table.add_row("Chemical", "Metal Sites", "✓ Detected",
                        f"{features['chemical_features']['metal_sites_count']} sites")
        
        if features['chemical_features']['disulfide_bonds']:
            table.add_row("Chemical", "Disulfide Bonds", "✓ Detected",
                        f"{features['chemical_features']['disulfide_count']} bonds")
        
        # Structural issues
        if features['structural_issues']['missing_residues']:
            table.add_row("Issues", "Missing Residues", "⚠ Detected", "")
        if features['structural_issues']['missing_atoms']:
            table.add_row("Issues", "Missing Atoms", "⚠ Detected", "")
        if features['structural_issues']['repaired_regions']:
            table.add_row("Issues", "Repaired Regions", "✓ Available", "")
        
        # Analysis results
        if features['analysis_results']['protonation_states']:
            table.add_row("Analysis", "Protonation States", "✓ Analyzed",
                        f"{features['analysis_results']['titratable_count']} residues")
        
        if features['analysis_results']['filtered_regions']:
            table.add_row("Analysis", "Filtered Regions", "✓ Available", "")
        
        self.console.print(table)

    def _select_visualization_mode(self) -> Optional[str]:
        """Select visualization mode based on available structures"""
        available = self.available_structures
        modes = []
        descriptions = []
        
        if available['original']:
            modes.append("original_only")
            descriptions.append("Show original structure only")
            
            if available['filtered']:
                modes.append("original_with_highlights")
                descriptions.append("Show original structure with filtered regions highlighted")
        
        if available['filtered']:
            modes.append("filtered_only")
            descriptions.append("Show filtered structure only")
        
        if available['repaired']:
            modes.append("repaired_only")
            descriptions.append("Show repaired structure only")
            
            if available['filtered']:
                modes.append("filtered_with_repairs")
                descriptions.append("Show filtered structure with repaired regions highlighted")
            
            if available['original']:
                modes.append("before_after_repair")
                descriptions.append("Compare original vs repaired structures")
        
        if len(modes) == 0:
            self.console.print("[red]No valid visualization modes available.[/red]")
            return None
        
        if len(modes) == 1:
            self.console.print(f"[green]Using available mode: {descriptions[0]}[/green]")
            return modes[0]
        
        # Multiple modes available, let user choose
        self.console.print("\n[bold]Available visualization modes:[/bold]")
        for i, (mode, desc) in enumerate(zip(modes, descriptions), 1):
            self.console.print(f"{i}. {desc}")
        
        mode_options_map = {str(i + 1): desc for i, desc in enumerate(descriptions)}
        choice = prompt_with_context(
            self.processor,
            "Select mode",
            choices=[str(i) for i in range(1, len(modes) + 1)],
            default="1",
            module="Structure Viewer",
            description="Select visualization mode",
            options_map=mode_options_map,
        )

        return modes[int(choice) - 1]

    def _select_features_to_highlight(self) -> List[str]:
        """Select which workspace features to highlight in visualization"""
        available_features = []
        feature_descriptions = {}
        
        # Check each feature type
        if self.workspace_features['mutations']['pending']:
            available_features.append('pending_mutations')
            feature_descriptions['pending_mutations'] = f"Pending mutations ({self.workspace_features['mutations']['pending_count']} residues)"
        
        if self.workspace_features['mutations']['applied']:
            available_features.append('applied_mutations')
            feature_descriptions['applied_mutations'] = f"Applied mutations ({self.workspace_features['mutations']['applied_count']} residues)"
        
        if self.workspace_features['chemical_features']['metal_sites']:
            available_features.append('metal_sites')
            feature_descriptions['metal_sites'] = f"Metal binding sites ({self.workspace_features['chemical_features']['metal_sites_count']} sites)"
        
        if self.workspace_features['chemical_features']['disulfide_bonds']:
            available_features.append('disulfide_bonds')
            feature_descriptions['disulfide_bonds'] = f"Disulfide bonds ({self.workspace_features['chemical_features']['disulfide_count']} bonds)"
        
        if self.workspace_features['analysis_results']['protonation_states']:
            available_features.append('titratable_residues')
            feature_descriptions['titratable_residues'] = f"Titratable residues ({self.workspace_features['analysis_results']['titratable_count']} residues)"
        
        if self.workspace_features['analysis_results']['filtered_regions']:
            available_features.append('filtered_regions')
            feature_descriptions['filtered_regions'] = "Filtered regions"
        
        if self.workspace_features['structural_issues']['repaired_regions']:
            available_features.append('repaired_regions')
            feature_descriptions['repaired_regions'] = "Repaired regions"
        
        if not available_features:
            self.console.print("[yellow]No special features detected for highlighting.[/yellow]")
            return []
        
        # Display available features
        self.console.print("\n[bold]Available features to highlight:[/bold]")
        for i, feature in enumerate(available_features, 1):
            self.console.print(f"{i}. {feature_descriptions[feature]}")
        
        self.console.print(f"{len(available_features) + 1}. Select all")
        self.console.print(f"{len(available_features) + 2}. None (just show basic structure)")
        
        # Get user selection
        if confirm_with_context(
            self.processor,
            "Select features interactively?",
            default=True,
            module="Structure Viewer",
            description="Select features to highlight interactively (per-feature yes/no)",
        ):
            selected = []
            for i, feature in enumerate(available_features, 1):
                if confirm_with_context(
                    self.processor,
                    f"Highlight {feature_descriptions[feature]}?",
                    default=True,
                    module="Structure Viewer",
                    description=f"Highlight {feature}",
                ):
                    selected.append(feature)
            return selected
        else:
            feature_opts = {str(i + 1): feature_descriptions[f] for i, f in enumerate(available_features)}
            feature_opts[str(len(available_features) + 1)] = "Select all"
            feature_opts[str(len(available_features) + 2)] = "None (just show basic structure)"
            choice = prompt_with_context(
                self.processor,
                "Select option",
                choices=[str(i) for i in range(1, len(available_features) + 3)],
                default=str(len(available_features) + 1),
                module="Structure Viewer",
                description="Select features to highlight",
                options_map=feature_opts,
            )
            
            choice_num = int(choice)
            if choice_num == len(available_features) + 1:  # Select all
                return available_features
            elif choice_num == len(available_features) + 2:  # None
                return []
            else:  # Individual feature
                return [available_features[choice_num - 1]]

    def _generate_feature_representations(self) -> List[Dict[str, Any]]:
        """Generate suggested representations for selected features"""
        representations = []
        
        # Always add a basic protein representation
        basic_rep = self._create_basic_protein_representation()
        if basic_rep:
            representations.append(basic_rep)
        
        # Add feature-specific representations
        for feature in self.selected_features:
            self.logger.info(f"Generating representation for feature: {feature}")
            
            if feature == 'pending_mutations':
                rep = self._create_mutation_representation('pending')
                if rep:
                    representations.append(rep)
                    self.logger.info(f"Created pending mutations representation with selection: {rep['selection']}")
                else:
                    self.logger.warning("Could not create pending mutations representation")
            
            elif feature == 'applied_mutations':
                rep = self._create_mutation_representation('applied')
                if rep:
                    representations.append(rep)
                    self.logger.info(f"Created applied mutations representation with selection: {rep['selection']}")
                else:
                    self.logger.warning("Could not create applied mutations representation")
            
            elif feature == 'metal_sites':
                reps = self._create_metal_site_representations()
                if reps:
                    representations.extend(reps)
                    for rep in reps:
                        self.logger.info(f"Created metal site representation '{rep['selection_name']}' with selection: {rep['selection']}")
                else:
                    self.logger.warning("Could not create metal site representations")
            
            elif feature == 'disulfide_bonds':
                rep = self._create_disulfide_representation()
                if rep:
                    representations.append(rep)
                    self.logger.info(f"Created disulfide bonds representation with selection: {rep['selection']}")
                else:
                    self.logger.warning("Could not create disulfide bonds representation")
            
            elif feature == 'titratable_residues':
                rep = self._create_protonation_representation()
                if rep:
                    representations.append(rep)
                    self.logger.info(f"Created titratable residues representation with selection: {rep['selection']}")
                else:
                    self.logger.warning("Could not create titratable residues representation")
            
            elif feature == 'filtered_regions':
                rep = self._create_filtered_regions_representation()
                if rep:
                    representations.append(rep)
                    self.logger.info(f"Created filtered regions representation with selection: {rep['selection']}")
                else:
                    self.logger.warning("Could not create filtered regions representation")
            
            elif feature == 'repaired_regions':
                rep = self._create_repaired_regions_representation()
                if rep:
                    representations.append(rep)
                    self.logger.info(f"Created repaired regions representation with selection: {rep['selection']}")
                else:
                    self.logger.info("Skipped repaired regions representation (no specific selection available)")
        
        self.logger.info(f"Generated {len(representations)} total representations")
        return representations

    def _create_basic_protein_representation(self) -> Optional[Dict[str, Any]]:
        """Create basic protein backbone representation"""
        # Get proper drawing parameters from RepresentationManager
        drawing_params = self.rep_manager.get_drawing_method_parameters('NewCartoon')
        
        return {
            'selection_name': 'Protein Backbone',
            'description': 'Main protein structure',
            'selection': 'protein',
            'drawing_method': 'NewCartoon',
            'drawing_parameters': drawing_params,
            'coloring_method': 'ColorID',
            'color_value': 'Silver',  # ColorID 6
            'material': 'Opaque',
            'suggested': True,
            'customizable': True
        }

    def _create_mutation_representation(self, mutation_type: str) -> Optional[Dict[str, Any]]:
        """Create representation for mutations"""
        mutations_data = self.workspace_features['mutations'][f'{mutation_type}_data']
        
        if not mutations_data:
            return None
        
        # Build VMD selection for mutations
        selection_parts = []
        for chain_id, res_num, from_aa, to_aa in mutations_data:
            selection_parts.append(f"(chain {chain_id} and resid {res_num})")
        
        if not selection_parts:
            return None
        
        color = 'Orange' if mutation_type == 'pending' else 'Lime'  # Orange for pending, Lime for applied
        
        # Get proper drawing parameters from RepresentationManager
        drawing_params = self.rep_manager.get_drawing_method_parameters('NewCartoon')
        
        return {
            'selection_name': f'{mutation_type.title()} Mutations',
            'description': f'Residues with {mutation_type} mutations ({len(mutations_data)} total)',
            'selection': ' or '.join(selection_parts),
            'drawing_method': 'NewCartoon',
            'drawing_parameters': drawing_params,
            'coloring_method': 'ColorID',
            'color_value': color,
            'material': 'Opaque',
            'suggested': True,
            'customizable': True
        }

    def _create_metal_site_representations(self) -> List[Dict[str, Any]]:
        """Create representations for metal sites"""
        metal_sites_data = self.workspace_features['chemical_features']['metal_sites_data']
        
        if not metal_sites_data:
            return []
        
        representations = []
        
        # Build selections for metal centers and ligands
        metal_selections = []
        ligand_selections = []
        
        for site in metal_sites_data:
            try:
                # Handle different metal site data structures
                if isinstance(site, dict):
                    metal_info = site.get('metal', {})
                    ligands = site.get('ligands', [])
                    
                    # Metal center - be more flexible with field names
                    if 'chain' in metal_info and ('resid' in metal_info or 'residue_number' in metal_info):
                        chain = metal_info['chain']
                        resid = metal_info.get('resid', metal_info.get('residue_number'))
                        if resid is not None:
                            metal_selections.append(f"(chain {chain} and resid {resid})")
                    
                    # Ligands
                    for ligand in ligands:
                        if isinstance(ligand, dict):
                            if 'chain' in ligand and ('resid' in ligand or 'residue_number' in ligand):
                                chain = ligand['chain']
                                resid = ligand.get('resid', ligand.get('residue_number'))
                                if resid is not None:
                                    ligand_selections.append(f"(chain {chain} and resid {resid})")
                
                elif hasattr(site, 'metal') and hasattr(site, 'ligands'):
                    # Object-style access
                    metal = site.metal
                    ligands = site.ligands
                    
                    # Metal center
                    if hasattr(metal, 'chain') and (hasattr(metal, 'resid') or hasattr(metal, 'residue_number')):
                        chain = metal.chain
                        resid = getattr(metal, 'resid', getattr(metal, 'residue_number', None))
                        if resid is not None:
                            metal_selections.append(f"(chain {chain} and resid {resid})")
                    
                    # Ligands
                    for ligand in ligands:
                        if hasattr(ligand, 'chain') and (hasattr(ligand, 'resid') or hasattr(ligand, 'residue_number')):
                            chain = ligand.chain
                            resid = getattr(ligand, 'resid', getattr(ligand, 'residue_number', None))
                            if resid is not None:
                                ligand_selections.append(f"(chain {ligand.chain} and resid {resid})")
            
            except Exception as e:
                self.logger.warning(f"Error processing metal site data: {e}")
                continue
        
        # Metal centers representation
        if metal_selections:
            drawing_params = self.rep_manager.get_drawing_method_parameters('VDW')
            representations.append({
                'selection_name': 'Metal Centers',
                'description': f'Metal ions ({len(metal_selections)} centers)',
                'selection': ' or '.join(metal_selections),
                'drawing_method': 'VDW',
                'drawing_parameters': drawing_params,
                'coloring_method': 'Element',
                'material': 'Glossy',
                'suggested': True,
                'customizable': True
            })
        else:
            self.logger.warning("No valid metal center selections could be built")
        
        # Metal ligands representation
        if ligand_selections:
            drawing_params = self.rep_manager.get_drawing_method_parameters('Licorice')
            representations.append({
                'selection_name': 'Metal Ligands',
                'description': f'Metal coordination residues ({len(ligand_selections)} residues)',
                'selection': ' or '.join(ligand_selections),
                'drawing_method': 'Licorice',
                'drawing_parameters': drawing_params,
                'coloring_method': 'Element',
                'material': 'Opaque',
                'suggested': True,
                'customizable': True
            })
        else:
            self.logger.warning("No valid metal ligand selections could be built")
        
        return representations

    def _create_disulfide_representation(self) -> Optional[Dict[str, Any]]:
        """Create representation for disulfide bonds"""
        disulfide_data = self.workspace_features['chemical_features']['disulfide_data']
        
        if not disulfide_data:
            return None
        
        # Build selection for cysteines involved in disulfide bonds
        selection_parts = []
        for chain1, res1, chain2, res2 in disulfide_data:
            selection_parts.extend([
                f"(chain {chain1} and resid {res1})",
                f"(chain {chain2} and resid {res2})"
            ])
        
        if not selection_parts:
            return None
        
        # Get proper drawing parameters from RepresentationManager
        drawing_params = self.rep_manager.get_drawing_method_parameters('NewCartoon')
        
        return {
            'selection_name': 'Disulfide Bonds',
            'description': f'Cysteine residues in disulfide bonds ({len(disulfide_data)} bonds)',
            'selection': ' or '.join(selection_parts),
            'drawing_method': 'NewCartoon',
            'drawing_parameters': drawing_params,
            'coloring_method': 'ColorID',
            'color_value': 'Yellow',  # Yellow for disulfide bonds
            'material': 'Opaque',
            'suggested': True,
            'customizable': True
        }

    def _create_protonation_representation(self) -> Optional[Dict[str, Any]]:
        """Create representation for titratable residues"""
        protonation_data = self.workspace_features['analysis_results']['protonation_data']
        
        if not protonation_data:
            return None
        
        # Extract titratable residues from results - handle different data structures
        selection_parts = []
        
        # Try different possible data structures for protonation results
        if isinstance(protonation_data, dict):
            for residue_key, residue_data in protonation_data.items():
                if isinstance(residue_data, dict):
                    # Format 1: {residue_key: {chain: X, resid: Y, ...}}
                    if 'chain' in residue_data and 'resid' in residue_data:
                        chain = residue_data['chain']
                        resid = residue_data['resid']
                        selection_parts.append(f"(chain {chain} and resid {resid})")
                    # Format 2: nested structure analysis
                    elif 'residue_info' in residue_data:
                        res_info = residue_data['residue_info']
                        if 'chain' in res_info and 'resid' in res_info:
                            chain = res_info['chain']
                            resid = res_info['resid']
                            selection_parts.append(f"(chain {chain} and resid {resid})")
                
                # Format 3: residue_key contains chain and resid (e.g., "A_123_HIS")
                elif isinstance(residue_key, str) and '_' in residue_key:
                    parts = residue_key.split('_')
                    if len(parts) >= 2:
                        chain = parts[0]
                        try:
                            resid = int(parts[1])
                            selection_parts.append(f"(chain {chain} and resid {resid})")
                        except ValueError:
                            continue
        
        if not selection_parts:
            # Fallback: if we can't parse specific residues, don't create the representation
            self.logger.warning("Could not parse specific residues from protonation results")
            return None
        
        # Get proper drawing parameters from RepresentationManager
        drawing_params = self.rep_manager.get_drawing_method_parameters('NewCartoon')
        
        return {
            'selection_name': 'Titratable Residues',
            'description': f'Residues with analyzed protonation states ({len(selection_parts)} residues)',
            'selection': ' or '.join(selection_parts),
            'drawing_method': 'NewCartoon',
            'drawing_parameters': drawing_params,
            'coloring_method': 'ColorID',
            'color_value': 'Purple',  # Purple for titratable residues
            'material': 'Opaque',
            'suggested': True,
            'customizable': True
        }

    def _create_filtered_regions_representation(self) -> Optional[Dict[str, Any]]:
        """Create representation for filtered regions"""
        filter_data = self.workspace_features['analysis_results']['filter_data']
        
        if not filter_data:
            return None
        
        # Build selection from filter data
        selection_parts = []
        for chain_id, chain_data in filter_data.items():
            for comp_type, residue_list in chain_data.items():
                for res_id in residue_list:
                    selection_parts.append(f"(chain {chain_id} and resid {res_id})")
        
        if not selection_parts:
            return None
        
        # Get proper drawing parameters from RepresentationManager
        drawing_params = self.rep_manager.get_drawing_method_parameters('NewCartoon')
        
        return {
            'selection_name': 'Filtered Regions',
            'description': 'Regions selected by PDB Filter',
            'selection': ' or '.join(selection_parts),
            'drawing_method': 'NewCartoon',
            'drawing_parameters': drawing_params,
            'coloring_method': 'ColorID',
            'color_value': 'Cyan',  # Cyan for filtered regions
            'material': 'Opaque',
            'suggested': True,
            'customizable': True
        }

    def _create_repaired_regions_representation(self) -> Optional[Dict[str, Any]]:
        """Create representation for repaired regions using proper residue mapping"""
        # Get the completeness results to understand what was repaired
        completeness_results = self.get_from_workspace('completeness_results', {})
        
        if not completeness_results:
            return None
        
        # Check if we have residue mapping information from MODELLER
        residue_mapping = self.get_from_workspace('residue_mapping', {})
        
        # Extract specific residues that were identified as missing and presumably repaired
        repaired_selections = []
        
        # Get missing residues that were repaired
        missing_residues = completeness_results.get('missing_residues', {})
        original_missing_residues = set()
        
        for method, chain_results in missing_residues.items():
            if isinstance(chain_results, dict):
                for chain_id, missing_list in chain_results.items():
                    if isinstance(missing_list, list):
                        for missing_res in missing_list:
                            # Handle different formats of missing residue data
                            if isinstance(missing_res, dict):
                                if 'residue_number' in missing_res:
                                    res_num = missing_res['residue_number']
                                    original_missing_residues.add((chain_id, res_num))
                                elif 'resid' in missing_res:
                                    res_num = missing_res['resid']
                                    original_missing_residues.add((chain_id, res_num))
                            elif isinstance(missing_res, (int, str)):
                                # Simple residue number
                                original_missing_residues.add((chain_id, int(missing_res)))
        
        # If we have residue mapping, use it to convert to new residue numbers
        if residue_mapping:
            self.logger.info("Using residue mapping to convert original to MODELLER-renumbered residues")
            for chain_id in residue_mapping:
                if 'original_to_new' in residue_mapping[chain_id]:
                    original_to_new = residue_mapping[chain_id]['original_to_new']
                    
                    # Convert original missing residues to new numbering
                    for orig_chain, orig_res in original_missing_residues:
                        if orig_chain == chain_id and orig_res in original_to_new:
                            new_res = original_to_new[orig_res]
                            repaired_selections.append(f"(chain {chain_id} and resid {new_res})")
                            self.logger.info(f"Mapped original residue {orig_chain}:{orig_res} -> {chain_id}:{new_res}")
        else:
            # Fallback: try to identify repaired regions by comparing structures
            self.logger.info("No residue mapping found, attempting structure comparison approach")
            repaired_selections = self._identify_repaired_regions_by_comparison()
        
        # Remove duplicates
        repaired_selections = list(set(repaired_selections))
        
        if not repaired_selections:
            self.logger.warning("No specific repaired regions could be identified")
            return None
        
        # Get proper drawing parameters from RepresentationManager
        drawing_params = self.rep_manager.get_drawing_method_parameters('NewCartoon')
        
        return {
            'selection_name': 'Repaired Regions',
            'description': f'Regions repaired by Structure Fixer ({len(repaired_selections)} regions)',
            'selection': ' or '.join(repaired_selections),
            'drawing_method': 'NewCartoon',
            'drawing_parameters': drawing_params,
            'coloring_method': 'ColorID',
            'color_value': 'Green',  # Green for repaired regions
            'material': 'Opaque',
            'suggested': True,
            'customizable': True
        }

    def _identify_repaired_regions_by_comparison(self) -> List[str]:
        """
        Fallback method to identify repaired regions by comparing structures
        """
        original_structure = self.get_from_workspace('structure')
        filtered_structure = self.get_from_workspace('filtered_structure') 
        repaired_structure = self.get_from_workspace('repaired_structure')
        
        if not all([original_structure, filtered_structure, repaired_structure]):
            return []
        
        repaired_selections = []
        
        try:
            # Compare filtered vs repaired to identify added residues
            # This is a simplified approach - in reality this would need more sophisticated analysis
            
            # Get residue counts per chain
            filtered_residues = {}
            repaired_residues = {}
            
            # Count residues in filtered structure
            for model in filtered_structure:
                for chain in model:
                    chain_id = chain.id
                    filtered_residues[chain_id] = len([r for r in chain if r.id[0] == ' '])
            
            # Count residues in repaired structure  
            for model in repaired_structure:
                for chain in model:
                    chain_id = chain.id
                    repaired_residues[chain_id] = len([r for r in chain if r.id[0] == ' '])
            
            # If repaired has more residues, assume the extra ones are repairs
            # This is crude but better than nothing
            for chain_id in repaired_residues:
                if chain_id in filtered_residues:
                    if repaired_residues[chain_id] > filtered_residues[chain_id]:
                        extra_residues = repaired_residues[chain_id] - filtered_residues[chain_id]
                        self.logger.info(f"Chain {chain_id}: {extra_residues} residues appear to be repaired")
                        
                        # Since MODELLER renumbers sequentially, we can't easily identify specific residues
                        # For now, skip this fallback as it would be inaccurate
                        pass
            
        except Exception as e:
            self.logger.warning(f"Error in structure comparison: {e}")
        
        return repaired_selections

    def _review_and_customize_representations(self) -> List[Dict[str, Any]]:
        """Review suggested representations and allow customization"""
        if not self.suggested_representations:
            self.console.print("[yellow]No representations were generated.[/yellow]")
            return []
        
        self.console.print(f"\n[bold]Generated {len(self.suggested_representations)} representations:[/bold]")
        
        # Display suggested representations
        table = Table(title="Suggested Representations")
        table.add_column("#", width=3)
        table.add_column("Name", width=20)
        table.add_column("Description", width=35)
        table.add_column("Style", width=25)
        
        for i, rep in enumerate(self.suggested_representations, 1):
            style_info = f"{rep['drawing_method']}, {rep['coloring_method']}"
            if 'color_value' in rep:
                style_info += f" ({rep['color_value']})"
            
            table.add_row(str(i), rep['selection_name'], rep['description'], style_info)
        
        self.console.print(table)
        
        # Ask user how to proceed
        self.console.print("\n[bold]How would you like to proceed?[/bold]")
        self.console.print("1. Accept all suggestions as-is", highlight=False)
        self.console.print("2. Select which representations to include", highlight=False)
        self.console.print("3. Use existing representations manager for full customization", highlight=False)
        self.console.print("4. Cancel", highlight=False)
        
        choice = prompt_with_context(
            self.processor,
            "Choose option",
            choices=["1", "2", "3", "4"],
            default="1",
            module="Structure Viewer",
            description="Suggested-representations action",
            options_map={
                "1": "Accept all suggestions as-is",
                "2": "Select which representations to include",
                "3": "Use existing representations manager for full customization",
                "4": "Cancel",
            },
        )

        if choice == "1":
            return self.suggested_representations
        elif choice == "2":
            return self._select_representations_subset()
        elif choice == "3":
            return self._use_representations_manager()
        else:
            return []

    def _select_representations_subset(self) -> List[Dict[str, Any]]:
        """Allow user to select which representations to include"""
        selected_reps = []

        for i, rep in enumerate(self.suggested_representations):
            if confirm_with_context(
                self.processor,
                f"Include '{rep['selection_name']}'?",
                default=True,
                module="Structure Viewer",
                description=f"Include representation '{rep['selection_name']}'",
            ):
                selected_reps.append(rep)
        
        return selected_reps

    def _use_representations_manager(self) -> List[Dict[str, Any]]:
        """Use the existing RepresentationsManager for full customization"""
        # Get the working structure for the representations manager
        working_structure = self._get_working_structure()
        
        if not working_structure:
            self.console.print("[red]Could not get structure for representations manager.[/red]")
            return self.suggested_representations
        
        # Create representations manager with existing infrastructure
        rep_manager = RepresentationsManager(
            console=self.console,
            component_classifier=self.component_classifier,
            structure=working_structure,
            processor=self.processor,
        )
        
        # Pre-populate with suggested representations
        for rep in self.suggested_representations:
            rep_manager.representations.append(rep)
        
        # Run the interactive manager
        self.console.print("\n[bold]Opening full representations manager...[/bold]")
        self.console.print("Pre-loaded with suggested representations. You can modify, add, or remove as needed.")
        
        # Run the representations manager
        final_reps = rep_manager.run_representations_manager()
        
        return final_reps if final_reps else self.suggested_representations

    def _generate_workspace_aware_visualization(self, representations: List[Dict[str, Any]]):
        """Generate VMD files with workspace-aware features"""
        # Determine which structure to use
        working_structure = self._get_working_structure()
        
        if not working_structure:
            self.console.print("[red]Could not determine structure for visualization.[/red]")
            return
        
        # Generate output filename
        output_base = prompt_with_context(
            self.processor,
            "Enter output filename (without extension)",
            default="workspace_visualization",
            module="Structure Viewer",
            description="VMD output filename (without extension)",
        )
        
        # Add workspace context to representations
        enhanced_reps = self._enhance_representations_with_context(representations)
        
        # Use existing VMDStateGenerator to create the visualization
        try:
            pdb_file, state_file = self.state_generator.generate_state_file(
                filepath=output_base,
                structure=working_structure,
                representations=enhanced_reps,
                view_settings={'background': 'white', 'axes': True}
            )
            
            self.console.print(f"\n[green]Generated VMD visualization files:[/green]")
            self.console.print(f"  PDB Structure: {pdb_file}")
            self.console.print(f"  VMD State: {state_file}")
            
            # Add workspace information as comments
            self._add_workspace_documentation(state_file)
            
            self.console.print(f"\n[bold cyan]To use:[/bold cyan]")
            self.console.print(f"  1. Open VMD")
            self.console.print(f"  2. Go to File > Load Visualization State")
            self.console.print(f"  3. Select: {state_file}")
            
        except Exception as e:
            self.console.print(f"[red]Error generating VMD files: {str(e)}[/red]")
            self.logger.error(f"VMD generation error: {str(e)}")

    def _get_working_structure(self):
        """Get the appropriate structure for the selected visualization mode.

        Uses StructureSelector's get_structure_by_key() for mode-based selection,
        which is appropriate when the user has chosen a specific pipeline stage to visualize.
        """
        from proprep.utils.structure_selector import StructureSelector

        # Create selector if processor has workspace
        if hasattr(self, 'processor') and self.processor and hasattr(self.processor, 'workspace'):
            selector = StructureSelector(
                self.processor.workspace, self.console, processor=self.processor
            )

            if self.visualization_mode in ['original_only', 'original_with_highlights']:
                # Get original/loaded structure
                return selector.get_structure_by_key(
                    'structure',
                    fallback_keys=['rcsb_structure', 'local_structure', 'alphafold_structure'],
                    require_exists=False  # BioPython objects don't need file exists check
                )

            elif self.visualization_mode in ['filtered_only', 'filtered_with_repairs']:
                # For filtered with repairs, prefer repaired structure
                if self.visualization_mode == 'filtered_with_repairs':
                    structure = selector.get_structure_by_key(
                        'repaired_structure',
                        require_exists=False
                    )
                    if structure:
                        return structure
                # Fallback to filtered structure
                return selector.get_structure_by_key(
                    'filtered_structure',
                    require_exists=False
                )

            elif self.visualization_mode == 'repaired_only':
                return selector.get_structure_by_key(
                    'repaired_structure',
                    require_exists=False
                )

        # Default: fallback to direct workspace access
        return self.get_from_workspace('structure')

    def _enhance_representations_with_context(self, representations: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Add workspace context information to representations"""
        enhanced = []
        
        for rep in representations:
            enhanced_rep = rep.copy()
            
            # Add feature source information
            if 'suggested' in rep:
                enhanced_rep['comment'] = f"Generated from workspace feature: {rep.get('description', 'Unknown')}"
            
            # Ensure all required fields are present for VMDStateGenerator
            required_fields = ['drawing_method', 'coloring_method', 'material', 'selection', 'drawing_parameters']
            for field in required_fields:
                if field not in enhanced_rep:
                    # Set reasonable defaults
                    if field == 'drawing_parameters':
                        # Get default parameters for the drawing method
                        drawing_method = enhanced_rep.get('drawing_method', 'Cartoon')
                        enhanced_rep[field] = self.rep_manager.get_drawing_method_parameters(drawing_method)
                    else:
                        defaults = {
                            'drawing_method': 'Cartoon',
                            'coloring_method': 'Structure',
                            'material': 'Opaque',
                            'selection': 'all'
                        }
                        enhanced_rep[field] = defaults[field]
            
            enhanced.append(enhanced_rep)
        
        return enhanced

    def _add_workspace_documentation(self, state_file: str):
        """Add workspace context as comments to the VMD state file"""
        try:
            # Read existing file
            with open(state_file, 'r') as f:
                content = f.read()
            
            # Create documentation header
            doc_header = f"""
# Workspace Analysis Context:
# Visualization Mode: {self.visualization_mode}
# Selected Features: {', '.join(self.selected_features) if self.selected_features else 'None'}
# 
# Available Structures:
"""
            
            for struct_type, available in self.available_structures.items():
                doc_header += f"#   {struct_type}: {'Yes' if available else 'No'}\n"
            
            doc_header += "#\n# Generated Features:\n"
            for rep in self.suggested_representations:
                doc_header += f"#   - {rep.get('selection_name', 'Unknown')}: {rep.get('description', 'No description')}\n"
            
            doc_header += "#\n"
            
            # Insert documentation after the existing header
            lines = content.split('\n')
            header_end = 0
            for i, line in enumerate(lines):
                if line.startswith('#') or line.strip() == '':
                    header_end = i + 1
                else:
                    break
            
            # Insert our documentation
            lines.insert(header_end, doc_header)
            
            # Write back to file
            with open(state_file, 'w') as f:
                f.write('\n'.join(lines))
                
        except Exception as e:
            self.logger.warning(f"Could not add workspace documentation: {str(e)}")

    # Helper methods
    def _get_workspace_structure(self, key: str):
        """Get structure from workspace with error handling"""
        try:
            return self.get_from_workspace(key)
        except Exception:
            return None

    def _has_structural_issues(self, completeness_results: Dict[str, Any], issue_type: str) -> bool:
        """Check if completeness results contain specific structural issues"""
        if not completeness_results or issue_type not in completeness_results:
            return False
        
        issue_data = completeness_results[issue_type]
        if isinstance(issue_data, dict):
            # Check if any method found issues
            for method, results in issue_data.items():
                if isinstance(results, dict) and any(results.values()):
                    return True
        
        return False

    def get_from_workspace(self, key: str, default=None):
        """Get value from processor workspace"""
        if hasattr(self, 'processor') and self.processor:
            if hasattr(self.processor, 'workspace'):
                if hasattr(self.processor.workspace, 'get'):
                    return self.processor.workspace.get(key, default)
                elif hasattr(self.processor.workspace, '_data'):
                    return self.processor.workspace._data.get(key, default)
        
        return default

    def _display_workspace_analysis(self):
        """Display detailed workspace analysis"""
        features = self._analyze_workspace_features()
        self._display_feature_summary(features)

    def _use_visualization_presets(self):
        """Use predefined visualization presets"""
        self.console.print("\n[bold]Available Visualization Presets:[/bold]")
        
        presets = [
            ("metal_analysis", "Metal Site Analysis", "Focus on metal binding sites and coordination"),
            ("mutation_analysis", "Mutation Analysis", "Highlight mutation sites and changes"),
            ("structure_quality", "Structure Quality", "Show structural issues and repairs"),
            ("comprehensive", "Comprehensive View", "Show all detected features")
        ]
        
        for i, (preset_id, name, desc) in enumerate(presets, 1):
            self.console.print(f"{i}. {name}: {desc}")
        
        preset_options_map = {str(i + 1): f"{name}: {desc}" for i, (_, name, desc) in enumerate(presets)}
        choice = prompt_with_context(
            self.processor,
            "Select preset",
            choices=[str(i) for i in range(1, len(presets) + 1)],
            default="1",
            module="Structure Viewer",
            description="Select visualization preset",
            options_map=preset_options_map,
        )
        
        preset_id = presets[int(choice) - 1][0]
        self._apply_visualization_preset(preset_id)

    def _apply_visualization_preset(self, preset_id: str):
        """Apply a visualization preset"""
        # Define preset configurations
        preset_configs = {
            'metal_analysis': {
                'mode': 'original_with_highlights',
                'features': ['metal_sites'],
                'description': 'Metal site analysis preset'
            },
            'mutation_analysis': {
                'mode': 'original_with_highlights', 
                'features': ['pending_mutations', 'applied_mutations'],
                'description': 'Mutation analysis preset'
            },
            'structure_quality': {
                'mode': 'repaired_only' if self.get_from_workspace('repaired_structure') else 'original_only',
                'features': ['repaired_regions'],
                'description': 'Structure quality preset'
            },
            'comprehensive': {
                'mode': 'original_with_highlights',
                'features': ['metal_sites', 'disulfide_bonds', 'mutations', 'protonation_states'],
                'description': 'Comprehensive analysis preset'
            }
        }
        
        config = preset_configs.get(preset_id, preset_configs['comprehensive'])
        
        self.console.print(f"\n[bold]Applying preset: {config['description']}[/bold]")
        
        # Set configuration
        self.visualization_mode = config['mode']
        self.selected_features = config['features']
        
        # Generate and apply representations
        self.suggested_representations = self._generate_feature_representations()
        
        if self.suggested_representations:
            self._generate_workspace_aware_visualization(self.suggested_representations)
        else:
            self.console.print("[yellow]No representations could be generated for this preset.[/yellow]")


# ============================================================================
# PDB DOWNLOAD AND STANDALONE FUNCTIONALITY (Keep existing functionality)
# ============================================================================

def download_pdb(pdb_id):
    """Download PDB file from RCSB"""
    try:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        filename = f"{pdb_id.lower()}.pdb"
        
        console = Console()
        with console.status(f"[bold green]Downloading {pdb_id.upper()}..."):
            urllib.request.urlretrieve(url, filename)
        
        if os.path.exists(filename):
            console.print(f"[green]Downloaded {pdb_id.upper()} to {filename}[/green]")
            return filename
        else:
            console.print(f"[red]Failed to download {pdb_id.upper()}[/red]")
            return None
    except Exception as e:
        console.print(f"[red]Error downloading {pdb_id}: {str(e)}[/red]")
        return None


def main():
    """Main function for standalone mode"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Enhanced VMD Visualization Creator - Standalone Mode"
    )
    parser.add_argument("--pdb", help="Path to local PDB file")
    parser.add_argument("--pdbid", help="PDB ID to download from RCSB")
    args = parser.parse_args()

    # Create visualizer module
    visualizer = EnhancedVMDVisualizationModule()
    visualizer.initialize()  # Initialize in standalone mode

    # Handle command line args
    if args.pdbid:
        # Download PDB file
        pdb_file = download_pdb(args.pdbid)
        if pdb_file:
            visualizer.workspace.set("pdb_id", args.pdbid)
            visualizer.workspace.set("pdb_file", pdb_file)
            # Parse the PDB file
            try:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("protein", pdb_file)
                visualizer.workspace.set("structure", structure)
                visualizer._run_visualization_creator()
                return
            except Exception as e:
                visualizer.console.print(f"[red]Error parsing PDB file: {str(e)}[/red]")

    elif args.pdb:
        if os.path.exists(args.pdb):
            visualizer.workspace.set("pdb_file", args.pdb)
            # Parse the PDB file
            try:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("protein", args.pdb)
                visualizer.workspace.set("structure", structure)
                visualizer._run_visualization_creator()
                return
            except Exception as e:
                visualizer.console.print(f"[red]Error parsing PDB file: {str(e)}[/red]")
        else:
            visualizer.console.print(f"[red]PDB file not found: {args.pdb}[/red]")

    # Interactive mode
    visualizer.console.print(
        "[bold cyan]Enhanced VMD Visualization Module - Standalone Mode[/bold cyan]"
    )
    visualizer.console.print(
        "\nThis module allows you to create sophisticated VMD visualizations with an easy-to-use interface."
    )


if __name__ == "__main__":
    main()
