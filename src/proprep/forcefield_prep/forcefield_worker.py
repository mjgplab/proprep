import glob
import json
import logging
import os
import sys
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.prompts import prompt_with_context, confirm_with_context

import proprep.structure_prep.chem_comp_dict_fetcher as ccd

logger = logging.getLogger(__name__)


class NonStandardResidue:
    """Class to hold information about a non-standard residue"""

    def __init__(self, name: str, chain_id: str, resid: int, category: str = "unknown"):
        """
        Initialize the non-standard residue information

        Args:
            name: The residue name (e.g., MSE)
            chain_id: The chain identifier
            resid: The residue number in the structure
            category: The category (modified_amino_acid, metal_site, unknown)
        """
        self.name = name
        self.chain_id = chain_id
        self.resid = resid
        self.category = category
        self.parent_residue = None
        self.ccd_data = None

    def get_location_str(self) -> str:
        """Get a formatted string of the residue location"""
        return f"{self.chain_id}:{self.resid}"

    def __str__(self) -> str:
        """String representation of the non-standard residue"""
        if self.category == "modified_amino_acid" and self.parent_residue:
            return f"{self.name} ({self.chain_id}:{self.resid}) - Modified {self.parent_residue}"
        else:
            return f"{self.name} ({self.chain_id}:{self.resid}) - {self.category.capitalize()}"


class ForcefieldWorker:
    """Worker class for forcefield parameterization functionality"""

    def __init__(self, workspace: Dict[str, Any], processor=None):
        self.workspace = workspace
        self.processor = processor
        self.console = Console()
        self.ccd_parser = ccd.CCDParser(use_cache=True)
        self.non_standard_residues = []

        self.standard_aa = {
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        }

        self.modified_aa_map = {
            "2AS": "ASP",
            "3AH": "HIS",
            "5HP": "GLU",
            "5OW": "LYS",
            "ACL": "ARG",
            "AGM": "ARG",
            "AIB": "ALA",
            "ALM": "ALA",
            "ALO": "THR",
            "ALY": "LYS",
            "ARM": "ARG",
            "ASA": "ASP",
            "ASB": "ASP",
            "ASK": "ASP",
            "ASL": "ASP",
            "ASQ": "ASP",
            "AYA": "ALA",
            "BCS": "CYS",
            "BHD": "ASP",
            "BMT": "THR",
            "BNN": "ALA",
            "BUC": "CYS",
            "BUG": "LEU",
            "C5C": "CYS",
            "C6C": "CYS",
            "CAS": "CYS",
            "CCS": "CYS",
            "CEA": "CYS",
            "CGU": "GLU",
            "CHG": "ALA",
            "CLE": "LEU",
            "CME": "CYS",
            "CSD": "ALA",
            "CSO": "CYS",
            "CSP": "CYS",
            "CSS": "CYS",
            "CSW": "CYS",
            "CSX": "CYS",
            "CXM": "MET",
            "CY1": "CYS",
            "CY3": "CYS",
            "CYG": "CYS",
            "CYM": "CYS",
            "CYQ": "CYS",
            "DAH": "PHE",
            "DAL": "ALA",
            "DAR": "ARG",
            "DAS": "ASP",
            "DCY": "CYS",
            "DGL": "GLU",
            "DGN": "GLN",
            "DHA": "ALA",
            "DHI": "HIS",
            "DIL": "ILE",
            "DIV": "VAL",
            "DLE": "LEU",
            "DLY": "LYS",
            "DNP": "ALA",
            "DPN": "PHE",
            "DPR": "PRO",
            "DSN": "SER",
            "DSP": "ASP",
            "DTH": "THR",
            "DTR": "TRP",
            "DTY": "TYR",
            "DVA": "VAL",
            "EFC": "CYS",
            "FLA": "ALA",
            "FME": "MET",
            "GGL": "GLU",
            "GL3": "GLY",
            "GLZ": "GLY",
            "GMA": "GLU",
            "GSC": "GLY",
            "HAC": "ALA",
            "HAR": "ARG",
            "HIC": "HIS",
            "HIP": "HIS",
            "HMR": "ARG",
            "HPQ": "PHE",
            "HTR": "TRP",
            "HYP": "PRO",
            "IAS": "ASP",
            "IIL": "ILE",
            "IYR": "TYR",
            "KCX": "LYS",
            "LLP": "LYS",
            "LLY": "LYS",
            "LTR": "TRP",
            "LYM": "LYS",
            "LYZ": "LYS",
            "MAA": "ALA",
            "MEN": "ASN",
            "MHS": "HIS",
            "MIS": "SER",
            "MK8": "LEU",
            "MLE": "LEU",
            "MPQ": "GLY",
            "MSA": "GLY",
            "MSE": "MET",
            "MSO": "MET",
            "NEP": "HIS",
            "NLE": "LEU",
            "NLN": "LEU",
            "NVA": "VAL",
            "ORN": "ARG",
            "PCA": "GLU",
            "PHI": "PHE",
            "PHL": "PHE",
            "PR3": "CYS",
            "PRR": "ALA",
            "PTR": "TYR",
            "SAC": "SER",
            "SAR": "GLY",
            "SEP": "SER",
            "SET": "SER",
            "SHC": "CYS",
            "SOC": "CYS",
            "STY": "TYR",
            "SVA": "SER",
            "TIH": "ALA",
            "TPL": "TRP",
            "TPO": "THR",
            "TPQ": "ALA",
            "TRG": "LYS",
            "TRO": "TRP",
            "TYB": "TYR",
            "TYI": "TYR",
            "TYQ": "TYR",
            "TYS": "TYR",
            "TYY": "TYR",
        }

        self._load_modified_aa_map()

    def analyze_nonstandard_residues(self) -> List[NonStandardResidue]:
        """Analyze the structure for non-standard residues using filter_selections"""
        filter_selections = self.workspace.get("filter_selections")

        if not filter_selections:
            self.console.print(
                "[red]No filter selections found in workspace. Run PDB Filter first.[/red]"
            )
            return []

        self.non_standard_residues = []

        for selection in filter_selections:
            if selection.category in ["modified_protein", "unknown"]:
                # Create non-standard residue object
                residue = NonStandardResidue(
                    name=selection.residue_name,
                    chain_id=selection.chain_id,
                    resid=selection.residue_number,
                    category=(
                        "modified_amino_acid"
                        if selection.category == "modified_protein"
                        else "unknown"
                    ),
                )

                # Check if this is a modified amino acid and set parent
                if residue.name.upper() in self.modified_aa_map:
                    residue.parent_residue = self.modified_aa_map[residue.name.upper()]
                    residue.category = "modified_amino_acid"

                self.non_standard_residues.append(residue)

        # Group by residue name for display
        residue_groups = defaultdict(list)
        for res in self.non_standard_residues:
            residue_groups[res.name].append(res)

        # Display summary
        if self.non_standard_residues:
            table = Table(title="Non-Standard Residues Found")
            table.add_column("Residue", style="cyan")
            table.add_column("Count", style="green")
            table.add_column("Category", style="yellow")
            table.add_column("Parent/Notes", style="blue")
            table.add_column("Locations", style="magenta")

            for res_name, instances in residue_groups.items():
                count = len(instances)
                category = instances[0].category.replace("_", " ").capitalize()
                parent_notes = (
                    instances[0].parent_residue
                    if instances[0].parent_residue
                    else "Unknown"
                )

                locations = [res.get_location_str() for res in instances[:3]]
                if count > 3:
                    locations.append("...")
                location_str = ", ".join(locations)

                table.add_row(
                    res_name, str(count), category, parent_notes, location_str
                )

            self.console.print(table)
            self.console.print(
                f"\nTotal: {len(self.non_standard_residues)} non-standard residues found"
            )
        else:
            self.console.print(
                "[green]No non-standard residues found that require parameterization[/green]"
            )

        return self.non_standard_residues

    def reclassify_residue(self):
        """Interactive function to change the classification of a non-standard residue"""
        if not self.non_standard_residues:
            self.non_standard_residues = self.workspace.get("non_standard_residues", [])

            if not self.non_standard_residues:
                self.console.print(
                    "[yellow]No non-standard residues analyzed yet. Running analysis...[/yellow]"
                )
                self.analyze_nonstandard_residues()

                if not self.non_standard_residues:
                    self.console.print(
                        "[red]No non-standard residues found to reclassify[/red]"
                    )
                    return

        self.console.print("\n[bold]Select a residue to reclassify:[/bold]")
        table = Table(title="Available Non-Standard Residues")
        table.add_column("No.", style="cyan", justify="right")
        table.add_column("Name", style="magenta")
        table.add_column("Location", style="green")
        table.add_column("Current Category", style="yellow")
        table.add_column("Parent/Notes", style="blue")

        displayed_residues = set()
        displayed_list = []

        for res in self.non_standard_residues:
            res_key = (res.name, res.chain_id, res.resid)

            if res_key in displayed_residues:
                continue

            displayed_residues.add(res_key)
            displayed_list.append(res)

            category_display = res.category.replace("_", " ").capitalize()
            location = res.get_location_str()

            notes = ""
            if res.category == "modified_amino_acid" and res.parent_residue:
                notes = f"Parent: {res.parent_residue}"
            elif res.category == "metal_site":
                notes = "Metal coordination site"

            table.add_row(
                str(len(displayed_list)), res.name, location, category_display, notes
            )

        self.console.print(table)

        reclass_options_map = {
            str(i + 1): f"{r.name} ({r.category})" for i, r in enumerate(displayed_list)
        }
        reclass_options_map["q"] = "Cancel"
        choice = prompt_with_context(
            self.processor,
            "Enter residue number to reclassify (or 'q' to cancel)",
            choices=[str(i) for i in range(1, len(displayed_list) + 1)] + ["q"],
            default="q",
            module="Forcefield Worker",
            description="Select non-standard residue to reclassify",
            options_map=reclass_options_map,
        )

        if choice == "q":
            self.console.print("[yellow]Reclassification cancelled[/yellow]")
            return

        selected_idx = int(choice) - 1
        selected_residue = displayed_list[selected_idx]

        self.console.print(
            f"\nSelected: {selected_residue.name} ({selected_residue.get_location_str()})"
        )
        self.console.print(
            f"Current classification: {selected_residue.category.replace('_', ' ').capitalize()}"
        )
        if (
            selected_residue.category == "modified_amino_acid"
            and selected_residue.parent_residue
        ):
            self.console.print(
                f"Current parent residue: {selected_residue.parent_residue}"
            )

        categories = ["modified_amino_acid", "metal_site", "unknown"]

        self.console.print("\nAvailable classifications:")
        for i, cat in enumerate(categories, 1):
            self.console.print(f"  {i}. {cat.replace('_', ' ').capitalize()}")

        cat_options_map = {str(i + 1): c for i, c in enumerate(categories)}
        cat_choice = prompt_with_context(
            self.processor,
            "Select new classification",
            choices=[str(i) for i in range(1, len(categories) + 1)],
            default="1",
            module="Forcefield Worker",
            description="Residue classification category",
            options_map=cat_options_map,
        )

        old_category = selected_residue.category
        new_category = categories[int(cat_choice) - 1]
        selected_residue.category = new_category

        self.console.print(
            f"[green]Successfully reclassified {selected_residue.name} from '{old_category.replace('_', ' ').capitalize()}' to '{new_category.replace('_', ' ').capitalize()}'[/green]"
        )

        if new_category == "modified_amino_acid":
            if selected_residue.name.upper() in self.modified_aa_map:
                mapped_parent = self.modified_aa_map[selected_residue.name.upper()]
                self.console.print(
                    f"[green]Found in mapping table: {selected_residue.name} is derived from {mapped_parent}[/green]"
                )

                use_mapped = confirm_with_context(
                    self.processor,
                    f"Use mapped parent residue ({mapped_parent})?",
                    default=True,
                    module="Forcefield Worker",
                    description=f"Use mapped parent residue {mapped_parent}",
                )

                if use_mapped:
                    selected_residue.parent_residue = mapped_parent
                    self.console.print(
                        f"[green]Set parent residue to {mapped_parent}[/green]"
                    )

                    try_ccd = confirm_with_context(
                        self.processor,
                        "Also check Chemical Component Dictionary?",
                        default=False,
                        module="Forcefield Worker",
                        description="Also check CCD for residue info",
                    )

                    if not try_ccd:
                        if selected_residue.name not in self.modified_aa_map:
                            add_to_map = confirm_with_context(
                                self.processor,
                                f"Add {selected_residue.name} -> {selected_residue.parent_residue} to mapping table for future use?",
                                default=True,
                                module="Forcefield Worker",
                                description=f"Add {selected_residue.name} to residue mapping table",
                            )

                            if add_to_map:
                                self.modified_aa_map[selected_residue.name.upper()] = (
                                    selected_residue.parent_residue
                                )
                                self._save_modified_aa_map()
                                self.console.print(
                                    f"[green]Added {selected_residue.name} -> {selected_residue.parent_residue} to mapping table[/green]"
                                )

                        self._prompt_parameterization(selected_residue)
                        return

            if not selected_residue.parent_residue:
                standard_aa_list = list(self.standard_aa)
                standard_aa_list.sort()

                self.console.print("\n[bold]Select parent residue:[/bold]")
                for i, aa in enumerate(standard_aa_list, 1):
                    self.console.print(f"  {i}. {aa}")

                aa_choice_map = {str(i + 1): aa for i, aa in enumerate(standard_aa_list)}
                aa_choice_map["c"] = "Custom 3-letter code"
                aa_choice = prompt_with_context(
                    self.processor,
                    "Select parent amino acid (or enter 'c' for custom 3-letter code)",
                    choices=[str(i) for i in range(1, len(standard_aa_list) + 1)]
                    + ["c"],
                    default="c",
                    module="Forcefield Worker",
                    description="Select parent amino acid",
                    options_map=aa_choice_map,
                )

                if aa_choice == "c":
                    parent_aa = prompt_with_context(
                        self.processor,
                        "Enter 3-letter code for parent amino acid",
                        default="UNK",
                        module="Forcefield Worker",
                        description="Parent amino acid 3-letter code",
                    ).upper()
                else:
                    parent_aa = standard_aa_list[int(aa_choice) - 1]

                selected_residue.parent_residue = parent_aa
                self.console.print(f"[green]Set parent residue to {parent_aa}[/green]")

                add_to_map = confirm_with_context(
                    self.processor,
                    f"Add {selected_residue.name} -> {parent_aa} to mapping table for future use?",
                    default=True,
                    module="Forcefield Worker",
                    description=f"Add {selected_residue.name} mapping to table",
                )

                if add_to_map:
                    self.modified_aa_map[selected_residue.name.upper()] = parent_aa
                    self._save_modified_aa_map()
                    self.console.print(
                        f"[green]Added {selected_residue.name} -> {parent_aa} to mapping table[/green]"
                    )

            try_ccd = confirm_with_context(
                self.processor,
                "Try to retrieve data from Chemical Component Dictionary?",
                default=True,
                module="Forcefield Worker",
                description="Try CCD for residue data retrieval",
            )

            if try_ccd:
                try:
                    self.console.print(f"Querying CCD for {selected_residue.name}...")
                    ccd_data = self.ccd_parser.get_residue_data(selected_residue.name)
                    selected_residue.ccd_data = ccd_data

                    if (
                        "mon_nstd_parent_comp_id" in ccd_data
                        and ccd_data["mon_nstd_parent_comp_id"]
                    ):
                        ccd_parent = ccd_data["mon_nstd_parent_comp_id"]
                        self.console.print(
                            f"[green]CCD reports parent as {ccd_parent} (you selected {selected_residue.parent_residue})[/green]"
                        )

                        use_ccd_parent = confirm_with_context(
                            self.processor,
                            f"Use CCD parent ({ccd_parent}) instead of current parent ({selected_residue.parent_residue})?",
                            default=True,
                            module="Forcefield Worker",
                            description=f"Use CCD parent ({ccd_parent})",
                        )

                        if use_ccd_parent:
                            old_parent = selected_residue.parent_residue
                            selected_residue.parent_residue = ccd_parent
                            self.console.print(
                                f"[green]Updated parent residue from {old_parent} to {ccd_parent}[/green]"
                            )

                            if (
                                selected_residue.name.upper() in self.modified_aa_map
                                and self.modified_aa_map[selected_residue.name.upper()]
                                != ccd_parent
                            ):
                                update_map = confirm_with_context(
                                    self.processor,
                                    f"Update mapping table to use CCD parent ({ccd_parent}) for {selected_residue.name}?",
                                    default=True,
                                    module="Forcefield Worker",
                                    description=f"Update mapping table with CCD parent for {selected_residue.name}",
                                )

                                if update_map:
                                    self.modified_aa_map[
                                        selected_residue.name.upper()
                                    ] = ccd_parent
                                    self._save_modified_aa_map()
                                    self.console.print(
                                        f"[green]Updated mapping table entry: {selected_residue.name} -> {ccd_parent}[/green]"
                                    )
                    else:
                        self.console.print(
                            "[yellow]No parent residue found in CCD data[/yellow]"
                        )
                except Exception as e:
                    self.console.print(f"[red]Error querying CCD: {str(e)}[/red]")

        self._prompt_parameterization(selected_residue)

    def update_aa_mapping(self):
        """Function to add or update an entry in the modified amino acid mapping table"""
        self.console.print("\n[bold]Update Modified Amino Acid Mapping Table[/bold]")

        self._display_aa_mapping_table()

        mod_aa = prompt_with_context(
            self.processor,
            "Enter the 3-letter code for the modified amino acid",
            default="",
            module="Forcefield Worker",
            description="Modified amino acid 3-letter code",
        ).upper()

        if not mod_aa:
            self.console.print("[yellow]Operation cancelled[/yellow]")
            return

        if mod_aa in self.modified_aa_map:
            current_parent = self.modified_aa_map[mod_aa]
            self.console.print(
                f"[yellow]Note: {mod_aa} is already mapped to {current_parent}[/yellow]"
            )

        standard_aa_list = list(self.standard_aa)
        standard_aa_list.sort()

        self.console.print("\n[bold]Select parent residue:[/bold]")
        for i, aa in enumerate(standard_aa_list, 1):
            self.console.print(f"  {i}. {aa}")

        aa_choice = prompt_with_context(None,
            "Select parent amino acid (or enter 'c' for custom 3-letter code)",
            choices=[str(i) for i in range(1, len(standard_aa_list) + 1)] + ["c"],
            default="c",
        )

        if aa_choice == "c":
            parent_aa = prompt_with_context(
                self.processor,
                "Enter 3-letter code for parent amino acid",
                default="UNK",
                module="Forcefield Worker",
                description="Parent amino acid 3-letter code",
            ).upper()
        else:
            parent_aa = standard_aa_list[int(aa_choice) - 1]

        self.modified_aa_map[mod_aa] = parent_aa
        self._save_modified_aa_map()

        self.console.print(
            f"[green]Successfully added/updated mapping: {mod_aa} -> {parent_aa}[/green]"
        )

        self._display_aa_mapping_table()

    def parameterize_modified_amino_acid(
        self, residue_name: str, residues: List[NonStandardResidue]
    ):
        """Parameterize a modified amino acid residue"""
        if not residues:
            self.console.print("[red]No residues to parameterize[/red]")
            return

        parent_residue = residues[0].parent_residue
        ccd_data = residues[0].ccd_data

        if not parent_residue:
            self.console.print(
                f"[red]Error: No parent residue found for {residue_name}[/red]"
            )
            return

        self.console.print(
            f"\n[bold]Parameterizing Modified Amino Acid: {residue_name}[/bold]"
        )
        self.console.print(f"Parent Standard Residue: {parent_residue}")

        if ccd_data and "name" in ccd_data:
            self.console.print(f"Full Name: {ccd_data['name']}")

            # Show warning about CCD entries potentially being different compounds
            self.console.print(
                "\n[yellow]⚠ Note:[/yellow] The CCD may have an entry for this residue code that differs "
                "from your intended compound."
            )
            self.console.print(
                "[grey50]  For example, 'TYO' in the CCD is a ring-opened oxidative degradation product,[/grey50]"
            )
            self.console.print(
                "[grey50]  not a tyrosyl radical. Verify the Full Name above matches your intended molecule.[/grey50]"
            )

        confirm_parent = confirm_with_context(
            self.processor,
            f"Use {parent_residue} as parent residue?",
            default=True,
            module="Forcefield Worker",
            description=f"Use {parent_residue} as parent residue",
        )

        if not confirm_parent:
            standard_aa_list = list(self.standard_aa)
            standard_aa_list.sort()

            self.console.print("\n[bold]Select parent residue:[/bold]")
            for i, aa in enumerate(standard_aa_list, 1):
                self.console.print(f"  {i}. {aa}")

            aa_map_choice_map = {str(i + 1): aa for i, aa in enumerate(standard_aa_list)}
            aa_map_choice_map["c"] = "Custom 3-letter code"
            aa_choice = prompt_with_context(
                self.processor,
                "Select parent amino acid (or enter 'c' for custom 3-letter code)",
                choices=[str(i) for i in range(1, len(standard_aa_list) + 1)] + ["c"],
                default="c",
                module="Forcefield Worker",
                description="Select parent amino acid (override)",
                options_map=aa_map_choice_map,
            )

            if aa_choice == "c":
                parent_residue = prompt_with_context(
                    self.processor,
                    "Enter 3-letter code for parent amino acid",
                    default=parent_residue,
                    module="Forcefield Worker",
                    description="Custom parent residue 3-letter code",
                ).upper()
            else:
                parent_residue = standard_aa_list[int(aa_choice) - 1]

        self.console.print("[yellow]IMPORTANT WORKFLOW INFORMATION[/yellow]")
        self.console.print(
            f"The parameterization workflow will use the [bold]parent amino acid ({parent_residue})[/bold] structure as a starting point."
        )
        self.console.print(
            f"During the workflow, you will need to modify this structure to create the [bold]modified amino acid ({residue_name})[/bold]."
        )
        from proprep.forcefield_prep.forcefield_parameterizer import (
            _modaa_output_dirname)

        self.console.print(
            f"The output files will be saved in a directory named "
            f"[bold]{_modaa_output_dirname(residue_name)}[/bold]."
        )

        try:
            sys.path.append(os.path.dirname(os.path.abspath(__file__)))
            import modified_amino_acid_parameterizer

            output_dir = os.path.join(os.getcwd(), _modaa_output_dirname(residue_name))
            os.makedirs(output_dir, exist_ok=True)

            original_dir = os.getcwd()
            os.chdir(output_dir)

            self.console.print(f"[green]Output directory: {output_dir}[/green]")
            self.console.print(
                "[green]Launching Modified Amino Acid Parameterizer workflow...[/green]"
            )

            existing_prep = glob.glob(f"{residue_name.lower()}.prep")
            existing_frcmod = glob.glob(f"{residue_name.lower()}.frcmod")

            if existing_prep and existing_frcmod:
                self.console.print(
                    f"[yellow]Parameter files already exist for {residue_name}[/yellow]"
                )
                overwrite = prompt_with_context(
                    self.processor,
                    "Overwrite existing parameter files?",
                    choices=["y", "n"],
                    default="n",
                    module="Forcefield Worker",
                    description="Overwrite existing parameter files",
                    options_map={"y": "Yes, overwrite", "n": "No, keep existing"},
                )

                if overwrite.lower() != "y":
                    self.console.print(
                        f"Using existing parameter files for {residue_name}"
                    )

                    if not self.workspace.has("parameterized_residues"):
                        self.workspace.set("parameterized_residues", {})

                    param_residues = self.workspace.get("parameterized_residues", {})
                    param_residues[residue_name] = {
                        "prep_file": os.path.abspath(existing_prep[0]),
                        "frcmod_file": os.path.abspath(existing_frcmod[0]),
                        "parent_residue": parent_residue,
                        "reused_existing": True,
                        "success": True,
                    }
                    self.workspace.set("parameterized_residues", param_residues)

                    os.chdir(original_dir)
                    return

            workflow_result = modified_amino_acid_parameterizer.run_workflow(
                amino_acid=parent_residue, output_dir=output_dir, interactive=True
            )

            if workflow_result["success"]:
                prep_file = workflow_result["parameter_files"].get("prep_file")
                frcmod_file = workflow_result["parameter_files"].get("frcmod_file")

                if (
                    workflow_result.get("status") == "paused"
                    or workflow_result.get("status") == "pending_calculations"
                ):
                    missing_files = workflow_result.get("missing_files", [])
                    self.console.print(
                        f"[yellow]Parameterization requires Gaussian calculations to complete[/yellow]"
                    )
                    self.console.print(f"The following files need to be created:")
                    for file in missing_files:
                        self.console.print(f"  - {file}")

                    if not self.workspace.has("pending_parameterizations"):
                        self.workspace.set("pending_parameterizations", {})

                    pending_param = self.workspace.get("pending_parameterizations", {})
                    pending_param[residue_name] = {
                        "output_dir": output_dir,
                        "missing_files": missing_files,
                        "status": workflow_result.get("status"),
                        "parent_residue": parent_residue,
                    }
                    self.workspace.set("pending_parameterizations", pending_param)

                    self.console.print(
                        "[yellow]The parameterization process has been paused.[/yellow]"
                    )
                    self.console.print(
                        "You will need to perform the required Gaussian calculations"
                    )
                    self.console.print(
                        "and then restart the parameterization workflow."
                    )

                elif prep_file or frcmod_file:
                    # Store in workspace (usage instructions already shown by library generation)
                    if not self.workspace.has("parameterized_residues"):
                        self.workspace.set("parameterized_residues", {})

                    param_residues = self.workspace.get("parameterized_residues", {})
                    param_residues[residue_name] = {
                        "prep_file": prep_file,
                        "frcmod_file": frcmod_file,
                        "parent_residue": parent_residue,
                        "success": prep_file and frcmod_file,
                        "partial": bool(prep_file) != bool(frcmod_file),
                    }
                    self.workspace.set("parameterized_residues", param_residues)
                else:
                    self.console.print(
                        "[yellow]Workflow completed but no parameter files were generated[/yellow]"
                    )
            else:
                error = workflow_result.get(
                    "error", workflow_result.get("message", "Unknown error")
                )
                self.console.print(
                    f"[red]Error in parameterization workflow: {error}[/red]"
                )

            os.chdir(original_dir)

        except ImportError:
            self.console.print(
                "[yellow]Modified Amino Acid Parameterizer module not found[/yellow]"
            )
            self.console.print(
                f"[italic]This would launch the parameterization for {residue_name} ({parent_residue})[/italic]"
            )
        except Exception as e:
            import traceback

            traceback.print_exc()
            self.console.print(f"[red]Error during parameterization: {str(e)}[/red]")

            if "original_dir" in locals():
                os.chdir(original_dir)

    def display_help(self):
        """Display help information for forcefield parameterization"""
        help_text = """
[bold]Force Field Parameterizer Help[/bold]

This module helps parameterize non-standard residues for MD simulations.

[bold]Available Actions:[/bold]
1. [cyan]Analyze Non-Standard Residues[/cyan] - Scan structure for residues requiring parameterization
2. [cyan]Reclassify Residue[/cyan] - Change the classification of a non-standard residue
3. [cyan]Parameterize Residue[/cyan] - Generate force field parameters for a residue
4. [cyan]Update AA Mapping[/cyan] - Add/update modified amino acid to parent mappings

[bold]Workflow:[/bold]
1. First run "Analyze Non-Standard Residues" to identify problematic residues
2. Use "Reclassify Residue" to set correct categories and parent residues
3. Use "Parameterize Residue" to generate AMBER-compatible parameters
4. The generated .prep and .frcmod files can be used in tleap

[bold]Supported Categories:[/bold]
- [green]Modified Amino Acid[/green]: Variants of standard amino acids (e.g., MSE, SEP)
- [yellow]Metal Site[/yellow]: Metal coordination sites requiring special handling
- [red]Unknown[/red]: Residues that need manual classification

[bold]Notes:[/bold]
- Modified amino acids are automatically mapped to parent residues when possible
- CCD (Chemical Component Dictionary) data is used to validate parent assignments
- Parameterization generates AMBER-compatible .prep and .frcmod files
"""

        panel = Panel(
            help_text, title="Force Field Parameterizer Help", border_style="blue"
        )
        self.console.print(panel)

    def _load_modified_aa_map(self):
        """Load the modified amino acid mapping table from a config file"""
        try:
            config_file = os.path.join(
                os.path.expanduser("~"), ".proprep", "modified_aa_map.json"
            )
            if os.path.exists(config_file):
                with open(config_file, "r") as f:
                    user_map = json.load(f)

                self.modified_aa_map.update(user_map)
                # logger.info(
                #     f"Loaded user-defined mapping table with {len(user_map)} entries"
                # )
        except Exception as e:
            self.console.print(
                f"[yellow]Error loading mapping table: {str(e)}[/yellow]"
            )

    def _save_modified_aa_map(self):
        """Save the modified amino acid mapping table to a config file"""
        try:
            config_dir = os.path.join(os.path.expanduser("~"), ".proprep")
            os.makedirs(config_dir, exist_ok=True)

            config_file = os.path.join(config_dir, "modified_aa_map.json")
            with open(config_file, "w") as f:
                json.dump(self.modified_aa_map, f, indent=2)

            self.console.print(f"[green]Saved mapping table to {config_file}[/green]")
        except Exception as e:
            self.console.print(f"[red]Error saving mapping table: {str(e)}[/red]")

    def _display_aa_mapping_table(self):
        """Display the current modified amino acid mapping table"""
        if not self.modified_aa_map:
            self.console.print(
                "[yellow]The modified amino acid mapping table is empty[/yellow]"
            )
            return

        table = Table(title="Modified Amino Acid Mapping")
        table.add_column("Modified AA", style="cyan")
        table.add_column("Parent AA", style="green")

        sorted_keys = sorted(self.modified_aa_map.keys())

        for mod_aa in sorted_keys:
            parent_aa = self.modified_aa_map[mod_aa]
            table.add_row(mod_aa, parent_aa)

        self.console.print(table)
        self.console.print(f"Total entries: {len(self.modified_aa_map)}")

    def _prompt_parameterization(self, residue):
        """Helper method to prompt for parameterization of a residue"""
        if residue.category == "modified_amino_acid" and residue.parent_residue:
            parameterize = confirm_with_context(
                self.processor,
                f"Would you like to parameterize {residue.name} now?",
                default=False,
                module="Forcefield Worker",
                description=f"Parameterize {residue.name} now",
            )

            if parameterize:
                residues_for_param = [
                    r for r in self.non_standard_residues if r.name == residue.name
                ]
                self.parameterize_modified_amino_acid(residue.name, residues_for_param)
