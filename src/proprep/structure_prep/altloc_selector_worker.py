"""
AltLoc Selector Module

Identifies and processes alternate locations (conformations) in protein structures.
"""

import logging
import os
from collections import defaultdict

from proprep.utils.prompts import prompt_with_context


logger = logging.getLogger(__name__)


class AltLocSelector:
    """Identifies and processes alternate locations in protein structures."""

    def __init__(self, processor=None):
        """Initialize the alternate location selector."""
        self.processor = processor
        self.structure = None
        self.input_file = None
        self.alt_residues = defaultdict(set)
        self.occupancy_stats = {}
        self.atom_counts = {}
        self.selections = {}
        self.counted_atoms = {}  # NEW: Track which atoms we've already counted

    def setup(self, structure=None, input_file=None, selected_chains=None):
        """
        Set up the selector with structure data.

        Args:
            structure: BioPython Structure object
            input_file: Path to PDB file to process
            selected_chains: List of selected chain IDs
        """
        self.structure = structure
        self.input_file = input_file
        self.selected_chains = selected_chains or []

        # Reset data
        self.alt_residues = defaultdict(set)
        self.occupancy_stats = {}
        self.atom_counts = {}
        self.selections = {}
        self.counted_atoms = {}  # NEW: Reset the atom counting tracker

        if self.structure is None and self.input_file is None:
            raise ValueError("Either structure or input_file must be provided")

    def identify_alt_locs(self):
        """
        Identify residues with alternate locations in the structure.

        Returns:
            Dictionary of alternate locations
        """
        logger.info("Identifying alternate locations")

        # Parse the PDB file first
        atoms = self._parse_pdb()

        # Initialize data structures
        self.alt_residues = defaultdict(set)
        self.occupancy_stats = {}
        self.atom_counts = defaultdict(int)
        self.counted_atoms = {}  # NEW: Initialize atom counting tracker

        # First pass: identify residues with alternate locations and their atoms
        residue_atoms = defaultdict(set)  # Track atom names for each residue

        for line in atoms:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            alt_loc = line[16:17].strip()
            chain_id = line[21:22]

            # Skip chains that weren't selected by the user
            if self.selected_chains and chain_id not in self.selected_chains:
                continue

            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            ins_code = line[26:27].strip()
            atom_name = line[12:16].strip()

            # Create a unique residue identifier
            res_id = (chain_id, res_name, res_num, ins_code)

            # Track all atom names in this residue (regardless of alt loc)
            residue_atoms[res_id].add(atom_name)

            # Track alternate locations
            if alt_loc:
                self.alt_residues[res_id].add(alt_loc)

                # Initialize occupancy stats for this residue if needed
                if res_id not in self.occupancy_stats:
                    self.occupancy_stats[res_id] = {
                        "conformations": {},
                        "occupancy_sum_issues": False,
                        "atoms": defaultdict(dict),
                    }

                # Store occupancy data
                occupancy = float(line[54:60].strip())
                self.occupancy_stats[res_id]["atoms"][atom_name][alt_loc] = occupancy

        # Second pass: count atoms correctly for each residue
        for res_id, atom_names in residue_atoms.items():
            # Count unique atom names for this residue
            self.atom_counts[res_id] = len(atom_names)

        # Calculate occupancy statistics
        self._calculate_occupancy_statistics()

        # Return results
        if not self.alt_residues:
            logger.info("No alternate locations found in the structure")
            return {}

        logger.info(f"Found {len(self.alt_residues)} residues with alternate locations")
        self._print_alt_locs_summary()

        return self.alt_residues

    def select_conformations(self, auto_select=False, interactive=True):
        """
        Select which conformation to keep for each residue.

        Args:
            auto_select: Whether to automatically select highest occupancy
            interactive: Whether to prompt for user input

        Returns:
            Dictionary of selected conformations
        """
        if not self.alt_residues:
            logger.warning("No alternate locations to select from")
            return {}

        # Reset selections
        self.selections = {}

        if auto_select:
            # Automatically select highest occupancy conformation for each residue
            for res_id, alt_locs in self.alt_residues.items():
                highest_occ_alt = self._get_highest_occupancy_alt(res_id)
                self.selections[res_id] = highest_occ_alt

            logger.info("Automatically selected highest occupancy conformations")
            self._print_selections_summary()
            return self.selections

        if not interactive:
            # Default to highest occupancy in non-interactive mode
            for res_id, alt_locs in self.alt_residues.items():
                highest_occ_alt = self._get_highest_occupancy_alt(res_id)
                self.selections[res_id] = highest_occ_alt

            logger.info(
                "Selected highest occupancy conformations (non-interactive mode)"
            )
            return self.selections

        # First, ask if user wants to auto-select highest occupancy for all residues
        print(
            f"\nFound {len(self.alt_residues)} residues with alternate conformations."
        )
        auto_choice = (
            prompt_with_context(
                self.processor,
                "Automatically select highest occupancy conformation for all residues? (y/n)",
                module="AltLoc Selector",
                description="Auto-select highest-occupancy conformations for all residues",
            )
            .strip()
            .lower()
        )

        if auto_choice in ["y", "yes"]:
            # Apply highest occupancy selection to all residues
            for res_id, alt_locs in self.alt_residues.items():
                highest_occ_alt = self._get_highest_occupancy_alt(res_id)
                self.selections[res_id] = highest_occ_alt

            print("\nSelected highest occupancy conformation for all residues.")
            self._print_selections_summary()
            return self.selections

        # If user doesn't want auto-selection, proceed with residue-by-residue selection
        print("\nSelect which alternate conformation to keep for each residue:")
        print("-----------------------------------------------------------")

        for res_id, alt_locs in self.alt_residues.items():
            try:
                # Safely handle tuple unpacking in case res_id is malformed
                if isinstance(res_id, (tuple, list)) and len(res_id) >= 4:
                    chain_id, res_name, res_num, ins_code = res_id[0], res_id[1], res_id[2], res_id[3]
                elif isinstance(res_id, str):
                    # Handle string representation like "('A', 'SER', '37', '')"
                    print(f"Warning: Received string residue ID: {res_id}")
                    try:
                        import ast
                        parsed_res_id = ast.literal_eval(res_id)
                        if len(parsed_res_id) >= 4:
                            chain_id, res_name, res_num, ins_code = parsed_res_id[0], parsed_res_id[1], parsed_res_id[2], parsed_res_id[3]
                        else:
                            print(f"Error: Invalid residue ID format: {res_id}")
                            continue
                    except (ValueError, SyntaxError) as e:
                        print(f"Error: Could not parse residue ID: {res_id} - {e}")
                        continue
                else:
                    print(f"Error: Invalid residue ID type or length: {res_id}")
                    continue
                    
                ins_str = f":{ins_code}" if ins_code else ""

                print(
                    f"\nResidue: {res_name} {res_num}{ins_str} (Chain {chain_id}) - {self.atom_counts.get(res_id, 0)} atom(s)"
                )
            except Exception as e:
                print(f"Error processing residue {res_id}: {e}")
                continue

            # Display conformations with their occupancies
            print("Available conformations:")
            has_issues = self.occupancy_stats[res_id]["occupancy_sum_issues"]

            # Calculate total average occupancy
            conf_stats = self.occupancy_stats[res_id]["conformations"]
            total_avg_occ = sum(stats["average"] for stats in conf_stats.values())

            for alt_loc in sorted(alt_locs):
                stats = conf_stats[alt_loc]
                avg_occ = stats["average"]

                # Show normalized occupancy if there are issues
                if has_issues and "normalized_avg" in stats:
                    print(
                        f"  {alt_loc}: Occupancy = {avg_occ:.2f} (Normalized: {stats['normalized_avg']:.2f})"
                    )
                else:
                    print(f"  {alt_loc}: Occupancy = {avg_occ:.2f}")

            # Report on occupancy sum issues
            if abs(total_avg_occ - 1.0) > 0.05:  # More than 5% off from 1.0
                print(
                    f"  Note: Total occupancy ({total_avg_occ:.2f}) doesn't sum to 1.0"
                )
                if has_issues:
                    print("  Some atoms have inconsistent occupancy values")

            # Get the highest occupancy conformation
            highest_occ_alt = self._get_highest_occupancy_alt(res_id)
            print(
                f"Suggestion: Conformation {highest_occ_alt} has the highest occupancy"
            )

            # Get user input
            while True:
                choice = (
                    prompt_with_context(
                        self.processor,
                        f"Select conformation [{'|'.join(sorted(alt_locs))}]",
                        module="AltLoc Selector",
                        description=f"Select conformation for residue {res_id}",
                    )
                    .strip()
                    .upper()
                )
                if choice in alt_locs:
                    self.selections[res_id] = choice
                    break
                else:
                    print(
                        f"Invalid choice. Please select from: {', '.join(sorted(alt_locs))}"
                    )

        print("\nConformation selections have been recorded.")
        self._print_selections_summary()
        return self.selections

    def process_and_save_structure(
        self, output_file, method=1, save_selected_chains_only=False
    ):
        """
        Process and save the structure with selected conformations.

        Args:
            output_file: Path to save the processed structure
            method: Processing method:
                    1 = Remove alternate conformations (keep only selected)
                    2 = Set occupancy to 1.0 for selected and 0.0 for others
                    3 = Normalize occupancies
            save_selected_chains_only: Whether to save only the selected chains

        Returns:
            Path to the processed structure file
        """
        if not self.alt_residues:
            logger.warning("No alternate locations to process")
            return None

        if not self.selections:
            logger.warning("No selections made. Select conformations first.")
            return None

        # Generate output filename if not provided
        if not output_file:
            base = os.path.splitext(os.path.basename(self.input_file))[0]
            chains_suffix = "_selected_chains" if save_selected_chains_only else ""
            output_file = f"{base}_altloc_processed{chains_suffix}.pdb"

        # Parse PDB file if not done already
        atoms = self._parse_pdb()

        # Process the structure according to the selected method
        if method == 1:
            self._create_new_pdb(atoms, output_file, save_selected_chains_only)
        elif method == 2:
            self._create_occupancy_pdb(atoms, output_file, save_selected_chains_only)
        elif method == 3:
            self._create_normalized_pdb(atoms, output_file, save_selected_chains_only)
        else:
            logger.error(f"Invalid method: {method}")
            return None

        return output_file

    def process_structure(self, output_file=None, method=1):
        """
        Create a new PDB file with the selected conformations.

        Args:
            output_file: Path to save the processed structure
            method: Processing method:
                    1 = Remove alternate conformations (keep only selected)
                    3 = Normalize occupancies (keep relative proportions but ensure sum = 1.0)

        Returns:
            Path to the processed structure file
        """
        # Ask if user wants to save only selected chains
        save_selected_only = False
        if (
            len(self.selected_chains) > 0
            and self.selected_chains != self._get_all_chains()
        ):
            print("\nChain output options:")
            print("1. Save all chains (processed + unmodified)")
            print("2. Save only selected chains that were processed")

            while True:
                choice = prompt_with_context(
                    self.processor, "Select chain output option [1/2]",
                    module="AltLoc Selector",
                    description="Chain output option",
                    options_map={
                        "1": "Save all chains (processed + unmodified)",
                        "2": "Save only selected chains that were processed",
                    },
                ).strip()
                if choice == "1":
                    save_selected_only = False
                    break
                elif choice == "2":
                    save_selected_only = True
                    break
                else:
                    print("Invalid choice. Please enter 1 or 2.")

        # We no longer support option 2, so ensure method is either 1 or 3
        if method != 1 and method != 3:
            logger.warning(
                f"Unsupported method {method}. Using method 1 (remove alternate conformations)."
            )
            method = 1

        return self.process_and_save_structure(output_file, method, save_selected_only)

    def get_output_options(self):
        """
        Ask the user for the desired output processing method.

        Returns:
            int: Selected method (1 or 3)
        """
        print("\nOutput options:")
        print("1. Remove alternate conformations (keep only selected)")
        print(
            "2. Normalize occupancies (keep relative proportions but ensure sum = 1.0)"
        )

        while True:
            choice = prompt_with_context(
                self.processor, "Select output method [1/2]",
                module="AltLoc Selector",
                description="Output processing method",
                options_map={
                    "1": "Remove alternate conformations (keep only selected)",
                    "2": "Normalize occupancies (keep relative proportions, sum = 1.0)",
                },
            ).strip()
            if choice == "1":
                return 1
            elif choice == "2":
                return 3  # Map choice 2 to method 3 (normalization)
            else:
                print("Invalid choice. Please enter 1 or 2.")

    def _parse_pdb(self):
        """
        Parse a PDB file and extract ATOM/HETATM records.

        Returns:
            List of PDB lines
        """
        if self.input_file is None:
            raise ValueError("No input file specified")

        atoms = []
        try:
            with open(self.input_file, "r") as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atoms.append(line)
                    elif line.startswith("ANISOU"):
                        # Skip ANISOU records as they're tied to atoms and will be filtered later
                        continue
                    else:
                        # Keep other records (like HEADER, REMARK, etc.)
                        atoms.append(line)
        except Exception as e:
            logger.error(f"Error parsing PDB file: {str(e)}")
            raise

        return atoms

    def _calculate_occupancy_statistics(self):
        """Calculate occupancy statistics for each residue and conformation."""
        for res_id in self.alt_residues:
            if res_id not in self.occupancy_stats:
                continue

            # Get atoms and their occupancies
            atoms = self.occupancy_stats[res_id]["atoms"]
            alt_locs = self.alt_residues[res_id]

            # Initialize conformations data
            self.occupancy_stats[res_id]["conformations"] = {
                alt_loc: {"sum": 0.0, "count": 0, "average": 0.0}
                for alt_loc in alt_locs
            }

            # Check for occupancy sum issues
            occupancy_sum_issues = False
            for atom_name, occ_dict in atoms.items():
                occ_sum = sum(occ_dict.values())
                if abs(occ_sum - 1.0) > 0.05:  # More than 5% off from 1.0
                    occupancy_sum_issues = True
                    break

            self.occupancy_stats[res_id]["occupancy_sum_issues"] = occupancy_sum_issues

            # Calculate statistics for each conformation
            for atom_name, occ_dict in atoms.items():
                for alt_loc, occ in occ_dict.items():
                    if alt_loc in self.occupancy_stats[res_id]["conformations"]:
                        conf_stats = self.occupancy_stats[res_id]["conformations"][
                            alt_loc
                        ]
                        conf_stats["sum"] += occ
                        conf_stats["count"] += 1

            # Calculate average occupancy for each conformation
            for alt_loc, stats in self.occupancy_stats[res_id]["conformations"].items():
                if stats["count"] > 0:
                    stats["average"] = stats["sum"] / stats["count"]

            # Calculate normalized occupancies if there are issues
            if occupancy_sum_issues:
                for atom_name, occ_dict in atoms.items():
                    occ_sum = sum(occ_dict.values())
                    if occ_sum > 0:
                        for alt_loc in occ_dict:
                            if alt_loc in self.occupancy_stats[res_id]["conformations"]:
                                norm_value = occ_dict[alt_loc] / occ_sum
                                self.occupancy_stats[res_id]["conformations"][alt_loc][
                                    "normalized_avg"
                                ] = norm_value

    def _get_all_chains(self):
        """Get all chain IDs in the structure"""
        chains = set()
        for line in self._parse_pdb():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21:22]
                chains.add(chain_id)
        return list(chains)

    def _get_highest_occupancy_alt(self, res_id):
        """
        Get the alternate location with the highest occupancy for a residue.

        Args:
            res_id: Residue identifier tuple

        Returns:
            String: Alternate location identifier
        """
        if res_id not in self.occupancy_stats:
            return None

        conf_stats = self.occupancy_stats[res_id]["conformations"]
        if not conf_stats:
            return None

        return max(conf_stats.items(), key=lambda x: x[1]["average"])[0]

    def _print_alt_locs_summary(self):
        """Print a summary of alternate locations."""
        if not self.alt_residues:
            logger.info("No alternate conformations found in this PDB file.")
            return

        logger.info(
            f"Found {len(self.alt_residues)} residues with alternate conformations."
        )

        print("\nSummary of alternate conformations:")
        print("----------------------------------")
        for res_id, alt_locs in self.alt_residues.items():
            chain_id, res_name, res_num, ins_code = res_id
            ins_str = f":{ins_code}" if ins_code else ""

            print(
                f"Residue: {res_name} {res_num}{ins_str} (Chain {chain_id}) - {self.atom_counts[res_id]} atom(s)"
            )

            # Calculate total occupancy
            conf_stats = self.occupancy_stats[res_id]["conformations"]
            total_occ = sum(stats["average"] for stats in conf_stats.values())

            # Display each conformation with its occupancy
            for alt_loc in sorted(alt_locs):
                stats = conf_stats[alt_loc]
                avg_occ = stats["average"]

                # Show percentage of total
                pct_of_total = (avg_occ / total_occ * 100) if total_occ > 0 else 0

                print(
                    f"  Conformation {alt_loc}: Occupancy = {avg_occ:.2f} ({pct_of_total:.1f}% of total)"
                )

            # Report any occupancy issues
            if abs(total_occ - 1.0) > 0.05:  # More than 5% off from 1.0
                print(f"  Note: Total occupancy ({total_occ:.2f}) doesn't sum to 1.0")
                if self.occupancy_stats[res_id]["occupancy_sum_issues"]:
                    print(
                        "  Some atoms in this residue have inconsistent occupancy values"
                    )

    def _print_selections_summary(self):
        """Print a summary of selected conformations."""
        if not self.selections:
            logger.info("No selections made")
            return

        print("\nSummary of selections:")
        print("---------------------")
        for res_id, alt_loc in self.selections.items():
            chain_id, res_name, res_num, ins_code = res_id
            ins_str = f":{ins_code}" if ins_code else ""

            # Get occupancy for the selected conformation
            avg_occ = None
            if res_id in self.occupancy_stats:
                conf_stats = self.occupancy_stats[res_id]["conformations"]
                if alt_loc in conf_stats:
                    avg_occ = conf_stats[alt_loc]["average"]

            if avg_occ is not None:
                print(
                    f"Residue: {res_name} {res_num}{ins_str} (Chain {chain_id}) - Selected: {alt_loc} (Occupancy: {avg_occ:.2f})"
                )
            else:
                print(
                    f"Residue: {res_name} {res_num}{ins_str} (Chain {chain_id}) - Selected: {alt_loc}"
                )

    def _create_new_pdb(self, atoms, output_file, save_selected_chains_only=False):
        """
        Create a new PDB file with only the selected conformations.

        Args:
            atoms: List of PDB lines
            output_file: Path to save the processed structure
            save_selected_chains_only: Whether to save only the selected chains
        """
        with open(output_file, "w") as out:
            for line in atoms:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    # Write non-atom lines directly
                    out.write(line)
                    continue

                alt_loc = line[16:17].strip()
                chain_id = line[21:22]

                # Check if this chain should be included
                if save_selected_chains_only and chain_id not in self.selected_chains:
                    # Skip chains not selected when saving only selected chains
                    continue

                # If chain wasn't selected, write the line unchanged
                if chain_id not in self.selected_chains:
                    out.write(line)
                    continue

                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                ins_code = line[26:27].strip()

                res_id = (chain_id, res_name, res_num, ins_code)

                # Write the atom if:
                # 1. It has no alternate location, or
                # 2. It has the selected alternate location for its residue
                if not alt_loc or (
                    res_id in self.selections and self.selections[res_id] == alt_loc
                ):
                    # Remove the alternate location designation
                    new_line = line[:16] + " " + line[17:]
                    out.write(new_line)
                # Skip atoms with non-selected alternate locations in selected chains

    def _create_normalized_pdb(
        self, atoms, output_file, save_selected_chains_only=False
    ):
        """
        Create a new PDB file with normalized occupancies for alternate conformations.

        Args:
            atoms: List of PDB lines
            output_file: Path to save the processed structure
            save_selected_chains_only: Whether to save only the selected chains
        """
        # First, collect all atoms by residue and altloc to calculate normalization factors
        residue_altloc_atoms = defaultdict(lambda: defaultdict(list))

        for i, line in enumerate(atoms):
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            alt_loc = line[16:17].strip()
            if not alt_loc:
                continue

            chain_id = line[21:22]

            # Skip chains that weren't selected
            if chain_id not in self.selected_chains:
                continue

            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            ins_code = line[26:27].strip()

            # Create a unique residue identifier
            res_id = (chain_id, res_name, res_num, ins_code)

            # Store the atom index and occupancy
            occupancy = float(line[54:60].strip())
            residue_altloc_atoms[res_id][alt_loc].append((i, occupancy))

        # Calculate normalization factors for each residue
        normalization_factors = {}
        for res_id, altloc_atoms in residue_altloc_atoms.items():
            # Get all unique atom positions in this residue
            atom_positions = defaultdict(dict)

            for alt_loc, atom_list in altloc_atoms.items():
                for idx, occ in atom_list:
                    atom_name = atoms[idx][12:16].strip()
                    atom_positions[atom_name][alt_loc] = occ

            # For each atom position, calculate the sum of occupancies
            sum_by_position = {}
            for atom_name, alt_dict in atom_positions.items():
                sum_by_position[atom_name] = sum(alt_dict.values())

            # Use the average normalization factor for the residue
            if sum_by_position:
                avg_sum = sum(sum_by_position.values()) / len(sum_by_position)
                normalization_factors[res_id] = 1.0 / avg_sum if avg_sum > 0 else 1.0
            else:
                normalization_factors[res_id] = 1.0

        # Create the new PDB with normalized occupancies
        with open(output_file, "w") as out:
            for i, line in enumerate(atoms):
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    # Write non-atom lines directly
                    out.write(line)
                    continue

                alt_loc = line[16:17].strip()
                chain_id = line[21:22]

                # Check if this chain should be included
                if save_selected_chains_only and chain_id not in self.selected_chains:
                    # Skip chains not selected when saving only selected chains
                    continue

                # If chain wasn't selected, write the line unchanged
                if chain_id not in self.selected_chains:
                    out.write(line)
                    continue

                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                ins_code = line[26:27].strip()

                res_id = (chain_id, res_name, res_num, ins_code)

                if alt_loc and res_id in normalization_factors:
                    # Get current occupancy
                    occupancy = float(line[54:60].strip())

                    # Calculate normalized occupancy
                    norm_factor = normalization_factors[res_id]
                    new_occupancy = occupancy * norm_factor

                    # Format to ensure exactly 2 decimal places in a 6-character field
                    occupancy_str = f"{new_occupancy:6.2f}"

                    # Update the occupancy field (columns 55-60)
                    new_line = line[:54] + occupancy_str + line[60:]
                    out.write(new_line)
                else:
                    # Keep as is if no alternate location or not in our list
                    out.write(line)

    # For backward compatibility with method 2 (keeping it as a no-op)
    def _create_occupancy_pdb(
        self, atoms, output_file, save_selected_chains_only=False
    ):
        """
        Method 2 is deprecated, redirect to method 1.
        """
        logger.warning("Occupancy method (2) is deprecated. Using method 1 instead.")
        return self._create_new_pdb(atoms, output_file, save_selected_chains_only)
