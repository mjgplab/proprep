"""
Structure Selector - Standardized Workspace Structure Access

Provides centralized, configurable access to structures stored in the workspace.
Supports both automatic priority-based selection and interactive user selection.

This module standardizes structure retrieval across ProPrep, making it easy to:
- Add new structure types to the workspace
- Change priority orderings
- Allow user choice when needed
- Maintain backward compatibility

Author: ProPrep Developer
Date: 2025-11-08
"""

import os
import logging
from typing import Optional, List, Dict, Any, Tuple, Union
from dataclasses import dataclass, field
from rich.console import Console
from proprep.utils.prompts import prompt_with_context
from rich.table import Table

logger = logging.getLogger(__name__)


# ============================================================================
# Source Filters - Categories of structure types
# ============================================================================

SOURCE_FILTERS = {
    "preprocessing": [
        "preprocessing_protein_input",  # Extracted protein for protonation/tleap
        "preprocessing_ligand_input",   # Extracted ligand for small_mol_parameterizer
    ],
    "experimental": [
        "rcsb_pdb_file",
        "rcsb_pdb_files",  # Batch-downloaded RCSB structures (list-type key)
        "local_pdb_file",
    ],
    "predicted": [
        "alphafold_pdb_file",
        "alphafill_pdb_file",
        "alphafold_homolog_pdb_file",
    ],
    "processed": [
        "protonation_pdb_file",  # After key rename from structure_with_prot_resnames
        "structure_with_prot_resnames",  # Legacy support during migration
        "transformed_pdb_file",
        "repaired_pdb_file",
        "filtered_pdb_file",
        "hstripped_pdb_file",  # Hydrogen atoms removed for MD preparation
        "topology_extracted_pdb",  # Extracted from prmtop/rst7 via cpptraj
    ],
    "aligned": [
        "aligned_target_pdb_file",
        "aligned_ref_pdb_file",
    ],
}


@dataclass
class StructureType:
    """
    Definition of a structure type that can be stored in workspace.

    Attributes:
        workspace_key: Key used in workspace dict (e.g., "filtered_pdb_file")
        display_name: Human-readable name (e.g., "Filtered")
        priority: Priority level (1=highest, higher numbers=lower priority)
        description: Optional detailed description
    """
    workspace_key: str
    display_name: str
    priority: int
    description: str = ""

    def __post_init__(self):
        """Register this structure type in the global registry"""
        StructureRegistry.register(self)


@dataclass
class StructureInfo:
    """Information about a structure found in workspace"""
    structure_type: StructureType
    file_path: str
    exists: bool = field(init=False)
    file_size: Optional[int] = field(init=False, default=None)

    def __post_init__(self):
        """Validate file existence and get size"""
        self.exists = os.path.exists(self.file_path) if self.file_path else False
        if self.exists:
            try:
                self.file_size = os.path.getsize(self.file_path)
            except:
                self.file_size = None


class StructureRegistry:
    """
    Global registry of structure types.

    Uses a class-level dict to maintain registered structure types.
    Allows dynamic registration at runtime.
    """
    _registry: Dict[str, StructureType] = {}
    _initialized = False

    @classmethod
    def register(cls, structure_type: StructureType):
        """
        Register a structure type.

        Args:
            structure_type: StructureType instance to register
        """
        cls._registry[structure_type.workspace_key] = structure_type
        logger.debug(f"Registered structure type: {structure_type.workspace_key} (display: {structure_type.display_name}, priority: {structure_type.priority})")

    @classmethod
    def get(cls, workspace_key: str) -> Optional[StructureType]:
        """Get structure type by workspace key"""
        cls._ensure_initialized()
        return cls._registry.get(workspace_key)

    @classmethod
    def get_all(cls) -> List[StructureType]:
        """Get all registered structure types sorted by priority"""
        cls._ensure_initialized()
        return sorted(cls._registry.values(), key=lambda x: x.priority)

    @classmethod
    def get_by_priority_range(cls, min_priority: int, max_priority: int) -> List[StructureType]:
        """Get structure types within priority range"""
        cls._ensure_initialized()
        return [
            st for st in cls.get_all()
            if min_priority <= st.priority <= max_priority
        ]

    @classmethod
    def clear(cls):
        """Clear all registered types (mainly for testing)"""
        cls._registry.clear()
        cls._initialized = False

    @classmethod
    def _ensure_initialized(cls):
        """Ensure default structure types are registered"""
        if not cls._initialized:
            cls._register_defaults()
            cls._initialized = True

    @classmethod
    def _register_defaults(cls):
        """Register default ProPrep structure types"""
        # These are created but auto-register via __post_init__
        _ = [
            # Final preprocessed structure - highest priority (the output of tLEaP + metal reinsertion)
            StructureType(
                workspace_key="prepared_pdb",
                display_name="Prepared (tLEaP Output)",
                priority=-1,
                description="Final structure after tLEaP processing with hydrogens and reinserted metals"
            ),
            # Preprocessing inputs - high priority but temporary (used DURING preprocessing)
            # These are temporary structures used during preprocessing workflows
            StructureType(
                workspace_key="preprocessing_protein_input",
                display_name="Preprocessing: Protein",
                priority=0,
                description="Extracted protein residues for protonation/tleap preprocessing"
            ),
            StructureType(
                workspace_key="preprocessing_ligand_input",
                display_name="Preprocessing: Ligand",
                priority=0,
                description="Extracted ligand for small molecule parameterization preprocessing"
            ),
            # Primary key for protonation-updated structures (new standardized name)
            StructureType(
                workspace_key="protonation_pdb_file",
                display_name="Protonation-Updated",
                priority=1,
                description="After H++ or protonation state assignment with renamed residues"
            ),
            # Legacy alias - same priority so either key works during migration
            # TODO: Remove this after all modules are migrated to use protonation_pdb_file
            StructureType(
                workspace_key="structure_with_prot_resnames",
                display_name="Protonation-Updated (Legacy)",
                priority=1,
                description="Legacy key - use protonation_pdb_file instead"
            ),
            StructureType(
                workspace_key="transformed_pdb_file",
                display_name="Redox-Transformed",
                priority=2,
                description="After redox site transformations (Metallo integration)"
            ),
            StructureType(
                workspace_key="repaired_pdb_file",
                display_name="Completeness-Repaired",
                priority=3,
                description="After structure completeness fixes (MODELLER)"
            ),
            StructureType(
                workspace_key="filtered_pdb_file",
                display_name="Filtered",
                priority=4,
                description="After PDB filtering (chain/residue selection)"
            ),
            StructureType(
                workspace_key="oriented_pdb_file",
                display_name="Oriented",
                priority=4,
                description="Structure oriented along Cartesian axes"
            ),
            StructureType(
                workspace_key="topology_extracted_pdb",
                display_name="Topology-Extracted",
                priority=5,
                description="PDB extracted from prmtop/rst7 topology via cpptraj"
            ),
            StructureType(
                workspace_key="aligned_target_pdb_file",
                display_name="Aligned Target",
                priority=5,
                description="Target structure after structural alignment (aligned to reference)"
            ),
            StructureType(
                workspace_key="aligned_ref_pdb_file",
                display_name="Aligned Reference",
                priority=6,
                description="Reference structure from structural alignment"
            ),
            StructureType(
                workspace_key="hstripped_pdb_file",
                display_name="H-Stripped",
                priority=7,
                description="Hydrogen atoms removed for MD preparation"
            ),
            StructureType(
                workspace_key="biological_assembly_pdb_file",
                display_name="Biological Assembly",
                priority=8,
                description="Biological assembly generated from asymmetric unit transformations"
            ),
            StructureType(
                workspace_key="homology_model_pdb_file",
                display_name="Homology Model",
                priority=8,
                description="MODELLER-built homology model derived from a BLAST hit"
            ),
            StructureType(
                workspace_key="alphafold_homolog_pdb_file",
                display_name="AlphaFold Homolog",
                priority=9,
                description="Selected homolog from BLAST search (via AlphaFold Database)"
            ),
            StructureType(
                workspace_key="rcsb_pdb_file",
                display_name="RCSB PDB",
                priority=10,
                description="Downloaded from RCSB Protein Data Bank (experimental data)"
            ),
            # List-type key for multiple RCSB downloads (expanded in get_available_structures)
            StructureType(
                workspace_key="rcsb_pdb_files",
                display_name="RCSB PDB (Multiple)",
                priority=10,
                description="Multiple downloaded RCSB PDB structures (batch download)"
            ),
            StructureType(
                workspace_key="local_pdb_file",
                display_name="Local PDB",
                priority=11,
                description="Loaded from local PDB file on disk"
            ),
            StructureType(
                workspace_key="alphafill_pdb_file",
                display_name="AlphaFill",
                priority=12,
                description="AlphaFold structure enriched with transplanted ligands/cofactors"
            ),
            StructureType(
                workspace_key="alphafold_pdb_file",
                display_name="AlphaFold",
                priority=13,
                description="AlphaFold predicted structure from AlphaFold Database"
            ),
            StructureType(
                workspace_key="pdb_file",
                display_name="Legacy/Loaded",
                priority=14,
                description="Backward compatibility key"
            )
        ]


class StructureSelector:
    """
    Centralized structure selection from workspace.

    Usage:
        # Automatic (priority-based) selection:
        selector = StructureSelector(workspace, console)
        pdb_file = selector.get_structure()

        # Interactive (user choice):
        selector = StructureSelector(workspace, console)
        pdb_file = selector.get_structure(interactive=True)

        # Get all available structures:
        available = selector.get_available_structures()

        # Custom priority (exclude certain types):
        pdb_file = selector.get_structure(
            exclude_keys=["pdb_file"]
        )
    """

    def __init__(self, workspace, console: Optional[Console] = None, processor=None):
        """
        Initialize structure selector.

        Args:
            workspace: ProPrep workspace object
            console: Optional Rich console for output
            processor: Optional processor for session recording context
        """
        self.workspace = workspace
        self.console = console or Console()
        self.processor = processor
        self._cache = {}  # Cache structure info to avoid repeated file checks

    def get_structure(
        self,
        interactive: bool = False,
        exclude_keys: Optional[List[str]] = None,
        include_legacy: bool = True,
        silent: bool = False,
        min_priority: Optional[int] = None,
        max_priority: Optional[int] = None,
        return_key: bool = False,
        # New parameters for flexible selection
        source_filter: Optional[str] = None,
        priority_override: Optional[List[str]] = None,
        require_file_path: bool = True,
        prefer_file_over_object: bool = True,
    ) -> Optional[Union[str, Tuple[str, str]]]:
        """
        Get structure file path from workspace.

        Args:
            interactive: If True, prompt user to choose from available structures
            exclude_keys: List of workspace keys to exclude from selection
            include_legacy: If False, exclude "pdb_file" key (backward compat)
            silent: If True, suppress console output
            min_priority: Minimum priority level to consider (1=highest)
            max_priority: Maximum priority level to consider
            return_key: If True, return tuple of (file_path, workspace_key)
            source_filter: Filter by source category. Options:
                - "experimental": Only RCSB and local PDB files
                - "predicted": Only AlphaFold, AlphaFill, homologs
                - "processed": Only filtered, repaired, transformed, protonation
                - "aligned": Only aligned structures
                - None: No filtering (default)
            priority_override: Custom priority order as list of workspace keys.
                If provided, ONLY these keys are considered, in this order.
                Overrides StructureRegistry priority entirely.
            require_file_path: If True (default), only return if value is a
                valid file path that exists on disk.
            prefer_file_over_object: If True (default), for keys that might have
                both *_pdb_file and *_structure variants, prefer file paths.
                Currently not implemented - reserved for future use.

        Returns:
            Path to selected PDB file, or tuple of (path, key) if return_key=True,
            or None if no valid file found
        """
        available = self.get_available_structures(
            exclude_keys=exclude_keys,
            include_legacy=include_legacy,
            min_priority=min_priority,
            max_priority=max_priority,
            source_filter=source_filter,
            priority_override=priority_override,
            require_file_path=require_file_path,
        )

        if not available:
            if not silent:
                self.console.print("[red]No valid PDB file found in workspace[/red]")
                self.console.print("[yellow]Please load a PDB file first[/yellow]")
                logger.error("No valid PDB file found in workspace")
            return None

        if interactive:
            return self._interactive_selection(available, silent=silent, return_key=return_key)
        else:
            return self._priority_selection(available, silent=silent, return_key=return_key)

    def get_structure_object(
        self,
        interactive: bool = False,
        exclude_keys: Optional[List[str]] = None,
        silent: bool = False,
        return_key: bool = False,
        source_filter: Optional[str] = None,
        priority_override: Optional[List[str]] = None,
    ) -> Optional[Union[Any, Tuple[Any, str]]]:
        """
        Get BioPython Structure object from workspace.

        Similar to get_structure() but returns BioPython Structure objects
        instead of file paths. Automatically maps *_pdb_file keys to their
        corresponding *_structure keys.

        Args:
            interactive: If True, prompt user to choose from available structures
            exclude_keys: List of workspace keys to exclude from selection
            silent: If True, suppress console output
            return_key: If True, return tuple of (structure, workspace_key)
            source_filter: Filter by source category (experimental, predicted, etc.)
            priority_override: Custom priority order as list of workspace keys

        Returns:
            BioPython Structure object, or tuple of (structure, key) if return_key=True,
            or None if no valid structure found
        """
        # Build list of structure object keys based on file key priorities
        available = self.get_available_structures(
            exclude_keys=exclude_keys,
            source_filter=source_filter,
            priority_override=priority_override,
            require_file_path=False,  # We'll check structure objects separately
        )

        # Map file keys to structure object keys and check availability
        structure_candidates = []
        for info in available:
            file_key = info.structure_type.workspace_key
            # Convert file key to structure object key
            if file_key.endswith("_pdb_file"):
                structure_key = file_key.replace("_pdb_file", "_structure")
            elif file_key.endswith("_file"):
                structure_key = file_key.replace("_file", "_structure")
            else:
                structure_key = file_key + "_structure"

            # Check if structure object exists in workspace
            structure_obj = self.workspace.get(structure_key)
            if structure_obj is not None:
                # Enforce the *_structure contract: these keys must hold a
                # parsed Structure object, not a path. Some writers store a
                # path string under a *_structure key for other consumers
                # (e.g. redox_transformation_manager sets transformed_structure
                # to str(output_file) for the modAA parameterizer). Returning a
                # string here would let callers iterate it and crash on `.id`;
                # skip it so they fall through to the file-path selector, which
                # is the correct way to consume a path.
                if isinstance(structure_obj, (str, os.PathLike)):
                    logger.debug(
                        "Skipping workspace key %s: holds a path, not a "
                        "Structure object", structure_key
                    )
                    continue
                structure_candidates.append((structure_obj, structure_key, info))

        if not structure_candidates:
            if not silent:
                self.console.print("[red]No valid Structure objects found in workspace[/red]")
                self.console.print("[yellow]Try loading a PDB file first[/yellow]")
            return None

        # For single candidate or non-interactive, use first (highest priority)
        if not interactive or len(structure_candidates) == 1:
            structure_obj, structure_key, info = structure_candidates[0]
            if not silent:
                self.console.print(
                    f"[green]Using {info.structure_type.display_name} structure object[/green]"
                )
            if return_key:
                return structure_obj, structure_key
            return structure_obj

        # Interactive selection - show table of available structures
        table = Table(title="Available Structure Objects in Workspace")
        table.add_column("Option", style="cyan", justify="center")
        table.add_column("Type", style="green")
        table.add_column("Workspace Key", style="yellow")

        for idx, (_, key, info) in enumerate(structure_candidates, 1):
            table.add_row(
                str(idx),
                info.structure_type.display_name,
                key
            )

        self.console.print(table)

        # Prompt for selection
        from proprep.utils.prompts import prompt_with_context

        session_options_map = {
            str(idx): f"{info.structure_type.display_name} ({key})"
            for idx, (_, key, info) in enumerate(structure_candidates, 1)
        }

        if self.processor:
            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select structure to use",
                choices=[str(i) for i in range(1, len(structure_candidates) + 1)],
                default="1",
                module="Structure Selector",
                description="Select structure object from workspace",
                options_map=session_options_map,
            )
        else:
            choice = prompt_with_context(
                None,
                "Select structure",
                choices=[str(i) for i in range(1, len(structure_candidates) + 1)],
                default="1",
            )

        try:
            idx = int(choice) - 1
            if 0 <= idx < len(structure_candidates):
                structure_obj, structure_key, info = structure_candidates[idx]
                if return_key:
                    return structure_obj, structure_key
                return structure_obj
        except ValueError:
            pass

        # Invalid selection, use first candidate
        structure_obj, structure_key, _ = structure_candidates[0]
        if return_key:
            return structure_obj, structure_key
        return structure_obj

    def _expand_list_key(
        self,
        structure_type: StructureType,
        file_paths: List[str],
        require_file_path: bool = True
    ) -> List[StructureInfo]:
        """
        Expand a list-type workspace key into individual StructureInfo objects.

        For example, rcsb_pdb_files = ["1ABC.pdb", "2XYZ.pdb"] becomes two
        StructureInfo entries with display names "RCSB PDB (1ABC)" and "RCSB PDB (2XYZ)".

        Args:
            structure_type: The StructureType for this list key
            file_paths: List of file paths from the workspace
            require_file_path: If True, only include paths that exist

        Returns:
            List of StructureInfo objects, one per file in the list
        """
        expanded = []

        # Get download info for labeling (if available)
        # e.g., rcsb_pdb_files -> rcsb_download_info
        info_key = structure_type.workspace_key.replace("_files", "_download_info")
        download_info = self.workspace.get(info_key, [])

        for idx, file_path in enumerate(file_paths):
            if not isinstance(file_path, str):
                continue

            # Extract PDB ID for display name
            pdb_id = ""
            if idx < len(download_info) and isinstance(download_info[idx], dict):
                pdb_id = download_info[idx].get("pdb_id", "")

            # Extract from filename if not in download_info
            if not pdb_id:
                basename = os.path.basename(file_path)
                pdb_id = basename.split(".")[0].upper()

            # Create display name with PDB ID
            # "RCSB PDB (Multiple)" -> "RCSB PDB (1ABC)"
            base_display = structure_type.display_name.replace(" (Multiple)", "")
            display_name = f"{base_display} ({pdb_id})"

            # Create a unique workspace key for this item (for identification)
            item_key = f"{structure_type.workspace_key}[{idx}]"

            # Create a temporary structure type for this specific item
            # Note: We don't register this - it's just for display
            item_type = StructureType.__new__(StructureType)
            item_type.workspace_key = item_key
            item_type.display_name = display_name
            item_type.priority = structure_type.priority
            item_type.description = f"Item {idx + 1} from {structure_type.workspace_key}"
            # Skip __post_init__ registration by using __new__

            info = StructureInfo(
                structure_type=item_type,
                file_path=file_path
            )

            if not require_file_path or info.exists:
                expanded.append(info)

        return expanded

    def get_available_structures(
        self,
        exclude_keys: Optional[List[str]] = None,
        include_legacy: bool = True,
        min_priority: Optional[int] = None,
        max_priority: Optional[int] = None,
        source_filter: Optional[str] = None,
        priority_override: Optional[List[str]] = None,
        require_file_path: bool = True,
        expand_lists: bool = True,
    ) -> List[StructureInfo]:
        """
        Get all available structures from workspace.

        Args:
            exclude_keys: List of workspace keys to exclude
            include_legacy: If False, exclude "pdb_file" key
            min_priority: Minimum priority level to consider
            max_priority: Maximum priority level to consider
            source_filter: Filter by source category (experimental, predicted, processed, aligned)
            priority_override: Custom priority order as list of workspace keys
            require_file_path: If True, only include entries where file exists
            expand_lists: If True (default), expand list-type keys (e.g., rcsb_pdb_files)
                         into individual StructureInfo objects for each item

        Returns:
            List of StructureInfo objects for available structures (sorted by priority)
        """
        exclude_keys = list(exclude_keys) if exclude_keys else []
        if not include_legacy and "pdb_file" not in exclude_keys:
            exclude_keys.append("pdb_file")

        # Apply source filter to exclude_keys
        if source_filter is not None:
            if source_filter not in SOURCE_FILTERS:
                raise ValueError(
                    f"Unknown source_filter: '{source_filter}'. "
                    f"Valid options: {list(SOURCE_FILTERS.keys())}"
                )
            # Get allowed keys for this filter
            allowed_keys = SOURCE_FILTERS[source_filter]
            # All registered types
            all_types = StructureRegistry.get_all()
            # Exclude any key not in the allowed list
            for st in all_types:
                if st.workspace_key not in allowed_keys and st.workspace_key not in exclude_keys:
                    exclude_keys.append(st.workspace_key)

        # Handle priority_override - completely custom ordering
        if priority_override is not None:
            available = []
            for idx, key in enumerate(priority_override):
                if key in exclude_keys:
                    continue

                # Get or create structure type for this key
                structure_type = StructureRegistry.get(key)
                if structure_type is None:
                    # Create a temporary structure type for unregistered keys
                    structure_type = StructureType(
                        workspace_key=key,
                        display_name=key.replace("_", " ").title(),
                        priority=idx + 1,  # Use position in list as priority
                        description="Custom priority override"
                    )

                # Check cache first
                cache_key = key
                if cache_key in self._cache:
                    info = self._cache[cache_key]
                    if not require_file_path or info.exists:
                        available.append(info)
                    continue

                # Get from workspace
                file_path = self.workspace.get(key)
                if file_path:
                    info = StructureInfo(
                        structure_type=structure_type,
                        file_path=file_path
                    )
                    self._cache[cache_key] = info

                    if not require_file_path or info.exists:
                        available.append(info)

            return available

        # Standard priority-based selection from registry
        if min_priority is not None or max_priority is not None:
            min_p = min_priority or 1
            max_p = max_priority or 999
            structure_types = StructureRegistry.get_by_priority_range(min_p, max_p)
        else:
            structure_types = StructureRegistry.get_all()

        available = []
        for structure_type in structure_types:
            if structure_type.workspace_key in exclude_keys:
                continue

            # Check cache first
            cache_key = structure_type.workspace_key
            if cache_key in self._cache:
                info = self._cache[cache_key]
                if not require_file_path or info.exists:
                    available.append(info)
                continue

            # Get from workspace
            value = self.workspace.get(structure_type.workspace_key)
            if value is None:
                continue

            # Handle list-type values (e.g., rcsb_pdb_files)
            if isinstance(value, list) and expand_lists:
                # Expand list into individual StructureInfo objects
                expanded = self._expand_list_key(structure_type, value, require_file_path)
                available.extend(expanded)
            elif isinstance(value, str):
                # Standard single-file handling
                info = StructureInfo(
                    structure_type=structure_type,
                    file_path=value
                )
                self._cache[cache_key] = info

                if not require_file_path or info.exists:
                    available.append(info)

        return available

    def _priority_selection(
        self,
        available: List[StructureInfo],
        silent: bool = False,
        return_key: bool = False
    ) -> Optional[Union[str, Tuple[str, str]]]:
        """
        Select structure based on priority (highest first).

        Args:
            available: List of available structures
            silent: If True, suppress console output
            return_key: If True, return tuple of (file_path, workspace_key)

        Returns:
            Path to highest-priority structure, or tuple of (path, key) if return_key=True
        """
        if not available:
            return None

        # Already sorted by priority in get_available_structures
        selected = available[0]

        if not silent:
            self.console.print(
                f"[green]Using {selected.structure_type.display_name} PDB file: "
                f"{selected.file_path}[/green]"
            )

        logger.debug(
            f"Selected {selected.structure_type.display_name} PDB: "
            f"{selected.file_path}"
        )

        if return_key:
            return selected.file_path, selected.structure_type.workspace_key
        return selected.file_path

    def _interactive_selection(
        self,
        available: List[StructureInfo],
        silent: bool = False,
        return_key: bool = False
    ) -> Optional[Union[str, Tuple[str, str]]]:
        """
        Prompt user to select from available structures.

        Args:
            available: List of available structures
            silent: If True, suppress console output (forces priority selection)
            return_key: If True, return tuple of (file_path, workspace_key)

        Returns:
            Path to user-selected structure, or tuple of (path, key) if return_key=True
        """
        if silent:
            # Can't do interactive selection silently, fall back to priority
            return self._priority_selection(available, silent=True, return_key=return_key)

        if len(available) == 1:
            # Only one option, no need to ask
            return self._priority_selection(available, silent=False, return_key=return_key)

        # Display table of available structures
        table = Table(title="Available Structures in Workspace")
        table.add_column("Option", style="cyan", justify="center")
        table.add_column("Type", style="green")
        table.add_column("File Path", style="yellow")
        table.add_column("Size", style="magenta", justify="right")

        options_map = {}
        for idx, info in enumerate(available, 1):
            # Format file size
            if info.file_size:
                if info.file_size > 1024 * 1024:
                    size_str = f"{info.file_size / (1024*1024):.2f} MB"
                elif info.file_size > 1024:
                    size_str = f"{info.file_size / 1024:.2f} KB"
                else:
                    size_str = f"{info.file_size} bytes"
            else:
                size_str = "Unknown"

            table.add_row(
                str(idx),
                info.structure_type.display_name,
                info.file_path,
                size_str
            )
            options_map[str(idx)] = info

        self.console.print(table)

        # Prompt for selection with session recording context
        if self.processor:
            # Use session-recording-aware prompt
            from proprep.application.menu_commands import prompt_with_context

            # Build options map for session recording
            session_options_map = {
                str(idx): info.structure_type.display_name
                for idx, info in enumerate(available, 1)
            }

            choice = prompt_with_context(
                processor=self.processor,
                prompt="Select structure to use",
                choices=[str(i) for i in range(1, len(available) + 1)],
                default="1",
                module="Structure Selector",
                description="Select structure from workspace",
                options_map=session_options_map
            )
        else:
            # Simple prompt without session recording
            choice = prompt_with_context(None,
                "\nSelect structure",
                choices=[str(i) for i in range(1, len(available) + 1)],
                default="1"
            )

        selected = options_map[choice]
        self.console.print(
            f"[green]✓ Selected {selected.structure_type.display_name}: "
            f"{selected.file_path}[/green]"
        )

        if return_key:
            return selected.file_path, selected.structure_type.workspace_key
        return selected.file_path

    def _interactive_multi_selection(self, available) -> List[str]:
        """
        Prompt user to select one or multiple structures from available options.

        Supports multiple input formats:
        - "1" - single selection
        - "1,3,5" - comma-separated
        - "1-3" - range notation
        - "all" - select all structures

        Args:
            available: List of available structures

        Returns:
            List of selected file paths (can be single or multiple)
        """
        if len(available) == 1:
            # Only one option, auto-select it
            self.console.print(
                f"[green]✓ Selected {available[0].structure_type.display_name}: "
                f"{available[0].file_path}[/green]"
            )
            return [available[0].file_path]

        # Display table of available structures
        table = Table(title="Available Structures in Workspace")
        table.add_column("Option", style="cyan", justify="center")
        table.add_column("Type", style="green")
        table.add_column("File Path", style="yellow")
        table.add_column("Size", style="magenta", justify="right")

        for idx, info in enumerate(available, 1):
            # Format file size
            if info.file_size:
                if info.file_size > 1024 * 1024:
                    size_str = f"{info.file_size / (1024*1024):.2f} MB"
                elif info.file_size > 1024:
                    size_str = f"{info.file_size / 1024:.2f} KB"
                else:
                    size_str = f"{info.file_size} bytes"
            else:
                size_str = "Unknown"

            table.add_row(
                str(idx),
                info.structure_type.display_name,
                info.file_path,
                size_str
            )

        self.console.print(table)
        self.console.print("\n[grey50]Enter selection: number (e.g., 1), range (e.g., 1-3), comma-separated (e.g., 1,3,5), or 'all'[/grey50]")

        # Get selection with session recording - retry loop for invalid input
        from proprep.utils.prompts import prompt_with_context

        multi_options_map = {
            str(idx): f"{info.structure_type.display_name} ({info.file_path})"
            for idx, info in enumerate(available, 1)
        }
        multi_options_map["all"] = "All available structures"

        selected_files = []
        while not selected_files:
            if self.processor:
                from proprep.application.menu_commands import prompt_with_context

                selection = prompt_with_context(
                    processor=self.processor,
                    prompt="Select structure(s)",
                    default="1",
                    module="Structure Selector",
                    description="Select structure(s) from workspace (number, range, comma list, or 'all')",
                    options_map=multi_options_map,
                )
            else:
                selection = prompt_with_context(None, "Select structure(s)", default="1")

            # Parse the selection
            if selection.lower() == "all":
                # Select all structures
                for info in available:
                    selected_files.append(info.file_path)
                self.console.print(f"[green]✓ Selected all {len(selected_files)} structures[/green]")
            else:
                # Parse comma-separated and/or ranges
                try:
                    indices = set()  # Use set to avoid duplicates
                    parts = selection.split(",")

                    for part in parts:
                        part = part.strip()
                        if "-" in part:
                            # Range notation (e.g., "1-3")
                            start, end = part.split("-", 1)
                            start_idx = int(start.strip())
                            end_idx = int(end.strip())
                            for idx in range(start_idx, end_idx + 1):
                                if 1 <= idx <= len(available):
                                    indices.add(idx)
                        else:
                            # Single number
                            idx = int(part)
                            if 1 <= idx <= len(available):
                                indices.add(idx)

                    if not indices:
                        # No valid indices found
                        self.console.print(f"[yellow]Invalid selection. Please enter numbers between 1 and {len(available)}[/yellow]")
                        continue

                    # Convert indices to file paths
                    for idx in sorted(indices):
                        selected_files.append(available[idx - 1].file_path)

                    if len(selected_files) == 1:
                        self.console.print(
                            f"[green]✓ Selected {available[list(indices)[0] - 1].structure_type.display_name}: "
                            f"{selected_files[0]}[/green]"
                        )
                    else:
                        self.console.print(f"[green]✓ Selected {len(selected_files)} structures[/green]")

                except (ValueError, IndexError) as e:
                    self.console.print(f"[yellow]Invalid input format. Please try again.[/yellow]")
                    continue

        return selected_files

    def display_available_structures(self) -> None:
        """Display all available structures in workspace (for debugging/info)"""
        available = self.get_available_structures()

        if not available:
            self.console.print("[yellow]No structures found in workspace[/yellow]")
            return

        table = Table(title="Workspace Structures")
        table.add_column("Priority", style="cyan", justify="center")
        table.add_column("Type", style="green")
        table.add_column("Workspace Key", style="blue")
        table.add_column("File Path", style="yellow")
        table.add_column("Exists", style="magenta", justify="center")

        for info in available:
            table.add_row(
                str(info.structure_type.priority),
                info.structure_type.display_name,
                info.structure_type.workspace_key,
                info.file_path,
                "✓" if info.exists else "✗"
            )

        self.console.print(table)

    def clear_cache(self):
        """Clear the structure info cache (call if workspace changes)"""
        self._cache.clear()

    def get_structure_by_key(
        self,
        key: str,
        fallback_keys: Optional[List[str]] = None,
        require_exists: bool = True,
    ) -> Optional[Union[str, Any]]:
        """
        Get structure by specific workspace key name.

        Use this when you need a SPECIFIC structure type, not priority-based selection.
        This is useful for mode-based selection where the user has chosen which
        pipeline stage to view/use (e.g., VMD visualizer modes).

        Supports indexed access for list keys:
        - "rcsb_pdb_files[0]" - first item from list
        - "rcsb_pdb_files[-1]" - last item from list

        Args:
            key: Primary workspace key to retrieve (supports indexed access)
            fallback_keys: List of keys to try if primary not found (in order)
            require_exists: If True and value is a file path string, verify file exists

        Returns:
            Value from workspace (file path string or BioPython Structure object),
            or None if not found

        Example:
            # VMD visualizer - get specific structure for visualization mode
            selector = StructureSelector(workspace, console)
            if mode == "repaired_only":
                structure = selector.get_structure_by_key(
                    "repaired_structure",
                    fallback_keys=["filtered_structure", "structure"]
                )

            # Access specific item from list
            first_rcsb = selector.get_structure_by_key("rcsb_pdb_files[0]")
        """
        # Check for indexed access pattern (e.g., "rcsb_pdb_files[0]")
        if "[" in key and key.endswith("]"):
            base_key, index_str = key.rsplit("[", 1)
            index_str = index_str.rstrip("]")

            try:
                index = int(index_str)
                value_list = self.workspace.get(base_key)

                if isinstance(value_list, list) and abs(index) < len(value_list):
                    value = value_list[index]
                    if require_exists and isinstance(value, str):
                        if os.path.exists(value):
                            return value
                        # Fall through to try fallbacks
                    else:
                        return value
            except (ValueError, IndexError):
                pass
            # Fall through to try fallbacks
        else:
            # Standard key lookup
            value = self.workspace.get(key)
            if value is not None:
                if require_exists and isinstance(value, str):
                    if os.path.exists(value):
                        return value
                    # File path but doesn't exist - try fallbacks
                else:
                    # Not a string (BioPython object) or exists check disabled
                    return value

        # Try fallback keys in order
        if fallback_keys:
            for fallback_key in fallback_keys:
                value = self.workspace.get(fallback_key)
                if value is not None:
                    if require_exists and isinstance(value, str):
                        if os.path.exists(value):
                            return value
                        # Continue to next fallback
                    else:
                        return value

        return None

    def get_structure_status(self) -> Dict[str, Any]:
        """
        Get comprehensive status of available structures in workspace.

        This is useful for menu display and determining which modules/options
        should be available based on current workspace state.

        Returns:
            Dict with:
                - has_any: bool - Any structure available
                - has_experimental: bool - RCSB or local PDB available
                - has_predicted: bool - AlphaFold/AlphaFill available
                - has_filtered: bool - Filtered structure available
                - has_repaired: bool - Repaired structure available
                - has_transformed: bool - Transformed structure available
                - has_protonation: bool - Protonation-updated structure available
                - available_keys: List[str] - All available workspace keys
                - highest_priority_key: Optional[str] - Key of highest priority structure
                - count: int - Total number of available structures

        Example:
            selector = StructureSelector(workspace)
            status = selector.get_structure_status()

            if status["has_any"]:
                # Enable structure processing menu options
                if status["has_filtered"]:
                    # Show repair option
                    ...
        """
        available = self.get_available_structures()
        available_keys = [info.structure_type.workspace_key for info in available]

        # Check each category
        experimental_keys = SOURCE_FILTERS["experimental"]
        predicted_keys = SOURCE_FILTERS["predicted"]

        status = {
            "has_any": len(available) > 0,
            "has_experimental": any(k in available_keys for k in experimental_keys),
            "has_predicted": any(k in available_keys for k in predicted_keys),
            "has_filtered": "filtered_pdb_file" in available_keys,
            "has_repaired": "repaired_pdb_file" in available_keys,
            "has_transformed": "transformed_pdb_file" in available_keys,
            "has_protonation": (
                "protonation_pdb_file" in available_keys or
                "structure_with_prot_resnames" in available_keys  # Legacy support
            ),
            "available_keys": available_keys,
            "highest_priority_key": available[0].structure_type.workspace_key if available else None,
            "count": len(available),
        }

        return status


# ============================================================================
# Backward-Compatible Helper Functions
# ============================================================================

def get_priority_pdb_file(
    workspace,
    console: Optional[Console] = None,
    silent: bool = False
) -> Optional[str]:
    """
    Get PDB file using priority-based selection.

    This is a backward-compatible wrapper around StructureSelector
    that maintains the same API as the old scattered implementations.

    Args:
        workspace: ProPrep workspace object
        console: Optional Rich console for output
        silent: If True, suppress console output

    Returns:
        Path to selected PDB file, or None if no valid file found
    """
    selector = StructureSelector(workspace, console)
    return selector.get_structure(interactive=False, silent=silent)


def get_interactive_pdb_file(
    workspace,
    console: Optional[Console] = None,
    processor=None,
    exclude_keys: Optional[List[str]] = None,
    return_key: bool = False
) -> Optional[Union[str, Tuple[str, str]]]:
    """
    Get PDB file with interactive user selection.

    Prompts user to choose from available structures in workspace.
    Falls back to priority selection if only one structure available.

    Args:
        workspace: ProPrep workspace object
        console: Optional Rich console for output
        processor: Optional processor for session recording context
        exclude_keys: Optional list of workspace keys to exclude from selection
        return_key: If True, return tuple of (file_path, workspace_key)

    Returns:
        Path to selected PDB file, or tuple of (path, key) if return_key=True, or None if no valid file found
    """
    selector = StructureSelector(workspace, console, processor)
    return selector.get_structure(interactive=True, exclude_keys=exclude_keys, return_key=return_key)


def get_available_structure_keys(workspace) -> List[str]:
    """
    Get list of available structure workspace keys.

    Args:
        workspace: ProPrep workspace object

    Returns:
        List of workspace keys for structures present in workspace
    """
    selector = StructureSelector(workspace)
    available = selector.get_available_structures()
    return [info.structure_type.workspace_key for info in available]


# ============================================================================
# Structure Object Selection (not just file paths)
# ============================================================================

def get_priority_structure_object(
    workspace,
    console: Optional[Console] = None,
    silent: bool = False,
    structure_key_suffix: str = "_structure"
) -> Optional[Any]:
    """
    Get Structure object (BioPython Structure) using priority-based selection.

    This function looks for Structure objects in workspace, not file paths.
    It checks keys like: original_structure, alphafold_structure, etc.

    Args:
        workspace: ProPrep workspace object
        console: Optional Rich console for output
        silent: If True, suppress console output
        structure_key_suffix: Suffix for structure object keys (default "_structure")

    Returns:
        BioPython Structure object, or None if no valid structure found
    """
    if console is None:
        console = Console()

    # Get all registered structure types sorted by priority
    structure_types = StructureRegistry.get_all()

    for structure_type in structure_types:
        # Convert file path key to structure object key
        # e.g., "original_pdb_file" -> "original_structure"
        base_key = structure_type.workspace_key.replace("_pdb_file", "").replace("_file", "")
        structure_obj_key = f"{base_key}{structure_key_suffix}"

        structure_obj = workspace.get(structure_obj_key)
        if structure_obj is not None:
            if not silent:
                console.print(
                    f"[green]Using {structure_type.display_name} Structure object "
                    f"(workspace key: {structure_obj_key})[/green]"
                )
            logger.info(f"Selected {structure_type.display_name} Structure: {structure_obj_key}")
            return structure_obj

    # No valid structure found
    if not silent:
        console.print("[red]No valid Structure object found in workspace[/red]")
        console.print("[yellow]Please load a structure first[/yellow]")
    logger.error("No valid Structure object found in workspace")
    return None


def get_interactive_structure_object(
    workspace,
    console: Optional[Console] = None,
    structure_key_suffix: str = "_structure",
    processor=None
) -> tuple[Optional[Any], Optional[str]]:
    """
    Get Structure object with interactive user selection.

    Prompts user to choose from available Structure objects in workspace.
    Falls back to priority selection if only one structure available.

    Args:
        workspace: ProPrep workspace object
        console: Optional Rich console for output
        structure_key_suffix: Suffix for structure object keys (default "_structure")
        processor: Optional processor for session recording context

    Returns:
        Tuple of (BioPython Structure object, corresponding file path key) or (None, None)
        The file path key can be used to get the PDB file: workspace.get(file_path_key)
    """
    if console is None:
        console = Console()

    # Find all available structure objects
    structure_types = StructureRegistry.get_all()
    available = []

    for structure_type in structure_types:
        base_key = structure_type.workspace_key.replace("_pdb_file", "").replace("_file", "")
        structure_obj_key = f"{base_key}{structure_key_suffix}"

        structure_obj = workspace.get(structure_obj_key)
        if structure_obj is not None:
            available.append({
                'type': structure_type,
                'key': structure_obj_key,
                'object': structure_obj,
                'file_key': structure_type.workspace_key  # The original file path key
            })

    if not available:
        console.print("[red]No Structure objects found in workspace[/red]")
        console.print("[yellow]Please load a structure first[/yellow]")
        return None, None

    if len(available) == 1:
        # Only one structure, use it automatically
        selected = available[0]
        console.print(
            f"[green]Using {selected['type'].display_name} Structure object "
            f"(workspace key: {selected['key']})[/green]"
        )
        return selected['object'], selected['file_key']

    # Multiple structures - let user choose
    table = Table(title="Available Structures in Workspace")
    table.add_column("Option", style="cyan", justify="center")
    table.add_column("Type", style="green")
    table.add_column("Workspace Key", style="yellow")
    table.add_column("Description", style="white")

    options_map = {}
    for idx, item in enumerate(available, 1):
        table.add_row(
            str(idx),
            item['type'].display_name,
            item['key'],
            item['type'].description or "N/A"
        )
        options_map[str(idx)] = item

    console.print(table)

    # Build options map for session recording
    session_options_map = {
        str(idx): f"{item['type'].display_name} ({item['key']})"
        for idx, item in enumerate(available, 1)
    }

    choice = prompt_with_context(
        processor,
        "\nSelect structure to use",
        choices=[str(i) for i in range(1, len(available) + 1)],
        default="1",
        module="Structure Selector",
        description="Select structure to use for processing",
        options_map=session_options_map
    )

    selected = options_map[choice]
    console.print(
        f"[green]✓ Selected {selected['type'].display_name} Structure "
        f"(workspace key: {selected['key']})[/green]"
    )

    return selected['object'], selected['file_key']


# ============================================================================
# Dynamic Registration API
# ============================================================================

def register_structure_type(
    workspace_key: str,
    display_name: str,
    priority: int,
    description: str = ""
) -> StructureType:
    """
    Register a new custom structure type dynamically.

    This allows modules to add new structure types without modifying this file.
    The new type will be integrated into the priority ordering.

    Args:
        workspace_key: Workspace key name (e.g., "qmmm_optimized_structure")
        display_name: Human-readable name (e.g., "QM/MM Optimized")
        priority: Priority level (1=highest, higher numbers=lower priority)
        description: Optional detailed description

    Returns:
        The created StructureType instance

    Example:
        # Register in your module's __init__.py or module class:
        from proprep.utils.structure_selector import register_structure_type

        register_structure_type(
            workspace_key="qmmm_optimized_structure",
            display_name="QM/MM Optimized",
            priority=2,  # Between PROTONATION_UPDATED (1) and TRANSFORMED (2)
            description="After QM/MM geometry optimization"
        )

        # Now StructureSelector will automatically find and prioritize it!
    """
    # Check if already registered
    existing = StructureRegistry.get(workspace_key)
    if existing:
        logger.warning(
            f"Structure type '{workspace_key}' already registered with priority "
            f"{existing.priority}. Overwriting with new priority {priority}."
        )

    # Create and register (auto-registers via __post_init__)
    structure_type = StructureType(
        workspace_key=workspace_key,
        display_name=display_name,
        priority=priority,
        description=description
    )

    logger.debug(
        f"Registered structure type: {workspace_key} "
        f"(display: {display_name}, priority: {priority})"
    )

    return structure_type


def update_structure_priority(workspace_key: str, new_priority: int):
    """
    Update the priority of an existing structure type.

    Args:
        workspace_key: The workspace key of the structure type
        new_priority: New priority value

    Raises:
        KeyError: If structure type not found
    """
    structure_type = StructureRegistry.get(workspace_key)
    if not structure_type:
        raise KeyError(f"Structure type '{workspace_key}' not found in registry")

    old_priority = structure_type.priority
    structure_type.priority = new_priority

    logger.info(
        f"Updated priority for '{workspace_key}': {old_priority} -> {new_priority}"
    )
