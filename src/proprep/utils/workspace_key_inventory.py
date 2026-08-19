"""
Workspace Key Inventory

Collects and displays a complete inventory of workspace keys used across
all ProPrep modules. Useful for documentation (SI tables) and maintenance.

Usage:
    python -m proprep.utils.workspace_key_inventory [--format table|markdown|json]
"""

import argparse
import importlib
import json
import sys
from typing import Dict, List, Optional

from proprep.utils.module_registry import registry


# ── Module import list (mirrors pdbprocessor._import_required_modules) ──────

MODULES_TO_IMPORT = [
    # ── Structure loading & preparation ──
    "proprep.structure_prep.pdb_loader",
    "proprep.structure_prep.structure_loader",
    "proprep.structure_prep.biological_assembly",
    "proprep.structure_prep.structure_orientation",
    "proprep.structure_prep.interactive_structure_viewer",
    # ── Analysis & annotation ──
    "proprep.structure_prep.redox_detector_module",
    "proprep.emboss.emboss_module",
    "proprep.structure_prep.homology_searcher",
    "proprep.structure_prep.pdb_filter",
    "proprep.structure_prep.structure_alignment",
    "proprep.structure_prep.vmd_visualizer",
    "proprep.structure_prep.amino_acid_mutator",
    "proprep.structure_prep.structure_completeness",
    "proprep.structure_prep.protonation_state_analyzer",
    "proprep.structure_prep.md_restraint_manager",
    # ── Redox site preparation ──
    "proprep.redoxsite_prep.redoxsite_integration",
    # ── Force field parameterization ──
    "proprep.forcefield_prep.forcefield_parameterizer",
    "proprep.forcefield_prep.small_molecule_parameterizer",
    # ── Topology & solvation ──
    "proprep.tleap_prep.tleap_input_generator",
    "proprep.membrane_prep.membrane_builder",
    # ── Simulation setup ──
    "proprep.md_prep.molecular_dynamics_manager",
    # ── QM/MM ──
    "proprep.oniom_prep.oniom_qmmm_preparator",
    "proprep.orca_prep.orca_qmmm_preparator",
]


# ── Modules to exclude from inventory (dead code, deprecated) ────────────────

MODULES_TO_EXCLUDE = [
    "AltLoc Selector",
    "Disulfide Bond Detector",
]


# ── StructurePreprocessor keys (not a ProcessingModule, can't self-register) ─
# Only list keys that the preprocessor ORIGINATES, not keys it merely updates.

PREPROCESSOR_KEYS = [
    "structure_pdb_file",
    "preprocessing_output_dir",
    "preprocessing_protein_ff",
    "preprocessing_water_model",
    "preprocessing_organic_ff",
    "preprocessing_organometallic_ff",
    "preprocessing_isolated_metals",
    "preprocessing_metal_free_pdb",
    "preprocessing_residue_sequence_map",
    "prepared_pdb",
    "preprocessing_atom_data",
    "preprocessing_triage",
    "preprocessing_protein_input",
    "preprocessing_lib_files",
    "preprocessing_frcmod_files",
    "preprocessing_metal_reinsertion_map",
    "preprocessing_nonstandard_ff",
    "preprocessing_metal_types",
    "preprocessing_modified_aa",
    "redox_sites",
    "ligand_pdb_file",
    "ligand_resname",
]

PREPROCESSOR_NAME = "Structure Preprocessor (pipeline)"


# ── Key descriptions: type hint + one-line description ───────────────────────

KEY_DESCRIPTIONS: Dict[str, Dict[str, str]] = {
    # ── Structure loading ──
    "rcsb_pdb_file": {"type": "str (path)", "description": "Path to PDB file downloaded from RCSB"},
    "rcsb_structure": {"type": "Structure", "description": "BioPython Structure object from RCSB"},
    "rcsb_metadata": {"type": "dict", "description": "PDB metadata from RCSB"},
    "rcsb_pdb_files": {"type": "list[str]", "description": "Paths to all downloaded RCSB PDB files (batch)"},
    "rcsb_structures": {"type": "list[Structure]", "description": "All downloaded Structure objects (batch)"},
    "rcsb_metadata_list": {"type": "list[dict]", "description": "Metadata for all downloaded structures (batch)"},
    "rcsb_download_info": {"type": "list[dict]", "description": "Download info dicts (pdb_id, path, timestamp)"},
    "local_pdb_file": {"type": "str (path)", "description": "Path to locally loaded PDB file"},
    "local_structure": {"type": "Structure", "description": "BioPython Structure object from local file"},
    "local_metadata": {"type": "dict", "description": "Metadata extracted from local PDB file"},
    "alphafold_pdb_file": {"type": "str (path)", "description": "Path to AlphaFold Database structure"},
    "alphafold_structure": {"type": "Structure", "description": "BioPython Structure object from AlphaFold"},
    "alphafold_uniprot_id": {"type": "str", "description": "UniProt ID used for AlphaFold retrieval"},
    "alphafold_confidence": {"type": "dict", "description": "pLDDT scores and confidence metrics"},
    "alphafill_pdb_file": {"type": "str (path)", "description": "Path to AlphaFill structure with transplanted ligands"},
    "alphafill_structure": {"type": "Structure", "description": "BioPython Structure object from AlphaFill"},
    "alphafill_uniprot_id": {"type": "str", "description": "UniProt ID used for AlphaFill retrieval"},
    "alphafill_transplants": {"type": "list", "description": "List of transplanted ligands/cofactors"},
    "alphafill_metadata": {"type": "dict", "description": "AlphaFill metadata"},
    "biological_assembly_pdb_file": {"type": "str (path)", "description": "Path to biological assembly PDB file"},
    "biological_assembly_structure": {"type": "Structure", "description": "BioPython Structure for biological assembly"},
    "biological_assembly_metadata": {"type": "dict", "description": "Biological assembly metadata"},
    "hstripped_pdb_file": {"type": "str (path)", "description": "Path to PDB file with hydrogens stripped (for re-preprocessing)"},
    "hstripped_structure": {"type": "Structure", "description": "BioPython Structure with hydrogens stripped"},
    # ── PDB Loader (legacy) ──
    "original_pdb_file": {"type": "str (path)", "description": "Original PDB file path (canonical name)"},
    "original_structure": {"type": "Structure", "description": "Original BioPython Structure object"},
    "original_metadata": {"type": "dict", "description": "Original PDB metadata"},
    "original_residue_numbering_info": {"type": "dict", "description": "Original residue numbering information"},
    "pdb_file": {"type": "str (path)", "description": "Current PDB file path (deprecated, backward compat)"},
    "metadata": {"type": "dict", "description": "Current PDB metadata (deprecated, backward compat)"},
    "structure": {"type": "Structure", "description": "Current BioPython Structure (deprecated, backward compat)"},
    "residue_numbering_info": {"type": "dict", "description": "Residue numbering info (deprecated, backward compat)"},
    # ── Homology search ──
    "blast_results": {"type": "list", "description": "BLAST search results"},
    "blast_raw_results": {"type": "dict", "description": "Raw BLAST XML results"},
    "alphafold_homologs": {"type": "list", "description": "AlphaFold homolog entries"},
    "alphafold_homolog_pdb_file": {"type": "str (path)", "description": "Path to selected AlphaFold homolog PDB"},
    "alphafold_homolog_structure": {"type": "Structure", "description": "BioPython Structure for AlphaFold homolog"},
    "homology_model_pdb_file": {"type": "str (path)", "description": "Path to MODELLER-built homology model PDB"},
    "homology_model_structure": {"type": "Structure", "description": "BioPython Structure for homology model"},
    "homology_model_alignment": {"type": "dict", "description": "Alignment info (identity, gaps, score per chain, template)"},
    # ── Structure alignment ──
    "alignment_results": {"type": "dict", "description": "Structural alignment results (RMSD, etc.)"},
    "alignment_residues": {"type": "list", "description": "Residue-level alignment data"},
    "aligned_target_pdb_file": {"type": "str (path)", "description": "Path to aligned target PDB file"},
    "aligned_target_structure": {"type": "Structure", "description": "Aligned target Structure object"},
    "aligned_ref_pdb_file": {"type": "str (path)", "description": "Path to aligned reference PDB file"},
    "aligned_ref_structure": {"type": "Structure", "description": "Aligned reference Structure object"},
    "aligned_structures": {"type": "list", "description": "All aligned structures"},
    # ── Structure filtering ──
    "filtered_structure": {"type": "Structure", "description": "Structure after chain/model/residue filtering"},
    "filtered_pdb_file": {"type": "str (path)", "description": "Path to filtered PDB file"},
    "filter_selections": {"type": "dict", "description": "User's filter selections (chains, models, residues)"},
    "is_hplusplus_structure": {"type": "bool", "description": "Whether structure was protonated by H++"},
    # ── Structure orientation ──
    "oriented_pdb_file": {"type": "str (path)", "description": "Path to oriented PDB file"},
    "orientation_record": {"type": "dict", "description": "Transformation record (rotation matrix, parameters)"},
    # ── AltLoc / legacy ──
    "altloc_processing_results": {"type": "dict", "description": "AltLoc processing results (legacy)"},
    "processed_structure": {"type": "Structure", "description": "Structure after AltLoc selection (legacy)"},
    # ── Disulfide bonds (legacy) ──
    "disulfide_bonds": {"type": "list", "description": "Detected disulfide bonds (legacy)"},
    "disulfide_structure": {"type": "Structure", "description": "Structure with disulfide bonds (legacy)"},
    "disulfide_tleap_commands": {"type": "list[str]", "description": "TLEaP bond commands for disulfides (legacy)"},
    # ── Structure completeness / repair ──
    "completeness_results": {"type": "dict", "description": "Structure completeness analysis results"},
    "repaired_structure": {"type": "Structure", "description": "Structure after missing residue/loop repair"},
    "repaired_pdb_file": {"type": "str (path)", "description": "Path to repaired PDB file"},
    "mutations_applied": {"type": "list", "description": "Standard mutations applied during repair"},
    "nonstandard_mutations_applied": {"type": "list", "description": "Non-standard mutations applied during repair"},
    # ── Amino acid mutations ──
    "pending_mutations": {"type": "list", "description": "Queued standard amino acid mutations"},
    "pending_nonstandard_mutations": {"type": "list", "description": "Queued non-standard amino acid mutations"},
    # ── Redox site detection ──
    "detected_redox_sites": {"type": "list[dict]", "description": "Detected redox-active sites with metadata"},
    "redox_transformer_mappings": {"type": "dict", "description": "Mappings for redox state transformers"},
    "remove_hydrogens_for_md": {"type": "bool", "description": "Flag to remove hydrogens before MD preprocessing"},
    "redox_sites": {"type": "list[dict]", "description": "Synced redox sites after preprocessing"},
    # ── Redox site preparation ──
    "untransformed_structure": {"type": "Structure", "description": "Structure before redox transformation"},
    "transformed_pdb_file": {"type": "str (path)", "description": "Path to redox-transformed PDB file"},
    "transformed_structure": {"type": "Structure", "description": "Structure after redox transformation"},
    "generated_microstate_pdbs": {"type": "list[str]", "description": "Paths to generated microstate PDB files"},
    "metal_sites": {"type": "list[dict]", "description": "Metal site definitions from redox preparation"},
    "excluded_residues": {"type": "list", "description": "Residues excluded during redox preparation"},
    # ── Protonation state analysis ──
    "protonation_results": {"type": "dict", "description": "Protonation state analysis results"},
    "microstate_protonation_results": {"type": "dict", "description": "Microstate-specific protonation results"},
    "protonation_method": {"type": "str", "description": "Method used for protonation (propka, pdb2pqr, etc.)"},
    "md_residue_names": {"type": "dict", "description": "MD-compatible residue name mappings (HIE/HID/HIP, etc.)"},
    "net_charge": {"type": "float", "description": "Net charge of the system"},
    "textbook_pkas": {"type": "dict", "description": "Textbook pKa values used as reference"},
    "propka_pka_values": {"type": "dict", "description": "PROPKA-calculated pKa values"},
    "pdb2pqr_output_file": {"type": "str (path)", "description": "Path to PDB2PQR output file"},
    "pdb2pqr_pqr_file": {"type": "str (path)", "description": "Path to PQR file from PDB2PQR"},
    "protonation_pdb_file": {"type": "str (path)", "description": "Path to protonated PDB file"},
    "pdb2pqr_target_ph": {"type": "float", "description": "Target pH used for PDB2PQR protonation"},
    "constant_ph_residues": {"type": "list", "description": "Residues selected for constant-pH MD"},
    "constant_ph_data": {"type": "dict", "description": "Constant-pH simulation configuration data"},
    "ideal_capacitance": {"type": "float", "description": "Ideal capacitance from protonation analysis"},
    "specific_capacitance": {"type": "float", "description": "Specific capacitance from protonation analysis"},
    "capacitance_per_group": {"type": "dict", "description": "Per-group capacitance values"},
    "capacitance_profile": {"type": "dict", "description": "Capacitance profile data"},
    # ── Forcefield parameterization ──
    "parameterized_residues": {"type": "dict", "description": "Parameterized non-standard residue data"},
    "non_standard_residues": {"type": "list", "description": "Detected non-standard residues"},
    "user_residue_classifications": {"type": "dict", "description": "User classifications for non-standard residues"},
    "pending_parameterizations": {"type": "list", "description": "Queued parameterization tasks"},
    "global_atom_registry_data": {"type": "dict", "description": "Global atom type registry for consistency"},
    "small_molecules": {"type": "list", "description": "Detected small molecules for parameterization"},
    # ── Structure preprocessing pipeline ──
    "structure_pdb_file": {"type": "str (path)", "description": "Input PDB file for preprocessing pipeline"},
    "preprocessing_output_dir": {"type": "str (path)", "description": "Output directory for preprocessing"},
    "preprocessing_protein_ff": {"type": "str", "description": "Selected protein force field (leaprc)"},
    "preprocessing_water_model": {"type": "str", "description": "Selected water model (leaprc)"},
    "preprocessing_organic_ff": {"type": "str", "description": "Selected organic molecule force field (leaprc)"},
    "preprocessing_organometallic_ff": {"type": "str", "description": "Selected organometallic force field"},
    "preprocessing_isolated_metals": {"type": "list", "description": "Isolated metal ions found in structure"},
    "preprocessing_metal_free_pdb": {"type": "str (path)", "description": "PDB file with metals removed"},
    "preprocessing_residue_sequence_map": {"type": "dict", "description": "Residue-to-sequence mapping after preprocessing"},
    "prepared_pdb": {"type": "str (path)", "description": "Final preprocessed PDB file"},
    "preprocessing_atom_data": {"type": "dict", "description": "Atom data from preprocessing steps"},
    "preprocessing_triage": {"type": "dict", "description": "Residue triage results (standard/nonstandard/metal)"},
    "preprocessing_protein_input": {"type": "str (path)", "description": "Protein PDB input for tLEaP"},
    "selected_standard_forcefields": {"type": "dict", "description": "Selected standard force fields for tLEaP"},
    "preprocessing_lib_files": {"type": "list[str]", "description": "Library files from parameterization"},
    "preprocessing_frcmod_files": {"type": "list[str]", "description": "Frcmod files from parameterization"},
    "preprocessing_metal_reinsertion_map": {"type": "dict", "description": "Map for reinserting metals after preprocessing"},
    "preprocessing_nonstandard_ff": {"type": "dict", "description": "Non-standard residue force field selections"},
    "preprocessing_metal_types": {"type": "dict", "description": "Metal type classifications for parameterization"},
    "preprocessing_modified_aa": {"type": "list", "description": "Selected modified amino acid leaprcs"},
    "solvation_parameters": {"type": "dict", "description": "Solvation model parameters (explicit/implicit)"},
    "ligand_pdb_file": {"type": "str (path)", "description": "Extracted ligand PDB file for parameterization"},
    "ligand_resname": {"type": "str", "description": "Residue name of extracted ligand"},
    # ── TLEaP topology generation ──
    "combined_tleap_commands": {"type": "list[str]", "description": "Combined tLEaP commands for topology generation"},
    "single_state_ff_requirements": {"type": "dict", "description": "Force field requirements for single-state topology"},
    "single_state_selected_forcefields": {"type": "dict", "description": "Selected force fields for single-state topology"},
    "tleap_parameters": {"type": "dict", "description": "tLEaP execution parameters"},
    "tleap_input_file": {"type": "str (path)", "description": "Path to generated tLEaP input file"},
    "tleap_template": {"type": "dict", "description": "tLEaP input template for single state"},
    "microstate_tleap_template": {"type": "dict", "description": "tLEaP input template for microstates"},
    "generated_microstate_tleap_files": {"type": "list[str]", "description": "Paths to generated microstate tLEaP files"},
    "parm7_file": {"type": "str (path)", "description": "Path to AMBER topology (.parm7) file"},
    "rst7_file": {"type": "str (path)", "description": "Path to AMBER coordinate (.rst7) file"},
    "_active_tleap_input_file": {"type": "str (path)", "description": "Currently active tLEaP input file (internal)"},
    "cpin_config": {"type": "dict", "description": "Constant-pH input configuration"},
    "cpin_file": {"type": "str (path)", "description": "Path to cpinutil-generated CPin file"},
    # ── MD restraints ──
    "disang_file": {"type": "str (path)", "description": "Path to exported DISANG restraint file"},
    "disang_export_results": {"type": "dict", "description": "DISANG export metadata (path, counts, status)"},
    "redox_restraint_mask": {"type": "str", "description": "AMBER mask for redox-site restraints"},
    "redox_restraint_info": {"type": "dict", "description": "Redox restraint metadata"},
    "restraint_structure_source": {"type": "str (path)", "description": "Source structure for restraint generation"},
    "restraint_mask_generated": {"type": "bool", "description": "Whether restraint mask was generated"},
    "restraints": {"type": "list[dict]", "description": "List of restraint definitions"},
    # ── MD simulation setup ──
    "md_simulation_queue": {"type": "list", "description": "Queue of MD simulation configurations"},
    "md_workflows": {"type": "list", "description": "Defined MD workflow sequences"},
    "md_template_assignments": {"type": "dict", "description": "Template-to-structure assignments for MD"},
    "restraint_integration_config": {"type": "dict", "description": "Restraint integration configuration for MD"},
    "topology_extracted_pdb": {"type": "str (path)", "description": "PDB extracted from topology for visualization"},
    "md_structure_pairs": {"type": "dict", "description": "Structure pairs (prmtop/rst7) for MD"},
    "preferred_amber_engine": {"type": "str", "description": "Selected AMBER engine (pmemd.MPI or pmemd.cuda)"},
    "mpi_tasks": {"type": "int", "description": "Number of MPI tasks for pmemd.MPI"},
    "gpu_ids": {"type": "str", "description": "GPU device IDs for pmemd.cuda"},
    # ── Membrane builder ──
    "membrane_packed_pdb": {"type": "str (path)", "description": "Path to packmol-assembled membrane-protein PDB"},
    "membrane_config": {"type": "dict", "description": "Membrane builder configuration (lipid composition, dimensions, etc.)"},
    "membrane_leaprc_requirements": {"type": "list[str]", "description": "leaprc files required by the membrane system (lipid FF)"},
    "membrane_box_dimensions": {"type": "dict", "description": "Box dimensions of the packed membrane system"},
    "membrane_ion_summary": {"type": "dict", "description": "Ion counts and types used to neutralize/salt the system"},
    "membrane_solutes": {"type": "list", "description": "Solute molecules embedded in the membrane system"},
    "is_membrane_system": {"type": "bool", "description": "Flag indicating the system includes a lipid bilayer"},
    # ── ONIOM QM/MM ──
    "oniom_setup": {"type": "dict", "description": "ONIOM QM/MM layer and calculation configuration"},
    "oniom_input_file": {"type": "str (path)", "description": "Path to Gaussian ONIOM input file"},
    # ── ORCA QM/MM ──
    "orca_qmmm_setup": {"type": "object", "description": "ORCA QM/MM setup configuration object"},
    "orca_input_file": {"type": "str (path)", "description": "Path to generated ORCA QM/MM input file"},
    "orca_ff_file": {"type": "str (path)", "description": "Path to ORCA force field parameter file (.prms)"},
    "orca_pdb_file": {"type": "str (path)", "description": "Path to PDB file prepared for ORCA QM/MM"},
    # ── EMBOSS analysis ──
    "emboss_sequence_analysis": {"type": "dict", "description": "EMBOSS sequence analysis results (pepstats, etc.)"},
    "emboss_alignments": {"type": "dict", "description": "EMBOSS pairwise alignment results"},
    "emboss_motifs": {"type": "dict", "description": "EMBOSS motif/pattern search results"},
    "emboss_batch_analysis": {"type": "dict", "description": "EMBOSS batch analysis results"},
    # ── Miscellaneous ──
    "topology": {"type": "str (path)", "description": "Path to topology file (required by MD Manager)"},
    "coordinates": {"type": "str (path)", "description": "Path to coordinate file (required by MD Manager)"},

    # ── PB Titrate: Poisson-Boltzmann titration / constant-pH ──
    "pb_titrate_params": {"type": "dict", "description": "PB solver settings (epsin, istrng, space, nfocus, bcopt), worker count, target pH, cutout radius"},
    "pb_titrate_sites": {"type": "list[dict]", "description": "Titratable sites after user and metal-coordination exclusions, with residue envelopes"},
    "pb_titrate_prmtop": {"type": "Path", "description": "Solvent-stripped topology used for all PB calculations (implicit solvent)"},
    "pb_titrate_rst7": {"type": "Path", "description": "Coordinates matching the solvent-stripped PB topology"},
    "pb_titrate_prmtop_solvated": {"type": "Path", "description": "Original un-stripped topology, retained for state application and ion rebalancing"},
    "pb_titrate_rst7_solvated": {"type": "Path", "description": "Coordinates matching the un-stripped topology"},
    "pb_titrate_self_energies": {"type": "str (path)", "description": "Pickle of per-site multistate energies and effective pKa values; rewritten per site to allow resume"},
    "pb_titrate_pka": {"type": "dict", "description": "Map of (resname, resnum) to single-site effective pKa"},
    "pb_titrate_pka_csv": {"type": "str (path)", "description": "Human-readable table of per-site effective pKa and dominant state"},
    "pb_titrate_coupling": {"type": "str (path)", "description": "Pickle of site self-energies and the pairwise site-site interaction matrix"},
    "pb_titrate_state_map": {"type": "dict", "description": "Working per-site protonation assignment, mapping (resname, resnum) to a constant-pH state name"},
    "pb_titrate_solver_result": {"type": "str", "description": "Name of the solver that produced the current state assignment"},
    "pb_titrate_unconverged_sites": {"type": "list", "description": "Sites still flipping after mean-field sweeps; seeds the targeted-coupling subgraph"},
    "pb_titrate_metal_fixed": {"type": "list[dict]", "description": "Titratable side chains held at fixed charge because they coordinate a metal ion"},
    "pb_titrate_terminal_excluded": {"type": "list[dict]", "description": "Chain-terminal residues excluded from titration (no multi-state library entries)"},
    "pb_titrate_minimized": {"type": "dict", "description": "Map of topology path to the minimized coordinates produced for it"},
    "pb_rename_pdb_file": {"type": "str (path)", "description": "PDB with titratable residues renamed to modern force-field codes for their PB-recommended states"},
    "pb_titrate_pka_vs_propka_csv": {"type": "str (path)", "description": "Comparison table of PB versus PROPKA pKa values"},
    "titrate_recommendations": {"type": "dict", "description": "Per-site recommended constant-pH state, proton count, pKa correction and net charge"},
    "titrate_report_csv": {"type": "str (path)", "description": "Human-readable per-site recommendation report"},
    "prmtop_titrated": {"type": "str (path)", "description": "Topology with recommended protonation-state charges applied and ions rebalanced"},
    "rst7_titrated": {"type": "str (path)", "description": "Coordinates written alongside the titrated topology"},
    "simulation_pH": {"type": "float", "description": "Intended simulation pH, read as the default target pH (externally seeded)"},
    "ionic_strength_mM": {"type": "float", "description": "Monovalent salt concentration seeding the PB ionic-strength prompt (externally seeded)"},
    "prmtop": {"type": "str (path)", "description": "Legacy alias for the input topology (read-only back-compatibility)"},
    "prmtop_file": {"type": "str (path)", "description": "Legacy alias for the input topology (read-only back-compatibility)"},
    "rst7": {"type": "str (path)", "description": "Legacy alias for the input coordinates (read-only back-compatibility)"},

    # ── Preprocessing (additional) ──
    "preprocessing_atom_types": {"type": "list[str]", "description": "tLEaP addAtomTypes entries for custom MCPB atom types"},
    "preprocessing_bond_commands": {"type": "list[str]", "description": "Explicit tLEaP bond commands for metal coordination, generated while the coordinate mapping is valid"},
    "preprocessing_metal_clusters": {"type": "dict", "description": "Pure inorganic clusters (triage category F) withheld whole from the standard force-field pass"},
    "_preprocessing_tleap_active": {"type": "bool", "description": "Re-entrancy flag marking the preprocessor's internal, boxless tLEaP build"},
    "redox_sites_pristine": {"type": "list[RedoxSite]", "description": "Untouched snapshot of redox sites, re-synchronized from on every rerun"},
    "mcpb_step_results": {"type": "dict", "description": "Per-SITE, per-step MCPB results (paths, QM charge, RESP intermediates) enabling resume; keyed by site_id"},
    "ligand_mol2_file": {"type": "str (path)", "description": "Charged mol2 for a single extracted ligand (externally seeded)"},
    "ligand_frcmod_file": {"type": "str (path)", "description": "GAFF frcmod for a single extracted ligand (externally seeded)"},
    "ff_resolved_atom_types": {"type": "list[str]", "description": "addAtomTypes entries emitted after resolving atom-type collisions between parameter sets"},
    "transformer_info": {"type": "list[dict]", "description": "Per-redox-site transformer assignments, states, atom types and resolved force-field set"},
    "microstate_metadata_path": {"type": "str (path)", "description": "Path to the JSON describing every generated redox microstate and its states"},
    "standard_ff_leaprcs_sourced": {"type": "set[str]", "description": "leaprc files the selected force fields will source, to avoid duplicate sourcing"},
    "reordered_pdb_file": {"type": "str (path)", "description": "PDB with chains reordered into tLEaP-acceptable order; highest-priority structure for topology building"},
    "reordering_skipped": {"type": "bool", "description": "Set when the user declines chain reordering, so later stages fall back"},
    "structure_with_prot_resnames": {"type": "str (path)", "description": "Legacy alias for the protonation-state-named PDB (superseded by protonation_pdb_file)"},

    # ── Structure sources and analysis (additional) ──
    "alphafold_all_models": {"type": "list[str]", "description": "Paths to every ranked model from a local prediction run"},
    "alphafold_config": {"type": "dict", "description": "User-configured AlphaFold prediction parameters"},
    "alphafold_output_dir": {"type": "str (path)", "description": "Directory holding local AlphaFold/ColabFold run outputs"},
    "alphafold_pae_matrix": {"type": "dict", "description": "Predicted aligned error data from the top-ranked model"},
    "alphafold_sequence_source": {"type": "str", "description": "Provenance label for the sequence submitted to prediction"},
    "alphafold_used_template": {"type": "str (path)", "description": "Custom template supplied to a template-based prediction"},
    "aligned_pdb_file": {"type": "str (path)", "description": "Legacy single-reference alignment result (read-only)"},
    "trajectory_files": {"type": "list[str]", "description": "Ordered AMBER trajectory segment paths loaded alongside a topology"},
    "propka_problematic_residues": {"type": "set[str]", "description": "Residues PROPKA could not parse, accumulated so they can be excluded on retry"},
    "residue_mapping": {"type": "dict", "description": "Original-to-repaired residue numbering map (read-only; no in-tree writer)"},
    "processor": {"type": "object", "description": "Defensive fallback handle to the application object (no in-tree writer)"},
    "_metadata_info_shown": {"type": "bool", "description": "One-shot flag so the PDB metadata hint is shown only once per session"},

    # ── MD and workflow orchestration (additional) ──
    "md_restraints": {"type": "list[MDRestraint]", "description": "Restraint objects carrying AMBER flat-bottom or parabolic parameters"},
    "group_restraints": {"type": "list[dict]", "description": "AMBER GROUP-format positional restraint specifications"},
    "slurm_mode": {"type": "bool", "description": "Set once SLURM batch scripts have been generated or submitted"},
    "slurm_output_dir": {"type": "str (path)", "description": "Directory holding generated SLURM job and submit scripts"},
    "workflow_preset": {"type": "str", "description": "Name of the MD workflow preset expanded into concrete input files"},
    "workflow_files": {"type": "list[str]", "description": "Generated .mdin input files, in step order"},
    "workflow_sequence": {"type": "str (path)", "description": "JSON defining step order and inter-step dependencies"},
    "workflow_dir": {"type": "str (path)", "description": "Directory into which workflow inputs and scripts were written"},
    "workflow_assignments": {"type": "list[dict]", "description": "Per-structure MD protocol assignments"},
    "workflow_current_stage": {"type": "str", "description": "Guided-workflow stage the session is currently on"},
    "workflow_completed_stages": {"type": "list[str]", "description": "Guided-workflow stages the user has completed"},
    "workflow_skipped_stages": {"type": "list[str]", "description": "Guided-workflow stages the user explicitly skipped"},
    "active_workflow": {"type": "dict", "description": "Descriptor of the MD protocol in effect (externally seeded)"},
    "md_metadata": {"type": "dict", "description": "Auxiliary MD bookkeeping, including custom workflow metadata (externally seeded)"},
    "md_save_json_backup": {"type": "bool", "description": "Opt-in flag to also write the simulation queue to disk (externally seeded)"},
    "files": {"type": "list[str]", "description": "Registered file names, scanned for restraint files (externally seeded)"},

    # ── Session, CLI and batch-mode flags ──
    "project_directory": {"type": "str (path)", "description": "Run's project directory, fixed before session recording starts"},
    "menu_mode": {"type": "str", "description": "Active menu interface mode (guided workflow or full menu)"},
    "menu_layout": {"type": "str", "description": "Per-run override of full-menu presentation (grid or list)"},
    "jump_to_analysis": {"type": "bool", "description": "One-shot flag opening the analysis menu at startup"},
    "jump_to_pdbview": {"type": "bool", "description": "One-shot flag opening the structure viewer at startup"},
    "pdbview_target": {"type": "str", "description": "PDB ID or path to open when the viewer is launched directly"},
    "pdb_id": {"type": "str", "description": "Four-character RCSB accession for the structure under preparation"},
    "output_dir": {"type": "str (path)", "description": "Base output directory (externally seeded; defaults to the current directory)"},
    "working_directory": {"type": "str (path)", "description": "Session working directory used to resolve output paths (externally seeded)"},
    "default_pH": {"type": "float", "description": "Default pH used when one is needed without prompting (externally seeded)"},
    "auto_filter": {"type": "bool", "description": "Batch-mode flag to run filtering non-interactively (externally seeded)"},
    "auto_repair": {"type": "bool", "description": "Batch-mode flag to run structure repair non-interactively (externally seeded)"},
    "auto_detect_metals": {"type": "bool", "description": "Batch-mode flag to run metal-site detection non-interactively (externally seeded)"},
    "analyze_protonation": {"type": "bool", "description": "Batch-mode flag to run protonation analysis (externally seeded)"},
    "run_propka": {"type": "bool", "description": "Batch-mode flag selecting PROPKA for pKa prediction (externally seeded)"},
    "set_md_names": {"type": "bool", "description": "Batch-mode flag to apply MD residue names (externally seeded)"},
    "constant_pH": {"type": "bool", "description": "Batch-mode flag selecting constant-pH residue naming (externally seeded)"},
    "update_structure": {"type": "bool", "description": "Batch-mode flag to write protonation changes into the structure (externally seeded)"},
    "custom_sequences": {"type": "dict", "description": "User-supplied sequences stored for reuse in EMBOSS analyses"},
    "mutated_sequence": {"type": "str", "description": "Mutated sequence offered as an EMBOSS analysis input (externally seeded)"},
    "filtered_structure_pdb_content": {"type": "str", "description": "Raw PDB text used as a last-resort source for restraint-mask generation (externally seeded)"},
}


# ── Component abbreviations (for the compact SI table) ───────────────────────
# Three-letter codes keep the Produced By / Required By columns narrow enough
# that the description column stays readable. Emitted with a legend by
# --abbrev; every component appearing in the table must have an entry here,
# which _abbrev() enforces rather than silently falling back to the full name.

_ABBREV = {
    "AlphaFold Predictor": "AFP",
    "AMBER Input Generator": "AIG",
    "AMBER Workflow Manager": "AWM",
    "Amino Acid Mutator": "AAM",
    "Application Core": "APP",
    "Batch Processor": "BAT",
    "Biological Assembly Generator": "BAG",
    "EMBOSS Analysis": "EMB",
    "Force Field Parameterizer": "FFP",
    "Frame Extractor": "FRX",
    "Homology Searcher": "HOM",
    "MD Restraint Manager": "MRM",
    "Membrane Builder": "MEM",
    "Metal Site Parameterizer": "MSP",
    "Modified Amino Acid Parameterizer": "MAP",
    "Molecular Dynamics Manager": "MDM",
    "ONIOM QM/MM Preparator": "ONI",
    "ORCA QM/MM Preparator": "ORC",
    "PB Titrate": "PBT",
    "PDB Filter": "FIL",
    "PDB Loader": "PDL",
    "Protonation State Analyzer": "PSA",
    "QM/MM Preparator": "QMM",
    "Redox Site Detector": "RSD",
    "Redox Site Preparer": "RSP",
    "Small Molecule Parameterizer": "SMP",
    "Structure Aligner": "ALN",
    "Structure Fixer": "FIX",
    "Structure Loader": "STL",
    "Structure Orientator": "ORI",
    PREPROCESSOR_NAME: "PRE",
    "Structure Viewer": "VIS",
    "Structure Viewer (VMD)": "VMD",
    "Topology Generator": "TOP",
    "Workflow State Manager": "WSM",
}


def _abbrev(name: str) -> str:
    """Abbreviation for a component name, or the name itself if unmapped."""
    return _ABBREV.get(name, name)


def _legend_for(inventory: Dict[str, Dict[str, List[str]]]) -> List[tuple]:
    """(code, full name) pairs for components actually present, code-sorted."""
    used = set()
    for info in inventory.values():
        used.update(info["produced_by"])
        used.update(info["required_by"])
    return sorted(((_abbrev(n), n) for n in used), key=lambda p: p[0])


# ── Inventory collection ─────────────────────────────────────────────────────

def import_all_modules() -> List[str]:
    """Import all ProPrep modules to trigger registration. Returns list of failures."""
    failures = []
    for module_name in MODULES_TO_IMPORT:
        try:
            importlib.import_module(module_name)
        except ImportError as e:
            failures.append(f"{module_name}: {e}")
    return failures


# ── Source-scan completeness pass ────────────────────────────────────────────
# The declaration pass below only sees keys that modules formally declare via
# get_workspace_outputs()/get_workspace_requirements() AND that appear in
# MODULES_TO_IMPORT. To GUARANTEE no key is silently omitted, also scan the
# package source for every string-literal workspace key actually set or read,
# and merge in any the declaration pass missed (e.g. undeclared module outputs,
# and keys set by non-module code: workflow state, settings, navigation flags).

import io
import os
import re
import tokenize

_ACRONYMS = {
    "pdb": "PDB", "tleap": "tLEaP", "md": "MD", "qm": "QM", "mm": "MM",
    "ff": "FF", "rcsb": "RCSB", "mcpb": "MCPB", "qmmm": "QM/MM", "vmd": "VMD",
    "emboss": "EMBOSS", "resp": "RESP", "oniom": "ONIOM", "orca": "ORCA",
    "blast": "BLAST", "altloc": "AltLoc", "pdb2pqr": "PDB2PQR", "ph": "pH",
    "gpu": "GPU", "mpi": "MPI", "id": "ID", "aa": "AA",
}


# Files whose humanized stem would be uninformative or would not match the
# component name used elsewhere in the table. Keyed by filename stem.
_FILE_LABELS = {
    "workflow": "PB Titrate",
    "workflow_state_manager": "Workflow State Manager",
    "structure_preprocessor": PREPROCESSOR_NAME,
    "pdbprocessor": "Application Core",
    "main": "Application Core",
    "menu_commands": "Application Core",
    "module_commands": "Application Core",
    "workspace_commands": "Application Core",
    "processor_command": "Application Core",
    "feedback_command": "Application Core",
    "batch_processor": "Batch Processor",
    "metallo_worker": "Redox Site Preparer",
    "metallo_integration": "Redox Site Preparer",
    "redox_transformation_manager": "Redox Site Preparer",
    "oniom_qmmm_preparator": "ONIOM QM/MM Preparator",
    "orca_qmmm_preparator": "ORCA QM/MM Preparator",
    "qmmm_preparator": "QM/MM Preparator",
    "amber_controller": "AMBER Workflow Manager",
    "workflow_centric_step1": "AMBER Workflow Manager",
    "modified_amino_acid_parameterizer": "Modified Amino Acid Parameterizer",
    "metal_site_parameterizer": "Metal Site Parameterizer",
    "tleap_input_generator": "Topology Generator",
    "membrane_builder": "Membrane Builder",
    "structure_completeness": "Structure Fixer",
    "sequence_input": "EMBOSS Analysis",
    # Collapse filename-derived labels onto the component's registered NAME,
    # so one component does not appear in the table under two spellings.
    "structure_alignment": "Structure Aligner",
    "structure_alignment_commands": "Structure Aligner",
    "redox_detector_module": "Redox Site Detector",
    "comprehensive_redox_detector": "Redox Site Detector",
    "alphafold_predictor": "AlphaFold Predictor",
    "vmd_visualizer": "Structure Viewer (VMD)",
    "interactive_structure_viewer": "Structure Viewer",
    "pdb_filter": "PDB Filter",
    "pdb_loader": "PDB Loader",
    "structure_loader": "Structure Loader",
    "amino_acid_mutator": "Amino Acid Mutator",
    "homology_searcher": "Homology Searcher",
    "biological_assembly": "Biological Assembly Generator",
    "structure_orientation": "Structure Orientator",
    "md_restraint_manager": "MD Restraint Manager",
    "molecular_dynamics_manager": "Molecular Dynamics Manager",
    "protonation_state_analyzer": "Protonation State Analyzer",
    "small_molecule_parameterizer": "Small Molecule Parameterizer",
    "forcefield_parameterizer": "Force Field Parameterizer",
    "emboss_module": "EMBOSS Analysis",
    "redoxsite_integration": "Redox Site Preparer",
}


def _humanize(stem: str) -> str:
    """Readable producer/consumer label from a source filename stem."""
    if stem in _FILE_LABELS:
        return _FILE_LABELS[stem]
    return " ".join(_ACRONYMS.get(w, w.capitalize()) for w in stem.split("_"))


# Any receiver expression ending in a `workspace` or `ws` attribute/name, so
# that aliases are covered as well as the canonical spelling: `workspace.set`,
# `self.workspace.set`, `self.processor.workspace.set`, and the bare `ws.set`
# used throughout pb_titrate. `\b(?:workspace|ws)\b` cannot match inside a
# longer identifier such as `workspace_files`, so those are still skipped.
_WS_RECV = r'(?:[A-Za-z_][\w\.]*\.)?\b(?:workspace|ws)\b'

# The wrapper helpers come in a two-argument form (key, value) and a
# three-argument form (workspace, key, value); the optional leading argument
# below accepts both.
_WS_ARG = r'(?:[A-Za-z_][\w\.]*\s*,\s*)?'

_WRITE_PATS = [
    re.compile(r'(?:update_workspace|set_in_workspace|store_in_workspace)\s*\(\s*'
               + _WS_ARG + r'[\'"]([A-Za-z_]\w+)[\'"]'),
    re.compile(_WS_RECV + r'\.set\s*\(\s*[\'"]([A-Za-z_]\w+)[\'"]'),
    re.compile(r'workspace\[\s*[\'"]([A-Za-z_]\w+)[\'"]\s*\]\s*='),
]
_READ_PATS = [
    re.compile(_WS_RECV + r'\.(?:get|has)\s*\(\s*[\'"]([A-Za-z_]\w+)[\'"]'),
    re.compile(r'(?:get_from_workspace|get_from_workspace_obj)\s*\(\s*'
               + _WS_ARG + r'[\'"]([A-Za-z_]\w+)[\'"]'),
    re.compile(r'[\'"]([A-Za-z_]\w+)[\'"]\s*in\s+(?:self\.)?workspace'),
    re.compile(r'workspace\[\s*[\'"]([A-Za-z_]\w+)[\'"]\s*\](?!\s*=)'),
]
_SCAN_IGNORE = {"debug", "key"}

# Dead-code source files (the modules in MODULES_TO_EXCLUDE — AltLoc Selector
# and Disulfide Bond Detector — whose functionality was absorbed into
# structure_completeness.py). Skip them in the scan so their keys are not
# resurrected; keys they share with live code are still picked up elsewhere.
_SCAN_EXCLUDE_FILES = {
    "altloc_selector.py", "altloc_selector_worker.py", "altloc_selector_commands.py",
    "disulfide_bond_detector.py", "disulfide_bond_detector_worker.py", "disulfide_commands.py",
    # This module documents the keys rather than participating in the
    # workspace; scanning it would turn its own examples into keys.
    "workspace_key_inventory.py",
    # Legacy metallo_prep worker and its integration shim, superseded by the
    # RedoxSite system in redoxsite_integration.py. Both import
    # proprep.metallo_prep, a package no longer in the tree, so neither can be
    # imported and nothing they appear to write is ever written. Keys they
    # share with live code (metal_sites, excluded_residues) are still picked
    # up from their live call sites.
    "metallo_worker.py", "metallo_integration.py",
}


def _strip_comments(text: str) -> str:
    """Blank out comment text so commented-out calls are not scanned.

    Without this, an illustrative comment showing a commented-out getter
    call is indistinguishable from a real call site, and invents a key
    that no code ever uses.
    Falls back to the raw text if the file does not tokenize.
    """
    try:
        out, prev_end, tokens = [], (1, 0), tokenize.generate_tokens(
            io.StringIO(text).readline)
        for tok in tokens:
            if tok.start[0] != prev_end[0]:
                prev_end = (tok.start[0], 0)
            if tok.start[1] > prev_end[1]:
                out.append(" " * (tok.start[1] - prev_end[1]))
            out.append("" if tok.type == tokenize.COMMENT else tok.string)
            if tok.end[0] != tok.start[0]:
                out.append("\n")
                prev_end = (tok.end[0], tok.end[1])
            else:
                prev_end = tok.end
        return "".join(out)
    except (tokenize.TokenError, IndentationError, SyntaxError):
        return text


def scan_source_keys():
    """Scan the proprep package source for string-literal workspace keys.

    Returns (producers, consumers): dicts mapping key -> set of human-readable
    module names (from the source filename) that set / read that key.
    """
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # proprep/
    producers: Dict[str, set] = {}
    consumers: Dict[str, set] = {}
    for dirpath, _dirs, files in os.walk(pkg_root):
        if "tests" in dirpath.split(os.sep):
            continue
        for fn in files:
            if not fn.endswith(".py") or fn in _SCAN_EXCLUDE_FILES:
                continue
            try:
                with open(os.path.join(dirpath, fn), encoding="utf-8", errors="ignore") as fh:
                    text = _strip_comments(fh.read())
            except OSError:
                continue
            mod = _humanize(os.path.splitext(fn)[0])
            for pat in _WRITE_PATS:
                for k in pat.findall(text):
                    if k not in _SCAN_IGNORE:
                        producers.setdefault(k, set()).add(mod)
            for pat in _READ_PATS:
                for k in pat.findall(text):
                    if k not in _SCAN_IGNORE:
                        consumers.setdefault(k, set()).add(mod)
    return producers, consumers


def collect_full_inventory() -> Dict[str, Dict[str, List[str]]]:
    """
    Collect workspace key inventory from all sources.

    Returns dict mapping key names to:
        {"produced_by": [...], "required_by": [...]}
    """
    failures = import_all_modules()
    if failures:
        print(f"Warning: {len(failures)} module(s) failed to import:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)

    # Collect from registered ProcessingModules
    inventory = registry.collect_workspace_key_inventory()

    # Remove excluded (dead/deprecated) modules from all entries
    for key_info in inventory.values():
        for field in ("produced_by", "required_by"):
            key_info[field] = [m for m in key_info[field] if m not in MODULES_TO_EXCLUDE]

    # Drop keys that became empty after exclusion (only had dead-code producers)
    inventory = {
        k: v for k, v in inventory.items()
        if v["produced_by"] or v["required_by"]
    }

    # Merge StructurePreprocessor keys
    for key in PREPROCESSOR_KEYS:
        if key not in inventory:
            inventory[key] = {"produced_by": [], "required_by": []}
        if PREPROCESSOR_NAME not in inventory[key]["produced_by"]:
            inventory[key]["produced_by"].append(PREPROCESSOR_NAME)

    # Completeness pass: fold in any string-literal key the declaration pass
    # missed (undeclared module outputs, plus non-module setters such as the
    # workflow state manager, settings, and navigation flags). Declared keys
    # keep their canonical attribution; only genuinely-missing keys are added,
    # with producers/consumers derived from the source scan.
    scan_prod, scan_cons = scan_source_keys()
    for key in sorted(set(scan_prod) | set(scan_cons)):
        if key in inventory:
            continue
        inventory[key] = {
            "produced_by": sorted(scan_prod.get(key, set())),
            "required_by": sorted(scan_cons.get(key, set())),
        }

    return inventory


# ── Output formatters ────────────────────────────────────────────────────────

def format_table(inventory: Dict[str, Dict[str, List[str]]]) -> None:
    """Print inventory as a Rich table."""
    try:
        from rich.console import Console
        from rich.table import Table
    except ImportError:
        print("Rich library required for table output. Use --format markdown or json.")
        sys.exit(1)

    console = Console()
    table = Table(title="ProPrep Workspace Key Inventory", show_lines=True)
    table.add_column("Key", style="cyan", no_wrap=True)
    table.add_column("Produced By", style="green")
    table.add_column("Required By", style="yellow")
    table.add_column("Type", style="grey50")
    table.add_column("Description", style="grey50")

    for key in sorted(inventory.keys()):
        info = inventory[key]
        desc = KEY_DESCRIPTIONS.get(key, {})
        table.add_row(
            key,
            ", ".join(info["produced_by"]) or "-",
            ", ".join(info["required_by"]) or "-",
            desc.get("type", ""),
            desc.get("description", ""),
        )

    console.print(table)
    console.print(f"\nTotal keys: {len(inventory)}")


def format_markdown(inventory: Dict[str, Dict[str, List[str]]],
                    abbrev: bool = False) -> None:
    """Print inventory as Markdown table (for SI).

    With abbrev=True, component names are replaced by three-letter codes and
    a legend is printed first, which keeps the table narrow enough to set in
    a portrait page.
    """
    if abbrev:
        legend = _legend_for(inventory)
        print("**Component abbreviations.** "
              + "; ".join(f"{code}, {name}" for code, name in legend)
              + ".")
        print()

    label = _abbrev if abbrev else (lambda n: n)
    print("| Key | Produced By | Required By | Type | Description |")
    print("|-----|------------|-------------|------|-------------|")
    for key in sorted(inventory.keys()):
        info = inventory[key]
        desc = KEY_DESCRIPTIONS.get(key, {})
        producers = ", ".join(label(m) for m in info["produced_by"]) or "-"
        consumers = ", ".join(label(m) for m in info["required_by"]) or "-"
        typ = desc.get("type", "")
        description = desc.get("description", "")
        print(f"| {key} | {producers} | {consumers} | {typ} | {description} |")
    print(f"\nTotal keys: {len(inventory)}")


def format_json(inventory: Dict[str, Dict[str, List[str]]]) -> None:
    """Print inventory as JSON."""
    output = {}
    for key in sorted(inventory.keys()):
        info = inventory[key]
        desc = KEY_DESCRIPTIONS.get(key, {})
        output[key] = {
            "produced_by": info["produced_by"],
            "required_by": info["required_by"],
            **desc,
        }
    print(json.dumps(output, indent=2))


# ── CLI entry point ──────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="ProPrep Workspace Key Inventory"
    )
    parser.add_argument(
        "--format",
        choices=["table", "markdown", "json"],
        default="table",
        help="Output format (default: table)",
    )
    parser.add_argument(
        "--abbrev",
        action="store_true",
        help="Abbreviate component names to three-letter codes and print a "
             "legend first (markdown format only; yields a more compact table)",
    )
    args = parser.parse_args()

    inventory = collect_full_inventory()

    if args.format == "table":
        format_table(inventory)
    elif args.format == "markdown":
        format_markdown(inventory, abbrev=args.abbrev)
    elif args.format == "json":
        format_json(inventory)


if __name__ == "__main__":
    main()
