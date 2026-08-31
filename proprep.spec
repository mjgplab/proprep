# proprep.spec
# -*- mode: python ; coding: utf-8 -*-
# Updated: 2026-08-27

import sys
import os
from pathlib import Path

# Add your source directory to Python path
sys.path.insert(0, str(Path.cwd() / "src"))

block_cipher = None

# Locate parmed data files dynamically
_parmed_data = []
try:
    import parmed
    _parmed_dir = os.path.dirname(parmed.__file__)
    _modeller_data = os.path.join(_parmed_dir, 'modeller', 'data')
    if os.path.isdir(_modeller_data):
        _parmed_data.append((_modeller_data, 'parmed/modeller/data'))
except ImportError:
    pass

# Locate MODELLER runtime data (residue topologies, parameter libraries, etc.)
# These are ~38 MB and required for MODELLER to function.
_modeller_datas = []
try:
    import glob
    _conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if _conda_prefix:
        _mod_dirs = glob.glob(os.path.join(_conda_prefix, 'lib', 'modeller-*'))
        if _mod_dirs:
            _mod_root = _mod_dirs[0]
            _mod_ver = os.path.basename(_mod_root)  # e.g. "modeller-10.8"
            # Bundle modlib/ (parameter files, topology libraries)
            _modlib = os.path.join(_mod_root, 'modlib')
            if os.path.isdir(_modlib):
                _modeller_datas.append((_modlib, f'{_mod_ver}/modlib'))
            # Bundle lib/ (shared libraries MODELLER needs at runtime)
            _modlib_lib = os.path.join(_mod_root, 'lib')
            if os.path.isdir(_modlib_lib):
                _modeller_datas.append((_modlib_lib, f'{_mod_ver}/lib'))
except Exception:
    pass

# Modules excluded from BOTH executables to reduce size.
# NOTE: scipy is deliberately NOT excluded - it is a lazy runtime dependency.
_EXCLUDES = [
    # The conda MODELLER package's config.py has the BUILD machine's license
    # key baked in. Never ship it: the runtime hook supplies a stand-in
    # (see hook-modeller-config.py / proprep.utils.modeller_license).
    'modeller.config',
        # Exclude unnecessary modules to reduce size
        'tkinter',
        'matplotlib',
        'PIL',
        'Pillow',
        'pytest',
        'IPython',
        'jupyter',
        'notebook',
        'pandas',
        'tensorflow',
        'torch',
        'cv2',
        'wx',
        'PyQt5',
        'PyQt6',
        'PySide2',
        'PySide6',
]

a = Analysis(
    ['src/proprep/main.py'],
    pathex=[str(Path.cwd() / "src")],
    binaries=[],
    datas=[
        # Force field parameter files (includes redox_centers/)
        ('src/proprep/forcefield_params', 'proprep/forcefield_params'),

        # MD workflow definitions (builtin + user)
        ('src/proprep/md_workflows', 'proprep/md_workflows'),

        # MD MDIN templates (builtin + user)
        ('src/proprep/md_templates', 'proprep/md_templates'),

        # MD prep user data (templates, workflow presets)
        ('src/proprep/md_prep/user_data', 'proprep/md_prep/user_data'),

        # MD prep annotated templates
        ('src/proprep/md_prep/templates', 'proprep/md_prep/templates'),

        # Structure viewer HTML template
        ('src/proprep/structure_prep/templates', 'proprep/structure_prep/templates'),

        # PB-titrate model compounds (prmtop/rst7 pairs read by
        # pb_titrate/model_compounds.py via Path(__file__).parent / "models")
        ('src/proprep/pb_titrate/models', 'proprep/pb_titrate/models'),

        # HPC cluster profiles (SLURM generator)
        ('src/proprep/md_prep/cluster_profiles', 'proprep/md_prep/cluster_profiles'),

    ] + _parmed_data + _modeller_datas,
    hiddenimports=[
        # =====================================================================
        # Core ProPrep modules
        # =====================================================================
        'proprep.main',

        # --- Application layer ---
        'proprep.application.pdbprocessor',
        'proprep.application.processor',
        'proprep.application.command',
        'proprep.application.command_history',
        'proprep.application.feedback_command',
        'proprep.application.menu_commands',
        'proprep.application.module_commands',
        'proprep.application.processor_command',
        'proprep.application.workspace_commands',
        'proprep.application.forcefield_explorer',

        # --- Utilities ---
        'proprep.utils.module_registry',
        'proprep.utils.workspace',
        'proprep.utils.mpsa_workspace',
        'proprep.utils.debug_utils',
        'proprep.utils.debug_workspace',
        'proprep.utils.session_recorder',
        'proprep.utils.session_editor',
        'proprep.utils.prompts',
        'proprep.utils.structure_selector',
        'proprep.utils.enhanced_menu',
        'proprep.utils.settings_manager',
        'proprep.utils.batch_processor',
        'proprep.utils.workflow_state_manager',
        'proprep.utils.workflow_checklist',
        'proprep.utils.workspace_key_inventory',
        'proprep.utils.paths',
        'proprep.utils.template_converter',

        # =====================================================================
        # Structure preparation
        # =====================================================================
        'proprep.structure_prep.pdb_loader',
        'proprep.structure_prep.structure_loader',
        'proprep.structure_prep.homology_searcher',
        'proprep.structure_prep.homology_modeler',
        'proprep.structure_prep.pdb_filter',
        'proprep.structure_prep.pdb_filter_commands',
        'proprep.structure_prep.pdb_filter_worker',
        'proprep.structure_prep.vmd_visualizer',
        'proprep.structure_prep.structure_completeness',
        'proprep.structure_prep.amino_acid_mutator',
        'proprep.structure_prep.amino_acid_mutator_commands',
        'proprep.structure_prep.disulfide_bond_detector',
        'proprep.structure_prep.disulfide_bond_detector_worker',
        'proprep.structure_prep.disulfide_commands',
        'proprep.structure_prep.protonation_state_analyzer',
        'proprep.structure_prep.protonation_commands',
        'proprep.structure_prep.protonation_worker',
        'proprep.structure_prep.altloc_selector',
        'proprep.structure_prep.altloc_selector_commands',
        'proprep.structure_prep.altloc_selector_worker',
        'proprep.structure_prep.blast_commands',
        'proprep.structure_prep.blast_worker',
        'proprep.structure_prep.chem_comp_dict_fetcher',
        'proprep.structure_prep.md_restraint_manager',
        'proprep.structure_prep.md_restraint_commands',
        'proprep.structure_prep.comprehensive_redox_detector',
        'proprep.structure_prep.redox_detector_module',
        'proprep.structure_prep.redox_site_editor',
        'proprep.structure_prep.structure_orientation',
        'proprep.structure_prep.biological_assembly',
        'proprep.structure_prep.interactive_structure_viewer',
        'proprep.structure_prep.viewer_commands',
        'proprep.structure_prep.viewer_server',
        'proprep.structure_prep.pdb_fetcher',
        'proprep.structure_prep.pdb_searcher',
        'proprep.structure_prep.uniprot_searcher',
        'proprep.structure_prep.alphafold_fetcher',
        'proprep.structure_prep.alphafold_commands',
        'proprep.structure_prep.alphafold_worker',
        'proprep.structure_prep.alphafold_predictor',
        'proprep.structure_prep.alphafill_fetcher',
        'proprep.structure_prep.structure_alignment',
        'proprep.structure_prep.structure_alignment_commands',

        # =====================================================================
        # Redox site preparation
        # =====================================================================
        'proprep.redoxsite_prep.redoxsite_integration',
        'proprep.redoxsite_prep.metallo_integration',
        'proprep.redoxsite_prep.redox_commands',
        'proprep.redoxsite_prep.metallo_commands',
        'proprep.redoxsite_prep.metallo_worker',
        'proprep.redoxsite_prep.restraint_mask_generator',

        # --- Transformation framework ---
        'proprep.redoxsite_prep.transformation.redox_transformer_framework',
        'proprep.redoxsite_prep.transformation.redox_transformation_manager',
        'proprep.redoxsite_prep.transformation.redox_id_mapper',
        'proprep.redoxsite_prep.transformation.transformer_creator',

        # --- Transformer creator package ---
        'proprep.redoxsite_prep.transformation.transformer_creator_package.creator',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.code_generator',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.data_models',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.site_analyzer',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.template_extractor',

        # --- Transformer creator v3 ---
        'proprep.redoxsite_prep.transformation.transformer_creator_v3',
        'proprep.redoxsite_prep.transformation.transformer_creator_v3.v3_creator',
        'proprep.redoxsite_prep.transformation.transformer_creator_v3.code_generator',
        'proprep.redoxsite_prep.transformation.transformer_creator_v3.command_parser',
        'proprep.redoxsite_prep.transformation.transformer_creator_v3.constraint_resolver',
        'proprep.redoxsite_prep.transformation.transformer_creator_v3.data_models',
        'proprep.redoxsite_prep.transformation.transformer_creator_v3.phrase_parser',

        # --- Transformer creator collectors ---
        'proprep.redoxsite_prep.transformation.transformer_creator_package.collectors.component_matcher_collector',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.collectors.parameter_collector',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.collectors.site_requirements_collector',
        'proprep.redoxsite_prep.transformation.transformer_creator_package.collectors.transformation_collector',

        # --- Built-in transformers ---
        'proprep.redoxsite_prep.transformation.transformers.bis_his_c_type_heme',
        'proprep.redoxsite_prep.transformation.transformers.disulfide',
        'proprep.redoxsite_prep.transformation.transformers.no_transformation',

        # =====================================================================
        # EMBOSS module
        # =====================================================================
        'proprep.emboss.emboss_module',
        'proprep.emboss.emboss_commands',
        'proprep.emboss.emboss_core',
        'proprep.emboss.sequence_input',

        # =====================================================================
        # Force field preparation
        # =====================================================================
        'proprep.forcefield_prep.forcefield_parameterizer',
        'proprep.forcefield_prep.forcefield_commands',
        'proprep.forcefield_prep.forcefield_worker',
        'proprep.forcefield_prep.forcefield_data',
        'proprep.forcefield_prep.forcefield_data_collector',
        'proprep.forcefield_prep.forcefield_loader',
        'proprep.forcefield_prep.ff_parsers',
        'proprep.forcefield_prep.ff_types',
        'proprep.forcefield_prep.atom_type_counter',
        'proprep.forcefield_prep.modified_amino_acid_parameterizer',
        'proprep.forcefield_prep.metal_site_parameterizer',
        'proprep.forcefield_prep.small_molecule_parameterizer',
        'proprep.forcefield_prep.small_molecule_commands',
        'proprep.forcefield_prep.model_builder',
        'proprep.forcefield_prep.nonstandard_residue_handler',
        'proprep.forcefield_prep.structure_preprocessor',
        'proprep.forcefield_prep.paramfit_refinement',
        'proprep.forcefield_prep.pes_scan_refinement',
        'proprep.forcefield_prep.seminario_refinement',
        'proprep.forcefield_prep.pdb_writer',
        'proprep.forcefield_prep.xyz_writer',
        'proprep.forcefield_prep.terminal_classifier',

        # --- MCPB submodules ---
        'proprep.forcefield_prep.mcpb.antechamber_runner',
        'proprep.forcefield_prep.mcpb.atom_typer',
        'proprep.forcefield_prep.mcpb.esp_extractor',
        'proprep.forcefield_prep.mcpb.ff_library',
        'proprep.forcefield_prep.mcpb.ff_parameter_reader',
        'proprep.forcefield_prep.mcpb.fingerprint_generator',
        'proprep.forcefield_prep.mcpb.frcmod_builder',
        'proprep.forcefield_prep.mcpb.hessian_parser',
        'proprep.forcefield_prep.mcpb.ligand_grouping',
        'proprep.forcefield_prep.mcpb.metal_ion_database',
        'proprep.forcefield_prep.mcpb.mol2_writer',
        'proprep.forcefield_prep.mcpb.pre_frcmod_generator',
        'proprep.forcefield_prep.mcpb.qm_interface',
        'proprep.forcefield_prep.mcpb.resp_input_generator',
        'proprep.forcefield_prep.mcpb.resp_runner',
        'proprep.forcefield_prep.mcpb.seminario',
        'proprep.forcefield_prep.mcpb.terminal_classifier',
        'proprep.forcefield_prep.mcpb.global_atom_registry',
        'proprep.forcefield_prep.mcpb.integration_utils',

        # --- Force field params loader ---
        'proprep.forcefield_params.loader',

        # =====================================================================
        # MD preparation
        # =====================================================================
        'proprep.md_prep.molecular_dynamics_manager',
        'proprep.md_prep.molecular_dynamics_commands',
        'proprep.md_prep.layout_helpers',
        'proprep.md_prep.workflow_centric_step1',
        'proprep.md_prep.workflow_source_selectors',
        'proprep.md_prep.workflow_editor',
        'proprep.md_prep.workflow_loader',
        'proprep.md_prep.workflow_presets',
        'proprep.md_prep.workflow_commands',
        'proprep.md_prep.custom_workflow_creator',
        'proprep.md_prep.user_data_manager',
        'proprep.md_prep.template_workflow_manager',
        'proprep.md_prep.amber_input_generator',
        'proprep.md_prep.amber_input_generator_v2',
        'proprep.md_prep.amber_annotated_templates',
        'proprep.md_prep.amber_parameter_database',
        'proprep.md_prep.amber_parameter_display',
        'proprep.md_prep.amber_commands',
        'proprep.md_prep.amber_controller',
        'proprep.md_prep.amber_template_system',
        'proprep.md_prep.amber_wizard',
        'proprep.md_prep.amber_workflow_components',
        'proprep.md_prep.amber_workflow_manager',
        'proprep.md_prep.slurm_generator',
        'proprep.md_prep.trajectory_analyzer',
        'proprep.md_prep.step0_prototype',

        # --- AMBER MDIN keyword/parser subpackage ---
        'proprep.md_prep.amber._amber_keywords',
        'proprep.md_prep.amber.exceptions',
        'proprep.md_prep.amber.keyword_db',
        'proprep.md_prep.amber.mdin_builder',
        'proprep.md_prep.amber.mdin_parser',
        'proprep.md_prep.amber.mdin_validator',
        'proprep.md_prep.amber.mdin_writer',

        # =====================================================================
        # Membrane preparation
        # =====================================================================
        'proprep.membrane_prep.membrane_builder',
        'proprep.membrane_prep.membrane_config',
        'proprep.membrane_prep.membrane_presets',
        'proprep.membrane_prep.lipid_library',
        'proprep.membrane_prep.packmol_runner',

        # =====================================================================
        # TLeap preparation
        # =====================================================================
        'proprep.tleap_prep.tleap_input_generator',
        'proprep.tleap_prep.tleap_commands',
        'proprep.tleap_prep.bond_management_commands',
        'proprep.tleap_prep.bond_management_base',
        'proprep.tleap_prep.pdb_molecule_configurator',

        # =====================================================================
        # ONIOM preparation
        # =====================================================================
        'proprep.oniom_prep.oniom_commands',
        'proprep.oniom_prep.oniom_preparator',
        'proprep.oniom_prep.oniom_qmmm_preparator',
        'proprep.oniom_prep.oniom_atom_typer',
        'proprep.oniom_prep.oniom_writer',
        'proprep.oniom_prep.oniom_validator',
        'proprep.oniom_prep.connectivity_builder',
        'proprep.oniom_prep.layer_classifier',
        'proprep.oniom_prep.link_atom_placer',
        'proprep.oniom_prep.amber_lib_parser',
        'proprep.oniom_prep.data_structures',
        'proprep.oniom_prep.structure_utils',

        # =====================================================================
        # ORCA preparation
        # =====================================================================
        'proprep.orca_prep.orca_data_structures',
        'proprep.orca_prep.orca_ff_converter',
        'proprep.orca_prep.orca_qmmm_preparator',
        'proprep.orca_prep.orca_writer',

        # =====================================================================
        # QM/MM preparation (imported by name in pdbprocessor.py)
        # =====================================================================
        'proprep.qmmm_prep.qmmm_preparator',
        'proprep.qmmm_prep.frame_extractor',

        # =====================================================================
        # Third-party dependencies
        # =====================================================================

        # BioPython
        'Bio',
        'Bio.PDB',
        'Bio.PDB.PDBParser',
        'Bio.PDB.PDBIO',
        'Bio.SeqUtils',
        'Bio.Seq',

        # Rich
        'rich',
        'rich.console',
        'rich.table',
        'rich.panel',
        'rich.prompt',
        'rich.status',
        'rich.progress',
        'rich.markdown',
        'rich.syntax',

        # Other
        'requests',
        'freesasa',
        'yaml',
        'pyyaml',
        'numpy',
        'networkx',
        # scipy is imported lazily (pb_titrate cKDTree/curve_fit,
        # trajectory_analyzer cdist, molecular_dynamics_manager ndimage)
        'scipy',
        'scipy.spatial',
        'scipy.optimize',
        'scipy.ndimage',
        'sympy',
        'pydantic',
        'pydantic_core',
        'urllib3',
        'certifi',
        'charset_normalizer',
        'idna',
        'typing_extensions',
        'json',
        'dataclasses',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=['hook-modeller-config.py'],
    excludes=_EXCLUDES,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# ---------------------------------------------------------------------------
# proprep-web: the browser shell (entry point proprep.web.__main__:main).
# It only imports proprep.web.* + fastapi/uvicorn/httpx and execs the sibling
# ``proprep`` binary inside a PTY (see web/pty_session.py, frozen branch), so
# it gets its own small Analysis rather than a copy of the full one.
# ---------------------------------------------------------------------------
a_web = Analysis(
    ['src/proprep/web/__main__.py'],
    pathex=[str(Path.cwd() / "src")],
    binaries=[],
    datas=[
        # Served by web/server.py via Path(__file__).parent / "static"
        ('src/proprep/web/static', 'proprep/web/static'),
    ],
    hiddenimports=[
        'proprep.web',
        'proprep.web.server',
        'proprep.web.pty_session',
        'fastapi',
        'fastapi.responses',
        'fastapi.staticfiles',
        'uvicorn',
        'httpx',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=_EXCLUDES,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz_web = PYZ(a_web.pure, a_web.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='proprep',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)

exe_web = EXE(
    pyz_web,
    a_web.scripts,
    a_web.binaries,
    a_web.zipfiles,
    a_web.datas,
    [],
    name='proprep-web',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
)
