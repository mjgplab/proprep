"""
Docking preparation for AutoDock Vina and AutoDock4 scoring.

Deliberately NOT registered as a ``ProcessingModule`` yet -- importing this
package has no effect on the menus. Wiring it in means adding it to the three
registration lists (``application/pdbprocessor.py`` import block,
``menu_commands.py`` main-menu groups, ``menu_commands.py`` STAGE_TOOLS).

Docking needs no force-field parameters. Ligands and receptors are prepared
with Meeko, the AutoDock developers' own library, from chemistry ProPrep
checks and the user can see and edit:

* ``ccd_chemistry``      -- CCD entries as RDKit molecules, with every valence
                            problem reported per atom (never repaired silently)
* ``chemistry_edits``    -- named, replayable edits to bond orders, charges and H
* ``ligand_sources``     -- ligands from SMILES, SDF/mol2, or a HETATM residue
* ``ligand_prep``        -- ligand PDBQT with exact per-bond torsion choices
* ``receptor_decisions`` -- proposals: pdb2pqr/PROPKA protonation, metal
                            binding, termini from REMARK 465, LINK records
* ``receptor_templates`` -- Meeko templates for cofactors, ions and links
* ``receptor_prep``      -- receptor PDBQT from explicit decisions
* ``docking_run``        -- search box, settings, one Vina/Vinardo/AD4 run
* ``docking_results``    -- poses back to named molecules, RMSD, SDF
* ``dependencies``       -- which of RDKit, gemmi, Meeko, Vina and autogrid4
                            are installed, and what each is needed for

Nothing is imported here, so the package imports even where RDKit or Meeko
are missing; ``dependencies`` says what is.

The design evidence (real-structure probes of Meeko, Vina and AutoGrid) is
in ``tools/docking_phase0``.
"""
