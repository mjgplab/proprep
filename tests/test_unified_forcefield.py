"""
Tests for Unified Force Field System

Tests the unified force field parsing and data management:
- ff_parsers.py: Unified parsing functions
- forcefield_data.py: ForceFieldData class

These tests require AMBERHOME to be set to a valid AMBER installation.

Author: ProPrep Development Team
"""

import os
import pytest
from pathlib import Path

# Skip all tests if AMBERHOME not set
pytestmark = pytest.mark.skipif(
    not os.environ.get('AMBERHOME'),
    reason="AMBERHOME not set"
)


@pytest.fixture
def amberhome():
    """Get AMBERHOME path."""
    return Path(os.environ.get('AMBERHOME', ''))


@pytest.fixture
def ff_data():
    """Create ForceFieldData with common force fields loaded."""
    from proprep.forcefield_prep.forcefield_data import ForceFieldData
    from rich.console import Console

    data = ForceFieldData(console=Console())
    data.load_leaprc('leaprc.protein.ff14SB')
    data.load_leaprc('leaprc.water.tip3p')
    data.load_leaprc('leaprc.gaff2')
    return data


# =============================================================================
# Parser Tests
# =============================================================================

class TestParsers:
    """Tests for ff_parsers.py parsing functions."""

    def test_parse_leaprc(self, amberhome):
        """Test parse_leaprc extracts correct directives."""
        from proprep.forcefield_prep.ff_parsers import parse_leaprc

        leaprc_path = amberhome / 'dat' / 'leap' / 'cmd' / 'leaprc.protein.ff14SB'
        if not leaprc_path.exists():
            pytest.skip(f"leaprc.protein.ff14SB not found at {leaprc_path}")

        contents = parse_leaprc(leaprc_path)

        assert len(contents.lib_files) > 0, "Should find lib files"
        assert 'amino12.lib' in contents.lib_files, "Should include amino12.lib"

    def test_parse_lib_file(self, amberhome):
        """Test parse_lib_file extracts atom definitions."""
        from proprep.forcefield_prep.ff_parsers import parse_lib_file

        lib_path = amberhome / 'dat' / 'leap' / 'lib' / 'amino12.lib'
        if not lib_path.exists():
            pytest.skip(f"amino12.lib not found at {lib_path}")

        atoms = parse_lib_file(lib_path, 'amino12.lib')

        assert len(atoms) > 0, "Should extract atom definitions"

        # Check for ALA residue
        ala_atoms = [a for a in atoms if a.residue_name == 'ALA']
        assert len(ala_atoms) > 0, "Should include ALA atoms"

        # Check ALA CA atom
        ala_ca = [a for a in ala_atoms if a.atom_name == 'CA']
        assert len(ala_ca) == 1, "Should have one ALA CA"
        assert ala_ca[0].atom_type == 'CX', "ALA CA should be type CX"

    def test_parse_parm_dat(self, amberhome):
        """Test parse_parm_dat extracts parameters."""
        from proprep.forcefield_prep.ff_parsers import parse_parm_dat

        parm_path = amberhome / 'dat' / 'leap' / 'parm' / 'parm10.dat'
        if not parm_path.exists():
            pytest.skip(f"parm10.dat not found at {parm_path}")

        contents = parse_parm_dat(parm_path, 'parm10.dat')

        assert len(contents.mass_parameters) > 0, "Should extract mass parameters"
        assert len(contents.bond_parameters) > 0, "Should extract bond parameters"
        assert len(contents.angle_parameters) > 0, "Should extract angle parameters"
        assert len(contents.nonbonded_parameters) > 0, "Should extract VDW parameters"

        # Check a known bond
        ct_n = contents.bond_parameters.get('CT-N')
        assert ct_n is not None, "Should have CT-N bond"
        assert ct_n.force_constant > 0, "Bond should have positive force constant"

    def test_parse_frcmod(self, amberhome):
        """Test parse_frcmod extracts parameter modifications."""
        from proprep.forcefield_prep.ff_parsers import parse_frcmod

        frcmod_path = amberhome / 'dat' / 'leap' / 'parm' / 'frcmod.ff14SB'
        if not frcmod_path.exists():
            pytest.skip(f"frcmod.ff14SB not found at {frcmod_path}")

        contents = parse_frcmod(frcmod_path, 'frcmod.ff14SB')

        # ff14SB frcmod should have some modifications
        total_params = (
            len(contents.bond_parameters) +
            len(contents.angle_parameters) +
            len(contents.dihedral_parameters)
        )
        assert total_params > 0, "Should extract some parameters from frcmod"


# =============================================================================
# ForceFieldData Tests
# =============================================================================

class TestForceFieldData:
    """Tests for ForceFieldData class."""

    def test_load_leaprc(self):
        """Test loading a leaprc file."""
        from proprep.forcefield_prep.forcefield_data import ForceFieldData
        from rich.console import Console

        ff_data = ForceFieldData(console=Console())
        atoms_loaded = ff_data.load_leaprc('leaprc.protein.ff14SB')

        assert atoms_loaded > 0, "Should load atom definitions"
        assert len(ff_data.loaded_leaprcs) >= 1, "Should track loaded leaprcs"
        assert len(ff_data.loaded_libs) > 0, "Should load lib files"

    def test_get_atom_type(self, ff_data):
        """Test atom type lookup."""
        # Standard amino acid
        assert ff_data.get_atom_type('ALA', 'CA') == 'CX'
        assert ff_data.get_atom_type('ALA', 'N') == 'N'
        assert ff_data.get_atom_type('ALA', 'C') == 'C'

        # Unknown residue
        assert ff_data.get_atom_type('XYZ', 'CA') is None

    def test_get_atom_charge(self, ff_data):
        """Test atom charge lookup."""
        charge = ff_data.get_atom_charge('ALA', 'CA')
        assert charge is not None
        assert isinstance(charge, float)

    def test_has_residue(self, ff_data):
        """Test residue existence check."""
        assert ff_data.has_residue('ALA') is True
        assert ff_data.has_residue('GLY') is True
        # Note: AMBER uses HID/HIE/HIP for histidine variants, not HIS
        assert ff_data.has_residue('HID') is True
        assert ff_data.has_residue('XYZ') is False

    def test_get_bond_parameter(self, ff_data):
        """Test bond parameter lookup."""
        bond = ff_data.get_bond_parameter('CT', 'N')
        assert bond is not None
        assert bond.force_constant > 0
        assert bond.eq_length > 0

        # Test reverse ordering
        bond_rev = ff_data.get_bond_parameter('N', 'CT')
        assert bond_rev is not None
        assert bond_rev.force_constant == bond.force_constant

    def test_get_angle_parameter(self, ff_data):
        """Test angle parameter lookup."""
        angle = ff_data.get_angle_parameter('N', 'CT', 'C')
        assert angle is not None
        assert angle.force_constant > 0
        assert angle.eq_angle > 0

    def test_get_mass_parameter(self, ff_data):
        """Test mass parameter lookup."""
        mass = ff_data.get_mass_parameter('C')
        assert mass is not None
        assert mass.mass == pytest.approx(12.01, rel=0.01)

    def test_get_nonbonded_parameter(self, ff_data):
        """Test VDW parameter lookup."""
        vdw = ff_data.get_nonbonded_parameter('CT')
        assert vdw is not None
        assert vdw.radius > 0
        assert vdw.well_depth > 0

    def test_add_residue_alias(self, ff_data):
        """Test adding a residue alias."""
        # Add alias
        copied = ff_data.add_residue_alias('MYRES', 'ALA')
        assert copied > 0, "Should copy atoms"

        # Verify alias works
        assert ff_data.has_residue('MYRES') is True
        assert ff_data.get_atom_type('MYRES', 'CA') == 'CX'

    def test_get_available_residues(self, ff_data):
        """Test getting list of available residues."""
        residues = ff_data.get_available_residues()
        assert len(residues) > 0
        assert 'ALA' in residues
        assert 'GLY' in residues

    def test_get_statistics(self, ff_data):
        """Test statistics generation."""
        stats = ff_data.get_statistics()

        assert 'n_atom_definitions' in stats
        assert 'n_residues' in stats
        assert 'n_bond_params' in stats
        assert stats['n_atom_definitions'] > 0
        assert stats['n_residues'] > 0


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests with structure_preprocessor."""

    def test_select_force_field_non_interactive(self):
        """Test _select_force_field in non-interactive mode."""
        from proprep.forcefield_prep.structure_preprocessor import StructurePreprocessor
        from rich.console import Console

        preprocessor = StructurePreprocessor(processor=None, console=Console())
        ff_data = preprocessor._select_force_field(interactive=False)

        # Verify ForceFieldData was returned
        assert ff_data is not None
        assert len(ff_data.atom_definitions) > 0

        # Verify lookups work
        assert ff_data.get_atom_type('ALA', 'CA') == 'CX'
        assert ff_data.get_bond_parameter('CT', 'N') is not None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
