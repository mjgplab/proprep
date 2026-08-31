"""Batch microstate generation stamps fixed_E for cofactors with a
fixed_E/constant_E fork and threads it into the per-state parameters and
transformer_info, so the Topology Generator never offers the constant-E HEH
library for a microstate whose PDB carries HCO/HCR. Regression for the
HEH-vs-HCO mismatch that put every microstate at the same charge."""
import io

from rich.console import Console

from proprep.redoxsite_prep.transformation.redox_transformation_manager import RedoxTransformationManager
from proprep.redoxsite_prep.transformation.redox_transformer_framework import redox_transformer_registry
from proprep.tleap_prep.tleap_input_generator import TLeapInputGenerator


class _Center:
    chain, resname, resid = 'A', 'HEC', 92


class _Site:
    def __init__(self, site_id):
        self.site_id = site_id
        self.centers = [_Center()]


def _manager(assignments):
    m = object.__new__(RedoxTransformationManager)
    m.console = Console(file=io.StringIO(), force_terminal=False)
    m.transformation_parameters = {}
    m.site_state_protonation = {}
    m.site_redox_treatment = {}
    m.redox_sites = [_Site(s) for s in assignments]
    m.site_transformer_assignments = dict(assignments)
    return m


def test_batch_stamps_fixed_E_only_for_fork_transformers():
    m = _manager({'site_1': 'heme_bis_his_c_type', 'site_2': 'no_transformation'})
    m._stamp_batch_redox_treatment()
    assert m.site_redox_treatment == {'site_1': 'fixed_E'}


def test_fixed_E_reaches_params_transformer_info_and_names():
    m = _manager({'site_1': 'heme_bis_his_c_type'})
    m._stamp_batch_redox_treatment()
    params = m._resolve_state_params('site_1', 'oxidized', 'low_spin')
    assert params['redox_treatment'] == 'fixed_E'

    tclass = redox_transformer_registry.get_transformer('heme_bis_his_c_type')
    assert tclass.get_parameter_mappings({**params, 'ph_treatment': 'constant_pH'})['heme_name'] == 'HCO'
    assert tclass.get_parameter_mappings({'redox_state': 'reduced', 'spin_state': 'low_spin',
                                          'redox_treatment': 'fixed_E', 'ph_treatment': 'constant_pH'})['heme_name'] == 'HCR'

    info = m._build_batch_transformer_info([
        {'site_1': {'redox_state': 'oxidized', 'spin_state': 'low_spin'}},
        {'site_1': {'redox_state': 'reduced', 'spin_state': 'low_spin'}},
    ])
    assert {e['redox_treatment'] for e in info} == {'fixed_E'}
    assert all('HEH' not in (e['forcefield_set'] or '') for e in info)


def test_tleap_unknown_residue_lines_are_collected():
    lines = [
        "Loading PDB file: ./transformed_microstate_001.pdb",
        "Unknown residue: HCO   number: 91   type: Terminal/last",
        "..relaxing end constraints to try for a dbase match",
        "Unknown residue: HCO   number: 94   type: Nonterminal",
        "Unknown residue: PRN   number: 96   type: Nonterminal",
        "Total unperturbed charge:  -15.000000",
    ]
    assert TLeapInputGenerator._unknown_residues_in_tleap_output(lines) == ['HCO', 'PRN']
    assert TLeapInputGenerator._unknown_residues_in_tleap_output(["Total unperturbed charge: -11.0"]) == []
