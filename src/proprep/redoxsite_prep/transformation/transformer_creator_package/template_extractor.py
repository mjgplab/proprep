#!/usr/bin/env python3
"""
Template Extractor

Extracts transformer specifications from existing registered transformers
to use as starting templates for new transformer creation.
"""

import logging
from typing import Dict, List, Optional, Any
from pathlib import Path

from .data_models import (
    TransformerSpec,
    SiteRequirements,
    CenterRequirements,
    AtomRequirements,
    ResidueRequirement,
    BondRequirements,
    BondGroup,
    ComponentMatching,
    ComponentMatchingRule,
    TransformationSequence,
    Transformation,
    TransformationSelector,
    TransformationAction,
    Parameters,
    ParameterDefinition,
    ParameterType,
    CenterTypeEnum
)

logger = logging.getLogger(__name__)


class TemplateExtractor:
    """Extract transformer specifications from existing transformers"""

    def __init__(self):
        self.available_transformers = {}
        self._load_registry()

    def _load_registry(self):
        """Load all registered transformers"""
        try:
            from proprep.redoxsite_prep.transformation.redox_transformer_framework import (
                redox_transformer_registry
            )

            all_transformers = redox_transformer_registry.get_all_transformers()
            self.available_transformers = all_transformers

            logger.info(f"Loaded {len(all_transformers)} registered transformers")

        except Exception as e:
            logger.warning(f"Could not load transformer registry: {e}")
            self.available_transformers = {}

    def list_available_templates(self) -> List[Dict[str, str]]:
        """Get list of available templates with basic info"""
        templates = []

        for name, transformer_class in self.available_transformers.items():
            templates.append({
                'name': name,
                'class_name': transformer_class.__name__,
                'description': getattr(
                    transformer_class,
                    'DESCRIPTION',
                    transformer_class.__doc__ or 'No description'
                )
            })

        return templates

    def extract_template(self, transformer_name: str) -> Optional[TransformerSpec]:
        """
        Extract a complete TransformerSpec from an existing transformer

        Args:
            transformer_name: Name of registered transformer

        Returns:
            TransformerSpec if successful, None otherwise
        """
        if transformer_name not in self.available_transformers:
            logger.error(f"Transformer '{transformer_name}' not found in registry")
            return None

        transformer_class = self.available_transformers[transformer_name]

        try:
            # Extract all components
            site_reqs = self._extract_site_requirements(transformer_class)
            component_matching = self._extract_component_matching(transformer_class)
            transformations = self._extract_transformations(transformer_class)
            parameters = self._extract_parameters(transformer_class)

            # Get basic metadata
            name = getattr(transformer_class, 'TRANSFORMER_NAME', transformer_name)
            description = getattr(
                transformer_class,
                'DESCRIPTION',
                transformer_class.__doc__ or ''
            )
            supported_types = getattr(
                transformer_class,
                'SUPPORTED_SITE_TYPES',
                [name]
            )
            forcefield_path = getattr(transformer_class, 'FORCEFIELD_PATH', None)
            required_count = transformer_class.get_required_residue_count() if hasattr(
                transformer_class, 'get_required_residue_count'
            ) else 1

            # Create spec
            spec = TransformerSpec(
                name=name,
                class_name=transformer_class.__name__,
                description=description,
                supported_site_types=supported_types,
                site_requirements=site_reqs,
                component_matching=component_matching,
                transformation_sequence=transformations,
                parameters=parameters,
                forcefield_path=forcefield_path,
                template_source=transformer_name,
                required_residue_count=required_count
            )

            logger.info(f"Successfully extracted template from '{transformer_name}'")
            return spec

        except Exception as e:
            logger.error(f"Failed to extract template from '{transformer_name}': {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None

    def _extract_site_requirements(self, transformer_class) -> SiteRequirements:
        """Extract site requirements from get_site_requirements()"""
        try:
            reqs_dict = transformer_class.get_site_requirements()

            # Extract center requirements
            center_data = reqs_dict.get('centers', {})
            centers = CenterRequirements(
                required_count=center_data.get('required_count', 1),
                center_types=[
                    CenterTypeEnum(ct.value if hasattr(ct, 'value') else ct)
                    for ct in center_data.get('center_types', [])
                ],
                elements=center_data.get('elements')
            )

            # Extract atom requirements
            atoms = None
            if 'atoms' in reqs_dict:
                atom_data = reqs_dict['atoms']
                required_residues = {}

                for resname, counts in atom_data.get('required_residues', {}).items():
                    required_residues[resname] = ResidueRequirement(
                        min_count=counts.get('min_count', 0),
                        max_count=counts.get('max_count')
                    )

                atoms = AtomRequirements(
                    required_residues=required_residues,
                    alternative_groups=atom_data.get('alternative_groups')
                )

            # Extract bond requirements (simplified - this is complex)
            bonds = None
            if 'bonds' in reqs_dict:
                bond_data = reqs_dict['bonds']
                bond_groups = []

                for group in bond_data.get('required_bond_groups', []):
                    bond_groups.append(BondGroup(
                        description=group.get('description', ''),
                        min_count=group.get('min_count', 1),
                        bond_types=group.get('bond_types', {})
                    ))

                bonds = BondRequirements(
                    required_bond_groups=bond_groups,
                    require_one_group=bond_data.get('require_one_group', False)
                )

            return SiteRequirements(
                centers=centers,
                atoms=atoms,
                bonds=bonds
            )

        except Exception as e:
            logger.warning(f"Could not fully extract site requirements: {e}")
            # Return minimal valid requirements
            return SiteRequirements(
                centers=CenterRequirements(
                    required_count=1,
                    center_types=[CenterTypeEnum.METAL_ION]
                )
            )

    def _extract_component_matching(self, transformer_class) -> ComponentMatching:
        """Extract component matching logic from match_components()"""
        # This is challenging because match_components is imperative code
        # For now, return empty matching rules
        # User will need to define these interactively
        logger.info("Component matching rules must be defined interactively")
        return ComponentMatching(rules={})

    def _extract_transformations(self, transformer_class) -> TransformationSequence:
        """Extract transformations from get_transformation_sequence()"""
        try:
            # Call with dummy components and parameters to get structure
            dummy_components = {}
            dummy_parameters = {}

            trans_list = transformer_class.get_transformation_sequence(
                dummy_components,
                dummy_parameters
            )

            transformations = []
            for t in trans_list:
                transformations.append(Transformation(
                    id=t.get('id', ''),
                    description=t.get('description', ''),
                    selector=TransformationSelector(**t.get('selector', {})),
                    action=TransformationAction(**t.get('action', {}))
                ))

            return TransformationSequence(transformations=transformations)

        except Exception as e:
            logger.warning(f"Could not extract transformations: {e}")
            return TransformationSequence(transformations=[])

    def _extract_parameters(self, transformer_class) -> Optional[Parameters]:
        """Extract parameters from get_parameter_definitions() and get_parameter_mappings()"""
        try:
            # Get parameter definitions
            param_defs_dict = {}
            if hasattr(transformer_class, 'get_parameter_definitions'):
                defs = transformer_class.get_parameter_definitions()

                for param_name, param_def in defs.items():
                    param_type = param_def.get('type', 'choice')
                    param_defs_dict[param_name] = ParameterDefinition(
                        description=param_def.get('description', ''),
                        type=ParameterType(param_type),
                        options=param_def.get('options'),
                        default=param_def.get('default'),
                        value=param_def.get('value'),
                        note=param_def.get('note')
                    )

            # Get parameter mappings
            mappings_dict = {}
            if hasattr(transformer_class, 'get_parameter_mappings'):
                # This requires dummy parameters - skip for now
                pass

            if param_defs_dict or mappings_dict:
                return Parameters(
                    definitions=param_defs_dict,
                    mappings=mappings_dict
                )

            return None

        except Exception as e:
            logger.warning(f"Could not extract parameters: {e}")
            return None
