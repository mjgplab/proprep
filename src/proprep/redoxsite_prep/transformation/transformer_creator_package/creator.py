#!/usr/bin/env python3
"""
Transformer Creator - Main Orchestrator

Simple orchestrator that ties together all collectors and generates the transformer.
"""

import logging
from pathlib import Path
from typing import Optional

from .data_models import (
    TransformerSpec,
    create_empty_spec,
    ComponentMatching,
    TransformationSequence,
    SiteRequirements
)
from .template_extractor import TemplateExtractor
from .site_analyzer import SiteAnalyzer
from .collectors.site_requirements_collector import SiteRequirementsCollector
from .collectors.component_matcher_collector import ComponentMatcherCollector
from .collectors.transformation_collector import TransformationCollector
from .collectors.parameter_collector import ParameterCollector
from .code_generator import CodeGenerator

logger = logging.getLogger(__name__)


class TransformerCreator:
    """Main orchestrator for creating transformers"""

    def __init__(self, output_dir: Optional[Path] = None, workspace=None):
        """
        Initialize creator

        Args:
            output_dir: Directory to save generated transformers
                       (default: transformers/ subdirectory)
            workspace: ProPrep workspace (optional, enables RedoxSite analysis)
        """
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            # Default to transformers subdirectory in user's home
            self.output_dir = Path.home() / ".proprep" / "transformers"

        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Store workspace for RedoxSite analysis
        self.workspace = workspace

        # Initialize components
        self.template_extractor = TemplateExtractor()
        self.site_analyzer = SiteAnalyzer()
        self.site_req_collector = SiteRequirementsCollector()
        self.component_matcher_collector = ComponentMatcherCollector()
        self.transformation_collector = TransformationCollector()
        self.parameter_collector = ParameterCollector()

    def create_transformer(self) -> Optional[Path]:
        """
        Main method to create a transformer interactively

        Returns:
            Path to generated transformer file, or None if cancelled
        """
        print("\n" + "="*70)
        print("TRANSFORMER CREATOR v2")
        print("="*70)
        print("\nCreate a custom RedoxSite transformer with guided steps.\n")

        # Step 0: Template selection (optional)
        spec = self._template_selection_step()

        if spec is None:
            # User wants to start from scratch
            name = input("\nTransformer name (snake_case, e.g., my_cofactor): ").strip()
            description = input("Description: ").strip()

            if not name:
                print("Name required. Exiting.")
                return None

            spec = create_empty_spec(name, description)

        # Update basic info if needed
        print(f"\nTransformer name: {spec.name}")
        print(f"Class name: {spec.class_name}")
        print(f"Description: {spec.description}")

        # Step 1: Site Requirements
        # Skip if requirements were extracted from workspace (already complete)
        if not getattr(spec, '_from_workspace', False):
            spec.site_requirements = self.site_req_collector.collect(spec.site_requirements)
        else:
            print("\n" + "="*70)
            print("STEP 1: SITE REQUIREMENTS")
            print("="*70)
            print("\n[grey50]Requirements already extracted from workspace RedoxSite.[/grey50]")
            print("Do you want to modify these requirements?")
            modify = input("Edit requirements? [y/N]: ").strip().lower()
            if modify == 'y':
                spec.site_requirements = self.site_req_collector.collect(spec.site_requirements)
            else:
                print("✓ Using extracted requirements as-is")

        # Step 2: Component Matching
        # Always offer component matching, even if starting fresh
        component_names = self._extract_component_names(spec)

        # Ask user if they want to define components
        print("\n" + "="*70)
        print("STEP 2: COMPONENT MATCHING")
        print("="*70)
        print("\nComponents identify specific residues/atoms beyond the center.")
        print("(e.g., 'axial_his_1', 'axial_his_2', 'cys_ligand_1', etc.)\n")

        if component_names:
            print(f"Found {len(component_names)} component(s) from existing data: {', '.join(component_names)}")
            define = input("Edit component matching rules? [Y/n]: ").strip().lower()
            if define != 'n':
                spec.component_matching = self.component_matcher_collector.collect(
                    component_names,
                    spec.component_matching
                )
        else:
            print("No components defined yet.")
            define = input("Define component matching rules? [Y/n]: ").strip().lower()
            if define != 'n':
                # Ask how many components
                print("\nHow many components do you need to identify?")
                print("(e.g., 4 for four CYS ligands, 2 for two HIS ligands, etc.)")
                count_str = input("Number of components [0]: ").strip()
                try:
                    count = int(count_str) if count_str else 0
                    if count > 0:
                        # Collect component names first
                        print("\nEnter component names (snake_case):")
                        for i in range(count):
                            comp_name = input(f"  Component {i+1} name: ").strip()
                            if comp_name:
                                component_names.append(comp_name)

                        if component_names:
                            spec.component_matching = self.component_matcher_collector.collect(
                                component_names,
                                spec.component_matching
                            )
                except ValueError:
                    print("Invalid number, skipping component matching")

        # Step 3: Transformation Sequence
        all_components = ['center'] + component_names
        spec.transformation_sequence = self.transformation_collector.collect(
            all_components,
            spec.transformation_sequence
        )

        # Step 4: Parameters (optional)
        spec.parameters = self.parameter_collector.collect(spec.parameters)

        # Get required residue count
        print("\nHow many residue IDs does this transformer need?")
        print("  (e.g., 1 for simple sites, 3 for hemes with propionates)")
        count_str = input(f"Count [{spec.required_residue_count}]: ").strip()
        if count_str:
            spec.required_residue_count = int(count_str)

        # Forcefield path (if applicable)
        if spec.parameters:
            print("\nForcefield path (relative, e.g., 'heme/my_heme'):")
            ff_path = input(f"Path [{spec.forcefield_path or 'none'}]: ").strip()
            if ff_path and ff_path.lower() != 'none':
                spec.forcefield_path = ff_path

        # Summary
        self._display_summary(spec)

        # Confirm
        print("\nGenerate transformer code?")
        confirm = input("Proceed? [Y/n]: ").strip().lower()
        if confirm == 'n':
            print("Cancelled.")
            return None

        # Generate code
        generator = CodeGenerator(spec)
        code = generator.generate()

        # Save to file
        output_file = self.output_dir / f"{spec.name}.py"
        with open(output_file, 'w') as f:
            f.write(code)

        print(f"\n✓ Transformer generated: {output_file}")

        # Offer to save spec
        spec_file = self.output_dir / f"{spec.name}_spec.json"
        print(f"\nSave specification to {spec_file}?")
        save_spec = input("[Y/n]: ").strip().lower()
        if save_spec != 'n':
            with open(spec_file, 'w') as f:
                f.write(spec.model_dump_json(indent=2))
            print(f"✓ Specification saved: {spec_file}")

        return output_file

    def _template_selection_step(self) -> Optional[TransformerSpec]:
        """Step 0: Optional template selection or workspace analysis"""
        print("="*70)
        print("STEP 0: STARTING POINT")
        print("="*70)
        print("\nChoose how to start building your transformer:\n")

        templates = self.template_extractor.list_available_templates()

        # Check if workspace has detected RedoxSites
        has_workspace_sites = False
        if self.workspace:
            detected_sites = self.workspace.get('detected_redox_sites', [])
            has_workspace_sites = len(detected_sites) > 0

        print("Options:")
        print("  0. Start from scratch")

        if templates:
            for i, template in enumerate(templates, 1):
                print(f"  {i}. Use template: {template['name']} - {template['description'][:50]}")

        if has_workspace_sites:
            workspace_choice = len(templates) + 1
            print(f"  W. Analyze workspace RedoxSite ({len(detected_sites)} available)")

        choice_input = input("\nChoice [0]: ").strip() or "0"

        # Handle workspace choice
        if choice_input.upper() == 'W' and has_workspace_sites:
            return self._analyze_workspace_site(detected_sites)

        # Handle numeric choices
        try:
            choice_idx = int(choice_input)
            if choice_idx == 0:
                return None
            elif 1 <= choice_idx <= len(templates):
                template_name = templates[choice_idx - 1]['name']
                print(f"\nExtracting template: {template_name}")
                spec = self.template_extractor.extract_template(template_name)

                if spec:
                    print(f"✓ Template loaded: {template_name}")
                    print("  You can now modify any section.\n")
                    return spec
                else:
                    print("Failed to extract template. Starting from scratch.")
                    return None
            else:
                print("Invalid choice. Starting from scratch.")
                return None

        except ValueError:
            print("Invalid input. Starting from scratch.")
            return None

    def _analyze_workspace_site(self, detected_sites: list) -> Optional[TransformerSpec]:
        """Analyze a RedoxSite from workspace to generate exact requirements"""
        print("\n" + "="*70)
        print("ANALYZE WORKSPACE REDOX SITE")
        print("="*70)
        print("\nSelect a RedoxSite to use as the basis for your transformer.")
        print("The system will extract EXACT requirements (centers, residues, bonds).\n")

        # Display detected sites
        print("Detected RedoxSites:")
        for i, site in enumerate(detected_sites, 1):
            desc = self.site_analyzer.describe_site(site)
            print(f"  {i}. Site {site.site_id}: {desc}")

        # User selects site
        choice = input("\nSelect site number (or blank to cancel): ").strip()
        if not choice:
            print("Cancelled. Starting from scratch.")
            return None

        try:
            site_idx = int(choice) - 1
            if 0 <= site_idx < len(detected_sites):
                site = detected_sites[site_idx]
                print(f"\nAnalyzing site: {site.site_id}")

                # Extract exact requirements
                site_reqs = self.site_analyzer.analyze_redox_site(site)

                # Display what was extracted
                print("\n" + "-"*70)
                print("EXTRACTED REQUIREMENTS (EXACT)")
                print("-"*70)
                self._display_extracted_requirements(site_reqs)

                # Confirm
                use_it = input("\nUse these requirements? [Y/n]: ").strip().lower()
                if use_it == 'n':
                    print("Cancelled. Starting from scratch.")
                    return None

                # Get transformer metadata
                print("\nTransformer Metadata:")
                name = input("  Name (snake_case): ").strip()
                description = input("  Description: ").strip()

                if not name:
                    print("Name required. Cancelled.")
                    return None

                # Create partial spec with exact requirements
                spec = TransformerSpec(
                    name=name,
                    description=description or f"Transformer based on {site.site_id}",
                    site_requirements=site_reqs,
                    component_matching=ComponentMatching(),
                    transformation_sequence=TransformationSequence(),
                    supported_site_types=[name]
                )

                # Mark as workspace-derived so we can skip re-collection
                spec._from_workspace = True

                print(f"\n✓ Requirements extracted from {site.site_id}")
                print("  Proceeding to component matching...\n")

                return spec
            else:
                print("Invalid selection. Starting from scratch.")
                return None

        except ValueError as e:
            logger.error(f"ValueError in workspace site analysis: {e}")
            import traceback
            traceback.print_exc()
            print(f"Error: {e}")
            print("Starting from scratch.")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in workspace site analysis: {e}")
            import traceback
            traceback.print_exc()
            print(f"Unexpected error: {e}")
            print("Starting from scratch.")
            return None

    def _display_extracted_requirements(self, site_reqs: SiteRequirements):
        """Display extracted site requirements for user review"""
        # Centers
        print(f"\nCenters: {site_reqs.centers.required_count} required")
        print(f"  Types: {[ct.value for ct in site_reqs.centers.center_types]}")
        if site_reqs.centers.elements:
            print(f"  Elements: {site_reqs.centers.elements}")

        # Residues
        if site_reqs.atoms and site_reqs.atoms.required_residues:
            print(f"\nResidues: {len(site_reqs.atoms.required_residues)} types")
            for resname, req in sorted(site_reqs.atoms.required_residues.items()):
                print(f"  {resname}: exactly {req.min_count}")

        # Bonds
        if site_reqs.bonds and site_reqs.bonds.required_bond_groups:
            for group in site_reqs.bonds.required_bond_groups:
                print(f"\nBonds: {group.min_count} required")
                for bond_type, bonds in group.bond_types.items():
                    print(f"  {bond_type}: {len(bonds)} bonds")
                    for bond in bonds[:3]:  # Show first 3
                        print(f"    {bond[0][0]}:{bond[0][1]} - {bond[1][0]}:{bond[1][1]}")
                    if len(bonds) > 3:
                        print(f"    ... and {len(bonds)-3} more")

    def _extract_component_names(self, spec: TransformerSpec) -> list:
        """Extract component names from transformations or matching rules"""
        component_names = set()

        # From transformations
        for trans in spec.transformation_sequence.transformations:
            selector = trans.selector.model_dump(exclude_none=True)
            for value in selector.values():
                if isinstance(value, str) and '_' in value and value[0].islower():
                    # Looks like a component reference (e.g., "center_chain")
                    comp_name = value.rsplit('_', 1)[0]  # Remove _chain or _id suffix
                    if comp_name != 'center':  # Center is implicit
                        component_names.add(comp_name)

        # From existing matching rules
        for comp_name in spec.component_matching.rules.keys():
            component_names.add(comp_name)

        return sorted(list(component_names))

    def _display_summary(self, spec: TransformerSpec):
        """Display specification summary"""
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)

        print(f"\nName: {spec.name}")
        print(f"Class: {spec.class_name}")
        print(f"Description: {spec.description}")

        print(f"\nSite Requirements:")
        print(f"  Centers: {spec.site_requirements.centers.required_count} "
              f"{[ct.value for ct in spec.site_requirements.centers.center_types]}")
        if spec.site_requirements.atoms:
            print(f"  Required residues: {len(spec.site_requirements.atoms.required_residues)}")
        if spec.site_requirements.bonds:
            print(f"  Bond groups: {len(spec.site_requirements.bonds.required_bond_groups)}")

        print(f"\nComponent Matching:")
        print(f"  Components: {len(spec.component_matching.rules)}")
        for comp_name in spec.component_matching.rules.keys():
            print(f"    - {comp_name}")

        print(f"\nTransformations:")
        print(f"  Steps: {len(spec.transformation_sequence.transformations)}")
        for trans in spec.transformation_sequence.transformations:
            print(f"    - {trans.id}")

        if spec.parameters:
            print(f"\nParameters:")
            print(f"  Definitions: {len(spec.parameters.definitions)}")
            print(f"  Combinations: {len(spec.parameters.mappings)}")

        print(f"\nResidue IDs needed: {spec.required_residue_count}")
        if spec.forcefield_path:
            print(f"Forcefield path: {spec.forcefield_path}")


def main():
    """Command-line entry point"""
    import sys

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s'
    )

    # Create transformer
    creator = TransformerCreator()

    try:
        output_file = creator.create_transformer()

        if output_file:
            print("\n" + "="*70)
            print("SUCCESS!")
            print("="*70)
            print(f"\nYour transformer has been created: {output_file}")
            print("\nNext steps:")
            print("  1. Review the generated code")
            print("  2. Test with actual redox sites")
            print("  3. Refine matching rules and transformations as needed")
            print("  4. Add to your transformer registry\n")
            sys.exit(0)
        else:
            print("\nTransformer creation cancelled.")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n\nInterrupted by user.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error creating transformer: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
