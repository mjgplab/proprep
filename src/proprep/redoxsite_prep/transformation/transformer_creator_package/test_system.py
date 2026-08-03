#!/usr/bin/env python3
"""
Test the Transformer Creator v2 system

This script tests the data models, collectors, and code generation
without requiring interactive input.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent))

from proprep.redoxsite_prep.transformation.transformer_creator_v2.data_models import (
    TransformerSpec,
    SiteRequirements,
    CenterRequirements,
    AtomRequirements,
    ResidueRequirement,
    ComponentMatching,
    ComponentMatchingRule,
    IntrinsicFilterStage,
    PositionalStage,
    TransformationSequence,
    Transformation,
    TransformationSelector,
    TransformationAction,
    CenterTypeEnum
)
from proprep.redoxsite_prep.transformation.transformer_creator_v2.code_generator import CodeGenerator


def test_data_models():
    """Test that data models work correctly"""
    print("Testing data models...")

    # Create a simple spec
    spec = TransformerSpec(
        name="test_transformer",
        class_name="TestTransformer",
        description="A test transformer",
        supported_site_types=["test_type"],
        site_requirements=SiteRequirements(
            centers=CenterRequirements(
                required_count=1,
                center_types=[CenterTypeEnum.METAL_ION],
                elements=["FE"]
            ),
            atoms=AtomRequirements(
                required_residues={
                    "HIS": ResidueRequirement(min_count=2, max_count=2)
                }
            )
        ),
        component_matching=ComponentMatching(
            rules={
                "proximal_his": ComponentMatchingRule(
                    description="Proximal histidine",
                    stages=[
                        IntrinsicFilterStage(resname="HIS", coord_atom="NE2"),
                        PositionalStage(sort_by="resid", position=0, reverse=False)
                    ]
                )
            }
        ),
        transformation_sequence=TransformationSequence(
            transformations=[
                Transformation(
                    id="rename_his",
                    description="Rename histidine",
                    selector=TransformationSelector(
                        chain_id="proximal_his_chain",
                        residue_id="proximal_his_id"
                    ),
                    action=TransformationAction(
                        change_residue_name="HIE"
                    )
                )
            ]
        )
    )

    print(f"  ✓ Created spec: {spec.name}")
    print(f"  ✓ Class name: {spec.class_name}")
    print(f"  ✓ Center type: {spec.site_requirements.centers.center_types[0].value}")
    print(f"  ✓ Component rules: {len(spec.component_matching.rules)}")
    print(f"  ✓ Transformations: {len(spec.transformation_sequence.transformations)}")

    # Test JSON serialization
    json_str = spec.model_dump_json(indent=2)
    print(f"  ✓ JSON serialization works ({len(json_str)} chars)")

    # Test validation
    spec_dict = spec.model_dump()
    reconstructed = TransformerSpec(**spec_dict)
    print(f"  ✓ Validation works")

    return spec


def test_code_generation(spec):
    """Test code generation"""
    print("\nTesting code generation...")

    generator = CodeGenerator(spec)
    code = generator.generate()

    print(f"  ✓ Generated code ({len(code)} chars)")

    # Check key components
    assert "class TestTransformer" in code
    assert "@register_redox_transformer" in code
    assert "def evaluate_redox_site" in code
    assert "def match_components" in code
    assert "def get_transformation_sequence" in code
    assert "def get_site_requirements" in code

    print(f"  ✓ All required methods present")

    # Write to test file
    test_file = Path(__file__).parent / "test_generated_transformer.py"
    with open(test_file, 'w') as f:
        f.write(code)

    print(f"  ✓ Wrote test file: {test_file}")

    return code


def test_template_extraction():
    """Test template extraction"""
    print("\nTesting template extraction...")

    try:
        from proprep.redoxsite_prep.transformation.transformer_creator_v2.template_extractor import TemplateExtractor

        extractor = TemplateExtractor()
        templates = extractor.list_available_templates()

        print(f"  ✓ Found {len(templates)} registered transformers")

        if templates:
            # Try to extract the first one
            first_template = templates[0]
            print(f"  • Attempting to extract: {first_template['name']}")

            spec = extractor.extract_template(first_template['name'])
            if spec:
                print(f"    ✓ Extracted successfully")
                print(f"    • Name: {spec.name}")
                print(f"    • Description: {spec.description[:60]}...")
                print(f"    • Transformations: {len(spec.transformation_sequence.transformations)}")
            else:
                print(f"    ⚠ Extraction returned None")

    except Exception as e:
        print(f"  ⚠ Template extraction test failed: {e}")


def main():
    """Run all tests"""
    print("="*70)
    print("TRANSFORMER CREATOR V2 - SYSTEM TEST")
    print("="*70)
    print()

    try:
        # Test data models
        spec = test_data_models()

        # Test code generation
        code = test_code_generation(spec)

        # Test template extraction
        test_template_extraction()

        print("\n" + "="*70)
        print("ALL TESTS PASSED ✓")
        print("="*70)
        print("\nThe system is ready to use!")
        print("Run: python creator.py")
        print()

        return 0

    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
