#!/usr/bin/env python3
"""
Simple test of Transformer Creator v2 components

Tests data models and code generation without requiring full ProPrep imports
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

# Test imports directly
try:
    import data_models
    import code_generator

    from data_models import (
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
    from code_generator import CodeGenerator
    print("✓ Imports successful\n")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)


def test_data_models():
    """Test that data models work correctly"""
    print("="*70)
    print("TEST 1: Data Models")
    print("="*70)

    # Create a simple spec
    spec = TransformerSpec(
        name="test_transformer",
        class_name="TestTransformer",
        description="A test transformer for validation",
        supported_site_types=["test_type"],
        site_requirements=SiteRequirements(
            centers=CenterRequirements(
                required_count=1,
                center_types=[CenterTypeEnum.METAL_ION],
                elements=["FE"]
            ),
            atoms=AtomRequirements(
                required_residues={
                    "HIS": ResidueRequirement(min_count=2, max_count=2),
                    "CYS": ResidueRequirement(min_count=1, max_count=None)
                }
            )
        ),
        component_matching=ComponentMatching(
            rules={
                "proximal_his": ComponentMatchingRule(
                    description="Proximal histidine ligand",
                    stages=[
                        IntrinsicFilterStage(resname="HIS", coord_atom="NE2"),
                        PositionalStage(sort_by="resid", position=0, reverse=False)
                    ]
                ),
                "distal_his": ComponentMatchingRule(
                    description="Distal histidine ligand",
                    stages=[
                        IntrinsicFilterStage(resname="HIS", coord_atom="NE2"),
                        PositionalStage(sort_by="resid", position=1, reverse=False)
                    ]
                )
            }
        ),
        transformation_sequence=TransformationSequence(
            transformations=[
                Transformation(
                    id="rename_proximal_his",
                    description="Rename proximal histidine to HIE",
                    selector=TransformationSelector(
                        chain_id="proximal_his_chain",
                        residue_id="proximal_his_id"
                    ),
                    action=TransformationAction(
                        change_residue_name="HIE"
                    )
                ),
                Transformation(
                    id="rename_distal_his",
                    description="Rename distal histidine to HID",
                    selector=TransformationSelector(
                        chain_id="distal_his_chain",
                        residue_id="distal_his_id"
                    ),
                    action=TransformationAction(
                        change_residue_name="HID"
                    )
                )
            ]
        ),
        required_residue_count=1
    )

    print(f"✓ Created spec: {spec.name}")
    print(f"  - Class name: {spec.class_name}")
    print(f"  - Center type: {spec.site_requirements.centers.center_types[0].value}")
    print(f"  - Required residues: {len(spec.site_requirements.atoms.required_residues)}")
    print(f"  - Component rules: {len(spec.component_matching.rules)}")
    print(f"  - Transformations: {len(spec.transformation_sequence.transformations)}")

    # Test JSON serialization (Pydantic v2)
    json_str = spec.model_dump_json(indent=2)
    print(f"\n✓ JSON serialization works ({len(json_str)} chars)")

    # Test validation (Pydantic v2)
    spec_dict = spec.model_dump()
    reconstructed = TransformerSpec(**spec_dict)
    print(f"✓ Validation works (reconstructed from dict)")

    return spec


def test_code_generation(spec):
    """Test code generation"""
    print("\n" + "="*70)
    print("TEST 2: Code Generation")
    print("="*70)

    generator = CodeGenerator(spec)
    code = generator.generate()

    print(f"✓ Generated code ({len(code)} chars, {len(code.splitlines())} lines)")

    # Check key components
    checks = [
        ("Class definition", "class TestTransformer"),
        ("Decorator", "@register_redox_transformer"),
        ("evaluate_redox_site method", "def evaluate_redox_site"),
        ("match_components method", "def match_components"),
        ("get_transformation_sequence", "def get_transformation_sequence"),
        ("get_site_requirements", "def get_site_requirements"),
        ("Proximal HIS matching", "proximal_his"),
        ("Distal HIS matching", "distal_his"),
        ("Transformation rename", "rename_proximal_his")
    ]

    print("\nChecking generated code content:")
    for check_name, check_string in checks:
        if check_string in code:
            print(f"  ✓ {check_name}")
        else:
            print(f"  ❌ {check_name} - NOT FOUND!")
            return None

    # Write to test file
    test_file = Path(__file__).parent / "test_generated_transformer.py"
    with open(test_file, 'w') as f:
        f.write(code)

    print(f"\n✓ Wrote test file: {test_file}")
    print(f"  You can review the generated code at this location")

    return code


def show_code_snippet(code):
    """Show a snippet of generated code"""
    print("\n" + "="*70)
    print("CODE SNIPPET (first 50 lines)")
    print("="*70)

    lines = code.splitlines()
    for i, line in enumerate(lines[:50], 1):
        print(f"{i:3d} | {line}")

    if len(lines) > 50:
        print(f"... ({len(lines) - 50} more lines)")


def main():
    """Run all tests"""
    print("\n" + "="*70)
    print("TRANSFORMER CREATOR V2 - COMPONENT TEST")
    print("="*70)
    print()

    try:
        # Test data models
        spec = test_data_models()

        # Test code generation
        code = test_code_generation(spec)

        if code:
            # Show snippet
            show_code_snippet(code)

            print("\n" + "="*70)
            print("ALL TESTS PASSED ✓")
            print("="*70)
            print("\nThe core components are working correctly!")
            print("\nNext steps:")
            print("  1. Review test_generated_transformer.py")
            print("  2. Test collectors interactively with: python creator.py")
            print("  3. Try creating a real transformer\n")

            return 0
        else:
            print("\n❌ Code generation check failed")
            return 1

    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
