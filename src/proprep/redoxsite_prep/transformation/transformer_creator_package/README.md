# Transformer Creator

A simplified, modular system for creating RedoxSite transformers.

**Note**: This is the new transformer creator system, replacing the previous monolithic implementation (backed up as `transformer_creator_v1_backup.py`).

## Overview

This tool helps you create custom transformers for ProPrep's RedoxSite system through an interactive CLI interface. It uses a declarative approach where you specify WHAT you want rather than HOW to implement it.

## Architecture

```
transformer_creator_v2/
├── data_models.py           # Pydantic models for specifications
├── template_extractor.py    # Extract from existing transformers
├── code_generator.py        # Generate Python code
├── collectors/              # Interactive data collectors
│   ├── site_requirements_collector.py
│   ├── component_matcher_collector.py
│   ├── transformation_collector.py
│   └── parameter_collector.py
└── creator.py               # Main orchestrator
```

## Usage

### Command Line

```bash
cd /path/to/proprep/src/proprep/redoxsite_prep/transformation/transformer_creator_v2
python creator.py
```

### Programmatic

```python
from proprep.redoxsite_prep.transformation.transformer_creator_v2 import TransformerCreator

creator = TransformerCreator()
output_file = creator.create_transformer()
```

## Workflow

The creator guides you through 4-5 steps:

### Step 0: Template Selection (Optional)
- Choose an existing transformer to use as a starting point
- Or start from scratch

### Step 1: Site Requirements
Define what makes a site match this transformer:
- **Centers**: Metal ions, organic cofactors, or redox amino acids
- **Atoms**: Required residues (e.g., 2 HIS, 1 HEM)
- **Bonds**: Optional bond patterns (advanced)

### Step 2: Component Matching
Define HOW to find components using multi-stage rules:
- **Filter by properties**: Residue name, coordinating atom, distance
- **Sort and select**: By resid, distance, or chain
- **Calculate**: From other components (e.g., proximal_his = c_ring_cys + 1)
- **Exclude**: Already-matched components
- **Custom**: Your own logic

### Step 3: Transformation Sequence
Define transformations to apply:
- **Selector**: Which atoms (by chain, resid, resname, atom names)
- **Action**: What to do (rename residue, move atoms, change IDs)

### Step 4: Parameters (Optional)
Define user-configurable parameters:
- **Definitions**: What parameters exist (redox_state, spin_state, etc.)
- **Mappings**: How they map to residue names

## Example

Creating a simple zinc finger transformer:

```
Name: zinc_finger_c2h2
Description: Classic C2H2 zinc finger

Site Requirements:
  - 1 metal center (ZN)
  - 2 CYS, 2 HIS

Component Matching:
  - cys1: Filter CYS coordinating through SG, select first by resid
  - cys2: Filter CYS coordinating through SG, select second by resid
  - his1: Filter HIS coordinating through NE2, select first by resid
  - his2: Filter HIS coordinating through NE2, select second by resid

Transformations:
  - Rename CYS to CYX (coordinating cysteines)
  - Rename HIS to HIE/HID based on position

Parameters: None
```

## Generated Output

The creator generates:
1. **Python transformer file**: Complete transformer class ready to use
2. **JSON specification file**: Complete specification for future reference

## Key Features

- **Declarative**: Describe what you want, not how to implement it
- **Template-based**: Learn from existing transformers
- **Modular**: Each step is independent and can be modified
- **No hidden logic**: Everything is explicit and transparent
- **Simple CLI**: No complex UI dependencies

## Data Model

The core is the `TransformerSpec` Pydantic model:

```python
TransformerSpec(
    name="my_transformer",
    description="Description",
    site_requirements=SiteRequirements(...),
    component_matching=ComponentMatching(...),
    transformation_sequence=TransformationSequence(...),
    parameters=Parameters(...)  # Optional
)
```

All components are validated and can be serialized to JSON.

## Design Principles

1. **Declarative over imperative**: Users describe WHAT, not HOW
2. **Minimal interaction**: Only essential information
3. **Template-optional**: Can start from scratch
4. **Separation of concerns**: UI, data collection, code generation are separate
5. **Full transparency**: Nothing hidden from the user

## Future Enhancements

Potential improvements:
- Web UI for visual transformer creation
- Validation against actual RedoxSite objects
- Auto-testing of generated transformers
- Library of common matching patterns
- Visual site inspection during creation

## Comparison to v1

**Transformer Creator v1**:
- 11 steps with sub-menus
- ~140 methods, 8500+ lines
- Complex navigation, auto-save, pattern learning
- Heavy Rich UI dependency

**Transformer Creator v2**:
- 4-5 clear steps
- Modular architecture
- Simple CLI prompts
- Clean data models
- ~3000 lines, easy to maintain

## Support

For issues or questions, see the main ProPrep documentation.
