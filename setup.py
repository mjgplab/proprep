"""
ProPrep - Protein Preparation Assistant
A modular workflow manager for processing PDB structures for MD simulations.

This setup.py is designed for integration with the AmberTools build system.
It can also be used standalone: python setup.py install
"""

from setuptools import setup, find_packages

setup(
    name="proprep",
    version="1.15.0",
    description=(
        "Protein Preparation Assistant - A modular workflow manager "
        "for processing PDB structures for MD simulations."
    ),
    author="Matthew J. Guberman-Pfeffer",
    author_email="mjgp@nmsu.edu",
    license="MIT",
    python_requires=">=3.12",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    # IMPORTANT: keep this block in sync with pyproject.toml. AmberTools'
    # CMake build invokes `./setup.py install` directly (see
    # cmake/PythonBuildConfig.cmake → install_python_library()), so this
    # file is authoritative for the AmberTools install path. conda+pip
    # use pyproject.toml.
    package_data={
        "proprep": [
            # Force field parameter files
            "forcefield_params/**/*.frcmod",
            "forcefield_params/**/*.lib",
            "forcefield_params/**/*.mol2",
            "forcefield_params/**/*.json",
            # MD templates and workflows (builtin only)
            "md_templates/builtin/**/*.json",
            "md_templates/builtin/**/*.mdin",
            "md_workflows/builtin/*.json",
            # MD prep builtin data
            "md_prep/templates/annotated/builtin/*.json",
            "md_prep/user_data/templates/builtin/*.json",
            "md_prep/user_data/workflow_presets/builtin/*.json",
            "md_prep/cluster_profiles/*.json",
            # PB titrate prebuilt model compounds (ACE-X-NME tripeptides etc.)
            "pb_titrate/models/*.prmtop",
            "pb_titrate/models/*.rst7",
            # Structure prep templates
            "structure_prep/templates/*",
            # Web shell assets (static frontend + Jinja templates)
            "web/static/*",
            "web/templates/*",
        ],
    },
    entry_points={
        "console_scripts": [
            "proprep = proprep.main:main",
            "mpsa = proprep.main:main",
            "proprep-web = proprep.web.__main__:main",
        ],
    },
    install_requires=[
        "biopython>=1.83,<1.86",  # 1.86 changed PDBIO format string, producing 79-char lines
        "rich>=13.0",
        "requests>=2.28",
        "freesasa>=2.2",
        "pyyaml>=6.0",
        "numpy>=1.26,<2.0",
        "networkx",
        "scikit-learn>=1.3.0",
        "sympy>=1.12",
        "pydantic>=2.0,<3.0",
        "psutil>=5.9.0",
        "pdb2pqr>=3.5",
    ],
    extras_require={
        "tm-align": ["tmtools>=0.3.0"],  # not on conda-forge; optional feature
        "web": [
            "fastapi>=0.110",      # web shell server (proprep.web)
            "uvicorn[standard]>=0.27",
            "httpx>=0.27",         # async reverse proxy to the viewer server
        ],
        "feedback": [
            "pynacl>=1.5",         # encrypt session context attached to feedback issues
        ],
        "test": [],
        "dev": [
            "black==25.1.0",
            "coverage==7.8.0",
            "isort==6.0.1",
            "mypy==1.15.0",
            "pylint==3.3.4",
        ],
    },
)
