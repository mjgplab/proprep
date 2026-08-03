#!/usr/bin/env python3
"""
Build script for ProPrep using PyInstaller

This script builds ProPrep into a standalone executable that can be distributed
without requiring Python or dependencies to be installed.
"""

import os
import sys
import shutil
import subprocess
import platform
from pathlib import Path
import time

def print_header():
    """Print build script header"""
    print("=" * 60)
    print("ProPrep PyInstaller Build Script")
    print("=" * 60)
    print(f"Platform: {platform.system()} {platform.machine()}")
    print(f"Python: {sys.version}")
    print("=" * 60)

def check_requirements():
    """Check if PyInstaller is installed"""
    try:
        import PyInstaller
        print(f"✓ PyInstaller {PyInstaller.__version__} found")
        return True
    except ImportError:
        print("✗ PyInstaller not found. Installing...")
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pyinstaller>=6.0.0'])
            print("✓ PyInstaller installed successfully")
            return True
        except subprocess.CalledProcessError:
            print("✗ Failed to install PyInstaller")
            return False

def clean_build():
    """Clean previous build artifacts"""
    print("\nCleaning previous builds...")
    dirs_to_clean = ['build', 'dist', '__pycache__']
    
    for dir_name in dirs_to_clean:
        if os.path.exists(dir_name):
            print(f"  Removing {dir_name}...")
            shutil.rmtree(dir_name)
    
    # Clean .pyc files recursively
    pyc_count = 0
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file.endswith('.pyc'):
                os.remove(os.path.join(root, file))
                pyc_count += 1
    
    if pyc_count > 0:
        print(f"  Removed {pyc_count} .pyc files")
    
    print("✓ Cleanup completed")

def test_imports():
    """Test that all required modules can be imported"""
    print("\nTesting module imports...")
    
    critical_modules = [
        'proprep.main',
        'Bio.PDB',
        'rich.console',
        'yaml',
        'requests',
        'numpy',
    ]
    
    optional_modules = [
        'freesasa',
        'proprep.structure_prep.pdb_loader',
        'proprep.utils.module_registry',
        'proprep.application.pdbprocessor',
    ]
    
    failed_critical = []
    failed_optional = []
    
    # Test critical modules
    for module in critical_modules:
        try:
            __import__(module)
            print(f"  ✓ {module}")
        except ImportError as e:
            print(f"  ✗ {module}: {e}")
            failed_critical.append(module)
    
    # Test optional modules
    for module in optional_modules:
        try:
            __import__(module)
            print(f"  ✓ {module} (optional)")
        except ImportError as e:
            print(f"  ⚠ {module} (optional): {e}")
            failed_optional.append(module)
    
    if failed_critical:
        print(f"\n✗ Critical import failures: {failed_critical}")
        return False
    
    if failed_optional:
        print(f"\n⚠ Optional import failures: {failed_optional}")
        print("  Build will continue, but some features may not work")
    
    print("✓ Import test completed")
    return True

def build_executable():
    """Build the executable using PyInstaller"""
    print("\nBuilding ProPrep executable...")
    
    # Check if spec file exists
    spec_file = Path('proprep.spec')
    if not spec_file.exists():
        print("✗ proprep.spec file not found!")
        print("  Please create the spec file first (see documentation)")
        return False
    
    # Run PyInstaller
    cmd = [sys.executable, '-m', 'PyInstaller', 'proprep.spec', '--clean', '--noconfirm']
    
    print(f"Running: {' '.join(cmd)}")
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        build_time = time.time() - start_time
        
        print(f"✓ Build completed in {build_time:.1f} seconds")
        
        # Show build results
        dist_dir = Path('dist')
        if dist_dir.exists():
            files = list(dist_dir.glob('*'))
            print(f"\nBuild output in dist/:")
            for file in files:
                if file.is_file():
                    size = file.stat().st_size / (1024 * 1024)  # MB
                    print(f"  {file.name} ({size:.1f} MB)")
                else:
                    print(f"  {file.name}/ (directory)")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Build failed after {time.time() - start_time:.1f} seconds")
        print(f"\nError output:")
        print(e.stderr)
        if e.stdout:
            print(f"\nStandard output:")
            print(e.stdout)
        return False

def test_executable():
    """Test the built executable"""
    print("\nTesting built executable...")
    
    system = platform.system().lower()
    exe_name = 'proprep.exe' if system == 'windows' else 'proprep'
    exe_path = Path('dist') / exe_name
    
    if not exe_path.exists():
        print(f"✗ Executable not found: {exe_path}")
        return False
    
    # Test help command
    try:
        result = subprocess.run([str(exe_path), '--help'], 
                              capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print("✓ Executable runs successfully")
            return True
        else:
            print(f"✗ Executable failed with return code {result.returncode}")
            print(f"stderr: {result.stderr}")
            return False
    except subprocess.TimeoutExpired:
        print("✗ Executable test timed out")
        return False
    except Exception as e:
        print(f"✗ Error testing executable: {e}")
        return False

def create_distribution():
    """Create distribution package"""
    print("\nCreating distribution package...")
    
    system = platform.system().lower()
    arch = platform.machine().lower()
    
    # Normalize architecture names
    if arch in ['x86_64', 'amd64']:
        arch = 'x64'
    elif arch in ['arm64', 'aarch64']:
        arch = 'arm64'
    
    # Create distribution directory
    dist_name = f"proprep-{system}-{arch}"
    releases_dir = Path('releases')
    dist_path = releases_dir / dist_name
    
    if dist_path.exists():
        shutil.rmtree(dist_path)
    
    dist_path.mkdir(parents=True)
    
    # Copy executable
    exe_name = 'proprep.exe' if system == 'windows' else 'proprep'
    exe_path = Path('dist') / exe_name
    
    if exe_path.exists():
        shutil.copy2(exe_path, dist_path / exe_name)
        print(f"  ✓ Copied {exe_name}")
        
        # Copy additional files
        additional_files = ['README.md', 'LICENSE']
        for file_name in additional_files:
            if Path(file_name).exists():
                shutil.copy2(file_name, dist_path / file_name)
                print(f"  ✓ Copied {file_name}")
        
        # Create user guide
        user_guide = dist_path / 'GETTING_STARTED.txt'
        exe_cmd = 'proprep.exe' if system == 'windows' else './proprep'
        with open(user_guide, 'w') as f:
            f.write("ProPrep: Interactive Protein Preparation for AMBER\n")
            f.write("=" * 52 + "\n\n")
            f.write("QUICK START\n")
            f.write("-" * 11 + "\n\n")
            f.write(f"  {exe_cmd}                              Launch interactive interface\n")
            f.write(f"  {exe_cmd} --pdbid 1UBQ                 Load a PDB structure directly\n")
            f.write(f"  {exe_cmd} --pdbfile protein.pdb        Load a local PDB file\n")
            f.write(f"  {exe_cmd} --analysis                   Jump to simulation analysis\n")
            f.write(f"  {exe_cmd} --help                       Show all options\n\n")
            f.write("SETUP\n")
            f.write("-" * 5 + "\n\n")
            f.write("This executable bundles all Python dependencies and MODELLER.\n")
            f.write("AmberTools must be installed separately via conda.\n\n")
            f.write("If you do not have conda, install Miniconda first:\n")
            f.write("  https://docs.conda.io/en/latest/miniconda.html\n\n")
            f.write("1. Install AmberTools:\n\n")
            f.write("     conda create -n proprep_env -c conda-forge ambertools\n")
            f.write("     conda activate proprep_env\n")
            f.write("     export AMBERHOME=$CONDA_PREFIX\n\n")
            f.write("2. Configure MODELLER license key (optional but recommended):\n\n")
            f.write("   MODELLER is bundled but requires a free academic license key.\n")
            f.write("   Register at: https://salilab.org/modeller/registration.html\n\n")
            f.write("   Save your key (one-time setup):\n\n")
            f.write("     mkdir -p ~/.proprep\n")
            f.write("     echo 'YOUR_KEY' > ~/.proprep/modeller_key\n\n")
            f.write("   Without a key, ProPrep runs normally but structure repair\n")
            f.write("   and mutagenesis features are disabled.\n\n")
            f.write("RUNNING PROPREP\n")
            f.write("-" * 15 + "\n\n")
            f.write("Activate the conda environment, then run the executable:\n\n")
            f.write("     conda activate proprep_env\n")
            f.write(f"     {exe_cmd}\n\n")
            f.write("Common options:\n\n")
            f.write(f"  {exe_cmd} --pdbid 1UBQ                 Load a PDB structure directly\n")
            f.write(f"  {exe_cmd} --pdbfile protein.pdb        Load a local PDB file\n")
            f.write(f"  {exe_cmd} --analysis                   Jump to simulation analysis\n")
            f.write(f"  {exe_cmd} --help                       Show all options\n\n")
            f.write("USER DATA\n")
            f.write("-" * 9 + "\n\n")
            f.write("User-created protocols, templates, and settings are stored in\n")
            f.write("~/.proprep/. Session logs and workspace files are saved in the\n")
            f.write("project working directory.\n\n")
            f.write("See README.md for full documentation.\n")
        
        print(f"  ✓ Created GETTING_STARTED.txt")
        
        # Create run script for Unix systems
        if system != 'windows':
            run_script = dist_path / 'run_proprep.sh'
            with open(run_script, 'w') as f:
                f.write('#!/bin/bash\n')
                f.write('# ProPrep launcher script\n')
                f.write('cd "$(dirname "$0")"\n')
                f.write('./proprep "$@"\n')
            run_script.chmod(0o755)
            print(f"  ✓ Created run_proprep.sh")
        
        # Create archive
        try:
            archive_name = f"{dist_name}.zip"
            archive_path = releases_dir / archive_name
            
            shutil.make_archive(str(archive_path.with_suffix('')), 'zip', str(dist_path))
            print(f"  ✓ Created archive: {archive_name}")
            
            # Show final size
            archive_size = archive_path.stat().st_size / (1024 * 1024)
            print(f"  Archive size: {archive_size:.1f} MB")
            
        except Exception as e:
            print(f"  ⚠ Could not create archive: {e}")
        
        print(f"✓ Distribution created: {dist_path}")
        return True
    
    print(f"✗ Executable not found: {exe_path}")
    return False

def main():
    """Main build process"""
    print_header()
    
    # Check requirements
    if not check_requirements():
        print("\n✗ Build failed: PyInstaller not available")
        return 1
    
    # Test imports before building
    if not test_imports():
        print("\n✗ Build failed: Import test failed")
        return 1
    
    # Clean previous builds
    clean_build()
    
    # Build executable
    if not build_executable():
        print("\n✗ Build failed during compilation")
        return 1
    
    # Test the executable
    if not test_executable():
        print("\n⚠ Build completed but executable test failed")
        print("  The executable may still work, but please test manually")
    
    # Create distribution
    if create_distribution():
        print("\n" + "=" * 60)
        print("🎉 BUILD COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print("\nYour executable is ready for distribution:")
        print("  Check the 'releases/' directory for the packaged files")
        print("\nNext steps:")
        print("  1. Test the executable on target systems")
        print("  2. Create a public repository for distribution")
        print("  3. Upload binaries to GitHub releases")
        return 0
    else:
        print("\n✗ Build failed during distribution packaging")
        return 1

if __name__ == '__main__':
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n✗ Build cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
