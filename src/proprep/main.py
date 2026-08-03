#!/usr/bin/env python3
"""
MPSA - Main Entry Point with Automatic Session Recording

This enhanced version automatically records all sessions and checks for
existing session files to offer replay options.
"""

import argparse
import logging
import os
import sys
import traceback
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
import glob
import json
# Prevent readline from loading — it conflicts with Rich's ANSI prompt handling.
# ParmEd's fortranformat/_input.py has a leftover `import pdb` which transitively
# loads readline (pdb → rlcompleter → readline). Once loaded, readline corrupts
# Rich prompts: they vanish when the user types a character and presses backspace.
# Placing a dummy module here blocks the real readline before parmed is imported.
# The dummy must stub the functions that rlcompleter and pdb call at import time.
# See docs/AMBERTOOLS_INTEGRATION.md for the full import chain.
if 'readline' not in sys.modules:
    import types as _types
    _dummy_readline = _types.ModuleType('readline')
    _dummy_readline.set_completer = lambda *a, **kw: None
    _dummy_readline.set_completer_delims = lambda *a, **kw: None
    _dummy_readline.insert_text = lambda *a, **kw: None
    _dummy_readline.redisplay = lambda *a, **kw: None
    _dummy_readline.parse_and_bind = lambda *a, **kw: None
    _dummy_readline.get_completer = lambda: None
    _dummy_readline.get_completer_delims = lambda: ''
    sys.modules['readline'] = _dummy_readline
    del _types, _dummy_readline

from rich.console import Console
from rich.panel import Panel

from proprep.application.pdbprocessor import PDBProcessor
from proprep.utils import integrate_session_manager
from proprep.utils.session_recorder import safe_load_session_file
from proprep.utils.prompts import prompt_with_context, confirm_with_context

logger = logging.getLogger(__name__)


def show_welcome_banner(console: Console) -> None:
    """Display the welcome banner with ASCII art."""
    # ProPrep ASCII art banner
    banner = """
  ____            ____
 |  _ \\ _ __ ___  |  _ \\ _ __ ___ _ __
 | |_) | '__/ _ \\ | |_) | '__/ _ \\ '_ \\
 |  __/| | | (_) ||  __/| | |  __/ |_) |
 |_|   |_|  \\___/ |_|   |_|  \\___| .__/
                                 |_|
"""

    ver = get_version()
    banner_text = (
        f"[bright_blue]{banner}[/bright_blue]\n"
        f"[bold cyan]Perform Protein Preparation like a Pro[/bold cyan]  [grey50]v{ver}[/grey50]\n"
        "Learn the process, skip the pain\n\n"
        "[bright_magenta]Guberman-Pfeffer Lab, New Mexico State University[/bright_magenta]\n"
        "[bright_magenta]With contributions from Adara Walker and John Cairns[/bright_magenta]"
    )

    console.print(Panel(
        banner_text,
        border_style="bright_blue",
        expand=False
    ))

def get_version():
    try:
        return version("proprep")
    except Exception:
        return "unknown"

def setup_logging(verbose=False):
    """Configure logging for the application"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s",
    )

def find_recent_session_files(directory=".", pattern="proprep_session_*.json", max_files=5):
    """
    Find recent session files in the specified directory.
    
    Args:
        directory: Directory to search for session files
        pattern: Glob pattern for session files
        max_files: Maximum number of recent files to return
        
    Returns:
        List of (filepath, metadata) tuples sorted by modification time (newest first)
    """
    session_files = []
    
    # Find all matching files
    for filepath in glob.glob(os.path.join(directory, pattern)):
        try:
            # Get file modification time
            mtime = os.path.getmtime(filepath)
            
            # Try to read metadata from file
            metadata = {}
            try:
                data = safe_load_session_file(filepath, auto_recover=True)
                if data:
                    metadata = data.get('metadata', {})
                    metadata['start_time'] = data.get('start_time', '')
                    metadata['end_time'] = data.get('end_time', '')
                    metadata['interaction_count'] = len(data.get('interactions', []))
            except:
                pass
            
            session_files.append((filepath, mtime, metadata))
        except:
            continue
    
    # Sort by modification time (newest first)
    session_files.sort(key=lambda x: x[1], reverse=True)
    
    # Return only the requested number of files
    return [(f[0], f[2]) for f in session_files[:max_files]]

def browse_session_interactions(session_file):
    """
    Browse and optionally edit a session's interactions.

    Args:
        session_file: Path to session file

    Returns:
        Tuple of (mode, truncate_at, keep_following, new_value, edited_file) where:
        - mode: 'continue', 'edit', 'editor', or 'cancel'
        - truncate_at: Index to truncate at (None for continue/editor mode)
        - keep_following: Whether to keep interactions after truncate_at
        - new_value: New value for edited interaction (None unless mode='edit')
        - edited_file: Path to editor-saved file (only for mode='editor')
    """
    from rich.console import Console
    from rich.panel import Panel

    console = Console()

    # Load session data
    session_data = safe_load_session_file(session_file, auto_recover=True)
    if session_data is None:
        console.print(f"[red]Error: Could not load or recover session file: {session_file}[/red]")
        return ('cancel', None, False, None, None)

    interactions = session_data.get("interactions", [])
    metadata = session_data.get("metadata", {})

    if not interactions:
        console.print("[yellow]Session has no interactions to browse[/yellow]")
        return ('cancel', None, False, None, None)

    # Count display interactions (filter INPUT)
    display_count = sum(1 for i in interactions if i.get("type", "").lower() != "input")

    # Display session info
    console.print(f"\n[bold bright_blue]Session: {os.path.basename(session_file)}[/bold bright_blue]")

    info_lines = []
    if metadata.get('pdb_id'):
        info_lines.append(f"PDB ID: {metadata['pdb_id']}")
    if metadata.get('pdb_file'):
        info_lines.append(f"PDB File: {os.path.basename(metadata['pdb_file'])}")
    info_lines.append(f"Interactions: {display_count}")

    if info_lines:
        console.print(Panel("\n".join(info_lines), title="Session Info", border_style="bright_blue", expand=False))

    # Prompt for action
    console.print("\n[bold]Options:[/bold]")
    console.print("  1. Continue from end (append new interactions)")
    console.print("  2. Open session editor (edit/delete multiple interactions)")
    console.print("  3. Cancel")

    choice = prompt_with_context(None, "\nSelect option", choices=["1", "2", "3"], default="1")

    if choice == "1":
        return ('continue', None, False, None, None)

    elif choice == "2":
        # Launch interactive editor
        from proprep.utils.session_editor import InteractiveSessionEditor

        try:
            editor = InteractiveSessionEditor(session_file, mode="edit")
            edited_file = editor.run()

            if edited_file:
                return ('editor', None, False, None, edited_file)
            else:
                return ('cancel', None, False, None, None)
        except Exception as e:
            console.print(f"[red]Editor error: {e}[/red]")
            return ('cancel', None, False, None, None)

    else:
        return ('cancel', None, False, None, None)


def prompt_for_session_action(session_files):
    """
    Prompt user to select what to do with previous sessions.

    Args:
        session_files: List of (filepath, metadata) tuples

    Returns:
        Tuple of (action, filepath, mode, truncate_at, keep_following, new_value) where:
        - action: 'fresh', 'resume', or 'none'
        - filepath: Selected session file path (None for 'fresh' or 'none')
        - mode: 'continue' or 'edit' (None for fresh/none)
        - truncate_at: Index to truncate at (None unless mode='edit')
        - keep_following: Whether to keep interactions after truncate_at (False for fresh/none)
        - new_value: New value for edited interaction (None unless mode='edit')
    """
    from rich.console import Console
    from rich.table import Table
    
    console = Console()

    if not session_files:
        return ('fresh', None, None, None, False, None)

    # Ask what the user wants to do
    console.print("\n[bold bright_blue]Previous session files found![/bold bright_blue]")
    console.print("\n[bold]What would you like to do?[/bold]")
    console.print("  1. Start fresh session (record new interactions)")
    console.print("  2. Resume a previous session from any point")
    console.print("  3. Continue without session recording")

    action_choice = prompt_with_context(None,
        "\nSelect option ",
        choices=["1", "2", "3"],
        default="1"
    )

    if action_choice == "1":
        # Start fresh
        return ('fresh', None, None, None, False, None)

    elif action_choice == "3":
        # No sessions
        return ('none', None, None, None, False, None)

    elif action_choice == "2":
        # Resume - show session selection
        table = Table(title="Available Session Files")
        table.add_column("#", style="bright_blue", width=3)
        table.add_column("Filename", style="green")
        table.add_column("Date/Time", style="yellow")
        table.add_column("Interactions", style="magenta", width=12)
        table.add_column("Notes", style="blue")

        for idx, (filepath, metadata) in enumerate(session_files, 1):
            filename = os.path.basename(filepath)
            start_time = metadata.get('start_time', 'Unknown')
            if start_time != 'Unknown':
                try:
                    dt = datetime.fromisoformat(start_time)
                    start_time = dt.strftime("%Y-%m-%d %H:%M")
                except:
                    pass

            interactions = str(metadata.get('interaction_count', '?'))

            # Build notes from metadata
            notes = []
            if metadata.get('pdb_id'):
                notes.append(f"PDB: {metadata['pdb_id']}")
            if metadata.get('pdb_file'):
                notes.append(f"File: {os.path.basename(metadata['pdb_file'])}")
            if metadata.get('description'):
                notes.append(metadata['description'])

            notes_str = ", ".join(notes[:2]) if notes else "-"

            table.add_row(str(idx), filename, start_time, interactions, notes_str)

        console.print(table)

        # Get user selection
        choices = [str(i) for i in range(1, len(session_files) + 1)]
        choices.append("0")  # Add option to cancel

        choice = prompt_with_context(None,
            "\nSelect session to resume (0 to cancel) ",
            choices=choices,
            default="0"
        )

        if choice == "0":
            return ('fresh', None, None, None, False, None)

        selected_idx = int(choice) - 1
        selected_file = session_files[selected_idx][0]

        # Browse the session and get user's choice
        result = browse_session_interactions(selected_file)
        mode = result[0]
        edited_file = result[4] if len(result) > 4 else None

        if mode == 'cancel':
            return ('fresh', None, None, None, False, None)

        if mode == 'editor' and edited_file:
            # Editor produced a modified file — replay that instead
            return ('resume', edited_file, 'continue', None, False, None)

        return ('resume', selected_file, mode, None, False, None)

    # Default to fresh
    return ('fresh', None, None, None, False, None)

def backup_session_file(session_file):
    """
    Create a backup of a session file before editing.

    Args:
        session_file: Path to session file to backup

    Returns:
        Path to backup file
    """
    import shutil

    # Generate backup filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file = f"{session_file}.backup_{timestamp}"

    # Copy the file
    shutil.copy2(session_file, backup_file)

    return backup_file


def generate_session_filename(base_dir="."):
    """
    Generate a unique session filename with timestamp.

    Args:
        base_dir: Directory where session file will be saved

    Returns:
        Full path to session file
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"proprep_session_{timestamp}.json"
    return os.path.join(base_dir, filename)

def determine_menu_mode(args):
    """
    Determine which menu mode to use based on command-line flags and saved preferences.

    Priority order:
    1. Command-line flags (--workflow or --full-menu)
    2. Saved user preference
    3. Prompt user to choose

    Args:
        args: Parsed command-line arguments

    Returns:
        str: "workflow" or "full-menu"
    """
    from rich.console import Console
    from rich.panel import Panel
    from proprep.utils.settings_manager import SettingsManager

    console = Console()
    settings_mgr = SettingsManager()

    # Priority 1: Command-line flags
    if args.workflow_mode:
        return "workflow"
    if args.full_menu:
        return "full-menu"

    # Priority 2: Saved preference
    saved_mode = settings_mgr.get_menu_mode()

    if saved_mode == "workflow":
        return "workflow"
    elif saved_mode == "full-menu":
        return "full-menu"

    # Priority 3: Prompt user (saved_mode == "ask")
    console.print()
    console.print(Panel.fit(
        "[bold bright_blue]Welcome to ProPrep![/bold bright_blue]\n\n"
        "ProPrep offers two menu modes:\n\n"
        "[bold yellow]1. Workflow Mode (Recommended for new users)[/bold yellow]\n"
        "   • Guided step-by-step workflow\n"
        "   • Shows only relevant tools for each stage\n"
        "   • Progressive disclosure to avoid overwhelming options\n\n"
        "[bold yellow]2. Full Menu Mode (For experienced users)[/bold yellow]\n"
        "   • All tools visible at once\n"
        "   • Maximum flexibility and control\n"
        "   • No guided workflow progression",
        title="Choose Your Experience",
        border_style="bright_blue"
    ))

    choice = prompt_with_context(None,
        "\n[bold]Select menu mode[/bold] ",
        choices=["1", "2"],
        default=None,  # No default - force explicit choice
        show_choices=True
    )

    # Map choice to mode
    mode = "workflow" if choice == "1" else "full-menu"

    # Ask if they want to save this preference
    save_pref = prompt_with_context(None,
        "\n[grey50]Save this as your default preference?[/grey50] ",
        choices=["y", "n"],
        default="n"
    )

    if save_pref == "y":
        settings_mgr.set_menu_mode(mode)
        console.print(f"[green]✓ Saved '{mode}' as default[/green]")
        console.print(f"[grey50]Tip: Use --set-default to change this later[/grey50]")
    else:
        console.print("[grey50]Tip: Use --workflow or --full-menu flags to skip this prompt[/grey50]")

    console.print()  # Add spacing
    return mode

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ProPrep - Proper Protein Preparation workflow manager"
    )
    parser.add_argument("--pdbid", help="PDB ID to download and process")
    parser.add_argument("--pdbfile", help="Path to local PDB file to process")
    parser.add_argument("--output-dir", help="Directory to store output files")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose output"
    )
    parser.add_argument("--version", action="store_true", help="Show version and exit")

    # Maintainer feedback tools
    feedback_group = parser.add_argument_group('maintainer feedback tools')
    feedback_group.add_argument(
        "--generate-feedback-keypair", action="store_true",
        help="Generate a feedback encryption keypair (maintainer setup): saves the "
             "private key to ~/.proprep and prints the public key to embed."
    )
    feedback_group.add_argument(
        "--decrypt-feedback", metavar="PATH",
        help="Decrypt an encrypted session-context block copied from a feedback "
             "issue. PATH is a file containing the block, or '-' to read stdin."
    )

    # Session recording arguments
    session_group = parser.add_argument_group('session options')
    session_group.add_argument(
        "--no-session", "--no-auto-record",
        action="store_true",
        dest="no_auto_record",
        help="Disable session recording and session file discovery"
    )
    session_group.add_argument(
        "--session-dir",
        default=".",
        help="Directory for session files (default: current directory)"
    )
    session_group.add_argument(
        "--session-description",
        help="Description to include in session metadata"
    )
    session_group.add_argument(
        "--force-replay",
        metavar="FILE",
        help="Force replay of specific session file (skips prompt)"
    )
    session_group.add_argument(
        "--demo-delay",
        type=float,
        default=0.0,
        help="Delay in seconds between replay actions for demo mode (default: 0.0)"
    )
    session_group.add_argument(
        "--resume-session",
        metavar="FILE",
        help="Resume a previous session (allows continuation or editing)"
    )
    
    # Menu mode arguments
    menu_group = parser.add_argument_group('menu mode options')
    menu_mode_mutex = menu_group.add_mutually_exclusive_group()
    menu_mode_mutex.add_argument(
        "--workflow",
        action="store_true",
        dest="workflow_mode",
        help="Start in workflow mode (guided progressive disclosure)"
    )
    menu_mode_mutex.add_argument(
        "--full-menu",
        action="store_true",
        help="Start in full menu mode (show all tools)"
    )
    menu_group.add_argument(
        "--set-default",
        choices=["workflow", "full-menu", "ask"],
        help="Set the default menu mode preference and exit"
    )
    menu_layout_mutex = menu_group.add_mutually_exclusive_group()
    menu_layout_mutex.add_argument(
        "--menu-grid",
        action="store_true",
        help="Render the full menu as a 2-column panel grid (good for screenshots)"
    )
    menu_layout_mutex.add_argument(
        "--menu-list",
        action="store_true",
        help="Render the full menu as a single-column list (default)"
    )

    # Workflow shortcuts
    shortcut_group = parser.add_argument_group('workflow shortcuts')
    shortcut_group.add_argument(
        "--analysis",
        action="store_true",
        help="Jump directly to simulation analysis browser"
    )
    shortcut_group.add_argument(
        "--pdbview",
        metavar="PDB_ID_OR_FILE",
        help="Quick structure viewer: 4-letter PDB ID (e.g., 1ABC) or local .pdb file"
    )

    # Template and batch processing
    batch_group = parser.add_argument_group('template and batch processing')
    batch_group.add_argument(
        "--create-template",
        metavar="SESSION_FILE",
        help="Convert a recorded session to a reusable template"
    )
    batch_group.add_argument(
        "--validate-template",
        metavar="TEMPLATE_FILE",
        help="Validate a template file without executing. Add --batch-list to "
             "also check that those runs supply every variable the template needs"
    )
    batch_group.add_argument(
        "--template-info",
        metavar="TEMPLATE_FILE",
        help="Display information about a template file"
    )
    batch_group.add_argument(
        "--batch-replay",
        metavar="TEMPLATE_FILE",
        help="Process multiple proteins using a template"
    )
    batch_group.add_argument(
        "--batch-list",
        metavar="INPUT_FILE",
        help="File containing protein identifiers (one per line)"
    )
    batch_group.add_argument(
        "--batch-continue",
        action="store_true",
        help="Continue batch processing even if one protein fails"
    )
    batch_group.add_argument(
        "--web",
        dest="batch_web",
        action="store_const",
        const=True,
        default=None,
        help="Replay the batch in proprep-web (browser) mode so mode-gated "
             "prompts match a template recorded under proprep-web. Overrides "
             "the mode recorded in the template."
    )
    batch_group.add_argument(
        "--no-web",
        dest="batch_web",
        action="store_const",
        const=False,
        help="Replay the batch in plain-terminal mode (overrides the template)."
    )

    return parser.parse_args()

def setup_project_directory(project_name=None, interactive=True):
    """
    Create and change to a project directory for ProPrep session.

    All structures and outputs will be saved in this directory.

    Args:
        project_name: Name for the project directory. If None and interactive=True, prompts user.
        interactive: If True, prompts user for project name if not provided.

    Returns:
        str: Path to the created project directory
    """
    import os
    from pathlib import Path
    from rich.console import Console
    
    console = Console()

    # Get project name
    if project_name is None and interactive:
        console.print("\n[bold bright_blue]═══ Project Setup ═══[/bold bright_blue]\n")
        console.print("Enter a name for your project. All structures and outputs")
        console.print("will be saved in a directory with this name.")
        console.print("Press Enter to use the current directory.\n")

        project_name = prompt_with_context(None,
            "Project name",
            default="."
        ).strip()

    if not project_name:
        project_name = "."

    # Create the directory if it doesn't exist
    project_dir = Path(project_name)
    project_dir.mkdir(exist_ok=True)

    # Get absolute path BEFORE changing directory
    abs_path = project_dir.absolute()

    # Change to the project directory
    os.chdir(project_dir)

    console.print(f"[bold green]✓ Project directory: {abs_path}[/bold green]\n")

    return str(abs_path)



def main():
    """Main function - entry point for the application"""

    # Parse command line arguments
    args = parse_arguments()

    # Resolve file paths to absolute before any os.chdir() can invalidate them
    if args.pdbfile:
        args.pdbfile = str(Path(args.pdbfile).resolve())
    if args.pdbview and not (len(args.pdbview) == 4 and args.pdbview.isalnum()):
        args.pdbview = str(Path(args.pdbview).resolve())

    # Validate flag combinations
    if args.analysis and (args.pdbid or args.pdbfile):
        print("Error: --analysis cannot be combined with --pdbid or --pdbfile", file=sys.stderr)
        sys.exit(1)

    # Warn about orphan flags
    if args.demo_delay > 0 and not args.force_replay and not args.resume_session:
        print("Warning: --demo-delay has no effect without --force-replay or --resume-session", file=sys.stderr)
    if args.batch_list and not args.batch_replay and not args.validate_template:
        print("Warning: --batch-list has no effect without --batch-replay or "
              "--validate-template", file=sys.stderr)
    if args.batch_continue and not args.batch_replay:
        print("Warning: --batch-continue has no effect without --batch-replay", file=sys.stderr)
    if args.batch_web is not None and not args.batch_replay:
        print("Warning: --web/--no-web has no effect without --batch-replay", file=sys.stderr)

    # --pdbview is view-only: suppress session recording/discovery
    if args.pdbview:
        args.no_auto_record = True

    # Set up logging
    setup_logging(args.verbose)

    # Quiet noisy loggers unless --verbose is set
    if not args.verbose:
        logging.getLogger('proprep.structure_prep.pdb_loader').setLevel(logging.WARNING)
        logging.getLogger('proprep.structure_prep.pdb_filter_worker').setLevel(logging.WARNING)
        logging.getLogger('proprep.structure_prep.protonation_worker').setLevel(logging.WARNING)

    # Show version and exit if requested
    if args.version:
        from rich.console import Console; console = Console()
        console.print(f"ProPrep v{get_version()} - Proper Protein Preparation workflow manager", highlight=False)
        sys.exit(0)

    # Maintainer feedback-key tools (early exit, no app startup needed)
    if args.generate_feedback_keypair:
        from proprep.utils.feedback_crypto import cli_generate_keypair
        sys.exit(cli_generate_keypair())
    if args.decrypt_feedback:
        from proprep.utils.feedback_crypto import cli_decrypt
        sys.exit(cli_decrypt(args.decrypt_feedback))

    # Handle set-default flag
    if args.set_default:
        from rich.console import Console
        from proprep.utils.settings_manager import SettingsManager

        console = Console()
        settings_mgr = SettingsManager()
        settings_mgr.set_menu_mode(args.set_default)
        console.print(f"[bold green]✓ Default menu mode set to: {args.set_default}[/bold green]")
        sys.exit(0)

    # Handle template operations
    if args.create_template:
        from rich.console import Console
        from proprep.utils.session_editor import InteractiveSessionEditor

        console = Console()
        try:
            editor = InteractiveSessionEditor(args.create_template, mode="template")
            template_file = editor.run()
            sys.exit(0 if template_file else 1)
        except Exception as e:
            console.print(f"[red]Error creating template: {e}[/red]")
            if args.debug:
                traceback.print_exc()
            sys.exit(1)

    if args.validate_template:
        from rich.console import Console
        from proprep.utils.template_converter import (
            TemplateConverter, validate_variable_coverage
        )

        console = Console()
        converter = TemplateConverter(console)
        is_valid, errors = converter.validate_template(args.validate_template)
        warnings = []

        # Given an input list, also check the runs against the template. The
        # template's own shape can be perfectly valid while still being
        # unsatisfiable by these particular runs.
        if is_valid and args.batch_list:
            from proprep.utils.batch_processor import (
                load_input_list, normalize_input_list
            )
            try:
                with open(args.validate_template, 'r') as f:
                    template_data = json.load(f)
                rows = normalize_input_list(load_input_list(args.batch_list))
                cov_errors, warnings = validate_variable_coverage(
                    template_data, rows
                )
                errors.extend(cov_errors)
                is_valid = not errors
            except Exception as e:
                errors.append(f"Could not check {args.batch_list}: {e}")
                is_valid = False

        for warning in warnings:
            console.print(f"[yellow]! {warning}[/yellow]")

        if is_valid:
            console.print(f"[green]✓ Template is valid: {args.validate_template}[/green]")
            sys.exit(0)
        else:
            console.print(f"[red]✗ Template validation failed: {args.validate_template}[/red]")
            for error in errors:
                console.print(f"  [red]• {error}[/red]")
            sys.exit(1)

    if args.template_info:
        from rich.console import Console
        from proprep.utils.template_converter import TemplateConverter

        console = Console()
        converter = TemplateConverter(console)
        converter.display_template_info(args.template_info)
        sys.exit(0)

    # Handle batch processing
    if args.batch_replay:
        if not args.batch_list:
            from rich.console import Console; console = Console()
            console.print("[red]Error: --batch-replay requires --batch-list[/red]")
            console.print("Usage: proprep --batch-replay TEMPLATE_FILE --batch-list INPUT_FILE")
            sys.exit(1)

        from proprep.utils.batch_processor import run_batch_processing

        exit_code = run_batch_processing(
            template_file=args.batch_replay,
            input_list_file=args.batch_list,
            base_dir=args.output_dir or ".",
            continue_on_error=args.batch_continue,
            web_override=args.batch_web
        )
        sys.exit(exit_code)

    # Initialize variables that may be referenced in finally block
    session_file_to_replay = None
    session_file = None

    try:
        # Create PDB processor instance
        processor = PDBProcessor()

        # Show welcome banner (once at startup)
        show_welcome_banner(processor.console)

        # Integrate session recording/replay functionality
        integrate_session_manager(processor)

        # Determine menu mode BEFORE session recording starts
        # This allows the mode choice to be deterministic for session replay
        menu_mode = determine_menu_mode(args)
        processor.workspace.set("menu_mode", menu_mode)
        logger.debug(f"Menu mode set to: {menu_mode}")

        # Optional per-run full-menu layout override (else the persisted
        # SettingsManager preference applies inside the menu command).
        if getattr(args, "menu_grid", False):
            processor.workspace.set("menu_layout", "grid")
        elif getattr(args, "menu_list", False):
            processor.workspace.set("menu_layout", "list")

        # Set workflow shortcut flags if specified
        if args.analysis:
            processor.workspace.set("jump_to_analysis", True)
            logger.debug("Analysis shortcut flag set")

        if args.pdbview:
            processor.workspace.set("jump_to_pdbview", True)
            processor.workspace.set("pdbview_target", args.pdbview)
            logger.debug(f"PDB view shortcut flag set: {args.pdbview}")
        
        # Determine session mode (variables already initialized above)
        hybrid_mode = False
        truncate_at = None
        keep_following = False
        new_value = None

        # Check for resume session first
        if args.resume_session:
            if not os.path.exists(args.resume_session):
                logger.error(f"Session file not found: {args.resume_session}")
                sys.exit(1)

            # Browse the session and get user's choice
            result = browse_session_interactions(args.resume_session)
            mode = result[0]

            if mode == 'continue':
                # Continue from end in hybrid mode
                session_file_to_replay = args.resume_session
                hybrid_mode = True
                truncate_at = None

            elif mode == 'editor':
                # Interactive editor produced a modified session file
                edited_file = result[4]
                session_file_to_replay = edited_file
                hybrid_mode = True
                truncate_at = None

            else:
                # Cancel - start fresh
                session_file_to_replay = None

        # Check for forced replay
        elif args.force_replay:
            session_file_to_replay = args.force_replay
            if not os.path.exists(session_file_to_replay):
                logger.error(f"Session file not found: {session_file_to_replay}")
                sys.exit(1)

        # Otherwise, check for existing sessions and prompt
        # Skip session discovery for --analysis (analyzing existing data, not a new prep)
        elif not args.no_auto_record and not args.analysis:
            # Look for recent session files
            recent_sessions = find_recent_session_files(
                directory=args.session_dir,
                max_files=10
            )

            if recent_sessions:
                # Prompt user to select what to do with sessions
                action, selected_file, mode, truncate_at, keep_following, new_value = prompt_for_session_action(recent_sessions)

                if action == 'resume':
                    # Resume mode - hybrid (replay then record)
                    session_file_to_replay = selected_file
                    hybrid_mode = True
                    # For 'continue' mode (including editor-produced files),
                    # truncate_at is None — just replay the full session
                    # then record new interactions

                elif action == 'none':
                    # User chose to continue without session recording
                    args.no_auto_record = True
                    session_file_to_replay = None

                # For 'fresh' action, session_file_to_replay remains None (starts fresh recording)

        # Set debug mode if requested
        if args.debug:
            processor.workspace.set("debug", True)

        # Setup working directory BEFORE session recording
        # This ensures session files are saved to the correct directory
        if args.output_dir:
            # --output-dir sets the project directory explicitly
            output_dir = Path(args.output_dir)
            if not output_dir.exists():
                output_dir.mkdir(parents=True)
            os.chdir(output_dir)
            project_dir = str(output_dir.absolute())
            processor.workspace.set("project_directory", project_dir)
            logger.info(f"Working directory set to: {project_dir}")
        elif session_file_to_replay:
            # When replaying, we're already in the project directory (where session file is)
            project_dir = os.getcwd()
            processor.workspace.set("project_directory", project_dir)
            if args.verbose or (hasattr(args, 'demo_delay') and args.demo_delay > 0):
                print(f"[Session replaying from directory: {project_dir}]")
        elif args.analysis or args.pdbid or args.pdbfile or args.pdbview:
            # These flags all use CWD as project directory (no subdirectory, no prompt)
            project_dir = os.getcwd()
            processor.workspace.set("project_directory", project_dir)
            logger.info(f"Using current directory: {project_dir}")
        else:
            # Interactive mode: prompt for project directory
            project_dir = setup_project_directory(interactive=True)
            processor.workspace.set("project_directory", project_dir)
            logger.debug(f"Project directory: {project_dir}")

        # NOW start session recording (after directory setup)
        if session_file_to_replay:
            from rich.console import Console; console = Console()
            demo_delay = args.demo_delay if hasattr(args, 'demo_delay') else 0.0

            if hybrid_mode:
                # Hybrid mode: replay then record
                processor.start_session_hybrid(
                    session_file_to_replay,
                    truncate_at=truncate_at,
                    delay=demo_delay,
                    keep_following=keep_following,
                    new_value=new_value
                )
                logger.info(f"Hybrid session started: {session_file_to_replay}")

            else:
                # Replay-only mode
                processor.start_session_replay(session_file_to_replay, demo_delay)
                logger.info(f"Session replay started: {session_file_to_replay}")

                if demo_delay > 0:
                    console.print(f"\n[bold green]Replaying session from: {session_file_to_replay}[/bold green]")
                    console.print(f"[yellow]Demo mode: {demo_delay}s delay between actions[/yellow]")
                else:
                    console.print(f"\n[bold green]Replaying session from: {session_file_to_replay}[/bold green]")

        elif not args.no_auto_record:
            # Automatic recording mode
            session_file = generate_session_filename(args.session_dir)

            # Build metadata
            metadata = {
                "proprep_version": get_version(),
                "command_line": " ".join(sys.argv),
                "working_directory": os.getcwd(),
                "auto_recorded": True,
                # Record whether this session was captured under proprep-web
                # (PROPREP_WEB_SHELL set). Some prompts are mode-gated (e.g. the
                # PDB Filter "view structure?" Y/N is skipped in web mode), so a
                # faithful replay must run in the SAME mode it was recorded in.
                # The batch runner reads this back to match the mode automatically.
                "web_shell_mode": bool(os.environ.get("PROPREP_WEB_SHELL")),
            }

            # Add optional metadata
            if args.session_description:
                metadata["description"] = args.session_description
            if args.pdbid:
                metadata["pdb_id"] = args.pdbid
            if args.pdbfile:
                metadata["pdb_file"] = args.pdbfile

            # Start recording
            from rich.console import Console; console = Console()
            processor.start_session_recording(session_file, metadata)
            # Session recording info shown to user via console.print below
            console.print(f"\nSession recording to: {session_file}")

            # Inform user how to disable recording
            console.print("Tip: Use --no-session to disable automatic session recording\n")

        # Load PDB structure if specified via command line
        if args.pdbid:
            loader = processor.get_module_instance("Structure Loader")
            workspace = processor._get_workspace()
            loader._download_and_load_pdb(args.pdbid, workspace)
            logger.info(f"Downloaded and loaded PDB: {args.pdbid}")

        elif args.pdbfile:
            if not Path(args.pdbfile).exists():
                logger.error(f"PDB file not found: {args.pdbfile}")
                sys.exit(1)
            loader = processor.get_module_instance("Structure Loader")
            workspace = processor._get_workspace()
            loader._load_local_file_by_path(args.pdbfile, workspace)
            logger.info(f"Loaded PDB file: {args.pdbfile}")

        # Run the main menu
        processor.run_main_menu()

    except KeyboardInterrupt:
        from rich.console import Console; console = Console()
        logger.info("Application interrupted by user")
        console.print("\n[yellow]Application interrupted by user[/yellow]")
    except Exception as e:
        logger.error(f"Unhandled exception: {e}")
        if args.debug:
            traceback.print_exc()
        else:
            from rich.console import Console
            from rich.markup import escape
            console = Console()
            # Escape the error message to prevent Rich markup parsing issues
            escaped_error = escape(str(e))
            console.print(f"[red]Error: {escaped_error}[/red]")
            console.print("Run with --debug flag for more details")
        sys.exit(1)
    finally:
        # Clean up processor (this will also stop session recording/replay)
        if "processor" in locals():
            processor.cleanup()
            logger.debug("Application cleanup completed")
            
            # Show session file location if recording
            if (not args.no_auto_record and
                not session_file_to_replay and
                hasattr(processor, 'session_manager') and
                processor.session_manager.is_recording()):
                console.print(f"\n[bold green]Session saved to: {session_file}[/bold green]")

            # Citation and feedback message
            from rich.console import Console
            from rich.panel import Panel
            console = Console()
            console.print()
            console.print(Panel(
                "[bold]If you found ProPrep useful, please cite:[/bold]\n\n"
                "Walker, A. & Guberman-Pfeffer, M.J. (2026). ProPrep: An Interactive\n"
                "and Instructional Interface for Proper Protein Preparation with AMBER.\n"
                "[italic]bioRxiv[/italic]. https://doi.org/10.64898/2026.02.26.708365\n\n"
                "[green]Feedback (positive or negative) is welcome at: mjgp@nmsu.edu[/green]",
                title="[bold white]Thank you for using ProPrep[/bold white]",
                border_style="bright_blue",
                expand=False,
            ))

if __name__ == "__main__":
    main()
