"""
Layout Helpers for Molecular Dynamics Manager
Provides consistent, clean layouts for the 5-step workflow.
"""

from typing import List, Dict, Optional, Tuple, Any
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from proprep.utils.prompts import prompt_with_context


class WorkflowLayoutFormatter:
    """Handles consistent layout formatting across all workflow steps."""

    STEP_NAMES = [
        "Structure Files (Topology + Coordinates)",
        "Protocol Selection",
        "Restraint Integration",
        "Queue Management",
        "Hardware Config",
        "Review & Save"
    ]

    def __init__(self, console: Optional[Console] = None):
        self.console = console or Console()

    def progress_breadcrumb(self, current_step: int) -> str:
        """
        Generate compact progress breadcrumb.

        Args:
            current_step: 0-5 (Step 0 through Step 5)

        Returns:
            Progress string like: "Files ✓ → Workflow ✓ → [Restraints] → Queue → Hardware → Review"
        """
        step_abbrev = ["Files", "Protocol", "Restraints", "Queue", "Hardware", "Review"]
        parts = []

        for i, name in enumerate(step_abbrev):
            if i < current_step:
                parts.append(f"{name} ✓")
            elif i == current_step:
                parts.append(f"[{name}]")
            else:
                parts.append(name)

        return " → ".join(parts)

    def step_header(self, step_num: int, step_name: Optional[str] = None) -> Panel:
        """
        Create compact step header with progress.

        Args:
            step_num: 0-5 (Step 0 through Step 5)
            step_name: Optional override for step name

        Returns:
            Rich Panel with header
        """
        if step_name is None:
            step_name = self.STEP_NAMES[step_num]

        breadcrumb = self.progress_breadcrumb(step_num)

        content = Text()
        content.append(f"Step {step_num + 1} of 6: {step_name}\n", style="bold cyan")
        content.append(breadcrumb, style="grey50")

        return Panel(content, border_style="cyan", padding=(0, 1), expand=False)

    def simple_prompt(
        self,
        step_num: int,
        status_lines: List[str],
        actions: List[Tuple[str, str]],
        step_name: Optional[str] = None,
        help_text: Optional[str] = None
    ) -> str:
        """
        Create a clean, simple prompt layout.

        Args:
            step_num: 0-5
            status_lines: List of status info lines (e.g., ["Structure Pairs: None"])
            actions: List of (key, description) tuples
            step_name: Optional step name override
            help_text: Optional help text to show (usually hidden, shown on ?)

        Returns:
            Formatted prompt string
        """
        self.console.print(self.step_header(step_num, step_name))
        self.console.print()

        # Status section
        for line in status_lines:
            self.console.print(f"  {line}")

        if status_lines:
            self.console.print()

        # Actions
        self.console.print("What would you like to do?", style="bold")
        for key, description in actions:
            # Determine style based on key
            if key in ['x', 'q', 'exit']:
                style = "dim"
            elif key in ['n', 'next']:
                style = "green"
            elif key in ['b', 'back']:
                style = "yellow"
            elif key == 'h':
                style = "cyan"
            else:
                style = "white"

            key_str = f"[{style}]{key:>2}[/{style}]"
            self.console.print(f"  {key_str}. {description}")

        self.console.print()

        # Build the input choices string
        choice_keys = [k for k, _ in actions]
        choices_str = "/".join(choice_keys)

        return f"[{choices_str}]"

    def config_table(
        self,
        step_num: int,
        headers: List[str],
        rows: List[List[str]],
        actions: List[Tuple[str, str]],
        step_name: Optional[str] = None,
        title: Optional[str] = None
    ):
        """
        Create layout with table preview + actions.

        Args:
            step_num: 0-5
            headers: Table column headers
            rows: Table rows (list of lists)
            actions: Available actions
            step_name: Optional step name
            title: Optional table title
        """
        self.console.print(self.step_header(step_num, step_name))
        self.console.print()

        # Create table
        table = Table(title=title, show_header=True, header_style="bold cyan")
        for header in headers:
            table.add_column(header)

        for row in rows:
            table.add_row(*row)

        self.console.print(table)
        self.console.print()

        # Separate navigation from regular actions
        navigation_keys = {'b', 'n', 'x', 'back', 'next', 'exit'}
        regular_actions = [(k, d) for k, d in actions if k not in navigation_keys]
        nav_actions = [(k, d) for k, d in actions if k in navigation_keys]

        # Regular actions
        if regular_actions:
            self.console.print("What would you like to do?", style="bold")
            for key, description in regular_actions:
                self.console.print(f"   {key}. {description}")
            self.console.print()

        # Navigation section with cyan header
        if nav_actions:
            self.console.print("[bold cyan]Navigation:[/bold cyan]")
            for key, description in nav_actions:
                if key in ['x', 'exit']:
                    style = "dim"
                elif key in ['n', 'next']:
                    style = "green"
                elif key in ['b', 'back']:
                    style = "yellow"
                else:
                    style = "white"

                key_str = f"[{style}]{key}[/{style}]"
                self.console.print(f"   {key_str}. {description}")
            self.console.print()

        choice_keys = [k for k, _ in actions]
        return f"[{'/'.join(choice_keys)}]"

    def compact_list(
        self,
        items: List[Dict[str, str]],
        show_numbers: bool = True
    ) -> str:
        """
        Format a compact list of items.

        Args:
            items: List of dicts with 'title' and optional 'description'
            show_numbers: Whether to show item numbers

        Returns:
            Formatted list string
        """
        lines = []
        for i, item in enumerate(items, 1):
            if show_numbers:
                prefix = f"  {i}."
            else:
                prefix = "  •"

            title = item.get('title', '')
            desc = item.get('description', '')

            if desc:
                lines.append(f"{prefix} {title}")
                lines.append(f"     {desc}")
            else:
                lines.append(f"{prefix} {title}")

        return "\n".join(lines)

    def status_line(self, label: str, value: str, style: str = "white") -> str:
        """
        Create a compact status line.

        Args:
            label: Status label (e.g., "Structure Pairs")
            value: Status value (e.g., "1 selected")
            style: Rich style

        Returns:
            Formatted status line
        """
        return f"[bold]{label}:[/bold] [{style}]{value}[/{style}]"

    def inline_navigation_hint(self) -> str:
        """Standard inline navigation hint."""
        return "[grey50](h=help, n=next, b=back, x=exit)[/grey50]"


def format_simulation_name(name: str, max_length: int = 40) -> str:
    """
    Format long simulation names for display.

    Args:
        name: Full simulation name
        max_length: Maximum display length

    Returns:
        Truncated name with ellipsis if needed
    """
    if len(name) <= max_length:
        return name

    # Try to truncate at a sensible point
    parts = name.split('_')
    if len(parts) > 2:
        # Keep first and last parts
        return f"{parts[0]}...{parts[-1]}"
    else:
        return name[:max_length - 3] + "..."


def format_file_path(path: str, max_length: int = 60) -> str:
    """
    Format file paths for compact display.

    Args:
        path: Full file path
        max_length: Maximum display length

    Returns:
        Formatted path
    """
    if len(path) <= max_length:
        return path

    # Show start and end of path
    from pathlib import Path
    p = Path(path)

    if len(p.name) < max_length - 10:
        # Show .../<last_dir>/<filename>
        parent_name = p.parent.name
        return f".../{parent_name}/{p.name}"
    else:
        # Just truncate
        return f"...{path[-(max_length-3):]}"


class TemplatePreviewFormatter:
    """Handles template preview and comparison displays."""

    def __init__(self, console: Optional[Console] = None):
        self.console = console or Console()

    def preview_template(self, template_data: Dict[str, Any]):
        """
        Display template preview showing full annotated MDIN content.

        Args:
            template_data: Template data dictionary
        """
        from rich.panel import Panel
        from rich.syntax import Syntax

        # Header
        name = template_data.get('name', 'Unknown Template')
        description = template_data.get('description', '')

        header = Text()
        header.append(name, style="bold cyan")
        header.append("\n")
        header.append(description, style="grey50")

        self.console.print(Panel(header, border_style="cyan", padding=(0, 1), expand=False))
        self.console.print()

        # Load MDIN file content
        mdin_content = self._load_mdin_content(template_data)

        if mdin_content:
            # Show computed simulation time if applicable
            params = self._parse_mdin_parameters(mdin_content)
            if params:
                sim_time = self._compute_simulation_time(params)
                if sim_time:
                    self.console.print(f"[bold]Simulation time:[/bold] {sim_time}")
                    self.console.print()

            # Show full MDIN with syntax highlighting
            syntax = Syntax(
                mdin_content,
                "fortran",
                theme="monokai",
                line_numbers=True,
                word_wrap=True
            )
            self.console.print(Panel(syntax, border_style="green", padding=(0, 1), expand=False))
        else:
            self.console.print("[yellow]Could not load template content[/yellow]")

    def _compute_simulation_time(self, params: Dict[str, str]) -> Optional[str]:
        """Compute total simulation time from nstlim and dt."""
        try:
            imin = int(params.get('imin', '0'))
            if imin == 1:
                # Minimization — report cycles instead
                maxcyc = params.get('maxcyc')
                if maxcyc:
                    return f"{maxcyc} minimization cycles"
                return None

            nstlim = params.get('nstlim')
            dt = params.get('dt')
            if nstlim and dt:
                total_ps = int(nstlim) * float(dt)
                if total_ps >= 1000:
                    return f"{total_ps / 1000:.1f} ns ({nstlim} steps x {dt} ps)"
                else:
                    return f"{total_ps:.1f} ps ({nstlim} steps x {dt} ps)"
        except (ValueError, TypeError):
            pass
        return None

    def _load_mdin_content(self, template_data: Dict[str, Any]) -> Optional[str]:
        """Load MDIN file content from template data."""
        from pathlib import Path

        # Check if content is already in template_data
        if 'mdin_content' in template_data:
            return template_data['mdin_content']

        # Try to load from template_path
        template_path = template_data.get('template_path')
        if template_path:
            # Resolve builtin templates
            if template_path.startswith('builtin/'):
                # Get the builtin template directory
                from proprep.utils.paths import get_package_dir
                base_dir = get_package_dir() / 'md_templates' / 'builtin'
                template_path = template_path.replace('builtin/', '')
                full_path = base_dir / template_path
            else:
                full_path = Path(template_path)

            try:
                if full_path.exists():
                    with open(full_path, 'r') as f:
                        return f.read()
            except Exception:
                pass

        return None

    def _parse_mdin_parameters(self, mdin_content: str) -> Dict[str, str]:
        """Parse parameters from MDIN file content.

        Handles quoted string values that contain commas, '!', or '#'
        (e.g., AMBER restraintmask values like '@CA,C,O,N&!:WAT').
        """
        import re

        params = {}

        for line in mdin_content.split('\n'):
            line = line.strip()

            # Skip empty lines and full-line comments
            if not line or line.startswith('!') or line.startswith('#'):
                continue

            # Find "name =" prefix; everything after '=' is parsed with quote awareness
            m = re.match(r'(\w+)\s*=\s*(.*)$', line)
            if not m:
                continue

            param_name = m.group(1)
            rest = m.group(2)

            # Walk the value, respecting single/double quotes. Stop at an
            # unquoted ',', '!', or '#' (end-of-record or comment).
            quote = None
            end = len(rest)
            for i, ch in enumerate(rest):
                if quote:
                    if ch == quote:
                        quote = None
                elif ch in ("'", '"'):
                    quote = ch
                elif ch in (',', '!', '#'):
                    end = i
                    break

            param_value = rest[:end].strip().rstrip(',')
            params[param_name] = param_value

        return params

    def _extract_important_mdin_params(self, params: Dict[str, str]) -> List[Tuple[str, str, str]]:
        """Extract most important parameters from parsed MDIN parameters."""
        # Parameter descriptions
        param_descriptions = {
            'imin': 'Run minimization (1) or MD (0)',
            'nstlim': 'Number of MD steps',
            'maxcyc': 'Maximum minimization cycles',
            'ncyc': 'Steepest descent steps',
            'dt': 'Time step (ps)',
            'ntt': 'Temperature control method',
            'temp0': 'Target temperature (K)',
            'tempi': 'Initial temperature (K)',
            'ntp': 'Pressure control',
            'ntb': 'Periodic boundaries',
            'ntr': 'Use position restraints',
            'restraint_wt': 'Restraint force (kcal/mol·Å²)',
            'restraintmask': 'Atoms to restrain',
            'cut': 'Nonbonded cutoff (Å)',
            'ntc': 'SHAKE constraints',
            'ntf': 'Force evaluation',
            'ntwx': 'Trajectory write frequency',
            'ntpr': 'Energy print frequency',
            'ntwr': 'Restart write frequency',
            'ioutfm': 'Binary trajectory format',
            'iwrap': 'Wrap coordinates',
        }

        # Priority order
        priority_params = [
            'imin', 'nstlim', 'maxcyc', 'ncyc', 'dt',
            'ntt', 'temp0', 'tempi', 'ntp', 'ntb',
            'ntr', 'restraint_wt', 'restraintmask', 'cut'
        ]

        important = []

        for param_key in priority_params:
            if param_key in params:
                value = params[param_key]
                desc = param_descriptions.get(param_key, '')
                important.append((param_key, value, desc))

        # Limit to top 10 parameters
        return important[:10]

    def _extract_important_params(self, params: Dict) -> List[Tuple[str, Any, str]]:
        """
        Extract most important parameters for preview.

        Returns:
            List of (param_name, value, description) tuples
        """
        # Priority order for parameters
        priority_params = [
            'imin',      # Simulation type
            'nstlim',    # Number of steps
            'maxcyc',    # Max minimization cycles
            'dt',        # Time step
            'ntt',       # Temperature control
            'temp0',     # Target temperature
            'ntp',       # Pressure control
            'ntb',       # Periodic boundaries
            'ntr',       # Restraints
            'restraint_wt',  # Restraint weight
            'cut',       # Cutoff
        ]

        important = []

        for param_key in priority_params:
            if param_key in params:
                param = params[param_key]
                if isinstance(param, dict):
                    value = param.get('value', param.get('default', ''))
                    desc = param.get('description', '')
                else:
                    value = param
                    desc = ''

                # Only show if value is set
                if value is not None and value != '':
                    important.append((param_key, value, desc[:47] + '...' if len(desc) > 50 else desc))

        # Limit to top 8 parameters
        return important[:8]

    def compare_templates(self, template1: Dict, template2: Dict):
        """
        Display side-by-side comparison of two templates.

        Args:
            template1: First template data
            template2: Second template data
        """
        from rich.columns import Columns
        from rich.panel import Panel

        name1 = template1.get('name', 'Template 1')
        name2 = template2.get('name', 'Template 2')

        self.console.print(f"\n[bold cyan]Comparing: {name1} vs {name2}[/bold cyan]\n")

        # Load and parse MDIN files
        mdin1 = self._load_mdin_content(template1)
        mdin2 = self._load_mdin_content(template2)

        if not mdin1 or not mdin2:
            self.console.print("[yellow]Could not load one or both template files[/yellow]")
            return

        params1 = self._parse_mdin_parameters(mdin1)
        params2 = self._parse_mdin_parameters(mdin2)

        # Find all unique parameter keys
        all_keys = set(params1.keys()) | set(params2.keys())

        # Filter to important parameters for cleaner comparison
        important_keys = [
            'imin', 'nstlim', 'maxcyc', 'ncyc', 'dt',
            'ntt', 'temp0', 'tempi', 'ntp', 'ntb',
            'ntr', 'restraint_wt', 'restraintmask', 'cut',
            'ntc', 'ntf'
        ]

        # Use important keys if they exist, otherwise show all
        display_keys = [k for k in important_keys if k in all_keys]
        if not display_keys:
            display_keys = sorted(all_keys)

        # Build comparison table
        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("Parameter", style="white", width=15)
        table.add_column(name1[:20], style="green", width=25)
        table.add_column(name2[:20], style="blue", width=25)
        table.add_column("Difference", style="yellow", width=15)

        for key in display_keys:
            val1 = params1.get(key, "—")
            val2 = params2.get(key, "—")

            # Highlight differences
            if val1 != val2:
                diff = "✗ Different"
                style1 = "bold green"
                style2 = "bold blue"
            else:
                diff = "✓ Same"
                style1 = "dim green"
                style2 = "dim blue"

            table.add_row(
                key,
                f"[{style1}]{val1}[/{style1}]",
                f"[{style2}]{val2}[/{style2}]",
                diff
            )

        self.console.print(table)

    def _get_param_value(self, param) -> str:
        """Extract value from parameter (handles both dict and direct values)."""
        if param is None:
            return "—"
        if isinstance(param, dict):
            return str(param.get('value', param.get('default', '—')))
        return str(param)


def show_template_selection_with_preview(templates: List[Dict], console: Optional[Console] = None, processor=None):
    """
    Interactive template selection with preview capability.

    Args:
        templates: List of template dictionaries
        console: Rich console instance

    Returns:
        Selected template indices or None
    """
    console = console or Console()
    previewer = TemplatePreviewFormatter(console)

    while True:
        # Show template list
        console.print("\n[bold]Available Templates:[/bold]")
        for i, template in enumerate(templates, 1):
            name = template.get('name', 'Unknown')
            desc = template.get('description', '')
            console.print(f"  {i}. {name}")
            if desc:
                console.print(f"     [grey50]{desc[:70]}...[/grey50]" if len(desc) > 70 else f"     [grey50]{desc}[/grey50]")

        console.print("\n[bold]What would you like to do?[/bold]")
        console.print("  • Enter numbers to select (e.g., '1' or '1,3,5')")
        console.print("  • Type 'p <number>' to preview (e.g., 'p 1')")
        console.print("  • Type 'c <n1> <n2>' to compare (e.g., 'c 1 2')")
        console.print("  • Type 'all' to select all")
        console.print("  • Type 'b' to go back")
        console.print()

        choice = prompt_with_context(
            processor,
            "Your selection",
            module="MD Manager - Template Preview",
            description="Template selection (numbers, 'p N' to preview, 'c' to compare, 'all', 'b')",
        ).strip()

        # Handle preview
        if choice.startswith('p '):
            try:
                idx = int(choice.split()[1]) - 1
                if 0 <= idx < len(templates):
                    console.print()
                    previewer.preview_template(templates[idx])

                    input("\nPress Enter to continue...")
                else:
                    console.print("[red]Invalid template number[/red]")
            except (ValueError, IndexError):
                console.print("[red]Invalid preview command. Use 'p <number>'[/red]")
            continue

        # Handle comparison
        if choice.startswith('c '):
            try:
                parts = choice.split()
                idx1 = int(parts[1]) - 1
                idx2 = int(parts[2]) - 1
                if 0 <= idx1 < len(templates) and 0 <= idx2 < len(templates):
                    previewer.compare_templates(templates[idx1], templates[idx2])
                    input("\nPress Enter to continue...")
                else:
                    console.print("[red]Invalid template numbers[/red]")
            except (ValueError, IndexError):
                console.print("[red]Invalid compare command. Use 'c <n1> <n2>'[/red]")
            continue

        # Handle selection
        if choice == 'b':
            return None

        if choice == 'all':
            return list(range(len(templates)))

        # Parse number selection
        try:
            if ',' in choice:
                indices = [int(x.strip()) - 1 for x in choice.split(',')]
            else:
                indices = [int(choice) - 1]

            # Validate indices
            if all(0 <= idx < len(templates) for idx in indices):
                return indices
            else:
                console.print("[red]Invalid template number(s)[/red]")
        except ValueError:
            console.print("[red]Invalid input. Enter numbers, 'p <n>', 'c <n1> <n2>', 'all', or 'b'[/red]")
