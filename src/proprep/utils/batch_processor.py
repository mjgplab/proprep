"""
Batch Processor for ProPrep Templates

Enables batch processing of multiple proteins using template sessions.
Automates the execution of the same workflow on many structures.

Supports two input formats:
  - Plain text: one protein identifier per line (backward compatible)
  - CSV with headers: variable names in header row, one run per data row
"""

import csv
import json
import os
import sys
import traceback
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Union

from rich.console import Console
from rich.progress import Progress, TaskID, SpinnerColumn, BarColumn, TextColumn
from rich.panel import Panel
from rich.table import Table


class BatchProcessor:
    """Process multiple proteins using a template session."""

    def __init__(self, console: Optional[Console] = None):
        """Initialize the batch processor."""
        # Pin the console to the real stderr so it keeps working
        # even when sys.stderr is temporarily redirected during runs.
        self.console = console or Console(stderr=True, file=sys.stderr)
        self._real_stderr = sys.stderr
        self.results = []
        self.start_time = None
        self.end_time = None

    def process_batch(
        self,
        template_file: str,
        input_list: List[Union[str, Dict[str, str]]],
        base_dir: str = ".",
        continue_on_error: bool = True,
        web_override: Optional[bool] = None
    ) -> Tuple[int, int]:
        """
        Process a batch of proteins using a template.

        Args:
            template_file: Path to template session file
            input_list: List of protein identifiers (strings) or variable dicts (from CSV)
            base_dir: Base directory for output
            continue_on_error: Continue processing if one protein fails
            web_override: Force proprep-web (True) or plain-terminal (False)
                replay mode, overriding the template's recorded mode. None
                honors the template's recorded ``web_shell_mode``.

        Returns:
            Tuple of (successful_count, failed_count)
        """
        self.start_time = datetime.now()
        self.results = []

        # Validate template
        if not self._validate_template(template_file):
            return (0, len(input_list))

        # Load template to get variable info
        with open(template_file, 'r') as f:
            template_data = json.load(f)

        template_vars = template_data.get("template_variables", {})

        # Resolve the replay mode. Replays run IN-PROCESS and mode is detected
        # live via PROPREP_WEB_SHELL, so a template recorded under proprep-web
        # (where some prompts are skipped) only replays cleanly if we recreate
        # that same env. An explicit --web/--no-web overrides the recording.
        recorded_web = bool(template_data.get("metadata", {}).get("web_shell_mode", False))
        web_mode = recorded_web if web_override is None else web_override

        normalized_list = normalize_input_list(input_list)

        if not normalized_list:
            self.console.print("[yellow]No inputs to process[/yellow]")
            return (0, 0)

        # Check the runs against the template before any of them start. Left
        # to run time, a variable with no value becomes a prompt replay cannot
        # answer, which strict replay can only turn into a failed run.
        from proprep.utils.template_converter import validate_variable_coverage
        cov_errors, cov_warnings = validate_variable_coverage(
            template_data, normalized_list
        )
        for warning in cov_warnings:
            self.console.print(f"[yellow]! {warning}[/yellow]")
        if cov_errors:
            self.console.print(
                f"[red]Input list does not satisfy {template_file}:[/red]"
            )
            for error in cov_errors:
                self.console.print(f"  [red]• {error}[/red]")
            return (0, len(normalized_list))

        # Display batch info
        self._display_batch_header(template_file, normalized_list, template_vars)
        mode_label = "proprep-web (browser)" if web_mode else "plain terminal"
        override_note = "" if web_override is None else " [override]"
        self.console.print(f"[bold]Replay mode:[/bold] {mode_label}{override_note}\n")

        # Recreate the recorded launch environment for the whole batch so
        # mode-gated prompts fire (or don't) exactly as they did at capture
        # time. PROPREP_BATCH additionally tells the viewer layer to stay
        # fully headless — web mode would otherwise try to launch a real
        # localhost viewer server per protein (see _launch_viewer).
        prev_web_shell = os.environ.get("PROPREP_WEB_SHELL")
        prev_batch = os.environ.get("PROPREP_BATCH")
        if web_mode:
            os.environ["PROPREP_WEB_SHELL"] = "1"
        else:
            os.environ.pop("PROPREP_WEB_SHELL", None)
        os.environ["PROPREP_BATCH"] = "1"

        try:
            return self._run_batch_loop(
                template_file, normalized_list, template_vars,
                base_dir, continue_on_error
            )
        finally:
            # Restore the caller's environment verbatim.
            if prev_web_shell is None:
                os.environ.pop("PROPREP_WEB_SHELL", None)
            else:
                os.environ["PROPREP_WEB_SHELL"] = prev_web_shell
            if prev_batch is None:
                os.environ.pop("PROPREP_BATCH", None)
            else:
                os.environ["PROPREP_BATCH"] = prev_batch

    def _run_batch_loop(
        self,
        template_file: str,
        normalized_list: List[Dict[str, str]],
        template_vars: Dict[str, Any],
        base_dir: str,
        continue_on_error: bool
    ) -> Tuple[int, int]:
        """Run the per-protein replay loop (env already configured by caller)."""
        # Process each entry
        successful = 0
        failed = 0

        # Use a progress console pinned to the real stderr file object
        # so it keeps rendering even when sys.stderr is redirected.
        progress_console = Console(file=self._real_stderr, force_terminal=True)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=progress_console
        ) as progress:

            task = progress.add_task(
                f"Processing {len(normalized_list)} proteins...",
                total=len(normalized_list)
            )

            for i, variables in enumerate(normalized_list, 1):
                # Use input_protein as the display label
                display_label = variables.get("input_protein", str(variables))

                progress.update(
                    task,
                    description=f"Processing {display_label} ({i}/{len(normalized_list)})"
                )

                # Process this entry
                success, output_dir, error_msg = self._process_single_protein(
                    template_file=template_file,
                    variables=variables,
                    base_dir=base_dir,
                    template_vars=template_vars,
                    run_index=i
                )

                # Record result
                result = {
                    "protein": display_label,
                    "variables": variables,
                    "success": success,
                    "output_dir": output_dir,
                    "error": error_msg,
                    "timestamp": datetime.now().isoformat()
                }
                self.results.append(result)

                if success:
                    successful += 1
                else:
                    failed += 1
                    if not continue_on_error:
                        self.console.print(f"\n[red]Stopping batch due to error[/red]")
                        break

                progress.advance(task)

        self.end_time = datetime.now()

        # Generate summary
        self._display_batch_summary(template_file, normalized_list)

        # Save batch report
        self._save_batch_report(base_dir, template_file, normalized_list)

        return (successful, failed)

    def _validate_template(self, template_file: str) -> bool:
        """Validate template file exists and is properly formatted."""
        if not os.path.exists(template_file):
            self.console.print(f"[red]Template file not found: {template_file}[/red]")
            return False

        try:
            with open(template_file, 'r') as f:
                template_data = json.load(f)

            if not template_data.get("template"):
                self.console.print(f"[red]File is not a template: {template_file}[/red]")
                return False

            if "template_variables" not in template_data:
                self.console.print(f"[red]Template has no variables defined[/red]")
                return False

            return True

        except json.JSONDecodeError as e:
            self.console.print(f"[red]Invalid template file format: {e}[/red]")
            return False

    def _display_batch_header(
        self,
        template_file: str,
        input_list: List[Union[str, Dict[str, str]]],
        template_vars: Dict[str, Any]
    ):
        """Display batch processing header."""
        self.console.print("\n[bold cyan]═══ Batch Processing Mode ═══[/bold cyan]\n")
        self.console.print(f"[bold]Template:[/bold] {template_file}")
        self.console.print(f"[bold]Runs to process:[/bold] {len(input_list)}")

        # Show template variables
        if template_vars:
            var_info = []
            for var_name, var_def in template_vars.items():
                var_info.append(
                    f"  - {var_name}: {var_def.get('description', 'No description')}"
                )
            self.console.print(f"[bold]Template variables:[/bold]")
            for info in var_info:
                self.console.print(info)

        self.console.print()

    def _process_single_protein(
        self,
        template_file: str,
        variables: Dict[str, str],
        base_dir: str,
        template_vars: Dict[str, Any],
        run_index: int = 1
    ) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        Process a single run using the template with variable substitution.

        Args:
            template_file: Path to template session file
            variables: Dict of variable name -> value for substitution
            base_dir: Base directory for output
            template_vars: Template variable definitions

        Returns:
            Tuple of (success, output_dir, error_message)
        """
        # Generate output directory name
        output_dir_name = f"run_{run_index:02d}"
        output_dir = os.path.join(base_dir, output_dir_name)

        try:
            # Create output directory
            Path(output_dir).mkdir(exist_ok=True)

            # Save current directory
            original_dir = os.getcwd()

            # Change to output directory
            os.chdir(output_dir)

            # Redirect all console output to a log file so the terminal
            # only shows the batch progress bar.
            log_file_path = "run.log"  # already cd'd into output_dir
            log_file = open(log_file_path, "w")
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            sys.stdout = log_file
            sys.stderr = log_file

            # Everything from here on must run inside the try below: once
            # stdout is redirected and cwd has moved, any escape that skips the
            # restoring finally leaves those settings pointing at this run's
            # directory, corrupting every run after it.
            processor = None
            replayer = None
            interceptor = None

            try:
                # Import and run ProPrep with template replay
                from proprep.application.pdbprocessor import PDBProcessor
                from proprep.utils.session_recorder import SessionReplayer, InterceptedPrompt

                # Create processor with output directed to the log file. Do NOT
                # force_terminal here: run.log is a plain file, so forcing terminal
                # mode injects raw ANSI color escapes (e.g. \x1b[38;5;244m) that make
                # the log unreadable — and run.log is the ONLY way to diagnose a
                # stuck/failed batch run. A fixed width keeps panels/tables laid out
                # nicely without color codes. (The on-screen progress bar uses a
                # separate console pinned to the real stderr, so it stays colored.)
                processor = PDBProcessor()
                processor.console = Console(file=log_file, force_terminal=False, width=120)

                # Detect menu mode from template interactions and set it
                # before replay starts, so the correct menu is shown.
                with open(os.path.join(original_dir, template_file), 'r') as tf:
                    tpl_data = json.load(tf)
                for interaction in tpl_data.get("interactions", []):
                    module = interaction.get("context", {}).get("module", "")
                    if "Workflow Mode" in module:
                        processor.workspace.set("menu_mode", "workflow")
                        break
                    elif "All Tools" in module:
                        processor.workspace.set("menu_mode", "all")
                        break

                # Set up session replayer with variable substitution
                replayer = SessionReplayer(
                    replay_file=os.path.join(original_dir, template_file),
                    replay_delay=0.0,
                    variables=variables
                )

                # Install interceptor. Strict replay is what keeps a divergent
                # prompt from blocking on a terminal nobody is watching: stdout
                # is redirected into run.log, so a live prompt would print there
                # and then wait on stdin behind a frozen progress bar.
                interceptor = InterceptedPrompt(replayer=replayer, strict_replay=True)
                interceptor.install()

                # Start replay
                replayer.start_replay()

                # Run the main menu (which will be automated by replay)
                processor.run_main_menu()

                success = True
                error_msg = None

            except SystemExit:
                # The menu exit path calls sys.exit(0) — catch it so the
                # batch loop can continue to the next run.
                success = True
                error_msg = None

            except Exception as e:
                success = False
                error_msg = str(e)

                # Save error log. cwd is already output_dir at this point.
                with open("error.log", 'w') as f:
                    f.write(f"Error processing {variables}\n")
                    f.write(f"Timestamp: {datetime.now().isoformat()}\n")
                    f.write(f"Error: {error_msg}\n\n")
                    f.write("Traceback:\n")
                    f.write(traceback.format_exc())

            finally:
                # Restore stdout/stderr before any cleanup output
                sys.stdout = old_stdout
                sys.stderr = old_stderr
                log_file.close()

                # Clean up. These may be unset if setup raised part way.
                if interceptor is not None:
                    interceptor.uninstall()
                if replayer is not None:
                    replayer.stop_replay()
                if hasattr(processor, 'cleanup'):
                    processor.cleanup()

                # Reset the global module registry so the next run gets
                # fresh module instances with no leftover state (workspace
                # refs, cached paths, accumulated mutations, etc.).
                from proprep.utils.module_registry import registry
                registry.cleanup()
                registry.instances.clear()

                # Return to original directory
                os.chdir(original_dir)

            return (success, output_dir, error_msg)

        except Exception as e:
            # Directory creation or other setup error
            return (False, output_dir, f"Setup error: {str(e)}")

    def _display_batch_summary(self, template_file: str,
                               input_list: List[Dict[str, str]]):
        """Display batch processing summary."""
        planned = len(input_list)
        attempted = len(self.results)
        successful = sum(1 for r in self.results if r["success"])
        failed = attempted - successful
        # The loop appends one result per run in order, so anything past the
        # last result is a run the early stop never reached.
        not_attempted = planned - attempted

        def pct(count):
            return f" ({count / planned * 100:.1f}%)" if planned else ""

        duration = (self.end_time - self.start_time).total_seconds()
        minutes, seconds = divmod(int(duration), 60)

        self.console.print("\n[bold cyan]═══ Batch Processing Summary ═══[/bold cyan]\n")
        self.console.print(f"[bold]Duration:[/bold] {minutes}m {seconds}s")
        self.console.print(f"[bold]Planned:[/bold] {planned}")
        self.console.print(f"[bold green]Successful:[/bold green] {successful}{pct(successful)}")

        if failed > 0:
            self.console.print(f"[bold red]Failed:[/bold red] {failed}{pct(failed)}")

        # Say outright that these runs were never started. Reported as
        # failures they would look like 497 broken structures rather than a
        # batch that stopped after one.
        if not_attempted > 0:
            self.console.print(
                f"[bold yellow]Not attempted:[/bold yellow] {not_attempted}"
                f"{pct(not_attempted)}"
            )
            self.console.print(
                f"[yellow]Batch stopped after run {attempted} of {planned}. "
                f"Pass --batch-continue to work through failures instead of "
                f"stopping at the first one.[/yellow]"
            )

        # Show results table
        table = Table(title="Processing Results", show_header=True, header_style="bold magenta")
        table.add_column("Protein", style="cyan")
        table.add_column("Status", style="white")
        table.add_column("Output Directory", style="yellow")

        for result in self.results:
            status = "[green]✓ Success[/green]" if result["success"] else "[red]✗ Failed[/red]"
            output_dir = result.get("output_dir", "N/A")
            if result.get("error"):
                output_dir += f"\n[grey50]{result['error']}[/grey50]"

            table.add_row(
                result["protein"],
                status,
                output_dir
            )

        self.console.print("\n")
        self.console.print(table)

    def _write_retry_list(self, base_dir: str, timestamp: str,
                          rows: List[Dict[str, str]]) -> Optional[str]:
        """
        Write the runs still owed as an input list a later batch can consume.

        Re-running only what did not succeed is the normal follow-up, and
        rebuilding that by hand from a summary table does not scale to the
        run counts batch processing exists for. This is an ordinary input
        list rather than a bespoke manifest, so it feeds straight back into
        --batch-list with no conversion step.

        Args:
            base_dir: Base directory for outputs
            timestamp: Shared stamp so this file pairs with its JSON report
            rows: Runs to retry, in their original order

        Returns:
            Path to the retry list, or None if every run succeeded
        """
        if not rows:
            return None

        # Union of keys in first-seen order: rows come from the same input
        # list, but a hand-edited CSV may not have filled every column.
        fieldnames: List[str] = []
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)

        retry_file = os.path.join(base_dir, f"batch_retry_{timestamp}.csv")
        with open(retry_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        return retry_file

    def _save_batch_report(self, base_dir: str, template_file: str,
                           input_list: List[Dict[str, str]]):
        """Save batch processing report and a re-runnable list of what is owed."""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        report_file = os.path.join(base_dir, f"batch_report_{timestamp}.json")

        attempted = len(self.results)
        not_attempted = list(input_list[attempted:])
        failed_rows = [r["variables"] for r in self.results if not r["success"]]

        report_data = {
            "batch_info": {
                "template": template_file,
                "start_time": self.start_time.isoformat(),
                "end_time": self.end_time.isoformat(),
                "duration_seconds": (self.end_time - self.start_time).total_seconds(),
                "total_proteins": len(input_list),
                "attempted": attempted,
                "successful": sum(1 for r in self.results if r["success"]),
                "failed": len(failed_rows),
                "not_attempted": len(not_attempted),
                "stopped_early": bool(not_attempted)
            },
            "results": self.results,
            "not_attempted": not_attempted
        }

        with open(report_file, 'w') as f:
            json.dump(report_data, f, indent=2)

        self.console.print(f"\n[grey50]Batch report saved to: {report_file}[/grey50]")

        # Failures first, then anything the stop never reached: original order.
        retry_file = self._write_retry_list(
            base_dir, timestamp, failed_rows + not_attempted
        )
        if retry_file:
            self.console.print(
                f"[grey50]Runs still owed ({len(failed_rows) + len(not_attempted)}): "
                f"{retry_file}[/grey50]"
            )
            self.console.print(
                f"[grey50]Re-run them with: proprep --batch-replay "
                f"{template_file} --batch-list {retry_file}[/grey50]"
            )


def normalize_input_list(
    input_list: List[Union[str, Dict[str, str]]]
) -> List[Dict[str, str]]:
    """
    Convert every entry to a variable dict so runs have one uniform shape.

    A plain-text list carries only an identifier, which binds to the
    input_protein variable that TemplateConverter creates for the structure
    source. Validation and execution must normalize identically, or the two
    disagree about what a run supplies.

    Args:
        input_list: Entries from load_input_list (strings and/or dicts)

    Returns:
        List of variable dicts, blank entries dropped
    """
    normalized: List[Dict[str, str]] = []
    for item in input_list:
        if isinstance(item, str):
            item = item.strip()
            if item:
                normalized.append({"input_protein": item})
        elif isinstance(item, dict):
            normalized.append(item)
    return normalized


def load_input_list(input_file: str) -> List[Union[str, Dict[str, str]]]:
    """
    Load protein input list from file.

    Auto-detects format:
      - CSV with headers: first non-comment line contains commas -> returns list of dicts
      - Plain text: one protein identifier per line -> returns list of strings

    Args:
        input_file: Path to input file

    Returns:
        List of protein identifiers (strings) or variable dicts (from CSV)
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Filter out empty lines and comments, preserving order
    data_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith('#'):
            data_lines.append(stripped)

    if not data_lines:
        return []

    # Auto-detect CSV: commas in the first line, OR the first line looks
    # like a header (contains only identifier-safe chars, no spaces that
    # would indicate a PDB ID or file path).  Single-column CSVs produced
    # by the template system have a header like "enter_mutation_specification"
    # with no commas, so we also check for valid Python-identifier headers.
    first_line = data_lines[0]
    is_csv = ',' in first_line or (
        first_line.replace('_', '').isalnum()
        and not first_line[0].isdigit()
        and len(data_lines) > 1
        and first_line != data_lines[1]  # header differs from first data row
    )

    if is_csv:
        # Parse as CSV with headers
        csv_text = '\n'.join(data_lines)
        reader = csv.DictReader(StringIO(csv_text))
        return [dict(row) for row in reader]

    # Plain text: one protein per line (backward compatible)
    return data_lines


def run_batch_processing(
    template_file: str,
    input_list_file: str,
    base_dir: str = ".",
    continue_on_error: bool = True,
    web_override: Optional[bool] = None
) -> int:
    """
    Run batch processing from command line.

    Args:
        template_file: Path to template session file
        input_list_file: Path to file with protein identifiers
        base_dir: Base directory for outputs
        continue_on_error: Continue if one protein fails
        web_override: Force proprep-web (True) or plain-terminal (False) replay
            mode, overriding the mode recorded in the template. None means
            honor the template's recorded ``web_shell_mode``.

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    console = Console()

    try:
        # Load input list
        input_list = load_input_list(input_list_file)

        if not input_list:
            console.print(f"[red]No proteins found in input list: {input_list_file}[/red]")
            return 1

        # Create batch processor
        processor = BatchProcessor(console)

        # Process batch
        successful, failed = processor.process_batch(
            template_file=template_file,
            input_list=input_list,
            base_dir=base_dir,
            continue_on_error=continue_on_error,
            web_override=web_override
        )

        # Return exit code
        return 0 if failed == 0 else 1

    except Exception as e:
        console.print(f"[red]Batch processing error: {e}[/red]")
        traceback.print_exc()
        return 1
