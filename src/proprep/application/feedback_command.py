"""
Feedback Command for ProPrep

Allows users to submit feedback, bug reports, and feature requests
directly to the developer via GitHub Issues.
"""

import json
import platform
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any

import requests
from rich.panel import Panel
from rich.prompt import Prompt, Confirm
from rich.table import Table

from .processor_command import MenuCommand
from proprep.utils.prompts import prompt_with_context, confirm_with_context


class FeedbackCommand(MenuCommand):
    """Command to collect and submit user feedback."""

    def __init__(self, processor):
        super().__init__(processor, "User Feedback", "Submit feedback or report issues")

    def execute(self) -> None:
        """Display feedback form and collect user input."""
        self.enter_menu()

        try:
            self._show_feedback_intro()

            # Collect feedback information
            feedback_data = self._collect_feedback()

            if feedback_data is None:
                # User cancelled
                self.console.print("[yellow]Feedback cancelled.[/yellow]")
                return

            # Show preview and confirm
            if not self._confirm_submission(feedback_data):
                self.console.print("[yellow]Feedback not submitted.[/yellow]")
                return

            # Submit the feedback
            success, message = self._submit_feedback(feedback_data)

            if success:
                self.console.print(f"\n[green]✓ {message}[/green]")
            else:
                self.console.print(f"\n[red]✗ {message}[/red]")
                self._save_feedback_locally(feedback_data)

        finally:
            self.exit_menu()

    def _show_feedback_intro(self):
        """Display introduction panel."""
        intro_text = """[bold]Share Your Experience with ProPrep[/bold]

Your feedback helps improve ProPrep for everyone! You can:

  • Report bugs or issues
  • Request new features
  • Share suggestions for improvements
  • Ask questions about functionality

All feedback is submitted as a GitHub Issue (anonymously if you prefer).

[bold yellow]First Time Setup:[/bold yellow]
Before submitting feedback, you need to configure GitHub integration:
  1. Go to: [bright_blue]https://github.com/settings/tokens[/bright_blue]
  2. Click "Generate new token (classic)"
  3. Give it a name (e.g., "ProPrep Feedback")
  4. Check ONLY: [bold]public_repo[/bold] scope
  5. Generate and copy the token
  6. In ProPrep: Main Menu → [bold]p. Preferences[/bold] → [bold]2. Configure GitHub[/bold]
  7. Paste your token when prompted

[grey50]This is a one-time setup. Your token is stored securely in ~/.proprep/settings.json[/grey50]
"""
        self.console.print(Panel(intro_text, border_style="bright_blue", padding=(1, 2)))
        self.console.print()

    def _collect_feedback(self) -> Optional[Dict[str, Any]]:
        """
        Collect feedback information from user.

        Returns:
            Dictionary with feedback data, or None if cancelled
        """
        feedback_data = {}

        # 1. Feedback Type
        self.console.print("[bold]1. Feedback Type[/bold]")
        feedback_types = {
            "1": "Bug Report",
            "2": "Feature Request",
            "3": "Suggestion",
            "4": "Question",
            "5": "Other"
        }

        # Show options
        for key, label in feedback_types.items():
            self.console.print(f"  {key}. {label}")

        type_choice = prompt_with_context(
            self.processor,
            "Select feedback type",
            choices=list(feedback_types.keys()),
            default="3",
            module="Feedback",
            description="Feedback type selection",
            options_map=feedback_types
        )

        feedback_data["type"] = feedback_types[type_choice]
        self.console.print()

        # 2. Related Module (Optional)
        self.console.print("[bold]2. Related Module (Optional)[/bold]")
        self.console.print("[grey50]Which ProPrep module is this feedback about? (or 'none' for general)[/grey50]")

        # Use the same module list as main menu (curated, non-deprecated)
        curated_modules = [
            "PDB Loader",
            "Redox Site Detector",
            "Homology Searcher",
            "Structure Aligner",
            "PDB Filter",
            "Amino Acid Mutator",
            "Structure Fixer",
            "Protonation State Analyzer",
            "Redox Site Preparer",
            "MD Restraint Manager",
            "Force Field Parameterizer",
            "Topology Generator",
            "Molecular Dynamics Manager",
            "QM/MM Preparator",
            "Electron Transfer Analysis",
            "Structure Viewer",
        ]

        # Build module options (only include modules that exist in registry)
        module_options = {"0": "General / Not module-specific"}
        module_counter = 1
        for module_name in curated_modules:
            if module_name in self.processor.registry.modules:
                module_options[str(module_counter)] = module_name
                module_counter += 1

        # Show all options (sorted by number)
        for key in sorted(module_options.keys(), key=lambda x: int(x)):
            self.console.print(f"  {key}. {module_options[key]}")

        module_choice = prompt_with_context(
            self.processor,
            "Module number (or type module name)",
            default="0",
            module="Feedback",
            description="Related module",
            options_map=module_options
        )

        # Handle both number and name input
        if module_choice in module_options:
            feedback_data["module"] = module_options[module_choice]
        elif module_choice in curated_modules:
            feedback_data["module"] = module_choice
        else:
            feedback_data["module"] = "General"

        self.console.print()

        # 3. Email (Optional)
        self.console.print("[bold]3. Your Email (Optional)[/bold]")
        self.console.print("[grey50]Provide your email if you'd like a response (leave blank for anonymous)[/grey50]")

        email = prompt_with_context(
            self.processor,
            "Email address",
            default="",
            module="Feedback",
            description="Enter email address"
        )
        feedback_data["email"] = email if email.strip() else "Anonymous"
        self.console.print()

        # 4. Summary
        self.console.print("[bold]4. Brief Summary[/bold]")
        self.console.print("[grey50]One-line description of your feedback[/grey50]")

        summary = prompt_with_context(
            self.processor,
            "Summary",
            module="Feedback",
            description="Enter feedback summary"
        )
        if not summary.strip():
            self.console.print("[red]Summary is required.[/red]")
            return None
        feedback_data["summary"] = summary
        self.console.print()

        # 5. Detailed Description
        self.console.print("[bold]5. Detailed Description[/bold]")
        self.console.print("[grey50]Describe your feedback in detail. For bugs, include:[/grey50]")
        self.console.print("[grey50]  - What you were trying to do[/grey50]")
        self.console.print("[grey50]  - What happened instead[/grey50]")
        self.console.print("[grey50]  - Any error messages[/grey50]")
        self.console.print("[grey50]Press Enter twice when done (empty line to finish):[/grey50]")
        self.console.print()

        description_lines = []
        while True:
            line = prompt_with_context(
                self.processor, "",
                default="",
                module="Feedback",
                description="Feedback description line (blank to finish)",
            )
            if not line.strip():
                break
            description_lines.append(line)

        feedback_data["description"] = "\n".join(description_lines)
        if not feedback_data["description"].strip():
            self.console.print("[red]Description is required.[/red]")
            return None
        self.console.print()

        # 6. Session Context (Optional, encrypted to the maintainer)
        from proprep.utils.feedback_crypto import crypto_available, public_key_configured
        can_encrypt = crypto_available() and public_key_configured()

        if can_encrypt:
            self.console.print("[bold]6. Include Session Context?[/bold]")
            self.console.print(
                "[grey50]Attaches your full session history — every prompt and the answer\n"
                "you typed (PDB IDs, file paths, menu selections) — to help diagnose issues.[/grey50]"
            )
            self.console.print(
                "[grey50]It is encrypted to the maintainer's key before it is attached, so even\n"
                "though the issue itself is public, only the maintainer can read it.[/grey50]"
            )
            # Off by default: attaching the transcript is a deliberate choice.
            include_session = confirm_with_context(
                self.processor,
                "Include session context?",
                default=False,
                module="Feedback",
                description="Include session context in feedback"
            )
        else:
            # No encryption available/configured: never attach the transcript in
            # the clear, so the option is simply not offered.
            include_session = False

        if include_session:
            feedback_data["session_context"] = self._get_session_context()
        else:
            feedback_data["session_context"] = None

        self.console.print()

        # Add system info
        feedback_data["system_info"] = self._get_system_info()
        feedback_data["timestamp"] = datetime.now().isoformat()

        return feedback_data

    def _get_session_context(self) -> Dict[str, Any]:
        """
        Extract relevant session context for debugging.

        Returns:
            Dictionary with session context
        """
        context = {}

        # Get recent interactions if session recording is active
        if hasattr(self.processor, 'session_manager'):
            recorder = self.processor.session_manager.recorder
            if recorder and recorder.recording:
                # Include full session history for diagnosis
                interactions = recorder.session_data.get("interactions", [])

                # Simplify interaction data for readability
                context["session_interactions"] = [
                    {
                        "type": i.get("type"),
                        "prompt": i.get("prompt"),
                        "response": i.get("response"),
                        "context": i.get("context", {})
                    }
                    for i in interactions
                ]

        # Get current workspace state (keys only, not full data)
        workspace_keys = list(self.processor.workspace._data.keys())
        context["workspace_keys"] = workspace_keys

        # Get menu mode
        menu_mode = self.processor.workspace.get("menu_mode", "full-menu")
        context["menu_mode"] = menu_mode

        return context

    def _get_system_info(self) -> Dict[str, str]:
        """
        Get system information for debugging.

        Returns:
            Dictionary with system information
        """
        try:
            from importlib.metadata import version as _pkg_version
            proprep_version = _pkg_version("proprep")
        except Exception:
            proprep_version = "unknown"

        return {
            "proprep_version": proprep_version,
            "python_version": sys.version,
            "platform": platform.platform(),
            "os": platform.system(),
        }

    def _confirm_submission(self, feedback_data: Dict[str, Any]) -> bool:
        """
        Show preview and confirm submission.

        Args:
            feedback_data: Collected feedback data

        Returns:
            True if user confirms, False otherwise
        """
        self.console.print("\n[bold bright_blue]═══ Feedback Preview ═══[/bold bright_blue]\n")

        table = Table(show_header=False, box=None, padding=(0, 2))
        table.add_column("Field", style="bright_blue", width=20)
        table.add_column("Value")

        table.add_row("Type:", f"[yellow]{feedback_data['type']}[/yellow]")
        table.add_row("Module:", feedback_data['module'])
        table.add_row("Email:", feedback_data['email'])
        table.add_row("Summary:", f"[bold]{feedback_data['summary']}[/bold]")
        table.add_row("", "")  # Spacer
        table.add_row("Description:", "")

        self.console.print(table)

        # Print description with indentation
        for line in feedback_data['description'].split('\n'):
            self.console.print(f"  {line}")

        self.console.print()

        if feedback_data.get('session_context'):
            self.console.print("[grey50]✓ Session context will be included[/grey50]")

        self.console.print()

        return confirm_with_context(
            self.processor,
            "[bold]Submit this feedback to GitHub Issues?[/bold]",
            default=True,
            module="Feedback",
            description="Confirm submit feedback"
        )

    def _submit_feedback(self, feedback_data: Dict[str, Any]) -> tuple[bool, str]:
        """
        Submit feedback to GitHub Issues.

        Args:
            feedback_data: Collected feedback data

        Returns:
            Tuple of (success: bool, message: str)
        """
        # Get GitHub settings
        from proprep.utils.settings_manager import SettingsManager

        settings = SettingsManager()
        github_token = settings.get_github_token()
        github_repo = settings.get_github_repo()

        if not github_token or not github_repo:
            return (False,
                   "GitHub integration not configured. Please run Preferences menu to set up GitHub token.")

        # Format the issue
        issue_title = f"[{feedback_data['type']}] {feedback_data['summary']}"

        issue_body = self._format_issue_body(feedback_data)

        # Determine labels based on type
        labels = self._get_labels_for_type(feedback_data['type'])

        # Create GitHub issue
        try:
            url = f"https://api.github.com/repos/{github_repo}/issues"
            headers = {
                "Authorization": f"token {github_token}",
                "Accept": "application/vnd.github.v3+json"
            }
            payload = {
                "title": issue_title,
                "body": issue_body,
                "labels": labels
            }

            response = requests.post(url, headers=headers, json=payload, timeout=10)

            if response.status_code == 201:
                issue_data = response.json()
                issue_url = issue_data.get("html_url", "")
                return (True, f"Feedback submitted successfully! Issue: {issue_url}")
            else:
                error_msg = response.json().get("message", "Unknown error")
                return (False, f"Failed to create issue: {error_msg}")

        except requests.exceptions.Timeout:
            return (False, "Request timed out. Please check your internet connection.")
        except requests.exceptions.RequestException as e:
            return (False, f"Network error: {str(e)}")
        except Exception as e:
            return (False, f"Unexpected error: {str(e)}")

    def _format_issue_body(self, feedback_data: Dict[str, Any]) -> str:
        """
        Format feedback data as GitHub issue body.

        Args:
            feedback_data: Collected feedback data

        Returns:
            Formatted markdown string
        """
        body_parts = []

        # Header
        body_parts.append(f"**Feedback Type:** {feedback_data['type']}")
        body_parts.append(f"**Module:** {feedback_data['module']}")
        body_parts.append(f"**Submitted:** {feedback_data['timestamp']}")

        if feedback_data['email'] != "Anonymous":
            body_parts.append(f"**Contact:** {feedback_data['email']}")

        body_parts.append("")  # Blank line

        # Description
        body_parts.append("## Description")
        body_parts.append("")
        body_parts.append(feedback_data['description'])
        body_parts.append("")

        # System Info
        body_parts.append("## System Information")
        body_parts.append("")
        body_parts.append("```")
        for key, value in feedback_data['system_info'].items():
            body_parts.append(f"{key}: {value}")
        body_parts.append("```")
        body_parts.append("")

        # Session Context, encrypted to the maintainer key so it can sit inside
        # a public issue without exposing the transcript to anyone else.
        if feedback_data.get('session_context'):
            from proprep.utils.feedback_crypto import encrypt_for_maintainer
            plaintext = json.dumps(feedback_data['session_context'], indent=2)
            armored = encrypt_for_maintainer(plaintext)
            body_parts.append("## Session Context (encrypted)")
            body_parts.append("")
            if armored:
                body_parts.append(
                    "The reporter attached session context, encrypted to the "
                    "maintainer's key. Decrypt with `proprep --decrypt-feedback`."
                )
                body_parts.append("")
                body_parts.append("```")
                body_parts.append(armored)
                body_parts.append("```")
            else:
                # Encryption unexpectedly unavailable at format time: omit rather
                # than fall back to posting the transcript in the clear.
                body_parts.append(
                    "_Session context was requested but encryption was "
                    "unavailable, so it was omitted rather than posted in the clear._"
                )
            body_parts.append("")

        # Footer
        body_parts.append("---")
        body_parts.append("*Submitted via ProPrep Feedback Command*")

        return "\n".join(body_parts)

    def _get_labels_for_type(self, feedback_type: str) -> list[str]:
        """
        Get GitHub labels based on feedback type.

        Args:
            feedback_type: Type of feedback

        Returns:
            List of label strings
        """
        label_map = {
            "Bug Report": ["bug", "user-reported"],
            "Feature Request": ["enhancement", "user-requested"],
            "Suggestion": ["enhancement", "user-feedback"],
            "Question": ["question", "user-feedback"],
            "Other": ["user-feedback"]
        }

        return label_map.get(feedback_type, ["user-feedback"])

    def _save_feedback_locally(self, feedback_data: Dict[str, Any]):
        """
        Save feedback to local file as backup.

        Args:
            feedback_data: Collected feedback data
        """
        # Create feedback directory in ~/.proprep/feedback
        feedback_dir = Path.home() / ".proprep" / "feedback"
        feedback_dir.mkdir(parents=True, exist_ok=True)

        # Generate filename with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = feedback_dir / f"feedback_{timestamp}.json"

        # Save feedback
        try:
            with open(filename, 'w') as f:
                json.dump(feedback_data, f, indent=2)

            self.console.print(f"\n[bright_blue]Feedback saved locally: {filename}[/bright_blue]")
            self.console.print("[grey50]You can manually submit this later or email it to the developer.[/grey50]")
        except Exception as e:
            self.console.print(f"\n[red]Failed to save feedback locally: {e}[/red]")
