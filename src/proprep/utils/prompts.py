"""
Centralized Rich Prompt Wrappers with Session Recording Support

This module provides unified wrapper functions for Rich's prompt classes
that integrate with ProPrep's session recording system.

IMPORTANT - Readline Conflict:
    The `readline` module conflicts with Rich's ANSI escape sequences, causing
    prompts to vanish when the user types a character and presses backspace.
    Do NOT import readline anywhere in ProPrep. See the note in
    md_prep/amber_workflow_components.py for details.

Usage:
    from proprep.utils.prompts import (
        prompt_with_context,
        confirm_with_context,
        int_prompt_with_context,
        float_prompt_with_context,
    )

    # Basic usage
    response = prompt_with_context(processor, "Enter value", default="x")

    # With full context for session recording
    response = prompt_with_context(
        processor,
        "Select option",
        choices=["1", "2", "3"],
        default="1",
        module="My Module",
        description="Select processing option",
        options_map={"1": "Option A", "2": "Option B", "3": "Option C"}
    )
"""

from typing import Any, Dict, List, Optional, Union

from rich import get_console
from rich.console import Console
from rich.prompt import Confirm, FloatPrompt, IntPrompt, Prompt
from rich.theme import Theme


# Rich renders a prompt's default value — the "(1a)" in
# "Enter your choice (1a):" — with its built-in ``prompt.default`` style,
# which is cyan. Cyan is low-contrast on a white background and clashes
# with the menu's blue option keys. Remap ``prompt.default`` to the
# menu's ``bold blue`` so the parenthesized default matches the option
# numbering it refers to. This is pushed onto Rich's GLOBAL console (the
# one ``Prompt.ask`` uses when no ``console=`` is passed), so it recolors
# every prompt default app-wide. It overrides only RENDERING, not input,
# so session recording (which hooks input, not styling) is unaffected.
get_console().push_theme(Theme({"prompt.default": "bold blue"}, inherit=True))


class NavigationException(Exception):
    """Raised when user requests to navigate backwards in a workflow."""
    pass


def prompt_with_context(
    processor,
    prompt: str,
    choices: Optional[List[str]] = None,
    default: Optional[str] = None,
    module: Optional[str] = None,
    description: Optional[str] = None,
    options_map: Optional[Dict[str, str]] = None,
    show_choices: Optional[bool] = None,
) -> str:
    """
    Ask the user a question with rich context for session recording.

    This wrapper adds a trailing space to prevent the backspace vanishing issue
    and enriches session recordings with contextual metadata.

    Args:
        processor: The PDBProcessor object (has session_manager attribute)
        prompt: The question to ask (e.g., "Enter your choice")
        choices: Valid answers (e.g., ["1", "2", "3"])
        default: Default answer if user just presses Enter
        module: Name of this part of ProPrep (e.g., "Main Menu", "Structure Loader")
        description: What this question is asking (e.g., "Select menu option")
        options_map: Dictionary mapping choices to human-readable descriptions
                    (e.g., {"1": "Structure Preparation", "2": "Analysis"})
        show_choices: Whether to show choices inline. If None, auto-determines:
                     False if options_map provided (assumes choices printed separately),
                     True otherwise.

    Returns:
        The user's answer as a string
    """
    # Auto-determine show_choices if not specified
    if show_choices is None:
        show_choices = (options_map is None)

    # Prompt.prompt_suffix has been modified at module level to prevent backspace issue.
    # Rich treats its `default` parameter's sentinel as `...` (Ellipsis), not None —
    # passing default=None makes None itself the return value on empty input, which
    # then breaks callers that do option_map[choice]. Omit the kwarg when not set.
    #
    # case_sensitive=False: Rich's Prompt.ask matches `choices` case-sensitively by
    # default, so a prompt offering ["a","b","n"] (or ["Y","N"]) silently re-asks
    # when the user types the other case. Callers already normalize with .lower()
    # AFTER the fact, but that runs too late — Rich rejects first. Matching
    # case-insensitively (and Rich then returns the canonical choice, e.g. "A"->"a")
    # makes every menu accept either case, which is what users expect.
    ask_kwargs = {"choices": choices, "show_choices": show_choices,
                  "case_sensitive": False}
    if default is not None:
        ask_kwargs["default"] = default
    response = Prompt.ask(prompt, **ask_kwargs)

    # Add rich context to the session recording
    if processor and hasattr(processor, 'session_manager'):
        recorder = processor.session_manager.recorder
        if recorder and recorder.recording and recorder.session_data.get("interactions"):
            # Get the last recorded interaction
            last_interaction = recorder.session_data["interactions"][-1]

            # Add contextual metadata
            if module:
                last_interaction["context"]["module"] = module
            if description:
                last_interaction["context"]["description"] = description
            if options_map:
                last_interaction["context"]["options_map"] = options_map
                if response in options_map:
                    last_interaction["context"]["option_label"] = options_map[response]

    return response


def confirm_with_context(
    processor,
    prompt: str,
    default: Optional[bool] = None,
    module: Optional[str] = None,
    description: Optional[str] = None,
    allow_back: bool = False,
) -> Union[bool, str]:
    """
    Ask the user a yes/no question with rich context for session recording.

    This wrapper adds a trailing space to prevent the backspace vanishing issue
    and enriches session recordings with contextual metadata.

    Args:
        processor: The PDBProcessor object (has session_manager attribute)
        prompt: The question to ask
        default: Default answer (True for yes, False for no, None for no default)
        module: Name of this part of ProPrep (e.g., "Structure Loader")
        description: What this question is asking
        allow_back: If True, user can type 'back' or 'b' to navigate backwards.
                   In this case, raises NavigationException.

    Returns:
        The user's answer (True/False)

    Raises:
        NavigationException: If user types 'back' and allow_back=True
    """
    if allow_back:
        # Use prompt-based approach to support back navigation
        enhanced_prompt = f"{prompt}\n[grey50]  (Type 'back' or 'b' to return to previous question)[/grey50]"
        console = processor.console if processor and hasattr(processor, 'console') else Console()

        while True:
            response = prompt_with_context(
                processor=processor,
                prompt=enhanced_prompt,
                choices=['y', 'n', 'yes', 'no', 'back', 'b'],
                default='y' if default else 'n',
                module=module,
                description=description
            )

            # Check for back navigation
            if response.lower() in ['back', 'b']:
                raise NavigationException("User requested to go back")

            # Convert to boolean
            if response.lower() in ['y', 'yes']:
                return True
            elif response.lower() in ['n', 'no']:
                return False
            else:
                console.print("[red]Please enter 'y' or 'n'[/red]")
    else:
        # Confirm.prompt_suffix has been modified at module level to prevent backspace issue.
        # See prompt_with_context note: Rich's no-default sentinel is `...`, not None.
        ask_kwargs = {}
        if default is not None:
            ask_kwargs["default"] = default
        response = Confirm.ask(prompt, **ask_kwargs)

        # Add rich context to the session recording
        if processor and hasattr(processor, 'session_manager'):
            recorder = processor.session_manager.recorder
            if recorder and recorder.recording and recorder.session_data.get("interactions"):
                last_interaction = recorder.session_data["interactions"][-1]

                if module:
                    last_interaction["context"]["module"] = module
                if description:
                    last_interaction["context"]["description"] = description
                last_interaction["context"]["option_label"] = "Yes" if response else "No"

        return response


def int_prompt_with_context(
    processor,
    prompt: str,
    default: Optional[int] = None,
    module: Optional[str] = None,
    description: Optional[str] = None,
    show_default: bool = True,
    options_map: Optional[Dict[str, str]] = None,
) -> int:
    """
    Ask the user for an integer value with rich context for session recording.

    This wrapper adds a trailing space to prevent the backspace vanishing issue
    and enriches session recordings with contextual metadata.

    Args:
        processor: The PDBProcessor object (has session_manager attribute)
        prompt: The question to ask
        default: Default value if user just presses Enter
        module: Name of this part of ProPrep
        description: What this question is asking
        show_default: Whether to show the default value in the prompt
        options_map: Dictionary mapping string indices to human-readable descriptions
                    (e.g., {"0": "structure1.pdb", "1": "structure2.pdb"})

    Returns:
        The user's answer as an integer
    """
    # IntPrompt.prompt_suffix has been modified at module level to prevent backspace issue.
    # See prompt_with_context note: Rich's no-default sentinel is `...`, not None.
    ask_kwargs = {"show_default": show_default}
    if default is not None:
        ask_kwargs["default"] = default
    response = IntPrompt.ask(prompt, **ask_kwargs)

    # Add rich context to the session recording
    if processor and hasattr(processor, 'session_manager'):
        recorder = processor.session_manager.recorder
        if recorder and recorder.recording and recorder.session_data.get("interactions"):
            last_interaction = recorder.session_data["interactions"][-1]

            if module:
                last_interaction["context"]["module"] = module
            if description:
                last_interaction["context"]["description"] = description
            if options_map:
                last_interaction["context"]["options_map"] = options_map
                response_str = str(response)
                if response_str in options_map:
                    last_interaction["context"]["option_label"] = options_map[response_str]

    return response


def float_prompt_with_context(
    processor,
    prompt: str,
    default: Optional[float] = None,
    module: Optional[str] = None,
    description: Optional[str] = None,
    show_default: bool = True,
) -> float:
    """
    Ask the user for a float value with rich context for session recording.

    This wrapper adds a trailing space to prevent the backspace vanishing issue
    and enriches session recordings with contextual metadata.

    Args:
        processor: The PDBProcessor object (has session_manager attribute)
        prompt: The question to ask
        default: Default value if user just presses Enter
        module: Name of this part of ProPrep
        description: What this question is asking
        show_default: Whether to show the default value in the prompt

    Returns:
        The user's answer as a float
    """
    # FloatPrompt.prompt_suffix has been modified at module level to prevent backspace issue.
    # See prompt_with_context note: Rich's no-default sentinel is `...`, not None.
    ask_kwargs = {"show_default": show_default}
    if default is not None:
        ask_kwargs["default"] = default
    response = FloatPrompt.ask(prompt, **ask_kwargs)

    # Add rich context to the session recording
    if processor and hasattr(processor, 'session_manager'):
        recorder = processor.session_manager.recorder
        if recorder and recorder.recording and recorder.session_data.get("interactions"):
            last_interaction = recorder.session_data["interactions"][-1]

            if module:
                last_interaction["context"]["module"] = module
            if description:
                last_interaction["context"]["description"] = description

    return response


def prompt_float_with_retry(
    processor,
    prompt: str,
    default: float,
    module: Optional[str] = None,
    description: Optional[str] = None,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
    max_retries: int = 3,
    allow_back: bool = False,
) -> float:
    """
    Prompt for a float value with validation and retry logic.

    This is a higher-level utility that wraps prompt_with_context with
    additional validation for numeric ranges.

    Args:
        processor: The PDBProcessor object
        prompt: The question to ask
        default: Default value if user presses Enter or exceeds retries
        module: Name of this part of ProPrep
        description: What this question is asking
        min_value: Optional minimum allowed value
        max_value: Optional maximum allowed value
        max_retries: Maximum number of retry attempts (default: 3)
        allow_back: If True, user can type 'back' to navigate backwards

    Returns:
        Float value from user or default

    Raises:
        NavigationException: If user types 'back' and allow_back=True
    """
    console = processor.console if processor and hasattr(processor, 'console') else Console()
    attempts = 0

    while attempts < max_retries:
        try:
            # Build prompt text with range info
            prompt_text = f"{prompt} ({default})"
            if min_value is not None or max_value is not None:
                if min_value is not None and max_value is not None:
                    prompt_text = f"{prompt} (range: {min_value}-{max_value}, default: {default})"
                elif min_value is not None:
                    prompt_text = f"{prompt} (min: {min_value}, default: {default})"
                else:
                    prompt_text = f"{prompt} (max: {max_value}, default: {default})"

            # Add back hint if enabled
            if allow_back:
                prompt_text += "\n[grey50]  (Type 'back' or 'b' to return to previous question)[/grey50]"

            response = prompt_with_context(
                processor=processor,
                prompt=prompt_text,
                default=str(default),
                module=module,
                description=description
            ).strip()

            # Empty response = use default
            if not response:
                return default

            # Check for back navigation
            if allow_back and response.lower() in ['back', 'b']:
                raise NavigationException("User requested to go back")

            # Special commands
            if response.lower() in ['default', 'd']:
                console.print(f"[grey50]Using default: {default}[/grey50]")
                return default

            # Parse float
            value = float(response)

            # Validate range
            if min_value is not None and value < min_value:
                console.print(f"[red]Value must be at least {min_value}[/red]")
                attempts += 1
                continue
            if max_value is not None and value > max_value:
                console.print(f"[red]Value must be at most {max_value}[/red]")
                attempts += 1
                continue

            return value

        except ValueError:
            console.print("[red]Please enter a valid number[/red]")
            attempts += 1

    # Exceeded retries, use default
    console.print(f"[yellow]Using default value: {default}[/yellow]")
    return default


def prompt_int_with_retry(
    processor,
    prompt: str,
    default: int,
    module: Optional[str] = None,
    description: Optional[str] = None,
    min_value: Optional[int] = None,
    max_value: Optional[int] = None,
    max_retries: int = 3,
    allow_back: bool = False,
) -> int:
    """
    Prompt for an integer value with validation and retry logic.

    This is a higher-level utility that wraps prompt_with_context with
    additional validation for numeric ranges.

    Args:
        processor: The PDBProcessor object
        prompt: The question to ask
        default: Default value if user presses Enter or exceeds retries
        module: Name of this part of ProPrep
        description: What this question is asking
        min_value: Optional minimum allowed value
        max_value: Optional maximum allowed value
        max_retries: Maximum number of retry attempts (default: 3)
        allow_back: If True, user can type 'back' to navigate backwards

    Returns:
        Integer value from user or default

    Raises:
        NavigationException: If user types 'back' and allow_back=True
    """
    console = processor.console if processor and hasattr(processor, 'console') else Console()
    attempts = 0

    while attempts < max_retries:
        try:
            # Build prompt text with range info
            prompt_text = f"{prompt} ({default})"
            if min_value is not None or max_value is not None:
                if min_value is not None and max_value is not None:
                    prompt_text = f"{prompt} (range: {min_value}-{max_value}, default: {default})"
                elif min_value is not None:
                    prompt_text = f"{prompt} (min: {min_value}, default: {default})"
                else:
                    prompt_text = f"{prompt} (max: {max_value}, default: {default})"

            # Add back hint if enabled
            if allow_back:
                prompt_text += "\n[grey50]  (Type 'back' or 'b' to return to previous question)[/grey50]"

            response = prompt_with_context(
                processor=processor,
                prompt=prompt_text,
                default=str(default),
                module=module,
                description=description
            ).strip()

            # Empty response = use default
            if not response:
                return default

            # Check for back navigation
            if allow_back and response.lower() in ['back', 'b']:
                raise NavigationException("User requested to go back")

            # Special commands
            if response.lower() in ['default', 'd']:
                console.print(f"[grey50]Using default: {default}[/grey50]")
                return default

            # Parse integer
            value = int(response)

            # Validate range
            if min_value is not None and value < min_value:
                console.print(f"[red]Value must be at least {min_value}[/red]")
                attempts += 1
                continue
            if max_value is not None and value > max_value:
                console.print(f"[red]Value must be at most {max_value}[/red]")
                attempts += 1
                continue

            return value

        except ValueError:
            console.print("[red]Please enter a valid integer[/red]")
            attempts += 1

    # Exceeded retries, use default
    console.print(f"[yellow]Using default value: {default}[/yellow]")
    return default
