"""
Session Recording and Replay System for MPSA

This module implements a comprehensive session recording system that captures all user
interactions (menu selections and inputs) and allows for automated replay of sessions.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, TextIO, Union
from io import StringIO
import threading
import queue

from rich.prompt import Prompt, Confirm, IntPrompt, FloatPrompt


def detect_and_recover_json_corruption(file_path: str, backup: bool = True) -> Optional[Dict[str, Any]]:
    """
    Detect and automatically recover corrupted JSON session files.

    This function handles common corruption patterns:
    - Extra data after valid JSON (missing truncate() bug)
    - Partial JSON from interrupted writes
    - Malformed closing braces

    Args:
        file_path: Path to potentially corrupted JSON file
        backup: If True, create backup of original corrupted file

    Returns:
        Recovered session data dict, or None if unrecoverable

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Session file not found: {file_path}")

    # Try normal load first
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"[yellow]⚠ JSON corruption detected in {os.path.basename(file_path)}[/yellow]")
        print(f"[yellow]  Error: {e.msg} at line {e.lineno}, column {e.colno}[/yellow]")
        print(f"[yellow]  Attempting automatic recovery...[/yellow]")

    # Create backup if requested
    if backup:
        backup_path = file_path + '.corrupted_backup'
        try:
            import shutil
            shutil.copy2(file_path, backup_path)
            print(f"[grey50]  Backup created: {os.path.basename(backup_path)}[/grey50]")
        except Exception as e:
            print(f"[yellow]  Warning: Could not create backup: {e}[/yellow]")

    # Attempt recovery strategies
    recovered_data = None

    # Strategy 1: Read up to the error position and try to find valid JSON
    try:
        with open(file_path, 'r') as f:
            content = f.read()

        # Try to parse JSON by progressively removing content from the end
        # This handles "extra data" errors where valid JSON is followed by garbage
        lines = content.split('\n')

        for end_line in range(len(lines), 0, -1):
            test_content = '\n'.join(lines[:end_line]).rstrip()

            # Ensure it ends with a proper closing brace
            if not test_content.endswith('}'):
                continue

            try:
                recovered_data = json.loads(test_content)
                print(f"[green]✓ Recovery successful! Found valid JSON ending at line {end_line}/{len(lines)}[/green]")
                break
            except json.JSONDecodeError:
                continue

        if recovered_data is None:
            # Strategy 2: Try to find the last complete '}' and parse up to there
            brace_positions = [i for i, char in enumerate(content) if char == '}']

            for pos in reversed(brace_positions):
                test_content = content[:pos+1]
                try:
                    recovered_data = json.loads(test_content)
                    print(f"[green]✓ Recovery successful! Used content up to character {pos}[/green]")
                    break
                except json.JSONDecodeError:
                    continue

    except Exception as e:
        print(f"[red]✗ Recovery failed: {e}[/red]")
        return None

    if recovered_data:
        # Write recovered data back to file
        try:
            with open(file_path, 'w') as f:
                json.dump(recovered_data, f, indent=2)
            print(f"[green]✓ Recovered file saved successfully[/green]")

            # Validate the recovery
            interactions = recovered_data.get('interactions', [])
            metadata = recovered_data.get('metadata', {})
            print(f"[green]  Recovered {len(interactions)} interactions[/green]")
            if metadata:
                print(f"[grey50]  Metadata: {', '.join(metadata.keys())}[/grey50]")

        except Exception as e:
            print(f"[red]✗ Could not save recovered file: {e}[/red]")
            return None
    else:
        print(f"[red]✗ Could not recover session file[/red]")
        print(f"[yellow]  You may need to manually edit the file or restore from backup[/yellow]")

    return recovered_data


def safe_load_session_file(file_path: str, auto_recover: bool = True) -> Optional[Dict[str, Any]]:
    """
    Safely load a session file with automatic corruption detection and recovery.

    Args:
        file_path: Path to session file
        auto_recover: If True, automatically attempt recovery on corruption

    Returns:
        Session data dict, or None if file cannot be loaded/recovered
    """
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError:
        if auto_recover:
            return detect_and_recover_json_corruption(file_path, backup=True)
        else:
            raise
    except FileNotFoundError:
        return None


class SessionRecorder:
    """
    Records all user interactions during an MPSA session for later replay.

    This class intercepts user inputs and menu selections, recording them
    to a session log file that can be used for automated replay.

    Enhanced features:
    - Load and continue existing sessions
    - Truncate sessions at specific interaction indices
    - Support for rich context (module names, descriptions, option labels)
    """

    def __init__(self, record_file: Optional[str] = None):
        """
        Initialize the session recorder.

        Args:
            record_file: Path to the file where session will be recorded
        """
        self.record_file = record_file
        self.recording = False
        self.session_data = {
            "version": "1.1",  # Bumped version for rich context support
            "start_time": None,
            "end_time": None,
            "interactions": [],
            "metadata": {}
        }
        self._interaction_count = 0
        self._file_handle: Optional[TextIO] = None
        
    def start_recording(self, metadata: Dict[str, Any] = None):
        """Start recording a session."""
        if self.recording:
            return
            
        self.recording = True
        self.session_data["start_time"] = datetime.now().isoformat()
        
        if metadata:
            self.session_data["metadata"] = metadata
            
        # Create or open the recording file
        if self.record_file:
            self._file_handle = open(self.record_file, 'w')
            # Write initial session data
            self._write_header()
    
    def stop_recording(self):
        """Stop recording and finalize the session file."""
        if not self.recording:
            return

        self.recording = False
        self.session_data["end_time"] = datetime.now().isoformat()

        if self._file_handle:
            # Write the complete session data
            self._file_handle.seek(0)
            json.dump(self.session_data, self._file_handle, indent=2)
            self._file_handle.truncate()  # CRITICAL: Remove any trailing data
            self._file_handle.close()
            self._file_handle = None

    def load_existing_session(self, session_file: str):
        """
        Load an existing session to continue recording from it.

        This allows resuming a session and appending new interactions.

        Args:
            session_file: Path to existing session file
        """
        existing_data = safe_load_session_file(session_file, auto_recover=True)
        if existing_data is None:
            raise ValueError(f"Could not load or recover session file: {session_file}")

        # Inherit existing session data
        self.session_data = existing_data.copy()

        # Update interaction count to continue from where it left off
        self._interaction_count = len(self.session_data.get("interactions", []))

        # Update metadata to indicate continuation
        if "metadata" not in self.session_data:
            self.session_data["metadata"] = {}

        self.session_data["metadata"]["continued_from"] = session_file
        self.session_data["metadata"]["continuation_time"] = datetime.now().isoformat()

    def truncate_at(self, index: int):
        """
        Truncate the session at a specific interaction index.

        Removes all interactions after the specified index.
        Useful for editing sessions and discarding invalid future interactions.

        Args:
            index: Interaction index to truncate at (keeps 0 to index-1, removes index onwards)
        """
        interactions = self.session_data.get("interactions", [])

        if index < 0 or index > len(interactions):
            raise ValueError(f"Invalid truncation index: {index} (session has {len(interactions)} interactions)")

        # Keep only interactions before the index
        self.session_data["interactions"] = interactions[:index]
        self._interaction_count = len(self.session_data["interactions"])

        # Add metadata about truncation
        if "metadata" not in self.session_data:
            self.session_data["metadata"] = {}

        self.session_data["metadata"]["truncated_at"] = index
        self.session_data["metadata"]["truncation_time"] = datetime.now().isoformat()
        self.session_data["metadata"]["removed_interactions"] = len(interactions) - index
            
    def record_interaction(self, interaction_type: str, prompt: str, 
                          response: str, choices: List[str] = None,
                          context: Dict[str, Any] = None):
        """
        Record a single user interaction.
        
        Args:
            interaction_type: Type of interaction (menu_choice, text_input, confirm)
            prompt: The prompt shown to the user
            response: The user's response
            choices: Available choices (for menu selections)
            context: Additional context about the interaction
        """
        if not self.recording:
            return
            
        interaction = {
            "index": self._interaction_count,
            "timestamp": datetime.now().isoformat(),
            "type": interaction_type,
            "prompt": prompt,
            "response": response,
            "choices": choices,
            "context": context or {}
        }
        
        self._interaction_count += 1
        self.session_data["interactions"].append(interaction)
        
        # Write incrementally if file is open
        if self._file_handle:
            self._write_incremental()
            
    def annotate_last_interaction(self, extra_context: Dict[str, Any]):
        """
        Merge extra_context into the most recently recorded interaction's
        context dict, and flush to disk.

        Used by callers (e.g. file-picker menus) to attach resolver hints —
        the basename of the selected file, an external identifier, etc. — so
        replay can re-resolve responses that would otherwise drift when the
        underlying option list changes between recording and replay.
        """
        interactions = self.session_data.get("interactions", [])
        if not interactions:
            return
        last = interactions[-1]
        last.setdefault("context", {}).update(extra_context)
        if self._file_handle:
            self._write_incremental()

    def _write_header(self):
        """Write the initial header to the file."""
        if self._file_handle:
            # Write a placeholder that we'll update later
            json.dump(self.session_data, self._file_handle, indent=2)
            self._file_handle.flush()
            
    def _write_incremental(self):
        """Update the file with the latest session data."""
        if self._file_handle:
            self._file_handle.seek(0)
            json.dump(self.session_data, self._file_handle, indent=2)
            self._file_handle.truncate()
            self._file_handle.flush()


class SessionReplayer:
    """
    Replays recorded MPSA sessions by providing recorded responses.

    This class reads a session log file and provides the recorded
    responses when the application requests user input.

    Supports variable substitution for template-based batch processing.
    """

    def __init__(self, replay_file: str, replay_delay: float = 0.0, variables: Optional[Dict[str, str]] = None,
                 strict_variables: bool = True):
        """
        Initialize the session replayer.

        Args:
            replay_file: Path to the session log file to replay
            replay_delay: Delay in seconds between replay actions (for demos)
            variables: Dictionary of variable values for template substitution
            strict_variables: If True (batch replay), a required template
                variable with no provided value is a hard error. If False
                (interactive resume/continue), leave such placeholders
                unresolved so replay falls through to live input at that step
                instead of crashing.
        """
        self.replay_file = replay_file
        self.replay_delay = replay_delay
        self.variables = variables or {}
        self.strict_variables = strict_variables
        self.session_data = None
        self.interaction_index = 0
        self.replaying = False
        self.is_template = False
        # The most recently replayed interaction object (or None if the last
        # get_next_response call did not return a recorded response). Callers
        # use this to read recorded context (e.g. the originally selected
        # filename) so they can re-resolve indices that have shifted since
        # the log was recorded.
        self.last_returned_interaction: Optional[Dict[str, Any]] = None
        # True once the log and the run have parted company, so the divergence
        # is announced once rather than at every subsequent prompt.
        self._diverged = False
        self._load_session()
        
    def _load_session(self):
        """Load the session data from the replay file."""
        self.session_data = safe_load_session_file(self.replay_file, auto_recover=True)
        if self.session_data is None:
            raise ValueError(f"Could not load or recover session file: {self.replay_file}")

        # Check if this is a template
        self.is_template = self.session_data.get("template", False)

        if self.is_template:
            # Apply variable substitution
            self._substitute_variables()

    def _substitute_variables(self):
        """
        Substitute template variables with actual values.

        Replaces {{ variable_name }} in interaction responses with actual values.
        """
        import re

        # Get template variable definitions
        template_vars = self.session_data.get("template_variables", {})

        # Check that all required variables are provided. In strict (batch)
        # mode this is a hard error; in lenient (interactive resume) mode we
        # leave the placeholder unresolved so replay falls through to live
        # input when it reaches that interaction (see get_next_response).
        for var_name, var_def in template_vars.items():
            if var_def.get("required", True) and var_name not in self.variables:
                if self.strict_variables:
                    raise ValueError(
                        f"Required template variable '{var_name}' not provided. "
                        f"Description: {var_def.get('description', 'N/A')}"
                    )

        # Substitute variables in all interactions
        interactions = self.session_data.get("interactions", [])
        for interaction in interactions:
            response = interaction.get("response", "")
            if isinstance(response, str):
                # Find all {{ variable }} patterns
                pattern = r'\{\{\s*(\w+)\s*\}\}'
                matches = re.findall(pattern, response)

                # Substitute each variable
                for var_name in matches:
                    if var_name in self.variables:
                        placeholder = f"{{{{ {var_name} }}}}"
                        interaction["response"] = interaction["response"].replace(
                            placeholder,
                            str(self.variables[var_name])
                        )

    def start_replay(self):
        """Start replaying the session."""
        self.replaying = True
        self.interaction_index = 0
        
    def stop_replay(self):
        """Stop replaying the session."""
        self.replaying = False

    def truncate_at(self, index: int):
        """
        Truncate the session replay at a specific interaction index.

        This limits replay to only interactions before the specified index.

        Args:
            index: Interaction index to truncate at
        """
        if not self.session_data:
            return

        interactions = self.session_data.get("interactions", [])

        if index < 0 or index > len(interactions):
            raise ValueError(f"Invalid truncation index: {index} (session has {len(interactions)} interactions)")

        # Truncate interactions
        self.session_data["interactions"] = interactions[:index]

    def get_next_response(self, interaction_type: str, prompt: str) -> Optional[str]:
        """
        Get the next recorded response for the given prompt.

        Args:
            interaction_type: Type of interaction expected (prompt, confirm, input)
            prompt: The prompt being shown

        Returns:
            The recorded response or None if not found
        """
        if not self.replaying or not self.session_data:
            self.last_returned_interaction = None
            return None

        # Add delay for demo mode if configured
        if self.replay_delay > 0:
            import time
            time.sleep(self.replay_delay)

        interactions = self.session_data.get("interactions", [])

        # STRICT by default: a recorded answer is used only if it is the very
        # NEXT unconsumed interaction and it matches this prompt exactly.
        #
        # The alternative — scanning forward for a match — was meant to tolerate
        # a recorded prompt the current run does not ask. It cannot tell that
        # apart from the same question asked at a DIFFERENT point in the
        # workflow, and then it leaps. In the run that prompted this, a hydrogen
        # editor prompt now asked during step 8 matched its recording from step
        # 12, sixty-five interactions later; the scan consumed three "Select
        # action" checklist decisions on the way and replay went on to run step
        # 13 while the checklist was sitting at step 9.
        #
        # Strict matching also gives divergence recovery for free. A prompt that
        # does not match leaves the position untouched and falls through to live
        # input, so once the run asks something the head DOES match — the
        # checklist's next action, say — replay simply resumes there.
        if self.interaction_index >= len(interactions):
            self.last_returned_interaction = None
            return None

        interaction = interactions[self.interaction_index]
        if (interaction["type"] != interaction_type or
                interaction["prompt"] != prompt):
            if not self._diverged:
                self._diverged = True
                self._report_divergence(interaction, interaction_type, prompt)
            self.last_returned_interaction = None
            return None

        self.interaction_index += 1
        response = interaction["response"]
        # A response that still holds an unresolved {{ variable }} placeholder
        # (lenient resume of a template with no value for that variable) is not
        # a usable recorded answer. Fall through to live input so the user
        # supplies it in context; in hybrid mode the value they type is then
        # recorded into the session.
        import re
        if isinstance(response, str) and re.search(r'\{\{\s*\w+\s*\}\}', response):
            self.last_returned_interaction = None
            return None

        if self._diverged:
            self._diverged = False
            print("[Replay resynchronised — continuing from the recorded session]")

        self.last_returned_interaction = interaction
        return response

    def _report_divergence(self, expected: Dict[str, Any],
                           asked_type: str, asked_prompt: str) -> None:
        """Say once that the log and the run have parted company."""
        exp_prompt = (expected.get("prompt") or "").splitlines()
        exp_prompt = exp_prompt[0] if exp_prompt else ""
        asked = (asked_prompt or "").splitlines()
        asked = asked[0] if asked else ""
        print(
            f"\n[Replay diverged at interaction {self.interaction_index}]\n"
            f"    recorded next: ({expected.get('type')}) {exp_prompt!r}"
            f" -> {expected.get('response')!r}\n"
            f"    asked now:     ({asked_type}) {asked!r}\n"
            f"    Answer this one yourself; replay resumes when the recorded "
            f"question comes up again."
        )
        
    def has_more_interactions(self) -> bool:
        """Check if there are more interactions to replay."""
        if not self.session_data:
            return False
        return self.interaction_index < len(self.session_data.get("interactions", []))

    def peek_pending(self) -> Optional[Dict[str, Any]]:
        """
        Return the interaction replay is positioned at, without consuming it.

        get_next_response advances past every non-matching interaction while
        it searches, so the position is gone by the time it reports a miss.
        Callers that want to say what the recording expected must peek first.
        """
        if not self.session_data:
            return None
        interactions = self.session_data.get("interactions", [])
        if self.interaction_index < len(interactions):
            return interactions[self.interaction_index]
        return None


class ReplayDivergenceError(RuntimeError):
    """
    Raised when a strict replay reaches a prompt the recording cannot answer.

    Batch runs have nobody at the keyboard, so falling through to live input
    blocks on stdin until the run is killed -- and because the failed search
    consumes the rest of the recording, every later prompt in that run blocks
    too. Strict replay turns that silent stall into a failure the batch loop
    can report against the run that caused it.
    """

    def __init__(self, interaction_type: str, prompt: str,
                 pending: Optional[Dict[str, Any]] = None):
        self.interaction_type = interaction_type
        self.prompt = prompt
        self.pending = pending

        if pending is None:
            expected = "the recording has no interactions left"
        else:
            expected = (
                f"the recording expects a {pending.get('type', 'prompt')} here: "
                f"{pending.get('prompt', '')!r}"
            )
        super().__init__(
            f"Replay diverged: no recorded answer for {interaction_type} "
            f"{prompt!r}. At this point {expected}."
        )


class InterceptedPrompt:
    """
    Intercepts Rich Prompt calls and built-in input() to record/replay user interactions.
    
    This class wraps both Rich Prompt functionality and Python's built-in input()
    to integrate with the session recording and replay system.
    """
    
    def __init__(self, recorder: Optional[SessionRecorder] = None,
                 replayer: Optional[SessionReplayer] = None,
                 strict_replay: bool = False):
        """
        Initialize the intercepted prompt.

        Args:
            recorder: Session recorder instance
            replayer: Session replayer instance
            strict_replay: If True (batch), a prompt the recording cannot
                answer raises ReplayDivergenceError instead of falling
                through to live input, which no one is there to type.
        """
        self.recorder = recorder
        self.replayer = replayer
        self.strict_replay = strict_replay
        self._original_prompt_ask = Prompt.ask
        self._original_confirm_ask = Confirm.ask
        # IntPrompt/FloatPrompt are NOT subclasses of Prompt -- they inherit
        # PromptBase.ask -- so patching Prompt.ask never covered them. They fell
        # through to the builtin input() interception, where Rich has already
        # printed the question itself and passes nothing on, so every numeric
        # answer was recorded as type='input' with an EMPTY prompt string and
        # replayed positionally against any other numeric answer.
        self._original_int_ask = IntPrompt.ask
        self._original_float_ask = FloatPrompt.ask
        self._original_input = input
        self._in_rich_prompt = False  # Guard: True while inside Prompt/Confirm.ask

    def install(self):
        """Install the prompt interceptors."""
        # Monkey-patch the Rich prompt methods
        Prompt.ask = self._intercepted_prompt_ask
        Confirm.ask = self._intercepted_confirm_ask
        IntPrompt.ask = self._intercepted_int_ask
        FloatPrompt.ask = self._intercepted_float_ask
        
        # Monkey-patch the built-in input function
        import builtins
        builtins.input = self._intercepted_input
        
    def uninstall(self):
        """Restore the original prompt methods."""
        Prompt.ask = self._original_prompt_ask
        Confirm.ask = self._original_confirm_ask
        IntPrompt.ask = self._original_int_ask
        FloatPrompt.ask = self._original_float_ask
        
        # Restore the built-in input function
        import builtins
        builtins.input = self._original_input
        
    def _replayed_response(self, interaction_type: str,
                           prompt: str) -> Optional[str]:
        """
        Get the recorded answer for a prompt, or None to ask the user.

        None covers both "this prompt diverged from the recording" and
        "the recording ran out" -- indistinguishable here. Interactively that
        rightly means hand control back to the user; under strict replay it
        means the template cannot drive this run, so raise instead.
        """
        if not (self.replayer and self.replayer.replaying):
            return None

        pending = self.replayer.peek_pending()
        response = self.replayer.get_next_response(interaction_type, prompt)
        if response is None and self.strict_replay:
            raise ReplayDivergenceError(interaction_type, prompt, pending)
        return response

    def _intercepted_prompt_ask(self, prompt: str, **kwargs) -> str:
        """Intercepted version of Prompt.ask."""
        choices = kwargs.get('choices', None)
        default = kwargs.get('default', None)

        # Check if we're replaying
        response = self._replayed_response("prompt", prompt)
        if response is not None:
            print(f"{prompt} [{response}]")
            return response

        # Normal prompt — guard prevents _intercepted_input from
        # double-recording the same interaction
        self._in_rich_prompt = True
        try:
            response = self._original_prompt_ask(prompt, **kwargs)
        finally:
            self._in_rich_prompt = False

        # Record the interaction
        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "prompt",
                prompt,
                response,
                choices=choices,
                context={"default": default}
            )

        return response

    def _intercepted_numeric_ask(self, kind, original, cast, prompt: str, **kwargs):
        """Shared path for IntPrompt.ask / FloatPrompt.ask.

        Recorded as a ``"prompt"`` interaction carrying the real question text
        so it is matched by TEXT. Reaching the builtin-input interception
        instead recorded an empty prompt string (Rich prints the question
        itself), which made every numeric answer in a log interchangeable.
        """
        default = kwargs.get('default', None)

        response = self._replayed_response("prompt", prompt)
        if response is not None:
            try:
                value = cast(response)
                print(f"{prompt} [{response}]")
                return value
            except (TypeError, ValueError):
                # A recorded answer that is not a number belongs to another
                # question; fall through and ask rather than guess.
                print(f"{prompt} [{response!r} is not a {kind} — asking]")

        self._in_rich_prompt = True
        try:
            value = original(prompt, **kwargs)
        finally:
            self._in_rich_prompt = False

        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "prompt", prompt, str(value), context={"default": default},
            )
        return value

    def _intercepted_int_ask(self, prompt: str, **kwargs) -> int:
        """Intercepted version of IntPrompt.ask."""
        return self._intercepted_numeric_ask(
            "integer", self._original_int_ask, int, prompt, **kwargs)

    def _intercepted_float_ask(self, prompt: str, **kwargs) -> float:
        """Intercepted version of FloatPrompt.ask."""
        return self._intercepted_numeric_ask(
            "number", self._original_float_ask, float, prompt, **kwargs)

    def _intercepted_confirm_ask(self, prompt: str, **kwargs) -> bool:
        """Intercepted version of Confirm.ask."""
        default = kwargs.get('default', None)

        # Check if we're replaying
        response = self._replayed_response("confirm", prompt)
        if response is not None:
            bool_response = response.lower() in ('yes', 'y', 'true', '1')
            print(f"{prompt} [{'Yes' if bool_response else 'No'}]")
            return bool_response

        # Normal confirm — guard prevents double-recording
        self._in_rich_prompt = True
        try:
            response = self._original_confirm_ask(prompt, **kwargs)
        finally:
            self._in_rich_prompt = False

        # Record the interaction
        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "confirm",
                prompt,
                "yes" if response else "no",
                context={"default": default}
            )

        return response

    def _intercepted_input(self, prompt: str = "") -> str:
        """Intercepted version of built-in input()."""
        # Skip if called from inside a Rich prompt — that interaction
        # is already recorded by _intercepted_prompt_ask/_intercepted_confirm_ask
        if self._in_rich_prompt:
            return self._original_input(prompt)

        # Check if we're replaying
        response = self._replayed_response("input", prompt)
        if response is not None:
            print(f"{prompt}{response}")
            return response

        # Normal input
        response = self._original_input(prompt)

        # Record the interaction
        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "input",
                prompt,
                response,
                context={"builtin": True}
            )

        return response


class HybridInterceptor:
    """
    Hybrid interceptor that replays existing interactions, then switches to recording.

    This enables the workflow: replay existing session -> continue with new interactions.
    Once replay is exhausted, it automatically switches to recording mode.
    """

    def __init__(self, recorder: Optional[SessionRecorder] = None,
                 replayer: Optional[SessionReplayer] = None):
        """
        Initialize the hybrid interceptor.

        Args:
            recorder: Session recorder instance (for new interactions)
            replayer: Session replayer instance (for existing interactions)
        """
        self.recorder = recorder
        self.replayer = replayer
        self._original_prompt_ask = Prompt.ask
        self._original_confirm_ask = Confirm.ask
        # IntPrompt/FloatPrompt are NOT subclasses of Prompt -- they inherit
        # PromptBase.ask -- so patching Prompt.ask never covered them. They fell
        # through to the builtin input() interception, where Rich has already
        # printed the question itself and passes nothing on, so every numeric
        # answer was recorded as type='input' with an EMPTY prompt string and
        # replayed positionally against any other numeric answer.
        self._original_int_ask = IntPrompt.ask
        self._original_float_ask = FloatPrompt.ask
        self._original_input = input
        self._replay_exhausted = False
        self._in_rich_prompt = False  # Guard: True while inside Prompt/Confirm.ask

    def install(self):
        """Install the prompt interceptors."""
        # Monkey-patch the Rich prompt methods
        Prompt.ask = self._intercepted_prompt_ask
        Confirm.ask = self._intercepted_confirm_ask
        IntPrompt.ask = self._intercepted_int_ask
        FloatPrompt.ask = self._intercepted_float_ask

        # Monkey-patch the built-in input function
        import builtins
        builtins.input = self._intercepted_input

    def uninstall(self):
        """Restore the original prompt methods."""
        Prompt.ask = self._original_prompt_ask
        Confirm.ask = self._original_confirm_ask
        IntPrompt.ask = self._original_int_ask
        FloatPrompt.ask = self._original_float_ask

        # Restore the built-in input function
        import builtins
        builtins.input = self._original_input

    def _intercepted_prompt_ask(self, prompt: str, **kwargs) -> str:
        """Intercepted version of Prompt.ask with hybrid replay-record."""
        choices = kwargs.get('choices', None)
        default = kwargs.get('default', None)

        # Clear so callers can detect whether THIS prompt was answered by replay
        # vs live input (matters after _replay_exhausted, when get_next_response
        # isn't called and the flag would otherwise carry over from the last match).
        if self.replayer is not None:
            self.replayer.last_returned_interaction = None

        # Try replay first (if available and not exhausted)
        if self.replayer and self.replayer.replaying and not self._replay_exhausted:
            response = self.replayer.get_next_response("prompt", prompt)
            if response is not None:
                print(f"{prompt} [REPLAY: {response}]")
                return response
            else:
                # No match found for this prompt
                # Check if we've reached the end of all interactions
                if not self.replayer.has_more_interactions():
                    # Truly exhausted - no more interactions to try
                    self._replay_exhausted = True
                    print("\n[Replay complete - now recording new interactions]")
                # else: No match for THIS prompt, but continue trying for future prompts

        # Normal prompt — guard prevents _intercepted_input from
        # double-recording the same interaction
        self._in_rich_prompt = True
        try:
            response = self._original_prompt_ask(prompt, **kwargs)
        finally:
            self._in_rich_prompt = False

        # Record the interaction
        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "prompt",
                prompt,
                response,
                choices=choices,
                context={"default": default}
            )

        return response

    def _intercepted_numeric_ask(self, kind, original, cast, prompt: str, **kwargs):
        """Shared replay-record path for IntPrompt.ask / FloatPrompt.ask.

        Recorded as an ordinary ``"prompt"`` interaction carrying the real
        question text, so it is matched by TEXT like every other prompt. These
        used to reach the builtin-input interception instead, which sees an
        empty prompt string (Rich prints the question itself), so every numeric
        answer in a log was interchangeable with every other one — a newly
        added numeric question silently consumed the answer meant for a
        different one.
        """
        default = kwargs.get('default', None)

        if self.replayer is not None:
            self.replayer.last_returned_interaction = None

        if self.replayer and self.replayer.replaying and not self._replay_exhausted:
            response = self.replayer.get_next_response("prompt", prompt)
            if response is not None:
                try:
                    value = cast(response)
                    print(f"{prompt} [REPLAY: {response}]")
                    return value
                except (TypeError, ValueError):
                    # A recorded answer that is not a number belongs to some
                    # other question; ask rather than guess.
                    print(f"{prompt} [REPLAY: {response!r} is not a {kind} — asking]")
            elif not self.replayer.has_more_interactions():
                self._replay_exhausted = True
                print("\n[Replay complete - now recording new interactions]")

        self._in_rich_prompt = True
        try:
            value = original(prompt, **kwargs)
        finally:
            self._in_rich_prompt = False

        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "prompt", prompt, str(value), context={"default": default},
            )
        return value

    def _intercepted_int_ask(self, prompt: str, **kwargs) -> int:
        """Intercepted version of IntPrompt.ask with hybrid replay-record."""
        return self._intercepted_numeric_ask(
            "integer", self._original_int_ask, int, prompt, **kwargs)

    def _intercepted_float_ask(self, prompt: str, **kwargs) -> float:
        """Intercepted version of FloatPrompt.ask with hybrid replay-record."""
        return self._intercepted_numeric_ask(
            "number", self._original_float_ask, float, prompt, **kwargs)

    def _intercepted_confirm_ask(self, prompt: str, **kwargs) -> bool:
        """Intercepted version of Confirm.ask with hybrid replay-record."""
        default = kwargs.get('default', None)

        if self.replayer is not None:
            self.replayer.last_returned_interaction = None

        # Try replay first (if available and not exhausted)
        if self.replayer and self.replayer.replaying and not self._replay_exhausted:
            response = self.replayer.get_next_response("confirm", prompt)
            if response is not None:
                bool_response = response.lower() in ('yes', 'y', 'true', '1')
                print(f"{prompt} [REPLAY: {'Yes' if bool_response else 'No'}]")
                return bool_response
            else:
                # No match found for this prompt
                # Check if we've reached the end of all interactions
                if not self.replayer.has_more_interactions():
                    # Truly exhausted - no more interactions to try
                    self._replay_exhausted = True
                    print("\n[Replay complete - now recording new interactions]")
                # else: No match for THIS prompt, but continue trying for future prompts

        # Normal confirm — guard prevents double-recording
        self._in_rich_prompt = True
        try:
            response = self._original_confirm_ask(prompt, **kwargs)
        finally:
            self._in_rich_prompt = False

        # Record the interaction
        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "confirm",
                prompt,
                "yes" if response else "no",
                context={"default": default}
            )

        return response

    def _intercepted_input(self, prompt: str = "") -> str:
        """Intercepted version of built-in input() with hybrid replay-record."""
        # Skip if called from inside a Rich prompt — that interaction
        # is already handled by _intercepted_prompt_ask/_intercepted_confirm_ask
        if self._in_rich_prompt:
            return self._original_input(prompt)

        if self.replayer is not None:
            self.replayer.last_returned_interaction = None

        # Try replay first (if available and not exhausted)
        if self.replayer and self.replayer.replaying and not self._replay_exhausted:
            response = self.replayer.get_next_response("input", prompt)
            if response is not None:
                print(f"{prompt}[REPLAY: {response}]")
                return response
            else:
                # No match found for this prompt
                # Check if we've reached the end of all interactions
                if not self.replayer.has_more_interactions():
                    # Truly exhausted - no more interactions to try
                    self._replay_exhausted = True
                    print("\n[Replay complete - now recording new interactions]")
                # else: No match for THIS prompt, but continue trying for future prompts

        # Normal input (record if recorder is active)
        response = self._original_input(prompt)

        # Record the interaction
        if self.recorder and self.recorder.recording:
            self.recorder.record_interaction(
                "input",
                prompt,
                response,
                context={"builtin": True}
            )

        return response


class SessionManager:
    """
    High-level manager for session recording and replay in MPSA.
    
    This class provides a simple interface for integrating session
    recording/replay into the MPSA application.
    """
    
    def __init__(self):
        """Initialize the session manager."""
        self.recorder: Optional[SessionRecorder] = None
        self.replayer: Optional[SessionReplayer] = None
        self.interceptor: Optional[InterceptedPrompt] = None
        self.active = False
        
    def start_recording(self, filename: str, metadata: Dict[str, Any] = None):
        """
        Start recording a new session.
        
        Args:
            filename: Path to save the session log
            metadata: Additional metadata about the session
        """
        if self.active:
            raise RuntimeError("Session manager is already active")
            
        self.recorder = SessionRecorder(filename)
        self.recorder.start_recording(metadata)
        
        self.interceptor = InterceptedPrompt(recorder=self.recorder)
        self.interceptor.install()
        
        self.active = True
        # Session recording start message handled by main.py
        
    def start_replay(self, filename: str, replay_delay: float = 0.0):
        """
        Start replaying a recorded session.

        Args:
            filename: Path to the session log to replay
            replay_delay: Delay in seconds between interactions (for demos)
        """
        if self.active:
            raise RuntimeError("Session manager is already active")

        if not os.path.exists(filename):
            raise FileNotFoundError(f"Session file not found: {filename}")

        self.replayer = SessionReplayer(filename, replay_delay=replay_delay,
                                        strict_variables=False)
        self.replayer.start_replay()

        self.interceptor = InterceptedPrompt(replayer=self.replayer)
        self.interceptor.install()

        self.active = True
        print(f"[Session replay started: {filename}]")

    def start_hybrid_mode(self, filename: str, truncate_at: Optional[int] = None,
                          replay_delay: float = 0.0, keep_following_interactions: bool = False,
                          new_value: Optional[str] = None):
        """
        Start hybrid mode: replay existing session, then record new interactions.

        This enables resuming a session and continuing with new steps.
        All new interactions will be appended to the session file.

        Args:
            filename: Path to the session log to replay and continue
            truncate_at: Optional index to truncate/edit at
            replay_delay: Delay in seconds between replay interactions (for demos)
            keep_following_interactions: If True, when editing at truncate_at, keep
                                        interactions after that point and try to replay them.
                                        If False (default), discard all following interactions.
            new_value: If provided, replace the response at truncate_at with this value
        """
        if self.active:
            raise RuntimeError("Session manager is already active")

        if not os.path.exists(filename):
            raise FileNotFoundError(f"Session file not found: {filename}")

        # Back up the original file before any in-place rewrite.
        # The recorder opens this same path in 'w' mode below, which clobbers the
        # on-disk copy. If replay diverges or the user aborts mid-session, the
        # backup is the only way to recover the pre-edit log.
        backup_path = f"{filename}.bak-{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        try:
            import shutil
            shutil.copy2(filename, backup_path)
            print(f"[Session backup: {os.path.basename(backup_path)}]")
        except OSError as e:
            print(f"[yellow]Warning: could not back up session log before resume: {e}[/yellow]")
            print("[yellow]Aborting hybrid mode to avoid clobbering the original log.[/yellow]")
            raise

        # Set up replayer. Interactive resume is lenient about template
        # variables: unresolved placeholders fall through to live input rather
        # than crashing (batch replay stays strict via its own construction).
        self.replayer = SessionReplayer(filename, replay_delay=replay_delay,
                                        strict_variables=False)

        # Handle edit mode
        if truncate_at is not None and new_value is not None:
            # Edit mode: Update the response value at truncate_at
            if (self.replayer.session_data and
                'interactions' in self.replayer.session_data and
                truncate_at < len(self.replayer.session_data['interactions'])):

                # Update the response value
                self.replayer.session_data['interactions'][truncate_at]['response'] = new_value

                if keep_following_interactions:
                    print(f"[Session edit at interaction {truncate_at} - updated value and will replay following interactions]")
                else:
                    print(f"[Session edit at interaction {truncate_at} - updated value]")

        # Handle truncation if requested
        if truncate_at is not None and not keep_following_interactions:
            # Discard all interactions after truncate_at
            self.replayer.truncate_at(truncate_at + 1)  # +1 to keep the edited interaction
            print(f"[Session truncated after interaction {truncate_at} - following interactions discarded]")

        self.replayer.start_replay()

        # Set up recorder to continue from existing session
        self.recorder = SessionRecorder(filename)
        self.recorder.load_existing_session(filename)

        # Update recorder's session data with the new value if provided
        if truncate_at is not None and new_value is not None:
            if (truncate_at < len(self.recorder.session_data.get('interactions', []))):
                # Update the response in the recorder's session data too
                self.recorder.session_data['interactions'][truncate_at]['response'] = new_value

        # Truncate recorder's session data to match behavior
        if truncate_at is not None and not keep_following_interactions:
            # Truncate after the edited interaction (keep the edit, discard following)
            self.recorder.truncate_at(truncate_at + 1)

        self.recorder.start_recording()

        # Install hybrid interceptor
        self.interceptor = HybridInterceptor(
            recorder=self.recorder,
            replayer=self.replayer
        )
        self.interceptor.install()

        self.active = True

        # Show status
        from rich.console import Console
        console = Console()
        console.print(f"\n[bold green]Hybrid mode started:[/bold green]")
        console.print(f"  Session file: {filename}")
        if truncate_at is not None:
            mode_desc = "smart replay" if keep_following_interactions else "truncated"
            console.print(f"  Edit at: interaction {truncate_at} ({mode_desc})")
        console.print(f"  Mode: Replay existing → Record new")
        console.print()
        
    def stop(self):
        """Stop recording or replaying."""
        if not self.active:
            return
            
        if self.interceptor:
            self.interceptor.uninstall()
            
        if self.recorder:
            self.recorder.stop_recording()
            print("[Session recording stopped]")
            
        if self.replayer:
            self.replayer.stop_replay()
            print("[Session replay stopped]")
            
        self.recorder = None
        self.replayer = None
        self.interceptor = None
        self.active = False
        
    def is_recording(self) -> bool:
        """Check if currently recording."""
        return self.recorder is not None and self.recorder.recording

    def is_replaying(self) -> bool:
        """Check if currently replaying."""
        return self.replayer is not None and self.replayer.replaying

    def get_last_replayed_interaction(self) -> Optional[Dict[str, Any]]:
        """
        Return the most recently replayed interaction dict, or None.

        None means the last prompt was answered by live input (record path)
        rather than by a recorded response. Used by callers to decide whether
        to re-resolve a recorded selection or to annotate the just-recorded
        interaction with resolver hints.
        """
        if self.replayer is None:
            return None
        return self.replayer.last_returned_interaction

    def annotate_last_recorded(self, extra_context: Dict[str, Any]):
        """
        Attach extra context to the most recently recorded interaction.

        No-op if no recorder is active or no interactions have been recorded.
        """
        if self.recorder is None or not self.recorder.recording:
            return
        self.recorder.annotate_last_interaction(extra_context)


# Utility functions for session display

def format_interaction_display(interaction: Dict[str, Any]) -> str:
    """
    Format an interaction for display with smart context handling.

    This function intelligently displays interactions based on available context.
    It gracefully handles both simple prompts and rich context.

    Args:
        interaction: Interaction dict from session data

    Returns:
        Formatted string for display
    """
    context = interaction.get("context", {})
    prompt = interaction.get("prompt", "")
    response = interaction.get("response", "")
    interaction_type = interaction.get("type", "")

    # Try to get rich context first
    description = context.get("description")
    module_name = context.get("module")
    option_label = context.get("option_label")
    options_map = context.get("options_map")

    # Build display string based on available context
    parts = []

    # Add module name if available
    if module_name:
        parts.append(f"[{module_name}]")

    # Add description or prompt
    if description:
        parts.append(description)
    elif prompt:
        # Trim long prompts
        display_prompt = prompt if len(prompt) <= 60 else prompt[:57] + "..."
        parts.append(display_prompt)

    # Add response with context
    if option_label:
        # Rich context: show the label
        response_str = f"→ {option_label}"
    elif options_map and response in options_map:
        # Lookup response in options map
        response_str = f"→ {options_map[response]}"
    else:
        # Simple response
        if interaction_type == "confirm":
            response_str = f"→ {'Yes' if response.lower() in ('yes', 'y', 'true', '1') else 'No'}"
        else:
            response_str = f"→ {response}"

    parts.append(response_str)

    return " ".join(parts)


def get_interaction_summary(interaction: Dict[str, Any]) -> str:
    """
    Get a one-line summary of an interaction for list display.

    Args:
        interaction: Interaction dict from session data

    Returns:
        One-line summary string
    """
    timestamp = interaction.get("timestamp", "")
    interaction_type = interaction.get("type", "").upper()

    # Extract time from ISO timestamp
    time_str = "??:??"
    if timestamp:
        try:
            dt = datetime.fromisoformat(timestamp)
            time_str = dt.strftime("%H:%M:%S")
        except:
            pass

    # Get formatted display
    display = format_interaction_display(interaction)

    return f"{time_str}  {interaction_type:8} │ {display}"


# Integration functions for MPSA

def integrate_session_manager(processor):
    """
    Integrate the session manager into the PDBProcessor.

    This function adds session recording/replay capabilities to
    an existing PDBProcessor instance.

    Args:
        processor: PDBProcessor instance
    """
    # Add session manager as an attribute
    processor.session_manager = SessionManager()

    # Add methods to the processor
    processor.start_session_recording = lambda filename, metadata=None: \
        processor.session_manager.start_recording(filename, metadata)

    processor.start_session_replay = lambda filename, delay=0.0: \
        processor.session_manager.start_replay(filename, delay)

    processor.start_session_hybrid = lambda filename, truncate_at=None, delay=0.0, keep_following=False, new_value=None: \
        processor.session_manager.start_hybrid_mode(filename, truncate_at, delay, keep_following, new_value)

    processor.stop_session = lambda: processor.session_manager.stop()

    # Modify the cleanup method to stop session
    original_cleanup = processor.cleanup

    def enhanced_cleanup():
        processor.session_manager.stop()
        original_cleanup()

    processor.cleanup = enhanced_cleanup