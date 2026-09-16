"""
Rewind: undo for a prompt-driven session.

Typing ``undo`` at any prompt rewinds the session to an earlier answer. The
user picks the answer to revisit; ProPrep then unwinds to ``main()``, builds a
fresh processor, and replays the session log up to that answer in hybrid mode
before handing control back live. Nothing is reversed in place: ProPrep's state
includes the position in the Python call stack, and only re-executing the
program from its recorded inputs can put the interpreter back there.

Three pieces live here:

``RewindPrompt`` and friends
    Thin subclasses of Rich's prompt classes whose ``process_response`` accepts
    the keyword before Rich's own choice validation runs. Rich re-asks until
    ``process_response`` returns, so a prompt with ``choices`` would otherwise
    reject the word as an invalid choice and the interceptor would never see it.

``request_rewind``
    The picker. It lists the most recent recorded answers, asks which to
    rewind to and how, and raises :class:`RewindRequested`. Its own prompts go
    through the unpatched Rich functions and are never recorded.

``RewindRequested``
    The signal. It derives from ``BaseException`` so the many
    ``except Exception`` handlers between a prompt and ``main()`` let it pass,
    the same way ``SystemExit`` already unwinds a quit from any depth.
"""

from typing import Any, Dict, List, Optional

from rich.prompt import Confirm, FloatPrompt, IntPrompt, Prompt

REWIND_KEYWORD = "undo"
PICKER_PAGE = 12


def is_rewind_keyword(text: Any) -> bool:
    return isinstance(text, str) and text.strip().lower() == REWIND_KEYWORD


class RewindKeyword(Exception):
    """Raised inside Rich's validation loop when the keyword is typed.

    Internal: it travels only from ``process_response`` to the interceptor
    that invoked the prompt, which turns it into a picker run.
    """


class RewindRequested(BaseException):
    """Unwind to ``main()`` and replay the session to ``index``.

    Attributes:
        index: Position in the session log of the answer to revisit.
        mode: ``"reanswer"`` asks the question at ``index`` again live and
            discards what followed; ``"change"`` replaces the answer and replays
            what followed.
        new_value: Replacement answer for ``"change"`` mode.
        session_file: Path of the log being rewound.
    """

    def __init__(self, index: int, mode: str, session_file: str,
                 new_value: Optional[str] = None):
        self.index = index
        self.mode = mode
        self.new_value = new_value
        self.session_file = session_file
        super().__init__(f"rewind to interaction {index} ({mode})")

    def hybrid_args(self) -> Dict[str, Any]:
        """Arguments for ``SessionManager.start_hybrid_mode``.

        Hybrid mode keeps the interaction at ``truncate_at`` and drops what
        follows, so re-answering ``index`` means truncating after ``index - 1``.
        Index 0 therefore gives -1, which hybrid mode treats as "replay
        nothing": the first question is asked live.
        """
        if self.mode == "change":
            return {"truncate_at": self.index, "keep_following": True,
                    "new_value": self.new_value}
        return {"truncate_at": self.index - 1, "keep_following": False,
                "new_value": None}


# ---------------------------------------------------------------------------
# Keyword-aware Rich prompt classes
# ---------------------------------------------------------------------------

class RewindPrompt(Prompt):
    def process_response(self, value: str) -> str:
        if is_rewind_keyword(value):
            raise RewindKeyword()
        return super().process_response(value)


class RewindConfirm(Confirm):
    def process_response(self, value: str) -> bool:
        if is_rewind_keyword(value):
            raise RewindKeyword()
        return super().process_response(value)


class RewindIntPrompt(IntPrompt):
    def process_response(self, value: str) -> int:
        if is_rewind_keyword(value):
            raise RewindKeyword()
        return super().process_response(value)


class RewindFloatPrompt(FloatPrompt):
    def process_response(self, value: str) -> float:
        if is_rewind_keyword(value):
            raise RewindKeyword()
        return super().process_response(value)


def _call_with_class(original, rewind_cls, prompt: str, **kwargs):
    """Invoke a captured ``ask`` classmethod with ``rewind_cls`` as ``cls``.

    ``original`` is the bound classmethod the interceptor saved before
    patching (``Prompt.ask`` bound to ``Prompt``). Calling its underlying
    function with the keyword-aware subclass runs Rich's real prompt loop
    without re-entering the patched attribute. Anything else, such as a test
    stub, is called as it is.
    """
    func = getattr(original, "__func__", None)
    owner = getattr(original, "__self__", None)
    if func is not None and isinstance(owner, type) and issubclass(rewind_cls, owner):
        return func(rewind_cls, prompt, **kwargs)
    return original(prompt, **kwargs)


def ask_original(interceptor, original, rewind_cls, prompt: str, **kwargs):
    """Ask through the original Rich prompt, honouring the rewind keyword.

    Returns the answer. If the user typed the keyword, the picker runs; it
    either raises :class:`RewindRequested` or, when a rewind is not possible
    or was cancelled, returns so the question is asked again.
    """
    while True:
        try:
            return _call_with_class(original, rewind_cls, prompt, **kwargs)
        except RewindKeyword:
            request_rewind(interceptor)


def input_original(interceptor, original_input, prompt: str = "") -> str:
    """Builtin ``input`` counterpart of :func:`ask_original`."""
    while True:
        response = original_input(prompt)
        if not is_rewind_keyword(response):
            return response
        request_rewind(interceptor)


# ---------------------------------------------------------------------------
# Picker
# ---------------------------------------------------------------------------

def _question_text(interaction: Dict[str, Any]) -> str:
    context = interaction.get("context") or {}
    text = context.get("description") or interaction.get("prompt") or ""
    text = " ".join(str(text).split())
    return text if len(text) <= 44 else text[:41] + "..."


def _answer_text(interaction: Dict[str, Any]) -> str:
    context = interaction.get("context") or {}
    response = interaction.get("response", "")
    label = context.get("option_label")
    if not label:
        options_map = context.get("options_map") or {}
        label = options_map.get(response)
    if label:
        return f"{response}  {label}"
    if interaction.get("type") == "confirm":
        return "yes" if str(response).lower() in ("yes", "y", "true", "1") else "no"
    return str(response)


def _print_page(console, interactions: List[Dict[str, Any]], start: int, end: int):
    from rich.table import Table
    table = Table(title="Recent answers in this session", title_justify="left",
                  show_edge=False, pad_edge=False, expand=False)
    table.add_column("#", justify="right", style="bold blue")
    table.add_column("module", style="grey50")
    table.add_column("question")
    table.add_column("answer", style="green")
    for i in range(start, end):
        it = interactions[i]
        module = (it.get("context") or {}).get("module") or it.get("type", "")
        table.add_row(str(i), str(module), _question_text(it), _answer_text(it))
    console.print()
    console.print(table)


def _validate_new_value(interaction: Dict[str, Any], value: str) -> Optional[str]:
    """Return the canonical value, or None with a reason printed by the caller."""
    if interaction.get("type") == "confirm":
        low = value.strip().lower()
        if low in ("y", "yes", "true", "1"):
            return "yes"
        if low in ("n", "no", "false", "0"):
            return "no"
        return None
    choices = interaction.get("choices")
    if choices:
        lowered = [str(c).lower() for c in choices]
        if value.strip().lower() not in lowered:
            return None
        return str(choices[lowered.index(value.strip().lower())])
    return value


def request_rewind(interceptor) -> None:
    """Run the picker and raise :class:`RewindRequested`.

    Returns without raising when no rewind is possible (no active recorder,
    nothing recorded yet) or when the user cancels; the caller then re-asks
    its question.
    """
    from rich.console import Console
    console = Console()

    recorder = getattr(interceptor, "recorder", None)
    if recorder is None or not getattr(recorder, "recording", False) \
            or not getattr(recorder, "record_file", None):
        console.print(
            "[yellow]Rewind needs an active session recording, and this run has "
            "none (replay-only or --no-session). Answer the question to continue.[/yellow]")
        return

    interactions = recorder.session_data.get("interactions", [])
    if not interactions:
        console.print("[yellow]Nothing to rewind yet: no answers have been recorded.[/yellow]")
        return

    # The picker's own prompts must not be recorded or replayed. Rich's
    # console.input calls the patched builtin, which records unless the
    # interceptor believes it is already inside a Rich prompt.
    guard_before = getattr(interceptor, "_in_rich_prompt", False)
    interceptor._in_rich_prompt = True
    try:
        ask = getattr(interceptor, "_original_prompt_ask")
        total = len(interactions)
        end = total
        while True:
            start = max(0, end - PICKER_PAGE)
            _print_page(console, interactions, start, end)
            hint = "[grey50]number to rewind to, \\[m]ore, or \\[x] cancel[/grey50]"
            raw = _call_with_class(
                ask, Prompt, f"Rewind to which answer? {hint}",
                default=str(total - 1), show_default=True)
            text = str(raw).strip().lower()
            if text == "x":
                console.print("[grey50]Rewind cancelled.[/grey50]")
                return
            if text == "m":
                end = start if start > 0 else total
                continue
            try:
                index = int(text)
            except ValueError:
                console.print(f"[red]Enter a number between 0 and {total - 1}.[/red]")
                continue
            if not 0 <= index < total:
                console.print(f"[red]Enter a number between 0 and {total - 1}.[/red]")
                continue
            break

        target = interactions[index]
        console.print()
        console.print(f"  [bold]{index}[/bold]  {_question_text(target)}  "
                      f"[green]{_answer_text(target)}[/green]")
        mode_raw = _call_with_class(
            ask, Prompt,
            "\\[r]e-answer it live, or \\[c]hange it and replay what followed? "
            "[grey50](\\[x] cancel)[/grey50]",
            choices=["r", "c", "x"], default="r", show_choices=False,
            case_sensitive=False)
        mode = str(mode_raw).strip().lower()
        if mode == "x":
            console.print("[grey50]Rewind cancelled.[/grey50]")
            return

        new_value = None
        if mode == "c":
            question = _question_text(target)
            choices = target.get("choices")
            while True:
                raw_value = _call_with_class(ask, Prompt, f'New answer for "{question}"')
                value = _validate_new_value(target, str(raw_value))
                if value is not None:
                    new_value = value
                    break
                if target.get("type") == "confirm":
                    console.print("[red]Answer yes or no.[/red]")
                else:
                    console.print(f"[red]Answer must be one of: {', '.join(map(str, choices))}[/red]")
    finally:
        interceptor._in_rich_prompt = guard_before

    if mode == "c":
        after = total - index - 1
        console.print(
            f"[grey50][Rewinding: answer {index} becomes {new_value!r}; "
            f"{after} later answer{'s' if after != 1 else ''} will be replayed][/grey50]")
    else:
        dropped = total - index
        console.print(
            f"[grey50][Rewinding: question {index} will be asked again; "
            f"{dropped} answer{'s' if dropped != 1 else ''} from there on are set aside "
            f"in a backup of the session log][/grey50]")
    raise RewindRequested(index, "change" if mode == "c" else "reanswer",
                          recorder.record_file, new_value)
