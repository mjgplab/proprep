"""
AMBER MDIN Configuration Interface - Layer 1: MDIN Parser

Parses AMBER MDIN (sander/pmemd input) files using Fortran namelist syntax.

Handles:
- Title lines (before first &)
- Namelist blocks: &name ... /
- Key=value pairs with commas/newlines
- String values (single-quoted, may contain commas)
- Array syntax: mask(1)='...' -> mask_1='...'
- Logical values: .true./.false., T/F, 1/0
- Comments: ! to end of line (outside quotes)
- &wt blocks with TYPE tracking
- File redirect commands after &wt END

Usage:
    from proprep.md_prep.amber import parse_mdin, parse_mdin_file

    # From text
    parsed = parse_mdin(mdin_text)
    nstlim = parsed.get_value("nstlim")

    # From file
    parsed = parse_mdin_file("production.mdin")
"""

import logging
import re
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Any

from .exceptions import MdinSyntaxError

logger = logging.getLogger(__name__)


# =============================================================================
# Data structures
# =============================================================================


@dataclass
class MdinValue:
    """A parsed key=value entry with source location."""

    key: str  # lowercase normalized name
    raw_key: str  # original case as written
    value: Any  # Python-typed: int, float, str, bool, list
    raw_value: str  # original text as written
    line_number: int  # for error reporting
    section: str  # which namelist block this was in


@dataclass
class MdinBlock:
    """A parsed namelist block."""

    section: str  # e.g. "cntrl", "ewald" (lowercase)
    raw_section: str  # original case
    entries: dict[str, MdinValue] = field(default_factory=dict)  # key -> MdinValue
    line_start: int = 0
    line_end: int = 0


@dataclass
class WtBlock:
    """A parsed &wt block with its TYPE value."""

    type_value: str | None  # The TYPE='...' value, None if not found
    entries: dict[str, MdinValue] = field(default_factory=dict)
    line_start: int = 0
    line_end: int = 0


@dataclass
class FileRedirectEntry:
    """A file redirection command (e.g., DISANG=restraints.RST)."""

    redirect_type: str  # e.g., "DISANG", "DUMPAVE"
    filename: str
    line_number: int


@dataclass
class ParsedMdin:
    """Complete parsed MDIN file."""

    title: str
    blocks: dict[str, MdinBlock] = field(default_factory=dict)  # section name -> block
    wt_blocks: list[WtBlock] = field(default_factory=list)  # &wt blocks in order
    file_redirects: list[FileRedirectEntry] = field(default_factory=list)
    source_path: str | None = None
    raw_text: str = ""

    def get(self, key: str, section: str | None = None) -> MdinValue | None:
        """
        Look up a value. If section omitted, search all blocks.

        Args:
            key: Keyword name (case-insensitive)
            section: Section to search in, or None for all

        Returns:
            MdinValue if found, None otherwise
        """
        key_lower = key.lower()
        if section is not None:
            block = self.blocks.get(section.lower())
            if block:
                return block.entries.get(key_lower)
            return None

        # Search all blocks
        for block in self.blocks.values():
            if key_lower in block.entries:
                return block.entries[key_lower]
        return None

    def get_value(self, key: str, section: str | None = None) -> Any:
        """
        Get just the value (not MdinValue) for a keyword.

        Args:
            key: Keyword name (case-insensitive)
            section: Section to search in, or None for all

        Returns:
            The Python value if found, None otherwise
        """
        entry = self.get(key, section)
        return entry.value if entry else None

    def has(self, key: str, section: str | None = None) -> bool:
        """Check if a key is set (not defaulted)."""
        return self.get(key, section) is not None

    def all_keys(self) -> list[MdinValue]:
        """All set key-value pairs across all blocks."""
        result: list[MdinValue] = []
        for block in self.blocks.values():
            result.extend(block.entries.values())
        return result

    def sections(self) -> list[str]:
        """Names of all namelist blocks present."""
        return list(self.blocks.keys())


# =============================================================================
# Parser implementation
# =============================================================================


class _ParserState(Enum):
    """Parser state machine states."""

    TITLE = auto()  # Before any & block
    IN_BLOCK = auto()  # Inside &name ... / block
    AFTER_BLOCKS = auto()  # After final / (file redirects region)


class _MdinParser:
    """Internal parser implementation."""

    # Patterns for tokenization
    _NAMELIST_START = re.compile(r"&(\w+)", re.IGNORECASE)
    _NAMELIST_END = re.compile(r"^\s*/\s*$|(?<!['\w])/(?!['\w])")
    _KEY_VALUE = re.compile(
        r"(\w+(?:\(\d+\))?)\s*=\s*"  # key or key(n) followed by =
    )
    _STRING_VALUE = re.compile(r"'([^']*)'")
    _LOGICAL_TRUE = re.compile(r"^\.true\.$|^t$|^\.t\.$", re.IGNORECASE)
    _LOGICAL_FALSE = re.compile(r"^\.false\.$|^f$|^\.f\.$", re.IGNORECASE)
    _FILE_REDIRECT = re.compile(r"(\w+)\s*=\s*(\S+)", re.IGNORECASE)

    def __init__(self, text: str, source_path: str | None = None):
        self.text = text
        self.source_path = source_path
        self.lines = text.splitlines()
        self.line_idx = 0
        self.state = _ParserState.TITLE

        # Results
        self.title_lines: list[str] = []
        self.blocks: dict[str, MdinBlock] = {}
        self.wt_blocks: list[WtBlock] = []
        self.file_redirects: list[FileRedirectEntry] = []

        # Current parsing context
        self.current_block: MdinBlock | None = None
        self.current_wt_block: WtBlock | None = None
        self.in_wt_section = False
        self.wt_ended = False

    def parse(self) -> ParsedMdin:
        """Parse the MDIN text and return structured result."""
        self._process_lines()

        return ParsedMdin(
            title="\n".join(self.title_lines).strip(),
            blocks=self.blocks,
            wt_blocks=self.wt_blocks,
            file_redirects=self.file_redirects,
            source_path=self.source_path,
            raw_text=self.text,
        )

    def _process_lines(self) -> None:
        """Process all lines in the input."""
        accumulated_text = ""

        while self.line_idx < len(self.lines):
            line = self.lines[self.line_idx]
            line_num = self.line_idx + 1  # 1-based line numbers

            # Remove comments (but not inside quotes)
            line = self._strip_comments(line)

            if self.state == _ParserState.TITLE:
                # Look for namelist start
                match = self._NAMELIST_START.search(line)
                if match:
                    # Everything before this is title
                    pre_namelist = line[: match.start()].strip()
                    if pre_namelist:
                        self.title_lines.append(pre_namelist)

                    self._start_namelist(match.group(1), line_num)
                    rest_of_line = line[match.end() :]
                    self.state = _ParserState.IN_BLOCK

                    # Check if namelist ends on the same line (e.g., &wt TYPE='END', /)
                    end_match = self._NAMELIST_END.search(rest_of_line)
                    if end_match:
                        # Process content before the /
                        pre_end = rest_of_line[: end_match.start()]
                        self._process_namelist_content(pre_end, line_num)
                        self._end_namelist(line_num)

                        # Check for another namelist after /
                        post_end = rest_of_line[end_match.end() :]
                        new_match = self._NAMELIST_START.search(post_end)
                        if new_match:
                            self._start_namelist(new_match.group(1), line_num)
                            accumulated_text = post_end[new_match.end() :]
                            self.state = _ParserState.IN_BLOCK
                        elif self.wt_ended:
                            self.state = _ParserState.AFTER_BLOCKS
                        else:
                            self.state = _ParserState.TITLE
                    else:
                        accumulated_text = rest_of_line
                else:
                    # Part of title
                    if line.strip():
                        self.title_lines.append(line.strip())

            elif self.state == _ParserState.IN_BLOCK:
                # Check for namelist end
                if self._NAMELIST_END.search(line):
                    # Process any remaining text before the /
                    end_match = self._NAMELIST_END.search(line)
                    pre_end = line[: end_match.start()]
                    accumulated_text += " " + pre_end
                    self._process_namelist_content(accumulated_text, line_num)
                    accumulated_text = ""

                    # End the current block
                    self._end_namelist(line_num)

                    # Check if there's another namelist on the same line after /
                    post_end = line[end_match.end() :]
                    new_match = self._NAMELIST_START.search(post_end)
                    if new_match:
                        self._start_namelist(new_match.group(1), line_num)
                        accumulated_text = post_end[new_match.end() :]
                        self.state = _ParserState.IN_BLOCK
                    elif self.wt_ended:
                        self.state = _ParserState.AFTER_BLOCKS
                    else:
                        self.state = _ParserState.TITLE
                else:
                    # Check if a new namelist starts (shouldn't happen without /)
                    new_match = self._NAMELIST_START.search(line)
                    if new_match:
                        # Process accumulated and start new block
                        self._process_namelist_content(accumulated_text, line_num)
                        accumulated_text = ""
                        self._end_namelist(line_num)
                        self._start_namelist(new_match.group(1), line_num)
                        accumulated_text = line[new_match.end() :]
                    else:
                        # Accumulate content
                        accumulated_text += " " + line

            elif self.state == _ParserState.AFTER_BLOCKS:
                # Look for file redirects
                self._parse_file_redirect(line, line_num)

            self.line_idx += 1

        # Handle any remaining content
        if accumulated_text.strip() and self.current_block:
            self._process_namelist_content(accumulated_text, self.line_idx)
            self._end_namelist(self.line_idx)

    def _strip_comments(self, line: str) -> str:
        """Remove comments (! to end of line) but not inside quotes."""
        result = []
        in_quote = False
        for i, char in enumerate(line):
            if char == "'" and (i == 0 or line[i - 1] != "\\"):
                in_quote = not in_quote
            elif char == "!" and not in_quote:
                break
            result.append(char)
        return "".join(result)

    def _start_namelist(self, name: str, line_num: int) -> None:
        """Start a new namelist block."""
        name_lower = name.lower()

        if name_lower == "wt":
            self.in_wt_section = True
            self.current_wt_block = WtBlock(
                type_value=None,
                line_start=line_num,
            )
        else:
            self.in_wt_section = False
            self.current_block = MdinBlock(
                section=name_lower,
                raw_section=name,
                line_start=line_num,
            )

    def _end_namelist(self, line_num: int) -> None:
        """End the current namelist block."""
        if self.in_wt_section and self.current_wt_block:
            self.current_wt_block.line_end = line_num
            self.wt_blocks.append(self.current_wt_block)

            # Check if this was TYPE='END'
            if self.current_wt_block.type_value:
                if self.current_wt_block.type_value.upper() == "END":
                    self.wt_ended = True

            self.current_wt_block = None
            self.in_wt_section = False
        elif self.current_block:
            self.current_block.line_end = line_num
            # Store block (merge if section already exists)
            section = self.current_block.section
            if section in self.blocks:
                # Merge entries
                self.blocks[section].entries.update(self.current_block.entries)
                self.blocks[section].line_end = line_num
            else:
                self.blocks[section] = self.current_block
            self.current_block = None

    def _process_namelist_content(self, content: str, line_num: int) -> None:
        """Parse key=value pairs from namelist content."""
        if not content.strip():
            return

        # Normalize whitespace
        content = " ".join(content.split())

        # Find all key=value pairs
        pos = 0
        while pos < len(content):
            match = self._KEY_VALUE.search(content, pos)
            if not match:
                break

            raw_key = match.group(1)
            key = self._normalize_key(raw_key)
            value_start = match.end()

            # Find the value
            value, value_end, raw_value = self._extract_value(content, value_start)

            # Create MdinValue
            mdin_value = MdinValue(
                key=key,
                raw_key=raw_key,
                value=value,
                raw_value=raw_value,
                line_number=line_num,
                section=self.current_block.section if self.current_block else "wt",
            )

            # Store in appropriate block
            if self.in_wt_section and self.current_wt_block:
                if key.lower() == "type":
                    self.current_wt_block.type_value = str(value)
                self.current_wt_block.entries[key] = mdin_value
            elif self.current_block:
                self.current_block.entries[key] = mdin_value

            pos = value_end

    def _normalize_key(self, key: str) -> str:
        """Normalize a key, converting array syntax like mask(1) to mask_1."""
        # Convert array notation: name(n) -> name_n
        array_match = re.match(r"(\w+)\((\d+)\)", key)
        if array_match:
            return f"{array_match.group(1).lower()}_{array_match.group(2)}"
        return key.lower()

    def _extract_value(
        self, content: str, start: int
    ) -> tuple[Any, int, str]:
        """Extract a value starting at position start."""
        content = content[start:].lstrip()
        if not content:
            return None, start, ""

        # Check for string value
        if content.startswith("'"):
            match = self._STRING_VALUE.match(content)
            if match:
                raw = match.group(0)
                return match.group(1), start + len(raw) + (len(content) - len(content.lstrip())), raw
            # Unterminated string - find next ' or end
            end = content.find("'", 1)
            if end == -1:
                end = len(content)
            raw = content[: end + 1] if end < len(content) else content
            return content[1:end], start + len(content) - len(content.lstrip()) + end + 1, raw

        # Find the end of this value (comma, whitespace before next key, or end)
        value_chars = []
        i = 0
        while i < len(content):
            char = content[i]
            if char == ",":
                # Skip the comma
                i += 1
                break
            if char in " \t\n":
                # Check if this is followed by key=
                rest = content[i:].lstrip()
                if self._KEY_VALUE.match(rest):
                    break
            value_chars.append(char)
            i += 1

        raw_value = "".join(value_chars).strip()
        parsed_value = self._parse_value(raw_value)

        return parsed_value, start + len(content) - len(content.lstrip()) + i, raw_value

    def _parse_value(self, raw: str) -> Any:
        """Parse a raw value string into appropriate Python type."""
        raw = raw.strip()
        if not raw:
            return None

        # Check for logical values
        if self._LOGICAL_TRUE.match(raw):
            return True
        if self._LOGICAL_FALSE.match(raw):
            return False

        # Check for integer
        try:
            return int(raw)
        except ValueError:
            pass

        # Check for float
        try:
            return float(raw)
        except ValueError:
            pass

        # Return as string
        return raw

    def _parse_file_redirect(self, line: str, line_num: int) -> None:
        """Parse a file redirect line (after &wt END)."""
        line = line.strip()
        if not line:
            return

        # Look for TYPE = filename
        match = self._FILE_REDIRECT.match(line)
        if match:
            redirect_type = match.group(1).upper()
            filename = match.group(2)

            # Handle quoted filenames
            if filename.startswith("'") and filename.endswith("'"):
                filename = filename[1:-1]
            elif filename.startswith('"') and filename.endswith('"'):
                filename = filename[1:-1]

            self.file_redirects.append(
                FileRedirectEntry(
                    redirect_type=redirect_type,
                    filename=filename,
                    line_number=line_num,
                )
            )


# =============================================================================
# Public API
# =============================================================================


def parse_mdin(text: str, source_path: str | None = None) -> ParsedMdin:
    """
    Parse an MDIN file from text.

    Args:
        text: MDIN file content
        source_path: Optional path for error reporting

    Returns:
        ParsedMdin structure with all parsed content

    Raises:
        MdinSyntaxError: On malformed input with line number
    """
    parser = _MdinParser(text, source_path)
    return parser.parse()


def parse_mdin_file(path: str | Path) -> ParsedMdin:
    """
    Parse an MDIN file from disk.

    Args:
        path: Path to the MDIN file

    Returns:
        ParsedMdin structure with all parsed content

    Raises:
        FileNotFoundError: If file doesn't exist
        MdinSyntaxError: On malformed input
    """
    path = Path(path)
    text = path.read_text()
    return parse_mdin(text, str(path))
