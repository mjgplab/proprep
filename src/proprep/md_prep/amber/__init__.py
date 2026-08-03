"""
AMBER MDIN Configuration Interface

A comprehensive interface for building, parsing, validating, and exploring
AMBER MDIN input files, backed by a keyword database parsed from the AMBER Manual.

Layers:
- Layer 0 (KeywordDB): Fast indexed access to keyword metadata
- Layer 1 (Parser/Writer): MDIN file I/O
- Layer 2 (Validator): Configuration validation with dependency checking
- Layer 3 (Builder): Interactive TUI and programmatic builder

Usage Examples:

    # Keyword lookup
    from proprep.md_prep.amber import get_keyword_db

    db = get_keyword_db()
    ntt = db.get("ntt")
    print(f"ntt: {ntt.description}")
    print(f"Options: {ntt.options}")

    # Parse an existing MDIN file
    from proprep.md_prep.amber import parse_mdin_file

    parsed = parse_mdin_file("production.mdin")
    print(f"nstlim = {parsed.get_value('nstlim')}")

    # Generate a new MDIN file
    from proprep.md_prep.amber import write_mdin

    text = write_mdin(
        {
            "cntrl": {
                "imin": 0,
                "nstlim": 500000,
                "dt": 0.002,
                "ntt": 3,
                "gamma_ln": 1.0,
                "temp0": 300.0,
            }
        },
        title="Production MD",
        comments=True,
    )

    # Validate configuration
    from proprep.md_prep.amber import validate_mdin

    result = validate_mdin(parsed)
    for issue in result.errors:
        print(f"ERROR: {issue.message}")

    # Interactive builder (TUI)
    from proprep.md_prep.amber import run_mdin_builder

    run_mdin_builder()

    # Programmatic builder with fluent API
    from proprep.md_prep.amber import MdinBuilder

    builder = MdinBuilder()
    builder.set("imin", 0).set("nstlim", 500000).set("dt", 0.002)
    builder.set("ntt", 3).set("gamma_ln", 1.0).set("temp0", 300.0)
    result = builder.validate()
    if result.is_valid:
        builder.write("production.mdin")
"""

# Exceptions
from .exceptions import AmbiguousKeywordError, MdinSyntaxError

# Data classes from database
from ._amber_keywords import (
    FileRedirect,
    Keyword,
    WtParameter,
    WtType,
)

# Layer 0: Keyword Database
from .keyword_db import KeywordDB, get_keyword_db

# Layer 1: Parser
from .mdin_parser import (
    FileRedirectEntry,
    MdinBlock,
    MdinValue,
    ParsedMdin,
    WtBlock,
    parse_mdin,
    parse_mdin_file,
)

# Layer 1: Writer
from .mdin_writer import MdinWriter, write_mdin

# Layer 2: Validator
from .mdin_validator import (
    DependencyRule,
    MdinValidator,
    Severity,
    ValidationIssue,
    ValidationResult,
    validate_mdin,
    validate_mdin_file,
)

# Layer 3: Interactive Builder
from .mdin_builder import (
    BuilderState,
    MdinBuilder,
    run_mdin_builder,
)

__all__ = [
    # Exceptions
    "MdinSyntaxError",
    "AmbiguousKeywordError",
    # Database types
    "Keyword",
    "WtType",
    "WtParameter",
    "FileRedirect",
    # Layer 0
    "KeywordDB",
    "get_keyword_db",
    # Layer 1: Parser types
    "MdinValue",
    "MdinBlock",
    "WtBlock",
    "FileRedirectEntry",
    "ParsedMdin",
    # Layer 1: Parser functions
    "parse_mdin",
    "parse_mdin_file",
    # Layer 1: Writer
    "MdinWriter",
    "write_mdin",
    # Layer 2: Validator
    "Severity",
    "ValidationIssue",
    "ValidationResult",
    "DependencyRule",
    "MdinValidator",
    "validate_mdin",
    "validate_mdin_file",
    # Layer 3: Interactive Builder
    "BuilderState",
    "MdinBuilder",
    "run_mdin_builder",
]
