"""Helpers for emitting tLEaP-safe tokens.

tLEaP's lexer parses a bare token like ``9E2`` as a number in scientific
notation (9x10^2). Two consequences bite ProPrep whenever a residue/molecule
code is used verbatim in a generated tLEaP script:

* As a **unit variable** -- ``9e2 = loadmol2 "..."`` is "assign to a number",
  a fatal syntax error.
* As an **unquoted filename** -- ``loadamberparams 9e2.frcmod`` emits a
  non-fatal "Decimal point follows exponent in NUMBER-like token" error that
  still bumps tLEaP's error count (which callers may treat as failure).

Digit-leading codes are common in the PDB Chemical Component Dictionary
(9E2, 0G6, 1N7) and in some modified/non-standard residues, so any code that
turns a residue name into a tLEaP token should route through here.
"""

from typing import Any


def tleap_safe_unit_var(name: str) -> str:
    """Return a tLEaP-safe unit variable name derived from ``name``.

    tLEaP identifiers must start with a letter; a leading digit makes the token
    number-like. When ``name`` isn't already a safe identifier we prefix it
    with ``m`` and replace any non-alphanumeric/underscore characters.

    CAUTION: use this only for a *throwaway* handle (e.g. the variable in a
    ``X = loadmol2 ...; saveAmberParm X ...`` step). Do NOT use it for the
    variable you ``saveOff`` into a reusable library that must later match a
    structure residue: tLEaP matches ``loadpdb`` residues to templates by the
    OFF/unit **entry name** (which is the variable name at ``saveOff`` time),
    not by the residue name inside the mol2. Mangling ``9E2`` to ``m9e2`` here
    makes the lib entry ``m9e2``, which then fails to match a structure residue
    ``9E2`` ("Unknown residue"). For a residue name that must survive into a
    library, validate it up front with :func:`is_tleap_safe_resname` instead.

    >>> tleap_safe_unit_var("9E2")
    'm9E2'
    >>> tleap_safe_unit_var("LIG")
    'LIG'
    """
    safe = name if name else "mol"
    if not safe[0].isalpha():
        safe = "m" + safe
    safe = "".join(c if (c.isalnum() or c == "_") else "_" for c in safe)
    return safe


def is_tleap_safe_resname(name: str) -> bool:
    """True if ``name`` works as BOTH a tLEaP unit/OFF-entry name and a PDB resName.

    tLEaP matches ``loadpdb`` residues to templates by the OFF/unit entry name
    (not the residue name recorded inside a mol2), and it lexes a digit-leading
    bare token as scientific notation. PDB ``resName`` occupies only three
    columns. So a name is safe only when it is 1-3 characters long, starts with
    a letter, and is otherwise alphanumeric.

    >>> is_tleap_safe_resname("LIG"), is_tleap_safe_resname("9E2"), is_tleap_safe_resname("ABCD")
    (True, False, False)
    """
    if not name or len(name) > 3:
        return False
    if not name[0].isalpha():
        return False
    return all(c.isalnum() for c in name)


def suggest_tleap_safe_resname(name: str) -> str:
    """Best-effort tLEaP/PDB-safe residue name derived from ``name``.

    Strips non-alphanumerics, upper-cases, guarantees a leading letter (prefixes
    ``X`` when the cleaned name starts with a digit), and truncates to three
    characters. Falls back to ``LIG`` when nothing usable remains.

    >>> suggest_tleap_safe_resname("9E2")
    'X9E'
    >>> suggest_tleap_safe_resname("lig-1")
    'LIG'
    """
    cleaned = "".join(c for c in (name or "") if c.isalnum()).upper()
    if not cleaned:
        return "LIG"
    if not cleaned[0].isalpha():
        cleaned = "X" + cleaned
    return cleaned[:3]


def tleap_quote(path: Any) -> str:
    """Wrap a filename/path in double quotes for tLEaP so a number-like name
    (e.g. ``9E2.frcmod``) isn't misparsed as a scientific-notation token.

    Idempotent: a value that is already double-quoted is returned unchanged.
    """
    s = str(path)
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s
    return f'"{s}"'
