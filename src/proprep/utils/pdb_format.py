"""PDB fixed-column conventions for atom names and elements.

The PDB atom-name field (columns 13-16) encodes the element by where the name
starts, not just by its letters:

* a **two-letter** element starts in column 13 — ``FE1 `` is iron, ``CA  `` a
  calcium ion, ``MO  `` molybdenum;
* a **one-letter** element is right-justified into column 14, leaving column 13
  blank so the remaining columns can carry a remoteness indicator — `` CA `` is
  an alpha carbon, `` PA `` a phosphorus, `` S1 `` a sulfur.

That single column is the only thing distinguishing the calcium ion ``CA`` from
the alpha carbon ``CA``. Stripping the name before looking at it throws the
distinction away, which is how FAD's phosphates PA/PB have been read as
protactinium and lead, and how iron written from column 14 has been read as
fluorine.

Both directions live here so the writing and reading sides cannot drift apart.
"""

from typing import Optional

__all__ = ["atom_name_field", "element_from_name_field"]


def atom_name_field(atom_name: str, element: str = "") -> str:
    """The 4-character columns-13-16 field for an atom name.

    Args:
        atom_name: Atom name, with or without surrounding spaces.
        element: Element symbol, if known. Its length decides the
            justification; when it is unknown the name is indented (the
            one-letter convention), which is the safe default because a
            stripped two-letter name is ambiguous.

    Returns:
        Exactly four characters, ready to place at column 13.
    """
    name = (atom_name or "").strip()
    symbol = (element or "").strip()

    # A four-character name has nowhere to indent to and always starts at 13.
    if len(name) >= 4 or len(symbol) == 2:
        return f"{name:<4.4s}"
    return f" {name:<3.3s}"


def element_from_name_field(name_field: str) -> str:
    """Element implied by a raw columns-13-16 field, in title case.

    Args:
        name_field: The four characters as they appear in the file — NOT a
            stripped name, which no longer carries the justification. With
            BioPython this is ``Atom.fullname``.

    Returns:
        Title-case element symbol, or "" if the field is empty.
    """
    field = name_field or ""
    if len(field) >= 4:
        if field[0] == " ":
            return field.strip()[:1].title()
        return field[:2].strip().title()
    return field.strip().title()
