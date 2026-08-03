"""
Phrase parser for relationship definitions.

Parses natural-language-like phrases into structured relationship objects.
Strict keyword-based parsing with immediate validation against the RedoxSite.

Three categories:
  Connectivity: [role] bonds [through [atom]] to [[atom] on] [role] [with [bond_type]]
  Spatial:      [[atom] of] [role] is [closest to/farthest from/within N Å of] [[atom] on] [role]
  Sequence:     [role] is +N/-N from [role]
                [role] has a lower/higher resid than [role]
                [role] is on the same/different chain as [role]
"""

import logging
import re
from typing import Dict, List, Optional, Set, Tuple

from .data_models import (
    ConnectivityPhrase,
    RelationshipPhrase,
    RoleDefinition,
    SequencePhrase,
    SpatialPhrase,
)

logger = logging.getLogger(__name__)

# Valid bond types
BOND_TYPES = {"coordinate", "covalent", "thioether", "disulfide"}


class ParseError(Exception):
    """Raised when a phrase cannot be parsed."""
    pass


class ValidationError(Exception):
    """Raised when a parsed phrase fails validation against the site data."""
    pass


class PhraseParser:
    """
    Parse and validate relationship phrases.

    Args:
        roles: dict of label -> RoleDefinition (from step 2)
        bonds: list of RedoxSiteBond objects (from the RedoxSite)
    """

    def __init__(self, roles: Dict[str, RoleDefinition], bonds: list = None):
        self.roles = roles
        self.bonds = bonds or []
        self._role_labels: Set[str] = set(roles.keys())

    def update_roles(self, roles: Dict[str, RoleDefinition]):
        """Update roles (e.g., after new roles created during editing)."""
        self.roles = roles
        self._role_labels = set(roles.keys())

    def parse(self, text: str) -> RelationshipPhrase:
        """
        Parse a phrase string into a RelationshipPhrase.

        Raises ParseError if the phrase doesn't match any template.
        Raises ValidationError if the phrase references invalid roles/atoms.
        """
        text = text.strip()
        if not text:
            raise ParseError("Empty phrase")

        # Detect category by keywords
        phrase = self._try_parse_connectivity(text)
        if phrase is None:
            phrase = self._try_parse_spatial(text)
        if phrase is None:
            phrase = self._try_parse_sequence(text)
        if phrase is None:
            raise ParseError(
                f"Could not parse phrase: '{text}'\n"
                "Expected one of:\n"
                "  [role] bonds [through [atom]] to [[atom] on] [role] [with [type]]\n"
                "  [[atom] of] [role] is closest to/farthest from/within N Å of [[atom] on] [role]\n"
                "  [role] is +N/-N from [role]\n"
                "  [role] has a lower/higher resid than [role]\n"
                "  [role] is on the same/different chain as [role]"
            )

        self._validate(phrase)
        return phrase

    def get_available_atoms(self, role_label: str) -> List[str]:
        """Get atom names available on a role (for autocomplete)."""
        role = self.roles.get(role_label)
        if role:
            return sorted(role.atom_names)
        return []

    def get_available_bond_types(self) -> List[str]:
        """Get valid bond types (for autocomplete)."""
        return sorted(BOND_TYPES)

    # ===== CONNECTIVITY PARSING =====

    def _try_parse_connectivity(self, text: str) -> Optional[ConnectivityPhrase]:
        """
        Try to parse as connectivity phrase.
        Pattern: [role] bonds [through [atom]] to [[atom] on] [role] [with [bond_type]]
        """
        # Must contain "bonds" and "to"
        if " bonds " not in text and not text.startswith("bonds "):
            return None

        # Split on "bonds"
        parts = text.split(" bonds ", 1)
        if len(parts) != 2:
            return None

        subject_role = parts[0].strip()
        remainder = parts[1].strip()

        source_atom = None
        target_atom = None
        bond_type = None

        # Check for "with [bond_type]" at end
        if " with " in remainder:
            remainder, bond_type_str = remainder.rsplit(" with ", 1)
            bond_type = bond_type_str.strip().lower()
            remainder = remainder.strip()

        # Check for "through [atom]" before "to"
        if remainder.startswith("through "):
            remainder = remainder[len("through "):]
            # Next token is the source atom, then "to"
            if " to " not in remainder:
                raise ParseError(
                    f"Expected 'to' after 'through [atom]' in: {text}"
                )
            source_atom_str, remainder = remainder.split(" to ", 1)
            source_atom = source_atom_str.strip().upper()
            remainder = remainder.strip()
        elif " to " in remainder:
            # No "through" — just "bonds to ..."
            remainder = remainder.split(" to ", 1)[1].strip()
        else:
            raise ParseError(f"Expected 'to' in connectivity phrase: {text}")

        # Remainder is: [[atom] on] [role]
        target_role, target_atom = self._parse_atom_on_role(remainder)

        return ConnectivityPhrase(
            subject_role=subject_role,
            target_role=target_role,
            source_atom=source_atom,
            target_atom=target_atom,
            bond_type=bond_type,
            raw_text=text,
        )

    # ===== SPATIAL PARSING =====

    def _try_parse_spatial(self, text: str) -> Optional[SpatialPhrase]:
        """
        Try to parse as spatial phrase.
        Pattern: [[atom] of] [role] is [closest to/farthest from/within N Å of] [[atom] on] [role]
        """
        # Must contain "is closest to", "is farthest from", or "is within"
        comparison = None
        distance = None

        for comp, keyword in [
            ("closest", " is closest to "),
            ("farthest", " is farthest from "),
            ("within", " is within "),
        ]:
            if keyword in text:
                comparison = comp
                left, right = text.split(keyword, 1)
                break
        else:
            return None

        left = left.strip()
        right = right.strip()

        # Parse left side: [[atom] of] [role]
        source_atom = None
        if " of " in left:
            source_atom_str, subject_role = left.rsplit(" of ", 1)
            source_atom = source_atom_str.strip().upper()
            subject_role = subject_role.strip()
        else:
            subject_role = left

        # For "within", parse "N Å of [target]"
        if comparison == "within":
            # Pattern: N Å of [[atom] on] [role]
            match = re.match(r'([\d.]+)\s*[Åå]?\s*of\s+(.*)', right)
            if not match:
                raise ParseError(
                    f"Expected 'N Å of [role]' after 'is within' in: {text}"
                )
            distance = float(match.group(1))
            right = match.group(2).strip()

        # Parse right side: [[atom] on] [role]
        target_role, target_atom = self._parse_atom_on_role(right)

        return SpatialPhrase(
            subject_role=subject_role,
            target_role=target_role,
            source_atom=source_atom,
            target_atom=target_atom,
            comparison=comparison,
            distance=distance,
            raw_text=text,
        )

    # ===== SEQUENCE PARSING =====

    def _try_parse_sequence(self, text: str) -> Optional[SequencePhrase]:
        """
        Try to parse as sequence phrase.
        Forms:
          [role] is +N/-N from [role]
          [role] has a lower/higher resid than [role]
          [role] is on the same/different chain as [role]
        """
        # Form 1: arithmetic offset
        match = re.match(
            r'(\S+)\s+is\s+([+-]\d+)\s+from\s+(\S+)$', text
        )
        if match:
            return SequencePhrase(
                subject_role=match.group(1),
                offset=int(match.group(2)),
                offset_target_role=match.group(3),
                raw_text=text,
            )

        # Form 2: resid comparison
        match = re.match(
            r'(\S+)\s+has\s+a\s+(lower|higher)\s+resid\s+than\s+(\S+)$', text
        )
        if match:
            return SequencePhrase(
                subject_role=match.group(1),
                resid_comparison=match.group(2),
                comparison_target_role=match.group(3),
                raw_text=text,
            )

        # Form 3: chain relation
        match = re.match(
            r'(\S+)\s+is\s+on\s+the\s+(same|different)\s+chain\s+as\s+(\S+)$',
            text,
        )
        if match:
            return SequencePhrase(
                subject_role=match.group(1),
                chain_relation=match.group(2),
                chain_target_role=match.group(3),
                raw_text=text,
            )

        return None

    # ===== HELPERS =====

    def _parse_atom_on_role(self, text: str) -> Tuple[str, Optional[str]]:
        """
        Parse '[[atom] on] [role]' pattern.

        Returns (role, atom_or_None).
        'CAB on center' -> ('center', 'CAB')
        'center'        -> ('center', None)
        """
        text = text.strip()
        if " on " in text:
            atom_str, role = text.rsplit(" on ", 1)
            return role.strip(), atom_str.strip().upper()
        return text, None

    # ===== VALIDATION =====

    def _validate(self, phrase: RelationshipPhrase):
        """Validate a parsed phrase against known roles and atoms."""
        if isinstance(phrase, ConnectivityPhrase):
            self._validate_connectivity(phrase)
        elif isinstance(phrase, SpatialPhrase):
            self._validate_spatial(phrase)
        elif isinstance(phrase, SequencePhrase):
            self._validate_sequence(phrase)

    def _validate_role(self, label: str, context: str):
        """Validate that a role label exists."""
        if label not in self._role_labels:
            available = ", ".join(sorted(self._role_labels))
            raise ValidationError(
                f"Unknown role '{label}' in {context}. "
                f"Available roles: {available}"
            )

    def _validate_atom_on_role(self, atom: Optional[str], role_label: str, context: str):
        """Validate that an atom exists on a role."""
        if atom is None:
            return
        role = self.roles.get(role_label)
        if role and atom not in role.atom_names:
            available = ", ".join(sorted(role.atom_names))
            raise ValidationError(
                f"Atom '{atom}' not found on role '{role_label}' in {context}. "
                f"Available atoms: {available}"
            )

    def _validate_connectivity(self, phrase: ConnectivityPhrase):
        """Validate connectivity phrase."""
        self._validate_role(phrase.subject_role, phrase.raw_text)
        self._validate_role(phrase.target_role, phrase.raw_text)
        self._validate_atom_on_role(
            phrase.source_atom, phrase.subject_role, phrase.raw_text
        )
        self._validate_atom_on_role(
            phrase.target_atom, phrase.target_role, phrase.raw_text
        )
        if phrase.bond_type and phrase.bond_type not in BOND_TYPES:
            raise ValidationError(
                f"Unknown bond type '{phrase.bond_type}' in {phrase.raw_text}. "
                f"Valid types: {', '.join(sorted(BOND_TYPES))}"
            )

    def _validate_spatial(self, phrase: SpatialPhrase):
        """Validate spatial phrase."""
        self._validate_role(phrase.subject_role, phrase.raw_text)
        self._validate_role(phrase.target_role, phrase.raw_text)
        self._validate_atom_on_role(
            phrase.source_atom, phrase.subject_role, phrase.raw_text
        )
        self._validate_atom_on_role(
            phrase.target_atom, phrase.target_role, phrase.raw_text
        )
        if phrase.comparison == "within" and phrase.distance is None:
            raise ValidationError(
                f"Missing distance for 'within' comparison in {phrase.raw_text}"
            )
        if phrase.distance is not None and phrase.distance <= 0:
            raise ValidationError(
                f"Distance must be positive in {phrase.raw_text}"
            )

    def _validate_sequence(self, phrase: SequencePhrase):
        """Validate sequence phrase."""
        self._validate_role(phrase.subject_role, phrase.raw_text)
        if phrase.offset_target_role:
            self._validate_role(phrase.offset_target_role, phrase.raw_text)
        if phrase.comparison_target_role:
            self._validate_role(phrase.comparison_target_role, phrase.raw_text)
        if phrase.chain_target_role:
            self._validate_role(phrase.chain_target_role, phrase.raw_text)
