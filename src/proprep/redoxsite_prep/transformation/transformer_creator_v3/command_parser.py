"""
Command parser for PDB editing operations.

Parses role-based editing commands:
  rename [role] [new_resname]                    # rename residue
  rename [role]:[old_atoms] [new_atoms]          # rename atoms (positional 1:1)
  move [role]:[atoms] to [role]                  # move atoms to existing role
  move [role]:[atoms] to new [resname] as [role] # move atoms, create new residue+role
  hetatm [role]                                  # whole residue → HETATM
  atom [role]                                    # whole residue → ATOM
"""

import logging
from typing import Dict, List, Optional, Set

from .data_models import CommandType, EditCommand, RoleDefinition

logger = logging.getLogger(__name__)


class CommandParseError(Exception):
    """Raised when a command cannot be parsed."""
    pass


class CommandValidationError(Exception):
    """Raised when a parsed command fails validation."""
    pass


class CommandParser:
    """
    Parse and validate PDB editing commands.

    Args:
        roles: dict of label -> RoleDefinition
        site_state: current atom state per role (tracks moves/renames across commands)
    """

    def __init__(self, roles: Dict[str, RoleDefinition]):
        self.roles = roles
        # Track current atom names per role (mutated as commands execute)
        self._role_atoms: Dict[str, Set[str]] = {}
        for label, role in roles.items():
            self._role_atoms[label] = set(role.atom_names)

    def register_new_role(self, label: str, resname: str, atom_names: List[str]):
        """Register a new role created via 'move ... to new ... as ...'."""
        self._role_atoms[label] = set(atom_names)

    def get_role_atoms(self, label: str) -> Set[str]:
        """Get current atom names for a role."""
        return self._role_atoms.get(label, set())

    def parse(self, text: str) -> EditCommand:
        """
        Parse a command string into an EditCommand.

        Raises CommandParseError if the command doesn't match any template.
        Raises CommandValidationError if references are invalid.
        """
        text = text.strip()
        if not text:
            raise CommandParseError("Empty command")

        tokens = text.split()
        verb = tokens[0].lower()

        if verb == "rename":
            cmd = self._parse_rename(text, tokens)
        elif verb == "move":
            cmd = self._parse_move(text, tokens)
        elif verb == "hetatm":
            cmd = self._parse_record_type(text, tokens, CommandType.HETATM)
        elif verb == "atom":
            cmd = self._parse_record_type(text, tokens, CommandType.ATOM)
        else:
            raise CommandParseError(
                f"Unknown command verb '{verb}'. "
                "Expected: rename, move, hetatm, atom"
            )

        self._validate(cmd)
        return cmd

    # ===== RENAME PARSING =====

    def _parse_rename(self, text: str, tokens: List[str]) -> EditCommand:
        """
        Parse rename command.
        Two forms:
          rename [role] [new_resname]
          rename [role]:[old_atoms] [new_atoms]
        """
        if len(tokens) < 3:
            raise CommandParseError(
                "rename requires: rename [role] [new_resname] "
                "or rename [role]:[atoms] [new_atoms]"
            )

        target = tokens[1]

        if ":" in target:
            # Atom rename: rename role:OLD1,OLD2 NEW1,NEW2
            role_label, old_atoms_str = target.split(":", 1)
            old_atoms = [a.strip().upper() for a in old_atoms_str.split(",")]
            new_atoms = [a.strip().upper() for a in tokens[2].split(",")]

            return EditCommand(
                command_type=CommandType.RENAME_ATOMS,
                raw_text=text,
                source_role=role_label,
                old_atom_names=old_atoms,
                new_atom_names=new_atoms,
            )
        else:
            # Residue rename: rename role NEW_RESNAME
            return EditCommand(
                command_type=CommandType.RENAME_RESIDUE,
                raw_text=text,
                source_role=target,
                new_resname=tokens[2].upper(),
            )

    # ===== MOVE PARSING =====

    def _parse_move(self, text: str, tokens: List[str]) -> EditCommand:
        """
        Parse move command.
        Two forms:
          move [role]:[atoms] to [role]
          move [role]:[atoms] to new [resname] as [role]
        """
        # Find "to" keyword
        if " to " not in text:
            raise CommandParseError("move command requires 'to': move [role]:[atoms] to [target]")

        left, right = text.split(" to ", 1)
        left = left.strip()
        right = right.strip()

        # Parse source: "move role:ATOM1,ATOM2"
        source_part = left.split(None, 1)[1] if " " in left else left  # skip "move"
        if ":" not in source_part:
            raise CommandParseError(
                "move source must specify atoms: move [role]:[atoms] to ..."
            )
        role_label, atoms_str = source_part.split(":", 1)
        atom_names = [a.strip().upper() for a in atoms_str.split(",")]

        # Parse target: either "[role]" or "new [resname] as [role]"
        right_tokens = right.split()

        if right_tokens[0].lower() == "new":
            # move ... to new RESNAME as ROLE
            if len(right_tokens) < 4 or right_tokens[2].lower() != "as":
                raise CommandParseError(
                    "Expected: move [role]:[atoms] to new [resname] as [new_role]"
                )
            return EditCommand(
                command_type=CommandType.MOVE_ATOMS_NEW,
                raw_text=text,
                source_role=role_label,
                atom_names=atom_names,
                new_residue_name=right_tokens[1].upper(),
                new_role_label=right_tokens[3],
            )
        else:
            # move ... to [role]
            target_role = right_tokens[0]
            return EditCommand(
                command_type=CommandType.MOVE_ATOMS,
                raw_text=text,
                source_role=role_label,
                atom_names=atom_names,
                target_role=target_role,
            )

    # ===== RECORD TYPE PARSING =====

    def _parse_record_type(
        self, text: str, tokens: List[str], cmd_type: CommandType
    ) -> EditCommand:
        """Parse hetatm/atom command: hetatm [role] or atom [role]."""
        if len(tokens) < 2:
            raise CommandParseError(f"{tokens[0]} requires a role: {tokens[0]} [role]")

        return EditCommand(
            command_type=cmd_type,
            raw_text=text,
            source_role=tokens[1],
        )

    # ===== VALIDATION =====

    def _validate(self, cmd: EditCommand):
        """Validate a parsed command against current state."""
        # Validate source role exists
        if cmd.source_role not in self._role_atoms:
            available = ", ".join(sorted(self._role_atoms.keys()))
            raise CommandValidationError(
                f"Unknown role '{cmd.source_role}'. Available: {available}"
            )

        source_atoms = self._role_atoms[cmd.source_role]

        if cmd.command_type == CommandType.RENAME_ATOMS:
            # Validate old atoms exist
            for atom in cmd.old_atom_names:
                if atom not in source_atoms:
                    available = ", ".join(sorted(source_atoms))
                    raise CommandValidationError(
                        f"Atom '{atom}' not found on role '{cmd.source_role}'. "
                        f"Available: {available}"
                    )
            # Validate list lengths match
            if len(cmd.old_atom_names) != len(cmd.new_atom_names):
                raise CommandValidationError(
                    f"Atom rename lists must have same length: "
                    f"{len(cmd.old_atom_names)} old vs {len(cmd.new_atom_names)} new"
                )
            # Check for duplicate new names
            # (after renaming, considering atoms that aren't being renamed)
            remaining = source_atoms - set(cmd.old_atom_names)
            all_new = remaining | set(cmd.new_atom_names)
            if len(all_new) < len(remaining) + len(cmd.new_atom_names):
                raise CommandValidationError(
                    f"Atom rename would create duplicate names on '{cmd.source_role}'"
                )

        elif cmd.command_type in (CommandType.MOVE_ATOMS, CommandType.MOVE_ATOMS_NEW):
            # Validate atoms exist on source
            for atom in cmd.atom_names:
                if atom not in source_atoms:
                    available = ", ".join(sorted(source_atoms))
                    raise CommandValidationError(
                        f"Atom '{atom}' not found on role '{cmd.source_role}'. "
                        f"Available: {available}"
                    )

            if cmd.command_type == CommandType.MOVE_ATOMS:
                # Validate target role exists
                if cmd.target_role not in self._role_atoms:
                    available = ", ".join(sorted(self._role_atoms.keys()))
                    raise CommandValidationError(
                        f"Unknown target role '{cmd.target_role}'. Available: {available}"
                    )
                # Check for atom name collision on target
                target_atoms = self._role_atoms[cmd.target_role]
                collisions = set(cmd.atom_names) & target_atoms
                if collisions:
                    raise CommandValidationError(
                        f"Moving atoms {sorted(collisions)} to '{cmd.target_role}' "
                        f"would create duplicate atom names. Rename them first."
                    )

            elif cmd.command_type == CommandType.MOVE_ATOMS_NEW:
                # Validate new role label doesn't already exist
                if cmd.new_role_label in self._role_atoms:
                    raise CommandValidationError(
                        f"Role '{cmd.new_role_label}' already exists. "
                        "Choose a different label."
                    )

    def apply(self, cmd: EditCommand):
        """
        Apply a command to the internal state tracking.
        Call this after parse() succeeds to keep state consistent.
        """
        if cmd.command_type == CommandType.RENAME_ATOMS:
            atoms = self._role_atoms[cmd.source_role]
            for old, new in zip(cmd.old_atom_names, cmd.new_atom_names):
                atoms.discard(old)
                atoms.add(new)

        elif cmd.command_type == CommandType.MOVE_ATOMS:
            source = self._role_atoms[cmd.source_role]
            target = self._role_atoms[cmd.target_role]
            for atom in cmd.atom_names:
                source.discard(atom)
                target.add(atom)

        elif cmd.command_type == CommandType.MOVE_ATOMS_NEW:
            source = self._role_atoms[cmd.source_role]
            new_atoms = set()
            for atom in cmd.atom_names:
                source.discard(atom)
                new_atoms.add(atom)
            self._role_atoms[cmd.new_role_label] = new_atoms

        # RENAME_RESIDUE, HETATM, ATOM don't change atom tracking
