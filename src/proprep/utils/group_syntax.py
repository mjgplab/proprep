"""
Shared parsers for the "1:1-4 2:5,6" group-assignment syntax.

Originally introduced for the redox-site detector's site-grouping prompt.
Reused at the MCPB RESP step to let users group residues into
cross-residue charge-equivalence groups.

Format:
    Each group is "<group_id>:<index_list>"; groups separated by whitespace.
    The index list is comma-separated; ranges use "-" (e.g., 1-4).
    Group IDs must start at 1 and be sequential when used for full
    partition assignments (the detector's case); for sparse-subset
    assignments (the RESP equivalence case) sequentiality is not enforced
    here — caller validates per-feature semantics.

Examples:
    "1:1-4"          -> {1: [0, 1, 2, 3]}
    "1:1,2 2:3,4"    -> {1: [0, 1], 2: [2, 3]}
    "1:1 2:2 3:3"    -> {1: [0], 2: [1], 3: [2]}

All returned indices are 0-based for direct list indexing; the caller
displays them 1-based in the UI.
"""

from typing import Dict, List, Optional


def parse_selection_input(input_str: str, max_num: int) -> Optional[List[int]]:
    """
    Parse a comma/space-separated selection with optional ranges into 0-based indices.

    Returns None if any token is malformed or out-of-range.
    """
    try:
        indices = set()
        parts = input_str.replace(",", " ").split()

        for part in parts:
            if "-" in part:
                range_parts = part.split("-")
                if len(range_parts) != 2:
                    return None
                start = int(range_parts[0])
                end = int(range_parts[1])
                if start > end:
                    return None
                for i in range(start, end + 1):
                    if 1 <= i <= max_num:
                        indices.add(i - 1)
                    else:
                        return None
            else:
                num = int(part)
                if 1 <= num <= max_num:
                    indices.add(num - 1)
                else:
                    return None

        return sorted(indices)

    except ValueError:
        return None


def parse_group_assignments(assignment_str: str, max_index: int) -> Optional[Dict[int, List[int]]]:
    """
    Parse "1:1,2 2:3-5" style group assignments into a dict.

    Args:
        assignment_str: Input string in the documented format.
        max_index: Largest valid 1-based index (e.g., number of residues).

    Returns:
        Dict mapping group_id (int >= 1) -> list of 0-based indices,
        or None if any token is malformed, out-of-range, or duplicate group_id.
    """
    try:
        groups: Dict[int, List[int]] = {}

        for assignment in assignment_str.split():
            if ":" not in assignment:
                return None

            head, tail = assignment.split(":", 1)

            try:
                group_id = int(head)
            except ValueError:
                return None
            if group_id < 1:
                return None
            if group_id in groups:
                return None

            indices = parse_selection_input(tail, max_index)
            if indices is None:
                return None

            groups[group_id] = indices

        return groups

    except Exception:
        return None
