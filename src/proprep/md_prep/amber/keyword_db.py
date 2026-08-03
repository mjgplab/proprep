"""
AMBER MDIN Configuration Interface - Layer 0: Keyword Database

Provides fast, indexed access to AMBER keyword metadata from the
comprehensive keyword database. All lookups are case-insensitive.

Usage:
    from proprep.md_prep.amber import get_keyword_db, KeywordDB

    db = get_keyword_db()
    kw = db.get("nstlim")  # Look up by name
    cntrl_keywords = db.by_section("cntrl")  # All keywords in &cntrl
    common = db.commonly_changed()  # Keywords marked commonly_changed=True
    results = db.search("thermostat")  # Full-text search
"""

import logging
from typing import Any

from ._amber_keywords import (
    KEYWORDS,
    WT_TYPES,
    WT_PARAMETERS,
    FILE_REDIRECTIONS,
    Keyword,
    WtType,
    WtParameter,
    FileRedirect,
)
from .exceptions import AmbiguousKeywordError

logger = logging.getLogger(__name__)


class KeywordDB:
    """
    Singleton-pattern wrapper around the AMBER keyword database.

    Builds indexes on initialization for O(1) lookups. All lookups are
    case-insensitive to match AMBER namelist behavior.
    """

    def __init__(
        self,
        keywords: list[Keyword] | None = None,
        wt_types: list[WtType] | None = None,
        wt_params: list[WtParameter] | None = None,
        file_redirects: list[FileRedirect] | None = None,
    ):
        """
        Initialize the keyword database with indexes.

        Args:
            keywords: List of Keyword objects. Defaults to KEYWORDS from database.
            wt_types: List of WtType objects. Defaults to WT_TYPES from database.
            wt_params: List of WtParameter objects. Defaults to WT_PARAMETERS.
            file_redirects: List of FileRedirect objects. Defaults to FILE_REDIRECTIONS.
        """
        self._keywords = keywords if keywords is not None else KEYWORDS
        self._wt_types_list = wt_types if wt_types is not None else WT_TYPES
        self._wt_params = wt_params if wt_params is not None else WT_PARAMETERS
        self._file_redirects_list = (
            file_redirects if file_redirects is not None else FILE_REDIRECTIONS
        )

        # Build indexes
        self._by_name_section: dict[tuple[str, str], Keyword] = {}
        self._by_name: dict[str, list[Keyword]] = {}
        self._by_section: dict[str, list[Keyword]] = {}
        self._by_category: dict[str, list[Keyword]] = {}
        self._categories_by_section: dict[str, list[str]] = {}
        self._commonly_changed_list: list[Keyword] = []
        self._related_graph: dict[str, set[str]] = {}

        # &wt and file redirect indexes
        self._wt_types_by_name: dict[str, WtType] = {}
        self._file_redirects_by_name: dict[str, FileRedirect] = {}

        self._build_indexes()

    def _build_indexes(self) -> None:
        """Build all lookup indexes from the keyword list."""
        for kw in self._keywords:
            name_lower = kw.name.lower()
            section_lower = kw.section.lower()

            # Primary index: (name, section) -> Keyword
            self._by_name_section[(name_lower, section_lower)] = kw

            # Name-only index
            if name_lower not in self._by_name:
                self._by_name[name_lower] = []
            self._by_name[name_lower].append(kw)

            # Section index
            if section_lower not in self._by_section:
                self._by_section[section_lower] = []
            self._by_section[section_lower].append(kw)

            # Category index
            if kw.category:
                cat_lower = kw.category.lower()
                if cat_lower not in self._by_category:
                    self._by_category[cat_lower] = []
                self._by_category[cat_lower].append(kw)

                # Track categories per section
                if section_lower not in self._categories_by_section:
                    self._categories_by_section[section_lower] = []
                if kw.category not in self._categories_by_section[section_lower]:
                    self._categories_by_section[section_lower].append(kw.category)

            # Commonly changed subset
            if kw.commonly_changed:
                self._commonly_changed_list.append(kw)

            # Build bidirectional related graph
            if name_lower not in self._related_graph:
                self._related_graph[name_lower] = set()
            for related_name in kw.related:
                related_lower = related_name.lower()
                self._related_graph[name_lower].add(related_lower)
                # Add reverse link
                if related_lower not in self._related_graph:
                    self._related_graph[related_lower] = set()
                self._related_graph[related_lower].add(name_lower)

        # Sort categories alphabetically within each section
        for section in self._categories_by_section:
            self._categories_by_section[section].sort()

        # Build &wt type index
        for wt in self._wt_types_list:
            self._wt_types_by_name[wt.name.lower()] = wt

        # Build file redirect index
        for fr in self._file_redirects_list:
            self._file_redirects_by_name[fr.name.lower()] = fr

    # -------------------------------------------------------------------------
    # Exact lookups
    # -------------------------------------------------------------------------

    def get(self, name: str, section: str | None = None) -> Keyword | None:
        """
        Look up a keyword by exact name (case-insensitive).

        Args:
            name: Keyword name
            section: Namelist section (e.g., "cntrl", "ewald"). Required if
                    the keyword exists in multiple sections.

        Returns:
            Keyword object if found, None otherwise.

        Raises:
            AmbiguousKeywordError: If name exists in multiple sections and
                                   section is not specified.
        """
        name_lower = name.lower()

        if section is not None:
            return self._by_name_section.get((name_lower, section.lower()))

        keywords = self._by_name.get(name_lower, [])
        if len(keywords) == 0:
            return None
        if len(keywords) == 1:
            return keywords[0]

        # Multiple sections contain this keyword
        sections = [kw.section for kw in keywords]
        raise AmbiguousKeywordError(name, sections)

    def get_all(self, name: str) -> list[Keyword]:
        """
        Return all keywords with this name across all sections.

        Args:
            name: Keyword name (case-insensitive)

        Returns:
            List of Keyword objects (may be empty)
        """
        return list(self._by_name.get(name.lower(), []))

    # -------------------------------------------------------------------------
    # Collection lookups
    # -------------------------------------------------------------------------

    def by_section(self, section: str) -> list[Keyword]:
        """
        All keywords in a namelist section, in definition order.

        Args:
            section: Namelist section name (e.g., "cntrl", "ewald")

        Returns:
            List of Keyword objects
        """
        return list(self._by_section.get(section.lower(), []))

    def by_category(self, category: str) -> list[Keyword]:
        """
        All keywords in a category, in definition order.

        Args:
            category: Category name (case-insensitive)

        Returns:
            List of Keyword objects
        """
        return list(self._by_category.get(category.lower(), []))

    def categories(self, section: str | None = None) -> list[str]:
        """
        All category names, sorted alphabetically.

        Args:
            section: If provided, only categories within that section.

        Returns:
            List of category names
        """
        if section is not None:
            return list(self._categories_by_section.get(section.lower(), []))

        # All unique categories across all sections
        all_categories = set()
        for cats in self._categories_by_section.values():
            all_categories.update(cats)
        return sorted(all_categories)

    def sections(self) -> list[str]:
        """
        All namelist section names, sorted.

        Returns:
            List of section names (e.g., ["cntrl", "ewald", "pb", ...])
        """
        return sorted(self._by_section.keys())

    def commonly_changed(self, section: str | None = None) -> list[Keyword]:
        """
        Keywords flagged commonly_changed=True.

        Args:
            section: If provided, filter by section.

        Returns:
            List of Keyword objects
        """
        if section is None:
            return list(self._commonly_changed_list)
        section_lower = section.lower()
        return [kw for kw in self._commonly_changed_list if kw.section.lower() == section_lower]

    # -------------------------------------------------------------------------
    # Relationship traversal
    # -------------------------------------------------------------------------

    def related(
        self, name: str, section: str | None = None, depth: int = 1
    ) -> list[Keyword]:
        """
        Follow the `related` links from a keyword.

        Args:
            name: Starting keyword name
            section: Section for the starting keyword (required if ambiguous)
            depth: 1 for direct relatives only, 2+ for transitive closure

        Returns:
            Deduplicated list of related Keywords, excluding the starting keyword.
        """
        # Get the starting keyword to verify it exists
        start_kw = self.get(name, section)
        if start_kw is None:
            return []

        name_lower = name.lower()
        visited: set[str] = {name_lower}
        current_level: set[str] = {name_lower}

        for _ in range(depth):
            next_level: set[str] = set()
            for kw_name in current_level:
                for related_name in self._related_graph.get(kw_name, set()):
                    if related_name not in visited:
                        next_level.add(related_name)
                        visited.add(related_name)
            current_level = next_level
            if not current_level:
                break

        # Convert to Keyword objects, excluding the starting keyword
        visited.discard(name_lower)
        result: list[Keyword] = []
        for kw_name in visited:
            keywords = self._by_name.get(kw_name, [])
            result.extend(keywords)

        return result

    # -------------------------------------------------------------------------
    # Search
    # -------------------------------------------------------------------------

    def search(self, query: str, max_results: int = 20) -> list[Keyword]:
        """
        Full-text fuzzy search across keyword names, descriptions, notes,
        and option descriptions.

        Search strategy:
        1. Exact name match (highest priority)
        2. Name substring/prefix match
        3. Token overlap in description + notes
        4. Fuzzy name match (edit distance <= 2, for typo correction)

        Args:
            query: Search query string
            max_results: Maximum number of results to return

        Returns:
            List of Keywords ranked by relevance
        """
        if not query:
            return []

        query_lower = query.lower()
        query_tokens = set(query_lower.split())

        scored_results: list[tuple[float, Keyword]] = []

        for kw in self._keywords:
            score = 0.0
            name_lower = kw.name.lower()

            # Exact name match (highest priority)
            if name_lower == query_lower:
                score = 100.0
            # Name prefix match
            elif name_lower.startswith(query_lower):
                score = 80.0
            # Name substring match
            elif query_lower in name_lower:
                score = 60.0
            else:
                # Token overlap in description and notes
                text = (kw.description or "").lower()
                if kw.notes:
                    text += " " + kw.notes.lower()
                if kw.options:
                    for opt_desc in kw.options.values():
                        if isinstance(opt_desc, str):
                            text += " " + opt_desc.lower()

                text_tokens = set(text.split())
                overlap = len(query_tokens & text_tokens)
                if overlap > 0:
                    score = 10.0 + (overlap * 5.0)

                # Fuzzy name match (simple Levenshtein approximation)
                if score == 0 and len(name_lower) <= 10 and len(query_lower) <= 10:
                    dist = self._edit_distance(name_lower, query_lower)
                    if dist <= 2:
                        score = 5.0 - dist

            if score > 0:
                scored_results.append((score, kw))

        # Sort by score descending
        scored_results.sort(key=lambda x: -x[0])

        return [kw for _, kw in scored_results[:max_results]]

    @staticmethod
    def _edit_distance(s1: str, s2: str) -> int:
        """Compute Levenshtein edit distance between two strings."""
        if len(s1) < len(s2):
            return KeywordDB._edit_distance(s2, s1)

        if len(s2) == 0:
            return len(s1)

        previous_row = list(range(len(s2) + 1))
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row

        return previous_row[-1]

    # -------------------------------------------------------------------------
    # &wt and File Redirects
    # -------------------------------------------------------------------------

    def wt_types(self) -> list[WtType]:
        """All &wt TYPE values."""
        return list(self._wt_types_list)

    def wt_type(self, name: str) -> WtType | None:
        """Look up a specific &wt TYPE by name (case-insensitive)."""
        return self._wt_types_by_name.get(name.lower())

    def wt_parameters(self) -> list[WtParameter]:
        """All &wt shared parameters (ISTEP1, ISTEP2, etc.)."""
        return list(self._wt_params)

    def file_redirects(self) -> list[FileRedirect]:
        """All file redirect definitions."""
        return list(self._file_redirects_list)

    def file_redirect(self, name: str) -> FileRedirect | None:
        """Look up a specific file redirect by name (case-insensitive)."""
        return self._file_redirects_by_name.get(name.lower())

    # -------------------------------------------------------------------------
    # Utility
    # -------------------------------------------------------------------------

    def __len__(self) -> int:
        """Total number of keywords in the database."""
        return len(self._keywords)

    def __repr__(self) -> str:
        return (
            f"KeywordDB({len(self._keywords)} keywords, "
            f"{len(self._by_section)} sections, "
            f"{len(self._wt_types_list)} wt_types, "
            f"{len(self._file_redirects_list)} file_redirects)"
        )


# Module-level singleton
_keyword_db: KeywordDB | None = None


def get_keyword_db() -> KeywordDB:
    """
    Get the singleton KeywordDB instance.

    Returns:
        The shared KeywordDB instance with all indexes built.
    """
    global _keyword_db
    if _keyword_db is None:
        _keyword_db = KeywordDB()
    return _keyword_db
