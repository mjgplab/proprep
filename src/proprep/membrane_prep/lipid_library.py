"""
Lipid Library

Parses packmol-memgen's memgen.parm file from $AMBERHOME to provide lipid
search, browse, and category lookup. All lipid data is read at runtime —
no static copy is bundled.
"""

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class LipidEntry:
    """A single lipid from the memgen.parm database."""
    name: str
    full_name: str
    charge: int
    apl_ff: Optional[float]  # Area per lipid (force field)
    apl_exp: Optional[float]  # Area per lipid (experimental)
    volume: Optional[float]
    charmm_available: bool
    is_extended: bool = False  # Requires leaprc.lipid_ext (after #EXTENDED LIPIDS in memgen.parm)
    category: str = ""  # Assigned during categorization

    @property
    def charge_str(self) -> str:
        if self.charge == 0:
            return "0"
        elif self.charge > 0:
            return f"+{self.charge}"
        else:
            return str(self.charge)

    @property
    def apl_display(self) -> str:
        if self.apl_ff is not None:
            return f"{self.apl_ff:.0f}"
        return "N/A"


# Lipid categories keyed on packmol-memgen's naming convention (head-group
# suffix or documented prefix). Only names present in memgen.parm are listed
# explicitly; nothing about biological occurrence is asserted here, and any
# category that matches no database entry is not shown.
LIPID_CATEGORIES = {
    "Phosphatidylcholines (PC)": {
        "pattern": r"^(?:si)?[A-Z]{2}PC$|^POPC$|^SOPC$|^OSPC$|^OPPC$|^YOPC$|^LYPC$",
        "suffixes": ["PC"],
        "names": [],
    },
    "Phosphatidylethanolamines (PE)": {
        "pattern": r"^(?:si)?[A-Z]{2}PE$|^POPE$|^SOPE$",
        "suffixes": ["PE"],
        "names": [],
    },
    "Phosphatidylserines (PS)": {
        "pattern": r"^[A-Z]{2}PS$|^POPS$|^SOPS$",
        "suffixes": ["PS"],
        "names": [],
    },
    "Phosphatidylglycerols (PG)": {
        "pattern": r"^[A-Z]{2}PG$|^POPG$|^SOPG$",
        "suffixes": ["PG"],
        "names": [],
    },
    "Phosphatidic acids (PA)": {
        "pattern": r"^[A-Z]{2}PA$|^POPA$|^SOPA$",
        "suffixes": ["PA"],
        "names": [],
    },
    "Phosphatidylinositols (PI)": {
        "pattern": r"^[A-Z]{2}PI$|^SAPI$|^SLPI$|^PIP[123]?$|^PI[34][5P]$",
        "suffixes": ["PI"],
        "names": ["SAPI", "SLPI"],
    },
    "Sphingomyelins (SM)": {
        "pattern": r"^[A-Z]?SM$|^PSM$|^SSM$|^LSM$",
        "suffixes": ["SM"],
        "names": ["PSM", "SSM", "LSM"],
    },
    "Sterols": {
        "pattern": r"^CHL1$|^ERG$|^SIT$|^STI$|^CAM$",
        "suffixes": [],
        "names": ["CHL1", "ERG", "SIT", "STI", "CAM"],
    },
    "Cardiolipins (CL)": {
        "pattern": r"^[A-Z]{2}CL2?$|^T[A-Z]CL2?$|^AR[0-9]*CL$",
        "suffixes": ["CL", "CL2"],
        "names": [],
    },
    "SIRAH coarse-grain (si prefix)": {
        "pattern": r"^si[A-Z]",
        "suffixes": [],
        "names": [],
    },
}


class LipidLibrary:
    """
    Reads and provides access to the packmol-memgen lipid database.

    Data is loaded lazily from $AMBERHOME on first access.
    """

    def __init__(self):
        self._entries: Dict[str, LipidEntry] = {}
        self._loaded = False
        self._load_error: Optional[str] = None

    @property
    def is_loaded(self) -> bool:
        return self._loaded

    @property
    def load_error(self) -> Optional[str]:
        return self._load_error

    def _find_memgen_parm(self) -> Optional[Path]:
        """Locate memgen.parm file."""
        candidates = []

        # Check if packmol_memgen is installed as a Python package (highest priority)
        try:
            import packmol_memgen
            pkg_dir = Path(packmol_memgen.__file__).parent
            candidates.append(pkg_dir / "data" / "memgen.parm")
        except ImportError:
            pass

        # Check AMBERHOME-based paths
        amberhome = os.environ.get("AMBERHOME")
        if amberhome:
            candidates.append(
                Path(amberhome) / "lib" / "python" / "packmol_memgen" / "data" / "memgen.parm"
            )
            candidates.append(
                Path(amberhome) / "AmberTools" / "src" / "packmol_memgen" / "packmol_memgen" / "data" / "memgen.parm"
            )

        for path in candidates:
            if path.exists():
                return path

        return None

    def load(self) -> bool:
        """
        Load lipid data from memgen.parm.

        Returns:
            True if loaded successfully, False otherwise.
        """
        if self._loaded:
            return True

        parm_path = self._find_memgen_parm()
        if parm_path is None:
            self._load_error = (
                "Could not find memgen.parm. Ensure $AMBERHOME is set and "
                "packmol-memgen is installed."
            )
            logger.warning(self._load_error)
            return False

        try:
            self._parse_memgen_parm(parm_path)
            self._categorize_lipids()
            self._loaded = True
            logger.debug(f"Loaded {len(self._entries)} lipids from {parm_path}")
            return True
        except Exception as e:
            self._load_error = f"Error parsing memgen.parm: {e}"
            logger.error(self._load_error)
            return False

    def _parse_memgen_parm(self, path: Path):
        """Parse the memgen.parm file."""
        is_extended = False
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    if line.startswith("#EXTENDED LIPIDS"):
                        is_extended = True
                    continue

                parts = line.split()
                if len(parts) < 12:
                    continue

                name = parts[0]

                # Parse APL (force field)
                try:
                    apl_ff = float(parts[3])
                except (ValueError, IndexError):
                    apl_ff = None

                # Parse APL (experimental)
                try:
                    apl_exp = float(parts[4]) if parts[4] != "XX" else None
                except (ValueError, IndexError):
                    apl_exp = None

                # Parse volume
                try:
                    volume = float(parts[5]) if parts[5] != "XXXX" else None
                except (ValueError, IndexError):
                    volume = None

                # Parse charge
                try:
                    charge = int(parts[6])
                except (ValueError, IndexError):
                    charge = 0

                # CHARMM availability
                charmm = parts[11].upper() == "Y" if len(parts) > 11 else False

                # Full name — everything from column 12 onward
                full_name = " ".join(parts[12:]) if len(parts) > 12 else name

                self._entries[name] = LipidEntry(
                    name=name,
                    full_name=full_name,
                    charge=charge,
                    apl_ff=apl_ff,
                    apl_exp=apl_exp,
                    volume=volume,
                    charmm_available=charmm,
                    is_extended=is_extended,
                )

    def _categorize_lipids(self):
        """Assign categories to lipids based on naming patterns."""
        for entry in self._entries.values():
            entry.category = self._determine_category(entry.name)

    def _determine_category(self, name: str) -> str:
        """Determine the category of a lipid by name and full name."""
        entry = self._entries.get(name)
        full = entry.full_name.lower() if entry else ""

        for category, info in LIPID_CATEGORIES.items():
            # Check explicit names first
            if name in info.get("names", []):
                return category
            # Check regex pattern on short name
            if re.match(info["pattern"], name):
                return category
            # Check suffixes on short name
            for suffix in info.get("suffixes", []):
                if name.endswith(suffix) and len(name) <= 5:
                    return category

        # Fall back to categorizing by full name keywords
        if full:
            if "phosphocholine" in full or "phosphatidylcholine" in full:
                return "Phosphatidylcholines (PC)"
            if "phosphoethanolamine" in full or "phopshoethanolamine" in full:
                return "Phosphatidylethanolamines (PE)"
            if "phosphoserine" in full or "phosphatidylserine" in full:
                return "Phosphatidylserines (PS)"
            if "phosphoglycerol" in full or "phosphatidylglycerol" in full:
                return "Phosphatidylglycerols (PG)"
            if "phosphat" in full and ("acid" in full or "phosphate" in full or "-phospho" not in full):
                if "inositol" in full:
                    return "Phosphatidylinositols (PI)"
                return "Phosphatidic acids (PA)"
            if "inositol" in full:
                return "Phosphatidylinositols (PI)"
            if "sphingomyelin" in full:
                return "Sphingomyelins (SM)"
            if "ceramid" in full:
                return "Ceramides (CER)"
            if "cholest" in full or "ergost" in full or "stigmast" in full or "sitost" in full or "campest" in full:
                return "Sterols"
            if "cardiolipin" in full or name.endswith("CL") or name.endswith("CL2"):
                return "Cardiolipins (CL)"
            if "galactosyl" in full or ("sulfo" in full and "lipid" in full):
                return "Glycolipids"

        return "Other"

    def get_all(self) -> List[LipidEntry]:
        """Get all lipid entries."""
        self.load()
        return sorted(self._entries.values(), key=lambda e: e.name)

    def get(self, name: str) -> Optional[LipidEntry]:
        """Get a specific lipid by name."""
        self.load()
        return self._entries.get(name)

    def search(self, query: str) -> List[LipidEntry]:
        """
        Search lipids by name or full name (case-insensitive).

        Args:
            query: Search string to match against name or full_name.

        Returns:
            List of matching LipidEntry objects.
        """
        self.load()
        query_lower = query.lower()
        results = []
        for entry in self._entries.values():
            if (query_lower in entry.name.lower() or
                    query_lower in entry.full_name.lower()):
                results.append(entry)
        return sorted(results, key=lambda e: e.name)

    def get_by_category(self, category: str) -> List[LipidEntry]:
        """Get all lipids in a category."""
        self.load()
        return sorted(
            [e for e in self._entries.values() if e.category == category],
            key=lambda e: e.name,
        )

    def get_categories(self) -> List[Tuple[str, str, int]]:
        """
        Get list of categories with a charge summary and counts.

        The summary is computed from the database entries in the category
        (the charge column of memgen.parm), so it states only what the
        database says.

        Returns:
            List of (category_name, charge_summary, count) tuples.
        """
        self.load()
        result = []
        # Categories assigned from the database's own full names (the keyword
        # fallback in _categorize) may not be in LIPID_CATEGORIES; list them
        # too so no entry is hidden from the browser.
        assigned = {e.category for e in self._entries.values()}
        extra = sorted(assigned - set(LIPID_CATEGORIES) - {"Other"})
        names = list(LIPID_CATEGORIES) + extra + ["Other"]
        for cat_name in names:
            members = [e for e in self._entries.values() if e.category == cat_name]
            if members:
                result.append((cat_name, self._charge_summary(members), len(members)))
        return result

    @staticmethod
    def _charge_summary(entries: List["LipidEntry"]) -> str:
        charges = sorted({e.charge for e in entries})
        if len(charges) == 1:
            return f"charge {charges[0]:+d}" if charges[0] else "charge 0"
        return f"charges {charges[0]:+d} to {charges[-1]:+d}"

    def validate_lipid_string(self, lipid_str: str) -> Tuple[bool, str]:
        """
        Validate a lipid composition string.

        Checks that all lipid names exist in the database.

        Args:
            lipid_str: e.g. "DOPC:CHL1" or "DOPC//DOPE"

        Returns:
            Tuple of (is_valid, error_message).
        """
        self.load()

        # Split by leaflet separators
        leaflets = re.split(r"//|///", lipid_str)

        for leaflet in leaflets:
            lipid_names = leaflet.split(":")
            for name in lipid_names:
                name = name.strip()
                if name and name != "None" and name not in self._entries:
                    return False, f"Unknown lipid: '{name}'"

        return True, ""

    def needs_extended(self, lipid_str: str) -> bool:
        """
        Check whether a lipid composition requires leaprc.lipid_ext.

        Mirrors packmol-memgen's logic: any lipid after the #EXTENDED LIPIDS
        marker in memgen.parm triggers the extended force field requirement.

        Args:
            lipid_str: Lipid composition string (e.g., "DOPC:CHL1" or "DOPC//DOPE")

        Returns:
            True if any lipid in the composition is an extended lipid.
        """
        self.load()
        leaflets = re.split(r"//|///", lipid_str)
        for leaflet in leaflets:
            for name in leaflet.split(":"):
                name = name.strip()
                entry = self._entries.get(name)
                if entry and entry.is_extended:
                    return True
        return False

    def get_lipid_ext_addpath(self) -> Optional[str]:
        """
        Get the addPath directory needed for leaprc.lipid_ext.

        Returns the packmol_memgen data directory path, or None if not found.
        """
        try:
            import packmol_memgen
            return str(Path(packmol_memgen.__file__).parent / "data")
        except ImportError:
            return None

    def get_total_charge(self, lipid_str: str, ratio_str: str = "") -> Optional[int]:
        """
        Estimate the charge contribution from lipids (per leaflet pair).

        This is a rough estimate — actual charge depends on lipid counts
        calculated by packmol-memgen from geometry.
        """
        self.load()

        # Just report charge per lipid molecule for each component
        # Actual ion counts depend on the number of lipids placed
        return None  # Calculation requires geometry; leave to packmol-memgen
