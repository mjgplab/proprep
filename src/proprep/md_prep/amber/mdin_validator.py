"""
AMBER MDIN Configuration Interface - Layer 2: Validation Engine

Validates MDIN files against the keyword database, checking for:
- Lexical errors (unknown keywords, wrong sections, duplicates)
- Type and range violations
- Dependency and consistency issues
- Best practice recommendations

Usage:
    from proprep.md_prep.amber import MdinValidator, validate_mdin

    # Validate text
    result = validate_mdin(mdin_text)
    if not result.is_valid:
        for error in result.errors:
            print(f"{error.code}: {error.message}")

    # Validate file
    validator = MdinValidator()
    result = validator.validate_file("production.mdin")
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable

from .keyword_db import KeywordDB, get_keyword_db
from .mdin_parser import ParsedMdin, parse_mdin, parse_mdin_file

logger = logging.getLogger(__name__)


# =============================================================================
# Validation Result Model
# =============================================================================


class Severity(Enum):
    """Severity levels for validation issues."""

    ERROR = "error"  # Will definitely cause AMBER to fail
    WARNING = "warning"  # Likely a mistake or suboptimal setting
    INFO = "info"  # Suggestion or informational note
    STYLE = "style"  # Non-functional formatting suggestion


@dataclass
class ValidationIssue:
    """A single validation issue found in an MDIN file."""

    severity: Severity
    code: str  # Machine-readable code, e.g. "UNKNOWN_KEYWORD"
    message: str  # Human-readable description
    keyword: str | None = None
    section: str | None = None
    line_number: int | None = None
    suggestion: str | None = None  # "Did you mean 'nstlim'?"
    related_keywords: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        parts = []
        if self.line_number:
            parts.append(f"Line {self.line_number}")
        parts.append(f"{self.severity.value.upper()}")
        parts.append(f"[{self.code}]")
        parts.append(self.message)
        result = " ".join(parts)
        if self.suggestion:
            result += f"\n    {self.suggestion}"
        return result


@dataclass
class ValidationResult:
    """Complete validation result for an MDIN file."""

    issues: list[ValidationIssue] = field(default_factory=list)
    parsed: ParsedMdin | None = None

    @property
    def errors(self) -> list[ValidationIssue]:
        """All ERROR severity issues."""
        return [i for i in self.issues if i.severity == Severity.ERROR]

    @property
    def warnings(self) -> list[ValidationIssue]:
        """All WARNING severity issues."""
        return [i for i in self.issues if i.severity == Severity.WARNING]

    @property
    def infos(self) -> list[ValidationIssue]:
        """All INFO severity issues."""
        return [i for i in self.issues if i.severity == Severity.INFO]

    @property
    def is_valid(self) -> bool:
        """True if no ERROR severity issues."""
        return len(self.errors) == 0

    def __str__(self) -> str:
        lines = []
        for issue in self.issues:
            lines.append(str(issue))
        lines.append("")
        lines.append(
            f"{len(self.errors)} error(s), {len(self.warnings)} warning(s), "
            f"{len(self.infos)} info(s)"
        )
        return "\n".join(lines)


# =============================================================================
# Dependency Rule Type
# =============================================================================

# A dependency rule is a callable that receives the parsed MDIN and KeywordDB,
# and returns a list of validation issues.
DependencyRule = Callable[[ParsedMdin, KeywordDB], list[ValidationIssue]]


# =============================================================================
# Lexical Rules (Category 1)
# =============================================================================


def _check_unknown_keywords(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check for keywords not found in the database."""
    issues: list[ValidationIssue] = []

    for section_name, block in parsed.blocks.items():
        # Get valid keywords for this section
        valid_keywords = {kw.name.lower() for kw in db.by_section(section_name)}

        for key, entry in block.entries.items():
            key_lower = key.lower()
            if key_lower not in valid_keywords:
                # Check if it exists in another section (WRONG_SECTION)
                all_matches = db.get_all(key_lower)
                if all_matches:
                    correct_sections = [kw.section for kw in all_matches]
                    issues.append(
                        ValidationIssue(
                            severity=Severity.ERROR,
                            code="WRONG_SECTION",
                            message=(
                                f"Keyword '{entry.raw_key}' exists but belongs in "
                                f"&{correct_sections[0]}, not &{section_name}."
                            ),
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                            suggestion=f"Move to &{correct_sections[0]} section.",
                            related_keywords=correct_sections,
                        )
                    )
                else:
                    # Truly unknown - try fuzzy match
                    similar = db.search(key_lower, max_results=3)
                    suggestion = None
                    if similar:
                        # Filter to same section if possible
                        same_section = [
                            kw for kw in similar if kw.section.lower() == section_name
                        ]
                        if same_section:
                            suggestion = f"Did you mean '{same_section[0].name}'?"
                        elif similar:
                            suggestion = f"Did you mean '{similar[0].name}' (in &{similar[0].section})?"

                    issues.append(
                        ValidationIssue(
                            severity=Severity.WARNING,
                            code="UNKNOWN_KEYWORD",
                            message=f"Unknown keyword '{entry.raw_key}' in &{section_name}.",
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                            suggestion=suggestion,
                        )
                    )

    return issues


def _check_unknown_sections(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check for unrecognized namelist sections."""
    issues: list[ValidationIssue] = []
    valid_sections = set(db.sections())

    for section_name, block in parsed.blocks.items():
        if section_name not in valid_sections:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="UNKNOWN_SECTION",
                    message=f"Unrecognized namelist section '&{block.raw_section}'.",
                    section=section_name,
                    line_number=block.line_start,
                )
            )

    return issues


# =============================================================================
# Type & Range Rules (Category 2)
# =============================================================================


def _check_type_and_range(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check value types and ranges against keyword definitions."""
    issues: list[ValidationIssue] = []

    for section_name, block in parsed.blocks.items():
        for key, entry in block.entries.items():
            kw = db.get(key, section_name)
            if kw is None:
                continue  # Already handled by unknown keyword check

            value = entry.value

            # Type check
            expected_type = kw.value_type
            type_ok = True

            if expected_type == "int":
                if not isinstance(value, int) or isinstance(value, bool):
                    type_ok = False
                    issues.append(
                        ValidationIssue(
                            severity=Severity.ERROR,
                            code="TYPE_MISMATCH",
                            message=(
                                f"'{key}' expects an integer, got "
                                f"{type(value).__name__} ({value})."
                            ),
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                        )
                    )
            elif expected_type == "float":
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    type_ok = False
                    issues.append(
                        ValidationIssue(
                            severity=Severity.ERROR,
                            code="TYPE_MISMATCH",
                            message=(
                                f"'{key}' expects a number, got "
                                f"{type(value).__name__} ({value})."
                            ),
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                        )
                    )
            elif expected_type == "str":
                # Strings are flexible - most values can be treated as strings
                pass

            # Range check (only if type is ok)
            if type_ok and isinstance(value, (int, float)) and not isinstance(value, bool):
                if kw.min_val is not None and value < kw.min_val:
                    issues.append(
                        ValidationIssue(
                            severity=Severity.ERROR,
                            code="OUT_OF_RANGE",
                            message=(
                                f"'{key}' value {value} is below minimum {kw.min_val}."
                            ),
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                        )
                    )
                if kw.max_val is not None and value > kw.max_val:
                    issues.append(
                        ValidationIssue(
                            severity=Severity.ERROR,
                            code="OUT_OF_RANGE",
                            message=(
                                f"'{key}' value {value} is above maximum {kw.max_val}."
                            ),
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                        )
                    )

            # Option check (only for keywords with exhaustive options)
            if kw.options and type_ok:
                if value not in kw.options:
                    valid_opts = list(kw.options.keys())
                    issues.append(
                        ValidationIssue(
                            severity=Severity.WARNING,
                            code="INVALID_OPTION",
                            message=(
                                f"'{key}' value {value} is not a recognized option. "
                                f"Valid options: {valid_opts}"
                            ),
                            keyword=key,
                            section=section_name,
                            line_number=entry.line_number,
                        )
                    )

    return issues


# =============================================================================
# Dependency Rules (Category 3)
# =============================================================================


def _get_value(parsed: ParsedMdin, key: str, default: Any = None) -> Any:
    """Helper to get a value with default."""
    val = parsed.get_value(key)
    return val if val is not None else default


def _rule_md_vs_minimization(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check MD vs minimization consistency."""
    issues: list[ValidationIssue] = []
    imin = _get_value(parsed, "imin", 0)

    if imin == 0:
        # MD run - nstlim and dt should be set
        if parsed.get_value("nstlim") is None:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="MISSING_DEPENDENCY",
                    message="imin=0 (MD) but nstlim is not set. Defaults to 1 step.",
                    keyword="nstlim",
                    related_keywords=["imin", "nstlim", "dt"],
                )
            )
        if parsed.get_value("dt") is None:
            issues.append(
                ValidationIssue(
                    severity=Severity.INFO,
                    code="MISSING_DEPENDENCY",
                    message="imin=0 (MD) but dt is not set. Defaults to 0.001 ps.",
                    keyword="dt",
                    related_keywords=["imin", "dt"],
                )
            )
    elif imin == 1:
        # Minimization - maxcyc should be set
        if parsed.get_value("maxcyc") is None:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="MISSING_DEPENDENCY",
                    message="imin=1 (minimization) but maxcyc is not set. Defaults to 1 cycle.",
                    keyword="maxcyc",
                    related_keywords=["imin", "maxcyc", "ncyc"],
                )
            )

    return issues


def _rule_restart_consistency(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check restart flag consistency."""
    issues: list[ValidationIssue] = []
    irest = _get_value(parsed, "irest", 0)
    ntx = _get_value(parsed, "ntx", 1)

    if irest == 1 and ntx != 5:
        issues.append(
            ValidationIssue(
                severity=Severity.ERROR,
                code="INCONSISTENT_SETTING",
                message=(
                    "irest=1 (restart) requires ntx=5 to read velocities. "
                    f"Currently ntx={ntx}."
                ),
                keyword="ntx",
                suggestion="Set ntx=5 when irest=1.",
                related_keywords=["irest", "ntx"],
            )
        )

    if irest == 0 and ntx == 5:
        issues.append(
            ValidationIssue(
                severity=Severity.INFO,
                code="INCONSISTENT_SETTING",
                message=(
                    "ntx=5 reads velocities but irest=0 will ignore them and "
                    "assign new velocities."
                ),
                keyword="irest",
                related_keywords=["irest", "ntx"],
            )
        )

    return issues


def _rule_thermostat_completeness(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check thermostat settings are complete."""
    issues: list[ValidationIssue] = []
    ntt = _get_value(parsed, "ntt", 0)
    imin = _get_value(parsed, "imin", 0)

    # Only check for MD runs
    if imin != 0:
        return issues

    if ntt == 1:
        # Berendsen - tautp should be considered
        tautp = _get_value(parsed, "tautp", 1.0)
        if tautp < 0.5:
            issues.append(
                ValidationIssue(
                    severity=Severity.INFO,
                    code="TIGHT_COUPLING",
                    message=f"tautp={tautp} is quite tight. Values 0.5-5.0 ps are typical.",
                    keyword="tautp",
                    related_keywords=["ntt", "tautp"],
                )
            )
    elif ntt == 2:
        # Andersen - ig should be set
        ig = _get_value(parsed, "ig", -1)
        if ig != -1:
            issues.append(
                ValidationIssue(
                    severity=Severity.INFO,
                    code="IG_NOT_RANDOM",
                    message=(
                        "ntt=2 (Andersen) with fixed ig. For independent trajectories, "
                        "consider ig=-1 for random seed."
                    ),
                    keyword="ig",
                    related_keywords=["ntt", "ig"],
                )
            )
    elif ntt == 3:
        # Langevin - gamma_ln must be set
        gamma_ln = _get_value(parsed, "gamma_ln", 0.0)
        if gamma_ln == 0.0:
            issues.append(
                ValidationIssue(
                    severity=Severity.ERROR,
                    code="MISSING_DEPENDENCY",
                    message=(
                        "ntt=3 (Langevin) requires gamma_ln > 0. "
                        "Currently defaults to 0.0 (no friction)."
                    ),
                    keyword="gamma_ln",
                    suggestion="Set gamma_ln to 1.0-5.0 ps^-1 for typical simulations.",
                    related_keywords=["ntt", "gamma_ln", "temp0"],
                )
            )

        # Check ig for Langevin
        ig = _get_value(parsed, "ig", -1)
        if ig != -1:
            issues.append(
                ValidationIssue(
                    severity=Severity.INFO,
                    code="IG_NOT_RANDOM",
                    message=(
                        "ntt=3 (Langevin) with fixed ig. For independent trajectories "
                        "and to avoid synchronization artifacts, use ig=-1."
                    ),
                    keyword="ig",
                    related_keywords=["ntt", "ig"],
                )
            )

    return issues


def _rule_barostat_requirements(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check barostat settings."""
    issues: list[ValidationIssue] = []
    ntp = _get_value(parsed, "ntp", 0)
    imin = _get_value(parsed, "imin", 0)

    # Only check for MD runs with pressure control
    if imin != 0 or ntp == 0:
        return issues

    ntb = _get_value(parsed, "ntb", 1)
    if ntb != 2:
        issues.append(
            ValidationIssue(
                severity=Severity.ERROR,
                code="INCONSISTENT_SETTING",
                message=f"ntp={ntp} (constant pressure) requires ntb=2. Currently ntb={ntb}.",
                keyword="ntb",
                suggestion="Set ntb=2 for constant pressure simulations.",
                related_keywords=["ntp", "ntb"],
            )
        )

    barostat = _get_value(parsed, "barostat", 1)
    if barostat == 1:
        # Berendsen
        taup = _get_value(parsed, "taup", 1.0)
        if taup < 1.0:
            issues.append(
                ValidationIssue(
                    severity=Severity.INFO,
                    code="TIGHT_COUPLING",
                    message=f"taup={taup} is quite tight. Values 1.0-5.0 ps are recommended.",
                    keyword="taup",
                    related_keywords=["ntp", "barostat", "taup"],
                )
            )
    elif barostat == 2:
        # Monte Carlo barostat - recommended for GPU
        pass  # No additional requirements

    return issues


def _rule_shake_force_eval(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check SHAKE and force evaluation consistency."""
    issues: list[ValidationIssue] = []
    ntc = _get_value(parsed, "ntc", 1)
    ntf = _get_value(parsed, "ntf", 1)
    dt = _get_value(parsed, "dt", 0.001)
    imin = _get_value(parsed, "imin", 0)

    # ntc and ntf should typically match
    if ntc != ntf and ntc > 1:
        issues.append(
            ValidationIssue(
                severity=Severity.WARNING,
                code="INCONSISTENT_SETTING",
                message=(
                    f"ntc={ntc} but ntf={ntf}. When using SHAKE (ntc>1), "
                    "ntf should typically equal ntc."
                ),
                keyword="ntf",
                suggestion=f"Set ntf={ntc} to match ntc.",
                related_keywords=["ntc", "ntf"],
            )
        )

    # Time step checks (only for MD)
    if imin == 0:
        if ntc == 1 and dt > 0.001:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="LARGE_TIMESTEP",
                    message=(
                        f"dt={dt} ps without SHAKE (ntc=1). "
                        "Without SHAKE, dt should be <= 0.001 ps."
                    ),
                    keyword="dt",
                    suggestion="Either reduce dt to 0.001 or enable SHAKE with ntc=2.",
                    related_keywords=["ntc", "dt"],
                )
            )
        elif ntc >= 2 and dt > 0.002:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="LARGE_TIMESTEP",
                    message=(
                        f"dt={dt} ps with SHAKE. Maximum recommended is 0.002 ps "
                        "(or 0.004 ps with hydrogen mass repartitioning)."
                    ),
                    keyword="dt",
                    related_keywords=["ntc", "dt"],
                )
            )

    return issues


def _rule_periodic_boundary(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check periodic boundary condition consistency."""
    issues: list[ValidationIssue] = []
    ntb = _get_value(parsed, "ntb", 1)
    igb = _get_value(parsed, "igb", 0)
    ntp = _get_value(parsed, "ntp", 0)

    # GB requires non-periodic
    if ntb > 0 and igb > 0 and igb != 10:
        issues.append(
            ValidationIssue(
                severity=Severity.ERROR,
                code="INCONSISTENT_SETTING",
                message=(
                    f"Generalized Born (igb={igb}) requires non-periodic (ntb=0). "
                    f"Currently ntb={ntb}."
                ),
                keyword="ntb",
                suggestion="Set ntb=0 for GB simulations.",
                related_keywords=["ntb", "igb"],
            )
        )

    # Pressure requires periodicity
    if ntb == 0 and ntp > 0:
        issues.append(
            ValidationIssue(
                severity=Severity.ERROR,
                code="INCONSISTENT_SETTING",
                message=(
                    f"Constant pressure (ntp={ntp}) requires periodic boundaries (ntb>0). "
                    f"Currently ntb={ntb}."
                ),
                keyword="ntb",
                suggestion="Set ntb=2 for constant pressure.",
                related_keywords=["ntb", "ntp"],
            )
        )

    return issues


def _rule_gb_settings(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check Generalized Born settings."""
    issues: list[ValidationIssue] = []
    igb = _get_value(parsed, "igb", 0)

    if igb == 0:
        return issues  # GB not in use

    cut = _get_value(parsed, "cut", 9999.0 if igb > 0 else 8.0)
    if cut < 25.0 and cut != 9999.0:
        issues.append(
            ValidationIssue(
                severity=Severity.WARNING,
                code="SHORT_CUTOFF",
                message=(
                    f"GB (igb={igb}) with cut={cut} Å. Short cutoffs can cause "
                    "artifacts in GB. Consider cut >= 25 Å or using no cutoff."
                ),
                keyword="cut",
                related_keywords=["igb", "cut"],
            )
        )

    return issues


def _rule_constant_ph(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check constant pH settings."""
    issues: list[ValidationIssue] = []
    icnstph = _get_value(parsed, "icnstph", 0)

    if icnstph == 0:
        return issues

    # solvph must be set
    if parsed.get_value("solvph") is None:
        issues.append(
            ValidationIssue(
                severity=Severity.ERROR,
                code="MISSING_DEPENDENCY",
                message="icnstph > 0 requires solvph to be set.",
                keyword="solvph",
                related_keywords=["icnstph", "solvph", "ntcnstph"],
            )
        )

    # ntcnstph should be set
    if parsed.get_value("ntcnstph") is None:
        issues.append(
            ValidationIssue(
                severity=Severity.WARNING,
                code="MISSING_DEPENDENCY",
                message=(
                    "icnstph > 0 but ntcnstph not set. "
                    "Defaults to 10 steps between MC attempts."
                ),
                keyword="ntcnstph",
                related_keywords=["icnstph", "ntcnstph"],
            )
        )

    # Explicit solvent CpHMD needs ntrelax
    if icnstph == 2:
        if parsed.get_value("ntrelax") is None:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="MISSING_DEPENDENCY",
                    message=(
                        "icnstph=2 (explicit solvent CpHMD) but ntrelax not set. "
                        "Defaults to 100 relaxation steps."
                    ),
                    keyword="ntrelax",
                    related_keywords=["icnstph", "ntrelax"],
                )
            )

    return issues


def _rule_constant_redox(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check constant redox potential settings."""
    issues: list[ValidationIssue] = []
    icnste = _get_value(parsed, "icnste", 0)

    if icnste == 0:
        return issues

    # solve must be set
    if parsed.get_value("solve") is None:
        issues.append(
            ValidationIssue(
                severity=Severity.ERROR,
                code="MISSING_DEPENDENCY",
                message="icnste > 0 requires solve (redox potential) to be set.",
                keyword="solve",
                related_keywords=["icnste", "solve", "ntcnste"],
            )
        )

    # ntcnste should be set
    if parsed.get_value("ntcnste") is None:
        issues.append(
            ValidationIssue(
                severity=Severity.WARNING,
                code="MISSING_DEPENDENCY",
                message=(
                    "icnste > 0 but ntcnste not set. "
                    "Defaults to 10 steps between MC attempts."
                ),
                keyword="ntcnste",
                related_keywords=["icnste", "ntcnste"],
            )
        )

    # Explicit solvent needs ntrelaxe
    if icnste == 2:
        if parsed.get_value("ntrelaxe") is None:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="MISSING_DEPENDENCY",
                    message=(
                        "icnste=2 (explicit solvent) but ntrelaxe not set. "
                        "Defaults to 200 relaxation steps."
                    ),
                    keyword="ntrelaxe",
                    related_keywords=["icnste", "ntrelaxe"],
                )
            )

    return issues


def _rule_output_frequency(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check output frequency settings."""
    issues: list[ValidationIssue] = []
    imin = _get_value(parsed, "imin", 0)

    if imin != 0:
        return issues  # Only for MD

    nstlim = _get_value(parsed, "nstlim", 1)
    ntwx = _get_value(parsed, "ntwx", 0)
    ntpr = _get_value(parsed, "ntpr", 50)
    ntwr = _get_value(parsed, "ntwr", nstlim)

    # No trajectory output
    if ntwx == 0:
        issues.append(
            ValidationIssue(
                severity=Severity.INFO,
                code="NO_TRAJECTORY",
                message="ntwx=0: No trajectory file will be written.",
                keyword="ntwx",
                related_keywords=["ntwx", "nstlim"],
            )
        )

    # Energy never printed during run
    if ntpr > nstlim:
        issues.append(
            ValidationIssue(
                severity=Severity.WARNING,
                code="NO_OUTPUT",
                message=f"ntpr={ntpr} > nstlim={nstlim}: Energy will never be printed during run.",
                keyword="ntpr",
                related_keywords=["ntpr", "nstlim"],
            )
        )

    # Restart only at end
    if ntwr > nstlim and ntwr != -1:
        issues.append(
            ValidationIssue(
                severity=Severity.INFO,
                code="INFREQUENT_RESTART",
                message=(
                    f"ntwr={ntwr} > nstlim={nstlim}: Restart file only written at end. "
                    "Consider more frequent restarts for crash recovery."
                ),
                keyword="ntwr",
                related_keywords=["ntwr", "nstlim"],
            )
        )

    return issues


def _rule_restraints(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check restraint settings."""
    issues: list[ValidationIssue] = []
    ntr = _get_value(parsed, "ntr", 0)

    if ntr == 1:
        # restraintmask should be set
        if not parsed.get_value("restraintmask"):
            issues.append(
                ValidationIssue(
                    severity=Severity.ERROR,
                    code="MISSING_DEPENDENCY",
                    message="ntr=1 requires restraintmask to specify restrained atoms.",
                    keyword="restraintmask",
                    related_keywords=["ntr", "restraintmask", "restraint_wt"],
                )
            )

        # restraint_wt should be set
        restraint_wt = _get_value(parsed, "restraint_wt", 0.0)
        if restraint_wt == 0.0:
            issues.append(
                ValidationIssue(
                    severity=Severity.WARNING,
                    code="ZERO_RESTRAINT",
                    message="ntr=1 but restraint_wt=0.0: Restraints have no effect.",
                    keyword="restraint_wt",
                    related_keywords=["ntr", "restraint_wt"],
                )
            )

        # netfrc should be 0 when using restraints
        netfrc = _get_value(parsed, "netfrc", 1)
        if netfrc != 0:
            issues.append(
                ValidationIssue(
                    severity=Severity.INFO,
                    code="RESTRAINT_NETFRC",
                    message="When ntr=1, consider setting netfrc=0 in &ewald.",
                    keyword="netfrc",
                    related_keywords=["ntr", "netfrc"],
                )
            )

    return issues


# =============================================================================
# Best Practice Rules (Category 4)
# =============================================================================


def _rule_best_practices(
    parsed: ParsedMdin, db: KeywordDB
) -> list[ValidationIssue]:
    """Check best practices and common recommendations."""
    issues: list[ValidationIssue] = []
    imin = _get_value(parsed, "imin", 0)

    if imin != 0:
        return issues  # Only for MD

    # iwrap not set
    ntb = _get_value(parsed, "ntb", 1)
    if ntb > 0 and parsed.get_value("iwrap") is None:
        issues.append(
            ValidationIssue(
                severity=Severity.INFO,
                code="MISSING_IWRAP",
                message=(
                    "iwrap not set for periodic simulation. Consider iwrap=1 for "
                    "long runs to prevent coordinate overflow."
                ),
                keyword="iwrap",
                related_keywords=["iwrap", "ntb"],
            )
        )

    # ioutfm not set to NetCDF
    ioutfm = _get_value(parsed, "ioutfm", 1)
    if ioutfm == 0:
        issues.append(
            ValidationIssue(
                severity=Severity.INFO,
                code="ASCII_TRAJECTORY",
                message=(
                    "ioutfm=0 writes ASCII trajectories. NetCDF (ioutfm=1) is "
                    "recommended for smaller files and better precision."
                ),
                keyword="ioutfm",
            )
        )

    return issues


# =============================================================================
# Default Rules Collection
# =============================================================================


DEFAULT_RULES: list[DependencyRule] = [
    # Lexical (these run on raw parsed data)
    _check_unknown_keywords,
    _check_unknown_sections,
    # Type & Range
    _check_type_and_range,
    # Dependency rules
    _rule_md_vs_minimization,
    _rule_restart_consistency,
    _rule_thermostat_completeness,
    _rule_barostat_requirements,
    _rule_shake_force_eval,
    _rule_periodic_boundary,
    _rule_gb_settings,
    _rule_constant_ph,
    _rule_constant_redox,
    _rule_output_frequency,
    _rule_restraints,
    # Best practices
    _rule_best_practices,
]


# =============================================================================
# Validator Class
# =============================================================================


class MdinValidator:
    """
    Validates MDIN files against the keyword database and dependency rules.

    Usage:
        validator = MdinValidator()
        result = validator.validate_file("production.mdin")
        if not result.is_valid:
            for error in result.errors:
                print(error)
    """

    def __init__(
        self,
        db: KeywordDB | None = None,
        rules: list[DependencyRule] | None = None,
    ):
        """
        Initialize the validator.

        Args:
            db: KeywordDB instance. Defaults to singleton.
            rules: List of validation rules. Defaults to DEFAULT_RULES.
        """
        self.db = db if db is not None else get_keyword_db()
        self.rules = rules if rules is not None else list(DEFAULT_RULES)

    def validate(self, mdin: ParsedMdin) -> ValidationResult:
        """
        Validate a parsed MDIN structure.

        Args:
            mdin: ParsedMdin object to validate

        Returns:
            ValidationResult with all found issues
        """
        all_issues: list[ValidationIssue] = []

        for rule in self.rules:
            try:
                issues = rule(mdin, self.db)
                all_issues.extend(issues)
            except Exception as e:
                logger.warning(f"Rule {rule.__name__} failed: {e}")

        # Sort by severity (errors first), then by line number
        severity_order = {
            Severity.ERROR: 0,
            Severity.WARNING: 1,
            Severity.INFO: 2,
            Severity.STYLE: 3,
        }
        all_issues.sort(
            key=lambda i: (severity_order[i.severity], i.line_number or 0)
        )

        return ValidationResult(issues=all_issues, parsed=mdin)

    def validate_file(self, path: str | Path) -> ValidationResult:
        """
        Parse and validate an MDIN file.

        Args:
            path: Path to the MDIN file

        Returns:
            ValidationResult with all found issues
        """
        parsed = parse_mdin_file(path)
        return self.validate(parsed)

    def validate_text(self, text: str) -> ValidationResult:
        """
        Parse and validate MDIN text.

        Args:
            text: MDIN file content as string

        Returns:
            ValidationResult with all found issues
        """
        parsed = parse_mdin(text)
        return self.validate(parsed)

    def add_rule(self, rule: DependencyRule) -> None:
        """
        Add a custom validation rule.

        Args:
            rule: A callable that takes (ParsedMdin, KeywordDB) and returns
                  a list of ValidationIssue objects.
        """
        self.rules.append(rule)


# =============================================================================
# Convenience Function
# =============================================================================


def validate_mdin(
    text: str,
    rules: list[DependencyRule] | None = None,
) -> ValidationResult:
    """
    Parse and validate MDIN text.

    Convenience function that creates a validator and validates text.

    Args:
        text: MDIN file content
        rules: Optional custom rules list

    Returns:
        ValidationResult with all found issues
    """
    validator = MdinValidator(rules=rules)
    return validator.validate_text(text)


def validate_mdin_file(
    path: str | Path,
    rules: list[DependencyRule] | None = None,
) -> ValidationResult:
    """
    Parse and validate an MDIN file.

    Convenience function that creates a validator and validates a file.

    Args:
        path: Path to MDIN file
        rules: Optional custom rules list

    Returns:
        ValidationResult with all found issues
    """
    validator = MdinValidator(rules=rules)
    return validator.validate_file(path)
