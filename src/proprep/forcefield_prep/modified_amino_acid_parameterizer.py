import glob
import logging
import os
import queue
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from proprep.utils.prompts import prompt_with_context, confirm_with_context, int_prompt_with_context

from proprep.utils.prompts import prompt_with_context, confirm_with_context
from proprep.utils.workflow_checklist import WorkflowChecklist, WorkflowStep

# Setup logging
logger = logging.getLogger(__name__)

# Module-level console for Rich output
_console = Console()


def create_combined_frcmod(temp_frcmod, gaff_frcmod, combined_frcmod):
    """
    Create a combined frcmod file by replacing ATTN lines in temp_frcmod
    with corresponding lines from gaff_frcmod.

    Parameters:
    -----------
    temp_frcmod : str
        Path to the temporary frcmod file with ATTN warnings
    gaff_frcmod : str
        Path to the GAFF-based frcmod file
    combined_frcmod : str
        Path to the output combined frcmod file

    Returns:
    --------
    dict
        Dictionary containing the result of the combination
    """
    try:
        # Read the GAFF frcmod file and create a dictionary of parameters
        gaff_params = {}
        current_section = None

        with open(gaff_frcmod, "r") as f:
            for line in f:
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Check for section headers
                if (
                    line.startswith("MASS")
                    or line.startswith("BOND")
                    or line.startswith("ANGLE")
                    or line.startswith("DIHE")
                    or line.startswith("IMPROPER")
                    or line.startswith("NONBON")
                ):
                    current_section = line
                    continue

                # Skip comments
                if line.startswith("#"):
                    continue

                # Store the parameter line by its key (first part of the line)
                if current_section:
                    # Extract the parameter key based on section type
                    if current_section == "MASS":
                        # For MASS, key is atom type
                        key = line.split()[0].strip()
                    elif current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                        # For bonded terms, key is the atom types involved
                        parts = line.split()
                        key = parts[0].strip()
                        # Handle multi-term dihedral parameters
                        if len(parts) >= 2 and parts[1].strip().startswith("-"):
                            # Include the divider in the key for multi-term dihedrals
                            key += " " + parts[1].strip()
                    else:  # NONBON
                        # For nonbonded terms, key is atom type
                        key = line.split()[0].strip()

                    # Store the full line under this key and section
                    if current_section not in gaff_params:
                        gaff_params[current_section] = {}
                    gaff_params[current_section][key] = line

        # Process the temp frcmod file and replace ATTN warnings
        temp_lines = []
        replacements = 0
        current_section = None

        with open(temp_frcmod, "r") as f:
            for line in f:
                original_line = line
                line = line.strip()

                # Keep track of current section
                if (
                    line.startswith("MASS")
                    or line.startswith("BOND")
                    or line.startswith("ANGLE")
                    or line.startswith("DIHE")
                    or line.startswith("IMPROPER")
                    or line.startswith("NONBON")
                ):
                    current_section = line
                    temp_lines.append(original_line)
                    continue

                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    temp_lines.append(original_line)
                    continue

                # Check if line contains ATTN warning
                if "ATTN" in line and current_section:
                    # Extract the parameter key
                    if current_section == "MASS":
                        key = line.split()[0].strip()
                    elif current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                        parts = line.split()
                        key = parts[0].strip()
                        # Handle multi-term dihedral parameters
                        if len(parts) >= 2 and parts[1].strip().startswith("-"):
                            key += " " + parts[1].strip()
                    else:  # NONBON
                        key = line.split()[0].strip()

                    # Look for replacement in GAFF parameters
                    if (
                        current_section in gaff_params
                        and key in gaff_params[current_section]
                    ):
                        # Replace with GAFF parameter
                        gaff_line = gaff_params[current_section][key]
                        # Add a comment indicating the replacement
                        if not gaff_line.endswith("\n"):
                            gaff_line += "\n"
                        temp_lines.append(gaff_line)
                        replacements += 1
                        _console.print(f"  [green]✓[/green] Replaced: {line.strip()} → {gaff_line.strip()}")
                    else:
                        # Keep original line if no replacement found
                        temp_lines.append(original_line)
                        _console.print(f"  [yellow]○[/yellow] No GAFF replacement found for: {line.strip()}")
                else:
                    # Keep original line
                    temp_lines.append(original_line)

        # Write the combined frcmod file
        with open(combined_frcmod, "w") as f:
            f.writelines(temp_lines)

        return {
            "success": True,
            "message": f"Successfully created combined frcmod file with {replacements} replacements",
            "replacements": replacements,
        }

    except Exception as e:
        import traceback

        traceback.print_exc()
        return {
            "success": False,
            "error": f"Error creating combined frcmod file: {str(e)}",
        }


def _frcmod_penalty(line):
    """Parse a frcmod parameter line's penalty score / ATTN status.

    Returns (penalty, attn): penalty is a float when the line carries
    ``penalty score=`` (None otherwise); attn is True when the line is an
    ``ATTN, need revision`` placeholder (parmchk2 could not assign the term).
    """
    attn = "ATTN" in line
    penalty = None
    if "penalty score=" in line:
        try:
            tok = line.split("penalty score=")[1].strip().split()[0]
            penalty = float(tok.rstrip(")]}>"))
        except (IndexError, ValueError):
            penalty = None
    return penalty, attn


def _index_frcmod_params(frcmod_file):
    """Index a frcmod's bonded/nonbonded parameter lines with penalty metadata.

    Returns a dict:
      ``lines``  — the raw file lines (newline-terminated), in order.
      ``params`` — list of {idx, section, key, penalty, attn} for each
                   parameter line in BOND/ANGLE/DIHE/IMPROPER/NONBON, where
                   ``key`` is the fixed-width atom-type identifier
                   (extract_parameter_name), shared across FFs for the same
                   molecule since the prep's atom types are identical.
      ``by_key`` — {(section, key): [idx, ...]} (multi-term dihedrals map a key
                   to several line indices).
    """
    from proprep.forcefield_prep.small_molecule_parameterizer import extract_parameter_name

    section_keywords = {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}
    param_sections = {"BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}

    with open(frcmod_file, "r") as f:
        lines = f.readlines()

    params = []
    by_key = {}
    section = None
    for idx, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped in section_keywords:
            section = stripped
            continue
        if section not in param_sections:
            continue
        key = extract_parameter_name(raw.rstrip("\n"), section)
        penalty, attn = _frcmod_penalty(raw)
        params.append({"idx": idx, "section": section, "key": key,
                       "penalty": penalty, "attn": attn})
        by_key.setdefault((section, key), []).append(idx)

    return {"lines": lines, "params": params, "by_key": by_key}


def _group_frcmod_candidates(index):
    """Collapse an index's params into per-(section,key) severity summaries.

    Returns {(section, key): {penalty, attn}} where ``penalty`` is the worst
    (max) penalty across the key's lines and ``attn`` is True if any line is an
    ATTN placeholder. Keys with neither a penalty nor ATTN are omitted (they are
    already well parameterized and need no GAFF2 consideration).
    """
    summary = {}
    for p in index["params"]:
        sk = (p["section"], p["key"])
        cur = summary.get(sk)
        pen = p["penalty"]
        if cur is None:
            summary[sk] = {"penalty": pen, "attn": p["attn"]}
        else:
            if pen is not None:
                cur["penalty"] = pen if cur["penalty"] is None else max(cur["penalty"], pen)
            cur["attn"] = cur["attn"] or p["attn"]
    return {sk: v for sk, v in summary.items() if v["attn"] or v["penalty"] is not None}


def _match_dihedral_penalty(candidates, type_quad):
    """Find the penalty summary for a dihedral given its four atom types.

    ``candidates`` is the ``_group_frcmod_candidates`` output. ``type_quad`` is a
    tuple of four atom-type strings. Matches the DIHE key in either atom order and
    honours ``X`` wildcards in the frcmod key. Prefers a fully specific match over
    a wildcard one. Returns the matching value dict ({penalty, attn}) or None.
    """
    tq = tuple(t.strip() for t in type_quad)
    tqr = tuple(reversed(tq))

    def _matches(key_types):
        for target in (tq, tqr):
            if all(k == "X" or k == t for k, t in zip(key_types, target)):
                return True
        return False

    wildcard_hit = None
    for (section, key), v in candidates.items():
        if section != "DIHE":
            continue
        kt = tuple(t.strip() for t in key.split("-"))
        if len(kt) != 4 or not _matches(kt):
            continue
        if "X" not in kt:
            return v  # specific match wins
        wildcard_hit = wildcard_hit or v
    return wildcard_hit


def select_gaff_replacements(temp_index, gaff_index, processor=None, interactive=True):
    """Compare protein-FF penalties against GAFF2 and choose which to replace.

    Shows every penalized/ATTN parameter in the protein-FF frcmod alongside the
    GAFF2 counterpart's penalty, recommends the ones GAFF2 clearly improves
    (protein-FF ATTN → GAFF2 assigned, or strictly lower GAFF2 penalty), and
    lets the user edit the selection. Returns a set of (section, key) to replace.
    """
    temp_cand = _group_frcmod_candidates(temp_index)
    gaff_cand = {}
    # Best (min) real penalty per key in GAFF2, and whether it is still ATTN.
    for p in gaff_index["params"]:
        sk = (p["section"], p["key"])
        cur = gaff_cand.get(sk)
        if cur is None:
            gaff_cand[sk] = {"penalty": p["penalty"], "attn": p["attn"]}
        else:
            if p["penalty"] is not None:
                cur["penalty"] = (p["penalty"] if cur["penalty"] is None
                                  else min(cur["penalty"], p["penalty"]))
            cur["attn"] = cur["attn"] and p["attn"]
    gaff_by_key = gaff_index["by_key"]

    if not temp_cand:
        return set()

    # Severity for ordering: ATTN worst, then by descending penalty.
    def severity(v):
        return float("inf") if v["attn"] else (v["penalty"] or 0.0)

    ordered = sorted(temp_cand.items(), key=lambda kv: severity(kv[1]), reverse=True)

    rows = []
    recommended = []
    for (section, key), tv in ordered:
        gv = gaff_cand.get((section, key))
        have_gaff = (section, key) in gaff_by_key
        if not have_gaff or gv is None:
            gaff_txt, rec = "—", False
        elif gv["attn"] and gv["penalty"] is None:
            gaff_txt, rec = "ATTN", False
        else:
            gp = gv["penalty"]
            gaff_txt = f"{gp:.1f}" if gp is not None else "assigned"
            if tv["attn"]:
                rec = True                      # FF had no params; GAFF2 supplies them
            elif tv["penalty"] is not None and gp is not None:
                rec = gp < tv["penalty"]        # GAFF2 strictly better
            elif tv["penalty"] is not None and gp is None:
                rec = True                      # GAFF2 assigned cleanly
            else:
                rec = False
        ff_txt = "ATTN" if tv["attn"] else f"{tv['penalty']:.1f}"
        rows.append((section, key, ff_txt, gaff_txt, rec))
        if rec:
            recommended.append((section, key))

    table = Table(title="Penalty comparison: protein FF vs GAFF2", expand=False)
    table.add_column("#", style="grey50", justify="right")
    table.add_column("Section", style="blue")
    table.add_column("Parameter", style="cyan")
    table.add_column("FF penalty", style="yellow", justify="right")
    table.add_column("GAFF2 penalty", style="green", justify="right")
    table.add_column("GAFF2 better?", justify="center")
    for i, (section, key, ff_txt, gaff_txt, rec) in enumerate(rows, 1):
        table.add_row(str(i), section, key, ff_txt, gaff_txt,
                      "[green]✓[/green]" if rec else "")
    _console.print(table)

    _console.print(
        "\n[grey50]A ✓ means GAFF2 assigns the term with a lower penalty (or supplies one "
        "where the protein FF could not).\n"
        "Replacing splices the GAFF2 line into the frcmod, keyed by the same atom types.[/grey50]"
    )

    rec_indices = [i for i, r in enumerate(rows, 1) if r[4]]
    default_sel = ",".join(str(i) for i in rec_indices) if rec_indices else "none"

    _console.print(Panel(
        "[bold]Which penalized parameters should be replaced with GAFF2?[/bold]\n\n"
        "  • Numbers/ranges (e.g. [cyan]1,2,5-8[/cyan]) — replace those rows\n"
        "  • [cyan]all[/cyan]  — replace every row that has a GAFF2 parameter\n"
        "  • [cyan]rec[/cyan]  — replace the recommended (✓) rows [grey50](default)[/grey50]\n"
        "  • [cyan]none[/cyan] — keep the protein-FF parameters as-is",
        title="GAFF2 Replacement Selection",
        border_style="blue",
        expand=False,
    ))

    if not interactive:
        selection = "rec"
    else:
        selection = prompt_with_context(
            processor, "Parameters to replace with GAFF2", default="rec",
            module="Modified Amino Acid Parameterizer",
            description="Select GAFF2 replacements",
        ).strip().lower()

    if selection in ("none", ""):
        return set()
    if selection in ("rec", "recommended", "default"):
        return set(recommended)
    if selection == "all":
        return set((s, k) for (s, k, _f, _g, _r) in rows if (s, k) in gaff_by_key)

    chosen = set()
    try:
        for part in selection.replace(" ", "").split(","):
            if "-" in part:
                a, b = part.split("-")
                for i in range(int(a), int(b) + 1):
                    if 1 <= i <= len(rows):
                        s, k = rows[i - 1][0], rows[i - 1][1]
                        chosen.add((s, k))
            else:
                i = int(part)
                if 1 <= i <= len(rows):
                    s, k = rows[i - 1][0], rows[i - 1][1]
                    chosen.add((s, k))
    except ValueError:
        _console.print("[yellow]Could not parse selection; using recommended set.[/yellow]")
        return set(recommended)
    # Only keep choices GAFF2 can actually satisfy.
    return set(sk for sk in chosen if sk in gaff_by_key)


def create_penalty_aware_combined_frcmod(temp_index, gaff_index, combined_frcmod, selected_keys):
    """Write a combined frcmod, replacing the selected (section,key) parameter
    blocks in the protein-FF frcmod with the GAFF2 lines for the same key.

    Handles multi-term dihedrals: every protein-FF line for a selected key is
    replaced by every GAFF2 line for that key (1↔N and N↔1 both supported).
    """
    try:
        lines = list(temp_index["lines"])
        gaff_lines = gaff_index["lines"]
        gaff_by_key = gaff_index["by_key"]
        temp_by_key = temp_index["by_key"]

        replace_at = {}   # first temp line index of a key -> [gaff line strings]
        drop = set()      # remaining temp line indices for that key
        replaced = []
        for (section, key) in selected_keys:
            if (section, key) not in gaff_by_key or (section, key) not in temp_by_key:
                continue
            temp_idxs = temp_by_key[(section, key)]
            g_lines = [gaff_lines[i] for i in gaff_by_key[(section, key)]]
            g_lines = [gl if gl.endswith("\n") else gl + "\n" for gl in g_lines]
            replace_at[temp_idxs[0]] = g_lines
            for j in temp_idxs[1:]:
                drop.add(j)
            replaced.append((section, key))
            for gl in g_lines:
                _console.print(f"  [green]✓[/green] {section} {key}: {gl.strip()}")

        out = []
        for i, ln in enumerate(lines):
            if i in replace_at:
                out.extend(replace_at[i])
            elif i in drop:
                continue
            else:
                out.append(ln)

        with open(combined_frcmod, "w") as f:
            f.writelines(out)

        return {"success": True, "replacements": len(replaced), "replaced": replaced}
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"success": False, "error": f"Error creating combined frcmod: {e}"}


def run_parmchk2_gaff(prep_file, gaff_frcmod):
    """Generate a GAFF2 frcmod for a prep file (used to fill ATTN gaps).

    The Amber protein-FF parmchk2 run leaves ATTN warnings for the modified /
    ligand atoms it has no parameters for; GAFF2 supplies those. This is a thin
    wrapper over run_parmchk2 with the gaff2 parameter set.
    """
    return run_parmchk2(prep_file, gaff_frcmod, parm_set="gaff2")


def improve_frcmod_parameters(residue_symbol, prep_file, temp_frcmod):
    """
    Improve parameters in the temporary frcmod file by creating a GAFF-based frcmod
    and then combining them to replace ATTN warnings.

    Parameters:
    -----------
    residue_symbol : str
        Residue symbol/name
    prep_file : str
        Path to the prep file
    temp_frcmod : str
        Path to the temporary frcmod file

    Returns:
    --------
    dict
        Dictionary containing the result of the improvement process
    """
    import os

    _console.print("\n[bold cyan]Improving Force Field Parameters[/bold cyan]")

    # Check if files exist
    if not os.path.exists(prep_file):
        return {"success": False, "error": f"Prep file not found: {prep_file}"}

    if not os.path.exists(temp_frcmod):
        return {
            "success": False,
            "error": f"Temporary frcmod file not found: {temp_frcmod}",
        }

    # Generate output filenames
    base_name = os.path.splitext(os.path.basename(temp_frcmod))[0]
    residue_name = residue_symbol.lower()
    gaff_frcmod = f"{residue_name}_gaff.frcmod"
    final_frcmod = f"{residue_name}.frcmod"

    # Check for ATTN warnings in temp_frcmod
    attn_count = 0
    with open(temp_frcmod, "r") as f:
        for line in f:
            if "ATTN" in line:
                attn_count += 1

    if attn_count == 0:
        _console.print(f"[green]✓[/green] No ATTN warnings found in {temp_frcmod}. No improvements needed.")
        # Create a copy of the temp_frcmod as the final frcmod
        import shutil

        shutil.copy(temp_frcmod, final_frcmod)

        return {
            "success": True,
            "message": f"No ATTN warnings found. Copied {temp_frcmod} to {final_frcmod}.",
            "final_frcmod": final_frcmod,
        }

    _console.print(f"[yellow]⚠[/yellow] Found {attn_count} ATTN warnings in {temp_frcmod}.")

    # Generate GAFF-based frcmod file
    _console.print(f"\n[cyan]→[/cyan] Generating GAFF-based parameters for {residue_symbol}...")
    gaff_result = run_parmchk2_gaff(prep_file, gaff_frcmod)

    if not gaff_result["success"]:
        _console.print(
            f"[red]✗ Error generating GAFF parameters: {gaff_result.get('error', 'Unknown error')}[/red]"
        )
        return gaff_result

    _console.print(f"[green]✓[/green] Successfully generated GAFF parameters in {gaff_frcmod}.")

    # Create combined frcmod file
    _console.print(f"\n[cyan]→[/cyan] Combining parameters from {temp_frcmod} and {gaff_frcmod}...")
    combine_result = create_combined_frcmod(temp_frcmod, gaff_frcmod, final_frcmod)

    if not combine_result["success"]:
        _console.print(
            f"[red]✗ Error creating combined frcmod file: {combine_result.get('error', 'Unknown error')}[/red]"
        )
        return combine_result

    replacements = combine_result.get("replacements", 0)
    _console.print(f"[green]✓[/green] Successfully created combined frcmod file {final_frcmod}.")
    _console.print(
        f"  Replaced {replacements} of {attn_count} ATTN warnings with GAFF parameters."
    )

    if replacements < attn_count:
        _console.print(
            f"\n[yellow]⚠ WARNING: {attn_count - replacements} ATTN warnings still remain in the final frcmod file.[/yellow]"
        )
        _console.print(
            "[grey50]  You should manually check and fix these parameters before using the force field.[/grey50]"
        )

    return {
        "success": True,
        "message": f"Successfully improved force field parameters",
        "temp_frcmod": temp_frcmod,
        "gaff_frcmod": gaff_frcmod,
        "final_frcmod": final_frcmod,
        "attn_warnings": attn_count,
        "replacements": replacements,
    }


def find_amber_parm_files():
    """
    Find and list all parm*.dat files in the AMBERHOME/dat/leap/parm directory.

    Returns:
    --------
    dict
        A dictionary containing paths to available parameter files
    """
    import glob
    import os

    # Try to get AMBERHOME from environment variable
    amberhome = os.environ.get("AMBERHOME")

    if not amberhome:
        _console.print("[yellow]⚠[/yellow] AMBERHOME environment variable not found.")
        # Try to find it in common locations
        common_locations = [
            "/usr/local/amber",
            "/opt/amber",
            "/usr/local/amber22",
            "/opt/amber22",
            "/usr/local/amber20",
            "/opt/amber20",
            "/usr/local/amber18",
            "/opt/amber18",
            "/usr/local/AmberTools",
            "/opt/AmberTools",
        ]

        for location in common_locations:
            if os.path.exists(location):
                amberhome = location
                _console.print(f"[cyan]ℹ[/cyan] Found potential AMBERHOME at {amberhome}")
                break

    if amberhome:
        # Search for parameter files
        parm_dir = os.path.join(amberhome, "dat", "leap", "parm")
        if os.path.exists(parm_dir):
            parm_files = sorted(glob.glob(os.path.join(parm_dir, "parm*.dat")))

            if parm_files:
                result = {
                    "success": True,
                    "amberhome": amberhome,
                    "parm_dir": parm_dir,
                    "parm_files": parm_files,
                }
                return result
            else:
                return {
                    "success": False,
                    "error": f"No parameter files found in {parm_dir}",
                }
        else:
            return {
                "success": False,
                "error": f"Parameter directory not found: {parm_dir}",
            }
    else:
        return {
            "success": False,
            "error": "AMBERHOME not found. Please set the AMBERHOME environment variable.",
        }


def find_prep_files():
    """
    Find all .prep files in the current directory.

    Returns:
    --------
    list
        A list of .prep files in the current directory
    """
    import glob

    prep_files = sorted(glob.glob("*.prep"))
    return prep_files


def _parse_prep_residue(prep_file):
    """Parse an AMBER prep file's single residue.

    Returns a dict: ``atoms`` (ordered list of {seq,name,type,topo,na,charge}),
    ``bonds`` (set of frozenset({seqA,seqB}) from the Z-matrix tree + LOOP
    ring-closures), ``head`` (main-chain N atom name), ``tail`` (main-chain C
    atom name), ``by_seq``, ``by_name``. Types/charges are residuegen's
    authoritative AMBER assignments (do NOT re-perceive them — antechamber
    mis-types the cap-removed backbone N).
    """
    lines = open(prep_file).read().splitlines()
    atoms = []
    for ln in lines:
        p = ln.split()
        if len(p) >= 11 and p[0].isdigit() and p[2] != "DU":
            atoms.append({"seq": int(p[0]), "name": p[1], "type": p[2],
                          "topo": p[3], "na": int(p[4]), "charge": float(p[10])})
    by_seq = {a["seq"]: a for a in atoms}
    by_name = {a["name"]: a for a in atoms}
    seqset = set(by_seq)
    bonds = set()
    for a in atoms:
        if a["na"] in seqset:                       # Z-matrix parent bond
            bonds.add(frozenset((a["seq"], a["na"])))
    if "LOOP" in lines:                             # ring-closure bonds (by name)
        for ln in lines[lines.index("LOOP") + 1:]:
            s = ln.strip()
            if s in ("IMPROPER", "DONE", "STOP"):
                break
            q = ln.split()
            if len(q) == 2 and q[0] in by_name and q[1] in by_name:
                bonds.add(frozenset((by_name[q[0]]["seq"], by_name[q[1]]["seq"])))
    heads = [a for a in atoms if a["topo"] == "M" and a["type"] == "N"]
    tails = [a for a in atoms if a["topo"] == "M" and a["type"] == "C"]
    return {"atoms": atoms, "bonds": bonds, "by_seq": by_seq, "by_name": by_name,
            "head": heads[0]["name"] if heads else None,
            "tail": tails[-1]["name"] if tails else None}


def _derive_conjugate_atom_map(opt_gjf, capped_pdb, ac_file):
    """Recover each library atom's original PDB name + source residue.

    The QM log strips atom names (Gaussian keeps only elements + coordinates),
    so antechamber assigns fresh names and residuegen reorders atoms — there is
    no name- or index-preserved link back to the input structure. But the
    optimization gjf carries the crystal coordinates in the pipeline's atom
    order (copied verbatim from the capped PDB), and the ``.ac`` has that same
    order, so matching gjf ↔ capped-PDB by coordinate recovers, for each library
    (antechamber) atom name, its original name and source residue.

    Returns ``{antechamber_name: {"orig": original_name, "resname": source}}``.
    """
    def _isf(x):
        try:
            float(x); return True
        except ValueError:
            return False

    gjf = []
    for ln in open(opt_gjf):
        p = ln.split()
        if len(p) == 5 and p[0][0].isalpha() and p[1] in ("0", "-1") and all(_isf(v) for v in p[2:]):
            gjf.append((float(p[2]), float(p[3]), float(p[4])))
    cap = []
    for ln in open(capped_pdb):
        if ln.startswith(("ATOM", "HETATM")):
            cap.append((ln[12:16].strip(), ln[17:20].strip(),
                        float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    acn = []
    for ln in open(ac_file):
        if ln.startswith(("ATOM", "HETATM")):
            acn.append(ln[12:16].strip())

    amap, used = {}, set()
    for i, g in enumerate(gjf):
        for j, c in enumerate(cap):
            if j in used:
                continue
            if abs(g[0] - c[2]) < 1e-2 and abs(g[1] - c[3]) < 1e-2 and abs(g[2] - c[4]) < 1e-2:
                used.add(j)
                if i < len(acn):
                    amap[acn[i]] = {"orig": c[0], "resname": c[1]}
                break
    return amap


def split_conjugate_residue(prep_file, frcmod_file, aa_code, cof_code,
                            output_dir=".", linkage=None, atom_source=None,
                            aa_restype="protein"):
    """Split a combined AA↔ligand adduct residue into two tleap libraries.

    Route B parameterizes the covalent adduct (e.g. flavocysteine) as ONE
    RESP-charged residue. For force-field integration we keep the amino acid and
    the cofactor as SEPARATE residues joined by an explicit ``bond`` command
    (the established ProPrep pattern for covalent linkages), so this cuts the
    single residue at the linkage bond into:

    * ``aa_code``  — the amino acid, a chain residue (head N / tail C set so
      tleap peptide-bonds it), with the linkage atom left as an external
      connection point.
    * ``cof_code`` — the cofactor, a standalone unit (no head/tail), linkage
      atom external.

    RESP charges are partitioned as-is (each fragment non-integer, summing to
    the integer total — valid because the two residues always co-occur). The
    frcmod is REUSED UNCHANGED: its parameters are atom-type-keyed, so the
    cross-boundary bond/angle/dihedral terms resolve from the one loaded frcmod.

    ``linkage`` may be given as ``(aa_side_atom_name, cof_side_atom_name)`` to
    pin the cut. Otherwise, if ``atom_source`` (``{atom_name: source_resname}``)
    is supplied, the cut is the graph bridge whose two atoms come from different
    source residues (general — any covalent adduct). Failing both, it falls back
    to the S–CA thioether bridge (sulfur kept on the amino-acid side).

    Returns a dict with ``aa_lib``, ``cof_lib``, ``aa_code``, ``cof_code``,
    ``aa_charge``, ``cof_charge``, ``linkage`` (atom-name pair), ``head``,
    ``tail``, ``aa_atoms`` / ``cof_atoms`` (ordered atom-name lists).
    """
    import collections

    prep = _parse_prep_residue(prep_file)
    atoms, by_seq, by_name = prep["atoms"], prep["by_seq"], prep["by_name"]
    head, tail = prep["head"], prep["tail"]
    if not head or not tail:
        return {"success": False, "error": "Could not identify backbone head/tail in prep"}

    adj = collections.defaultdict(set)
    for e in prep["bonds"]:
        a, b = tuple(e)
        adj[a].add(b); adj[b].add(a)

    # Resolve the linkage bond (a graph bridge between the two source residues).
    link = None
    if linkage:
        try:
            link = (by_name[linkage[0]]["seq"], by_name[linkage[1]]["seq"])
        except KeyError:
            return {"success": False, "error": f"Linkage atoms {linkage} not found in prep"}
    elif atom_source:
        # The cut is the bridge crossing the two source residues.
        for e in prep["bonds"]:
            a, b = tuple(e)
            sa, sb = atom_source.get(by_seq[a]["name"]), atom_source.get(by_seq[b]["name"])
            if sa and sb and sa != sb:
                link = (a, b)
                break
    if link is None:
        for a in atoms:                       # thioether S–CA fallback
            if a["type"] == "S":
                for nb in adj[a["seq"]]:
                    if by_seq[nb]["type"] == "CA":
                        link = (a["seq"], nb)
    if link is None:
        return {"success": False,
                "error": "Could not determine the linkage bond; pass linkage explicitly"}

    def component(start, cut):
        seen, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if {x, y} == set(cut):
                    continue
                if y not in seen:
                    seen.add(y); stack.append(y)
        return seen

    c0, c1 = component(link[0], link), component(link[1], link)
    head_seq = by_name[head]["seq"]
    aa_set = c0 if head_seq in c0 else c1
    cof_set = c1 if aa_set is c0 else c0

    # Cartesian coords (from the Z-matrix) via antechamber; keep prep types/charges.
    tmp_mol2 = os.path.join(output_dir, "_split_coords.mol2")
    subprocess.run(["antechamber", "-i", prep_file, "-fi", "prepi",
                    "-o", tmp_mol2, "-fo", "mol2", "-at", "amber", "-pf", "y"],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    coord = {}
    if os.path.exists(tmp_mol2):
        ml = open(tmp_mol2).read().splitlines()
        ai, bi = ml.index("@<TRIPOS>ATOM"), ml.index("@<TRIPOS>BOND")
        for ln in ml[ai + 1:bi]:
            p = ln.split()
            coord[p[1]] = (p[2], p[3], p[4])
    if len(coord) != len(atoms):
        return {"success": False, "error": "Failed to derive coordinates for residue split"}

    def write_fragment(ids, resname, fn):
        ids = sorted(ids)
        remap = {o: i + 1 for i, o in enumerate(ids)}
        frag_bonds = [tuple(e) for e in prep["bonds"]
                      if all(s in remap for s in e)]
        with open(fn, "w") as f:
            f.write("@<TRIPOS>MOLECULE\n%s\n %d %d 1 0 0\nSMALL\nUSER_CHARGES\n"
                    "@<TRIPOS>ATOM\n" % (resname, len(ids), len(frag_bonds)))
            for o in ids:
                a = by_seq[o]; x, y, z = coord[a["name"]]
                f.write("%7d %-8s%10s%10s%10s %-6s%4d %-8s%10.6f\n"
                        % (remap[o], a["name"], x, y, z, a["type"], 1, resname, a["charge"]))
            f.write("@<TRIPOS>BOND\n")
            for i, (x, y) in enumerate(frag_bonds, 1):
                f.write("%6d %4d %4d 1\n" % (i, remap[x], remap[y]))
            f.write("@<TRIPOS>SUBSTRUCTURE\n     1 %s         1 ****"
                    "               0 ****  **** \n" % resname)

    aa_mol2 = os.path.join(output_dir, f"{aa_code.lower()}.mol2")
    cof_mol2 = os.path.join(output_dir, f"{cof_code.lower()}.mol2")
    write_fragment(aa_set, aa_code, aa_mol2)
    write_fragment(cof_set, cof_code, cof_mol2)

    aa_lib = os.path.join(output_dir, f"{aa_code.lower()}.lib")
    cof_lib = os.path.join(output_dir, f"{cof_code.lower()}.lib")
    restype_line = f"set {aa_code}.1 restype {aa_restype}\n" if aa_restype else ""
    leap = (
        "source leaprc.protein.ff14SB\n"
        "source leaprc.gaff2\n"
        f"loadamberparams {frcmod_file}\n"
        f"{aa_code} = loadmol2 {aa_mol2}\n"
        f"set {aa_code} head {aa_code}.1.{head}\n"
        f"set {aa_code} tail {aa_code}.1.{tail}\n"
        f"set {aa_code}.1 connect0 {aa_code}.1.{head}\n"
        f"set {aa_code}.1 connect1 {aa_code}.1.{tail}\n"
        f"{restype_line}"
        f"saveoff {aa_code} {aa_lib}\n"
        f"{cof_code} = loadmol2 {cof_mol2}\n"
        f"saveoff {cof_code} {cof_lib}\n"
        "quit\n"
    )
    leap_in = os.path.join(output_dir, "_split_mklib.leap")
    open(leap_in, "w").write(leap)
    subprocess.run(["tleap", "-f", leap_in],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if not (os.path.exists(aa_lib) and os.path.exists(cof_lib)):
        return {"success": False, "error": "tleap failed to build split libraries"}

    return {
        "success": True,
        "aa_lib": aa_lib, "cof_lib": cof_lib,
        "aa_code": aa_code, "cof_code": cof_code,
        "aa_charge": round(sum(by_seq[i]["charge"] for i in aa_set), 6),
        "cof_charge": round(sum(by_seq[i]["charge"] for i in cof_set), 6),
        "linkage": (by_seq[link[0]]["name"], by_seq[link[1]]["name"]),
        "head": head, "tail": tail,
        "aa_atoms": [by_seq[i]["name"] for i in sorted(aa_set)],
        "cof_atoms": [by_seq[i]["name"] for i in sorted(cof_set)],
    }


def _rename_pdb_residues_atoms(pdb_path, per_residue, output_path):
    """Rename residue names AND atom names in a PDB (in place or to output).

    ``per_residue`` maps ``(chain, resid) -> {"resname": new, "atoms":
    {old_atom: new_atom}}``. Mirrors the metal-site immediate path
    (structure_preprocessor._apply_atom_name_mappings) so the prepared PDB's
    atoms match the deposited library and tLEaP's loadpdb resolves them.
    """
    out = []
    for line in open(pdb_path):
        # ANISOU/TER records share the atom-name/resname/chain/resid columns of
        # their partner ATOM record and MUST be renamed in lockstep, otherwise
        # the file is left internally inconsistent (ATOM says CF1/N1 while the
        # paired ANISOU still says CYF/N), which corrupts downstream parsing.
        if line.startswith(("ATOM  ", "HETATM", "ANISOU", "TER")):
            chain = line[21] if len(line) > 21 else " "
            try:
                resid = int(line[22:26])
            except (ValueError, IndexError):
                out.append(line); continue
            spec = per_residue.get((chain, resid))
            if spec:
                # TER records carry no atom name; only rename it where present.
                if not line.startswith("TER"):
                    atom = line[12:16].strip()
                    new_atom = spec.get("atoms", {}).get(atom)
                    if new_atom:
                        line = line[:12] + f"{new_atom:>4}" + line[16:]
                line = line[:17] + f"{spec['resname']:>3}" + line[20:]
            out.append(line)
        else:
            out.append(line)
    with open(output_path, "w") as f:
        f.writelines(out)
    return output_path


def _register_for_topology_generator(workspace, lib_file=None, frcmod_file=None,
                                     console=None):
    """Append lib/frcmod paths to the workspace keys the Topology Generator reads.

    ``preprocessing_lib_files`` and ``preprocessing_frcmod_files`` are what
    tleap_input_generator turns into loadoff / loadamberparams lines. A
    parameterizer that skips this deposit leaves its outputs on disk and out of
    the system build. Either argument accepts a single path or an iterable.

    Returns True if anything new was registered.
    """
    if workspace is None:
        return False

    def _as_list(value):
        if value is None:
            return []
        return list(value) if isinstance(value, (list, tuple, set)) else [value]

    def _extend(key, paths):
        paths = [p for p in paths if p]
        if not paths:
            return 0
        existing = workspace.get(key, []) or []
        if not isinstance(existing, list):
            existing = []
        added = 0
        for path in paths:
            path = os.path.abspath(path)
            if path not in existing:
                existing.append(path)
                added += 1
        workspace.set(key, existing)
        return added

    count = _extend("preprocessing_lib_files", _as_list(lib_file))
    count += _extend("preprocessing_frcmod_files", _as_list(frcmod_file))
    if count and console is not None:
        console.print(f"  [green]✓[/green] Registered {count} file(s) with the "
                      f"Topology Generator")
    return count > 0


def integrate_modaa_from_structure(console, workspace, residue_name,
                                   source_residues, capped_pdb, prep_file,
                                   frcmod_file, output_dir, conformer_label="xtal"):
    """Force-field integration for a Route B covalent adduct (metal-site style).

    Splits the combined adduct residue into an amino-acid + cofactor pair,
    deposits both libs + the (reused) frcmod into the user library, renames
    residues AND atoms in the prepared PDB so the Topology Generator can build
    immediately, populates the ``preprocessing_*`` workspace keys, and emits a
    reusable rename transformer. Returns a result dict; best-effort, so a
    failure is reported but the parameterization (already complete) still stands.

    ``source_residues`` is a list of ``{"name", "chain_id", "resid"}`` for the
    original RedoxSite members (their real chain/resid in the prepared PDB).
    ``workspace`` is the processor workspace (``.get``/``.set``) or None.
    """
    from pathlib import Path
    from proprep.forcefield_params.user_library import (
        promote_state, PromotionRequest, DEFAULT_REDOX_STATE, DEFAULT_SPIN_STATE,
    )
    from proprep.forcefield_prep.mcpb.integration_utils import generate_unique_residue_names
    from proprep.forcefield_prep.structure_preprocessor import STANDARD_RESIDUES
    from proprep.redoxsite_prep.transformation.auto_rename import emit_rename_transformer

    aa = residue_name.lower()
    # Topology-only map (identical across conformers) — use the reference
    # conformer's optimization input.
    opt_gjf = os.path.join(output_dir, f"{aa}_{conformer_label}_opt.gjf")
    ac_file = os.path.join(output_dir, f"{residue_name.upper()[:3]}.ac")
    for pth in (opt_gjf, ac_file, capped_pdb, prep_file, frcmod_file):
        if not pth or not os.path.exists(pth):
            return {"success": False,
                    "message": f"FF integration skipped (missing {os.path.basename(pth or '?')})"}

    # 1) atom-name/source map (antechamber name -> {orig, resname}).
    amap = _derive_conjugate_atom_map(opt_gjf, capped_pdb, ac_file)
    if not amap:
        return {"success": False, "message": "Could not derive the atom-name map"}
    atom_source = {k: v["resname"] for k, v in amap.items()}

    # 2) probe split to learn which SOURCE residue is the amino acid.
    probe = split_conjugate_residue(prep_file, frcmod_file, "AA0", "CO0",
                                    output_dir=output_dir, atom_source=atom_source)
    if not probe.get("success"):
        return {"success": False, "message": probe.get("error", "split failed")}
    aa_source = atom_source.get(probe["aa_atoms"][0])
    cof_source = atom_source.get(probe["cof_atoms"][0])

    def _find(nm):
        for r in (source_residues or []):
            if r.get("name") == nm:
                return r
        return None
    aa_res, cof_res = _find(aa_source), _find(cof_source)
    if not aa_res or not cof_res:
        return {"success": False,
                "message": "Could not match split fragments to source residues"}

    # 3) three-letter codes for the split pair.
    #
    #    The amino-acid fragment IS the residue this workflow has been called all
    #    along: the name chosen when the capped model was built, and the name on
    #    the working directory, every file prefix and the AC file. Reuse it
    #    rather than regenerating from the SOURCE resname — regenerating silently
    #    discarded a name the user deliberately picked (choose "K12" at capping
    #    and the deposited library unit would still come out "CS1").
    #
    #    Only the cofactor needs a fresh code, and it is generated against an
    #    explicit reserved set: the names tLEaP already resolves, the adduct's
    #    own code, and the source residues. (generate_unique_residue_names ADDS
    #    what it generates to the set it is handed, so pass a copy.)
    from proprep.forcefield_params.loader import get_registered_residue_names
    aa_code = (residue_name or aa_source).strip().upper()[:3]
    reserved = set(STANDARD_RESIDUES) | {"HOH", "WAT", aa_code}
    reserved |= {(r.get("name") or "").strip().upper()
                 for r in (source_residues or []) if r.get("name")}
    # Also avoid names already registered in the forcefield library, so the split
    # cofactor does not collide with an unrelated library residue at deposit. The
    # adduct's own aa_code is kept in `reserved` above even if a prior run of the
    # SAME residue registered it (on_collision="overwrite" handles that leaf).
    reserved |= set(get_registered_residue_names())
    cof_code = generate_unique_residue_names(
        [(cof_res["resid"], cof_source)], existing_names=set(reserved)
    ).get((cof_res["resid"], cof_source)) or f"{cof_source[:2].upper()}1"

    # 4) final split into the two libraries.
    split = split_conjugate_residue(prep_file, frcmod_file, aa_code, cof_code,
                                    output_dir=output_dir, atom_source=atom_source)
    if not split.get("success"):
        return {"success": False, "message": split.get("error", "split failed")}
    aa_lib, cof_lib = split["aa_lib"], split["cof_lib"]

    # 5) deposit both libs + the reused frcmod into the user library.
    req = PromotionRequest(
        family="modified_aa", type_name=residue_name,
        set_name=f"{residue_name}_RESP",
        frcmod_src=frcmod_file, lib_srcs=[aa_lib, cof_lib],
        redox_state=DEFAULT_REDOX_STATE, spin_state=DEFAULT_SPIN_STATE,
        set_meta={"description": f"Route B covalent adduct {residue_name} "
                                 f"({aa_source}->{aa_code} + {cof_source}->{cof_code})"},
        # residue_name + atom_types live at the spin level (see _splice_leaf).
        spin_meta={"residue_name": {aa_source: aa_code, cof_source: cof_code}},
        on_collision="overwrite",
    )
    dep = promote_state(req)
    library_path = dep.get("library_path")
    console.print(f"  [green]✓[/green] Library: [grey50]{library_path}[/grey50]")

    forcefield_path = None
    if library_path:
        parts = Path(library_path).parts
        if "specialized_residues" in parts:
            i = parts.index("specialized_residues")
            forcefield_path = "/".join(parts[i + 1:]) or None

    # inverse per-fragment atom maps: original PDB name -> library name.
    aa_atoms = {v["orig"]: k for k, v in amap.items() if v["resname"] == aa_source}
    cof_atoms = {v["orig"]: k for k, v in amap.items() if v["resname"] == cof_source}

    # 6) immediate use: rename residues + atoms in the prepared PDB.
    # Use the canonical priority-based selector (same source the Topology
    # Generator's loadpdb will resolve to) rather than guessing workspace keys.
    prepared = None
    if workspace is not None:
        try:
            from proprep.utils.structure_selector import get_priority_pdb_file
            prepared = get_priority_pdb_file(workspace, silent=True)
        except Exception:
            prepared = None
        if not prepared:
            prepared = (workspace.get("prepared_pdb", None)
                        or workspace.get("transformed_pdb_file", None)
                        or workspace.get("transformed_structure", None))
    if prepared and os.path.exists(prepared):
        per_res = {
            (aa_res["chain_id"], aa_res["resid"]): {"resname": aa_code, "atoms": aa_atoms},
            (cof_res["chain_id"], cof_res["resid"]): {"resname": cof_code, "atoms": cof_atoms},
        }
        # NEVER overwrite the source (it may be the user's original crystal
        # PDB). Write the renamed structure to a sibling file and register it
        # as the highest-priority prepared_pdb so the Topology Generator's
        # loadpdb resolves to it.
        stem = Path(prepared).stem
        renamed = os.path.join(os.path.dirname(prepared) or ".",
                               f"{stem}_modaa_renamed.pdb")
        _rename_pdb_residues_atoms(prepared, per_res, renamed)
        if workspace is not None:
            workspace.set("prepared_pdb", renamed)
        console.print(f"  [green]✓[/green] Renamed residues/atoms (new file, original preserved): "
                      f"[grey50]{renamed}[/grey50]")

        # 6b) Sync the workspace RedoxSite objects to the renamed PDB. The rename
        # changed residue names (CYF->CF1, RBF->RF1) and atom names (->lib
        # names) but NOT coordinates, so the detected_redox_sites objects are
        # now stale: PDB Filter / Topology Generator would emit tLEaP `bond`
        # commands against the old names (e.g. `bond mol.44.SG ...` when the lib
        # atom is now `S1`), which tLEaP rejects. Coordinate-based sync re-keys
        # every site's atoms/centers/coord_to_pdb/bonds/residue_groups off the
        # renamed structure. Mirrors the metal-site path.
        if workspace is not None:
            sites = workspace.get("detected_redox_sites") or []
            if sites:
                try:
                    from proprep.structure_prep.comprehensive_redox_detector import (
                        sync_redox_sites_from_pdb,
                    )
                    from proprep.forcefield_prep.structure_preprocessor import (
                        _ensure_redox_site_objects,
                    )
                    # detected_redox_sites round-trips through JSON workspace state,
                    # so on a RESUMED session the sites come back as serialization
                    # wrapper dicts ({"__type__": "RedoxSite", ...}), not RedoxSite
                    # objects — the checklist deserializer leaves complex types for
                    # the caller to reconstruct. sync_redox_sites_from_pdb needs
                    # objects (it reads site.atoms), so normalize first (a no-op for
                    # object-form sites) and store the objects back so downstream
                    # consumers match a fresh, non-resumed run.
                    sites = _ensure_redox_site_objects(sites)
                    changed = sync_redox_sites_from_pdb(renamed, sites)
                    workspace.set("detected_redox_sites", sites)
                    if changed:
                        console.print("  [green]✓[/green] Synced RedoxSite definitions to renamed "
                                      "residue/atom names (bond commands stay valid).")
                except Exception as exc:  # noqa: BLE001
                    console.print(f"  [yellow]Note: RedoxSite sync skipped ({exc}); verify tLEaP "
                                  f"bond commands reference the renamed atoms.[/yellow]")
    else:
        console.print("[grey50]FF integration: no prepared PDB in workspace to rename "
                      "(Topology Generator still loads the libs).[/grey50]")

    # 7) workspace keys the Topology Generator consumes.
    _register_for_topology_generator(
        workspace, lib_file=[aa_lib, cof_lib], frcmod_file=frcmod_file,
    )

    # 8) reusable rename transformer (Tier-1, with atom_renames).
    tpath = None
    try:
        tpath = emit_rename_transformer(
            [{"resname": aa_source, "target": aa_code, "atom_renames": aa_atoms},
             {"resname": cof_source, "target": cof_code, "atom_renames": cof_atoms}],
            name=f"modaa_{residue_name}_{aa_code}_{cof_code}".lower(),
            description=f"Reuse {residue_name} modified-AA adduct "
                        f"({aa_source}->{aa_code} + {cof_source}->{cof_code})",
            # Bake the SAME states the deposit used so the Topology Generator's
            # FF lookup resolves modified_aa/<name>/<redox>/<spin> to the
            # deposited set (was defaulting to 'unknown/unknown' → dir not found).
            redox_state=DEFAULT_REDOX_STATE,
            spin_state=DEFAULT_SPIN_STATE,
            forcefield_path=forcefield_path,
            provenance={"source": "modaa_from_structure", "residue": residue_name,
                        "library_path": library_path,
                        "metadata_path": dep.get("metadata_path")},
            site_types=["modified_aa"],
        )
        console.print(f"  [green]✓[/green] Reuse transformer: [grey50]{tpath}[/grey50]")
    except Exception as exc:  # noqa: BLE001
        console.print(f"  [yellow]Note: reuse transformer not created ({exc}); params still saved.[/yellow]")

    return {
        "success": True,
        "message": f"Integrated {aa_source}->{aa_code} (chain AA) + {cof_source}->{cof_code} (cofactor)",
        "library_path": library_path, "aa_code": aa_code, "cof_code": cof_code,
        "aa_lib": aa_lib, "cof_lib": cof_lib, "frcmod_file": frcmod_file,
        "transformer": tpath,
    }


def generate_amber_library(prep_file, frcmod_file, residue_name, show_next_steps=True):
    """
    Generate an AMBER library file (.lib) from prep and frcmod files using tleap.

    ``show_next_steps`` controls the "copy these files / load them in tleap"
    panel. It is correct for the de-novo (Route A) standalone library, but
    misleading for the from-structure (Route B) adduct, where this single lib is
    only an intermediate that step 10 splits + deposits + wires into a
    transformer — so Route B passes False.

    Parameters:
    -----------
    prep_file : str
        Path to the prep file
    frcmod_file : str
        Path to the frcmod file
    residue_name : str
        Name of the residue (used for library naming)

    Returns:
    --------
    dict
        Dictionary containing the result of the library generation
    """
    import os
    import subprocess

    # Normalize residue name for filenames
    res_lower = residue_name.lower()
    lib_file = f"{res_lower}.lib"
    tleap_input = f"{res_lower}_savelib.in"

    # Generate tleap input file
    tleap_content = f"""# tleap input to create AMBER library for {residue_name}
# Generated by ProPrep

# Load the frcmod file with bonded parameters
loadAmberParams "{frcmod_file}"

# Load the prep file with residue definition and charges
loadAmberPrep "{prep_file}"

# Save as AMBER library file
saveOff {residue_name} "{lib_file}"

quit
"""

    _console.print(f"\n[bold cyan]Step 10: Generating AMBER Library File[/bold cyan]")
    _console.print(f"Creating tleap input: {tleap_input}")

    try:
        with open(tleap_input, "w") as f:
            f.write(tleap_content)

        # Show what we're doing with a Panel
        _console.print(Panel(
            f"[bold]Command:[/bold] tleap -f {tleap_input}\n\n"
            f"Converting prep file to AMBER library (OFF) format.\n"
            f"The frcmod is loaded to validate parameters during conversion.\n\n"
            f"[bold]tleap commands:[/bold]\n"
            f'  [cyan]loadAmberParams[/cyan] "{frcmod_file}"\n'
            f'  [cyan]loadAmberPrep[/cyan]   "{prep_file}"\n'
            f'  [cyan]saveOff[/cyan]         {residue_name} "{lib_file}"',
            title="Running tleap",
            border_style="grey50",
            expand=False
        ))

        # Run tleap
        process = subprocess.run(
            ["tleap", "-f", tleap_input],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Check if library was created
        if os.path.exists(lib_file):
            _console.print(f"\n[green]✓ Successfully created library file: {lib_file}[/green]")
            if show_next_steps:
                _console.print(Panel(
                    f"[bold]To use this residue in your simulations:[/bold]\n\n"
                    f"Copy these files to your working directory:\n"
                    f"   • [cyan]{frcmod_file}[/cyan] - force field parameters (bonds, angles, dihedrals)\n"
                    f"   • [cyan]{lib_file}[/cyan] - residue topology and charges\n\n"
                    f"In tleap, load them with:\n"
                    f'   [cyan]loadAmberParams "{frcmod_file}"[/cyan]\n'
                    f'   [cyan]loadOff "{lib_file}"[/cyan]\n\n'
                    f"The residue '[bold]{residue_name}[/bold]' will then be available for use.\n\n"
                    f"[grey50]Note: The lib file contains the residue definition (same info as {prep_file}).\n"
                    f"You need both lib and frcmod because lib has topology/charges while\n"
                    f"frcmod has the actual force field parameters.[/grey50]",
                    title="Next Steps",
                    border_style="green",
                    expand=False
                ))

            return {
                "success": True,
                "lib_file": lib_file,
                "tleap_input": tleap_input,
                "message": f"Successfully created {lib_file}",
            }
        else:
            _console.print(f"[yellow]⚠ Warning: tleap completed but {lib_file} was not created[/yellow]")
            _console.print("[grey50]tleap output:[/grey50]")
            _console.print(process.stdout)
            if process.stderr:
                _console.print("[grey50]tleap errors:[/grey50]")
                _console.print(process.stderr)

            return {
                "success": False,
                "error": f"Library file {lib_file} was not created",
                "tleap_output": process.stdout,
                "tleap_stderr": process.stderr,
            }

    except FileNotFoundError:
        _console.print("[red]✗ Error: tleap command not found. Make sure AmberTools is installed and in your PATH.[/red]")
        return {
            "success": False,
            "error": "tleap command not found",
        }
    except Exception as e:
        _console.print(f"[red]✗ Error running tleap: {str(e)}[/red]")
        return {
            "success": False,
            "error": str(e),
        }


def run_parmchk2(prep_file, frcmod_file, parm_set=None, frc_file=None, parm_dat_file=None):
    """
    Run parmchk2 on a prep file to generate an frcmod file.

    Parameters:
    -----------
    prep_file : str
        Path to the prep file
    frcmod_file : str
        Path to the output frcmod file
    parm_set : str, optional
        Parameter set shortcut for -s flag (gaff, gaff2, parm99, parm10, lipid14)
    frc_file : str, optional
        Additional frcmod to load via -frc flag (ff99SB, ff14SB, ff03, etc.)
    parm_dat_file : str, optional
        Full path to parameter file for -p flag (used when -s shortcut unavailable)

    Returns:
    --------
    dict
        Dictionary containing the result of the parmchk2 command
    """
    import os
    import subprocess

    try:
        # Base command
        cmd = ["parmchk2", "-i", prep_file, "-f", "prepi", "-o", frcmod_file, "-a", "Y"]

        # Add parameter set shortcut if specified
        if parm_set:
            cmd.extend(["-s", parm_set])

        # Add additional frcmod if specified
        if frc_file:
            cmd.extend(["-frc", frc_file])

        # Add full parameter file path if specified (overrides -s)
        if parm_dat_file:
            cmd.extend(["-p", parm_dat_file])

        # Build explanation text
        parm_set_desc = {
            "gaff": "GAFF - General Amber Force Field",
            "gaff2": "GAFF2 - General Amber Force Field 2",
            "parm99": "parm99.dat - base for ff99SB",
            "parm10": "parm10.dat - base for ff14SB",
            "lipid14": "lipid14.dat - lipid force field",
        }
        frc_desc = {
            "ff99SB": "ff99SB protein corrections",
            "ff14SB": "ff14SB protein corrections (torsions, etc.)",
            "ff03": "ff03 protein corrections",
            "bsc1": "bsc1 DNA corrections",
            "ol15": "OL15 DNA corrections",
            "ol3": "OL3 RNA corrections",
            "yil": "YIL RNA corrections",
        }

        # Build flags explanation
        flags_text = (
            f"[cyan]-i[/cyan] {prep_file:<20} Input prep file with atom types and charges\n"
            f"[cyan]-f[/cyan] prepi{' '*15} Input format (AMBER prep)\n"
            f"[cyan]-o[/cyan] {frcmod_file:<20} Output frcmod file for missing/modified parameters\n"
            f"[cyan]-a[/cyan] Y{' '*19} Print all parameters including those already in parm file"
        )

        if parm_dat_file:
            flags_text += f"\n[cyan]-p[/cyan] {os.path.basename(parm_dat_file):<20} Base parameter file (full path, overrides -s)"
        elif parm_set:
            desc = parm_set_desc.get(parm_set, parm_set)
            flags_text += f"\n[cyan]-s[/cyan] {parm_set:<20} {desc}"

        if frc_file:
            desc = frc_desc.get(frc_file, frc_file)
            flags_text += f"\n[cyan]-frc[/cyan] {frc_file:<18} Load additional {desc}"

        # Display command with explanations using Panel
        _console.print(Panel(
            f"[bold]Command:[/bold] {' '.join(cmd)}\n\n"
            f"[bold]Flags explained:[/bold]\n{flags_text}",
            title="Running parmchk2",
            border_style="grey50",
            expand=False
        ))

        # Run the command
        process = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )

        # Check if output file was created
        if os.path.exists(frcmod_file):
            return {
                "success": True,
                "message": f"Successfully created {frcmod_file}",
                "output": process.stdout,
            }
        else:
            return {
                "success": False,
                "error": f"parmchk2 did not create the output file {frcmod_file}",
                "output": process.stdout,
                "stderr": process.stderr,
            }

    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "error": f"Error running parmchk2: {e.stderr}",
            "output": e.stdout,
        }
    except FileNotFoundError:
        return {
            "success": False,
            "error": "parmchk2 command not found. Make sure AmberTools is installed and in your PATH.",
        }


def _all_bond_angle_params(frcmod_file):
    """Every BOND and ANGLE term present in the frcmod.

    Returns ``(param_name, 0.0, 'full', section)`` tuples in the same shape
    ``analyze_frcmod_penalties`` produces, so they drop straight into the
    Seminario pipeline. Used by the Full-Seminario scope: refine the residue's
    COMPLETE custom bond/angle set from the QM Hessian, not just the penalty/
    ATTN-flagged subset. (Standard protein-FF terms are never written to a
    parmchk2 frcmod, so "all terms in the frcmod" is exactly the residue's
    nonstandard bond/angle set — the standard backbone terms stay ff14SB.)
    """
    from proprep.forcefield_prep.small_molecule_parameterizer import extract_parameter_name
    out, seen, section = [], set(), None
    try:
        with open(frcmod_file) as f:
            for raw in f:
                s = raw.strip()
                if not s or s.startswith('#'):
                    continue
                if s in ('MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON', 'END'):
                    section = s
                    continue
                if section in ('BOND', 'ANGLE'):
                    name = extract_parameter_name(s, section)
                    if name and name not in seen:
                        seen.add(name)
                        out.append((name, 0.0, 'full', section))
    except OSError:
        return []
    return out


def generate_bonded_parameters(residue_symbol=None, processor=None, standalone_use=True):
    """
    Find prep files, generate bonded parameters using parmchk2, and improve parameters
    by replacing ATTN warnings with GAFF parameters.

    Parameters:
    -----------
    residue_symbol : str, optional
        Three-letter code for the custom residue. If not provided, will try to detect from prep files.
    processor : optional
        Processor object with session_manager for session recording support.

    Returns:
    --------
    dict
        Dictionary containing the result of the parameter generation
    """
    import os

    _console.print("\n[bold cyan]Step 9: Generating Bonded Parameters Using parmchk2[/bold cyan]")

    # Find all prep files in the current directory
    prep_files = find_prep_files()

    if not prep_files:
        _console.print("[yellow]○[/yellow] No prep files found in the current directory.")
        prep_file = prompt_with_context(
            processor, "Please enter the path to the prep file", default="",
            module="Modified Amino Acid Parameterizer", description="Path to prep file"
        )
        if os.path.exists(prep_file) and prep_file.endswith(".prep"):
            prep_files = [prep_file]
        else:
            return {"success": False, "error": f"Invalid prep file: {prep_file}"}

    # If multiple prep files found, let user choose
    if len(prep_files) > 1:
        _console.print("\n[cyan]Multiple prep files found:[/cyan]")
        for i, file in enumerate(prep_files):
            _console.print(f"  [cyan]{i+1}.[/cyan] {file}")

        try:
            choice_str = prompt_with_context(
                processor, "Select a prep file by number", default="1",
                choices=[str(i) for i in range(1, len(prep_files) + 1)],
                module="Modified Amino Acid Parameterizer", description="Select prep file"
            )
            choice = int(choice_str)
            if 1 <= choice <= len(prep_files):
                selected_prep = prep_files[choice - 1]
            else:
                _console.print(f"[yellow]⚠[/yellow] Invalid choice. Using the first file: {prep_files[0]}")
                selected_prep = prep_files[0]
        except ValueError:
            _console.print(f"[yellow]⚠[/yellow] Invalid input. Using the first file: {prep_files[0]}")
            selected_prep = prep_files[0]
    else:
        selected_prep = prep_files[0]

    _console.print(f"\n[green]✓[/green] Selected prep file: {selected_prep}")

    # Extract base filename from prep file for output files
    prep_base_name = os.path.splitext(os.path.basename(selected_prep))[0]

    # Extract residue code from prep file if not provided
    if not residue_symbol:
        # Try to extract from filename
        residue_symbol = prep_base_name.upper()

        # Confirm with user
        residue_symbol = prompt_with_context(
            processor, "Enter residue symbol/name for the frcmod file",
            default=residue_symbol,
            module="Modified Amino Acid Parameterizer", description="Residue symbol"
        )

    # Generate output frcmod filename based on prep file name, not residue symbol
    frcmod_file = f"{prep_base_name}_temp.frcmod"

    # Show force field selection menu with Panel
    _console.print(Panel(
        "[bold]Force Field Options:[/bold]\n\n"
        "[cyan]1.[/cyan] ff19SB [green](Recommended)[/green]\n"
        "   Most modern protein force field with improved backbone parameters\n\n"
        "[cyan]2.[/cyan] ff14SB\n"
        "   Widely used, well-validated protein force field\n\n"
        "[cyan]3.[/cyan] ff99SB\n"
        "   Legacy force field, use for compatibility with older simulations\n\n"
        "[cyan]4.[/cyan] GAFF2\n"
        "   General Amber Force Field - for non-protein applications",
        title="Select Target Force Field",
        border_style="blue",
        expand=False
    ))

    try:
        choice = prompt_with_context(
            processor, "Choice", choices=["1", "2", "3", "4"], default="1",
            module="Modified Amino Acid Parameterizer", description="Select force field",
            options_map={"1": "ff19SB", "2": "ff14SB", "3": "ff99SB", "4": "GAFF2"}
        )

        # Map choice to parmchk2 parameters
        parm_set = None
        frc_file = None
        parm_dat_file = None

        if choice == "1":
            # ff19SB - not in -s shortcuts, need full path
            amberhome = os.environ.get("AMBERHOME", "")
            parm_dat_file = os.path.join(amberhome, "dat", "leap", "parm", "parm19.dat")
            if not os.path.exists(parm_dat_file):
                _console.print(f"[yellow]⚠ Warning: {parm_dat_file} not found. Using GAFF2 default.[/yellow]")
                parm_set = "gaff2"
                parm_dat_file = None
        elif choice == "2":
            # ff14SB - use parm10 base with ff14SB corrections
            parm_set = "parm10"
            frc_file = "ff14SB"
        elif choice == "3":
            # ff99SB - use parm99 base with ff99SB corrections
            parm_set = "parm99"
            frc_file = "ff99SB"
        elif choice == "4":
            # GAFF2 - just use the shortcut
            parm_set = "gaff2"
        else:
            _console.print(f"[yellow]⚠ Invalid choice '{choice}'. Using GAFF2 default.[/yellow]")
            parm_set = "gaff2"

    except ValueError:
        _console.print("[yellow]⚠ Invalid input. Using GAFF2 default.[/yellow]")
        parm_set = "gaff2"
        frc_file = None
        parm_dat_file = None

    # Build description of parameter set used
    if parm_dat_file:
        param_description = os.path.basename(parm_dat_file)
    elif parm_set and frc_file:
        param_description = f"{parm_set} + {frc_file}"
    elif parm_set:
        param_description = parm_set
    else:
        param_description = "default (gaff)"

    # Run parmchk2
    parmchk2_result = run_parmchk2(
        selected_prep, frcmod_file,
        parm_set=parm_set, frc_file=frc_file, parm_dat_file=parm_dat_file
    )

    if not parmchk2_result["success"]:
        _console.print(
            f"\n[red]✗ Error generating bonded parameters: {parmchk2_result.get('error', 'Unknown error')}[/red]"
        )
        return parmchk2_result

    _console.print(f"\n[green]✓ {parmchk2_result['message']}[/green]")

    import shutil

    final_frcmod = f"{prep_base_name}.frcmod"

    # Show EVERY parameter parmchk2 flagged — penalty scores AND ATTN — so the
    # user can see poorly-transferred terms (high penalty) that are not ATTN.
    from proprep.forcefield_prep.small_molecule_parameterizer import analyze_frcmod_penalties
    analyze_frcmod_penalties(frcmod_file, _console)

    temp_index = _index_frcmod_params(frcmod_file)
    temp_candidates = _group_frcmod_candidates(temp_index)

    if not temp_candidates:
        _console.print(f"[green]✓[/green] No penalty scores or ATTN warnings in {frcmod_file}. No improvements needed.")
        shutil.copy(frcmod_file, final_frcmod)
        lib_result = generate_amber_library(selected_prep, final_frcmod, residue_symbol, show_next_steps=standalone_use)
        return {
            "lib_file": lib_result.get("lib_file") if lib_result.get("success") else None,
            "success": True,
            "message": "Successfully generated bonded parameters (no penalties or ATTN warnings)",
            "prep_file": selected_prep,
            "frcmod_file": frcmod_file,
            "final_frcmod": final_frcmod,
            "improved": False,
            "param_file": param_description,
        }

    improve = confirm_with_context(
        processor,
        "\nWould you like to review these against GAFF2 and replace selected parameters?",
        default=True, module="Modified Amino Acid Parameterizer",
        description="Improve parameters with GAFF2"
    )

    if not improve:
        _console.print(
            "\n[yellow]⚠ Skipping parameter improvement. Penalty/ATTN parameters remain in the frcmod file.[/yellow]"
        )
        shutil.copy(frcmod_file, final_frcmod)
        lib_result = generate_amber_library(selected_prep, final_frcmod, residue_symbol, show_next_steps=standalone_use)
        return {
            "success": True,
            "message": "Successfully generated bonded parameters",
            "prep_file": selected_prep,
            "frcmod_file": frcmod_file,
            "final_frcmod": final_frcmod,
            "lib_file": lib_result.get("lib_file") if lib_result.get("success") else None,
            "param_file": param_description,
            "improved": False,
        }

    # Generate the GAFF2 frcmod (parmchk2 maps the AMBER atom types to GAFF by
    # analogy, so its parameter keys match the protein-FF frcmod's).
    gaff_frcmod = f"{prep_base_name}_gaff.frcmod"
    _console.print(f"\n[cyan]→[/cyan] Generating GAFF2 parameters for {residue_symbol}...")
    gaff_result = run_parmchk2_gaff(selected_prep, gaff_frcmod)

    if not gaff_result["success"]:
        _console.print(
            f"[red]✗ Error generating GAFF2 parameters: {gaff_result.get('error', 'Unknown error')}[/red]"
        )
        _console.print(f"\n[yellow]⚠ Keeping the protein-FF frcmod {frcmod_file} (penalty/ATTN terms remain).[/yellow]")
        shutil.copy(frcmod_file, final_frcmod)
        lib_result = generate_amber_library(selected_prep, final_frcmod, residue_symbol, show_next_steps=standalone_use)
        return {
            "success": True,
            "message": "Generated protein-FF parameters but failed to generate GAFF2 parameters",
            "prep_file": selected_prep,
            "frcmod_file": frcmod_file,
            "final_frcmod": final_frcmod,
            "lib_file": lib_result.get("lib_file") if lib_result.get("success") else None,
            "param_file": param_description,
            "has_penalty_warnings": True,
        }

    _console.print(f"[green]✓[/green] GAFF2 parameters written to {gaff_frcmod}.")

    gaff_index = _index_frcmod_params(gaff_frcmod)
    selected_keys = select_gaff_replacements(temp_index, gaff_index, processor, interactive=True)

    if not selected_keys:
        _console.print("\n[grey50]No parameters selected for replacement; keeping the protein-FF frcmod.[/grey50]")
        shutil.copy(frcmod_file, final_frcmod)
        lib_result = generate_amber_library(selected_prep, final_frcmod, residue_symbol, show_next_steps=standalone_use)
        return {
            "success": True,
            "message": "Successfully generated bonded parameters (no GAFF2 replacements selected)",
            "prep_file": selected_prep,
            "frcmod_file": frcmod_file,
            "gaff_frcmod": gaff_frcmod,
            "final_frcmod": final_frcmod,
            "lib_file": lib_result.get("lib_file") if lib_result.get("success") else None,
            "param_file": param_description,
            "improved": False,
        }

    _console.print(f"\n[cyan]→[/cyan] Splicing {len(selected_keys)} GAFF2 parameter(s) into {final_frcmod}...")
    combine_result = create_penalty_aware_combined_frcmod(
        temp_index, gaff_index, final_frcmod, selected_keys
    )

    if not combine_result["success"]:
        _console.print(
            f"[red]✗ Error creating combined frcmod file: {combine_result.get('error', 'Unknown error')}[/red]"
        )
        _console.print(f"\n[yellow]⚠ Using original frcmod file {frcmod_file} (penalty/ATTN terms remain).[/yellow]")
        shutil.copy(frcmod_file, final_frcmod)
        lib_result = generate_amber_library(selected_prep, final_frcmod, residue_symbol, show_next_steps=standalone_use)
        return {
            "success": True,
            "message": "Generated protein-FF and GAFF2 parameters but failed to combine them",
            "prep_file": selected_prep,
            "frcmod_file": frcmod_file,
            "gaff_frcmod": gaff_frcmod,
            "final_frcmod": final_frcmod,
            "lib_file": lib_result.get("lib_file") if lib_result.get("success") else None,
            "param_file": param_description,
        }

    replacements = combine_result.get("replacements", 0)
    remaining = len(temp_candidates) - replacements
    _console.print(f"[green]✓[/green] Wrote combined frcmod file {final_frcmod}.")
    _console.print(
        f"  Replaced {replacements} of {len(temp_candidates)} penalized/ATTN parameter(s) with GAFF2."
    )
    if remaining > 0:
        _console.print(
            f"[grey50]  {remaining} penalized/ATTN parameter(s) were left as the protein-FF value "
            "(not selected, or GAFF2 offered no improvement).[/grey50]"
        )

    # Show the penalty table of the COMBINED frcmod so the user sees the result of
    # the splice — GAFF2 lines carry their own penalty scores, so the remaining
    # high-penalty/ATTN terms (if any) are exactly what still needs attention
    # (e.g. Seminario). Without this the user only sees the pre-merge table and a
    # replacement count, never the improved state they actually ship.
    _console.print("\n[bold]Combined frcmod after GAFF2 replacement:[/bold]")
    analyze_frcmod_penalties(final_frcmod, _console)

    lib_result = generate_amber_library(selected_prep, final_frcmod, residue_symbol, show_next_steps=standalone_use)

    return {
        "success": True,
        "message": "Successfully generated and improved bonded parameters",
        "prep_file": selected_prep,
        "frcmod_file": frcmod_file,
        "gaff_frcmod": gaff_frcmod,
        "final_frcmod": final_frcmod,
        "lib_file": lib_result.get("lib_file") if lib_result.get("success") else None,
        "penalized_parameters": len(temp_candidates),
        "replacements": replacements,
        "improved": True,
        "param_file": param_description,
    }


def count_conformations_in_esp_file(esp_file):
    """
    Count the number of conformations in a concatenated ESP file.

    Parameters:
    -----------
    esp_file : str
        Path to the ESP file

    Returns:
    --------
    int
        Number of conformations in the ESP file
    """
    try:
        with open(esp_file, "r") as f:
            lines = f.readlines()

        # Count ESP file headers which are lines containing exactly three integers
        # For example: "24 7372    0" where the first integer is the number of atoms
        conf_count = 0

        for line in lines:
            # Strip any leading/trailing whitespace
            line = line.strip()

            # Split the line into parts
            parts = line.split()

            # Check if it has exactly three parts
            if len(parts) == 3:
                # Try to convert all three parts to integers
                try:
                    int(parts[0])
                    int(parts[1])
                    int(parts[2])
                    # This looks like an ESP file header
                    conf_count += 1
                except ValueError:
                    # Not all parts are integers
                    continue

        if conf_count > 0:
            return conf_count

    except Exception as e:
        _console.print(f"[red]✗ Error counting conformations in ESP file: {str(e)}[/red]")
        return 0


def analyze_ac_file_for_separator_bonds(ac_file):
    """
    Analyze the AC file to identify separator bonds between the residue and caps.
    The function assumes ACE-XXX-NME structure, where:
    - The first 6 atoms are the ACE cap
    - The final 8 atoms are the NME cap

    Parameters:
    -----------
    ac_file : str
        Path to the AC file

    Returns:
    --------
    dict
        Dictionary containing identified separator bonds
    """
    try:
        with open(ac_file, "r") as f:
            lines = f.readlines()

        # Parse atoms
        atoms = []
        for line in lines:
            if line.startswith("ATOM"):
                parts = line.split()
                atom_id = int(parts[1])
                atom_name = parts[2]
                atoms.append(
                    {
                        "id": atom_id,
                        "name": atom_name,
                        "index": len(atoms),  # Track position in file
                    }
                )

        # Define cap sizes
        ace_cap_size = 6   # ACE: HH31, CH3, HH32, HH33, C, O
        nme_cap_size = 6   # NME: N, H, CH3, HH31, HH32, HH33

        # Assign regions to atoms
        for i, atom in enumerate(atoms):
            if i < ace_cap_size:
                atom["region"] = "ACE"
            elif i >= len(atoms) - nme_cap_size:
                atom["region"] = "NME"
            else:
                atom["region"] = "RESIDUE"

        # Create a map from atom name to atom for each region
        ace_atoms = {atom["name"]: atom for atom in atoms if atom["region"] == "ACE"}
        residue_atoms = {
            atom["name"]: atom for atom in atoms if atom["region"] == "RESIDUE"
        }
        nme_atoms = {atom["name"]: atom for atom in atoms if atom["region"] == "NME"}

        # Parse bonds based on atom names, not IDs
        bonds = []
        for line in lines:
            if line.startswith("BOND"):
                parts = line.split()
                if len(parts) >= 7:
                    bond_id = int(parts[1])
                    atom1_name = parts[5]
                    atom2_name = parts[6]

                    # Try to find these atoms
                    atom1 = None
                    atom2 = None

                    # Look in all regions
                    for region in [ace_atoms, residue_atoms, nme_atoms]:
                        if atom1_name in region:
                            atom1 = region[atom1_name]
                        if atom2_name in region:
                            atom2 = region[atom2_name]

                    # Only add bond if we found both atoms
                    if atom1 and atom2:
                        bonds.append(
                            {
                                "id": bond_id,
                                "atom1": atom1,
                                "atom2": atom2,
                                "atom1_name": atom1_name,
                                "atom2_name": atom2_name,
                            }
                        )

        # Find bonds connecting different regions
        region_crossing_bonds = [
            bond for bond in bonds if bond["atom1"]["region"] != bond["atom2"]["region"]
        ]

        # Find C-N bonds connecting different regions
        cn_bonds = []
        for bond in region_crossing_bonds:
            atom1_name = bond["atom1_name"]
            atom2_name = bond["atom2_name"]

            if atom1_name.startswith("C") and atom2_name.startswith("N"):
                cn_bonds.append(
                    {
                        "id": bond["id"],
                        "carbon": bond["atom1"],
                        "nitrogen": bond["atom2"],
                        "carbon_name": atom1_name,
                        "nitrogen_name": atom2_name,
                    }
                )
            elif atom2_name.startswith("C") and atom1_name.startswith("N"):
                cn_bonds.append(
                    {
                        "id": bond["id"],
                        "carbon": bond["atom2"],
                        "nitrogen": bond["atom1"],
                        "carbon_name": atom2_name,
                        "nitrogen_name": atom1_name,
                    }
                )

        # Find ACE-to-residue bond (carbon in ACE, nitrogen in residue)
        ace_to_residue = [
            bond
            for bond in cn_bonds
            if bond["carbon"]["region"] == "ACE"
            and bond["nitrogen"]["region"] == "RESIDUE"
        ]

        # Find residue-to-NME bond (carbon in residue, nitrogen in NME)
        residue_to_nme = [
            bond
            for bond in cn_bonds
            if bond["carbon"]["region"] == "RESIDUE"
            and bond["nitrogen"]["region"] == "NME"
        ]

        # Format the bonds for return
        suggested_bonds = []

        if ace_to_residue:
            # Pick the bond with the carbon atom that has the highest index (last one in ACE cap)
            ace_to_residue.sort(key=lambda b: b["carbon"]["index"], reverse=True)
            ace_bond = ace_to_residue[0]

            # Format the bond for return - using atom names
            # Make the first atom from residue, second atom from cap
            suggested_bonds.append(
                {
                    "atom1_name": ace_bond["nitrogen_name"],
                    "atom2_name": ace_bond["carbon_name"],
                    "atom1_idx": ace_bond["nitrogen"]["id"],
                    "atom2_idx": ace_bond["carbon"]["id"],
                    "type": "ACE_TO_RESIDUE",
                }
            )

            _console.print(
                f"[green]✓[/green] Found N-terminus separator bond: [cyan]{ace_bond['nitrogen_name']} - {ace_bond['carbon_name']}[/cyan]"
            )

        if residue_to_nme:
            # Pick the bond with the carbon atom that has the highest index (last one in residue)
            residue_to_nme.sort(key=lambda b: b["carbon"]["index"], reverse=True)
            nme_bond = residue_to_nme[0]

            # Format the bond for return - using atom names
            # Make the first atom from residue, second atom from cap
            suggested_bonds.append(
                {
                    "atom1_name": nme_bond["carbon_name"],
                    "atom2_name": nme_bond["nitrogen_name"],
                    "atom1_idx": nme_bond["carbon"]["id"],
                    "atom2_idx": nme_bond["nitrogen"]["id"],
                    "type": "RESIDUE_TO_NME",
                }
            )

            _console.print(
                f"[green]✓[/green] Found C-terminus separator bond: [cyan]{nme_bond['carbon_name']} - {nme_bond['nitrogen_name']}[/cyan]"
            )

        return {
            "success": True,
            "atoms": atoms,
            "bonds": bonds,
            "suggested_sep_bonds": suggested_bonds,
        }

    except Exception as e:
        import traceback

        traceback.print_exc()
        return {"success": False, "error": f"Error analyzing AC file: {str(e)}"}


def identify_terminal_atoms(ac_file, sep_bonds):
    """
    Identify the terminal atoms in the AC file that need standard AMBER charges.

    Parameters:
    -----------
    ac_file : str
        Path to the AC file
    sep_bonds : list
        List of separator bonds between residue and caps

    Returns:
    --------
    dict
        Dictionary containing identified terminal atoms
    """
    try:
        with open(ac_file, "r") as f:
            lines = f.readlines()

        # Parse atoms
        atoms = []
        for line in lines:
            if line.startswith("ATOM"):
                parts = line.split()
                atom_id = int(parts[1])
                atom_name = parts[2]
                atoms.append(
                    {
                        "id": atom_id,
                        "name": atom_name,
                        "index": len(atoms),  # Track position in file
                    }
                )

        # Define cap sizes
        ace_cap_size = 6   # ACE: HH31, CH3, HH32, HH33, C, O
        nme_cap_size = 6   # NME: N, H, CH3, HH31, HH32, HH33

        # Assign regions to atoms
        for i, atom in enumerate(atoms):
            if i < ace_cap_size:
                atom["region"] = "ACE"
            elif i >= len(atoms) - nme_cap_size:
                atom["region"] = "NME"
            else:
                atom["region"] = "RESIDUE"

        # Create a map from atom ID to atom
        atom_by_id = {atom["id"]: atom for atom in atoms}

        # Create a map from atom name to atom for each region
        ace_atoms = {atom["name"]: atom for atom in atoms if atom["region"] == "ACE"}
        residue_atoms = {
            atom["name"]: atom for atom in atoms if atom["region"] == "RESIDUE"
        }
        nme_atoms = {atom["name"]: atom for atom in atoms if atom["region"] == "NME"}

        # Parse bonds and map atom names to atoms
        bonds = []
        for line in lines:
            if line.startswith("BOND"):
                parts = line.split()
                if len(parts) >= 7:
                    bond_id = int(parts[1])
                    atom1_name = parts[5]
                    atom2_name = parts[6]

                    # Try to find these atoms
                    atom1 = None
                    atom2 = None

                    # Look in all regions
                    for region in [ace_atoms, residue_atoms, nme_atoms]:
                        if atom1_name in region:
                            atom1 = region[atom1_name]
                        if atom2_name in region:
                            atom2 = region[atom2_name]

                    # Only add bond if we found both atoms
                    if atom1 and atom2:
                        bonds.append(
                            {
                                "id": bond_id,
                                "atom1": atom1,
                                "atom2": atom2,
                                "atom1_name": atom1_name,
                                "atom2_name": atom2_name,
                            }
                        )

        # Find backbone N and H atoms at N-terminus
        n_terminus_atoms = {}

        # First look through the separator bonds to identify the residue N atom
        for bond in sep_bonds:
            if isinstance(bond, dict) and "atom1_name" in bond and "atom2_name" in bond:
                # For N-terminus (ACE-residue bond), the first atom should be from residue (N)
                if bond["atom1_name"].startswith("N"):
                    n_atom_name = bond["atom1_name"]
                    n_atom = None

                    # Find the N atom in the residue atoms dictionary
                    if n_atom_name in residue_atoms:
                        n_atom = residue_atoms[n_atom_name]
                        n_terminus_atoms["N"] = n_atom_name

                        # Now find an H atom bonded to this N
                        for b in bonds:
                            if (
                                b["atom1"]["id"] == n_atom["id"]
                                and b["atom2_name"].startswith("H")
                            ) or (
                                b["atom2"]["id"] == n_atom["id"]
                                and b["atom1_name"].startswith("H")
                            ):
                                # Found an H atom bonded to N
                                h_atom_name = (
                                    b["atom2_name"]
                                    if b["atom1"]["id"] == n_atom["id"]
                                    else b["atom1_name"]
                                )
                                n_terminus_atoms["H"] = h_atom_name
                                break

        # Find backbone C and O atoms at C-terminus
        c_terminus_atoms = {}

        # Look through the separator bonds to identify the residue C atom
        for bond in sep_bonds:
            if isinstance(bond, dict) and "atom1_name" in bond and "atom2_name" in bond:
                # For C-terminus (residue-NME bond), the first atom should be from residue (C)
                if bond["atom1_name"].startswith("C"):
                    c_atom_name = bond["atom1_name"]
                    c_atom = None

                    # Find the C atom in the residue atoms dictionary
                    if c_atom_name in residue_atoms:
                        c_atom = residue_atoms[c_atom_name]
                        c_terminus_atoms["C"] = c_atom_name

                        # Now find an O atom bonded to this C
                        for b in bonds:
                            if (
                                b["atom1"]["id"] == c_atom["id"]
                                and b["atom2_name"].startswith("O")
                            ) or (
                                b["atom2"]["id"] == c_atom["id"]
                                and b["atom1_name"].startswith("O")
                            ):
                                # Found an O atom bonded to C
                                o_atom_name = (
                                    b["atom2_name"]
                                    if b["atom1"]["id"] == c_atom["id"]
                                    else b["atom1_name"]
                                )
                                c_terminus_atoms["O"] = o_atom_name
                                break

        # If we couldn't find the atoms using the separator bonds, try using residue atoms directly
        if not n_terminus_atoms.get("N") or not n_terminus_atoms.get("H"):
            # Find the first N atom in the residue region
            for atom in atoms:
                if atom["region"] == "RESIDUE" and atom["name"].startswith("N"):
                    n_atom = atom
                    n_terminus_atoms["N"] = n_atom["name"]

                    # Now find an H atom bonded to this N
                    for b in bonds:
                        if (
                            b["atom1"]["id"] == n_atom["id"]
                            and b["atom2_name"].startswith("H")
                        ) or (
                            b["atom2"]["id"] == n_atom["id"]
                            and b["atom1_name"].startswith("H")
                        ):
                            # Found an H atom bonded to N
                            h_atom_name = (
                                b["atom2_name"]
                                if b["atom1"]["id"] == n_atom["id"]
                                else b["atom1_name"]
                            )
                            n_terminus_atoms["H"] = h_atom_name
                            break
                    break

        if not c_terminus_atoms.get("C") or not c_terminus_atoms.get("O"):
            # Find the last C atom in the residue region that has an O bonded to it
            residue_c_atoms = [
                atom
                for atom in atoms
                if atom["region"] == "RESIDUE" and atom["name"].startswith("C")
            ]

            # Sort by index to get the last one first
            residue_c_atoms.sort(key=lambda a: a["index"], reverse=True)

            for c_atom in residue_c_atoms:
                # Check if this C atom has an O bonded to it
                has_o_bond = False
                o_atom_name = None

                for b in bonds:
                    if (
                        b["atom1"]["id"] == c_atom["id"]
                        and b["atom2_name"].startswith("O")
                    ) or (
                        b["atom2"]["id"] == c_atom["id"]
                        and b["atom1_name"].startswith("O")
                    ):
                        # Found an O atom bonded to C
                        has_o_bond = True
                        o_atom_name = (
                            b["atom2_name"]
                            if b["atom1"]["id"] == c_atom["id"]
                            else b["atom1_name"]
                        )
                        break

                if has_o_bond:
                    c_terminus_atoms["C"] = c_atom["name"]
                    c_terminus_atoms["O"] = o_atom_name
                    break

        # Print what we found
        _console.print("[bold]Terminal atom identification results:[/bold]")
        _console.print(f"  N-terminus atoms: [cyan]{n_terminus_atoms}[/cyan]")
        _console.print(f"  C-terminus atoms: [cyan]{c_terminus_atoms}[/cyan]")

        return {
            "success": True,
            "n_terminus": n_terminus_atoms,
            "c_terminus": c_terminus_atoms,
        }

    except Exception as e:
        import traceback

        traceback.print_exc()
        return {
            "success": False,
            "error": f"Error identifying terminal atoms: {str(e)}",
        }


def assign_standard_amber_backbone_charges(atom_identification):
    """
    Assign standard AMBER backbone charges to the identified terminal atoms.

    Parameters:
    -----------
    atom_identification : dict
        Dictionary containing identified terminal atoms

    Returns:
    --------
    list
        List of (atom_name, charge) tuples
    """
    # Standard AMBER backbone charges (ff14SB main-chain values, amino12.lib)
    standard_charges = {
        "N": -0.4157,  # backbone amide N
        "H": 0.2719,  # backbone amide H
        "C": 0.5973,  # backbone carbonyl C
        "O": -0.5679,  # backbone carbonyl O
    }

    atom_charges = []

    # Assign charges to N-terminus atoms
    if "N" in atom_identification["n_terminus"]:
        n_atom = atom_identification["n_terminus"]["N"]
        atom_charges.append((n_atom, standard_charges["N"]))

    if "H" in atom_identification["n_terminus"]:
        h_atom = atom_identification["n_terminus"]["H"]
        atom_charges.append((h_atom, standard_charges["H"]))

    # Assign charges to C-terminus atoms
    if "C" in atom_identification["c_terminus"]:
        c_atom = atom_identification["c_terminus"]["C"]
        atom_charges.append((c_atom, standard_charges["C"]))

    if "O" in atom_identification["c_terminus"]:
        o_atom = atom_identification["c_terminus"]["O"]
        atom_charges.append((o_atom, standard_charges["O"]))

    return atom_charges


def generate_residuegen_input(amino_acid, ac_file, esp_file, net_charge, processor=None):
    """
    Generate an input file for the residuegen program.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    ac_file : str
        Path to the AC file
    esp_file : str
        Path to the ESP file
    net_charge : int
        Net charge of the residue

    Returns:
    --------
    dict
        Dictionary containing the path to the generated input file
    """
    amino_acid = amino_acid.strip().lower()
    resname = amino_acid.upper()

    _console.print(f"\n[bold cyan]Step 8: Generating Residuegen Input File for {resname}[/bold cyan]")

    # Check if AC and ESP files exist
    if not os.path.exists(ac_file):
        return {"success": False, "error": f"AC file not found: {ac_file}"}

    if not os.path.exists(esp_file):
        return {"success": False, "error": f"ESP file not found: {esp_file}"}

    # Count conformations in ESP file
    conf_num = count_conformations_in_esp_file(esp_file)
    if conf_num == 0:
        _console.print(
            "[yellow]⚠ Warning: Could not determine the number of conformations in the ESP file.[/yellow]"
        )
        conf_num = int_prompt_with_context(
            processor,
            "Enter the number of conformations in the ESP file",
            module="Modified AA Parameterizer",
            description="Number of conformations in ESP file",
        )

    # Analyze AC file for separator bonds
    ac_analysis = analyze_ac_file_for_separator_bonds(ac_file)
    sep_bonds = []

    if ac_analysis["success"] and ac_analysis["suggested_sep_bonds"]:
        _console.print("\n[cyan]Suggested separator bonds based on AC file analysis:[/cyan]")
        for i, bond in enumerate(ac_analysis["suggested_sep_bonds"]):
            _console.print(f"  [cyan]{i+1}.[/cyan] {bond['atom1_name']} - {bond['atom2_name']}")
            # Add the bond type indicators for ACE-to-RESIDUE (1 N C) and RESIDUE-to-NME (-1 C N)
            if bond.get("type") == "ACE_TO_RESIDUE":
                sep_bonds.append((bond["atom1_name"], bond["atom2_name"], "1 N C"))
            elif bond.get("type") == "RESIDUE_TO_NME":
                sep_bonds.append((bond["atom1_name"], bond["atom2_name"], "-1 C N"))
            else:
                sep_bonds.append((bond["atom1_name"], bond["atom2_name"], ""))

    # If we couldn't find separator bonds, ask the user
    if not sep_bonds:
        _console.print("\n[yellow]○[/yellow] Could not automatically determine separator bonds.")
        _console.print(
            "[grey50]Please enter separator bonds in the format 'Atom_Name1 Atom_Name2 [bond_type]'[/grey50]"
        )
        _console.print(
            "[grey50]Example: N1 C2 1 N C (where N1 belongs to residue and C2 belongs to ACE cap)[/grey50]"
        )
        _console.print(
            "[grey50]Example: C5 N2 -1 C N (where C5 belongs to residue and N2 belongs to NME cap)[/grey50]"
        )

        sep_bond1 = prompt_with_context(
            processor,
            "Enter first separator bond",
            module="Modified AA Parameterizer",
            description="First separator bond (Atom1 Atom2 [bond_type])",
        )
        if sep_bond1:
            parts = sep_bond1.split()
            if len(parts) >= 2:
                atom1 = parts[0]
                atom2 = parts[1]
                # Check if bond type is specified
                bond_type = (
                    " ".join(parts[2:]) if len(parts) > 2 else "1 N C"
                )  # Default for first bond
                sep_bonds.append((atom1, atom2, bond_type))

        sep_bond2 = prompt_with_context(
            processor,
            "Enter second separator bond (optional)",
            default="",
            module="Modified AA Parameterizer",
            description="Second separator bond (optional)",
        )
        if sep_bond2:
            parts = sep_bond2.split()
            if len(parts) >= 2:
                atom1 = parts[0]
                atom2 = parts[1]
                # Check if bond type is specified
                bond_type = (
                    " ".join(parts[2:]) if len(parts) > 2 else "-1 C N"
                )  # Default for second bond
                sep_bonds.append((atom1, atom2, bond_type))
    else:
        # Confirm with user
        _console.print("\n[cyan]Do these separator bonds look correct?[/cyan]")
        for i, (atom1, atom2, bond_type) in enumerate(sep_bonds):
            _console.print(f"  [cyan]{i+1}.[/cyan] {atom1} - {atom2} {bond_type}")

        confirm = confirm_with_context(
            processor,
            "Use these bonds?",
            default=True,
            module="Modified AA Parameterizer",
            description="Use detected separator bonds",
        )

        if not confirm:
            # User wants to enter custom bonds
            sep_bonds = []
            _console.print(
                "\n[grey50]Please enter separator bonds in the format 'Atom_Name1 Atom_Name2 [bond_type]'[/grey50]"
            )
            _console.print(
                "[grey50]Example: N1 C2 1 N C (where N1 belongs to residue and C2 belongs to ACE cap)[/grey50]"
            )
            _console.print(
                "[grey50]Example: C5 N2 -1 C N (where C5 belongs to residue and N2 belongs to NME cap)[/grey50]"
            )

            sep_bond1 = prompt_with_context(
            processor,
            "Enter first separator bond",
            module="Modified AA Parameterizer",
            description="First separator bond (Atom1 Atom2 [bond_type])",
        )
            if sep_bond1:
                parts = sep_bond1.split()
                if len(parts) >= 2:
                    atom1 = parts[0]
                    atom2 = parts[1]
                    # Check if bond type is specified
                    bond_type = (
                        " ".join(parts[2:]) if len(parts) > 2 else "1 N C"
                    )  # Default for first bond
                    sep_bonds.append((atom1, atom2, bond_type))

            sep_bond2 = prompt_with_context(
            processor,
            "Enter second separator bond (optional)",
            default="",
            module="Modified AA Parameterizer",
            description="Second separator bond (optional)",
        )
            if sep_bond2:
                parts = sep_bond2.split()
                if len(parts) >= 2:
                    atom1 = parts[0]
                    atom2 = parts[1]
                    # Check if bond type is specified
                    bond_type = (
                        " ".join(parts[2:]) if len(parts) > 2 else "-1 C N"
                    )  # Default for second bond
                    sep_bonds.append((atom1, atom2, bond_type))

    # Identify terminal atoms and assign standard AMBER charges
    # Convert sep_bonds format for compatibility with identify_terminal_atoms
    bonds_for_terminal = [
        {"atom1_name": atom1, "atom2_name": atom2} for atom1, atom2, _ in sep_bonds
    ]
    terminal_atoms = identify_terminal_atoms(ac_file, bonds_for_terminal)

    atom_charges = []
    if terminal_atoms["success"]:
        atom_charges = assign_standard_amber_backbone_charges(terminal_atoms)

        _console.print("\n[cyan]Assigned standard AMBER backbone charges:[/cyan]")
        for atom_name, charge in atom_charges:
            _console.print(f"  {atom_name}: {charge:.4f}")

        # Inform user about standard charges
        _console.print(Panel(
            "[bold]Standard AMBER backbone charges (ff14SB):[/bold]\n"
            "  backbone N atom: -0.4157\n"
            "  backbone H atom: 0.2719\n"
            "  backbone C atom: 0.5973\n"
            "  backbone O atom: -0.5679",
            title="Note",
            border_style="grey50",
            expand=False
        ))

        confirm = confirm_with_context(
            processor,
            "\nUse these standard charges?",
            default=True,
            module="Modified AA Parameterizer",
            description="Use standard AMBER backbone charges",
        )

        if not confirm:
            # User wants to enter custom charges
            atom_charges = []
            _console.print("\n[grey50]Enter atom charges in the format 'Atom_Name Charge'[/grey50]")
            _console.print("[grey50]Enter a blank line when done[/grey50]")

            while True:
                charge_input = prompt_with_context(
                    processor,
                    "Enter atom charge (or blank to finish)",
                    default="",
                    module="Modified AA Parameterizer",
                    description="Custom atom charge (Atom_Name Charge, or blank to finish)",
                )
                if not charge_input:
                    break

                parts = charge_input.split()
                if len(parts) >= 2:
                    try:
                        atom_name = parts[0]
                        charge = float(parts[1])
                        atom_charges.append((atom_name, charge))
                    except ValueError:
                        _console.print("[red]Invalid charge value. Please enter a number.[/red]")
    else:
        # If we couldn't identify terminal atoms, ask the user
        _console.print("\n[yellow]○[/yellow] Could not automatically identify terminal atoms for standard charges.")
        _console.print("[grey50]Please enter the four key backbone atoms that need standard charges:[/grey50]")

        # Ask for N-terminus atoms
        n_atom = prompt_with_context(
            processor,
            "Enter N-terminus N atom name",
            module="Modified AA Parameterizer",
            description="N-terminus N atom name",
        )
        h_atom = prompt_with_context(
            processor,
            "Enter N-terminus H atom name",
            module="Modified AA Parameterizer",
            description="N-terminus H atom name",
        )

        # Ask for C-terminus atoms
        c_atom = prompt_with_context(
            processor,
            "Enter C-terminus C atom name",
            module="Modified AA Parameterizer",
            description="C-terminus C atom name",
        )
        o_atom = prompt_with_context(
            processor,
            "Enter C-terminus O atom name",
            module="Modified AA Parameterizer",
            description="C-terminus O atom name",
        )

        # Assign standard charges
        if n_atom:
            atom_charges.append((n_atom, -0.4157))
        if h_atom:
            atom_charges.append((h_atom, 0.2719))
        if c_atom:
            atom_charges.append((c_atom, 0.5973))
        if o_atom:
            atom_charges.append((o_atom, -0.5679))

    # residuegen output settings (the defaults are correct in almost all cases;
    # these prompts expose residuegen's three output fields).
    prep_file = f"{amino_acid.lower()}.prep"
    prep_file = prompt_with_context(
        processor,
        "Output prep file to write (AMBER .prep)",
        default=prep_file,
        module="Modified AA Parameterizer",
        description="Filename for the AMBER .prep that residuegen will write",
    )

    residue_file_name = f"{amino_acid.lower()}.res"
    residue_file_name = prompt_with_context(
        processor,
        "Internal residue-file label inside the prep",
        default=residue_file_name,
        module="Modified AA Parameterizer",
        description="residuegen RESIDUE_FILE_NAME field — a label recorded in the prep header, not a file you use",
    )

    residue_symbol = resname
    residue_symbol = prompt_with_context(
        processor,
        "3-letter residue code for this residue",
        default=residue_symbol,
        module="Modified AA Parameterizer",
        description="residuegen RESIDUE_SYMBOL — the residue name (e.g. CYR) that tLEaP/loadpdb will match",
    )

    # Generate the residuegen input file
    input_file = f"{amino_acid.lower()}_residuegen.in"

    with open(input_file, "w") as f:
        f.write("#residuegen input file\n\n")

        f.write(
            "#INPUT_FILE:    structure file in ac format, generated from a Gaussian output with 'antechamber'\n"
        )
        f.write(f"INPUT_FILE      {ac_file}\n\n")

        f.write("#CONF_NUM:  Number of conformations applied\n")
        f.write(f"CONF_NUM        {conf_num}\n\n")

        f.write(
            "#ESP_FILE:      esp file generated from gaussian output with 'espgen'\n"
        )
        f.write(f"ESP_FILE        {esp_file}\n\n")

        f.write(
            "#SEP_BOND:  bonds that separate residue and caps, input in a format of (Atom_Name1 Atom_Name2),\n"
        )
        f.write(
            "#       where Atom_Name1 belongs to the residue and Atom_Name2 belongs to a cap.\n"
        )
        for atom1, atom2, bond_type in sep_bonds:
            f.write(f"SEP_BOND        {atom1} {atom2} {bond_type}\n")
        f.write("\n")

        f.write("#NET_CHARGE:    net charge of the residue\n")
        f.write(f"NET_CHARGE      {net_charge}\n\n")

        if atom_charges:
            f.write(
                "#ATOM_CHARGE:   predefined atom charge, input in a format of (Atom_Name Partial_Charge)\n"
            )
            f.write(
                "#               Standard AMBER backbone charges applied to terminal atoms\n"
            )
            for atom_name, charge in atom_charges:
                f.write(f"ATOM_CHARGE     {atom_name} {charge:.4f}\n")
            f.write("\n")

        f.write("#PREP_FILE:     prep file name\n")
        f.write(f"PREP_FILE:      {prep_file}\n\n")

        f.write("#RESIDUE_FILE_NAME:     residue file name in PREP_FILE\n")
        f.write(f"RESIDUE_FILE_NAME:      {residue_file_name}\n\n")

        f.write("#RESIDUE_SYMBOL:    residue symbol in PREP_FILE\n")
        f.write(f"RESIDUE_SYMBOL:     {residue_symbol}\n")

    _console.print(f"\n[green]✓ Successfully generated residuegen input file: {input_file}[/green]")

    return {
        "success": True,
        "input_file": input_file,
        "prep_file": prep_file,
        "residue_symbol": residue_symbol,
    }


def run_residuegen(input_file):
    """
    Run the residuegen program with the provided input file.

    Parameters:
    -----------
    input_file : str
        Path to the residuegen input file

    Returns:
    --------
    dict
        Dictionary containing the result of the residuegen command
    """
    try:
        _console.print(f"\n[cyan]→[/cyan] Running residuegen with input file: {input_file}")

        cmd = ["residuegen", input_file]

        process = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )

        # Check if prep file was created
        input_dir = os.path.dirname(input_file) or "."

        # Extract prep file name from input file
        prep_file = None
        with open(input_file, "r") as f:
            for line in f:
                if line.startswith("PREP_FILE:"):
                    prep_file = line.split()[1].strip()
                    break

        if prep_file and os.path.exists(os.path.join(input_dir, prep_file)):
            return {
                "success": True,
                "message": f"Successfully created {prep_file}",
                "output": process.stdout,
                "prep_file": prep_file,
            }
        else:
            # Try to find any prep files that might have been created
            prep_files = glob.glob(os.path.join(input_dir, "*.prep"))

            if prep_files:
                return {
                    "success": True,
                    "message": f"Successfully created prep file(s): {', '.join(os.path.basename(f) for f in prep_files)}",
                    "output": process.stdout,
                    "prep_file": prep_files[0],  # Return the first one found
                }
            else:
                return {
                    "success": False,
                    "error": f"residuegen did not create the expected prep file",
                    "output": process.stdout,
                    "stderr": process.stderr,
                }

    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "error": f"Error running residuegen: {e.stderr}",
            "output": e.stdout,
        }
    except FileNotFoundError:
        return {
            "success": False,
            "error": "residuegen command not found. Make sure AmberTools is installed and in your PATH.",
        }


def generate_and_run_residuegen(amino_acid, ac_file, esp_file, net_charge, processor=None):
    """
    Generate a residuegen input file and run the residuegen program.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    ac_file : str
        Path to the AC file
    esp_file : str
        Path to the ESP file
    net_charge : int
        Net charge of the residue

    Returns:
    --------
    dict
        Dictionary containing the results of the residuegen process
    """
    _console.print(f"\n[bold cyan]Step 8: Residuegen Process for {amino_acid.upper()}[/bold cyan]")

    # Generate the residuegen input file
    input_result = generate_residuegen_input(amino_acid, ac_file, esp_file, net_charge, processor=processor)

    if not input_result["success"]:
        _console.print(
            f"\n[red]✗ Error generating residuegen input file: {input_result.get('error', 'Unknown error')}[/red]"
        )
        return input_result

    # Ask if the user wants to run residuegen now
    run_now = confirm_with_context(
        processor,
        "\nDo you want to run residuegen now?",
        default=True,
        module="Modified AA Parameterizer",
        description="Run residuegen now",
    )

    if run_now:
        # Run residuegen
        residuegen_result = run_residuegen(input_result["input_file"])

        if residuegen_result["success"]:
            _console.print(f"\n[green]✓ {residuegen_result['message']}[/green]")

            # Return the complete result
            return {
                "success": True,
                "message": f"Successfully completed residuegen process",
                "input_file": input_result["input_file"],
                "prep_file": residuegen_result.get("prep_file"),
                "residue_symbol": input_result["residue_symbol"],
            }
        else:
            _console.print(
                f"\n[red]✗ Error running residuegen: {residuegen_result.get('error', 'Unknown error')}[/red]"
            )
            return residuegen_result
    else:
        _console.print(f"\n[yellow]○[/yellow] Skipping residuegen execution. You can run it manually with:")
        _console.print(f"  [cyan]residuegen -i {input_result['input_file']}[/cyan]")

        return {
            "success": True,
            "message": f"Generated residuegen input file but skipped execution",
            "input_file": input_result["input_file"],
            "residue_symbol": input_result["residue_symbol"],
        }


def extract_scf_energy(log_file):
    """
    Extract the last SCF energy value from a Gaussian log file.

    Parameters:
    -----------
    log_file : str
        Path to the Gaussian log file

    Returns:
    --------
    float or None
        The SCF energy value in atomic units, or None if not found
    """
    try:
        with open(log_file, "r") as f:
            content = f.read()

        # Find all "SCF Done" lines
        scf_matches = re.findall(
            r"SCF Done:\s+E\([^)]+\)\s*=\s*([+-]?\d+\.\d+)", content
        )

        if scf_matches:
            # Get the last (most recent) SCF energy value
            last_energy = float(scf_matches[-1])
            return last_energy
        else:
            _console.print(f"[yellow]⚠[/yellow] No SCF energy found in {log_file}")
            return None

    except Exception as e:
        _console.print(f"[red]✗ Error reading energy from {log_file}: {str(e)}[/red]")
        return None


def find_lowest_energy_structure(log_files):
    """
    Find the log file with the lowest SCF energy among a list of log files.

    Parameters:
    -----------
    log_files : list
        List of Gaussian log files to check

    Returns:
    --------
    dict
        Dictionary containing the path to the lowest energy log file and its energy
    """
    if not log_files:
        return {"success": False, "error": "No log files provided"}

    _console.print(f"\n[cyan]→[/cyan] Analyzing energies of {len(log_files)} log files...")

    lowest_energy = None
    lowest_energy_file = None
    energies = {}

    for log_file in log_files:
        energy = extract_scf_energy(log_file)

        if energy is not None:
            energies[log_file] = energy

            # Check if this is the lowest energy so far
            if lowest_energy is None or energy < lowest_energy:
                lowest_energy = energy
                lowest_energy_file = log_file

    if lowest_energy_file:
        _console.print(f"\n[green]✓ Lowest energy structure:[/green]")
        _console.print(f"  File: {os.path.basename(lowest_energy_file)}")
        _console.print(f"  Energy: {lowest_energy} A.U.")

        # Also print a few other energies for comparison
        sorted_energies = sorted(energies.items(), key=lambda x: x[1])

        if len(sorted_energies) > 1:
            _console.print("\n[cyan]Structure energies (lowest to highest):[/cyan]")
            for i, (file, energy) in enumerate(
                sorted_energies[:5]
            ):  # Show top 5 lowest energies
                _console.print(f"  [cyan]{i+1}.[/cyan] {os.path.basename(file)}: {energy} A.U.")

            # Show energy range
            if len(sorted_energies) > 5:
                _console.print(f"  [grey50]... {len(sorted_energies)-5} more structures ...[/grey50]")

            highest_file, highest_energy = sorted_energies[-1]
            _console.print(f"  Highest: {os.path.basename(highest_file)}: {highest_energy} A.U.")
            _console.print(f"  Energy range: {highest_energy - lowest_energy} A.U.")

        return {
            "success": True,
            "lowest_energy_file": lowest_energy_file,
            "lowest_energy": lowest_energy,
            "energies": energies,
        }
    else:
        return {
            "success": False,
            "error": "Could not find energy values in any of the log files",
        }


def run_antechamber(log_file, output_ac, resname, net_charge):
    """
    Run antechamber to generate an ac file from a Gaussian log file.

    Parameters:
    -----------
    log_file : str
        Path to the Gaussian log file
    output_ac : str
        Path to the output ac file
    resname : str
        Custom residue name
    net_charge : int
        Net charge of the molecule

    Returns:
    --------
    dict
        Dictionary containing the result of the antechamber command
    """
    try:
        _console.print(Panel(
            f"[bold]Command:[/bold] antechamber\n\n"
            f"[bold]Parameters:[/bold]\n"
            f"  [cyan]Input:[/cyan]  {os.path.basename(log_file)}\n"
            f"  [cyan]Output:[/cyan] {output_ac}\n"
            f"  [cyan]Residue name:[/cyan] {resname}\n"
            f"  [cyan]Net charge:[/cyan] {net_charge}",
            title="Running antechamber",
            border_style="grey50",
            expand=False
        ))

        cmd = [
            "antechamber",
            "-fi",
            "gout",  # Input format is Gaussian output
            "-fo",
            "ac",  # Output format is AC
            "-i",
            log_file,  # Input file
            "-o",
            output_ac,  # Output file
            "-c",
            "resp",  # Charge method is RESP
            "-rn",
            resname,  # Residue name
            "-at",
            "amber",  # Atom types are AMBER
            "-s",
            "2",  # Status: 2 (for verbose output)
            "-nc",
            str(net_charge),  # Net charge
        ]

        process = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )

        # Check if output file was created
        if os.path.exists(output_ac):
            return {
                "success": True,
                "message": f"Successfully created {output_ac}",
                "output": process.stdout,
            }
        else:
            return {
                "success": False,
                "error": f"antechamber did not create the output file {output_ac}",
                "output": process.stdout,
                "stderr": process.stderr,
            }

    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "error": f"Error running antechamber: {e.stderr}",
            "output": e.stdout,
        }
    except FileNotFoundError:
        return {
            "success": False,
            "error": "antechamber command not found. Make sure AmberTools is installed and in your PATH.",
        }


def generate_ac_file(amino_acid, processor=None):
    """
    Find the lowest energy structure and generate an AC file using antechamber.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code

    Returns:
    --------
    dict
        Dictionary containing the result of the AC file generation
    """
    amino_acid = amino_acid.strip().lower()

    _console.print(f"\n[bold cyan]Step 7: Generating AC File for {amino_acid.upper()}[/bold cyan]")

    # Get list of all log files from extracted structure directories
    structure_check = check_extracted_structure_calculations(amino_acid)

    all_log_files = []
    for conf in ["ahelix", "bsheet"]:
        if structure_check[conf]["exists"]:
            all_log_files.extend(structure_check[conf]["log_files"])

    if not all_log_files:
        _console.print("\n[red]✗ No log files found for extracted structures. Cannot generate AC file.[/red]")
        return {
            "success": False,
            "error": "No log files found for extracted structures",
        }

    # Find the lowest energy structure
    energy_result = find_lowest_energy_structure(all_log_files)

    if not energy_result["success"]:
        _console.print(
            f"\n[red]✗ Error finding lowest energy structure: {energy_result.get('error', 'Unknown error')}[/red]"
        )
        return energy_result

    # Get the custom residue name and net charge from the user
    _console.print("\n[cyan]To generate the AC file, please provide the following information:[/cyan]")

    # Default residue name is the uppercase amino acid code
    default_resname = amino_acid.upper()
    custom_resname = prompt_with_context(
        processor,
        "Enter custom residue name",
        default=default_resname,
        module="Modified AA Parameterizer",
        description="Custom residue name for AC file",
    )

    # Default net charge is 0
    try:
        net_charge = int_prompt_with_context(
            processor,
            "Enter the net charge of the molecule",
            default=0,
            module="Modified AA Parameterizer",
            description="Net charge of molecule for AC file",
        )
    except ValueError:
        _console.print("[yellow]⚠ Invalid charge value. Using default value of 0.[/yellow]")
        net_charge = 0

    # Generate the output AC filename
    output_ac = f"{custom_resname}.ac"

    # Run antechamber
    antechamber_result = run_antechamber(
        energy_result["lowest_energy_file"], output_ac, custom_resname, net_charge
    )

    if antechamber_result["success"]:
        _console.print(f"\n[green]✓ Successfully generated AC file: {output_ac}[/green]")
        _console.print(
            f"[grey50]This file can be used with parmchk2 to generate FRCMOD file for the custom residue.[/grey50]"
        )

        # Return the complete result
        return {
            "success": True,
            "message": f"Successfully generated AC file: {output_ac}",
            "ac_file": output_ac,
            "residue_name": custom_resname,
            "charge": net_charge,
            "source_log": energy_result["lowest_energy_file"],
            "energy": energy_result["lowest_energy"],
        }
    else:
        _console.print(
            f"\n[red]✗ Error generating AC file: {antechamber_result.get('error', 'Unknown error')}[/red]"
        )
        return antechamber_result


def check_extracted_structure_calculations(amino_acid):
    """
    Check if Gaussian calculations for the extracted structures are complete.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code

    Returns:
    --------
    dict
        A dictionary containing information about the status of the calculations
    """
    amino_acid = amino_acid.strip().lower()

    # Define the structure directories
    ahelix_dir = f"{amino_acid}_ahelix_structures"
    bsheet_dir = f"{amino_acid}_bsheet_structures"

    result = {
        "success": True,
        "ahelix": {
            "exists": False,
            "gjf_files": [],
            "log_files": [],
            "missing_log_files": [],
            "extra_log_files": [],
        },
        "bsheet": {
            "exists": False,
            "gjf_files": [],
            "log_files": [],
            "missing_log_files": [],
            "extra_log_files": [],
        },
    }

    # Check alpha helix directory
    if os.path.exists(ahelix_dir):
        result["ahelix"]["exists"] = True

        # Get gjf and log files
        gjf_files = sorted(glob.glob(os.path.join(ahelix_dir, "*.gjf")))
        log_files = sorted(glob.glob(os.path.join(ahelix_dir, "*.log")))

        result["ahelix"]["gjf_files"] = gjf_files
        result["ahelix"]["log_files"] = log_files

        # Identify missing log files
        for gjf_file in gjf_files:
            expected_log = os.path.splitext(gjf_file)[0] + ".log"
            if expected_log not in log_files:
                result["ahelix"]["missing_log_files"].append(gjf_file)

        # Identify unexpected log files
        for log_file in log_files:
            expected_gjf = os.path.splitext(log_file)[0] + ".gjf"
            if expected_gjf not in gjf_files:
                result["ahelix"]["extra_log_files"].append(log_file)

    # Check beta sheet directory
    if os.path.exists(bsheet_dir):
        result["bsheet"]["exists"] = True

        # Get gjf and log files
        gjf_files = sorted(glob.glob(os.path.join(bsheet_dir, "*.gjf")))
        log_files = sorted(glob.glob(os.path.join(bsheet_dir, "*.log")))

        result["bsheet"]["gjf_files"] = gjf_files
        result["bsheet"]["log_files"] = log_files

        # Identify missing log files
        for gjf_file in gjf_files:
            expected_log = os.path.splitext(gjf_file)[0] + ".log"
            if expected_log not in log_files:
                result["bsheet"]["missing_log_files"].append(gjf_file)

        # Identify unexpected log files
        for log_file in log_files:
            expected_gjf = os.path.splitext(log_file)[0] + ".gjf"
            if expected_gjf not in gjf_files:
                result["bsheet"]["extra_log_files"].append(log_file)

    return result


def run_espgen_on_log_files(log_files):
    """
    Run the espgen program on each log file to generate ESP files.

    Parameters:
    -----------
    log_files : list
        List of Gaussian log files to process

    Returns:
    --------
    dict
        A dictionary containing results of the ESP generation
    """
    result = {"success": True, "esp_files": [], "failed_files": []}

    if not log_files:
        return {"success": False, "error": "No log files provided to process"}

    _console.print(f"\n[cyan]→[/cyan] Running espgen on {len(log_files)} log files...")

    for i, log_file in enumerate(log_files):
        _console.print(f"[{i+1}/{len(log_files)}] Processing {os.path.basename(log_file)}...")

        # Generate output ESP filename
        esp_file = os.path.splitext(log_file)[0] + ".esp"

        try:
            # Run espgen command
            cmd = ["espgen", "-i", log_file, "-o", esp_file]
            process = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

            # Check if the ESP file was created
            if os.path.exists(esp_file):
                result["esp_files"].append(esp_file)
                _console.print(f"  [green]✓[/green] Generated: {os.path.basename(esp_file)}")
            else:
                result["failed_files"].append(log_file)
                _console.print(f"  [red]✗[/red] Error: ESP file not created for {os.path.basename(log_file)}")

        except subprocess.CalledProcessError as e:
            _console.print(f"  [red]✗[/red] Error running espgen: {e.stderr}")
            result["failed_files"].append(log_file)
        except FileNotFoundError:
            _console.print(
                "[red]✗ Error: espgen command not found. Make sure AmberTools is installed and in your PATH.[/red]"
            )
            return {
                "success": False,
                "error": "espgen command not found. Make sure AmberTools is installed and in your PATH.",
            }

    # Update success status
    if result["failed_files"]:
        result["success"] = False
        result["error"] = (
            f"Failed to generate ESP files for {len(result['failed_files'])} log files"
        )

    return result


def concatenate_esp_files(esp_files, output_file):
    """
    Concatenate multiple ESP files into a single file.

    Parameters:
    -----------
    esp_files : list
        List of ESP files to concatenate
    output_file : str
        Path to the output concatenated file

    Returns:
    --------
    dict
        A dictionary containing results of the concatenation
    """
    if not esp_files:
        return {"success": False, "error": "No ESP files to concatenate"}

    _console.print(f"\n[cyan]→[/cyan] Concatenating {len(esp_files)} ESP files into {output_file}...")

    try:
        # Open output file for writing
        with open(output_file, "w") as outfile:
            # Process each ESP file
            for i, esp_file in enumerate(esp_files):
                _console.print(f"  [green]✓[/green] Adding: {os.path.basename(esp_file)}")

                with open(esp_file, "r") as infile:
                    # Read the content of the ESP file
                    content = infile.read()

                    # Write to the output file
                    #                   outfile.write(f"# File {i+1}: {os.path.basename(esp_file)}\n")
                    outfile.write(content)
        #                   outfile.write("\n\n")

        return {
            "success": True,
            "message": f"Successfully concatenated {len(esp_files)} ESP files to {output_file}",
        }

    except Exception as e:
        return {"success": False, "error": f"Error concatenating ESP files: {str(e)}"}


def process_extracted_structures_esp(amino_acid, processor=None):
    """
    Process the log files from extracted structures to generate and concatenate ESP files.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code

    Returns:
    --------
    dict
        A dictionary containing results of the processing
    """
    amino_acid = amino_acid.strip().lower()

    _console.print(
        f"\n[bold cyan]Step 6: Processing ESP Data for {amino_acid.upper()} Extracted Structures[/bold cyan]"
    )

    # Check status of extracted structure calculations
    check_result = check_extracted_structure_calculations(amino_acid)

    # Define the structure directories
    ahelix_dir = f"{amino_acid}_ahelix_structures"
    bsheet_dir = f"{amino_acid}_bsheet_structures"

    # Report status of calculations
    _console.print("\n[cyan]Status of Gaussian calculations for extracted structures:[/cyan]")

    missing_calculations = 0
    has_complete_calculations = False

    for conf in ["ahelix", "bsheet"]:
        if check_result[conf]["exists"]:
            gjf_count = len(check_result[conf]["gjf_files"])
            log_count = len(check_result[conf]["log_files"])
            missing_count = len(check_result[conf]["missing_log_files"])
            extra_count = len(check_result[conf]["extra_log_files"])

            _console.print(f"\n[bold]{conf.upper()}[/bold] conformation ({amino_acid}_{conf}_structures):")
            _console.print(f"  Found {gjf_count} input files (.gjf)")
            _console.print(f"  Found {log_count} output files (.log)")

            if gjf_count == 0:
                _console.print("  [yellow]○[/yellow] No structure files found. Skipping this conformation.")
                continue

            if missing_count > 0:
                _console.print(f"  [yellow]⚠[/yellow] Missing {missing_count} log files for existing input files")
                missing_calculations += missing_count
            else:
                _console.print("  [green]✓[/green] All input files have corresponding output files")
                has_complete_calculations = True

            if extra_count > 0:
                _console.print(
                    f"  [yellow]⚠ Warning: Found {extra_count} log files without corresponding input files[/yellow]"
                )
        else:
            _console.print(
                f"\n[bold]{conf.upper()}[/bold] conformation: [grey50]Directory not found ({amino_acid}_{conf}_structures)[/grey50]"
            )

    # If there are missing calculations, offer to run them
    if missing_calculations > 0:
        run_missing = confirm_with_context(
            processor,
            f"\nThere are {missing_calculations} missing Gaussian calculations. Would you like to run them now?",
            default=False,
            module="Modified AA Parameterizer",
            description="Run missing Gaussian calculations now",
        )

        if run_missing:
            # Collect all missing calculation files
            missing_files = []
            for conf in ["ahelix", "bsheet"]:
                if check_result[conf]["exists"]:
                    missing_files.extend(check_result[conf]["missing_log_files"])

            # Run batch calculations for missing files
            if missing_files:
                # Find Gaussian executable
                gaussian_exe = find_gaussian_executable()
                if not gaussian_exe:
                    _console.print(
                        "[yellow]○[/yellow] Could not find Gaussian executable. Please specify the path manually."
                    )
                    gaussian_exe = prompt_with_context(
                        processor,
                        "Enter path to Gaussian executable (g16, g09, etc.)",
                        module="Modified AA Parameterizer",
                        description="Path to Gaussian executable",
                    )

                    # Verify executable exists
                    if not os.path.exists(gaussian_exe):
                        return {
                            "success": False,
                            "error": f"Gaussian executable not found: {gaussian_exe}",
                        }

                # Ask for maximum concurrent calculations
                try:
                    max_concurrent = int_prompt_with_context(
                        processor,
                        "Enter maximum number of concurrent Gaussian calculations",
                        default=1,
                        module="Modified AA Parameterizer",
                        description="Maximum concurrent Gaussian calculations",
                    )
                    if max_concurrent < 1:
                        max_concurrent = 1
                except ValueError:
                    _console.print("[yellow]⚠ Invalid input. Using default value of 1.[/yellow]")
                    max_concurrent = 1

                # Run in batch mode
                batch_results = batch_run_gaussian_files(
                    missing_files, max_concurrent, gaussian_exe
                )

                # Update check results after batch run
                check_result = check_extracted_structure_calculations(amino_acid)

                # Check if we now have complete calculations
                for conf in ["ahelix", "bsheet"]:
                    if (
                        check_result[conf]["exists"]
                        and len(check_result[conf]["missing_log_files"]) == 0
                    ):
                        has_complete_calculations = True

    # If we still don't have complete calculations, ask if user wants to continue anyway
    if not has_complete_calculations:
        _console.print("\n[yellow]⚠ Warning: Not all Gaussian calculations are complete.[/yellow]")
        cont = confirm_with_context(
            processor,
            "Would you like to continue with the available calculations?",
            default=False,
            module="Modified AA Parameterizer",
            description="Continue with available calculations",
        )

        if not cont:
            _console.print(
                "[yellow]○[/yellow] Exiting. Please run the remaining Gaussian calculations and restart the program."
            )
            return {"success": False, "message": "Incomplete calculations"}

    # Collect all available log files
    all_log_files = []
    for conf in ["ahelix", "bsheet"]:
        if check_result[conf]["exists"]:
            all_log_files.extend(check_result[conf]["log_files"])

    if not all_log_files:
        _console.print("[red]✗[/red] No log files found. Cannot proceed with ESP generation.")
        return {"success": False, "error": "No log files found"}

    # Run espgen on all log files
    espgen_result = run_espgen_on_log_files(all_log_files)

    if not espgen_result["success"]:
        _console.print(
            f"[red]✗[/red] Error generating ESP files: {espgen_result.get('error', 'Unknown error')}"
        )
        return espgen_result

    # Concatenate ESP files
    if espgen_result["esp_files"]:
        # Define output file for concatenated ESP
        output_file = f"{amino_acid}_combined.esp"

        # Concatenate ESP files
        concat_result = concatenate_esp_files(espgen_result["esp_files"], output_file)

        if not concat_result["success"]:
            _console.print(
                f"[red]✗[/red] Error concatenating ESP files: {concat_result.get('error', 'Unknown error')}"
            )
            return concat_result

        _console.print(f"\n[green]✓[/green] ESP processing complete!")
        _console.print(f"  Generated {len(espgen_result['esp_files'])} individual ESP files")
        _console.print(f"  Created combined ESP file: [cyan]{output_file}[/cyan]")

        return {
            "success": True,
            "message": f"Successfully processed {len(all_log_files)} log files and created {output_file}",
            "esp_files": espgen_result["esp_files"],
            "combined_file": output_file,
        }
    else:
        _console.print("[red]✗[/red] No ESP files were generated. Cannot create combined file.")
        return {"success": False, "error": "No ESP files generated"}


def run_gaussian_calculation(input_file, output_file, gaussian_exe="g16"):
    """
    Run a single Gaussian calculation.

    Parameters:
    -----------
    input_file : str
        Path to the Gaussian input file
    output_file : str
        Path to the output log file
    gaussian_exe : str
        Gaussian executable (default: 'g16')

    Returns:
    --------
    dict
        A dictionary containing success status and other information
    """
    try:
        _console.print(f"[cyan]→[/cyan] Starting Gaussian calculation for [cyan]{input_file}[/cyan]...")
        start_time = time.time()

        # Run Gaussian with input file redirected to output file
        process = subprocess.run(
            [gaussian_exe, input_file],
            stdout=open(output_file, "w"),
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

        end_time = time.time()
        elapsed_time = end_time - start_time

        _console.print(
            f"[green]✓[/green] Completed Gaussian calculation for [cyan]{input_file}[/cyan] in {elapsed_time:.1f} seconds"
        )
        return {
            "success": True,
            "input_file": input_file,
            "output_file": output_file,
            "elapsed_time": elapsed_time,
        }
    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "input_file": input_file,
            "output_file": output_file,
            "error": f"Gaussian calculation failed: {e.stderr}",
        }
    except Exception as e:
        return {
            "success": False,
            "input_file": input_file,
            "output_file": output_file,
            "error": f"Error during Gaussian calculation: {str(e)}",
        }


def worker_thread(job_queue, result_queue, gaussian_exe):
    """
    Worker thread function to process Gaussian jobs from a queue.

    Parameters:
    -----------
    job_queue : Queue
        Queue containing input and output file pairs
    result_queue : Queue
        Queue to store results of calculations
    gaussian_exe : str
        Path to Gaussian executable
    """
    while True:
        try:
            # Get job from queue with a timeout
            job = job_queue.get(timeout=5)
            if job is None:  # Sentinel value to exit
                job_queue.task_done()
                break

            input_file, output_file = job
            result = run_gaussian_calculation(input_file, output_file, gaussian_exe)
            result_queue.put(result)
            job_queue.task_done()
        except queue.Empty:
            # If queue is empty for 5 seconds, check if we should exit
            continue
        except Exception as e:
            _console.print(f"[red]✗[/red] Worker thread error: {str(e)}")
            # Put error result in result queue
            result_queue.put(
                {
                    "success": False,
                    "input_file": "unknown",
                    "output_file": "unknown",
                    "error": f"Worker thread error: {str(e)}",
                }
            )
            job_queue.task_done()


def batch_run_gaussian_files(input_files, max_concurrent, gaussian_exe="g16"):
    """
    Run Gaussian calculations in batches with specified concurrency.

    Parameters:
    -----------
    input_files : list
        List of Gaussian input files to process
    max_concurrent : int
        Maximum number of concurrent calculations
    gaussian_exe : str
        Path to Gaussian executable

    Returns:
    --------
    dict
        A dictionary containing success status and results
    """
    results = {"total": len(input_files), "completed": 0, "failed": 0, "results": []}

    if not input_files:
        _console.print("[yellow]○[/yellow] No input files to process.")
        return results

    _console.print(f"\n[bold]Running Gaussian calculations in batch mode[/bold]")
    _console.print(f"  Total files to process: [cyan]{len(input_files)}[/cyan]")
    _console.print(f"  Maximum concurrent calculations: [cyan]{max_concurrent}[/cyan]")
    _console.print(f"  Using Gaussian executable: [cyan]{gaussian_exe}[/cyan]")

    # Create job queue and result queue
    job_queue = queue.Queue()
    result_queue = queue.Queue()

    # Create worker threads
    workers = []
    for i in range(max_concurrent):
        worker = threading.Thread(
            target=worker_thread,
            args=(job_queue, result_queue, gaussian_exe),
            daemon=True,
        )
        workers.append(worker)
        worker.start()

    # Add jobs to queue
    for input_file in input_files:
        # Generate output file name by replacing .gjf with .log
        output_file = os.path.splitext(input_file)[0] + ".log"
        job_queue.put((input_file, output_file))

    # Setup signal handler for graceful termination
    original_sigint_handler = signal.getsignal(signal.SIGINT)

    def signal_handler(sig, frame):
        _console.print("\n[yellow]⚠[/yellow] Received interrupt signal. Stopping batch processing...")
        # Put sentinel values to stop workers
        for _ in range(max_concurrent):
            job_queue.put(None)
        # Restore original handler so second Ctrl+C will terminate immediately
        signal.signal(signal.SIGINT, original_sigint_handler)

    signal.signal(signal.SIGINT, signal_handler)

    try:
        # Process results as they come in
        start_time = time.time()
        completed_count = 0
        failed_count = 0

        # Wait until we've processed all files
        while completed_count + failed_count < len(input_files):
            try:
                # Get result from queue with timeout
                result = result_queue.get(timeout=2)

                if result["success"]:
                    completed_count += 1
                    _console.print(
                        f"[green]✓[/green] [{completed_count + failed_count}/{len(input_files)}] "
                        f"Completed: [cyan]{os.path.basename(result['input_file'])}[/cyan] → "
                        f"[cyan]{os.path.basename(result['output_file'])}[/cyan] "
                        f"({result.get('elapsed_time', 0):.1f}s)"
                    )
                else:
                    failed_count += 1
                    _console.print(
                        f"[red]✗[/red] [{completed_count + failed_count}/{len(input_files)}] "
                        f"Failed: [cyan]{os.path.basename(result['input_file'])}[/cyan] → "
                        f"{result['error']}"
                    )

                # Store result
                results["results"].append(result)
                result_queue.task_done()

            except queue.Empty:
                # If no results for 2 seconds, just continue waiting
                continue

        end_time = time.time()
        total_time = end_time - start_time

        # Add sentinel values to stop workers
        for _ in range(max_concurrent):
            job_queue.put(None)

        # Wait for workers to terminate
        for worker in workers:
            worker.join(timeout=5)

        # Update results
        results["completed"] = completed_count
        results["failed"] = failed_count
        results["total_time"] = total_time

        _console.print("\n[bold]Batch processing completed:[/bold]")
        _console.print(f"  Total files processed: [cyan]{completed_count + failed_count}[/cyan]")
        _console.print(f"  Successfully completed: [green]{completed_count}[/green]")
        _console.print(f"  Failed: [red]{failed_count}[/red]")
        _console.print(f"  Total elapsed time: [cyan]{total_time:.1f}[/cyan] seconds")

        return results

    except KeyboardInterrupt:
        _console.print("\n[yellow]⚠[/yellow] Batch processing interrupted by user.")
        # Results will contain what was processed so far
        results["completed"] = completed_count
        results["failed"] = failed_count
        return results
    except Exception as e:
        _console.print(f"[red]✗[/red] Error during batch processing: {str(e)}")
        results["error"] = str(e)
        return results
    finally:
        # Restore original signal handler
        signal.signal(signal.SIGINT, original_sigint_handler)


def find_gaussian_executable():
    """
    Find the Gaussian executable by searching common names and locations.

    Returns:
    --------
    str
        Path to Gaussian executable, or None if not found
    """
    # Common Gaussian executable names
    common_names = ["g16", "g09", "gaussian", "g16.exe", "g09.exe"]

    # Check if any of the common names are in PATH
    for name in common_names:
        try:
            subprocess.run(
                [name, "--version"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            return name
        except FileNotFoundError:
            pass

    # If not found in PATH, check common installation directories
    common_paths = [
        "/opt/gaussian/g16",
        "/usr/local/gaussian/g16",
        "/opt/gaussian/g09",
        "/usr/local/gaussian/g09",
        "C:\\G16W\\g16.exe",
        "C:\\G09W\\g09.exe",
    ]

    for path in common_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            return path

    return None


def run_gaussian_batch_mode(amino_acid, processor=None):
    """
    Run Gaussian calculations in batch mode for a given amino acid.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code

    Returns:
    --------
    dict
        A dictionary containing success status and results
    """
    # Check for input files
    amino_acid = amino_acid.strip().lower()

    # Find files in ahelix and bsheet directories
    ahelix_dir = f"{amino_acid}_ahelix_structures"
    bsheet_dir = f"{amino_acid}_bsheet_structures"

    ahelix_files = sorted(glob.glob(os.path.join(ahelix_dir, "*.gjf")))
    bsheet_files = sorted(glob.glob(os.path.join(bsheet_dir, "*.gjf")))

    # Combine files with ahelix first, then bsheet
    all_files = ahelix_files + bsheet_files

    if not all_files:
        _console.print(
            f"[red]✗[/red] No Gaussian input files found for [cyan]{amino_acid}[/cyan] in expected directories:"
        )
        _console.print(f"  - {ahelix_dir}/*.gjf")
        _console.print(f"  - {bsheet_dir}/*.gjf")
        return {"success": False, "error": "No input files found"}

    # Find Gaussian executable
    gaussian_exe = find_gaussian_executable()
    if not gaussian_exe:
        _console.print("[yellow]⚠[/yellow] Could not find Gaussian executable. Please specify the path manually.")
        gaussian_exe = prompt_with_context(
            processor,
            "Enter path to Gaussian executable (g16, g09, etc.)",
            module="Modified AA Parameterizer",
            description="Path to Gaussian executable",
        )

        # Verify executable exists
        if not os.path.exists(gaussian_exe):
            return {
                "success": False,
                "error": f"Gaussian executable not found: {gaussian_exe}",
            }

    # Ask for maximum concurrent calculations
    try:
        max_concurrent = int_prompt_with_context(
            processor,
            "Enter maximum number of concurrent Gaussian calculations",
            default=1,
            module="Modified AA Parameterizer",
            description="Maximum concurrent Gaussian calculations",
        )
        if max_concurrent < 1:
            max_concurrent = 1
    except ValueError:
        _console.print("[yellow]⚠[/yellow] Invalid input. Using default value of 1.")
        max_concurrent = 1

    # Run in batch mode
    results = batch_run_gaussian_files(all_files, max_concurrent, gaussian_exe)

    return {"success": True, "results": results}


def run_batch_for_extracted_structures(amino_acid, structure_dirs, processor=None):
    """
    Run Gaussian calculations in batch mode for extracted structures.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    structure_dirs : list
        List of directories containing extracted structure files

    Returns:
    --------
    dict
        A dictionary containing success status and results
    """
    # Collect all .gjf files from the structure directories
    all_files = []
    for dir_path in structure_dirs:
        if os.path.exists(dir_path):
            gjf_files = sorted(glob.glob(os.path.join(dir_path, "*.gjf")))
            all_files.extend(gjf_files)

    if not all_files:
        _console.print(f"[red]✗[/red] No Gaussian input files found in the specified directories.")
        return {"success": False, "error": "No input files found"}

    _console.print(f"\n[green]✓[/green] Found [cyan]{len(all_files)}[/cyan] Gaussian input files to process.")
    for dir_path in structure_dirs:
        if os.path.exists(dir_path):
            file_count = len(glob.glob(os.path.join(dir_path, "*.gjf")))
            if file_count > 0:
                _console.print(f"  - {dir_path}: [cyan]{file_count}[/cyan] files")

    # Ask if the user wants to run these calculations
    run_gaussian = confirm_with_context(
        processor,
        "\nDo you want to run Gaussian calculations for these extracted structures?",
        default=False,
        module="Modified AA Parameterizer",
        description="Run Gaussian for extracted structures",
    )

    if not run_gaussian:
        _console.print("[yellow]○[/yellow] Skipping Gaussian calculations. You can run them manually later.")
        return {"success": False, "skipped": True}

    # Ask for running mode
    _console.print(Panel(
        "[bold]1[/bold]. Batch mode (run multiple calculations in parallel)\n"
        "[bold]2[/bold]. Manual mode (instructions for running manually)",
        title="How do you want to run Gaussian calculations?",
        expand=False
    ))
    mode = prompt_with_context(
        processor,
        "Enter option",
        choices=["1", "2"],
        default="1",
        module="Modified AA Parameterizer",
        description="Run mode for extracted structures",
        options_map={"1": "Run all in batch", "2": "Run manually"},
    )

    if mode == "1":
        # Find Gaussian executable
        gaussian_exe = find_gaussian_executable()
        if not gaussian_exe:
            _console.print(
                "[yellow]⚠[/yellow] Could not find Gaussian executable. Please specify the path manually."
            )
            gaussian_exe = prompt_with_context(
                processor,
                "Enter path to Gaussian executable (g16, g09, etc.)",
                module="Modified AA Parameterizer",
                description="Path to Gaussian executable",
            )

            # Verify executable exists
            if not os.path.exists(gaussian_exe):
                return {
                    "success": False,
                    "error": f"Gaussian executable not found: {gaussian_exe}",
                }

        # Ask for maximum concurrent calculations
        try:
            max_concurrent = int_prompt_with_context(
                processor,
                "Enter maximum number of concurrent Gaussian calculations",
                default=1,
                module="Modified AA Parameterizer",
                description="Maximum concurrent Gaussian calculations",
            )
            if max_concurrent < 1:
                max_concurrent = 1
        except ValueError:
            _console.print("[yellow]⚠[/yellow] Invalid input. Using default value of 1.")
            max_concurrent = 1

        # Run in batch mode
        results = batch_run_gaussian_files(all_files, max_concurrent, gaussian_exe)

        return {"success": True, "results": results}
    else:
        # Show manual mode instructions
        dir_list = "\n".join(f"     - {dir_path}/" for dir_path in structure_dirs if os.path.exists(dir_path))
        manual_instructions = f"""[bold]To run the Gaussian calculations manually:[/bold]

  1. Navigate to each of these directories:
{dir_list}

  2. For each .gjf file, run Gaussian with the command:
     [cyan]g16 input.gjf > input.log[/cyan]
     (Replace 'g16' with the appropriate Gaussian command for your system)

  3. After all calculations are complete, run this script again
     to analyze the results"""
        _console.print(Panel(
            manual_instructions,
            title="Manual Gaussian Calculation Instructions",
            expand=False
        ))

        return {"success": True, "manual_mode": True}


def extract_pes_structures(log_file, output_dir, gaussian_params, atoms_to_remove=None):
    """
    Extract structures from a Gaussian PES scan log file and save converged ones as Gaussian input files.
    Optionally removes specified atoms from the extracted structures and sets frozen atoms for optimization.

    Parameters:
    -----------
    log_file : str
        Path to the Gaussian log file
    output_dir : str
        Directory to save extracted structures
    gaussian_params : dict
        Parameters for Gaussian input files, including:
        - memory: Memory specification
        - procs: Number of processors
        - keywords: Gaussian keywords
        - charge: Molecular charge
        - multiplicity: Spin multiplicity
        - frozen_atoms: Can be a list of atom indices to freeze, "ALL" to freeze everything,
                       or a dict with {'mode': 'ALL_EXCEPT', 'except_list': [indices]}
    atoms_to_remove : list, optional
        List of atom indices to remove from the structures

    Returns:
    --------
    dict
        Information about extracted structures, convergence, etc.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Initialize result data
    result = {
        "total_structures": 0,
        "converged": 0,
        "failed": 0,
        "files_created": [],
        "failed_points": [],
    }

    _console.print(f"\n[cyan]Processing log file:[/cyan] {log_file}")
    _console.print(f"[cyan]Saving structures to:[/cyan] {output_dir}")

    # Print info about atoms to remove
    if atoms_to_remove:
        _console.print(
            f"[yellow]Removing {len(atoms_to_remove)} atoms with indices: {', '.join(map(str, atoms_to_remove))}[/yellow]"
        )

    # Print info about frozen atoms
    frozen_atoms = gaussian_params.get("frozen_atoms", [])
    if frozen_atoms == "ALL":
        _console.print("[cyan]ℹ All atoms will be frozen during optimization[/cyan]")
    elif isinstance(frozen_atoms, dict) and frozen_atoms.get("mode") == "ALL_EXCEPT":
        except_list = frozen_atoms.get("except_list", [])
        _console.print(
            f"[cyan]ℹ All atoms will be frozen EXCEPT indices: {', '.join(map(str, except_list))}[/cyan]"
        )
    elif frozen_atoms:
        _console.print(
            f"[cyan]ℹ Freezing {len(frozen_atoms)} atoms with indices: {', '.join(map(str, frozen_atoms))}[/cyan]"
        )

    try:
        # Read the entire log file
        with open(log_file, "r") as f:
            lines = f.readlines()

        # Extract base filename without extension
        base_name = os.path.splitext(os.path.basename(log_file))[0]

        # Find all scan point information
        scan_points = {}  # Will store scan_point -> [list of step numbers]
        current_scan_point = None
        current_step = None

        # Find optimized structures
        opt_param_lines = []
        for i, line in enumerate(lines):
            if "!   Optimized Parameters   !" in line:
                opt_param_lines.append(i)

        # Process each line to track scan points and steps
        for i, line in enumerate(lines):
            # Detect scan point lines with correct pattern matching
            if "on scan point" in line:
                parts = line.split()
                try:
                    # Format: "Step number X out of a maximum of Y on scan point Z out of W"
                    step_idx = parts.index("number") + 1
                    scan_idx = parts.index("point") + 1

                    step_num = int(parts[step_idx])
                    scan_num = int(parts[scan_idx])

                    current_scan_point = scan_num
                    current_step = step_num

                    if scan_num not in scan_points:
                        scan_points[scan_num] = []
                    scan_points[scan_num].append(step_num)
                except (ValueError, IndexError) as e:
                    _console.print(f"[red]✗ Error parsing line: {line.strip()}[/red]")
                    _console.print(f"[red]  Exception: {str(e)}[/red]")
                    continue

        _console.print(f"[green]✓ Found {len(scan_points)} scan points:[/green]")
        for scan_num, steps in sorted(scan_points.items()):
            _console.print(
                f"  Scan point {scan_num}: {len(steps)} steps (steps: {min(steps)}-{max(steps)})"
            )

        _console.print(f"[green]✓ Found {len(opt_param_lines)} optimized parameter sections[/green]")

        # Process each optimized parameter section
        extracted_structures = {}  # To prevent duplicates

        for i in opt_param_lines:
            # Look backward to find the most recent scan point info
            scan_num = None
            step_num = None

            for j in range(i, max(0, i - 500), -1):
                if "on scan point" in lines[j]:
                    parts = lines[j].split()
                    try:
                        step_idx = parts.index("number") + 1
                        scan_idx = parts.index("point") + 1

                        step_num = int(parts[step_idx])
                        scan_num = int(parts[scan_idx])
                        break
                    except (ValueError, IndexError):
                        continue

            if scan_num is None:
                _console.print(
                    f"[yellow]⚠ Warning: Could not determine scan point for optimized parameters at line {i}[/yellow]"
                )
                continue

            # Check if we've already extracted a structure for this scan point
            # If so, we only want the last one (final convergence)
            if scan_num in extracted_structures:
                old_step = extracted_structures[scan_num][0]
                if step_num <= old_step:
                    # We already have a more recent structure for this scan point
                    continue

            # Look for energy
            energy = "unknown"
            for j in range(i, max(0, i - 100), -1):
                if "SCF Done:" in lines[j]:
                    parts = lines[j].split()
                    try:
                        energy_idx = parts.index("=") + 1
                        energy = parts[energy_idx]
                    except (ValueError, IndexError):
                        pass
                    break

            # Look for coordinates
            coords_start = None
            for j in range(i, max(0, i - 500), -1):
                if "Standard orientation:" in lines[j]:
                    coords_start = j + 5  # Coordinates start 5 lines after this
                    break

            if coords_start is None:
                _console.print(
                    f"[yellow]⚠ Warning: Could not find standard orientation for scan point {scan_num}[/yellow]"
                )
                continue

            # Find the end of the coordinate section
            coords_end = None
            for j in range(coords_start, min(len(lines), coords_start + 200)):
                if (
                    "---------------------------------------------------------------------"
                    in lines[j]
                ):
                    coords_end = j
                    break

            if coords_end is None:
                _console.print(
                    f"[yellow]⚠ Warning: Could not find end of coordinates for scan point {scan_num}[/yellow]"
                )
                continue

            # Extract coordinates
            coords_block = lines[coords_start:coords_end]

            # Store that we've processed this scan point
            extracted_structures[scan_num] = (step_num, coords_block, energy)

        # Create files for each extracted structure
        for scan_num, (step_num, coords_block, energy) in sorted(
            extracted_structures.items()
        ):
            # Create Gaussian input file named by scan point
            output_file = os.path.join(output_dir, f"{base_name}_scan{scan_num}.gjf")

            with open(output_file, "w") as out_f:
                # Write Gaussian input file header
                out_f.write(f"%mem={gaussian_params['memory']}\n")
                out_f.write(f"%nprocshared={gaussian_params['procs']}\n")

                # Use keywords as provided - the -1/0 flags work with standard opt keyword
                # No need to add modredundant unless specifically requested
                keywords = gaussian_params["keywords"]

                # Ensure there's an opt keyword if we're freezing atoms
                if (
                    "frozen_atoms" in gaussian_params
                    and gaussian_params["frozen_atoms"]
                ):
                    if "opt" not in keywords:
                        # If opt is not specified at all, add it
                        keywords = f"{keywords} opt"

                out_f.write(f"#p {keywords}\n\n")
                out_f.write(
                    f"Structure from PES scan - Scan point {scan_num}, Step {step_num}, Energy {energy}\n\n"
                )
                out_f.write(
                    f"{gaussian_params['charge']} {gaussian_params['multiplicity']}\n"
                )

                # Process and write coordinates, skipping atoms to be removed
                processed_atoms = (
                    []
                )  # Keep track of processed atoms to update atom indices for constraints
                atom_index_map = (
                    {}
                )  # Maps original atom indices to new indices after removal

                for coord_line in coords_block:
                    parts = coord_line.split()
                    if len(parts) >= 6:
                        try:
                            atom_idx = int(parts[0])

                            # Skip this atom if it's in the removal list
                            if atoms_to_remove and atom_idx in atoms_to_remove:
                                continue

                            atom_type = int(parts[1])

                            # Convert atomic number to element symbol
                            element = get_element_from_atomic_number(atom_type)

                            x, y, z = parts[3:6]

                            # Check if this atom should be frozen
                            frozen_atoms = gaussian_params.get("frozen_atoms", [])
                            freeze_flag = 0  # Default to not frozen (0)

                            if frozen_atoms == "ALL":
                                # All atoms should be frozen
                                freeze_flag = -1
                            elif (
                                isinstance(frozen_atoms, dict)
                                and frozen_atoms.get("mode") == "ALL_EXCEPT"
                            ):
                                # All atoms except those in the except_list should be frozen
                                except_list = frozen_atoms.get("except_list", [])
                                if atom_idx not in except_list:
                                    freeze_flag = -1
                            elif atom_idx in frozen_atoms:
                                # Specific atoms are frozen
                                freeze_flag = -1

                            out_f.write(f"{element} {freeze_flag} {x} {y} {z}\n")

                            # Add to processed atoms and update index map
                            processed_atoms.append(atom_idx)
                            atom_index_map[atom_idx] = len(processed_atoms)

                        except (ValueError, IndexError) as e:
                            # Skip lines that can't be properly parsed
                            continue

                out_f.write("\n")  # End of file

                # If there are any constraints in the original calculation, we need to update them
                # This would be expanded for specific constraint types if needed
                if "constraints" in gaussian_params and gaussian_params["constraints"]:
                    for constraint in gaussian_params["constraints"]:
                        # Example of updating a dihedral constraint
                        # Format: "D atom1 atom2 atom3 atom4 F value"
                        if constraint.startswith("D"):
                            parts = constraint.split()
                            if len(parts) >= 5:
                                atoms = [
                                    int(parts[1]),
                                    int(parts[2]),
                                    int(parts[3]),
                                    int(parts[4]),
                                ]
                                # Check if all atoms in constraint still exist after removal
                                if all(atom in atom_index_map for atom in atoms):
                                    # Update the atom indices
                                    new_atoms = [atom_index_map[atom] for atom in atoms]
                                    # Write the updated constraint
                                    out_f.write(
                                        f"D {new_atoms[0]} {new_atoms[1]} {new_atoms[2]} {new_atoms[3]} {' '.join(parts[5:])}\n"
                                    )

            result["files_created"].append(output_file)
            result["converged"] += 1
            _console.print(
                f"[green]✓ Scan point {scan_num} (Step {step_num}):[/green] Created {os.path.basename(output_file)}"
            )

        # Total number of scan points
        result["total_structures"] = len(scan_points)

        _console.print(f"\n[bold]Summary for {log_file}:[/bold]")
        _console.print(f"  Total scan points found: {result['total_structures']}")
        _console.print(f"  Converged structures extracted: {result['converged']}")
        _console.print(f"  Files created: {len(result['files_created'])}")

        return result

    except Exception as e:
        import traceback

        traceback.print_exc()
        _console.print(f"[red]✗ Error processing log file {log_file}: {str(e)}[/red]")
        return {"success": False, "error": str(e)}


def get_element_from_atomic_number(atomic_number):
    """
    Convert atomic number to element symbol.

    Parameters:
    -----------
    atomic_number : int
        Atomic number of the element

    Returns:
    --------
    str
        Element symbol
    """
    elements = {
        1: "H",
        2: "He",
        3: "Li",
        4: "Be",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "Ne",
        11: "Na",
        12: "Mg",
        13: "Al",
        14: "Si",
        15: "P",
        16: "S",
        17: "Cl",
        18: "Ar",
        19: "K",
        20: "Ca",
        21: "Sc",
        22: "Ti",
        23: "V",
        24: "Cr",
        25: "Mn",
        26: "Fe",
        27: "Co",
        28: "Ni",
        29: "Cu",
        30: "Zn",
        31: "Ga",
        32: "Ge",
        33: "As",
        34: "Se",
        35: "Br",
        36: "Kr",
        37: "Rb",
        38: "Sr",
        39: "Y",
        40: "Zr",
        41: "Nb",
        42: "Mo",
        43: "Tc",
        44: "Ru",
        45: "Rh",
        46: "Pd",
        47: "Ag",
        48: "Cd",
        49: "In",
        50: "Sn",
        51: "Sb",
        52: "Te",
        53: "I",
        54: "Xe",
    }
    return elements.get(atomic_number, "X")  # Return 'X' if not found


def analyze_pes_log_files(amino_acid, log_files, processor=None):
    """
    Analyze PES scan log files and extract converged structures.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    log_files : list
        List of Gaussian log files to analyze
    processor : optional
        Processor object with session_manager for session recording support.

    Returns:
    --------
    dict
        Summary of analysis results
    """
    _console.print(f"\n[bold cyan]Step 5: Analyzing PES scan results for {amino_acid}[/bold cyan]")

    # Get Gaussian parameters for output files from user (asked only once)
    _console.print("\n[bold]Specify Gaussian settings for all extracted structures:[/bold]")
    memory = prompt_with_context(
        processor, "Enter memory for Gaussian job", default="8GB",
        module="Modified Amino Acid Parameterizer", description="Gaussian memory allocation"
    )
    procs = prompt_with_context(
        processor, "Enter number of processors for Gaussian job", default="4",
        module="Modified Amino Acid Parameterizer", description="Number of processors"
    )
    keywords = prompt_with_context(
        processor, "Enter Gaussian keywords",
        default="opt HF/6-31g(d) pop=mk iop(6/33=2,6/42=6)",
        module="Modified Amino Acid Parameterizer", description="Gaussian keywords"
    )

    # Get charge and multiplicity
    charge = prompt_with_context(
        processor, "Enter charge for the molecule", default="0",
        module="Modified Amino Acid Parameterizer", description="Molecular charge"
    )
    multiplicity = prompt_with_context(
        processor, "Enter multiplicity (singlet=1, doublet=2, etc.)", default="1",
        module="Modified Amino Acid Parameterizer", description="Spin multiplicity"
    )

    try:
        charge = int(charge)
        multiplicity = int(multiplicity)
    except ValueError:
        _console.print("[red]✗ Error: Charge and multiplicity must be integers. Using defaults (0 1).[/red]")
        charge = 0
        multiplicity = 1

    # Ask if any atoms should be removed
    remove_atoms = confirm_with_context(
        processor, "\nDo you want to remove any atoms from the extracted structures?",
        default=False, module="Modified Amino Acid Parameterizer",
        description="Remove atoms from structures"
    )
    atoms_to_remove = []

    if remove_atoms:
        atoms_input = prompt_with_context(
            processor, "Enter the indices of atoms to remove (space-separated, e.g., '1 5 9')",
            default="", module="Modified Amino Acid Parameterizer",
            description="Atom indices to remove"
        )
        try:
            # Parse the input into a list of integers
            if atoms_input.strip():
                atoms_to_remove = [int(idx) for idx in atoms_input.split()]
                _console.print(
                    f"[cyan]ℹ Will remove {len(atoms_to_remove)} atoms with indices: {', '.join(map(str, atoms_to_remove))}[/cyan]"
                )
        except ValueError:
            _console.print("[red]✗ Error parsing atom indices. No atoms will be removed.[/red]")
            atoms_to_remove = []

    # Ask about freezing atoms in the optimization
    frozen_atoms = []
    _console.print(Panel(
        "[bold]Atom Freezing Options:[/bold]\n\n"
        "[cyan]1.[/cyan] No atoms frozen (all atoms optimize)\n"
        "[cyan]2.[/cyan] Freeze all atoms\n"
        "[cyan]3.[/cyan] Freeze specific atoms\n"
        "[cyan]4.[/cyan] Freeze all EXCEPT specific atoms",
        title="Optimization Constraints",
        border_style="blue",
        expand=False
    ))
    freeze_option = prompt_with_context(
        processor, "Enter option", default="1", choices=["1", "2", "3", "4"],
        module="Modified Amino Acid Parameterizer", description="Atom freezing option",
        options_map={"1": "No atoms frozen", "2": "Freeze all", "3": "Freeze specific", "4": "Freeze all except"}
    )

    try:
        freeze_option = int(freeze_option)

        if freeze_option == 1:
            # No atoms frozen
            _console.print("[cyan]ℹ No atoms will be frozen during optimization.[/cyan]")

        elif freeze_option == 2:
            # Freeze all atoms
            _console.print("[cyan]ℹ All atoms will be frozen during optimization.[/cyan]")
            # We'll handle this by setting a flag and assigning the atoms later
            # after we know how many atoms there are in each structure
            frozen_atoms = "ALL"

        elif freeze_option == 3:
            # Freeze specific atoms
            atoms_input = prompt_with_context(
                processor, "Enter the indices of atoms to freeze (space-separated, e.g., '1 5 9')",
                default="", module="Modified Amino Acid Parameterizer",
                description="Atom indices to freeze"
            )
            if atoms_input.strip():
                frozen_atoms = [int(idx) for idx in atoms_input.split()]
                _console.print(
                    f"[cyan]ℹ Will freeze {len(frozen_atoms)} atoms with indices: {', '.join(map(str, frozen_atoms))}[/cyan]"
                )
            else:
                _console.print("[cyan]ℹ No atoms specified, no atoms will be frozen.[/cyan]")

        elif freeze_option == 4:
            # Freeze all EXCEPT specific atoms
            atoms_input = prompt_with_context(
                processor, "Enter the indices of atoms to NOT freeze (space-separated, e.g., '1 5 9')",
                default="", module="Modified Amino Acid Parameterizer",
                description="Atom indices to NOT freeze"
            )
            if atoms_input:
                unfreeze_atoms = [int(idx) for idx in atoms_input.split()]
                _console.print(
                    f"[cyan]ℹ Will allow {len(unfreeze_atoms)} atoms to optimize with indices: {', '.join(map(str, unfreeze_atoms))}[/cyan]"
                )
                frozen_atoms = {"mode": "ALL_EXCEPT", "except_list": unfreeze_atoms}
            else:
                _console.print("[cyan]ℹ No atoms specified, all atoms will be frozen.[/cyan]")
                frozen_atoms = "ALL"

        else:
            _console.print("[yellow]⚠ Invalid option. No atoms will be frozen during optimization.[/yellow]")

    except ValueError:
        _console.print("[red]✗ Error parsing option. No atoms will be frozen during optimization.[/red]")

    # Store parameters for passing to extraction function
    gaussian_params = {
        "memory": memory,
        "procs": procs,
        "keywords": keywords,
        "charge": charge,
        "multiplicity": multiplicity,
        "frozen_atoms": frozen_atoms,
    }

    results = {
        "success": True,
        "ahelix": {"structures": 0, "converged": 0, "failed": 0, "output_dir": ""},
        "bsheet": {"structures": 0, "converged": 0, "failed": 0, "output_dir": ""},
    }

    structure_dirs = []  # Track directories where structures are saved

    for log_file in log_files:
        # Determine conformation type from filename
        if "ahelix" in log_file:
            conformation = "ahelix"
        elif "bsheet" in log_file:
            conformation = "bsheet"
        else:
            _console.print(
                f"[yellow]⚠ Warning: Could not determine conformation type from filename: {log_file}[/yellow]"
            )
            continue

        # Create output directory for this conformation
        output_dir = f"{amino_acid.lower()}_{conformation}_structures"
        structure_dirs.append(output_dir)
        results[conformation]["output_dir"] = output_dir

        # Extract structures from log file, passing the atoms to remove
        extraction_result = extract_pes_structures(
            log_file, output_dir, gaussian_params, atoms_to_remove
        )

        if "success" in extraction_result and extraction_result["success"] == False:
            _console.print(f"[red]✗ Error analyzing {log_file}: {extraction_result['error']}[/red]")
            results["success"] = False
            continue

        # Update results
        results[conformation]["structures"] += extraction_result["total_structures"]
        results[conformation]["converged"] += extraction_result["converged"]
        results[conformation]["failed"] += extraction_result["failed"]

    # Print overall summary using a Table
    _console.print("\n")
    summary_table = Table(title=f"Analysis Summary for {amino_acid.upper()}")
    summary_table.add_column("Conformation", style="cyan")
    summary_table.add_column("Total Scan Points", style="white")
    summary_table.add_column("Converged", style="green")
    summary_table.add_column("Failed", style="red")
    summary_table.add_column("Output Directory", style="grey50")

    for conf in ["ahelix", "bsheet"]:
        if results[conf]["structures"] > 0:
            summary_table.add_row(
                conf.upper(),
                str(results[conf]["structures"]),
                str(results[conf]["converged"]),
                str(results[conf]["failed"]),
                f"{results[conf]['output_dir']}/"
            )
    _console.print(summary_table)

    # NEW: Ask if the user wants to run Gaussian calculations on the extracted structures
    if results["success"] and (
        results["ahelix"]["converged"] > 0 or results["bsheet"]["converged"] > 0
    ):
        batch_result = run_batch_for_extracted_structures(amino_acid, structure_dirs, processor=processor)

    return results


def process_opt_only_workflow(amino_acid, opt_log_files, processor=None):
    """
    Process optimization-only workflow by extracting final geometries from
    optimization logs and generating ESP input files.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    opt_log_files : list
        List of optimization log files (e.g., ['lys_ahelix.log', 'lys_bsheet.log'])
    processor : optional
        Processor object with session_manager for session recording support.

    Returns:
    --------
    dict
        Dictionary containing success status and any error messages
    """
    _console.print(f"\n[bold cyan]Step 5: Generating ESP inputs from optimized structures for {amino_acid.upper()}[/bold cyan]")
    _console.print(f"[cyan]ℹ Using {len(opt_log_files)} optimized structures[/cyan]")

    amino_acid = amino_acid.strip().lower()

    # Try to find mol2 files to get dihedral atoms for backbone constraints
    dihedral_info = {}
    mol2_files = glob.glob(f"{amino_acid}_*.mol2")
    for mol2_file in mol2_files:
        if "ahelix" in mol2_file:
            parse_result = parse_mol2_for_dihedral_atoms(mol2_file)
            if parse_result.get("success"):
                dihedral_info["ahelix"] = {"phi": parse_result["phi"], "psi": parse_result["psi"]}
        elif "bsheet" in mol2_file:
            parse_result = parse_mol2_for_dihedral_atoms(mol2_file)
            if parse_result.get("success"):
                dihedral_info["bsheet"] = {"phi": parse_result["phi"], "psi": parse_result["psi"]}

    # Ask about protonation state changes
    _console.print(Panel(
        "[bold]Protonation State[/bold]\n\n"
        "Will you modify the protonation state for the ESP calculation?\n"
        "[grey50](e.g., deprotonate a hydroxyl to create an oxide)[/grey50]",
        border_style="blue",
        expand=False
    ))
    protonation_change = confirm_with_context(
        processor, "Protonation state will change?", default=False,
        module="Modified Amino Acid Parameterizer", description="Protonation state change"
    )

    # Set default keywords based on protonation state
    if protonation_change:
        _console.print(Panel(
            "[bold]Geometry Constraints[/bold]\n\n"
            "Re-optimization is recommended to relax the modified functional group.\n"
            "However, without constraints, the α-helix and β-sheet conformations\n"
            "may converge to similar structures, losing conformational diversity.\n\n"
            "[cyan]d[/cyan] - Freeze backbone dihedrals (φ/ψ) only\n"
            "[cyan]f[/cyan] - Freeze all atoms, then you edit to unfreeze specific atoms",
            border_style="blue",
            expand=False
        ))

        constraint_choice = prompt_with_context(
            processor, "Constraint method?", default="d", choices=["d", "f"],
            module="Modified Amino Acid Parameterizer", description="Geometry constraint method",
            options_map={"d": "Freeze dihedrals", "f": "Freeze atoms"}
        )
        use_freeze_atoms = (constraint_choice == "f")

        if use_freeze_atoms:
            default_keywords = "opt HF/6-31G(d) Pop=mk IOp(6/33=2,6/42=6)"
            _console.print("\n[cyan]ℹ Using freeze atom approach:[/cyan]")
            _console.print("  All atoms will be marked frozen (-1). Edit the .gjf file to change")
            _console.print("  atoms you want to optimize from -1 to 0 (typically the modified group).")
        else:
            default_keywords = "opt=modredundant HF/6-31G(d) Pop=mk IOp(6/33=2,6/42=6)"
            if dihedral_info:
                _console.print("\n[cyan]ℹ Using dihedral freeze approach:[/cyan]")
                _console.print("  Backbone φ/ψ dihedrals will be frozen to preserve conformation.")
            else:
                _console.print("\n[yellow]⚠ Warning: Could not find mol2 files to extract dihedral atoms.[/yellow]")
                _console.print("You may need to manually add dihedral constraints to the .gjf files.")
    else:
        # No protonation change - single-point ESP (no optimization needed)
        default_keywords = "HF/6-31G(d) Pop=mk IOp(6/33=2,6/42=6)"
        use_freeze_atoms = False
        _console.print("\n[cyan]ℹ Using single-point ESP (no re-optimization).[/cyan]")

    # Get Gaussian parameters for ESP calculations
    _console.print(Panel("[bold]Gaussian Settings[/bold]", border_style="blue", expand=False))
    memory = prompt_with_context(
        processor, "Enter memory for Gaussian job", default="8GB",
        module="Modified Amino Acid Parameterizer", description="Gaussian memory allocation"
    )
    procs = prompt_with_context(
        processor, "Enter number of processors for Gaussian job", default="4",
        module="Modified Amino Acid Parameterizer", description="Number of processors"
    )
    keywords = prompt_with_context(
        processor, f"Enter Gaussian keywords", default=default_keywords,
        module="Modified Amino Acid Parameterizer", description="Gaussian keywords"
    )

    # Get charge and multiplicity
    charge = prompt_with_context(
        processor, "Enter charge for the molecule", default="0",
        module="Modified Amino Acid Parameterizer", description="Molecular charge"
    )
    multiplicity = prompt_with_context(
        processor, "Enter multiplicity (singlet=1, doublet=2, etc.)", default="1",
        module="Modified Amino Acid Parameterizer", description="Spin multiplicity"
    )

    try:
        charge = int(charge)
        multiplicity = int(multiplicity)
    except ValueError:
        _console.print("[red]✗ Error: Charge and multiplicity must be integers. Using defaults (0 1).[/red]")
        charge = 0
        multiplicity = 1

    # Create output directories for each conformation
    created_dirs = []
    extracted_count = 0

    for opt_log in opt_log_files:
        # Determine conformation type from filename
        if "ahelix" in opt_log.lower():
            conf_type = "ahelix"
        elif "bsheet" in opt_log.lower():
            conf_type = "bsheet"
        else:
            _console.print(f"[yellow]⚠ Warning: Could not determine conformation type from {opt_log}, skipping[/yellow]")
            continue

        # Create output directory
        output_dir = f"{amino_acid}_{conf_type}_structures"
        os.makedirs(output_dir, exist_ok=True)
        created_dirs.append(output_dir)

        _console.print(f"\n[cyan]→ Extracting optimized geometry from {opt_log}...[/cyan]")

        # Extract final optimized geometry
        geometry = extract_final_geometry_from_log(opt_log)

        if not geometry:
            _console.print(f"[red]✗ Error: Could not extract geometry from {opt_log}[/red]")
            continue

        # Generate Gaussian input file for ESP calculation
        output_gjf = os.path.join(output_dir, f"{amino_acid}_{conf_type}_opt.gjf")

        with open(output_gjf, 'w') as f:
            f.write(f"%mem={memory}\n")
            f.write(f"%nprocshared={procs}\n")
            f.write(f"# {keywords}\n\n")
            f.write(f"{amino_acid.upper()} {conf_type} ESP calculation\n\n")
            f.write(f"{charge} {multiplicity}\n")

            # Write coordinates (with freeze flags if using freeze atoms approach)
            if use_freeze_atoms and protonation_change:
                for atom_line in geometry:
                    # Insert -1 (frozen) flag after element symbol
                    parts = atom_line.split()
                    if len(parts) >= 4:
                        element = parts[0]
                        coords = " ".join(parts[1:])
                        f.write(f"{element:2s}  -1  {coords}\n")
                    else:
                        f.write(f"{atom_line}\n")
            else:
                for atom_line in geometry:
                    f.write(f"{atom_line}\n")

            f.write("\n")  # Blank line after coordinates

            # Add dihedral constraints if using dihedral freeze approach
            if protonation_change and not use_freeze_atoms and conf_type in dihedral_info:
                phi = dihedral_info[conf_type]["phi"]
                psi = dihedral_info[conf_type]["psi"]
                f.write(f"D {phi[0]} {phi[1]} {phi[2]} {phi[3]} F\n")
                f.write(f"D {psi[0]} {psi[1]} {psi[2]} {psi[3]} F\n")
                f.write("\n")

        _console.print(f"[green]✓ Created ESP input: {output_gjf}[/green]")
        if use_freeze_atoms and protonation_change:
            _console.print(f"  [grey50]Note: All atoms frozen (-1). Edit file to set atoms to optimize to 0.[/grey50]")
        elif protonation_change and not use_freeze_atoms and conf_type in dihedral_info:
            _console.print(f"  [grey50]Note: φ/ψ dihedral constraints added to preserve backbone conformation.[/grey50]")
        extracted_count += 1

    if extracted_count == 0:
        return {
            "success": False,
            "error": "Failed to extract geometries from optimization logs"
        }

    _console.print(f"\n[green]{'='*70}[/green]")
    _console.print(f"[green]✓ Successfully created {extracted_count} ESP input files[/green]")
    _console.print(f"[green]{'='*70}[/green]")
    _console.print("\n[bold]Next steps:[/bold]")
    _console.print("  1. Run these Gaussian ESP calculations:")
    for dir_name in created_dirs:
        gjf_files = glob.glob(os.path.join(dir_name, "*.gjf"))
        for gjf in gjf_files:
            _console.print(f"     - {gjf}")
    _console.print("\n  2. After calculations complete, the workflow will continue automatically")

    return {
        "success": True,
        "message": f"Generated {extracted_count} ESP input files from optimized structures",
        "directories": created_dirs,
        "file_count": extracted_count
    }


def extract_final_geometry_from_log(log_file):
    """
    Extract the final optimized geometry from a Gaussian optimization log file.

    Parameters:
    -----------
    log_file : str
        Path to the Gaussian optimization log file

    Returns:
    --------
    list
        List of strings containing atom coordinates in Gaussian format (e.g., "C 0.0 0.0 0.0")
        Returns None if extraction fails
    """
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()

        # Find the LAST "Standard orientation:" section (final optimized geometry)
        last_std_orient = None
        for i in range(len(lines) - 1, -1, -1):
            if "Standard orientation:" in lines[i]:
                last_std_orient = i
                break

        if last_std_orient is None:
            _console.print(f"[red]✗ Error: Could not find 'Standard orientation:' in {log_file}[/red]")
            return None

        # Coordinates start 5 lines after "Standard orientation:"
        coords_start = last_std_orient + 5

        # Find the end of coordinate section (marked by dashes)
        coords_end = None
        for i in range(coords_start, min(len(lines), coords_start + 200)):
            if "---------------------------------------------------------------------" in lines[i]:
                coords_end = i
                break

        if coords_end is None:
            _console.print(f"  [red]✗[/red] Error: Could not find end of coordinate section in {log_file}")
            return None

        # Parse coordinates
        # Format: Center Number Atomic Number Atomic Type X Y Z
        geometry = []
        element_symbols = {
            1: 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S', 15: 'P', 9: 'F', 17: 'Cl', 35: 'Br', 53: 'I'
        }

        for i in range(coords_start, coords_end):
            parts = lines[i].split()
            if len(parts) >= 6:
                atomic_num = int(parts[1])
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])

                # Get element symbol
                element = element_symbols.get(atomic_num, f"X{atomic_num}")

                geometry.append(f"{element:2s} {x:12.8f} {y:12.8f} {z:12.8f}")

        if not geometry:
            _console.print(f"  [red]✗[/red] Error: No coordinates extracted from {log_file}")
            return None

        _console.print(f"  [green]✓[/green] Extracted [cyan]{len(geometry)}[/cyan] atoms from final optimized geometry")
        return geometry

    except Exception as e:
        _console.print(f"  [red]✗[/red] Error reading {log_file}: {str(e)}")
        return None


def check_for_log_files(amino_acid):
    """
    Check if Gaussian log files exist for the given amino acid.
    Now also checks for optimization logs to support PES-skip workflow.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code

    Returns:
    --------
    dict
        A dictionary containing whether files exist and the paths to any found files
    """
    amino_acid = amino_acid.strip().lower()

    # Define expected PES log file patterns
    pes_patterns = [f"{amino_acid}_ahelix_pes.log", f"{amino_acid}_bsheet_pes.log"]

    # Define optimization log file patterns (for PES-skip workflow)
    opt_patterns = [f"{amino_acid}_ahelix.log", f"{amino_acid}_bsheet.log"]

    # Check for PES files
    found_pes_files = []
    for pattern in pes_patterns:
        matching_files = glob.glob(pattern)
        found_pes_files.extend(matching_files)

    # Check for optimization files
    found_opt_files = []
    for pattern in opt_patterns:
        matching_files = glob.glob(pattern)
        found_opt_files.extend(matching_files)

    # Determine which PES files are missing
    missing_pes_files = [
        pattern for pattern in pes_patterns if pattern not in found_pes_files
    ]

    return {
        "exist": len(found_pes_files) > 0,
        "complete": len(missing_pes_files) == 0,
        "found_files": found_pes_files,
        "missing_files": missing_pes_files,
        # New fields for opt-only workflow
        "opt_files_found": found_opt_files,
        "has_opt_files": len(found_opt_files) > 0,
        "opt_files_complete": len(found_opt_files) == 2,
    }


def run_tleap_for_amino_acid(amino_acid, output_prefix=None):
    """
    Generate and run a tleap input file for a specified amino acid.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code (e.g., ALA, PHE, GLY)
    output_prefix : str, optional
        Prefix for output files. If None, uses 'ace-{amino_acid}-nme'

    Returns:
    --------
    dict
        A dictionary containing:
        - 'success': Boolean indicating if all operations completed successfully
        - 'message': Status message
        - 'files': List of generated files if successful
        - 'error': Error message if unsuccessful
    """
    # Validate input
    amino_acid = amino_acid.strip().upper()
    if len(amino_acid) != 3:
        return {
            "success": False,
            "error": "Invalid amino acid code. Must be exactly 3 letters.",
        }

    # Set output prefix if not provided
    if output_prefix is None:
        output_prefix = f"ace-{amino_acid.lower()}-nme"

    # Create tleap input file
    tleap_input = f"""source leaprc.protein.ff14SB
mol = sequence {{ ACE {amino_acid} NME }}
saveamberparm mol {output_prefix}.parm7 {output_prefix}.rst7
savepdb mol {output_prefix}.pdb
quit
"""

    # Write tleap input to file
    tleap_file = "tleap.in"
    with open(tleap_file, "w") as f:
        f.write(tleap_input)

    try:
        # Run tleap with the input file
        result = subprocess.run(
            ["tleap", "-s", "-f", tleap_file],
            capture_output=True,
            text=True,
            check=True,
        )

        # Check if output files were created
        expected_files = [
            f"{output_prefix}.parm7",
            f"{output_prefix}.rst7",
            f"{output_prefix}.pdb",
        ]
        files_exist = all(os.path.exists(file) for file in expected_files)

        if files_exist:
            return {
                "success": True,
                "message": f"Successfully created files for ACE-{amino_acid}-NME",
                "files": expected_files,
                "stdout": result.stdout,
            }
        else:
            return {
                "success": False,
                "error": "Not all expected output files were created.",
                "stdout": result.stdout,
            }

    except subprocess.CalledProcessError as e:
        return {"success": False, "error": f"Error running tleap: {e.stderr}"}
    except FileNotFoundError:
        return {
            "success": False,
            "error": "tleap command not found. Make sure AmberTools is installed and in your PATH.",
        }


def run_cpptraj_set_angles(amino_acid, parm_file, rst_file):
    """
    Run CPPTRAJ to set phi/psi angles for alpha helix and beta sheet conformations.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    parm_file : str
        Path to the parameter file (.parm7)
    rst_file : str
        Path to the restart file (.rst7)

    Returns:
    --------
    dict
        A dictionary containing:
        - 'success': Boolean indicating if all operations completed successfully
        - 'message': Status message
        - 'files': List of generated files if successful
        - 'error': Error message if unsuccessful
    """
    amino_acid = amino_acid.strip().upper()

    # Create CPPTRAJ input script
    cpptraj_input = f"""parm {parm_file}
loadcrd {rst_file} name res
#alpha-helix
rotatedihedral crdset res value -60 res 2 type phi
rotatedihedral crdset res value -40 res 2 type psi
crdout res {amino_acid.lower()}.ahelix.mol2
#beta-sheet
rotatedihedral crdset res value -143 res 2 type phi
rotatedihedral crdset res value 160 res 2 type psi
crdout res {amino_acid.lower()}.bsheet.mol2
quit
"""

    # Write CPPTRAJ input to file
    cpptraj_file = "cpptraj.in"
    with open(cpptraj_file, "w") as f:
        f.write(cpptraj_input)

    try:
        # Run CPPTRAJ with the input file
        result = subprocess.run(
            ["cpptraj", "-i", cpptraj_file], capture_output=True, text=True, check=True
        )

        # Check if output files were created
        expected_files = [
            f"{amino_acid.lower()}.ahelix.mol2",
            f"{amino_acid.lower()}.bsheet.mol2",
        ]
        files_exist = all(os.path.exists(file) for file in expected_files)

        if files_exist:
            return {
                "success": True,
                "message": f"Successfully created alpha helix and beta sheet conformations for {amino_acid}",
                "files": expected_files,
                "stdout": result.stdout,
            }
        else:
            return {
                "success": False,
                "error": "Not all expected output files were created.",
                "stdout": result.stdout,
            }

    except subprocess.CalledProcessError as e:
        return {"success": False, "error": f"Error running CPPTRAJ: {e.stderr}"}
    except FileNotFoundError:
        return {
            "success": False,
            "error": "cpptraj command not found. Make sure AmberTools is installed and in your PATH.",
        }


def parse_mol2_for_dihedral_atoms(mol2_file):
    """
    Parse a mol2 file to identify atoms involved in phi and psi dihedral angles.

    Parameters:
    -----------
    mol2_file : str
        Path to the mol2 file

    Returns:
    --------
    dict
        A dictionary containing atom indices for phi and psi dihedrals
    """
    try:
        with open(mol2_file, "r") as f:
            content = f.read()

        # Extract the ATOM section
        atom_section_match = re.search(
            r"@<TRIPOS>ATOM\n(.*?)@<TRIPOS>", content, re.DOTALL
        )
        if not atom_section_match:
            return {
                "success": False,
                "error": f"Could not find ATOM section in {mol2_file}",
            }

        atom_section = atom_section_match.group(1)
        atom_lines = atom_section.strip().split("\n")

        # Dictionary to store atom index -> (atom_name, residue_name)
        atoms = {}
        # Also store a mapping of residue_name -> list of (atom_idx, atom_name) pairs
        residue_atoms = defaultdict(list)

        for line in atom_lines:
            parts = line.split()
            if len(parts) >= 8:
                atom_idx = int(parts[0])
                atom_name = parts[1]
                residue_name = parts[7]
                atoms[atom_idx] = (atom_name, residue_name)
                residue_atoms[residue_name].append((atom_idx, atom_name))

        # Debug output
        _console.print(f"[grey50]Found residues: {list(residue_atoms.keys())}[/grey50]")
        for res, atoms_list in residue_atoms.items():
            _console.print(f"[grey50]Residue {res} has {len(atoms_list)} atoms[/grey50]")
            _console.print(
                f"[grey50]  First 5 atoms: {atoms_list[:5] if len(atoms_list) > 5 else atoms_list}[/grey50]"
            )

        # Phi dihedral atoms: C(ACE) - N(XXX) - CA(XXX) - C(XXX)
        # Psi dihedral atoms: N(XXX) - CA(XXX) - C(XXX) - N(NME)

        # Initialize variables to store atom indices
        c_ace = None
        n_mid = None
        ca_mid = None
        c_mid = None  # Carbonyl carbon of middle residue
        n_nme = None

        # Find atoms in ACE residue
        for atom_idx, atom_name in residue_atoms.get("ACE", []):
            if atom_name == "C":  # Usually labeled as 'C' in ACE
                c_ace = atom_idx
                _console.print(f"[grey50]Found ACE carbonyl carbon: atom {atom_idx} named {atom_name}[/grey50]")

        # Find atoms in middle (amino acid) residue - look for backbone atoms
        for atom_idx, atom_name in residue_atoms.get(
            "", []
        ):  # Sometimes middle residue has empty name
            if atom_name == "N":
                n_mid = atom_idx
                _console.print(f"[grey50]Found middle residue N: atom {atom_idx}[/grey50]")
            elif atom_name == "CA":
                ca_mid = atom_idx
                _console.print(f"[grey50]Found middle residue CA: atom {atom_idx}[/grey50]")
            elif atom_name == "C":  # Carbonyl carbon
                c_mid = atom_idx
                _console.print(f"[grey50]Found middle residue C: atom {atom_idx}[/grey50]")

        # Check other potential residue names if we didn't find our atoms
        if not (n_mid and ca_mid and c_mid):
            # Sometimes the middle residue has the amino acid code as the residue name
            for residue_name, atoms_list in residue_atoms.items():
                if residue_name not in ["ACE", "NME"]:  # Skip N and C terminus caps
                    for atom_idx, atom_name in atoms_list:
                        if not n_mid and atom_name == "N":
                            n_mid = atom_idx
                            _console.print(
                                f"[grey50]Found middle residue N in {residue_name}: atom {atom_idx}[/grey50]"
                            )
                        elif not ca_mid and atom_name == "CA":
                            ca_mid = atom_idx
                            _console.print(
                                f"[grey50]Found middle residue CA in {residue_name}: atom {atom_idx}[/grey50]"
                            )
                        elif not c_mid and atom_name == "C":
                            c_mid = atom_idx
                            _console.print(
                                f"[grey50]Found middle residue C in {residue_name}: atom {atom_idx}[/grey50]"
                            )

        # Find atoms in NME residue
        for atom_idx, atom_name in residue_atoms.get("NME", []):
            if atom_name == "N":
                n_nme = atom_idx
                _console.print(f"[grey50]Found NME nitrogen: atom {atom_idx}[/grey50]")

        # Check if we're still missing atoms - look for backbone atoms by patterns
        if not c_ace:
            # Look for atoms named 'C' in any residue that isn't the middle or NME
            for atom_idx, (atom_name, residue_name) in atoms.items():
                if atom_name == "C" and residue_name == "ACE":
                    c_ace = atom_idx
                    _console.print(f"[grey50]Found ACE C by pattern matching: atom {atom_idx}[/grey50]")

        if not n_nme:
            # Look for atoms named 'N' in any residue that isn't the middle or ACE
            for atom_idx, (atom_name, residue_name) in atoms.items():
                if atom_name == "N" and residue_name == "NME":
                    n_nme = atom_idx
                    _console.print(f"[grey50]Found NME N by pattern matching: atom {atom_idx}[/grey50]")

        # If we still can't find the middle residue atoms, try looking by atom name patterns
        if not (n_mid and ca_mid and c_mid):
            # Some mol2 files might use different naming schemes
            for atom_idx, (atom_name, residue_name) in atoms.items():
                if residue_name not in ["ACE", "NME"]:
                    if not n_mid and atom_name.startswith("N"):
                        n_mid = atom_idx
                        _console.print(
                            f"[grey50]Found middle N by pattern: atom {atom_idx} named {atom_name}[/grey50]"
                        )
                    elif not ca_mid and atom_name == "CA":
                        ca_mid = atom_idx
                        _console.print(
                            f"[grey50]Found middle CA by pattern: atom {atom_idx} named {atom_name}[/grey50]"
                        )
                    elif (
                        not c_mid
                        and atom_name.startswith("C")
                        and atom_name != "CA"
                        and atom_name != "CB"
                    ):
                        c_mid = atom_idx
                        _console.print(
                            f"[grey50]Found middle C by pattern: atom {atom_idx} named {atom_name}[/grey50]"
                        )

        # Check if all atoms were found
        missing_atoms = []
        if not c_ace:
            missing_atoms.append("C(ACE)")
        if not n_mid:
            missing_atoms.append("N(middle)")
        if not ca_mid:
            missing_atoms.append("CA(middle)")
        if not c_mid:
            missing_atoms.append("C(middle)")
        if not n_nme:
            missing_atoms.append("N(NME)")

        if missing_atoms:
            return {
                "success": False,
                "error": f"Could not find all dihedral atoms. Missing: {', '.join(missing_atoms)}",
            }

        # Construct phi and psi dihedral atom lists
        phi_atoms = [c_ace, n_mid, ca_mid, c_mid]
        psi_atoms = [n_mid, ca_mid, c_mid, n_nme]

        _console.print(f"[green]✓[/green] Phi dihedral atoms: [cyan]{phi_atoms}[/cyan]")
        _console.print(f"[green]✓[/green] Psi dihedral atoms: [cyan]{psi_atoms}[/cyan]")

        return {"success": True, "phi": phi_atoms, "psi": psi_atoms}

    except Exception as e:
        import traceback

        traceback.print_exc()
        return {"success": False, "error": f"Error parsing mol2 file: {str(e)}"}


def _mol2_atom_names(mol2_file):
    """Ordered atom names from a mol2's @<TRIPOS>ATOM section (1-based position)."""
    names = []
    try:
        with open(mol2_file) as f:
            content = f.read()
        m = re.search(r"@<TRIPOS>ATOM\n(.*?)@<TRIPOS>", content, re.DOTALL)
        if not m:
            return names
        for line in m.group(1).strip().split("\n"):
            parts = line.split()
            if len(parts) >= 2:
                names.append(parts[1])
    except Exception:
        return []
    return names


def create_gaussian_input(
    mol2_file,
    output_file,
    memory,
    procs,
    keywords,
    title,
    phi_atoms,
    psi_atoms,
    phi_value,
    psi_value,
    charge=0,
    multiplicity=1,
    scan_dihedral=None,
):
    """
    Create a Gaussian input file from a mol2 file with fixed dihedral constraints.

    ``scan_dihedral``, when given, is ``(quad, n_steps, step_size)`` where ``quad``
    is a 4-tuple of 1-based mol2 atom indices. It appends a relaxed-scan line
    ``D i j k l S n step`` after the φ/ψ ``F`` restraints, so the de-novo route can
    AUTO-generate a sidechain torsion scan (φ/ψ held at the backbone conformer)
    instead of only scaffolding it for hand-editing — unifying it with the
    from-structure route, whose scan the parameterizer already writes.

    Parameters:
    -----------
    mol2_file : str
        Path to the mol2 file
    output_file : str
        Path to the output Gaussian input file
    memory : str
        Memory specification for Gaussian (e.g., "8GB")
    procs : int
        Number of processors to use
    keywords : str
        Gaussian keywords/options line
    title : str
        Title for the Gaussian job
    phi_atoms : list
        List of 4 atom indices for the phi dihedral
    psi_atoms : list
        List of 4 atom indices for the psi dihedral
    phi_value : float
        Value for the phi dihedral angle
    psi_value : float
        Value for the psi dihedral angle
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity

    Returns:
    --------
    dict
        A dictionary indicating success or failure
    """
    try:
        # Read coordinates from mol2 file
        with open(mol2_file, "r") as f:
            content = f.read()

        # Extract the ATOM section
        atom_section_match = re.search(
            r"@<TRIPOS>ATOM\n(.*?)@<TRIPOS>", content, re.DOTALL
        )
        if not atom_section_match:
            return {
                "success": False,
                "error": f"Could not find ATOM section in {mol2_file}",
            }

        atom_section = atom_section_match.group(1)
        atom_lines = atom_section.strip().split("\n")

        # Mapping from mol2 atom types to element symbols
        # This dictionary maps common Amber atom types to their corresponding elements
        atom_type_to_element = {
            # Carbon types
            "C": "C",
            "CA": "C",
            "CB": "C",
            "CC": "C",
            "CD": "C",
            "CE": "C",
            "CF": "C",
            "CG": "C",
            "CH": "C",
            "CI": "C",
            "CJ": "C",
            "CK": "C",
            "CL": "C",
            "CM": "C",
            "CN": "C",
            "CO": "C",
            "CP": "C",
            "CQ": "C",
            "CR": "C",
            "CS": "C",
            "CT": "C",
            "CU": "C",
            "CV": "C",
            "CW": "C",
            "CX": "C",
            "CY": "C",
            "CZ": "C",
            "C0": "C",
            # Hydrogen types
            "H": "H",
            "HA": "H",
            "HB": "H",
            "HC": "H",
            "HD": "H",
            "HE": "H",
            "HF": "H",
            "HG": "H",
            "HH": "H",
            "HI": "H",
            "HJ": "H",
            "HK": "H",
            "HL": "H",
            "HM": "H",
            "HN": "H",
            "HO": "H",
            "HP": "H",
            "HQ": "H",
            "HR": "H",
            "HS": "H",
            "HT": "H",
            "HU": "H",
            "HV": "H",
            "HW": "H",
            "HX": "H",
            "HY": "H",
            "HZ": "H",
            "H1": "H",
            "H2": "H",
            "H3": "H",
            "H4": "H",
            "H5": "H",
            # Nitrogen types
            "N": "N",
            "NA": "N",
            "NB": "N",
            "NC": "N",
            "ND": "N",
            "NE": "N",
            "NF": "N",
            "NG": "N",
            "NH": "N",
            "NI": "N",
            "NJ": "N",
            "NK": "N",
            "NL": "N",
            "NM": "N",
            "NN": "N",
            "NO": "N",
            # Oxygen types
            "O": "O",
            "OA": "O",
            "OB": "O",
            "OC": "O",
            "OD": "O",
            "OE": "O",
            "OF": "O",
            "OG": "O",
            "OH": "O",
            "OI": "O",
            "OJ": "O",
            "OK": "O",
            "OL": "O",
            "OM": "O",
            "ON": "O",
            "OO": "O",
            "OW": "O",
            "O2": "O",
            # Sulfur types
            "S": "S",
            "SA": "S",
            "SB": "S",
            "SC": "S",
            "SD": "S",
            "SE": "S",
            "SF": "S",
            "SG": "S",
            "SH": "S",
            "SI": "S",
            "SJ": "S",
            "SK": "S",
            "SL": "S",
            "SM": "S",
            "SN": "S",
            "SO": "S",
            # Phosphorus types
            "P": "P",
            "PA": "P",
            "PB": "P",
            "PC": "P",
            "PD": "P",
            "PE": "P",
            "PF": "P",
            "PG": "P",
            "PH": "P",
            "PI": "P",
            "PJ": "P",
            "PK": "P",
            "PL": "P",
            "PM": "P",
            "PN": "P",
            "PO": "P",
            # Other common elements
            "Br": "Br",
            "BR": "Br",
            "br": "Br",
            "Cl": "Cl",
            "CL": "Cl",
            "cl": "Cl",
            "F": "F",
            "FL": "F",
            "fl": "F",
            "I": "I",
            "ID": "I",
            "id": "I",
            # Special case for wildcards
            "X": "X",
        }

        # Function to determine element from atom_name and atom_type
        def get_element(atom_name, atom_type=None):
            # First, try using the atom_type if available
            if atom_type and atom_type in atom_type_to_element:
                return atom_type_to_element[atom_type]

            # Next, check if the first one or two characters indicate the element
            if len(atom_name) >= 1:
                first_char = atom_name[0]
                if first_char in ["C", "H", "N", "O", "S", "P", "F", "I"]:
                    return first_char
                if len(atom_name) >= 2:
                    first_two = atom_name[0:2]
                    if first_two in [
                        "Br",
                        "Cl",
                        "Si",
                        "Na",
                        "Mg",
                        "Al",
                        "Fe",
                        "Zn",
                        "Cu",
                    ]:
                        return first_two

            # Default to carbon if we can't determine (but we should log this issue)
            _console.print(
                f"[yellow]⚠[/yellow] Warning: Could not determine element for atom {atom_name}, defaulting to Carbon"
            )
            return "C"

        # Prepare Gaussian input format
        gaussian_coords = []
        for line in atom_lines:
            parts = line.split()
            if len(parts) >= 6:
                atom_idx = int(parts[0])
                atom_name = parts[1]
                x, y, z = parts[2:5]

                # Get atom type from mol2 file if available (typically column 6)
                atom_type = parts[5] if len(parts) > 5 else None

                # Determine element
                element = get_element(atom_name, atom_type)

                gaussian_coords.append(f"{element}  {x}  {y}  {z}")

        # Create Gaussian input file
        with open(output_file, "w") as f:
            f.write(f"%mem={memory}\n")
            f.write(f"%nprocshared={procs}\n")
            f.write(f"#p {keywords}\n\n")
            f.write(f"{title}\n\n")
            f.write(
                f"{charge} {multiplicity}\n"
            )  # Use specified charge and multiplicity

            # Write coordinates
            for coord in gaussian_coords:
                f.write(f"{coord}\n")

            # Add blank line and dihedral constraints
            f.write("\n")
            f.write(
                f"D {phi_atoms[0]} {phi_atoms[1]} {phi_atoms[2]} {phi_atoms[3]} F\n"
            )
            f.write(
                f"D {psi_atoms[0]} {psi_atoms[1]} {psi_atoms[2]} {psi_atoms[3]} F\n"
            )
            if scan_dihedral:
                quad, n_steps, step_size = scan_dihedral
                a1, a2, a3, a4 = quad
                f.write(f"D {a1} {a2} {a3} {a4} S {int(n_steps)} {float(step_size):.1f}\n")
            f.write("\n")

        return {"success": True, "file": output_file, "scan": bool(scan_dihedral)}

    except Exception as e:
        return {
            "success": False,
            "error": f"Error creating Gaussian input file: {str(e)}",
        }


def mol2_to_gaussian(amino_acid, mol2_files, gaussian_params, charge=0, multiplicity=1,
                     scan_dihedral=None):
    """
    Convert mol2 files to Gaussian input files with fixed dihedrals.

    ``scan_dihedral`` (``(quad, n_steps, step_size)``, 1-based mol2 indices), when
    given, auto-generates a relaxed sidechain-torsion scan per backbone conformer
    and names the outputs ``{aa}_{conf}_pes.gjf`` so the existing PES-scan
    detection (step 4/5) picks them up. Without it, the legacy φ/ψ-frozen opt
    scaffold (``{aa}_{conf}.gjf``) is written unchanged.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    mol2_files : list
        List of mol2 files to convert
    gaussian_params : dict
        Parameters for Gaussian input files
    charge : int
        Molecular charge
    multiplicity : int
        Spin multiplicity

    Returns:
    --------
    dict
        A dictionary containing success status and created files
    """
    created_files = []

    for mol2_file in mol2_files:
        # Determine conformation type from filename
        if "ahelix" in mol2_file:
            conformation = "ahelix"
            phi_value = -60
            psi_value = -40
        elif "bsheet" in mol2_file:
            conformation = "bsheet"
            phi_value = -143
            psi_value = 160
        else:
            return {
                "success": False,
                "error": f"Could not determine conformation type from filename: {mol2_file}",
            }

        # Parse mol2 file to get dihedral atoms
        parse_result = parse_mol2_for_dihedral_atoms(mol2_file)
        if not parse_result["success"]:
            return parse_result

        # Create Gaussian input file. A scan is named *_pes.gjf so the existing
        # PES detection routes it through conformer extraction.
        suffix = "_pes" if scan_dihedral else ""
        output_file = f"{amino_acid.lower()}_{conformation}{suffix}.gjf"
        gaussian_result = create_gaussian_input(
            mol2_file,
            output_file,
            gaussian_params["memory"],
            gaussian_params["procs"],
            gaussian_params["keywords"],
            f"{gaussian_params['title']} - {amino_acid} {conformation}",
            parse_result["phi"],
            parse_result["psi"],
            phi_value,
            psi_value,
            charge,
            multiplicity,
            scan_dihedral=scan_dihedral,
        )

        if not gaussian_result["success"]:
            return gaussian_result

        created_files.append(output_file)

    return {
        "success": True,
        "message": f"Successfully created Gaussian input files for {amino_acid}",
        "files": created_files,
    }


# Peptide backbone + cap atom names frozen during the from-structure QM opt: the
# backbone and ACE/NME caps stay at crystallographic geometry while the side
# chain and the covalently attached modification relax.
_MODAA_BACKBONE_FREEZE = {"N", "H", "HN", "CA", "HA", "HA2", "HA3", "C", "O"}


def capped_pdb_to_gaussian(pdb_file, output_gjf, gaussian_params, charge=0,
                           multiplicity=1, *, freeze_keys=None,
                           dihedral_restraints=None, chk=None,
                           scan_dihedral=None, coords_override=None):
    """Write a Gaussian opt input from an ACE/NME-capped model compound.

    Atoms are always emitted in ACE - middle - NME order (binned by residue
    name), so the sep-bond detector's first-6 / last-6 assumption holds no
    matter where reduce or interactive edits placed the hydrogens.

    Geometry constraints (Route B offers a choice; see the step-2 chooser):
      - freeze_keys: iterable of ``(resseq:int, name:str)`` atom keys to hold at
        crystal coordinates (Cartesian freeze flag -1). ``None`` keeps the legacy
        default (peptide backbone + ACE/NME caps frozen). An empty set freezes
        nothing (used when restraining dihedrals instead).
      - dihedral_restraints: list of 4-tuples of ``(resseq, name)`` atom keys;
        each becomes a ``D i j k l F`` line in an appended ModRedundant section
        (1-based indices in the emitted order), freezing that dihedral at its
        crystal value. The caller must put ``Opt=ModRedundant`` in the route.
      - chk: optional ``%chk`` filename (needed so the Hessian is saved for
        Seminario refinement).
      - scan_dihedral: optional ``(quad, n_steps, step_size)`` where ``quad`` is
        a 4-tuple of ``(resseq, name)`` keys. Emits a relaxed scan line
        ``D i j k l S n step`` so the torsion is driven through ``n_steps`` of
        ``step_size`` degrees while everything else (including any φ/ψ restraint
        lines above) relaxes. The scanned dihedral is automatically excluded from
        the ``F`` restraint set (cannot freeze and scan the same torsion). The
        caller must put ``Opt=ModRedundant`` in the route and omit ``Freq``.
      - coords_override: optional ``[(x, y, z), ...]`` replacing the PDB's
        coordinates, in the EMITTED (ACE - middle - NME) order. Lets a caller
        re-use this model's atom names, ordering and name-keyed restraint/scan
        machinery while supplying a different geometry — e.g. the step-3 scan
        starting from step 2's optimized geometry, read back from the opt log
        (which is in that same order, because the opt route carries NoSymm).
        Rejected if its length does not match the model's atom count.

    Returns a dict with success, file, frozen/total atom counts, the number of
    dihedral restraints written, and whether a scan line was emitted.
    """
    ace, mid, nme = [], [], []  # each: dict(key, resname, name, element, x, y, z)
    try:
        with open(pdb_file) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                name = line[12:16].strip()
                resname = line[17:20].strip()
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    resseq = 0
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                element = line[76:78].strip()
                if not element:
                    element = name[0] if name[:1].isalpha() else "C"
                rec = {"key": (resseq, name), "resname": resname, "name": name,
                       "element": element, "x": x, "y": y, "z": z}
                (ace if resname == "ACE" else nme if resname == "NME" else mid).append(rec)
    except Exception as e:
        return {"success": False, "error": f"Error reading {pdb_file}: {e}"}

    atoms = ace + mid + nme
    if not atoms:
        return {"success": False, "error": f"No atoms parsed from {pdb_file}"}

    if coords_override is not None:
        # Mismatched length means the override is not this model's geometry;
        # writing it would silently pair coordinates with the wrong atoms.
        if len(coords_override) != len(atoms):
            return {"success": False,
                    "error": (f"Geometry override has {len(coords_override)} atoms but "
                              f"{os.path.basename(pdb_file)} has {len(atoms)}")}
        for a, (x, y, z) in zip(atoms, coords_override):
            a["x"], a["y"], a["z"] = float(x), float(y), float(z)

    # 1-based emit-order index of each atom (first wins on duplicate keys) for
    # resolving freeze/dihedral atom keys.
    index_of = {}
    for i, a in enumerate(atoms, start=1):
        index_of.setdefault(a["key"], i)

    # Cartesian freeze flags.
    if freeze_keys is None:
        for a in atoms:
            a["flag"] = -1 if (a["resname"] in ("ACE", "NME")
                               or a["name"] in _MODAA_BACKBONE_FREEZE) else 0
    else:
        fk = set(freeze_keys)
        for a in atoms:
            a["flag"] = -1 if a["key"] in fk else 0
    frozen_n = sum(1 for a in atoms if a["flag"] == -1)

    # The scanned dihedral (if any) resolves to 1-based indices and is excluded
    # from the F restraint set — you cannot freeze and scan the same torsion.
    scan_line = None
    scan_quad_key = None
    if scan_dihedral:
        quad, n_steps, step_size = scan_dihedral
        scan_quad_key = tuple(quad)
        idxs = [index_of.get(k) for k in quad]
        if all(idxs):
            scan_line = ("D " + " ".join(str(i) for i in idxs)
                         + f" S {int(n_steps)} {float(step_size):.1f}")

    # ModRedundant dihedral-freeze lines (drop any that reference a missing atom
    # or that coincide with the scanned dihedral).
    dih_lines = []
    for quad in (dihedral_restraints or []):
        if scan_quad_key is not None and tuple(quad) == scan_quad_key:
            continue
        idxs = [index_of.get(k) for k in quad]
        if all(idxs):
            dih_lines.append("D " + " ".join(str(i) for i in idxs) + " F")

    try:
        with open(output_gjf, "w") as f:
            if chk:
                f.write(f"%chk={chk}\n")
            f.write(f"%mem={gaussian_params['memory']}\n")
            f.write(f"%nprocshared={gaussian_params['procs']}\n")
            f.write(f"#p {gaussian_params['keywords']}\n\n")
            f.write(f"{gaussian_params['title']}\n\n")
            f.write(f"{charge} {multiplicity}\n")
            for a in atoms:
                # Gaussian reads the integer freeze flag before the coordinates
                # (0 = optimize, -1 = frozen) when the route contains `opt`.
                f.write(f"{a['element']:<2s} {a['flag']:2d}  "
                        f"{a['x']:12.6f} {a['y']:12.6f} {a['z']:12.6f}\n")
            f.write("\n")
            if dih_lines or scan_line:
                for ln in dih_lines:
                    f.write(ln + "\n")
                if scan_line:
                    f.write(scan_line + "\n")
                f.write("\n")
    except Exception as e:
        return {"success": False, "error": f"Error writing {output_gjf}: {e}"}

    return {"success": True, "file": output_gjf, "frozen_atoms": frozen_n,
            "total_atoms": len(atoms), "dihedral_restraints": len(dih_lines),
            "scan_written": scan_line is not None}


def _parse_capped_model(pdb_file):
    """Parse an ACE/NME-capped model PDB into per-residue atom info.

    Returns a dict:
      ``residues`` – ordered list of {resseq, resname, names(set)}
      ``ace``/``nme``/``aa`` – the resseq of the ACE cap, NME cap, and the
      amino-acid residue (the non-cap residue carrying a full N/CA/C backbone),
      or None if absent.
    Used to derive backbone φ/ψ (and, later, linkage) dihedral atom keys and to
    present atom lists for the constraint chooser.
    """
    residues = {}
    order = []
    try:
        with open(pdb_file) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                name = line[12:16].strip()
                resname = line[17:20].strip()
                try:
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                if resseq not in residues:
                    residues[resseq] = {"resseq": resseq, "resname": resname, "names": set()}
                    order.append(resseq)
                residues[resseq]["names"].add(name)
    except Exception:
        return {"residues": [], "ace": None, "nme": None, "aa": None}

    ace = next((r for r in order if residues[r]["resname"] == "ACE"), None)
    nme = next((r for r in order if residues[r]["resname"] == "NME"), None)
    aa = next((r for r in order if residues[r]["resname"] not in ("ACE", "NME")
               and {"N", "CA", "C"} <= residues[r]["names"]), None)
    return {"residues": [residues[r] for r in order], "ace": ace, "nme": nme, "aa": aa}


def phi_psi_restraint_keys(pdb_file):
    """Return ``[phi, psi]`` dihedral atom-key quads for a capped model, or ``[]``.

    φ = ACE:C – AA:N – AA:CA – AA:C ; ψ = AA:N – AA:CA – AA:C – NME:N.
    Each atom key is ``(resseq, name)``.
    """
    m = _parse_capped_model(pdb_file)
    if not (m["ace"] and m["nme"] and m["aa"]):
        return []
    ac, nm, a = m["ace"], m["nme"], m["aa"]
    phi = [(ac, "C"), (a, "N"), (a, "CA"), (a, "C")]
    psi = [(a, "N"), (a, "CA"), (a, "C"), (nm, "N")]
    return [phi, psi]


# Covalent radii (Å) for distance-based bond perception in the capped model.
_COVALENT_RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05,
                   "P": 1.07, "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39,
                   "SE": 1.20, "ZN": 1.22, "FE": 1.32}


def _parse_capped_model_atoms(pdb_file):
    """Parse a capped model PDB into ordered atoms + CONECT pairs.

    Returns {"atoms": [{serial, resseq, name, resname, element, x, y, z}, ...],
             "conect": [(serial, serial), ...]}. reduce renumbers and re-emits
    the covalent CONECT, so the serials here match this file.
    """
    atoms, conect = [], []
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        serial = int(line[6:11])
                    except ValueError:
                        serial = len(atoms) + 1
                    name = line[12:16].strip()
                    resname = line[17:20].strip()
                    try:
                        resseq = int(line[22:26])
                    except ValueError:
                        resseq = 0
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    element = line[76:78].strip()
                    if not element:
                        element = name[0] if name[:1].isalpha() else "C"
                    atoms.append({"serial": serial, "resseq": resseq, "name": name,
                                  "resname": resname, "element": element.upper(),
                                  "x": x, "y": y, "z": z})
                elif line.startswith("CONECT"):
                    nums = line[6:].split()
                    try:
                        a = int(nums[0])
                        for b in nums[1:]:
                            conect.append((a, int(b)))
                    except (ValueError, IndexError):
                        pass
    except Exception:
        return {"atoms": [], "conect": []}
    return {"atoms": atoms, "conect": conect}


def _write_capped_model_pdb(atoms, conect_pairs, out_pdb):
    """Write a capped model PDB in the fixed-column layout the from-structure
    steps consume (blank altLoc column, chain A). ``atoms`` are dicts with
    name/resname/resseq/element/x/y/z; serials are re-assigned 1..N in order."""
    def _name_field(name):
        return name.ljust(4) if len(name) >= 4 else (" " + name).ljust(4)

    with open(out_pdb, "w") as f:
        for i, a in enumerate(atoms, start=1):
            f.write(
                "ATOM  {serial:5d} {name} {resname:>3s} A{resid:4d}    "
                "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n".format(
                    serial=i, name=_name_field(a["name"]),
                    resname=a["resname"], resid=a.get("resseq", a.get("resid", 0)),
                    x=a["x"], y=a["y"], z=a["z"],
                    el=(a.get("element") or a["name"][:1])[:2],
                )
            )
        for s1, s2 in conect_pairs:
            f.write(f"CONECT{s1:5d}{s2:5d}\n")
            f.write(f"CONECT{s2:5d}{s1:5d}\n")
        f.write("END\n")


def _reconcile_conformer_to_roster(ref_atoms, ref_conect, reduced_pdb, out_pdb,
                                   console=None):
    """Force a secondary conformer onto the reference atom roster.

    Multi-conformer RESP requires every conformer to expose the identical atom
    set in the identical order (one fitted charge per atom, shared across
    conformers). The reference conformer — the one the user curated — is the
    canonical roster. Here we take ``reduced_pdb`` (this conformer's heavy atoms
    plus reduce-placed hydrogens at its own geometry) and emit ``out_pdb`` with
    exactly the reference's atoms, in the reference's order and names, but with
    this conformer's coordinates (matched by atom name). Consequences:

    * A hydrogen the user removed on the reference is absent everywhere (it is
      simply not in the roster we iterate).
    * A hydrogen the reference kept is taken from this conformer's reduce output.
    * A reference atom missing from this conformer's reduce output (rare — a
      geometry-dependent tautomer that reduce placed differently) keeps the
      reference coordinate, with a warning, so the roster still matches exactly
      rather than silently diverging and breaking the fit.

    Returns (out_pdb, missing_names).
    """
    reduced = _parse_capped_model_atoms(reduced_pdb)["atoms"]
    # Key on (resname, resseq, name), not bare name: the ACE/NME caps reuse
    # standard backbone names (C, O, N, CH3) that also occur in the middle
    # residue, so a name-only match would cross-wire cap and body atoms.
    def _key(a):
        return (a.get("resname"), a.get("resseq"), a["name"])

    by_key = {}
    for a in reduced:
        by_key.setdefault(_key(a), a)

    out_atoms, missing = [], []
    for ra in ref_atoms:
        src = by_key.get(_key(ra))
        if src is None:
            missing.append(ra["name"])
            x, y, z = ra["x"], ra["y"], ra["z"]
        else:
            x, y, z = src["x"], src["y"], src["z"]
        out_atoms.append({"name": ra["name"], "resname": ra["resname"],
                          "resseq": ra.get("resseq", 0),
                          "element": ra.get("element", ra["name"][:1]),
                          "x": x, "y": y, "z": z})

    _write_capped_model_pdb(out_atoms, ref_conect, out_pdb)
    if missing and console is not None:
        console.print(
            f"[yellow]⚠ {len(missing)} reference atom(s) were not reproduced by reduce "
            f"on this conformer ({', '.join(missing)}); kept the reference coordinate so "
            f"the atom roster still matches for RESP. Verify placement in the viewer.[/yellow]")
    return out_pdb, missing


def _bonded(a, b):
    """Distance-based bond test using covalent radii (1.3x tolerance)."""
    d = ((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2) ** 0.5
    ra = _COVALENT_RADII.get(a["element"], 0.77)
    rb = _COVALENT_RADII.get(b["element"], 0.77)
    return 0.4 < d < (ra + rb) * 1.3


def _heavy_neighbors(atoms, idx):
    """Indices of heavy atoms bonded (by distance) to atoms[idx]."""
    return [j for j, other in enumerate(atoms)
            if j != idx and other["element"] != "H" and _bonded(atoms[idx], other)]


def _linkage_dihedral_candidates(model):
    """Candidate dihedrals spanning each covalent CONECT bond of the capped model.

    For a covalent bond B–C, a dihedral X–B–C–Y uses a heavy neighbour X of B
    (≠C) and a heavy neighbour Y of C (≠B). Returns a list of
    {"label": str, "keys": [(resseq, name) x4]} (mirror duplicates removed).
    """
    atoms = model["atoms"]
    by_serial = {a["serial"]: i for i, a in enumerate(atoms)}
    cands, seen = [], set()
    for s1, s2 in model["conect"]:
        i, j = by_serial.get(s1), by_serial.get(s2)
        if i is None or j is None:
            continue
        for x in _heavy_neighbors(atoms, i):
            if x == j:
                continue
            for y in _heavy_neighbors(atoms, j):
                if y == i or y == x:
                    continue
                quad = (x, i, j, y)
                canon = min(quad, tuple(reversed(quad)))
                if canon in seen:
                    continue
                seen.add(canon)
                cands.append({
                    "keys": [(atoms[k]["resseq"], atoms[k]["name"]) for k in quad],
                    "label": "–".join(atoms[k]["name"] for k in quad),
                })
    return cands


def _atom_label_map(atoms):
    """Display labels + a resolver for a capped model's atoms.

    Returns ``(labels, lookup)``: ``labels`` is the ordered display list (bare
    name when unique, else ``name@resseq``); ``lookup`` maps every label — and
    the always-unambiguous ``name@resseq`` form — to the ``(resseq, name)`` key.
    """
    from collections import Counter
    counts = Counter(a["name"] for a in atoms)
    labels, lookup = [], {}
    for a in atoms:
        qualified = f'{a["name"]}@{a["resseq"]}'
        lab = a["name"] if counts[a["name"]] == 1 else qualified
        labels.append(lab)
        lookup[lab] = (a["resseq"], a["name"])
        lookup[qualified] = (a["resseq"], a["name"])
    return labels, lookup


def _describe_constraint(mode, dihedrals, freeze_keys):
    """One-line human summary of the geometry constraint for the step-3 panel."""
    if mode == "restrain":
        return (f"φ/ψ torsion restraints to the crystal fold "
                f"({len(dihedrals or [])} dihedral(s) frozen, Opt=ModRedundant)")
    if mode == "freeze":
        return f"{len(freeze_keys or [])} atom(s) frozen at crystal coordinates (Cartesian)"
    return "none — full relaxation"


def _describe_dihedral(quad, model=None):
    """Render a dihedral (4 ``(resseq, name)`` keys) as ``A–B–C–D`` atom names."""
    return "–".join(str(k[1]) for k in (quad or []))


def _read_gjf_elements(gjf_file):
    """Element symbols, in order, from a Gaussian input's coordinate block."""
    elements = []
    try:
        with open(gjf_file) as f:
            lines = f.readlines()
    except Exception:
        return elements
    start = None
    for idx, ln in enumerate(lines):
        parts = ln.split()
        if len(parts) == 2 and all(p.lstrip("+-").isdigit() for p in parts):
            start = idx + 1  # line after "charge multiplicity"
            break
    if start is None:
        return elements
    for ln in lines[start:]:
        if not ln.strip():
            break
        elements.append(ln.split()[0])
    return elements


def _extract_last_orientation(log_file, n_atoms):
    """Coordinates [(x,y,z), ...] of the last orientation block in a Gaussian log."""
    blocks = []
    try:
        with open(log_file) as f:
            lines = f.readlines()
    except Exception:
        return None
    i = 0
    while i < len(lines):
        if ("Standard orientation:" in lines[i]) or ("Input orientation:" in lines[i]):
            j = i + 5  # skip the 4 header/underline lines to the first data row
            block = []
            while j < len(lines) and not lines[j].strip().startswith("---"):
                parts = lines[j].split()
                if len(parts) >= 6:
                    try:
                        block.append((float(parts[3]), float(parts[4]), float(parts[5])))
                    except ValueError:
                        break
                else:
                    break
                j += 1
            if block:
                blocks.append(block)
            i = j
        else:
            i += 1
    if not blocks:
        return None
    return blocks[-1]


def _extract_optimized_geometry(opt_log, opt_gjf):
    """Zip element symbols (from the opt input) with the last optimized coordinates.

    Returns a list of ``(element, x, y, z)`` in the opt input's atom order (ACE →
    middle → NME, so downstream stays canonical), or ``None`` on failure.
    """
    elements = _read_gjf_elements(opt_gjf)
    if not elements:
        return None
    coords = _extract_last_orientation(opt_log, len(elements))
    if not coords or len(coords) != len(elements):
        return None
    return [(elements[k], coords[k][0], coords[k][1], coords[k][2])
            for k in range(len(elements))]


def write_esp_single_point_gjf(atoms, output_gjf, route, memory, procs, title,
                               charge, multiplicity):
    """Write an ESP single-point Gaussian input from ``(element, x, y, z)`` atoms.

    No opt and no freeze flags — a plain single point at the supplied geometry.
    """
    try:
        with open(output_gjf, "w") as f:
            f.write(f"%mem={memory}\n")
            f.write(f"%nprocshared={procs}\n")
            f.write(f"#p {route}\n\n")
            f.write(f"{title}\n\n")
            f.write(f"{charge} {multiplicity}\n")
            for el, x, y, z in atoms:
                f.write(f"{el:<2s}  {x:12.6f} {y:12.6f} {z:12.6f}\n")
            f.write("\n")
    except Exception as e:
        return {"success": False, "error": f"Error writing {output_gjf}: {e}"}
    return {"success": True, "file": output_gjf, "total_atoms": len(atoms)}


def display_manual_steps_instructions(amino_acid, missing_files):
    """
    Display instructions for manual steps that the user needs to perform.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    missing_files : list
        List of missing log files
    """
    _console.print(Panel(
        f"[bold red]Manual Steps Required for {amino_acid.upper()}[/bold red]\n\n"
        "[bold]The following files need to be created:[/bold]\n" +
        "\n".join(f"  [yellow]•[/yellow] {file}" for file in missing_files) + "\n\n"
        "[bold]Please follow these steps:[/bold]\n"
        "  1. Open the Gaussian input files (.gjf) with GaussView or your preferred editor\n"
        f"     - {amino_acid.lower()}_ahelix.gjf\n"
        f"     - {amino_acid.lower()}_bsheet.gjf\n\n"
        "  2. Modify the structures as needed\n\n"
        "  3. Set up and run PES (Potential Energy Surface) scans in Gaussian\n\n"
        "  4. After Gaussian calculations are complete, run this program again\n"
        "     to proceed with analysis of the PES scan results\n\n"
        "[bold]Note:[/bold] The program will look for these specific log file names when restarted:\n" +
        "\n".join(f"  [cyan]•[/cyan] {file}" for file in missing_files),
        title="Action Required",
        border_style="yellow",
        expand=False
    ))


MODIFIED_AA_STEPS = [
    WorkflowStep(
        id="aa-1", name="Initial Structure",
        description="Generate structure with tleap",
        handler="_checklist_aa_1_structure",
        section="Structure Generation",
    ),
    WorkflowStep(
        id="aa-2", name="Backbone Angles",
        description="Set phi/psi angles with cpptraj",
        handler="_checklist_aa_2_angles",
        section="Structure Generation",
        dependencies=["aa-1"],
    ),
    WorkflowStep(
        id="aa-3", name="Gaussian Input",
        description="Create Gaussian input files for PES scans",
        handler="_checklist_aa_3_gaussian_input",
        section="Structure Generation",
        dependencies=["aa-2"],
    ),
    WorkflowStep(
        id="aa-4", name="PES Scan Results",
        description="Check/process Gaussian PES scan outputs",
        handler="_checklist_aa_4_pes_scan",
        section="QM Calculations",
        dependencies=["aa-3"],
        checkpoint=True,
        checkpoint_message="Run Gaussian PES scans on .gjf files, then resume.",
    ),
    WorkflowStep(
        id="aa-5", name="Extract Structures",
        description="Analyze PES, extract conformers",
        handler="_checklist_aa_5_extract",
        section="QM Calculations",
        dependencies=["aa-4"],
        checkpoint=True,
        checkpoint_message="Run Gaussian on extracted structures, then resume.",
    ),
    WorkflowStep(
        id="aa-6", name="ESP Generation",
        description="Generate ESP from optimized structures",
        handler="_checklist_aa_6_esp",
        section="Parameter Fitting",
        dependencies=["aa-5"],
    ),
    WorkflowStep(
        id="aa-7", name="AC File",
        description="Create AC file from lowest energy structure",
        handler="_checklist_aa_7_ac_file",
        section="Parameter Fitting",
        dependencies=["aa-6"],
    ),
    WorkflowStep(
        id="aa-8", name="RESP Charges",
        description="Run residuegen for RESP charge fitting",
        handler="_checklist_aa_8_resp",
        section="Parameter Fitting",
        dependencies=["aa-7"],
    ),
    WorkflowStep(
        id="aa-9", name="Bonded Parameters",
        description="Generate frcmod with parmchk2",
        handler="_checklist_aa_9_parmchk2",
        section="Parameter Fitting",
        dependencies=["aa-8"],
    ),
    WorkflowStep(
        id="aa-10", name="Force Field Integration",
        description="Create AMBER library and register it for the Topology Generator",
        handler="_checklist_aa_10_library",
        section="Parameter Fitting",
        dependencies=["aa-9"],
    ),
]


# The from-structure route (Route B) reuses the SAME executor handlers — each
# _run_step_* / _checklist_aa_* branches internally on conformer_mode — but the
# de-novo labels ("tleap", "phi/psi cpptraj", "PES scan", "extract conformers")
# describe work Route B never does, so steps 1-6 get accurate names here. Steps
# 7-9 (AC -> RESP -> parmchk2) are identical to Route A and are sliced in
# verbatim rather than duplicated; only the final step differs (FF integration
# vs a bare .lib write).
#
# One QM job per step, which is why steps 2 and 3 are separate:
#   step 2  B3LYP/6-31+G(d) Opt+Freq on the capped model, φ/ψ restrained by
#           default. ALWAYS run. Its checkpoint is the Hessian for the optional
#           Seminario refinement at step 9, so Freq is not negotiable here.
#   step 3  OPTIONAL relaxed torsion scan, started from step 2's OPTIMIZED
#           geometry (so it gates on the opt log). Samples the rotational
#           profile so RESP is not biased to the resolved pose. Skipping it is a
#           first-class choice: _is_scanned() keys off the scan input's presence
#           on disk, so steps 4-6 adapt with no extra state.
#   step 6  HF/6-31G* ESP single point(s) for RESP.
# Steps 3, 4 and 6 are QM gates (opt log, geometry logs, ESP logs).
#
# Route A must keep 2 and 3 separate for a different reason: its step 2 really
# does generate structures (cpptraj sets φ/ψ and writes the mol2 files step 3
# reads).
MODIFIED_AA_FROM_STRUCTURE_STEPS = [
    WorkflowStep(
        id="aa-1", name="Capped Model & Hydrogens",
        description="Add and curate hydrogens on the capped model built from the structure",
        handler="_checklist_aa_1_structure",
        section="Structure Generation",
    ),
    WorkflowStep(
        id="aa-2", name="Geometry Optimization",
        description="Constrain the crystal fold (φ/ψ by default) and write the opt+freq job",
        handler="_checklist_aa_2_angles",
        section="Structure Generation",
        dependencies=["aa-1"],
    ),
    WorkflowStep(
        id="aa-3", name="Dihedral Scan",
        description="Relaxed torsion scan from the optimized geometry — de-biases RESP",
        handler="_checklist_aa_3_gaussian_input",
        section="Structure Generation",
        dependencies=["aa-2"],
        optional=True,
    ),
    WorkflowStep(
        id="aa-4", name="QM Results",
        description="Ingest the opt/scan log(s); validate (imaginary modes / scan profile)",
        handler="_checklist_aa_4_pes_scan",
        section="QM Calculations",
        dependencies=["aa-3"],
        checkpoint=True,
        checkpoint_message="Run the geometry job(s) (opt+freq, and the scan if configured), then resume.",
    ),
    WorkflowStep(
        id="aa-5", name="Select Conformers",
        description="Choose which structures (conformers / scan points) feed the RESP fit",
        handler="_checklist_aa_5_extract",
        section="QM Calculations",
        dependencies=["aa-4"],
    ),
    WorkflowStep(
        id="aa-6", name="ESP Generation",
        description="Write + ingest the HF/6-31G* ESP single point(s) for the selected conformers",
        handler="_checklist_aa_6_esp",
        section="Parameter Fitting",
        dependencies=["aa-5"],
        checkpoint=True,
        checkpoint_message="Run the ESP single-point job(s), then resume.",
    ),
] + MODIFIED_AA_STEPS[6:9] + [
    # Route B's final step is FF integration, not a bare .lib write: split the
    # adduct into AA + cofactor libraries, deposit them into the user library,
    # rename the prepared PDB, wire the Topology-Generator keys, and emit a
    # reuse transformer — the modified-AA analogue of the metal-site integration
    # step. Reuses the shared _checklist_aa_10_library handler, which routes to
    # _run_step_10_from_structure on conformer_mode.
    WorkflowStep(
        id="aa-10", name="Force Field Integration",
        description="Split into AA+cofactor libs, deposit, prepare PDB, emit transformer",
        handler="_checklist_aa_10_library",
        section="Parameter Fitting",
        dependencies=["aa-9"],
    ),
]


class ModifiedAAWorkflowManager:
    """Workflow manager for modified amino acid parameterization.

    Step methods generate structures, run QM gates, and produce AMBER
    parameter files.  Orchestration is handled by WorkflowChecklist;
    the _checklist_aa_* handler methods bridge between the checklist
    and the existing _run_step_N() implementations.
    """

    WORKFLOW_NAME = "Modified Amino Acid Parameterization"
    PAUSE_STATUS = "paused_for_qm"

    def __init__(self, amino_acid, console=None, processor=None,
                 starting_pdb=None, conformer_mode="denovo_parent",
                 source_residues=None, conformer_pdbs=None):
        self.amino_acid = amino_acid
        self.console = console or _console
        self.processor = processor
        self.logger = logging.getLogger(self.__class__.__name__)
        self.step_results = {}
        # Original RedoxSite members ({name, chain_id, resid}) for the
        # from_structure route's FF-integration step (unique naming + prepared-PDB
        # rename + reuse transformer). Empty for the de-novo route.
        self.source_residues = source_residues or []
        # Route selector:
        #   "denovo_parent"  — build ACE-<parent>-NME from the ff14SB library and
        #                      sample alpha/beta backbone conformers (original path)
        #   "from_structure" — use capped model compound(s) extracted from the real
        #                      structure; one crystallographic conformer per
        #                      alternate conformation, no alpha/beta sampling
        #                      (starting_pdb / conformer_pdbs required)
        self.starting_pdb = starting_pdb
        self.conformer_mode = conformer_mode
        # Data-driven conformer list + per-conformer capped PDB map so the
        # from_structure steps loop over however many alternate conformations the
        # crystal supplied. conformer_pdbs is a list of (label, path) tuples; a
        # bare starting_pdb (or none) collapses to the single "xtal" conformer,
        # preserving the original single-conformer behaviour byte-for-byte.
        if conformer_mode == "from_structure":
            if conformer_pdbs:
                self.conformers = [label for label, _ in conformer_pdbs]
                self._conformer_pdb_map = {label: path for label, path in conformer_pdbs}
                self.starting_pdb = conformer_pdbs[0][1]
            else:
                self.conformers = ["xtal"]
                self._conformer_pdb_map = {"xtal": starting_pdb} if starting_pdb else {}
        else:
            self.conformers = ["ahelix", "bsheet"]
            self._conformer_pdb_map = {}

    def set_processor(self, processor):
        """Set processor for workspace access."""
        self.processor = processor

    # ── Steps 1-3: Structure generation ─────────────────────────────

    def _curated_conformer_pdb(self, label):
        """Resolve a conformer to its hydrogen-curated capped PDB.

        Step 1 records the curated path in-memory (``self._conformer_pdb_map``),
        but that link is not persisted: on a fresh session resumed after a QM
        pause the map is rebuilt from the ORIGINAL heavy-atom capped models. So a
        later step (e.g. re-running step 3) would otherwise write QM input from a
        hydrogen-less structure. Here we prefer the on-disk curated file — the
        reference conformer's ``*_H.pdb`` (written by reduce/HydrogenEditor) or a
        secondary's ``*_curated.pdb`` (written by the reconciler) — falling back
        to the in-memory path (already curated, or nothing curated yet).
        """
        raw = self._conformer_pdb_map.get(label)
        if not raw:
            return raw
        base = os.path.splitext(os.path.basename(raw))[0]
        if base.endswith("_H") or base.endswith("_curated"):
            return raw  # already the curated file in memory
        ref_label = self.conformers[0]
        candidate = os.path.abspath(
            f"{base}_H.pdb" if label == ref_label else f"{base}_curated.pdb")
        return candidate if os.path.exists(candidate) else raw

    def _is_scanned(self, label):
        """A conformer is scanned iff its relaxed-scan input exists on disk.

        Derived from disk (not in-memory state) so steps 4-6 survive a
        fresh-session resume after a QM pause without persisting the scan spec.
        """
        return os.path.exists(f"{self.amino_acid.lower()}_{label}_scan.gjf")

    def _parse_scan(self, label):
        """Parse (and cache) a scanned conformer's relaxed-scan log."""
        from proprep.forcefield_prep.pes_scan_refinement import parse_pes_scan_log
        cache = getattr(self, "_scan_cache", None)
        if cache is None:
            cache = self._scan_cache = {}
        if label not in cache:
            aa = self.amino_acid.lower()
            cache[label] = parse_pes_scan_log(f"{aa}_{label}_scan.log", self.console)
        return cache[label]

    def _representative_esp_log(self, label):
        """The single ESP log representing a conformer (for the AC file / topology).

        Unscanned: its one ESP log. Scanned: the lowest-energy scan point's ESP
        log (the resolved pose is itself one scan point, so nothing is lost).
        """
        aa = self.amino_acid.lower()
        if not self._is_scanned(label):
            return f"{aa}_{label}_esp.log"
        scan = self._parse_scan(label)
        energies = scan.get("energies") or []
        k = min(range(len(energies)), key=lambda i: energies[i]) if energies else 0
        return f"{aa}_{label}_p{k}_esp.log"

    # ── Conformer selection (unified step 5) ───────────────────────
    # A "candidate" is one structure that could feed the joint RESP fit. Its
    # key is (label, point): point is None for an unscanned conformer's optimized
    # geometry, or the 0-based scan-point index for a scanned conformer. This one
    # abstraction covers crystallographic altloc conformers and relaxed-scan
    # points uniformly, so multi-conformer RESP works over whatever subset the
    # user selects at step 5.

    def _structure_candidates(self):
        """Every candidate structure that could feed the RESP fit, ordered.

        Returns a list of dicts with keys: ``key`` ``(label, point)``, ``label``,
        ``point`` (None or scan-point index), ``rel_kcal`` (relative energy within
        that conformer's scan, or None), ``esp_gjf``/``esp_log`` (the ESP job
        names), and ``atoms`` (the source geometry as ``[(element,x,y,z), ...]``
        or None if its log can't be read yet).
        """
        aa = self.amino_acid.lower()
        cands = []
        for label in self.conformers:
            if self._is_scanned(label):
                scan = self._parse_scan(label)
                geoms = scan.get("geometries") or []
                els = scan.get("elements") or []
                energies = scan.get("energies") or []
                emin = min(energies) if energies else 0.0
                for k, g in enumerate(geoms):
                    atoms = [(els[j], g[j][1], g[j][2], g[j][3]) for j in range(len(g))]
                    rel = (energies[k] - emin) * 627.509 if k < len(energies) else None
                    cands.append({
                        "key": (label, k), "label": label, "point": k,
                        "rel_kcal": rel,
                        "esp_gjf": f"{aa}_{label}_p{k}_esp.gjf",
                        "esp_log": f"{aa}_{label}_p{k}_esp.log",
                        "atoms": atoms})
            else:
                opt_log, opt_gjf = f"{aa}_{label}_opt.log", f"{aa}_{label}_opt.gjf"
                atoms = (_extract_optimized_geometry(opt_log, opt_gjf)
                         if os.path.exists(opt_log) else None)
                cands.append({
                    "key": (label, None), "label": label, "point": None,
                    "rel_kcal": None,
                    "esp_gjf": f"{aa}_{label}_esp.gjf",
                    "esp_log": f"{aa}_{label}_esp.log",
                    "atoms": atoms})
        return cands

    def _default_selection_keys(self):
        """The default RESP selection: every unscanned conformer, plus the
        lowest-energy point of each scanned conformer.

        Deterministic so a fresh-session resume (which has no in-memory
        selection) reproduces exactly what step 5 defaulted to, without
        persisting any state.
        """
        keys = []
        by_label = {}
        for c in self._structure_candidates():
            by_label.setdefault(c["label"], []).append(c)
        for label in self.conformers:
            cs = by_label.get(label, [])
            if not cs:
                continue
            if len(cs) == 1 and cs[0]["point"] is None:
                keys.append(cs[0]["key"])
            else:
                lo = min(cs, key=lambda c: (c["rel_kcal"] is None, c["rel_kcal"] or 0.0))
                keys.append(lo["key"])
        return keys

    def _selected_candidates(self):
        """Candidates chosen to feed RESP, honoring an explicit step-5 selection.

        Resolution order: an in-memory selection (``self._selected_esp_keys``),
        then a selection recorded in the step-5 result, then the deterministic
        default. Falls back to all candidates if the selection matched nothing.
        """
        keys = getattr(self, "_selected_esp_keys", None)
        if keys is None:
            recorded = self.step_results.get("step_5", {}).get("selected_keys")
            keys = [tuple(k) for k in recorded] if recorded else None
        if keys is None:
            keys = self._default_selection_keys()
        keyset = {tuple(k) for k in keys}
        cands = self._structure_candidates()
        sel = [c for c in cands if c["key"] in keyset]
        return sel or cands

    def _representative_selected_esp_log(self):
        """One ESP log to source the AC file's topology, among SELECTED structures.

        Prefers the reference conformer's lowest-energy selected structure, else
        the first selected structure. Topology (atoms/types/connectivity) is
        identical across every geometry, so any works; this just keeps the choice
        deterministic and tied to the reference conformer when possible.
        """
        sel = self._selected_candidates()
        ref = self.conformers[0] if self.conformers else None
        pool = [c for c in sel if c["label"] == ref] or sel
        pool = sorted(pool, key=lambda c: (c["rel_kcal"] is None, c["rel_kcal"] or 0.0))
        return pool[0]["esp_log"] if pool else self._representative_esp_log(ref)

    def _resolve_from_structure_esp_on_disk(self):
        """The step-6 combined-ESP path, resolved from disk (Route B, for resume).

        Step 6 writes ``{aa}_combined.esp`` when more than one structure feeds the
        joint RESP fit, else leaves the single structure's own ``.esp``. On a
        resumed run its in-memory result is gone, so reconstruct the path. Prefer
        the existing combined file: it already encodes the ACTUAL step-5
        selection, which the disk-derived default may not reproduce. Returns the
        path or None.
        """
        aa = self.amino_acid.lower()
        combined = f"{aa}_combined.esp"
        if os.path.exists(combined):
            return combined
        # No combined file: a single selected structure's own espgen output
        # (esp_log with .log -> .esp).
        esps = []
        for c in self._selected_candidates():
            log = c.get("esp_log", "")
            esp = log[:-4] + ".esp" if log.endswith(".log") else ""
            if esp and os.path.exists(esp):
                esps.append(esp)
        return esps[0] if len(esps) == 1 else None

    def _from_structure_net_charge_on_disk(self):
        """Net charge for residuegen, read from the opt job's charge/mult line.

        On a full resume the charge chosen at step 2 is no longer in memory, and
        defaulting to 0 would give a charged adduct the wrong RESP total. The
        step-2 opt input records it deterministically as the first bare
        ``<charge> <mult>`` line (after the %-lines, route and title). Returns the
        int charge, or None if it can't be read.
        """
        ref = self.conformers[0] if self.conformers else "xtal"
        gjf = f"{self.amino_acid.lower()}_{ref}_opt.gjf"
        if not os.path.exists(gjf):
            return None
        try:
            with open(gjf) as f:
                for line in f:
                    parts = line.split()
                    if (len(parts) == 2 and parts[1].isdigit()
                            and (parts[0].lstrip("-").isdigit())):
                        return int(parts[0])
        except OSError:
            return None
        return None

    def _from_structure_prep_on_disk(self):
        """The residuegen prep path, resolved from disk (Route B, for resume).

        residuegen writes ``{aa}.prep``; recover it when step 8's in-memory
        result is gone. If the user renamed it at step 8, fall back to the sole
        ``.prep`` in the directory. Returns the path or None.
        """
        cand = f"{self.amino_acid.lower()}.prep"
        if os.path.exists(cand):
            return cand
        preps = [p for p in glob.glob("*.prep")]
        return preps[0] if len(preps) == 1 else None

    def _from_structure_frcmod_on_disk(self):
        """The combined frcmod path, resolved from disk (Route B, for resume).

        step 9 writes the final combined frcmod as ``{aa}.frcmod`` (distinct from
        the ``_temp``/``_gaff`` intermediates), recovered here when step 9's
        in-memory result is gone. Returns the path or None — never a glob, so an
        intermediate frcmod can't be mistaken for the final one.
        """
        cand = f"{self.amino_acid.lower()}.frcmod"
        return cand if os.path.exists(cand) else None

    def _from_structure_ac_on_disk(self):
        """The step-7 antechamber AC path, resolved from disk (Route B, resume).

        step 7 writes ``{AA}.ac`` (uppercased 3-char residue symbol). On a resumed
        run step 7's in-memory result is gone, so anything that reaches back to it
        — residuegen (step 8) and the Seminario connectivity MOL2 (inside step 9)
        — must recover the path by its canonical name. Returns the path or None.
        """
        cand = f"{self.amino_acid.upper()[:3]}.ac"
        return cand if os.path.exists(cand) else None

    def _run_step_1_from_structure(self, **kwargs):
        """Ingest the capped model(s) and curate hydrogens (Route B, step 1).

        One capped model per crystallographic alternate conformation was built
        upstream. Hydrogens are curated interactively on the reference conformer
        (the first) through the shared HydrogenEditor; the reference's final
        atom roster is then imposed on every other conformer by running reduce
        at that conformer's own geometry and reconciling names to the reference
        (see :func:`_reconcile_conformer_to_roster`). This keeps a single
        interactive curation while guaranteeing the identical atom set/order that
        multi-conformer RESP requires. The single-conformer case is unchanged.
        """
        from proprep.forcefield_prep.hydrogen_editor import (
            HydrogenEditor, check_reduce_availability,
            configure_reduce_options_aligned, run_reduce_aligned,
        )
        interactive = kwargs.get("interactive", True)
        labels = list(self.conformers)
        ref_label = labels[0]
        ref_pdb = self._conformer_pdb_map.get(ref_label)
        if not ref_pdb or not os.path.exists(ref_pdb):
            return {"success": False,
                    "message": f"from_structure route needs a capped PDB; missing: {ref_pdb}"}

        multi = len(labels) > 1
        if multi:
            self.console.print(
                f"[green]✓ Using {len(labels)} capped model conformers "
                f"({', '.join(labels)}) from structure.[/green]\n"
                f"[grey50]Hydrogens are curated interactively on reference conformer "
                f"{ref_label}; the other conformer(s) reuse those choices — reduce runs at "
                f"each conformer's own geometry, then names are reconciled to the reference "
                f"roster so every conformer carries the identical atom set for RESP.[/grey50]")
        else:
            self.console.print(
                f"[green]✓ Using capped model compound from structure:[/green] {ref_pdb}")

        # Reference conformer: full interactive curation.
        base_ref = os.path.splitext(os.path.basename(ref_pdb))[0]
        editor = HydrogenEditor(
            ref_pdb, base_ref, console=self.console, processor=self.processor,
            interactive=interactive, residue_name=self.amino_acid,
            module="Modified Amino Acid Parameterizer",
        )
        h_result = editor.run()
        ref_curated = os.path.abspath(h_result["pdb_file"])
        self.console.print(f"[grey50]{h_result['summary']}[/grey50]")
        self._conformer_pdb_map[ref_label] = ref_curated

        if not multi:
            self.starting_pdb = ref_curated
            return {"success": True,
                    "message": "Capped structure ingested and hydrogens curated",
                    "files": [ref_curated]}

        # Secondary conformers: reduce (non-interactively, same options for all)
        # then reconcile to the reference roster. reduce's H positions here are
        # only the starting geometry for that conformer's QM optimization, which
        # relaxes them — so default options are fine; what matters is the roster.
        parsed_ref = _parse_capped_model_atoms(ref_curated)
        ref_atoms, ref_conect = parsed_ref["atoms"], parsed_ref["conect"]
        if not ref_atoms:
            return {"success": False,
                    "message": f"Could not read the curated reference conformer: {ref_curated}"}

        reduce_ok, reduce_info = check_reduce_availability()
        reduce_opts = configure_reduce_options_aligned(
            False, self.console, self.processor,
            module="Modified Amino Acid Parameterizer")

        curated_files = [ref_curated]
        for label in labels[1:]:
            raw = self._conformer_pdb_map.get(label)
            if not raw or not os.path.exists(raw):
                return {"success": False,
                        "message": f"Missing capped PDB for conformer {label}: {raw}"}
            base_l = os.path.splitext(os.path.basename(raw))[0]
            out_pdb = os.path.abspath(f"{base_l}_curated.pdb")
            reduced_pdb = raw
            if reduce_ok:
                h_added = os.path.abspath(f"{base_l}_H.pdb")
                success, message = run_reduce_aligned(raw, h_added, reduce_opts, self.console)
                if success:
                    reduced_pdb = h_added
                else:
                    self.console.print(
                        f"[yellow]⚠ reduce failed on conformer {label} ({message}); "
                        f"reconciling heavy atoms only.[/yellow]")
            else:
                self.console.print(
                    f"[yellow]⚠ reduce not available ({reduce_info}); conformer {label} "
                    f"will carry only the reference hydrogens' names with heavy-atom "
                    f"coordinates.[/yellow]")
            _reconcile_conformer_to_roster(
                ref_atoms, ref_conect, reduced_pdb, out_pdb, console=self.console)
            self._conformer_pdb_map[label] = out_pdb
            curated_files.append(out_pdb)
            self.console.print(
                f"[green]✓ Conformer {label}: hydrogens placed and reconciled to the "
                f"reference roster[/green] [grey50]({out_pdb})[/grey50]")

        self.starting_pdb = ref_curated
        return {"success": True,
                "message": (f"{len(labels)} conformers ingested; hydrogens curated on "
                            f"{ref_label} and reconciled across {', '.join(labels[1:])}"),
                "files": curated_files}

    def _run_step_1(self, **kwargs):
        """Generate initial structure with tleap (or ingest a capped model)."""
        if self.conformer_mode == "from_structure":
            return self._run_step_1_from_structure(**kwargs)
        result = run_tleap_for_amino_acid(self.amino_acid)
        if result["success"]:
            self.console.print(f"[green]✓ Generated files:[/green]")
            for f in result.get("files", []):
                self.console.print(f"  [green]•[/green] {f}")
            return {"success": True, "message": "tleap completed", "files": result.get("files", [])}
        return {"success": False, "message": result.get("error", "tleap failed")}

    def _run_step_2_from_structure(self, **kwargs):
        """Constrain the fold and write the opt+freq input (Route B, step 2).

        Chooses how the crystal fold is held during optimization — by default
        φ/ψ torsion restraints, which leave every bond and angle free to relax —
        then writes one B3LYP/6-31+G(d) ``Opt Freq`` input per conformer.

        Freq is not optional here. Its checkpoint (``IOp(7/33=1)``) holds the
        Hessian that the OPTIONAL Seminario refinement reads at Bonded
        Parameters, so it is produced once here instead of as a second job
        later. A φ/ψ-restrained optimum is a legitimate Hessian source precisely
        because only a TORSION is constrained: every bond and angle still
        converges to its own equilibrium — the condition Seminario's projection
        needs — and the single coordinate left carrying a residual gradient is a
        dihedral, a term Seminario never fits. Frozen-ATOM mode is the opposite
        (it holds bonds/angles off equilibrium), which is why Seminario is
        disabled for it. The frequencies themselves are NOT a clean normal-mode
        analysis under restraint; see step 4's validation.

        The optional relaxed scan is step 3, which starts from the geometry
        THIS step's job produces.
        """
        interactive = kwargs.get("interactive", True)
        aa = self.amino_acid.lower()
        capped = self.starting_pdb
        if not capped or not os.path.exists(capped):
            return {"success": False,
                    "message": f"Capped structure not found: {capped}"}
        model = _parse_capped_model_atoms(capped)
        phipsi = phi_psi_restraint_keys(capped)

        if not interactive:
            self._geom_constraint = {"mode": "restrain", "dihedrals": list(phipsi),
                                     "freeze_keys": set()}
            self.console.print("[green]✓ Geometry constraint: φ/ψ torsion restraints (default)[/green]")
            return self._write_opt_inputs(interactive=False)

        self.console.print(Panel(
            "[bold]Geometry Constraint — how is the crystal fold held during optimization?[/bold]\n\n"
            "  [bold]1. Restrain dihedrals[/bold] (default)\n"
            "     Hold chosen torsions (φ/ψ, and optionally the covalent linkage or a\n"
            "     custom dihedral) at their crystal values via Opt=ModRedundant. Only\n"
            "     the torsion is constrained — every bond and angle still relaxes to its\n"
            "     own equilibrium, so the Hessian remains a valid Seminario source.\n\n"
            "  [bold]2. Freeze atoms[/bold]\n"
            "     Hold chosen atoms at their exact crystal coordinates (Cartesian). Any\n"
            "     bond/angle touching a frozen atom is held off equilibrium, so Seminario\n"
            "     is not applied for this mode (planned for a later update).\n\n"
            "  [bold]3. No constraint[/bold]\n"
            "     Full relaxation; the backbone fold may drift from the crystal geometry.",
            title="Geometry Constraint", border_style="blue", expand=False))

        choice = prompt_with_context(
            self.processor, "Constraint mode", choices=["1", "2", "3"], default="1",
            module="Modified Amino Acid Parameterizer", description="Geometry constraint mode")

        if choice == "3":
            self._geom_constraint = {"mode": "none", "dihedrals": [], "freeze_keys": set()}
        elif choice == "2":
            self._geom_constraint = {"mode": "freeze", "dihedrals": [],
                                     "freeze_keys": self._choose_freeze_atoms(model)}
            self.console.print(
                "[grey50]Seminario bond/angle refinement is disabled for frozen-atom mode "
                "(bonds/angles to frozen atoms are held off equilibrium).[/grey50]")
        else:  # restrain
            dihedrals = list(phipsi)
            self.console.print(
                f"[green]φ and ψ restrained by default ({len(dihedrals)} dihedral(s)).[/green]")
            dihedrals += self._choose_linkage_dihedrals(model)
            dihedrals += self._choose_custom_dihedrals(model)
            self._geom_constraint = {"mode": "restrain", "dihedrals": dihedrals, "freeze_keys": set()}

        c = self._geom_constraint
        self.console.print(
            f"[green]✓ Constraint:[/green] "
            f"{_describe_constraint(c['mode'], c['dihedrals'], c['freeze_keys'])}")
        return self._write_opt_inputs(interactive=True)

    def _maybe_configure_scan(self, model):
        """Offer a relaxed dihedral scan and populate ``self._scan_spec``.

        ``self._scan_spec = {dihedral, conformers, n_steps, step_size}`` where
        ``dihedral`` is a 4-tuple of ``(resseq, name)`` keys. The step-2
        constraint is left untouched: ``capped_pdb_to_gaussian`` already drops
        the scanned torsion from the ``F`` restraint set when it writes the scan
        (you cannot freeze and drive the same torsion), and step 2's opt has
        already been written against that constraint.
        """
        if not confirm_with_context(
                self.processor,
                "Scan a dihedral to sample its rotational profile "
                "(de-biases RESP charges; enables an optional torsion refit)?",
                default=False, module="Modified Amino Acid Parameterizer",
                description="Add a relaxed dihedral scan"):
            return
        quad = self._choose_scan_dihedral(model)
        if not quad:
            self.console.print("[grey50]No scan dihedral chosen; skipping the scan.[/grey50]")
            return

        labels = list(self.conformers)
        scan_conformers = [labels[0]]
        if len(labels) > 1:
            ans = prompt_with_context(
                self.processor,
                f"Scan which conformer(s)? ('all' or e.g. '{labels[0]}')",
                default=labels[0], module="Modified Amino Acid Parameterizer",
                description="Conformers to scan").strip().lower()
            if ans == "all":
                scan_conformers = labels
            else:
                picked = [l for l in labels
                          if l.lower() in ans.replace(",", " ").split()]
                scan_conformers = picked or [labels[0]]

        n_steps = int_prompt_with_context(
            self.processor, "Number of scan steps", default=24,
            module="Modified Amino Acid Parameterizer", description="Scan steps")
        step_str = prompt_with_context(
            self.processor, "Step size (degrees)", default="15.0",
            module="Modified Amino Acid Parameterizer", description="Scan step size")
        try:
            step_size = float(step_str)
        except ValueError:
            step_size = 15.0

        self._scan_spec = {"dihedral": quad, "conformers": scan_conformers,
                           "n_steps": n_steps, "step_size": step_size}
        self.console.print(
            f"[grey50]Note: this ADDS roughly {n_steps} relaxed optimizations (the scan) "
            f"plus up to {n_steps} ESP single points per scanned conformer, on top of "
            f"step 2's opt+freq. Skip this step to fit RESP on the resolved pose "
            f"alone.[/grey50]")

    def _choose_scan_dihedral(self, model):
        """Pick a single dihedral to scan (linkage candidate or custom labels).

        Returns a 4-tuple of ``(resseq, name)`` keys, or None.
        """
        atoms = model.get("atoms", [])
        if not atoms:
            return None
        cands = _linkage_dihedral_candidates(model)
        if cands:
            self.console.print("\n[cyan]Covalent-linkage dihedrals (across the covalent bond):[/cyan]")
            for k, c in enumerate(cands, 1):
                self.console.print(f"  {k}. {c['label']}")
            self.console.print("  c. Custom (enter four atom labels)")
            ans = prompt_with_context(
                self.processor, "Scan which dihedral? (number, 'c' for custom, or 'none')",
                default="1" if cands else "c", module="Modified Amino Acid Parameterizer",
                description="Dihedral to scan").strip().lower()
            if ans in ("none", ""):
                return None
            if ans.isdigit() and 1 <= int(ans) <= len(cands):
                return cands[int(ans) - 1]["keys"]
        # custom entry
        labels, lookup = _atom_label_map(atoms)
        self.console.print(f"[grey50]Atoms: {', '.join(labels)}[/grey50]")
        raw = prompt_with_context(
            self.processor, "Four atom labels for the dihedral to scan (space-separated), or 'none'",
            default="none", module="Modified Amino Acid Parameterizer",
            description="Custom scan dihedral atoms").strip()
        if raw.lower() in ("none", ""):
            return None
        toks = raw.split()
        if len(toks) != 4:
            self.console.print("[red]Enter exactly four atom labels.[/red]")
            return None
        keys = [lookup.get(t) for t in toks]
        if not all(keys):
            missing = [t for t, k in zip(toks, keys) if not k]
            self.console.print(f"[red]Unknown atom(s): {', '.join(missing)}[/red]")
            return None
        return keys

    def _choose_linkage_dihedrals(self, model):
        """Offer covalent-linkage dihedrals (derived from CONECT) to restrain."""
        cands = _linkage_dihedral_candidates(model)
        if not cands:
            return []
        self.console.print("\n[cyan]Covalent-linkage dihedrals (across the covalent bond):[/cyan]")
        for k, c in enumerate(cands, 1):
            self.console.print(f"  {k}. {c['label']}")
        ans = prompt_with_context(
            self.processor, "Restrain which linkage dihedral(s)? ('1,2', 'all', or 'none')",
            default="none", module="Modified Amino Acid Parameterizer",
            description="Linkage dihedrals to restrain").strip().lower()
        if ans in ("none", ""):
            return []
        if ans == "all":
            return [c["keys"] for c in cands]
        picked = []
        for tok in ans.replace(",", " ").split():
            if tok.isdigit() and 1 <= int(tok) <= len(cands):
                picked.append(cands[int(tok) - 1]["keys"])
        return picked

    def _choose_custom_dihedrals(self, model):
        """Let the user define extra dihedral restraints by four atom labels."""
        atoms = model["atoms"]
        if not atoms:
            return []
        if not confirm_with_context(
                self.processor, "Add a custom dihedral restraint?", default=False,
                module="Modified Amino Acid Parameterizer", description="Add custom dihedral"):
            return []
        labels, lookup = _atom_label_map(atoms)
        self.console.print(f"[grey50]Atoms: {', '.join(labels)}[/grey50]")
        picked = []
        while True:
            raw = prompt_with_context(
                self.processor, "Four atom labels for the dihedral (space-separated), or 'done'",
                default="done", module="Modified Amino Acid Parameterizer",
                description="Custom dihedral atoms").strip()
            if raw.lower() == "done":
                break
            toks = raw.split()
            if len(toks) != 4:
                self.console.print("[red]Enter exactly four atom labels.[/red]")
                continue
            keys = [lookup.get(t) for t in toks]
            if not all(keys):
                missing = [t for t, k in zip(toks, keys) if not k]
                self.console.print(f"[red]Unknown atom(s): {', '.join(missing)}[/red]")
                continue
            picked.append(keys)
            self.console.print(f"[green]✓ Added dihedral {'–'.join(toks)}[/green]")
        return picked

    def _choose_freeze_atoms(self, model):
        """Pick the Cartesian-frozen atom set; returns a set of (resseq, name)."""
        atoms = model["atoms"]
        self.console.print(
            "\n[cyan]Freeze which atoms?[/cyan]\n"
            "  1. Backbone + ACE/NME caps (default)\n"
            "  2. ACE/NME caps only\n"
            "  3. Custom atom list")
        sub = prompt_with_context(
            self.processor, "Freeze set", choices=["1", "2", "3"], default="1",
            module="Modified Amino Acid Parameterizer", description="Freeze atom set")
        if sub == "2":
            return {(a["resseq"], a["name"]) for a in atoms if a["resname"] in ("ACE", "NME")}
        if sub == "3":
            labels, lookup = _atom_label_map(atoms)
            self.console.print(f"[grey50]Atoms: {', '.join(labels)}[/grey50]")
            raw = prompt_with_context(
                self.processor, "Atom labels to freeze (comma/space-separated)", default="",
                module="Modified Amino Acid Parameterizer", description="Custom freeze atoms")
            keys = set()
            for t in raw.replace(",", " ").split():
                k = lookup.get(t)
                if k:
                    keys.add(k)
                else:
                    self.console.print(f"[yellow]Unknown atom: {t}[/yellow]")
            return keys
        # default: backbone + caps
        return {(a["resseq"], a["name"]) for a in atoms
                if a["resname"] in ("ACE", "NME") or a["name"] in _MODAA_BACKBONE_FREEZE}

    def _run_step_2(self, **kwargs):
        """Set backbone angles with cpptraj (skipped for the from-structure route)."""
        if self.conformer_mode == "from_structure":
            return self._run_step_2_from_structure(**kwargs)
        aa = self.amino_acid.lower()
        parm_file = f"ace-{aa}-nme.parm7"
        rst_file = f"ace-{aa}-nme.rst7"
        result = run_cpptraj_set_angles(self.amino_acid, parm_file, rst_file)
        if result["success"]:
            self.console.print(f"[green]✓ Generated files:[/green]")
            for f in result.get("files", []):
                self.console.print(f"  [green]•[/green] {f}")
            return {"success": True, "message": "cpptraj completed", "files": result.get("files", [])}
        return {"success": False, "message": result.get("error", "cpptraj failed")}

    def _prompt_scrf_solvent(self, interactive):
        """Prompt for optional implicit solvation; return the SCRF keyword or ''.

        Shared by the geometry (Job 1) and ESP (Job 2) from-structure steps.
        """
        if not interactive:
            return ""
        if not confirm_with_context(
                self.processor, "Add implicit solvation (SCRF)?", default=False,
                module="Modified Amino Acid Parameterizer",
                description="Add SCRF implicit solvation"):
            return ""
        self.console.print("\n[cyan]Common solvents (dielectric ε):[/cyan]")
        self.console.print("  Water (78.4), DiMethylSulfoxide (46.8), Acetonitrile (35.7), Methanol (32.6),")
        self.console.print("  Ethanol (24.9), Dichloromethane (8.9), TetraHydroFuran (7.4), Chloroform (4.7), Toluene (2.4)")
        self.console.print("[grey50]Full list: Gaussian SCRF documentation[/grey50]")
        solvent = prompt_with_context(
            self.processor, "Solvent name (Gaussian keyword)", default="Water",
            module="Modified Amino Acid Parameterizer", description="SCRF solvent")
        return f"SCRF=(Solvent={solvent})"

    def _write_opt_inputs(self, interactive=True):
        """Write one Opt+Freq Gaussian input per conformer (Route B, step 2).

        Called once the geometry constraint is chosen. Resource/charge/route
        answers are asked ONCE and reused for every conformer — only the
        coordinates differ between them — and are carried forward to the ESP
        job and step-4 validation.
        """
        aa = self.amino_acid.lower()
        capped_pdb = self.starting_pdb

        constraint = getattr(self, "_geom_constraint", None) or {
            "mode": "restrain", "dihedrals": phi_psi_restraint_keys(capped_pdb),
            "freeze_keys": set()}
        mode = constraint.get("mode", "restrain")
        dihedrals = constraint.get("dihedrals") or []
        freeze_keys = constraint.get("freeze_keys")
        if mode == "restrain":
            freeze_keys = set()
            opt_kw = "Opt=ModRedundant" if dihedrals else "Opt"
        elif mode == "freeze":
            dihedrals = []
            opt_kw = "Opt"
        else:  # none
            dihedrals, freeze_keys, opt_kw = [], set(), "Opt"

        functional, basis = "B3LYP", "6-31+G(d)"
        memory, procs = "8GB", "4"
        charge, multiplicity = 0, 1
        scrf = ""
        custom_route = None
        title = f"From-structure modified AA {self.amino_acid} (geometry opt+freq)"

        def _opt_route(func, bas):
            # NoSymm keeps Gaussian in the input orientation (no reorientation),
            # so the optimized geometry stays in the same atom order we wrote.
            # Step 3's scan and the ESP job both read coordinates back out of the
            # log positionally, so this is load-bearing, not cosmetic.
            sp = f" {scrf}" if scrf else ""
            return (f"{func}/{bas} {opt_kw} Freq NoSymm{sp} "
                    f"Geom=PrintInputOrient Integral=(Grid=UltraFine) IOp(7/33=1)")

        if interactive:
            self.console.print(Panel(
                "[bold]Geometry optimization + frequencies[/bold]\n"
                "Optimizes the capped model under the constraint above. The frequency\n"
                "job is not optional: its checkpoint carries the Hessian that the\n"
                "OPTIONAL Seminario refinement reads at the Bonded Parameters step, so\n"
                "it is produced here rather than as an extra job later.\n\n"
                f"[bold]Constraint:[/bold] {_describe_constraint(mode, dihedrals, freeze_keys)}\n\n"
                "[bold]Default route keywords:[/bold]\n"
                "  • [bold]B3LYP/6-31+G(d)[/bold]        Hybrid DFT + diffuse basis (covers anions)\n"
                "  • [bold]Opt / Opt=ModRedundant[/bold]  Geometry optimization (+ restraint section)\n"
                "  • [bold]Freq[/bold]                    Harmonic frequencies → the Hessian\n"
                "  • [bold]IOp(7/33=1)[/bold]             Save the Cartesian Hessian (for Seminario)\n"
                "  • [bold]NoSymm[/bold]                  No reorientation (keeps the input orientation)\n"
                "  • [bold]Geom=PrintInputOrient[/bold]   Coordinates printed in input orientation\n"
                "  • [bold]Integral=(Grid=UltraFine)[/bold] High-quality integration grid\n"
                "  • [bold]%chk[/bold]                    Checkpoint file holding the Hessian",
                title="Geometry Optimization Input", border_style="blue", expand=False))

            memory = prompt_with_context(
                self.processor, "Memory allocation", default="8GB",
                module="Modified Amino Acid Parameterizer", description="Gaussian memory")
            procs = prompt_with_context(
                self.processor, "Number of processors", default="4",
                module="Modified Amino Acid Parameterizer", description="Gaussian processors")
            self.console.print(
                "[grey50]Charge/multiplicity depend on the modification's protonation & redox "
                "state (e.g. oxidised flavin = 0 / 1; a radical semiquinone = 0 / 2).[/grey50]")
            charge = int_prompt_with_context(
                self.processor, "Net charge of the capped model", default=0,
                module="Modified Amino Acid Parameterizer", description="Molecular charge")
            multiplicity = int_prompt_with_context(
                self.processor, "Spin multiplicity (2S+1)", default=1,
                module="Modified Amino Acid Parameterizer", description="Spin multiplicity")

            scrf = self._prompt_scrf_solvent(interactive)

            self.console.print(f"\n[bold]Recommended route:[/bold] [cyan]#p {_opt_route(functional, basis)}[/cyan]")
            if not confirm_with_context(
                    self.processor, "Use this route?",
                    default=True, module="Modified Amino Acid Parameterizer",
                    description="Use recommended geometry route"):
                functional = prompt_with_context(
                    self.processor, "DFT functional / method", default="B3LYP",
                    module="Modified Amino Acid Parameterizer", description="Geometry functional")
                basis = prompt_with_context(
                    self.processor, "Basis set", default="6-31+G(d)",
                    module="Modified Amino Acid Parameterizer", description="Geometry basis set")
                self.console.print(
                    "[grey50]Full route line — include all keywords. Keep Freq and "
                    "IOp(7/33=1): they are what make the Hessian available to Seminario. "
                    "Keep Opt=ModRedundant if you are restraining dihedrals.[/grey50]")
                custom_route = prompt_with_context(
                    self.processor, "Full route line", default=_opt_route(functional, basis),
                    module="Modified Amino Acid Parameterizer", description="Geometry route line")

        route = custom_route if custom_route is not None else _opt_route(functional, basis)

        # Carry forward to step 3's scan, the ESP job (same molecule) and step-4
        # validation. Asked once, reused for every conformer.
        self._from_structure_charge = charge
        self._from_structure_mult = multiplicity
        self._from_structure_constraint_mode = mode
        self._from_structure_qm = {"memory": memory, "procs": procs,
                                   "functional": functional, "basis": basis,
                                   "scrf": scrf}

        output_gjfs = []
        for label in self.conformers:
            conf_pdb = self._curated_conformer_pdb(label) or self.starting_pdb
            if not conf_pdb or not os.path.exists(conf_pdb):
                return {"success": False,
                        "message": f"Capped structure not found for conformer {label}: {conf_pdb}"}
            output_gjf = f"{aa}_{label}_opt.gjf"
            conf_params = {"memory": memory, "procs": procs, "keywords": route,
                           "title": (title if len(self.conformers) == 1
                                     else f"{title} [{label}]")}
            result = capped_pdb_to_gaussian(
                conf_pdb, output_gjf, conf_params, charge, multiplicity,
                freeze_keys=freeze_keys, dihedral_restraints=dihedrals,
                chk=f"{aa}_{label}_opt.chk")
            if not result.get("success"):
                return {"success": False,
                        "message": result.get("error",
                                              f"Gaussian input generation failed for {label}")}
            output_gjfs.append(output_gjf)
            self.console.print(
                f"[green]✓ Wrote {output_gjf}[/green] [grey50](opt+freq; "
                f"{result['frozen_atoms']} frozen, {result['dihedral_restraints']} dihedral "
                f"restraint(s) / {result['total_atoms']} atoms; charge {charge}, "
                f"mult {multiplicity})[/grey50]")

        self.console.print(f"[grey50]Route:[/grey50] [cyan]#p {route}[/cyan]")
        msg = ("Geometry opt+freq input generated" if len(output_gjfs) == 1
               else f"Geometry opt+freq inputs generated for {len(output_gjfs)} conformers")
        return {"success": True, "message": msg,
                "files": output_gjfs, "charge": charge, "has_hessian": True}

    def _run_step_3_from_structure(self, **kwargs):
        """Optional relaxed dihedral scan from step 2's optimized geometry.

        A scan samples a torsion's rotational profile so the RESP fit is not
        biased to the single resolved pose, and it supplies the profile an
        optional paramfit torsion refit needs. It is an ADDITION to step 2's
        opt+freq, not a replacement: the Hessian still comes from step 2.

        The scan starts from step 2's OPTIMIZED geometry, so this step first
        gates on the opt log(s) — it pauses until they exist. Coordinates are
        read back positionally (the opt route carries NoSymm, so the log's atom
        order is the input's) and re-attached to the capped model's atom names,
        which is what lets the same name-keyed restraint machinery apply.

        Skipping this step is a first-class choice: ``_is_scanned()`` keys off
        the presence of the scan input on disk, so steps 4-6 adapt on their own.
        """
        interactive = kwargs.get("interactive", True)
        aa = self.amino_acid.lower()
        capped = self.starting_pdb
        if not capped or not os.path.exists(capped):
            return {"success": False, "message": f"Capped structure not found: {capped}"}

        constraint = getattr(self, "_geom_constraint", None) or {}
        mode = constraint.get("mode", "restrain")
        if mode == "freeze":
            # A Cartesian freeze cannot relax the scanned torsion. Not an error:
            # the step is optional, so report and continue from step 2's geometry.
            self.console.print(
                "[yellow]○[/yellow] A relaxed scan is not available in frozen-atom mode: "
                "the frozen Cartesian coordinates cannot relax as the torsion is driven.\n"
                "    Continuing from step 2's optimized geometry.")
            return {"success": True,
                    "message": "No scan: unavailable in frozen-atom geometry mode"}

        # ── QM gate: the scan needs step 2's optimized geometry ──
        missing = [f"{aa}_{label}_opt.log" for label in self.conformers
                   if not os.path.exists(f"{aa}_{label}_opt.log")]
        if missing:
            run_lines = "\n".join(
                f"Input:    [cyan]{aa}_{l}_opt.gjf[/cyan]  →  Expected: [cyan]{aa}_{l}_opt.log[/cyan]"
                for l in self.conformers
                if not os.path.exists(f"{aa}_{l}_opt.log"))
            self.console.print(Panel(
                "[bold]Run the geometry opt+freq job(s) from step 2, then resume.[/bold]\n\n"
                f"{run_lines}\n\n"
                "[grey50]The scan starts from the optimized geometry, so it cannot be\n"
                "written until the optimization has finished. To proceed without a scan,\n"
                "skip this step instead.[/grey50]",
                border_style="blue", expand=False))
            return {"success": True, "status": self.PAUSE_STATUS,
                    "message": f"Run the opt+freq job(s) ({len(missing)} pending), then resume.",
                    "missing_files": missing}

        model = _parse_capped_model_atoms(capped)
        self._scan_spec = None
        if interactive:
            self._maybe_configure_scan(model)
        if not self._scan_spec:
            # Declining the scan is a legitimate outcome, not a failure: no scan
            # input is written, so _is_scanned() stays False and steps 4-6 fit
            # the resolved pose alone.
            self.console.print(
                "[grey50]No scan configured — RESP will fit the resolved pose alone.[/grey50]")
            return {"success": True, "message": "No dihedral scan configured"}

        spec = self._scan_spec
        qm = getattr(self, "_from_structure_qm", None) or {}
        memory = qm.get("memory", "8GB")
        procs = qm.get("procs", "4")
        functional = qm.get("functional", "B3LYP")
        basis = qm.get("basis", "6-31+G(d)")
        scrf = qm.get("scrf", "")
        charge = getattr(self, "_from_structure_charge", 0)
        multiplicity = getattr(self, "_from_structure_mult", 1)
        dihedrals = [d for d in (constraint.get("dihedrals") or [])]

        # A relaxed scan has no single stationary point, so it carries neither
        # Freq nor the Hessian-saving IOp: those belong to step 2's job.
        sp = f" {scrf}" if scrf else ""
        route = (f"{functional}/{basis} Opt=ModRedundant NoSymm{sp} "
                 f"Geom=PrintInputOrient Integral=(Grid=UltraFine)")
        self.console.print(Panel(
            "[bold]Relaxed dihedral scan[/bold]\n"
            "Drives one torsion through its profile, relaxing everything else at each\n"
            "point. The points feed the joint RESP fit (and an optional torsion refit).\n"
            "Starting geometry: step 2's optimized structure.\n\n"
            f"[bold]Scan:[/bold] {_describe_dihedral(spec['dihedral'], model)} — "
            f"{spec['n_steps']} steps × {spec['step_size']}°\n"
            f"[bold]Conformer(s):[/bold] {', '.join(spec['conformers'])}\n\n"
            "[grey50]No Freq / IOp(7/33=1) here: a relaxed scan has no single stationary\n"
            "point, so the Hessian for any Seminario refinement comes from step 2.[/grey50]",
            title="Dihedral Scan Input", border_style="blue", expand=False))

        output_gjfs = []
        for label in spec["conformers"]:
            conf_pdb = self._curated_conformer_pdb(label) or self.starting_pdb
            opt_log, opt_gjf = f"{aa}_{label}_opt.log", f"{aa}_{label}_opt.gjf"
            opt_atoms = _extract_optimized_geometry(opt_log, opt_gjf)
            if not opt_atoms:
                return {"success": False,
                        "message": f"Could not read the optimized geometry from {opt_log}"}
            output_gjf = f"{aa}_{label}_scan.gjf"
            conf_params = {"memory": memory, "procs": procs, "keywords": route,
                           "title": f"From-structure modified AA {self.amino_acid} "
                                    f"[{label} scan {_describe_dihedral(spec['dihedral'])}]"}
            result = capped_pdb_to_gaussian(
                conf_pdb, output_gjf, conf_params, charge, multiplicity,
                freeze_keys=set(), dihedral_restraints=dihedrals,
                chk=f"{aa}_{label}_scan.chk",
                scan_dihedral=(spec["dihedral"], spec["n_steps"], spec["step_size"]),
                coords_override=[(x, y, z) for _el, x, y, z in opt_atoms])
            if not result.get("success"):
                return {"success": False,
                        "message": result.get("error", f"Scan input generation failed for {label}")}
            output_gjfs.append(output_gjf)
            self.console.print(
                f"[green]✓ Wrote {output_gjf}[/green] [grey50](scan "
                f"{spec['n_steps']}×{spec['step_size']}° from the optimized geometry; "
                f"{result['dihedral_restraints']} dihedral restraint(s) / "
                f"{result['total_atoms']} atoms; charge {charge}, mult {multiplicity})[/grey50]")

        self.console.print(f"[grey50]Route:[/grey50] [cyan]#p {route}[/cyan]")
        return {"success": True,
                "message": (f"Scan input generated for {len(output_gjfs)} conformer(s): "
                            f"{_describe_dihedral(spec['dihedral'], model)}"),
                "files": output_gjfs}

    def _run_step_3(self, **kwargs):
        """Create Gaussian input files."""
        if self.conformer_mode == "from_structure":
            return self._run_step_3_from_structure(**kwargs)

        interactive = kwargs.get("interactive", True)

        # Get mol2 files from step 2 results
        step2_result = self.step_results.get("step_2", {})
        mol2_files = step2_result.get("files", [])
        if not mol2_files:
            mol2_files = glob.glob(f"*{self.amino_acid.lower()}*.mol2")
        if not mol2_files:
            return {"success": False, "message": "No mol2 files found from step 2"}

        # Get Gaussian parameters
        if interactive:
            self.console.print(Panel("[bold]Gaussian Settings[/bold]", border_style="blue", expand=False))
            memory = prompt_with_context(
                self.processor, "Enter memory for Gaussian job", default="8GB",
                module="Modified Amino Acid Parameterizer", description="Gaussian memory")
            procs = prompt_with_context(
                self.processor, "Enter number of processors for Gaussian job", default="4",
                module="Modified Amino Acid Parameterizer", description="Number of processors")
            keywords = prompt_with_context(
                self.processor, "Enter Gaussian keywords",
                default="HF/6-31G* opt=(modredundant) freq iop(6/33=2,6/42=6) pop=mk",
                module="Modified Amino Acid Parameterizer", description="Gaussian keywords")
            title = prompt_with_context(
                self.processor, "Enter title for Gaussian job",
                default="Dihedral scan for amino acid",
                module="Modified Amino Acid Parameterizer", description="Job title")
            charge_str = prompt_with_context(
                self.processor, "Enter charge for the molecule", default="0",
                module="Modified Amino Acid Parameterizer", description="Molecular charge")
            mult_str = prompt_with_context(
                self.processor, "Enter multiplicity (singlet=1, doublet=2, etc.)", default="1",
                module="Modified Amino Acid Parameterizer", description="Spin multiplicity")
            try:
                charge = int(charge_str)
                multiplicity = int(mult_str)
            except ValueError:
                self.console.print("[red]✗ Invalid charge/multiplicity, using defaults (0, 1)[/red]")
                charge, multiplicity = 0, 1
        else:
            memory, procs = "8GB", "4"
            keywords = "HF/6-31G* opt=(modredundant) freq iop(6/33=2,6/42=6) pop=mk"
            title = "Dihedral scan for amino acid"
            charge, multiplicity = 0, 1

        # Optional auto-generated sidechain scan — the shared sampling mode, now
        # offered on the de-novo route too so it differs from the from-structure
        # route only in input source. Opt-in; declining keeps the legacy φ/ψ-frozen
        # opt scaffold and its {aa}_{conf}.gjf names exactly as before.
        scan_dihedral = None
        if interactive:
            scan_dihedral = self._maybe_configure_mol2_scan(mol2_files[0])

        gaussian_params = {"memory": memory, "procs": procs, "keywords": keywords, "title": title}
        result = mol2_to_gaussian(self.amino_acid, mol2_files, gaussian_params,
                                  charge, multiplicity, scan_dihedral=scan_dihedral)

        if not result["success"]:
            return {"success": False, "message": result.get("error", "Gaussian input generation failed")}

        self.console.print(f"[green]✓ Generated Gaussian input files:[/green]")
        for f in result.get("files", []):
            self.console.print(f"  [green]•[/green] {f}")
        if scan_dihedral:
            self.console.print(
                "[grey50]Auto-generated relaxed sidechain scan(s); the PES-scan steps (4-5) "
                "will extract conformers and feed them to the joint RESP fit.[/grey50]")
        return {"success": True, "message": "Gaussian input generated", "files": result.get("files", [])}

    def _maybe_configure_mol2_scan(self, mol2_file):
        """Offer an auto-generated sidechain torsion scan for the de-novo route.

        Returns ``(quad, n_steps, step_size)`` with 1-based mol2 atom indices, or
        None. The dihedral is chosen by atom name against the mol2's ATOM section;
        the same indices apply to every backbone conformer (identical topology).
        """
        if not confirm_with_context(
                self.processor,
                "Auto-generate a relaxed sidechain torsion scan (samples the rotamer "
                "profile for a de-biased, multi-conformer RESP fit)?",
                default=False, module="Modified Amino Acid Parameterizer",
                description="Add auto-generated sidechain scan"):
            return None
        names = _mol2_atom_names(mol2_file)  # ordered; 1-based index = position+1
        if not names:
            self.console.print("[yellow]⚠ Could not read atoms from the mol2; skipping scan.[/yellow]")
            return None
        lookup = {n: i + 1 for i, n in enumerate(names)}
        self.console.print(f"[grey50]Atoms: {', '.join(names)}[/grey50]")
        raw = prompt_with_context(
            self.processor, "Four atom names for the dihedral to scan (space-separated), or 'none'",
            default="none", module="Modified Amino Acid Parameterizer",
            description="Sidechain scan dihedral atoms").strip()
        if raw.lower() in ("none", ""):
            return None
        toks = raw.split()
        if len(toks) != 4 or not all(t in lookup for t in toks):
            bad = [t for t in toks if t not in lookup]
            self.console.print(f"[red]Need four known atom names; unknown/short: {bad or toks}[/red]")
            return None
        quad = tuple(lookup[t] for t in toks)
        n_steps = int_prompt_with_context(
            self.processor, "Number of scan steps", default=24,
            module="Modified Amino Acid Parameterizer", description="Scan steps")
        step_str = prompt_with_context(
            self.processor, "Step size (degrees)", default="15.0",
            module="Modified Amino Acid Parameterizer", description="Scan step size")
        try:
            step_size = float(step_str)
        except ValueError:
            step_size = 15.0
        return (quad, n_steps, step_size)

    # ── Step 4: QM Gate #1 ──────────────────────────────────────────

    def _run_step_4_from_structure(self, **kwargs):
        """Validate the geometry-optimization output (Route B, QM gate 1).

        ProPrep does not run Gaussian: it ingests the log the user produced. If
        the log is absent, pause. When present, parse the frequencies and check
        for imaginary modes. Interpretation is constraint-aware — under φ/ψ
        restraints or atom freezing the constrained coordinates are held off
        their free minimum, so low imaginary modes along them are expected; with
        no constraint any imaginary mode means the optimization did not reach a
        true minimum. (Future: re-optimize by perturbing along the most-negative
        mode — not implemented yet, so we surface the finding and let the user
        decide.)
        """
        interactive = kwargs.get("interactive", True)
        aa = self.amino_acid.lower()

        # Every conformer has a step-2 opt+freq log; a conformer that was also
        # scanned at step 3 has a second, scan log. Pause until every expected
        # log is present, listing all still missing so the user can run the whole
        # batch at once.
        expected = []
        for label in self.conformers:
            expected.append((label, f"{aa}_{label}_opt.gjf", f"{aa}_{label}_opt.log", "opt"))
            if self._is_scanned(label):
                expected.append((label, f"{aa}_{label}_scan.gjf",
                                 f"{aa}_{label}_scan.log", "scan"))
        missing = [e for e in expected if not os.path.exists(e[2])]
        if missing:
            lines = "\n".join(
                f"Input:    [cyan]{gjf}[/cyan]  →  Expected: [cyan]{log}[/cyan]"
                for _label, gjf, log, _kind in missing)
            self.console.print(Panel(
                f"[bold]Run the geometry job(s), then resume.[/bold]\n\n{lines}",
                border_style="blue", expand=False))
            return {"success": True, "status": self.PAUSE_STATUS,
                    "message": (f"Run Gaussian on {len(missing)} geometry input(s), "
                                f"then resume."),
                    "missing_files": [e[2] for e in missing]}

        from proprep.forcefield_prep.small_molecule_parameterizer import check_gaussian_frequencies
        mode = getattr(self, "_from_structure_constraint_mode", "restrain")

        # Validate each log by kind: an opt+freq log gets the imaginary-mode
        # check (only these can be "flagged"); a relaxed-scan log has no single
        # frequency calculation, so it is validated by parsing its points and
        # energies and reporting the profile.
        log_files = [log for _label, _gjf, log, _kind in expected]
        flagged = []
        for label, _gjf, log, kind in expected:
            if len(expected) > 1:
                self.console.print(f"\n[bold]Conformer {label}[/bold] ({log}):")
            if kind == "scan":
                scan = self._parse_scan(label)
                n_pts = len(scan.get("geometries") or [])
                if not scan.get("success") or n_pts == 0:
                    return {"success": False,
                            "message": f"Could not parse the relaxed scan log for {label}: {log}"}
                energies = scan.get("energies") or []
                span = ""
                if len(energies) >= 2:
                    span = (f", ΔE = {(max(energies) - min(energies)) * 627.509:.1f} kcal/mol "
                            f"across the profile")
                self.console.print(
                    f"[green]✓ Relaxed scan parsed: {n_pts} point(s){span}.[/green]")
                continue
            if not check_gaussian_frequencies(log, self.console):
                flagged.append((label, log))

        if not flagged:
            n_scans = sum(1 for e in expected if e[3] == "scan")
            scanned = f" + {n_scans} relaxed scan(s)" if n_scans else ""
            msg = ("Optimized geometry validated (no imaginary modes)" + scanned
                   if len(self.conformers) == 1
                   else f"All {len(self.conformers)} conformer geometries validated "
                        f"(no imaginary modes){scanned}")
            return {"success": True, "message": msg,
                    "log_files": log_files, "workflow_type": "from_structure"}

        # Imaginary frequencies present on at least one conformer — interpret
        # against the constraint.
        flagged_labels = ", ".join(label for label, _ in flagged)
        if mode in ("restrain", "freeze"):
            self.console.print(Panel(
                f"[yellow]Imaginary frequencies with an active geometry constraint "
                f"(conformer(s): {flagged_labels}).[/yellow]\n\n"
                "Under φ/ψ restraints or atom freezing the constrained coordinates are held\n"
                "off their free minimum, so small imaginary modes along those coordinates\n"
                "are expected. A large imaginary frequency, or one localised on a bond or\n"
                "angle, still indicates a real problem.\n\n"
                "[grey50]Automatic re-optimization along the most-negative mode is planned "
                "but not yet implemented.[/grey50]",
                border_style="yellow", expand=False))
            default_proceed = True
        else:
            self.console.print(Panel(
                f"[yellow]Imaginary frequencies with no geometry constraint "
                f"(conformer(s): {flagged_labels}) — the optimization did not reach a true\n"
                "minimum. Re-optimize before proceeding.[/yellow]\n\n"
                "[grey50]Automatic re-optimization along the most-negative mode is planned "
                "but not yet implemented.[/grey50]",
                border_style="yellow", expand=False))
            default_proceed = False

        if interactive:
            if not confirm_with_context(
                    self.processor, "Proceed to the ESP step anyway?", default=default_proceed,
                    module="Modified Amino Acid Parameterizer",
                    description="Proceed despite imaginary frequencies"):
                return {"success": False,
                        "message": "Stopped: imaginary frequencies. Re-optimize, then resume."}

        return {"success": True,
                "message": "Proceeding despite imaginary frequencies (constraint-aware)",
                "log_files": log_files, "workflow_type": "from_structure"}

    def _run_step_4(self, **kwargs):
        """Check for Gaussian log files or pause for external QM.

        Detects three cases:
        - PES scan logs complete -> proceed with PES workflow
        - Opt-only logs complete -> ask user if they want to skip PES
        - No logs -> pause for user to run Gaussian externally
        """
        interactive = kwargs.get("interactive", True)

        if self.conformer_mode == "from_structure":
            return self._run_step_4_from_structure(**kwargs)

        log_check = check_for_log_files(self.amino_acid)

        # Case 1: All PES logs found
        if log_check["complete"]:
            self.console.print("[green]✓ Found PES scan log files[/green]")
            return {"success": True, "message": "PES log files found",
                    "log_files": log_check["found_files"], "workflow_type": "pes"}

        # Case 2: Opt-only logs found (no PES)
        if log_check.get("opt_files_complete"):
            self.console.print("[green]✓ Found optimization log files[/green]")
            self.console.print("[yellow]○ PES logs not found[/yellow]")
            if interactive:
                skip_pes = confirm_with_context(
                    self.processor, "Skip PES and use 2 conformations (opt-only workflow)?",
                    default=False, module="Modified Amino Acid Parameterizer",
                    description="Skip PES scans")
                if skip_pes:
                    return {"success": True, "message": "Using opt-only workflow",
                            "log_files": log_check["opt_files_found"], "workflow_type": "opt_only"}
            # User declined opt-only — fall through to pause for PES

        # Case 3: No complete logs — show instructions and offer pause
        aa_lower = self.amino_acid.lower()
        self.console.print(Panel(
            f"[bold]Expected files:[/bold] {aa_lower}_ahelix_pes.log, {aa_lower}_bsheet_pes.log\n\n"
            "You can run these as-is (simple workflow) or add PES scans for more conformations.\n\n"
            f"[cyan]Simple:[/cyan] save output as {aa_lower}_ahelix.log, {aa_lower}_bsheet.log\n"
            f"[cyan]PES scan:[/cyan] save output as {aa_lower}_ahelix_pes.log, {aa_lower}_bsheet_pes.log\n\n"
            "[grey50]The workflow auto-detects which you chose based on filenames.[/grey50]",
            border_style="blue", expand=False))

        if not interactive:
            return {
                "success": True, "status": self.PAUSE_STATUS,
                "message": "Gaussian calculations required. Run externally and resume.",
                "missing_files": log_check["missing_files"],
            }

        choice = prompt_with_context(
            self.processor, "Proceed (y), Pause to run Gaussian (p), or Quit (q)?",
            choices=["y", "p", "q"], default="p",
            module="Modified Amino Acid Parameterizer", description="Continue workflow",
            options_map={"y": "Proceed", "p": "Pause", "q": "Quit"})

        if choice == "p":
            self.console.print("[cyan]Workflow paused. Run Gaussian calculations, then resume.[/cyan]")
            self.console.print("[bold]Expected log files:[/bold]")
            for f in log_check["missing_files"]:
                self.console.print(f"  [cyan]•[/cyan] {f}")
            return {
                "success": True, "status": self.PAUSE_STATUS,
                "message": "Paused for Gaussian PES calculations. Resume when log files are ready.",
                "missing_files": log_check["missing_files"],
            }
        elif choice == "q":
            return {"success": False, "message": "Terminated by user"}

        # choice == "y": user claims logs are ready, re-check
        log_check = check_for_log_files(self.amino_acid)
        if log_check["complete"]:
            return {"success": True, "message": "PES log files found",
                    "log_files": log_check["found_files"], "workflow_type": "pes"}
        if log_check.get("opt_files_complete"):
            return {"success": True, "message": "Opt log files found",
                    "log_files": log_check["opt_files_found"], "workflow_type": "opt_only"}
        return {"success": False, "message": "Log files still missing. Cannot proceed."}

    # ── Step 5: Analysis + QM Gate #2 ──────────────────────────────

    def _run_step_5_from_structure(self, **kwargs):
        """Select which structures feed the RESP fit (Route B, unified step 5).

        The QM geometry results are in hand (step 4). This step chooses which
        conformers / relaxed-scan points contribute their ESP to the joint RESP
        fit — the "extract lowest-energy structure(s)" decision. Selecting more
        than one geometry yields a multi-conformer RESP charge set (step 8's
        residuegen fits a single charge set jointly across all of them). No QM is
        launched here; the ESP inputs for the selected structures are written and
        run at step 6.

        Default selection: every unscanned conformer, plus the lowest-energy
        point of each scanned conformer. The scanned dihedral's full profile is
        still available to the torsion refit at step 9 regardless of what is
        selected here — selection governs the CHARGE fit only.
        """
        interactive = kwargs.get("interactive", True)
        cands = self._structure_candidates()
        if not cands:
            return {"success": False,
                    "message": "No candidate structures (missing geometry logs?)"}

        default_keys = self._default_selection_keys()

        # Single candidate (one conformer, unscanned): nothing to choose.
        if len(cands) == 1:
            self._selected_esp_keys = [cands[0]["key"]]
            return {"success": True,
                    "message": "Single structure selected for charge fitting",
                    "selected_keys": [list(cands[0]["key"])],
                    "workflow_type": "from_structure"}

        # Present the candidates with their relative energies, grouped by
        # conformer, and mark which are in the default (lowest-energy) selection.
        def _fmt(c, idx):
            tag = "" if c["point"] is None else f" · scan point {c['point']}"
            rel = "" if c["rel_kcal"] is None else f"  ΔE {c['rel_kcal']:+.1f} kcal/mol"
            star = " [green](default)[/green]" if c["key"] in default_keys else ""
            return f"  [cyan]{idx:>2}[/cyan]  {c['label']}{tag}{rel}{star}"

        if interactive:
            body = "\n".join(_fmt(c, i) for i, c in enumerate(cands, 1))
            self.console.print(Panel(
                "[bold]Which structures should feed the RESP charge fit?[/bold]\n\n"
                "Each selected geometry contributes one ESP; several geometries give a\n"
                "multi-conformer RESP fit (one charge set fitted jointly). The default\n"
                "keeps every distinct conformer at its lowest energy.\n\n"
                f"{body}\n\n"
                "[grey50]Charges are usually best fitted near energy minima; adding "
                "high-energy scan points can bias them.[/grey50]",
                title="Step 5 — Select Conformers", border_style="blue", expand=False))
            choice = prompt_with_context(
                self.processor,
                "Selection — \\[d]efault (lowest-energy), \\[a]ll structures, or \\[c]ustom list?",
                default="d", choices=["d", "a", "c"],
                module="Modified Amino Acid Parameterizer",
                description="Which structures feed RESP (default / all / custom)").lower()
            if choice == "a":
                selected = [c["key"] for c in cands]
            elif choice == "c":
                raw = prompt_with_context(
                    self.processor,
                    f"Structure numbers to include (comma/space separated, 1-{len(cands)})",
                    default=" ".join(str(i) for i, c in enumerate(cands, 1)
                                     if c["key"] in default_keys),
                    module="Modified Amino Acid Parameterizer",
                    description="Custom structure selection for RESP")
                idxs = []
                for tok in re.split(r"[,\s]+", raw.strip()):
                    if tok.isdigit() and 1 <= int(tok) <= len(cands):
                        idxs.append(int(tok) - 1)
                selected = [cands[i]["key"] for i in sorted(set(idxs))]
                if not selected:
                    self.console.print(
                        "[yellow]No valid numbers — using the default selection.[/yellow]")
                    selected = list(default_keys)
            else:
                selected = list(default_keys)
        else:
            selected = list(default_keys)

        self._selected_esp_keys = selected
        chosen = [c for c in cands if c["key"] in {tuple(k) for k in selected}]
        labels_shown = ", ".join(
            f"{c['label']}" + ("" if c["point"] is None else f"/p{c['point']}")
            for c in chosen)
        self.console.print(
            f"[green]✓ {len(chosen)} structure(s) selected for RESP:[/green] "
            f"[grey50]{labels_shown}[/grey50]")
        return {"success": True,
                "message": f"{len(chosen)} structure(s) selected for charge fitting",
                "selected_keys": [list(k) for k in selected],
                "workflow_type": "from_structure"}

    def _run_step_5(self, **kwargs):
        """Analyze PES/optimization results and extract structures.

        PES path: calls analyze_pes_log_files() which internally calls
        run_batch_for_extracted_structures() — this prompts the user to
        run Gaussian on extracted structures or does it in batch mode.

        Opt-only path: calls process_opt_only_workflow() which generates
        ESP input files. User may need to run Gaussian on them.

        After either path, checks if extracted structure calculations are
        complete. If not, pauses for QM Gate #2.
        """
        interactive = kwargs.get("interactive", True)

        if self.conformer_mode == "from_structure":
            return self._run_step_5_from_structure(**kwargs)

        # Get log files and workflow type from step 4, or re-detect
        step4 = self.step_results.get("step_4", {})
        log_files = step4.get("log_files", [])
        workflow_type = step4.get("workflow_type", "pes")

        if not log_files:
            log_check = check_for_log_files(self.amino_acid)
            if log_check["complete"]:
                log_files = log_check["found_files"]
                workflow_type = "pes"
            elif log_check.get("opt_files_complete"):
                log_files = log_check["opt_files_found"]
                workflow_type = "opt_only"
            else:
                return {"success": False, "message": "No log files found for analysis"}

        # Run the appropriate analysis
        if workflow_type == "pes":
            # analyze_pes_log_files extracts structures and internally calls
            # run_batch_for_extracted_structures, which either runs Gaussian
            # in batch mode or instructs the user to run manually.
            analyze_pes_log_files(self.amino_acid, log_files, self.processor)
        else:
            # process_opt_only_workflow extracts final geometries and generates
            # ESP input files for the user to run.
            opt_result = process_opt_only_workflow(self.amino_acid, log_files, self.processor)
            if not (opt_result and opt_result.get("success")):
                return {"success": False, "message": "Opt-only structure extraction failed"}

        # QM Gate #2: Check if extracted structure calculations are complete.
        # Both paths above create .gjf files in {aa}_{conf}_structures/ dirs.
        # We need the corresponding .log files to exist before proceeding.
        structure_check = check_extracted_structure_calculations(self.amino_acid.lower())
        any_structures = False
        all_complete = True

        for conf in ["ahelix", "bsheet"]:
            if structure_check[conf]["exists"] and structure_check[conf]["gjf_files"]:
                any_structures = True
                if structure_check[conf]["missing_log_files"]:
                    all_complete = False

        if not any_structures:
            return {"success": False, "message": "No extracted structure directories found after analysis"}

        if not all_complete:
            # Show what's missing
            self.console.print("\n[yellow]Extracted structures need Gaussian calculations:[/yellow]")
            for conf in ["ahelix", "bsheet"]:
                missing = structure_check[conf].get("missing_log_files", [])
                if missing:
                    self.console.print(f"  [yellow]{conf}:[/yellow] {len(missing)} calculations pending")
            return {
                "success": True, "status": self.PAUSE_STATUS,
                "message": "Extracted structures created. Run Gaussian on them, then resume.",
            }

        return {"success": True, "message": "Analysis completed, all structure calculations done",
                "workflow_type": workflow_type}

    # ── Steps 6-10: Post-QM processing (purely sequential) ─────────

    def _run_step_6(self, **kwargs):
        """Generate ESP from optimized structures.

        Runs espgen on the completed Gaussian log files in the extracted
        structure directories and concatenates them into a combined ESP file.
        """
        if self.conformer_mode == "from_structure":
            return self._run_step_6_from_structure(**kwargs)

        # Verify extracted structure calculations are complete before proceeding
        structure_check = check_extracted_structure_calculations(self.amino_acid.lower())
        for conf in ["ahelix", "bsheet"]:
            if structure_check[conf]["exists"] and structure_check[conf]["missing_log_files"]:
                return {"success": False,
                        "message": f"Extracted {conf} structure calculations incomplete. "
                                   f"Missing: {len(structure_check[conf]['missing_log_files'])} log files"}

        result = process_extracted_structures_esp(self.amino_acid, processor=self.processor)
        if result and result.get("success"):
            return {"success": True, "message": "ESP generated",
                    "combined_file": result.get("combined_file")}
        return {"success": False, "message": "ESP generation failed"}

    def _run_step_6_from_structure(self, **kwargs):
        """Write + ingest the ESP single point(s), then espgen (Route B, step 6).

        Operates on the structures SELECTED at step 5. On first entry it writes
        one HF/6-31G* ESP single-point input per selected geometry (gas phase by
        default) and pauses until every ESP log is present — the QM gate. On
        resume it runs espgen on those logs and concatenates them into one
        multi-conformer ESP so step 8's residuegen fits a single charge set
        jointly (it reads the geometry count and sets CONF_NUM). Method/basis,
        solvent, and the full route line are shown and editable; charge and
        multiplicity carry over from the geometry job.
        """
        interactive = kwargs.get("interactive", True)
        aa = self.amino_acid.lower()

        jobs = [(c["esp_gjf"], c["esp_log"], c["atoms"])
                for c in self._selected_candidates()]
        if not jobs:
            return {"success": False,
                    "message": "No ESP jobs could be enumerated (missing geometry logs?)"}

        # Any selected structure still needing an ESP run must have a source
        # geometry readable from its (opt or scan) log by now.
        for esp_gjf, esp_log, atoms in jobs:
            if not os.path.exists(esp_log) and atoms is None:
                return {"success": False,
                        "message": f"Could not read the source geometry for {esp_gjf}"}

        pending = [(g, l, a) for g, l, a in jobs if not os.path.exists(l)]
        if pending:
            try:
                self._write_esp_inputs(pending, interactive)
            except Exception as e:
                return {"success": False, "message": str(e)}
            run_lines = "\n".join(
                f"Input:    [cyan]{gjf}[/cyan]  →  Expected: [cyan]{log}[/cyan]"
                for gjf, log, _a in pending)
            self.console.print(Panel(
                f"[bold]Run the ESP single point(s), then resume.[/bold]\n\n{run_lines}",
                border_style="blue", expand=False))
            return {"success": True, "status": self.PAUSE_STATUS,
                    "message": f"Run Gaussian on {len(pending)} ESP input(s), then resume.",
                    "missing_files": [log for _g, log, _a in pending]}

        # Every selected ESP log present → espgen + (multi-conformer) concatenate.
        esp_logs = [log for _g, log, _a in jobs]
        esp_result = run_espgen_on_log_files(esp_logs)
        if not esp_result.get("success") or not esp_result.get("esp_files"):
            return {"success": False,
                    "message": esp_result.get("error", "espgen failed")}
        esp_files = esp_result["esp_files"]
        if len(esp_files) == 1:
            # Single geometry -> the one ESP file is the "combined" ESP.
            return {"success": True, "message": "ESP generated (single conformer)",
                    "combined_file": esp_files[0]}
        # Multiple geometries (conformers and/or scan points): concatenate
        # into one multi-conformer ESP so step 8's residuegen fits a single
        # charge set across all of them (it reads the count and sets CONF_NUM).
        combined = f"{aa}_combined.esp"
        concat = concatenate_esp_files(esp_files, combined)
        if not concat.get("success"):
            return {"success": False,
                    "message": concat.get("error", "ESP concatenation failed")}
        n_conf = count_conformations_in_esp_file(combined)
        self.console.print(
            f"[green]✓ Combined {len(esp_files)} ESP file(s) → {combined}[/green] "
            f"[grey50]({n_conf} geometries; RESP fits one charge set jointly across "
            f"all of them)[/grey50]")
        return {"success": True,
                "message": f"ESP generated and combined across {len(esp_files)} geometries",
                "combined_file": combined}

    def _write_esp_inputs(self, pending, interactive):
        """Write HF/6-31G* ESP single-point inputs for the pending structures.

        ``pending`` is a list of ``(esp_gjf, esp_log, atoms)`` for structures
        whose ESP log is not yet present. The route line (method/basis, solvent,
        IOps) is shown once and applied to every input; charge/multiplicity carry
        over from the geometry job. Returns the resolved route string.
        """
        charge = getattr(self, "_from_structure_charge", 0)
        multiplicity = getattr(self, "_from_structure_mult", 1)
        memory, procs = "8GB", "4"
        functional, basis = "HF", "6-31G*"
        esp_iops = "Pop=MK IOp(6/33=2,6/41=10,6/42=10) NoSymm"
        scrf = ""
        custom_route = None
        title_base = f"From-structure modified AA {self.amino_acid} (ESP single point)"

        def _esp_route(func, bas):
            sp = f" {scrf}" if scrf else ""
            return f"{func}/{bas}{sp} {esp_iops} Integral=(Grid=UltraFine)"

        if interactive:
            self.console.print(Panel(
                "[bold]ESP Single Point[/bold]\n"
                "Evaluates the electrostatic potential at the selected geometries for RESP\n"
                "charge fitting. HF/6-31G* is the charge model AMBER/ff14SB expect.\n\n"
                "[bold]Default route keywords:[/bold]\n"
                "  • [bold]HF/6-31G*[/bold]              Charge model RESP/ff14SB expect\n"
                "  • [bold]Pop=MK[/bold]                 Merz–Kollman ESP population analysis\n"
                "  • [bold]IOp(6/33=2)[/bold]            Write ESP points + potential (for RESP)\n"
                "  • [bold]IOp(6/41=10)[/bold]           10 ESP shells/layers\n"
                "  • [bold]IOp(6/42=10)[/bold]           10 points per unit area\n"
                "  • [bold]NoSymm[/bold]                 No reorientation (grid stays aligned to atoms)\n"
                "  • [bold]Integral=(Grid=UltraFine)[/bold] High-quality integration grid",
                title="Step 6 — ESP Generation", border_style="blue", expand=False))
            memory = prompt_with_context(
                self.processor, "Memory allocation", default="8GB",
                module="Modified Amino Acid Parameterizer", description="ESP memory")
            procs = prompt_with_context(
                self.processor, "Number of processors", default="4",
                module="Modified Amino Acid Parameterizer", description="ESP processors")
            scrf = self._prompt_scrf_solvent(interactive)
            self.console.print(f"\n[bold]Recommended route:[/bold] [cyan]#p {_esp_route(functional, basis)}[/cyan]")
            if not confirm_with_context(
                    self.processor, "Use this route?",
                    default=True, module="Modified Amino Acid Parameterizer",
                    description="Use recommended ESP route"):
                functional = prompt_with_context(
                    self.processor, "Method", default="HF",
                    module="Modified Amino Acid Parameterizer", description="ESP method")
                basis = prompt_with_context(
                    self.processor, "Basis set", default="6-31G*",
                    module="Modified Amino Acid Parameterizer", description="ESP basis set")
                self.console.print(
                    "[grey50]Full route line — Pop=MK and IOp(6/33=2) are required so espgen can "
                    "read the ESP for RESP fitting.[/grey50]")
                custom_route = prompt_with_context(
                    self.processor, "Full route line", default=_esp_route(functional, basis),
                    module="Modified Amino Acid Parameterizer", description="ESP route line")

        route = custom_route if custom_route is not None else _esp_route(functional, basis)

        for esp_gjf, esp_log, atoms in pending:
            title = f"{title_base} [{os.path.splitext(os.path.basename(esp_gjf))[0]}]"
            res = write_esp_single_point_gjf(
                atoms, esp_gjf, route, memory, procs, title, charge, multiplicity)
            if not res.get("success"):
                raise RuntimeError(
                    res.get("error", f"ESP input generation failed for {esp_gjf}"))
            self.console.print(
                f"[green]✓ Wrote {esp_gjf}[/green] [grey50]({res['total_atoms']} atoms; "
                f"charge {charge}, mult {multiplicity})[/grey50]")
        self.console.print(f"[grey50]Route:[/grey50] [cyan]#p {route}[/cyan]")
        return route

    def _run_step_7_from_structure(self, **kwargs):
        """Create the AC file from the ESP QM log (Route B).

        Reuses the net charge chosen at step 3 so the user isn't asked twice.
        The ESP single-point log carries the geometry. The AC file only encodes
        topology + atom types (identical across conformers), so it is built from
        one representative ESP log among the structures SELECTED at step 5; the
        charges it will later carry come from the joint multi-conformer RESP fit
        in step 8.
        """
        # A representative ESP log among the selected structures (reference
        # conformer's lowest-energy pick when available). Topology is identical
        # across every geometry, so this only fixes atom types + connectivity.
        log_file = self._representative_selected_esp_log()
        if not os.path.exists(log_file):
            return {"success": False, "message": f"ESP log not found: {log_file}"}

        resname = self.amino_acid.upper()[:3]
        net_charge = self.step_results.get("step_3", {}).get("charge", 0)
        output_ac = f"{resname}.ac"

        result = run_antechamber(log_file, output_ac, resname, net_charge)
        if not result.get("success"):
            return {"success": False,
                    "message": result.get("error", "antechamber failed")}
        return {"success": True, "message": f"AC file generated: {output_ac}",
                "ac_file": output_ac, "residue_name": resname, "charge": net_charge}

    def _run_step_7(self, **kwargs):
        """Create AC file from lowest energy structure."""
        if self.conformer_mode == "from_structure":
            return self._run_step_7_from_structure(**kwargs)
        result = generate_ac_file(self.amino_acid, processor=self.processor)
        if result and result.get("success"):
            return {"success": True, "message": "AC file generated",
                    "ac_file": result.get("ac_file"), "charge": result.get("charge", 0)}
        return {"success": False, "message": "AC file generation failed"}

    def _run_step_8(self, **kwargs):
        """Run residuegen for RESP charges."""
        step6 = self.step_results.get("step_6", {})
        step7 = self.step_results.get("step_7", {})
        esp_file = step6.get("combined_file")
        ac_file = step7.get("ac_file")
        net_charge = step7.get("charge", 0)

        # On resume, a step completed in a PRIOR session leaves no in-memory
        # result (step_results is per-process). Step 6 in particular pauses for
        # external QM, so its combined-ESP result is usually produced in an
        # earlier session and absent here — which is exactly what makes a
        # replayed run fail with "Missing ESP" even though the file is on disk.
        # Recover the artifacts by their canonical names. The from-structure
        # route names both deterministically; the de-novo AC name is user-chosen,
        # so that route keeps its step_results-only path.
        if self.conformer_mode == "from_structure":
            if not esp_file or not os.path.exists(esp_file):
                esp_file = self._resolve_from_structure_esp_on_disk() or esp_file
            if not ac_file or not os.path.exists(ac_file):
                ac_file = self._from_structure_ac_on_disk() or ac_file
            # Recover the net charge too when step 7's result is absent, so a
            # charged adduct doesn't silently fit to total charge 0 on resume.
            if "step_7" not in self.step_results:
                recovered = self._from_structure_net_charge_on_disk()
                if recovered is not None:
                    net_charge = recovered

        if not esp_file or not ac_file:
            return {"success": False, "message": "Missing ESP or AC file from previous steps"}

        result = generate_and_run_residuegen(self.amino_acid, ac_file, esp_file, net_charge, processor=self.processor)
        if result and result.get("success"):
            return {"success": True, "message": "residuegen completed",
                    "prep_file": result.get("prep_file"),
                    "residue_symbol": result.get("residue_symbol")}
        return {"success": False, "message": "residuegen failed"}

    def _run_step_9(self, **kwargs):
        """Generate bonded parameters with parmchk2."""
        if self.conformer_mode == "from_structure":
            return self._run_step_9_from_structure(**kwargs)
        step8 = self.step_results.get("step_8", {})
        residue_symbol = step8.get("residue_symbol")
        result = generate_bonded_parameters(residue_symbol, self.processor)
        if result and result.get("success"):
            frcmod = result.get("final_frcmod") or result.get("frcmod_file")
            prep = result.get("prep_file")
            return {"success": True, "message": "Bonded parameters generated",
                    "frcmod_file": frcmod, "prep_file": prep}
        return {"success": False, "message": "Parameter generation failed"}

    def _run_step_9_from_structure(self, **kwargs):
        """Bonded parameters for Route B: Amber+GAFF2 combine, optional Seminario.

        Retains Route A's procedure so the user sees exactly what the protein
        force field's atom-typing assigned and how GAFF2 patched the ATTN gaps
        for the adduct terms it has no parameters for
        (``generate_bonded_parameters`` = parmchk2 vs the chosen FF, then
        parmchk2 vs GAFF2, then ``create_combined_frcmod`` splices GAFF2 params
        into the ATTN lines). Seminario refinement of bond/angle force constants
        is then offered as an opt-in FALLBACK, derived from a QM Hessian: either
        an existing opt+freq checkpoint (an unscanned reference conformer), or a
        dedicated frequency / opt+freq job written here at the selected lowest-
        energy geometry (the unified step-9 choice). It is disabled when the
        geometry was optimized with frozen atoms (a frozen-atom stationary point
        holds bonds/angles off-equilibrium, contaminating any Seminario term that
        touches a frozen atom).
        """
        step8 = self.step_results.get("step_8", {})
        # On resume, step 8's result is gone; the residue symbol is the adduct
        # name this workflow has carried throughout (its prep is {aa}.prep, which
        # generate_bonded_parameters discovers on disk). Fall back to it so a
        # replayed run still typed-checks and finds the prep.
        residue_symbol = step8.get("residue_symbol") or self.amino_acid.upper()[:3]

        # Amber + GAFF2 combine, with the same transparency as Route A. The
        # single-residue lib it builds here is only an intermediate (step 10
        # splits + deposits it), so suppress the standalone "copy these files"
        # Next Steps panel that only applies to Route A.
        result = generate_bonded_parameters(residue_symbol, self.processor, standalone_use=False)
        if not (result and result.get("success")):
            return {"success": False, "message": "Parameter generation failed"}

        frcmod = result.get("final_frcmod") or result.get("frcmod_file")
        prep = result.get("prep_file")

        # Opt-in Seminario fallback for bond/angle force constants.
        refined = self._maybe_seminario(residue_symbol, frcmod, prep)
        if refined:
            frcmod = refined

        # Opt-in torsion refit of the scanned dihedral against the relaxed-scan
        # energies (only offered if a scan actually ran in this parameterization).
        refit = self._maybe_torsion_refit(
            residue_symbol, frcmod, prep, interactive=kwargs.get("interactive", True))
        if refit:
            frcmod = refit

        return {"success": True, "message": "Bonded parameters generated",
                "frcmod_file": frcmod, "prep_file": prep}

    def _recover_scan_indices(self, label):
        """1-based (i,j,k,l) of the scanned torsion, read from its scan input.

        Read from the ``D … S`` line in ``{aa}_{label}_scan.gjf`` so the scanned
        dihedral is recoverable even on a fresh-session resume (where the
        in-memory ``_scan_spec`` is gone). Returns a 4-tuple or None.
        """
        aa = self.amino_acid.lower()
        gjf = f"{aa}_{label}_scan.gjf"
        if not os.path.exists(gjf):
            return None
        try:
            with open(gjf) as f:
                for line in f:
                    m = re.match(r"\s*D\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+S\b", line)
                    if m:
                        return tuple(int(g) for g in m.groups())
        except Exception:
            return None
        return None

    def _maybe_torsion_refit(self, residue_symbol, frcmod_file, prep_file,
                             interactive=True):
        """Offer a paramfit refit of the scanned torsion to the scan energies.

        Only meaningful if a relaxed scan ran. Always writes the paramfit inputs
        (scan mdcrd + QM energies) so the fit can be done manually when paramfit
        is unavailable; runs the automated paramfit chain and splices the fitted
        DIHE term when paramfit + tleap are present. Returns the refit frcmod path
        if the automated fit succeeded, else None (caller keeps the current one).
        """
        if not interactive:
            return None  # opt-in feature; never auto-run in batch mode
        scanned = [l for l in self.conformers if self._is_scanned(l)]
        if not scanned:
            return None  # no scan → nothing to refit

        scan_label = scanned[0]
        idxs = self._recover_scan_indices(scan_label)
        if not idxs:
            self.console.print(
                "[yellow]⚠ A scan ran but its scanned-dihedral indices could not be "
                "recovered; skipping the torsion refit.[/yellow]")
            return None

        # Atom-type key of the scanned dihedral, via the AMBER-typed MOL2.
        mol2 = self._ac_to_mol2(residue_symbol)
        type_quad = None
        if mol2 and os.path.exists(mol2):
            try:
                from proprep.forcefield_prep.seminario_refinement import parse_mol2_connectivity
                atom_types = parse_mol2_connectivity(mol2)["atom_types"]
                if all(1 <= i <= len(atom_types) for i in idxs):
                    type_quad = tuple(atom_types[i - 1] for i in idxs)
            except Exception:
                type_quad = None

        dihe_name = "-".join(type_quad) if type_quad else None
        penalty_txt = ""
        if dihe_name:
            cand = _group_frcmod_candidates(_index_frcmod_params(frcmod_file))
            hit = _match_dihedral_penalty(cand, type_quad)
            if hit is not None:
                pen = hit.get("penalty")
                penalty_txt = (f" Its current parmchk2 term is an ATTN placeholder."
                               if hit.get("attn") else
                               f" Its current parmchk2 term carries a penalty of {pen:.1f}."
                               if pen is not None else "")
            else:
                penalty_txt = " It is well parameterized by parmchk2 (no penalty flagged)."

        self.console.print(Panel(
            f"[bold]Torsion refit (optional)[/bold]\n\n"
            f"A relaxed scan of dihedral "
            f"[cyan]{dihe_name or '-'.join(str(i) for i in idxs)}[/cyan] was run on conformer "
            f"[cyan]{scan_label}[/cyan].{penalty_txt}\n\n"
            "Fitting this single torsion to the scan energies (paramfit) replaces the\n"
            "analogy-based term with one derived from your own QM profile.",
            title="Bonded Parameters — Torsion Refit", border_style="blue", expand=False))

        if not confirm_with_context(
                self.processor, "Refit the scanned torsion to the scan energies?",
                default=bool(penalty_txt and "penalty" in penalty_txt),
                module="Modified Amino Acid Parameterizer",
                description="Refit scanned torsion"):
            return None

        # Always emit the paramfit inputs from the scan (usable for a manual fit).
        from proprep.forcefield_prep.pes_scan_refinement import (
            write_scan_mdcrd, write_scan_energy_file, merge_dihedral_into_frcmod)
        from proprep.forcefield_prep.paramfit_refinement import check_paramfit_availability
        scan = self._parse_scan(scan_label)
        if not scan.get("success") or not scan.get("geometries"):
            self.console.print("[yellow]⚠ Could not read the scan profile; skipping refit.[/yellow]")
            return None
        outdir = "torsion_refit"
        aa = self.amino_acid.lower()
        mdcrd = write_scan_mdcrd(scan, aa, outdir, self.console)
        energies = write_scan_energy_file(scan, aa, outdir, self.console)

        prmtop = self._build_scratch_prmtop(residue_symbol, prep_file, frcmod_file, outdir)
        paramfit_ok, paramfit_msg = check_paramfit_availability()
        if not (paramfit_ok and prmtop):
            reason = paramfit_msg if not paramfit_ok else "a scratch prmtop could not be built (tleap?)"
            self.console.print(Panel(
                f"[yellow]Automated fit unavailable ({reason}).[/yellow]\n\n"
                f"The scan data for a manual fit has been written:\n"
                f"  • trajectory: [cyan]{mdcrd}[/cyan]\n"
                f"  • QM energies (Hartree): [cyan]{energies}[/cyan]\n"
                f"  • dihedral to fit: [cyan]{dihe_name or '-'.join(str(i) for i in idxs)}[/cyan]\n\n"
                "Keeping the parmchk2/GAFF2 term for now.",
                border_style="yellow", expand=False))
            return None

        if not dihe_name:
            self.console.print("[yellow]⚠ Could not resolve the scanned dihedral's atom "
                               "types; skipping the automated fit.[/yellow]")
            return None
        fitted = self._run_paramfit_torsion(
            residue_symbol, prmtop, mdcrd, energies, type_quad, len(scan["geometries"]))
        if not fitted:
            return None
        merged = os.path.join(outdir, f"{aa}_torsionfit.frcmod")
        try:
            merge_dihedral_into_frcmod(frcmod_file, fitted, dihe_name, merged, self.console)
        except Exception as e:
            self.console.print(f"[yellow]⚠ Could not splice the fitted torsion: {e}[/yellow]")
            return None
        self.console.print(f"[green]✓ Torsion refit merged → {merged}[/green]")
        return merged

    def _build_scratch_prmtop(self, residue_symbol, prep_file, frcmod_file, outdir):
        """Build a throwaway prmtop of the capped model for paramfit (best-effort).

        Uses tleap: loadamberprep + loadamberparams + saveamberparm on the capped
        residue. Returns the prmtop path, or None if tleap is unavailable/failed.
        """
        import shutil, subprocess
        if not shutil.which("tleap") or not (prep_file and frcmod_file):
            return None
        os.makedirs(outdir, exist_ok=True)
        prmtop = os.path.join(outdir, f"{residue_symbol.lower()}_scratch.prmtop")
        inpcrd = os.path.join(outdir, f"{residue_symbol.lower()}_scratch.inpcrd")
        script = os.path.join(outdir, "scratch_leap.in")
        try:
            with open(script, "w") as f:
                f.write("source leaprc.protein.ff14SB\nsource leaprc.gaff2\n")
                f.write(f"loadamberprep {prep_file}\n")
                f.write(f"loadamberparams {frcmod_file}\n")
                f.write(f"m = loadpdb {self.starting_pdb}\n")
                f.write(f"saveamberparm m {prmtop} {inpcrd}\nquit\n")
            subprocess.run(["tleap", "-f", script], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except Exception:
            return None
        return prmtop if os.path.exists(prmtop) else None

    def _run_paramfit_torsion(self, residue_symbol, prmtop, mdcrd, energies,
                              type_quad, n_structures):
        """Run the paramfit K-fit → set-params → fit chain for one DIHE term.

        Returns a fitted-parameter frcmod path, or None on failure.
        """
        from proprep.forcefield_prep.paramfit_refinement import (
            run_paramfit_k_fitting, run_paramfit_set_params_automated,
            run_paramfit_parameter_fitting)
        try:
            k = run_paramfit_k_fitting(prmtop, mdcrd, energies, n_structures, self.console)
            if k is None:
                return None
            outdir = os.path.dirname(mdcrd) or "."
            params_file = os.path.join(outdir, f"{residue_symbol.lower()}_paramfit.params")
            selected = [("-".join(type_quad), float("inf"), "SCAN", "DIHE")]
            if not run_paramfit_set_params_automated(
                    prmtop, params_file, selected, self.console, force_constants_only=False):
                return None
            fitted = os.path.join(outdir, f"{residue_symbol.lower()}_fitted.frcmod")
            if not run_paramfit_parameter_fitting(
                    prmtop, mdcrd, energies, params_file, n_structures, k, fitted, self.console):
                return None
            return fitted if os.path.exists(fitted) else None
        except Exception as e:
            self.console.print(f"[yellow]⚠ paramfit torsion fit failed: {e}[/yellow]")
            return None

    def _ac_to_mol2(self, residue_symbol):
        """Convert the step-7 antechamber AC file to MOL2 for Seminario.

        Seminario maps its selected parameters onto atom indices by matching the
        frcmod's atom-type names (e.g. ``CA-S``) against the MOL2's atom-type
        column, so the MOL2 MUST carry the same AMBER types as the frcmod. The
        conversion therefore forces ``-at amber``; without it antechamber
        re-perceives GAFF types (``ca``/``ss``/``nc``…) and every match fails
        ("Could not find bond/angle for …"). The AC spans the full capped model,
        so the MOL2 keeps the same 65-atom order as the Job-1 Hessian.
        Returns the MOL2 path, or None on failure.
        """
        ac_file = self.step_results.get("step_7", {}).get("ac_file")
        # On resume, step 7's in-memory result is gone (step_results is
        # per-process). Seminario runs inside step 9, so this is the common case,
        # not the edge case — recover the AC by its canonical on-disk name rather
        # than aborting with "Could not build a MOL2 for Seminario connectivity".
        if not ac_file or not os.path.exists(ac_file):
            ac_file = self._from_structure_ac_on_disk()
        if not ac_file or not os.path.exists(ac_file):
            return None
        mol2_file = f"{residue_symbol.lower()}.mol2"
        cmd = ["antechamber", "-i", ac_file, "-fi", "ac",
               "-o", mol2_file, "-fo", "mol2", "-at", "amber", "-pf", "y"]
        try:
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError, OSError) as e:
            _console.print(f"[yellow]⚠[/yellow] Could not convert {ac_file} to MOL2: {e}")
            return None
        return mol2_file if os.path.exists(mol2_file) else None

    def _maybe_seminario(self, residue_symbol, frcmod_file, prep_file):
        """Offer Seminario refinement of bond/angle force constants (Route B).

        Returns the fitted frcmod path if refinement ran and succeeded, else
        None (caller keeps the combined frcmod). Disabled for the freeze
        constraint mode.
        """
        from proprep.forcefield_prep.small_molecule_parameterizer import analyze_frcmod_penalties

        mode = getattr(self, "_from_structure_constraint_mode", "restrain")
        aa = self.amino_acid.lower()
        ref_label = self.conformers[0]

        # Frozen-atom mode is disqualifying regardless of any checkpoint, so
        # rule it out before looking for one.
        if mode == "freeze":
            _console.print(
                "[yellow]○[/yellow] Seminario is unavailable for this residue: the geometry "
                "was optimized with [bold]frozen atoms[/bold].\n"
                "    A frozen-atom stationary point holds bonds and angles off their "
                "equilibrium values, so Hessian-derived\n"
                "    force constants for any term touching a frozen atom would be "
                "contaminated. Keeping the Amber/GAFF2 parameters."
            )
            # Future extension: refine only terms with no frozen atom.
            return None

        # The Hessian ALWAYS comes from step 2's opt+freq — that is the whole
        # point of running Freq there. Seminario derives LOCAL bond/angle
        # curvatures, so the reference conformer's optimized geometry is the
        # right source; the scan-point selection at step 5 is a RESP-charge
        # concern and irrelevant here. Prefer an already-formatted checkpoint,
        # else the binary one (convert_chk_to_fchk runs formchk). A
        # {aa}_seminario.* from a prior manual run is honored for back-compat.
        candidates = [f"{aa}_{ref_label}_opt.fchk", f"{aa}_{ref_label}_opt.chk",
                      f"{aa}_seminario.fchk", f"{aa}_seminario.chk"]
        chk_file = next((c for c in candidates if os.path.exists(c)), None)

        if not chk_file:
            # Step 2's opt+freq already computed the Hessian; its checkpoint just
            # isn't in this directory (a %chk binary is easy to leave behind on a
            # cluster). This is a "make the file available" situation, NOT a
            # "run another opt+freq" one — so say exactly that.
            opt_chk = f"{aa}_{ref_label}_opt.chk"
            _console.print(Panel(
                "[bold]Seminario needs step 2's Hessian checkpoint, which isn't here.[/bold]\n\n"
                f"Step 2's opt+freq already computed it (the [cyan]%chk[/cyan] line named\n"
                f"[cyan]{opt_chk}[/cyan], via IOp(7/33=1)) — it just isn't in this directory.\n\n"
                "To use Seminario, make that checkpoint available, then re-run this step:\n"
                f"  • copy [cyan]{opt_chk}[/cyan] back from wherever the opt job ran, or\n"
                f"  • if you kept the checkpoint elsewhere, run "
                f"[cyan]formchk {opt_chk} {aa}_{ref_label}_opt.fchk[/cyan].\n\n"
                "[grey50]No extra QM is needed — this is the same opt+freq you already ran.\n"
                "Keeping the Amber/GAFF2 parameters for now.[/grey50]",
                title="Bonded Parameters — Seminario", border_style="yellow",
                expand=False,
            ))
            return None

        _console.print(Panel(
            "[bold]Seminario refinement (optional)[/bold]\n\n"
            "parmchk2 could not read every bond/angle of this adduct straight from\n"
            "the force field, so it filled the missing ones [bold]by analogy[/bold] and tagged\n"
            "each with a [cyan]penalty score[/cyan] (how far it reached for the analogy; ATTN =\n"
            "no analogy found). You then swapped some for GAFF2. Those by-analogy\n"
            "terms — [bold]every[/bold] penalty-scored / ATTN term, including the ones now at\n"
            "penalty 0.0 — are empirical guesses. Seminario replaces a guess with a\n"
            "value derived directly from the step-2 QM Hessian of THIS adduct.\n\n"
            "You'll pick the scope next: just those by-analogy terms, or EVERY\n"
            "bond/angle term in the frcmod (a fully QM-derived custom set, including\n"
            "the terms parmchk2 matched directly). Either way Seminario sets both the\n"
            "force constant (from the Hessian) AND the equilibrium length/angle (from\n"
            "the QM geometry); it does not touch dihedrals or impropers, and standard\n"
            "backbone terms remain ff14SB.\n"
            f"[cyan]Source:[/cyan] {chk_file} (step 2 opt+freq Hessian).",
            title="Bonded Parameters — Seminario",
            border_style="blue",
            expand=False,
        ))

        penalties = analyze_frcmod_penalties(frcmod_file, _console)
        flagged = [p for p in penalties if p[3] in ("BOND", "ANGLE")]
        all_ba = _all_bond_angle_params(frcmod_file)
        if not all_ba:
            _console.print(
                "[grey50]No bond/angle terms in the combined frcmod; "
                "nothing to refine with Seminario.[/grey50]"
            )
            return None

        if flagged:
            table = Table(title="Bond/angle terms parmchk2 filled by analogy "
                                "(penalty-scored / ATTN) — the by-analogy set")
            table.add_column("#", style="cyan", justify="right")
            table.add_column("Parameter", style="white")
            table.add_column("Section", style="magenta")
            table.add_column("Status", style="yellow")
            for i, (name, score, status, section) in enumerate(flagged, 1):
                score_txt = "ATTN" if score == float("inf") else f"{score:.1f}"
                table.add_row(str(i), name, section, f"{status} ({score_txt})")
            _console.print(table)
            _console.print(
                f"[grey50]These {len(flagged)} term(s) are the ones parmchk2 assigned by "
                f"analogy (a penalty score — 0.0 just means it found a close analog — or "
                f"ATTN); they are the adduct's non-standard terms.\n"
                f"Scope: [bold]by-analogy[/bold] = refine those {len(flagged)} only (hybrid — "
                f"the {len(all_ba) - len(flagged)} directly-matched term(s) stay Amber/GAFF2); "
                f"[bold]all[/bold] = every bond/angle term in the frcmod ({len(all_ba)} total), "
                f"a fully QM-derived custom set (standard backbone terms remain ff14SB "
                f"either way).[/grey50]"
            )
            default_choice = "b"
        else:
            _console.print(
                f"[grey50]parmchk2 matched every bond/angle term directly (none needed an "
                f"analogy / carried a penalty). The combined frcmod still has {len(all_ba)} "
                f"custom bond/angle term(s) that Seminario can derive from the QM "
                f"Hessian.[/grey50]"
            )
            default_choice = "n"

        choice = prompt_with_context(
            self.processor,
            "\nSeminario scope — \\[b]y-analogy terms only, \\[a]ll frcmod bond/angle terms, or \\[n]o refinement?",
            default=default_choice,
            choices=["b", "a", "n"],
            module="Modified Amino Acid Parameterizer",
            description="Seminario refinement scope (Route B): by-analogy / all / none",
        ).lower()

        if choice == "n" or (choice == "b" and not flagged):
            _console.print("[grey50]Keeping the Amber/GAFF2 parameters.[/grey50]")
            return None

        ba_params = all_ba if choice == "a" else flagged
        _console.print(
            f"[cyan]Seminario will refine {len(ba_params)} bond/angle term(s) "
            f"({'all frcmod terms' if choice == 'a' else 'by-analogy terms only'}).[/cyan]"
        )

        mol2_file = self._ac_to_mol2(residue_symbol)
        if not mol2_file:
            _console.print(
                "[yellow]⚠[/yellow] Could not build a MOL2 for Seminario connectivity; "
                "keeping the Amber/GAFF2 parameters."
            )
            return None

        try:
            from proprep.forcefield_prep.seminario_refinement import run_seminario_refinement_workflow
            res = run_seminario_refinement_workflow(
                mol_name=residue_symbol.lower(), mol2_file=mol2_file,
                frcmod_file=frcmod_file, chk_file=chk_file,
                selected_params=ba_params, console=_console, interactive=True,
            )
        except ImportError as e:
            _console.print(f"[yellow]⚠[/yellow] Seminario module unavailable: {e}")
            return None
        except Exception as e:
            _console.print(f"[yellow]⚠[/yellow] Seminario refinement error: {e}")
            return None

        if res.get("refinement_success") and res.get("fitted_frcmod"):
            _console.print(f"[green]✓[/green] Seminario refinement complete: {res['fitted_frcmod']}")
            return res["fitted_frcmod"]
        _console.print(
            f"[yellow]⚠[/yellow] Seminario did not complete: {res.get('message', 'unknown error')}. "
            "Keeping the Amber/GAFF2 parameters."
        )
        return None

    def _run_step_10(self, **kwargs):
        """Create AMBER library file (de-novo) / FF integration (from-structure)."""
        if self.conformer_mode == "from_structure":
            return self._run_step_10_from_structure(**kwargs)
        step8 = self.step_results.get("step_8", {})
        step9 = self.step_results.get("step_9", {})
        prep_file = step8.get("prep_file") or step9.get("prep_file")
        frcmod_file = step9.get("frcmod_file")
        residue_symbol = step8.get("residue_symbol", self.amino_acid)

        if not prep_file or not frcmod_file:
            return {"success": False, "message": "Missing prep or frcmod file from previous steps"}

        result = generate_amber_library(prep_file, frcmod_file, residue_symbol)
        if not (result and result.get("success")):
            return {"success": False, "message": "Library generation failed"}

        # Register the lib + frcmod with the Topology Generator, the same way
        # the from-structure route and the metal-site step do. Without this the
        # de-novo library is written to disk but never loaded for the system
        # build, so tLEaP has no unit for the modified residue.
        lib_file = result.get("lib_file")
        workspace = getattr(self.processor, "workspace", None) if self.processor else None
        registered = _register_for_topology_generator(
            workspace, lib_file=lib_file, frcmod_file=frcmod_file,
            console=self.console,
        )

        # Reuse deposit, matching the from-structure route and the metal-site
        # step. No reuse transformer here: the de-novo route builds the residue
        # from its parent rather than renaming atoms off an input structure, so
        # there is no atom-name map to bake into one.
        from proprep.forcefield_prep.library_promotion import offer_library_promotion

        promo = offer_library_promotion(
            self.console, self.processor,
            category="modified_amino_acid",
            residue_name=residue_symbol,
            frcmod_file=frcmod_file,
            lib_search_dir=os.path.dirname(lib_file) if lib_file else os.getcwd(),
            prep_file=prep_file,
            lib_file=lib_file,
        )

        bits = []
        if registered:
            bits.append("registered")
        if (promo or {}).get("library_path"):
            bits.append("deposited")
        message = ("AMBER library generated"
                   + (f" and {' and '.join(bits)}" if bits else ""))
        return {"success": True, "message": message, "lib_file": lib_file}

    def _run_step_10_from_structure(self, **kwargs):
        """Force-field integration for a Route B adduct (final checklist step).

        Splits the combined residue into an amino-acid + cofactor pair, deposits
        both + the reused frcmod into the user library, renames residues/atoms in
        the prepared PDB, populates the Topology-Generator workspace keys, and
        emits a reuse transformer — mirroring the metal-site integration step.
        """
        step8 = self.step_results.get("step_8", {})
        step9 = self.step_results.get("step_9", {})
        prep_file = step8.get("prep_file") or step9.get("prep_file")
        frcmod_file = step9.get("frcmod_file")

        # On resume, steps 8-9 may have completed in a prior session, leaving no
        # in-memory result — recover the prep (step 8) and combined frcmod
        # (step 9) from their canonical on-disk names so integration still runs.
        if not prep_file or not os.path.exists(prep_file):
            prep_file = self._from_structure_prep_on_disk() or prep_file
        if not frcmod_file or not os.path.exists(frcmod_file):
            frcmod_file = self._from_structure_frcmod_on_disk() or frcmod_file

        if not prep_file or not frcmod_file:
            return {"success": False, "message": "Missing prep or frcmod file from previous steps"}

        workspace = getattr(self.processor, "workspace", None) if self.processor else None
        result = integrate_modaa_from_structure(
            console=self.console, workspace=workspace,
            residue_name=self.amino_acid, source_residues=self.source_residues,
            capped_pdb=self.starting_pdb, prep_file=prep_file,
            frcmod_file=frcmod_file, output_dir=os.getcwd(),
            conformer_label=self.conformers[0],
        )
        return result

    # ── WorkflowChecklist handler methods ─────────────────────────

    def _checklist_aa_1_structure(self):
        """Handler: generate initial structure with tleap."""
        result = self._run_step_1(interactive=True)
        if not result.get("success"):
            raise RuntimeError(result.get("message", "tleap failed"))
        self.step_results["step_1"] = result
        return {'summary': result.get("message", "Structure generated")}

    def _checklist_aa_2_angles(self):
        """Handler: set backbone angles with cpptraj."""
        result = self._run_step_2(interactive=True)
        if not result.get("success"):
            raise RuntimeError(result.get("message", "cpptraj failed"))
        self.step_results["step_2"] = result
        return {'summary': result.get("message", "Backbone angles set")}

    def _checklist_aa_3_gaussian_input(self):
        """Handler: create Gaussian input files.

        Route B's step 3 is the OPTIONAL relaxed scan, and it reads step 2's
        optimized geometry — so it can pause for QM the way steps 4-6 do.
        """
        result = self._run_step_3(interactive=True)
        self.step_results["step_3"] = result
        if result.get("status") == self.PAUSE_STATUS:
            return {'checkpoint': True}
        if not result.get("success"):
            raise RuntimeError(result.get("message", "Gaussian input generation failed"))
        return {'summary': result.get("message", "Gaussian input files created")}

    def _checklist_aa_4_pes_scan(self):
        """Handler: check for Gaussian PES outputs (QM Gate #1)."""
        result = self._run_step_4(interactive=True)
        self.step_results["step_4"] = result
        if result.get("status") == self.PAUSE_STATUS:
            return {'checkpoint': True}
        if not result.get("success"):
            raise RuntimeError(result.get("message", "PES scan check failed"))
        return {'summary': result.get("message", "PES scan results processed")}

    def _checklist_aa_5_extract(self):
        """Handler: analyze PES and extract structures (QM Gate #2)."""
        result = self._run_step_5(interactive=True)
        self.step_results["step_5"] = result
        if result.get("status") == self.PAUSE_STATUS:
            return {'checkpoint': True}
        if not result.get("success"):
            raise RuntimeError(result.get("message", "Structure extraction failed"))
        return {'summary': result.get("message", "Structures extracted")}

    def _checklist_aa_6_esp(self):
        """Handler: write + ingest the ESP single point(s), then run espgen.

        For the from-structure route this step now writes the ESP inputs for the
        selected conformers and pauses (QM gate) until their logs are present, so
        it can return a checkpoint just like the geometry gate at step 4.
        """
        result = self._run_step_6(interactive=True)
        self.step_results["step_6"] = result
        if result.get("status") == self.PAUSE_STATUS:
            return {'checkpoint': True}
        if not result.get("success"):
            raise RuntimeError(result.get("message", "ESP generation failed"))
        return {'summary': result.get("message", "ESP generated")}

    def _checklist_aa_7_ac_file(self):
        """Handler: create AC file from lowest energy structure."""
        result = self._run_step_7(interactive=True)
        if not result.get("success"):
            raise RuntimeError(result.get("message", "AC file generation failed"))
        self.step_results["step_7"] = result
        return {'summary': result.get("message", "AC file generated")}

    def _checklist_aa_8_resp(self):
        """Handler: run residuegen for RESP charges."""
        result = self._run_step_8(interactive=True)
        if not result.get("success"):
            raise RuntimeError(result.get("message", "residuegen failed"))
        self.step_results["step_8"] = result
        return {'summary': result.get("message", "RESP charges fitted")}

    def _checklist_aa_9_parmchk2(self):
        """Handler: generate bonded parameters with parmchk2."""
        result = self._run_step_9(interactive=True)
        if not result.get("success"):
            raise RuntimeError(result.get("message", "Parameter generation failed"))
        self.step_results["step_9"] = result
        return {'summary': result.get("message", "Bonded parameters generated")}

    def _checklist_aa_10_library(self):
        """Handler: AMBER library (de-novo) or FF integration (from-structure)."""
        result = self._run_step_10(interactive=True)
        if not result.get("success"):
            raise RuntimeError(result.get("message", "Final step failed"))
        self.step_results["step_10"] = result
        return {'summary': result.get("message", "Force field integrated")}


def run_workflow(amino_acid=None, output_dir=None, interactive=True, processor=None,
                 starting_pdb=None, conformer_mode="denovo_parent",
                 source_residues=None, conformer_pdbs=None):
    """
    Multi-step workflow to generate amino acid structures with different conformations,
    run PES scans, analyze the results, generate ESP files, create residue parameters,
    and generate bonded parameters.

    Parameters:
    -----------
    amino_acid : str, optional
        Three-letter amino acid code. If None, will prompt the user.
    output_dir : str, optional
        Directory to save output files. If None, uses current directory.
    interactive : bool, optional
        Whether to run in interactive mode, prompting the user for input.
    processor : optional
        Processor object with session_manager for session recording support.

    Returns:
    --------
    dict
        Dictionary containing results of the parameterization
    """
    # Resolve before any chdir so the capped structure(s) stay findable.
    if starting_pdb:
        starting_pdb = os.path.abspath(starting_pdb)
    if conformer_pdbs:
        conformer_pdbs = [(label, os.path.abspath(path)) for label, path in conformer_pdbs]

    original_dir = None
    if output_dir:
        original_dir = os.getcwd()
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)

    try:
        if amino_acid is None:
            amino_acid = prompt_with_context(
                processor, "Enter a three-letter amino acid code (e.g., ALA, PHE, GLY)",
                module="Modified Amino Acid Parameterizer", description="Amino acid code"
            )

        amino_acid = amino_acid.strip().upper()
        work_dir = Path(".")

        manager = ModifiedAAWorkflowManager(
            amino_acid, processor=processor,
            starting_pdb=starting_pdb, conformer_mode=conformer_mode,
            source_residues=source_residues, conformer_pdbs=conformer_pdbs)

        steps = (MODIFIED_AA_FROM_STRUCTURE_STEPS
                 if conformer_mode == "from_structure" else MODIFIED_AA_STEPS)
        checklist = WorkflowChecklist(
            steps=steps,
            executor=manager,
            processor=processor,
            workflow_name="Modified Amino Acid Parameterization",
            console=manager.console,
            state_dir=work_dir,
        )
        success = checklist.run()

        # Collect parameter files from step results
        results = {
            "success": success,
            "status": "completed" if success else "incomplete",
            "message": f"Modified amino acid parameterization for {amino_acid}",
            "parameter_files": {},
            "amino_acid": amino_acid,
        }

        step8 = manager.step_results.get("step_8", {})
        step9 = manager.step_results.get("step_9", {})
        step10 = manager.step_results.get("step_10", {})

        if step8.get("prep_file"):
            results["parameter_files"]["prep_file"] = os.path.abspath(step8["prep_file"])
        if step9.get("frcmod_file"):
            results["parameter_files"]["frcmod_file"] = os.path.abspath(step9["frcmod_file"])
        if step10.get("lib_file"):
            results["parameter_files"]["lib_file"] = os.path.abspath(step10["lib_file"])

        return results

    finally:
        if original_dir:
            os.chdir(original_dir)


def parameterize_modified_residues(component_analyzer, structure, output_dir=None, processor=None):
    """
    Interface between PDB processor and the amino acid parameterization workflow.
    This function identifies candidate residues for parameterization and passes them
    to the main workflow.

    Args:
        component_analyzer: ComponentAnalyzer instance with analyzed structure
        structure: BioPython structure object
        output_dir: Directory to store parameter files

    Returns:
        dict: Results of parameterization
    """
    logger.info("Starting identification of modified amino acids")

    # Set output directory
    if output_dir is None:
        output_dir = os.getcwd()

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Result tracking
    results = {
        "success": True,
        "parameterized_residues": {},
        "tleap_commands": [],
        "errors": [],
        "pending_workflows": {},  # Track workflows that need additional steps
    }

    # Get all residue types from the component analyzer
    all_residue_types = {}
    modified_residues = {}

    # First, collect all residue types and their occurrences
    for key, details in component_analyzer.residue_details.items():
        chain_id, res_num, icode, res_name = key

        # Add to all residues dictionary
        if res_name not in all_residue_types:
            all_residue_types[res_name] = []

        all_residue_types[res_name].append(
            {
                "chain_id": chain_id,
                "res_num": res_num,
                "icode": icode,
                "category": details["category"],
            }
        )

        # Separately collect modified residues
        if details["category"] == "modified_protein":
            if res_name not in modified_residues:
                modified_residues[res_name] = []

            modified_residues[res_name].append(
                {"chain_id": chain_id, "res_num": res_num, "icode": icode}
            )

    # Check if we found any residues
    if not all_residue_types:
        logger.info("No residues found in the structure")
        return {
            "success": True,
            "message": "No residues found in the structure",
            "parameterized_residues": {},
        }

    # Report findings to the user
    _console.print(f"\n[bold cyan]Residue Types in Structure[/bold cyan]")
    _console.print(f"Found {len(all_residue_types)} unique residue types in the structure:")
    _console.print(f"[grey50](* indicates residues automatically detected as modified amino acids)[/grey50]\n")

    # Sort residue types alphabetically and display in a table
    sorted_residue_list = sorted(all_residue_types.items())
    residue_table = Table(title="Available Residues")
    residue_table.add_column("#", style="cyan", width=4)
    residue_table.add_column("", width=2)  # Modified marker
    residue_table.add_column("Residue", style="bold")
    residue_table.add_column("Category", style="grey50")
    residue_table.add_column("Count", style="green")
    residue_table.add_column("Locations", style="grey50")

    for i, (res_name, instances) in enumerate(sorted_residue_list, 1):
        # Check if this is a modified residue
        is_modified = res_name in modified_residues
        marker = "[yellow]*[/yellow]" if is_modified else ""

        # Get locations for the first few instances
        locations = [f"{inst['chain_id']}:{inst['res_num']}" for inst in instances[:3]]
        if len(instances) > 3:
            locations.append("...")

        category = instances[0]["category"]
        residue_table.add_row(
            str(i), marker, res_name, category,
            str(len(instances)), ", ".join(locations)
        )
    _console.print(residue_table)

    # Ask if user wants to parameterize any residues
    if not confirm_with_context(
        processor,
        "\nWould you like to parameterize any of these residues?",
        default=True,
        module="Modified AA Parameterizer",
        description="Parameterize any modified residues",
    ):
        logger.info("User chose not to parameterize residues")
        return {
            "success": True,
            "message": "User chose not to parameterize residues",
            "parameterized_residues": {},
        }

    # Track parameterized residues
    parameterized_residues = {}
    pending_workflows = {}

    # Parameterization loop - continue until user chooses to stop
    while True:
        # User selection
        _console.print(
            "\n[bold]Select a residue to parameterize by number, or enter the residue name:[/bold]"
        )
        _console.print("[grey50](Enter 'q' to quit parameterization)[/grey50]")

        selection = prompt_with_context(
            processor,
            "Selection",
            module="Modified AA Parameterizer",
            description="Select residue to parameterize (number, name, or 'q' to quit)",
        )

        if selection.lower() == "q":
            break

        # Try to convert to integer for list selection
        selected_res = None
        try:
            idx = int(selection)
            if 1 <= idx <= len(sorted_residue_list):
                # Get the residue name by index (using sorted order)
                selected_res = sorted_residue_list[idx - 1][0]
            else:
                _console.print(
                    f"[red]✗ Invalid selection. Please enter a number between 1 and {len(sorted_residue_list)}[/red]"
                )
                continue
        except ValueError:
            # Not a number, treat as a residue name (convert to uppercase)
            selected_res = selection.upper()

            # Validate that the residue exists in the structure
            if selected_res not in all_residue_types:
                _console.print(
                    f"[red]✗ Residue '{selected_res}' not found in the structure. Please try again.[/red]"
                )
                continue

        # Check if this residue has already been parameterized
        if selected_res in parameterized_residues:
            _console.print(
                f"\n[yellow]⚠ Residue {selected_res} has already been parameterized in this session.[/yellow]"
            )
            if not confirm_with_context(
                processor,
                "Do you want to re-parameterize it?",
                default=True,
                module="Modified AA Parameterizer",
                description="Re-parameterize residue",
            ):
                continue

        # Prepare for parameterization of the selected residue
        _console.print(f"\n[cyan]→ Preparing to parameterize {selected_res}...[/cyan]")

        # Create residue-specific output directory
        residue_output_dir = os.path.join(output_dir, selected_res.lower())

        try:
            # Call the run_workflow function with the selected residue
            workflow_result = run_workflow(
                amino_acid=selected_res,
                output_dir=residue_output_dir,
                interactive=True,  # Ensure interactive mode for proper prompting
            )

            # Process the workflow results
            if workflow_result["success"]:
                # Check if parameter files were created
                prep_file = workflow_result["parameter_files"].get("prep_file")
                frcmod_file = workflow_result["parameter_files"].get("frcmod_file")

                # Check status of workflow
                if (
                    workflow_result.get("status") == "paused"
                    or workflow_result.get("status") == "pending_calculations"
                ):
                    # Workflow needs Gaussian calculations
                    missing_files = workflow_result.get("missing_files", [])
                    _console.print(
                        f"\n[yellow]○ Parameterization of {selected_res} requires Gaussian calculations to complete.[/yellow]"
                    )
                    _console.print(f"[bold]The following files need to be created:[/bold]")
                    for file in missing_files:
                        _console.print(f"  [yellow]•[/yellow] {file}")

                    # Track this workflow as pending
                    pending_workflows[selected_res] = {
                        "output_dir": residue_output_dir,
                        "missing_files": missing_files,
                        "status": workflow_result.get("status"),
                    }

                    # Ask if user wants to continue with other residues
                    _console.print(
                        "\n[grey50]You can continue parameterizing other residues, or quit to run Gaussian calculations.[/grey50]"
                    )

                elif prep_file or frcmod_file:
                    # Workflow completed successfully with parameter files
                    parameterized_residues[selected_res] = {}

                    if prep_file:
                        parameterized_residues[selected_res]["prep_file"] = prep_file
                        results["tleap_commands"].append(f"loadAmberPrep {prep_file}")

                    if frcmod_file:
                        parameterized_residues[selected_res][
                            "frcmod_file"
                        ] = frcmod_file
                        results["tleap_commands"].append(
                            f"loadAmberParams {frcmod_file}"
                        )

                    _console.print(f"\n[green]✓ Successfully parameterized {selected_res}[/green]")
                else:
                    # If workflow completed but didn't create parameter files
                    _console.print(
                        f"\n[yellow]⚠ Workflow completed for {selected_res}, but no parameter files were generated.[/yellow]"
                    )
                    msg = workflow_result.get(
                        "message",
                        "Manual steps may be required to complete parameterization",
                    )
                    _console.print(f"[grey50]Message: {msg}[/grey50]")
                    results["errors"].append(
                        f"Incomplete parameterization for {selected_res}: {msg}"
                    )
            else:
                # If workflow failed
                error = workflow_result.get(
                    "error", workflow_result.get("message", "Unknown error")
                )
                _console.print(
                    f"\n[red]✗ Error in parameterization workflow for {selected_res}: {error}[/red]"
                )
                results["errors"].append(
                    f"Error parameterizing {selected_res}: {error}"
                )

        except Exception as e:
            logger.error(f"Error during parameterization of {selected_res}: {str(e)}")
            _console.print(f"\n[red]✗ Error during parameterization: {str(e)}[/red]")
            results["errors"].append(
                f"Error during parameterization of {selected_res}: {str(e)}"
            )

        # Ask if user wants to parameterize another residue
        if not confirm_with_context(
            processor,
            "Would you like to parameterize another residue?",
            default=True,
            module="Modified AA Parameterizer",
            description="Parameterize another residue",
        ):
            break

    # Summary of parameterization
    if parameterized_residues:
        _console.print("\n[bold cyan]Parameterization Summary[/bold cyan]")
        _console.print(
            f"[green]✓ Successfully parameterized {len(parameterized_residues)} residue types:[/green]"
        )

        for res_name, files in parameterized_residues.items():
            prep_file = files.get("prep_file", "Not generated")
            frcmod_file = files.get("frcmod_file", "Not generated")

            if prep_file != "Not generated":
                prep_file = os.path.basename(prep_file)

            if frcmod_file != "Not generated":
                frcmod_file = os.path.basename(frcmod_file)

            _console.print(f"  [green]•[/green] {res_name}: {prep_file} and {frcmod_file}")

    # List pending workflows
    if pending_workflows:
        _console.print("\n[bold yellow]Pending Parameterizations[/bold yellow]")
        _console.print(
            f"[yellow]○ The following {len(pending_workflows)} residue types require Gaussian calculations:[/yellow]"
        )

        for res_name, info in pending_workflows.items():
            missing_files = info.get("missing_files", [])
            _console.print(f"  [yellow]•[/yellow] {res_name}: {len(missing_files)} missing files")
            for file in missing_files[:2]:  # Show just a couple of missing files
                _console.print(f"    [grey50]* {file}[/grey50]")
            if len(missing_files) > 2:
                _console.print(f"    [grey50]* ...and {len(missing_files)-2} more[/grey50]")

        _console.print("\n[bold]To complete these parameterizations:[/bold]")
        _console.print("1. Run the required Gaussian calculations")
        _console.print("2. Re-run the parameterization workflow for these residues")

    # Show tLEaP commands if any parameterizations completed
    if parameterized_residues:
        _console.print(
            "\n[bold]The following tLEaP commands will be used for the parameterized residues:[/bold]"
        )
        for cmd in results["tleap_commands"]:
            _console.print(f"  [cyan]{cmd}[/cyan]")

    # Update results
    results["parameterized_residues"] = parameterized_residues
    results["pending_workflows"] = pending_workflows

    if results["errors"]:
        # Only consider it unsuccessful if we had more errors than successes
        results["success"] = len(results["errors"]) < len(parameterized_residues)

    # Consider overall success if we either have parameterized residues or pending workflows
    if not parameterized_residues and not pending_workflows:
        results["success"] = False
        results["message"] = "No residues were successfully parameterized"

    logger.info(
        f"Completed parameterization of {len(parameterized_residues)} modified residues"
    )
    logger.info(f"Pending parameterization for {len(pending_workflows)} residues")

    return results


def integrate_parameterizer(processor, parameterize_modified=False):
    """
    Integrate parameterizer functionality with the PDB processor.

    Args:
        processor: PDBProcessor instance
        parameterize_modified: Whether to parameterize modified amino acids

    Returns:
        dict: Results of parameterization
    """
    logger.info("Initializing modified amino acid parameterization")

    # Initialize result tracking
    results = {
        "success": True,
        "parameterized_components": {"modified_residues": {}},
        "errors": [],
    }

    # Check if ComponentAnalyzer has been run
    if not hasattr(processor, "component_analyzer"):
        logger.warning(
            "Component analyzer not initialized. Cannot identify components for parameterization."
        )
        return {
            "success": False,
            "error": "Component analyzer must be run before parameterization",
        }

    # Create output directory for parameter files
    output_dir = processor._get_output_path("parameterized_components")
    os.makedirs(output_dir, exist_ok=True)

    # Initialize tLEaP commands dictionary if needed
    if not hasattr(processor, "tleap_commands"):
        processor.tleap_commands = {}

    # Handle modified amino acids parameterization
    if parameterize_modified:
        # Run modified amino acid parameterization
        try:
            modified_results = parameterize_modified_residues(
                processor.component_analyzer,
                processor.current_structure,
                output_dir=os.path.join(output_dir, "modified_residues"),
                processor=processor,
            )

            if modified_results["success"]:
                # Store results
                if "parameterized_residues" in modified_results:
                    results["parameterized_components"]["modified_residues"] = (
                        modified_results["parameterized_residues"]
                    )

                # Add tLEaP commands
                if (
                    "tleap_commands" in modified_results
                    and modified_results["tleap_commands"]
                ):
                    processor.tleap_commands["modified_residues"] = modified_results[
                        "tleap_commands"
                    ]
            else:
                results["errors"].append(
                    f"Modified residue parameterization error: {modified_results.get('error', 'Unknown error')}"
                )

        except Exception as e:
            logger.error(f"Error in modified amino acid parameterization: {str(e)}")
            results["errors"].append(
                f"Modified residue parameterization exception: {str(e)}"
            )
            results["success"] = False

    # Check if there were any errors
    if results["errors"]:
        results["success"] = False

    return results


def check_workflow_can_resume(amino_acid, output_dir):
    """
    Check if a paused workflow can be resumed by verifying if
    workflow state or required Gaussian log files exist.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    output_dir : str
        Directory where the workflow files are stored

    Returns:
    --------
    dict
        A dictionary indicating if the workflow can be resumed and status details
    """
    # Check new workflow_state.json first
    state_file = Path(output_dir) / "workflow_state.json"
    if state_file.exists():
        return {
            "can_resume": True,
            "message": "Workflow state found. Can resume from last completed step.",
            "source": "workflow_state",
        }

    # Fall back to legacy file-based detection
    original_dir = os.getcwd()

    try:
        os.chdir(output_dir)

        log_check = check_for_log_files(amino_acid)

        if log_check["complete"]:
            return {
                "can_resume": True,
                "message": "All required log files found. Workflow can be resumed.",
                "files": log_check["found_files"],
            }
        else:
            return {
                "can_resume": False,
                "message": "Some required log files are still missing.",
                "missing_files": log_check["missing_files"],
                "found_files": log_check.get("found_files", []),
            }
    finally:
        os.chdir(original_dir)


def resume_paused_workflow(amino_acid, output_dir):
    """
    Resume a previously paused workflow if the required log files are now available.

    Parameters:
    -----------
    amino_acid : str
        Three-letter amino acid code
    output_dir : str
        Directory where the workflow files are stored

    Returns:
    --------
    dict
        Results of the resumed workflow
    """
    # Check if workflow can be resumed
    check_result = check_workflow_can_resume(amino_acid, output_dir)

    if not check_result["can_resume"]:
        _console.print(f"\n[yellow]⚠[/yellow] Cannot resume parameterization of [cyan]{amino_acid}[/cyan] yet.")
        _console.print(f"[bold]The following files are still missing:[/bold]")
        for file in check_result["missing_files"]:
            _console.print(f"  [red]✗[/red] {file}")

        return {
            "success": False,
            "message": "Workflow cannot be resumed yet. Missing required log files.",
            "missing_files": check_result["missing_files"],
            "status": "pending_calculations",
        }

    _console.print(f"\n[green]✓[/green] Found all required log files for [cyan]{amino_acid}[/cyan].")
    _console.print(f"[cyan]→[/cyan] Resuming parameterization workflow...")

    # Store original directory
    original_dir = os.getcwd()

    try:
        # Change to the output directory
        os.chdir(output_dir)

        # Run the workflow, now that the log files are available
        # Skip the initial structure generation steps
        return run_workflow(amino_acid=amino_acid, interactive=True)
    finally:
        # Return to original directory
        os.chdir(original_dir)


def check_and_resume_pending_workflows(pending_workflows, processor=None):
    """
    Check all pending workflows and resume any that can be resumed.

    Parameters:
    -----------
    pending_workflows : dict
        Dictionary of pending workflows with their details

    Returns:
    --------
    dict
        Results of checking and resuming workflows
    """
    results = {
        "success": True,
        "resumed": {},
        "still_pending": {},
        "errors": [],
        "parameter_files": {},
        "tleap_commands": [],
    }

    if not pending_workflows:
        _console.print("[yellow]○[/yellow] No pending workflows to check.")
        return results

    _console.print(f"\n[bold]=== Checking Pending Workflows ===[/bold]")
    _console.print(f"Found [cyan]{len(pending_workflows)}[/cyan] pending parameterization workflows:")

    for res_name, info in pending_workflows.items():
        _console.print(f"\n[cyan]→[/cyan] Checking status of [cyan]{res_name}[/cyan] parameterization...")
        output_dir = info.get("output_dir")

        if not output_dir or not os.path.exists(output_dir):
            _console.print(f"  [red]✗[/red] Error: Output directory not found for {res_name}.")
            results["errors"].append(f"Output directory not found for {res_name}")
            continue

        # Check if this workflow can be resumed
        check_result = check_workflow_can_resume(res_name, output_dir)

        if check_result["can_resume"]:
            _console.print(f"  [green]✓[/green] All required files found for {res_name}. Resuming workflow...")

            # Ask user if they want to resume this workflow
            if confirm_with_context(
                processor,
                f"Resume parameterization of {res_name} now?",
                default=True,
                module="Modified AA Parameterizer",
                description=f"Resume parameterization of {res_name}",
            ):
                # Resume the workflow
                resume_result = resume_paused_workflow(res_name, output_dir)

                if resume_result["success"]:
                    # Process resumed workflow results
                    prep_file = resume_result["parameter_files"].get("prep_file")
                    frcmod_file = resume_result["parameter_files"].get("frcmod_file")

                    if prep_file or frcmod_file:
                        results["resumed"][res_name] = {}

                        if prep_file:
                            results["resumed"][res_name]["prep_file"] = prep_file
                            results["tleap_commands"].append(
                                f"loadAmberPrep {prep_file}"
                            )
                            results["parameter_files"][f"{res_name}_prep"] = prep_file

                        if frcmod_file:
                            results["resumed"][res_name]["frcmod_file"] = frcmod_file
                            results["tleap_commands"].append(
                                f"loadAmberParams {frcmod_file}"
                            )
                            results["parameter_files"][
                                f"{res_name}_frcmod"
                            ] = frcmod_file

                        _console.print(
                            f"  [green]✓[/green] Successfully resumed and completed parameterization of {res_name}"
                        )
                    elif resume_result.get("status") in [
                        "paused",
                        "pending_calculations",
                    ]:
                        # Still pending calculations
                        results["still_pending"][res_name] = info
                        _console.print(
                            f"  [yellow]○[/yellow] Parameterization of {res_name} still requires additional steps."
                        )
                    else:
                        # Workflow completed but no parameter files
                        _console.print(
                            f"  [yellow]⚠[/yellow] Workflow completed for {res_name}, but no parameter files were generated."
                        )
                        results["errors"].append(
                            f"Completed workflow for {res_name} did not generate parameter files"
                        )
                else:
                    # Resume failed
                    error = resume_result.get(
                        "error", resume_result.get("message", "Unknown error")
                    )
                    _console.print(f"  [red]✗[/red] Error resuming {res_name} workflow: {error}")
                    results["errors"].append(
                        f"Error resuming {res_name} workflow: {error}"
                    )
                    results["still_pending"][res_name] = info
            else:
                # User chose not to resume
                _console.print(f"  [yellow]○[/yellow] User chose not to resume {res_name} workflow.")
                results["still_pending"][res_name] = info
        else:
            # Cannot resume yet
            _console.print(f"  [yellow]⚠[/yellow] Cannot resume {res_name} workflow yet. Missing files:")
            for file in check_result["missing_files"]:
                _console.print(f"    [red]✗[/red] {file}")
            results["still_pending"][res_name] = info

    # Print summary
    _console.print("\n[bold]=== Pending Workflow Summary ===[/bold]")

    if results["resumed"]:
        _console.print(
            f"[green]✓[/green] Successfully resumed and completed [cyan]{len(results['resumed'])}[/cyan] workflows:"
        )
        for res_name, files in results["resumed"].items():
            prep_file = files.get("prep_file", "Not generated")
            frcmod_file = files.get("frcmod_file", "Not generated")

            if prep_file != "Not generated":
                prep_file = os.path.basename(prep_file)

            if frcmod_file != "Not generated":
                frcmod_file = os.path.basename(frcmod_file)

            _console.print(f"  [green]✓[/green] {res_name}: [cyan]{prep_file}[/cyan] and [cyan]{frcmod_file}[/cyan]")
    else:
        _console.print("[yellow]○[/yellow] No workflows were resumed and completed.")

    if results["still_pending"]:
        _console.print(f"\n[yellow]⚠[/yellow] [cyan]{len(results['still_pending'])}[/cyan] workflows are still pending:")
        for res_name in results["still_pending"].keys():
            _console.print(f"  [yellow]○[/yellow] {res_name}")

    if results["errors"]:
        _console.print(f"\n[red]✗[/red] [cyan]{len(results['errors'])}[/cyan] errors occurred:")
        for error in results["errors"][:3]:  # Show just a few errors to avoid clutter
            _console.print(f"  [red]✗[/red] {error}")
        if len(results["errors"]) > 3:
            _console.print(f"  [grey50]- ...and {len(results['errors'])-3} more[/grey50]")

    # Update success status
    if results["errors"] and not results["resumed"]:
        results["success"] = False

    return results


def main():
    """Entry point for command-line use of the workflow."""
    run_workflow()


if __name__ == "__main__":
    main()
