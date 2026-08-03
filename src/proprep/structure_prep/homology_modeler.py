"""
Homology Modeler

Self-contained worker for building homology models using MODELLER.
Handles sequence alignment (Bio.Align.PairwiseAligner), PIR file creation
with multi-chain and HETATM support, MODELLER execution, and quality assessment.
"""

import logging
import os
import shutil
import sys
import tempfile
from io import StringIO
from typing import Any, Dict, List, Optional, Tuple

from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.PDBIO import PDBIO
from Bio.Align import PairwiseAligner, substitution_matrices

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

logger = logging.getLogger(__name__)

# Optional MODELLER import — no warning here; structure_completeness.py already
# warns at startup, and the menu shows BLOCKED status when unavailable.
# Suppress C-level error output from _modeller.mod_start().
try:
    import ctypes as _ctypes
    _devnull = os.open(os.devnull, os.O_WRONLY)
    _old_stdout, _old_stderr = os.dup(1), os.dup(2)
    try:
        os.dup2(_devnull, 1)
        os.dup2(_devnull, 2)
        import modeller
        from modeller.automodel import AutoModel
        HAS_MODELLER = True
    finally:
        _ctypes.CDLL(None).fflush(None)
        os.dup2(_old_stdout, 1)
        os.dup2(_old_stderr, 2)
        os.close(_devnull)
        os.close(_old_stdout)
        os.close(_old_stderr)
except Exception:
    HAS_MODELLER = False


class HomologyModeler:
    """Worker for building homology models with MODELLER."""

    # Standard amino acid 3-to-1 mapping
    STD_AMINO_ACIDS = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # Common variants
        'MSE': 'M', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
        'HIE': 'H', 'HID': 'H', 'HIP': 'H',
        'CYX': 'C', 'CYM': 'C',
    }

    def __init__(self):
        pass

    @staticmethod
    def extract_template_sequences(pdb_file: str) -> Dict[str, str]:
        """
        Extract chain sequences from PDB file using PPBuilder.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dict mapping chain_id -> amino acid sequence (one-letter codes)
        """
        if not os.path.exists(pdb_file):
            logger.error(f"PDB file not found: {pdb_file}")
            return {}

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("template", pdb_file)
        ppb = PPBuilder()
        sequences = {}

        for model in structure:
            for chain in model:
                chain_id = chain.id
                pp_list = ppb.build_peptides(chain)
                if pp_list:
                    seq = "".join(str(pp.get_sequence()) for pp in pp_list)
                    if seq:
                        sequences[chain_id] = seq
            break  # Only first model

        return sequences

    @staticmethod
    def count_hetatm_per_chain(pdb_file: str) -> Dict[str, Dict[str, int]]:
        """
        Count non-water HETATMs and waters per chain.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Dict mapping chain_id -> {"hetatm": count, "water": count}
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("template", pdb_file)
        counts = {}

        for model in structure:
            for chain in model:
                water_count = 0
                other_hetatm_count = 0
                for residue in chain:
                    if residue.id[0] != " ":  # HETATM residue
                        if residue.resname in ['HOH', 'WAT']:
                            water_count += 1
                        else:
                            other_hetatm_count += 1
                counts[chain.id] = {
                    "hetatm": other_hetatm_count,
                    "water": water_count,
                }
            break  # Only first model

        return counts

    @staticmethod
    def align_sequences(template_seq: str, target_seq: str) -> Dict[str, Any]:
        """
        Global pairwise alignment using PairwiseAligner with BLOSUM62.

        Args:
            template_seq: Template protein sequence (one-letter codes)
            target_seq: Target protein sequence (one-letter codes)

        Returns:
            Dict with keys: template_aligned, target_aligned, score, identity,
            identity_pct, gaps, gaps_pct, length
        """
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        alignments = aligner.align(template_seq, target_seq)
        best = alignments[0]

        # Extract aligned strings
        aligned_strs = str(best).split('\n')
        # PairwiseAligner format: line 0 = seq1, line 1 = match, line 2 = seq2
        # But the format can vary; use the alignment object directly
        template_aligned = ""
        target_aligned = ""

        # Use alignment indices to build aligned strings
        t_idx = 0
        q_idx = 0
        for (t_start, t_end), (q_start, q_end) in zip(best.aligned[0], best.aligned[1]):
            # Add gaps for unaligned positions before this block
            while t_idx < t_start:
                template_aligned += template_seq[t_idx]
                target_aligned += "-"
                t_idx += 1
            while q_idx < q_start:
                template_aligned += "-"
                target_aligned += target_seq[q_idx]
                q_idx += 1
            # Add aligned block
            block_len = t_end - t_start
            template_aligned += template_seq[t_start:t_end]
            target_aligned += target_seq[q_start:q_end]
            t_idx = t_end
            q_idx = q_end

        # Trailing unaligned residues
        while t_idx < len(template_seq):
            template_aligned += template_seq[t_idx]
            target_aligned += "-"
            t_idx += 1
        while q_idx < len(target_seq):
            template_aligned += "-"
            target_aligned += target_seq[q_idx]
            q_idx += 1

        # Calculate stats
        aln_len = len(template_aligned)
        identity = sum(
            1 for a, b in zip(template_aligned, target_aligned)
            if a == b and a != '-'
        )
        gaps = sum(
            1 for a, b in zip(template_aligned, target_aligned)
            if a == '-' or b == '-'
        )

        return {
            "template_aligned": template_aligned,
            "target_aligned": target_aligned,
            "score": best.score,
            "identity": identity,
            "identity_pct": (identity / aln_len * 100) if aln_len > 0 else 0.0,
            "gaps": gaps,
            "gaps_pct": (gaps / aln_len * 100) if aln_len > 0 else 0.0,
            "length": aln_len,
        }

    @staticmethod
    def display_alignment(alignment_result: Dict[str, Any], console: Console,
                          chain_id: str = "") -> None:
        """
        Display alignment with Rich table + full wrapped alignment blocks.

        Args:
            alignment_result: Output from align_sequences()
            console: Rich Console for output
            chain_id: Optional chain label for display
        """
        chain_label = f" (Chain {chain_id})" if chain_id else ""

        # Summary table
        table = Table(title=f"Sequence Alignment{chain_label}")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="yellow")

        table.add_row("Alignment Length", str(alignment_result["length"]))
        table.add_row("Identity",
                       f"{alignment_result['identity']}/{alignment_result['length']} "
                       f"({alignment_result['identity_pct']:.1f}%)")
        table.add_row("Gaps",
                       f"{alignment_result['gaps']}/{alignment_result['length']} "
                       f"({alignment_result['gaps_pct']:.1f}%)")
        table.add_row("Score", f"{alignment_result['score']:.1f}")
        console.print(table)

        # Full alignment display in wrapped blocks (60 chars/line)
        template_aln = alignment_result["template_aligned"]
        target_aln = alignment_result["target_aligned"]
        block_size = 60
        aln_len = len(template_aln)
        lines = []

        for start in range(0, aln_len, block_size):
            end = min(start + block_size, aln_len)
            t_block = template_aln[start:end]
            q_block = target_aln[start:end]

            # Build match line
            match_line = ""
            for a, b in zip(t_block, q_block):
                if a == b and a != '-':
                    match_line += "|"
                elif a != '-' and b != '-':
                    # Check BLOSUM62 similarity
                    match_line += ":"
                else:
                    match_line += " "

            lines.append(f"Template {start+1:>5d}  {t_block}  {end}")
            lines.append(f"               {match_line}")
            lines.append(f"Target   {start+1:>5d}  {q_block}  {end}")
            lines.append("")

        console.print(Panel(
            "\n".join(lines),
            title=f"Alignment Detail{chain_label}",
            border_style="blue",
            expand=False,
        ))

    @staticmethod
    def create_pir_file(
        output_file: str,
        template_pdb: str,
        chain_ids: List[str],
        chain_alignments: Dict[str, Dict[str, str]],
        hetatm_counts: Dict[str, Dict[str, int]],
    ) -> None:
        """
        Create PIR alignment file for MODELLER with multi-chain and HETATM support.

        Args:
            output_file: Path to write PIR file
            template_pdb: Filename of template PDB (basename only)
            chain_ids: Ordered list of chain IDs to include
            chain_alignments: Dict chain_id -> {"template_aligned": str, "target_aligned": str}
            hetatm_counts: Dict chain_id -> {"hetatm": int, "water": int}
        """
        first_chain = chain_ids[0]
        last_chain = chain_ids[-1]

        with open(output_file, 'w') as f:
            # Template entry - use :: to auto-detect residue range from PDB
            f.write(">P1;template\n")
            f.write(f"structure:{template_pdb}::{first_chain}::{last_chain}::::\n")

            # Template sequences with chain separators and HETATM markers
            for i, chain_id in enumerate(chain_ids):
                aln = chain_alignments[chain_id]
                f.write(aln["template_aligned"])
                # Append HETATM dots and water markers
                counts = hetatm_counts.get(chain_id, {"hetatm": 0, "water": 0})
                if counts["hetatm"] > 0:
                    f.write("." * counts["hetatm"])
                if counts["water"] > 0:
                    f.write("w" * counts["water"])
                if i < len(chain_ids) - 1:
                    f.write("/")
            f.write("*\n\n")

            # Target entry
            f.write(">P1;target\n")
            f.write(f"sequence:target:::::::0.00:0.00\n")

            # Target sequences with chain separators and HETATM markers
            for i, chain_id in enumerate(chain_ids):
                aln = chain_alignments[chain_id]
                f.write(aln["target_aligned"])
                # Same HETATM/water markers so MODELLER preserves them
                counts = hetatm_counts.get(chain_id, {"hetatm": 0, "water": 0})
                if counts["hetatm"] > 0:
                    f.write("." * counts["hetatm"])
                if counts["water"] > 0:
                    f.write("w" * counts["water"])
                if i < len(chain_ids) - 1:
                    f.write("/")
            f.write("*\n")

    @staticmethod
    def run_modeller(
        work_dir: str,
        template_pdb: str,
        alignment_file: str,
        console: Optional[Console] = None,
    ) -> Tuple[bool, str, Optional[Any]]:
        """
        Execute MODELLER to build homology model.

        Args:
            work_dir: Working directory containing input files
            template_pdb: Path to template PDB file
            alignment_file: Path to PIR alignment file
            console: Optional console for output

        Returns:
            Tuple of (success, message, structure_or_None)
        """
        if not HAS_MODELLER:
            return False, "MODELLER not available", None

        import modeller
        from modeller.automodel import AutoModel, assess

        try:
            original_dir = os.getcwd()
            os.chdir(work_dir)

            # Set up environment
            env = modeller.Environ()
            env.io.atom_files_directory = ['.']
            env.io.hetatm = True
            env.io.water = True

            env.libs.topology.read(file="${LIB}/top_heav.lib")
            env.libs.parameters.read(file="${LIB}/par.lib")

            # Redirect stdout/stderr
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            captured_output = StringIO()

            try:
                sys.stdout = captured_output
                sys.stderr = captured_output

                mdl = AutoModel(
                    env,
                    alnfile=os.path.basename(alignment_file),
                    knowns="template",
                    sequence="target",
                )
                mdl.starting_model = 1
                mdl.ending_model = 1
                # Published, calibrated assessment scores computed via the
                # MODELLER API (no stdout scraping): normalized DOPE z-score
                # (Shen & Sali 2006) and GA341 (John & Sali 2003), with raw
                # DOPE / molpdf kept for provenance. Results land in mdl.outputs.
                mdl.assess_methods = (assess.normalized_dope, assess.GA341,
                                      assess.DOPE)
                mdl.make()

            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr
                os.chdir(original_dir)

            modeller_output = captured_output.getvalue()

            # Parse output PDB
            output_pdb = os.path.join(work_dir, "target.B99990001.pdb")
            if not os.path.exists(output_pdb):
                return False, "MODELLER output file not found", None

            parser = PDBParser(QUIET=True)
            model_structure = parser.get_structure("homology_model", output_pdb)

            # Calibrated quality assessment from the MODELLER API. The whole
            # model is built de novo here, so the global normalized-DOPE /
            # GA341 scores are the appropriate, calibrated metrics (unlike the
            # gap-filling path, there is no retained experimental region).
            if console:
                assessment = HomologyModeler._extract_assessment(mdl)
                HomologyModeler.display_quality_assessment(
                    assessment, model_structure, console
                )

            return True, output_pdb, model_structure

        except Exception as e:
            return False, f"MODELLER error: {str(e)}", None

    # ------------------------------------------------------------------
    # Quality assessment
    #
    # A homology model is built de novo by threading the target sequence onto
    # the template, so the WHOLE model is modelled. The decisive, calibrated
    # metrics are therefore MODELLER's own global scores, read straight from
    # the API (no stdout scraping, no invented MOLPDF-per-residue tiers):
    #   • Normalized DOPE z-score — Shen & Sali, Protein Sci. 2006, 15:2507
    #       (z <= -1 indicates a likely-correct, native-like fold).
    #   • GA341 — John & Sali, Nucleic Acids Res. 2003, 31:3982 (range 0-1;
    #       >= 0.7 indicates a reliable fold).
    # raw DOPE / molpdf are shown for provenance only (not comparable across
    # systems).
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_assessment(mdl) -> Dict:
        """Read calibrated scores from the MODELLER AutoModel outputs.

        Values come straight from mdl.outputs (populated by assess_methods),
        so a MODELLER output-format change cannot silently break scoring.
        Returns float-or-None for normalized_dope, ga341, dope and molpdf.
        """
        result = {
            'normalized_dope': None,
            'ga341': None,
            'dope': None,
            'molpdf': None,
        }
        try:
            output = mdl.outputs[0]
        except (AttributeError, IndexError, TypeError, KeyError):
            return result

        def _as_float(value):
            # GA341 is returned as a list whose first element is the score.
            if isinstance(value, (tuple, list)):
                value = value[0] if value else None
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        result['normalized_dope'] = _as_float(output.get('Normalized DOPE score'))
        result['ga341'] = _as_float(output.get('GA341 score'))
        result['dope'] = _as_float(output.get('DOPE score'))
        result['molpdf'] = _as_float(output.get('molpdf'))
        return result

    @staticmethod
    def _interpret_normalized_dope(z: float) -> Tuple[str, str]:
        """Return (label, rich_color) for a normalized DOPE z-score."""
        if z <= -1.0:
            return "Native-like fold (z ≤ −1)", "green"
        if z <= 0.0:
            return "Borderline (−1 < z ≤ 0)", "yellow"
        return "Likely incorrect fold (z > 0)", "red"

    @staticmethod
    def _interpret_ga341(score: float) -> Tuple[str, str]:
        """Return (label, rich_color) for a GA341 score."""
        if score >= 0.7:
            return "Reliable fold (≥ 0.7)", "green"
        return "Low reliability (< 0.7)", "yellow"

    @staticmethod
    def _count_protein(structure) -> Tuple[int, int]:
        """Count standard protein residues and their atoms in a structure."""
        n_res = 0
        n_atoms = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " ":  # standard (non-HETATM) residue
                        n_res += 1
                        n_atoms += len(list(residue.get_atoms()))
        return n_res, n_atoms

    @staticmethod
    def display_quality_assessment(assessment: Dict, model_structure,
                                   console: Console) -> None:
        """Display calibrated global scores for a homology model. No invented
        quality grade or MD-readiness verdict."""
        table = Table(title="Homology Model Assessment")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="yellow")
        table.add_column("Interpretation (literature-calibrated)", style="white")

        if model_structure is not None:
            n_res, n_atoms = HomologyModeler._count_protein(model_structure)
            if n_res:
                table.add_row("Protein residues", str(n_res), "")
            if n_atoms:
                table.add_row("Protein atoms", str(n_atoms), "")

        nd = assessment.get('normalized_dope')
        if nd is not None:
            label, color = HomologyModeler._interpret_normalized_dope(nd)
            table.add_row("Normalized DOPE", f"{nd:.2f}",
                          f"[{color}]{label}[/{color}]  (Shen & Sali 2006)")

        ga = assessment.get('ga341')
        if ga is not None:
            label, color = HomologyModeler._interpret_ga341(ga)
            table.add_row("GA341", f"{ga:.3f}",
                          f"[{color}]{label}[/{color}]  (John & Sali 2003)")

        dope = assessment.get('dope')
        if dope is not None:
            table.add_row("DOPE (raw)", f"{dope:.1f}",
                          "[grey50]non-normalized — provenance only[/grey50]")

        molpdf = assessment.get('molpdf')
        if molpdf is not None:
            table.add_row("molpdf", f"{molpdf:.1f}",
                          "[grey50]optimizer-internal — not comparable across "
                          "systems[/grey50]")

        console.print(table)

        if nd is None and ga is None:
            console.print(
                "[yellow]Calibrated assessment scores were not computed for "
                "this model.[/yellow]"
            )
        elif nd is not None and nd > 0:
            console.print(
                "[grey50]A normalized DOPE z-score above 0 suggests the fold "
                "may be unreliable; inspect the model and alignment.[/grey50]"
            )
