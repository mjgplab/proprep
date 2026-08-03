"""
RESP Runner

Executes AmberTools resp program for RESP charge fitting.
Runs two-stage RESP fitting and parses output charges.

Based on MCPB's resp execution implementation.
"""

from pathlib import Path
from typing import Optional, Dict, List, Tuple
import subprocess
import os
import logging


class RESPRunner:
    """
    Execute AmberTools resp program for RESP charge fitting.

    Runs two-stage RESP:
    - Stage 1: resp1.in → resp1.out, resp1.chg
    - Stage 2: resp2.in → resp2.out, resp2.chg (uses resp1.chg as qin)
    """

    def __init__(self, amberhome: Optional[str] = None, logger=None):
        """
        Initialize RESP runner.

        Args:
            amberhome: Path to AMBER installation (auto-detect if None)
            logger: Optional logger instance

        Raises:
            RuntimeError: If AMBER installation not found
        """
        self.amberhome = amberhome or self._detect_amberhome()
        self.logger = logger or logging.getLogger(__name__)
        self.resp_exe = self._find_resp_executable()

    def _detect_amberhome(self) -> str:
        """
        Detect AMBERHOME from environment.

        Returns:
            Path to AMBER installation

        Raises:
            RuntimeError: If AMBERHOME not found
        """
        amberhome = os.getenv('AMBERHOME')
        if not amberhome:
            raise RuntimeError(
                "AMBERHOME environment variable not set. "
                "Please source amber.sh or set AMBERHOME manually."
            )
        return amberhome

    def _find_resp_executable(self) -> str:
        """
        Find resp executable in AMBER installation.

        Returns:
            Path to resp executable

        Raises:
            RuntimeError: If resp executable not found
        """
        # Try AMBERHOME/bin/resp
        resp_path = Path(self.amberhome) / "bin" / "resp"
        if resp_path.exists() and os.access(resp_path, os.X_OK):
            return str(resp_path)

        # Try which resp
        try:
            result = subprocess.run(['which', 'resp'], capture_output=True, text=True)
            if result.returncode == 0:
                resp_path = result.stdout.strip()
                if resp_path:
                    return resp_path
        except Exception:
            pass

        raise RuntimeError(
            f"resp executable not found in {self.amberhome}/bin or PATH. "
            "Please ensure AmberTools is properly installed."
        )

    def run_resp_fitting(self,
                        esp_file: str,
                        respin1: str,
                        respin2: str,
                        output_dir: str,
                        timeout: int = 300) -> Dict[str, str]:
        """
        Run two-stage RESP fitting.

        MCPB command format (resp_fitting.py lines 446-449):
        Stage 1: resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -e esp.esp -s resp1_calc.esp
        Stage 2: resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg -t resp2.chg -e esp.esp -s resp2_calc.esp

        Args:
            esp_file: Path to .esp file
            respin1: Path to stage 1 input (resp1.in)
            respin2: Path to stage 2 input (resp2.in)
            output_dir: Directory for output files
            timeout: Timeout in seconds per stage

        Returns:
            Dict with paths to output files:
            {
                'stage1_output': 'resp1.out',
                'stage1_charges': 'resp1.chg',
                'stage2_output': 'resp2.out',
                'stage2_charges': 'resp2.chg'
            }

        Raises:
            RuntimeError: If resp execution fails
            TimeoutError: If execution exceeds timeout
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Run Stage 1
        self.logger.debug("Running RESP stage 1...")
        stage1_out, stage1_chg = self._run_resp_stage(
            stage=1,
            esp_file=esp_file,
            respin=respin1,
            qin_file=None,
            output_dir=output_dir,
            timeout=timeout
        )
        self.logger.debug(f"  Stage 1 complete: {Path(stage1_chg).name}")

        # Run Stage 2 (uses stage 1 charges as input)
        self.logger.debug("Running RESP stage 2...")
        stage2_out, stage2_chg = self._run_resp_stage(
            stage=2,
            esp_file=esp_file,
            respin=respin2,
            qin_file=stage1_chg,
            output_dir=output_dir,
            timeout=timeout
        )
        self.logger.debug(f"  Stage 2 complete: {Path(stage2_chg).name}")

        return {
            'stage1_output': stage1_out,
            'stage1_charges': stage1_chg,
            'stage2_output': stage2_out,
            'stage2_charges': stage2_chg
        }

    def _run_resp_stage(self,
                       stage: int,
                       esp_file: str,
                       respin: str,
                       qin_file: Optional[str],
                       output_dir: Path,
                       timeout: int) -> Tuple[str, str]:
        """
        Run single RESP stage.

        Args:
            stage: Stage number (1 or 2)
            esp_file: Path to .esp file
            respin: Path to respin file
            qin_file: Path to qin file (None for stage 1, resp1.chg for stage 2)
            output_dir: Output directory
            timeout: Timeout in seconds

        Returns:
            Tuple of (respout_path, respchg_path)

        Raises:
            RuntimeError: If resp execution fails
            TimeoutError: If execution exceeds timeout
        """
        respout = output_dir / f"resp{stage}.out"
        respchg = output_dir / f"resp{stage}.chg"
        resppch = output_dir / f"resp{stage}.pch"
        respcalc = output_dir / f"resp{stage}_calc.esp"

        # Build command using relative paths (resp.F has 80-char filename limit)
        def _rel(path, base):
            """Get relative path from base, or filename if in same directory."""
            try:
                return str(Path(path).relative_to(base))
            except ValueError:
                return str(path)

        cmd = [
            self.resp_exe,
            "-O",
            "-i", _rel(respin, output_dir),
            "-o", _rel(respout, output_dir),
            "-p", _rel(resppch, output_dir),
            "-t", _rel(respchg, output_dir),
            "-e", _rel(esp_file, output_dir),
            "-s", _rel(respcalc, output_dir)
        ]

        # Add qin for stage 2
        if qin_file:
            cmd.extend(["-q", _rel(qin_file, output_dir)])

        self.logger.debug(f"Running: {' '.join(cmd)}")

        # Execute
        try:
            result = subprocess.run(
                cmd,
                cwd=str(output_dir),
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                error_msg = f"resp stage {stage} failed with return code {result.returncode}\n"
                error_msg += f"STDOUT:\n{result.stdout}\n"
                error_msg += f"STDERR:\n{result.stderr}"
                raise RuntimeError(error_msg)

            # Check output files exist
            if not respout.exists() or not respchg.exists():
                raise RuntimeError(
                    f"resp stage {stage} did not generate expected output files. "
                    f"Check {respout} for errors."
                )

            return str(respout), str(respchg)

        except subprocess.TimeoutExpired:
            raise TimeoutError(
                f"resp stage {stage} exceeded timeout ({timeout}s). "
                "Consider increasing timeout or checking input files."
            )

    def parse_final_charges(self, punch_file: str) -> List[float]:
        """
        Parse final charges from resp2.chg file.

        Format: One or more charges per line, whitespace-separated

        Args:
            punch_file: Path to resp2.chg file

        Returns:
            List of RESP-fitted charges in atom order
        """
        charges = []

        with open(punch_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    # Split on whitespace and convert to float
                    for value in line.split():
                        try:
                            charges.append(float(value))
                        except ValueError:
                            pass  # Skip non-numeric values

        self.logger.debug(f"Parsed {len(charges)} charges from {Path(punch_file).name}")
        return charges
