"""
SLURM Job Script Generator

Generates SLURM job scripts for AMBER MD simulations on HPC clusters.
Supports single simulations and workflow dependencies.

Author: ProPrep Development Team
Date: 2025-01-20
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
)
import json


def safe_step_slug(name: str) -> str:
    """Sanitize a step name into something safe for directory + output paths.

    Parentheses, brackets, and other shell metachars in step names like
    'Equilibration (10 kcal)' or 'Production MD (100 ns)' caused unquoted
    sbatch invocations and downstream cd/source commands to fail. Strip
    them, collapse whitespace into underscores, and coalesce runs.
    """
    s = re.sub(r"\s+", "_", (name or "").strip())
    s = re.sub(r"[()\[\]{}<>|&;`$!*?\"'\\]", "", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "step"


@dataclass
class SlurmConfig:
    """SLURM job configuration settings."""
    # Job identification
    job_name: str = "amber_simulation"
    output_file: str = "slurm-%j.out"
    error_file: str = "slurm-%j.err"

    # Resource allocation
    partition: str = "gpu"
    nodes: int = 1
    ntasks: int = 1
    cpus_per_task: int = 4
    gpus: int = 0
    gpu_type: str = ""  # e.g., "v100", "a100"
    memory: str = "16G"  # e.g., "16G", "32000M"

    # Time limits
    time_limit: str = "02:00:00"  # HH:MM:SS

    # Email notifications
    email_notifications: bool = False
    email_address: str = ""
    email_types: str = "END,FAIL"  # BEGIN,END,FAIL,ALL

    # Environment
    modules: List[str] = None  # e.g., ["amber/22", "cuda/11.8"]
    conda_env: str = ""
    pre_commands: List[str] = None   # Bash lines emitted before the AMBER call
    post_commands: List[str] = None  # Bash lines emitted after the AMBER call
                                     # (runs regardless of exit code — typical
                                     # use is cleanup like `rm -rf "$TMPDIR"`)

    # SLURM directives
    account: str = ""  # Allocation account
    qos: str = ""  # Quality of Service
    constraint: str = ""  # Node constraints
    exclusive: bool = False
    extra_directives: List[str] = None  # Arbitrary "#SBATCH <line>" passthroughs

    # Explicit resource-class binary (pmemd.cuda / pmemd.MPI); overrides engine
    # inferred from other flags when generating from a cluster profile.
    binary: str = ""

    # Launcher prefix for MPI binaries (e.g. "srun", "srun --mpi=pmi2",
    # "srun --mpi=pmix"). Use srun, not mpirun: under Slurm with
    # task/affinity, mpirun pins all ranks to one core (see the fallback in
    # _build_amber_command). Empty string = no prefix (GPU / serial).
    mpi_launcher: str = ""

    # export KEY=VALUE lines emitted before the AMBER command
    env_vars: List[str] = None

    def __post_init__(self):
        """Initialize mutable defaults."""
        if self.modules is None:
            self.modules = []
        if self.pre_commands is None:
            self.pre_commands = []
        if self.post_commands is None:
            self.post_commands = []
        if self.extra_directives is None:
            self.extra_directives = []
        if self.env_vars is None:
            self.env_vars = []


class SlurmJobGenerator:
    """Generate SLURM job scripts for AMBER simulations."""

    def __init__(self, console: Optional[Console] = None, processor=None):
        """
        Initialize SLURM job generator.

        Args:
            console: Rich Console instance for output
            processor: Optional ProPrep processor for session recording context
        """
        self.console = console or Console()
        self.processor = processor
        self.config = SlurmConfig()

    def interactive_configure(self) -> SlurmConfig:
        """
        Interactive wizard to configure SLURM settings.

        Returns:
            Configured SlurmConfig object
        """
        self.console.print("\n[bold cyan]===== SLURM Job Configuration =====[/bold cyan]\n")
        self.console.print("[grey50]Configure SLURM settings for HPC cluster submission[/grey50]\n")

        # Job identification
        self.console.print("[bold]Job Identification:[/bold]")
        MOD = "SLURM Generator"
        self.config.job_name = prompt_with_context(
            self.processor, "Job name", default=self.config.job_name,
            module=MOD, description="SLURM job name",
        )

        # Resource allocation
        self.console.print("\n[bold]Resource Allocation:[/bold]")

        self.config.partition = prompt_with_context(
            self.processor, "Partition/Queue", default=self.config.partition,
            module=MOD, description="SLURM partition/queue",
        )

        self.config.nodes = int(prompt_with_context(
            self.processor, "Number of nodes", default=str(self.config.nodes),
            module=MOD, description="Number of nodes",
        ))

        self.config.ntasks = int(prompt_with_context(
            self.processor, "Tasks per node", default=str(self.config.ntasks),
            module=MOD, description="Tasks per node",
        ))

        self.config.cpus_per_task = int(prompt_with_context(
            self.processor, "CPUs per task", default=str(self.config.cpus_per_task),
            module=MOD, description="CPUs per task",
        ))

        # GPU configuration
        use_gpu = confirm_with_context(
            self.processor, "Use GPU?", default=True,
            module=MOD, description="Use GPU",
        )
        if use_gpu:
            self.config.gpus = int(prompt_with_context(
                self.processor, "Number of GPUs", default="1",
                module=MOD, description="Number of GPUs",
            ))

            self.config.gpu_type = prompt_with_context(
                self.processor, "GPU type (leave blank for any)", default="",
                module=MOD, description="GPU type (blank for any)",
            )
        else:
            self.config.gpus = 0
            self.config.gpu_type = ""

        self.config.memory = prompt_with_context(
            self.processor, "Memory (e.g., 16G, 32000M)", default=self.config.memory,
            module=MOD, description="SLURM memory allocation",
        )

        # Time limit
        self.console.print("\n[bold]Time Limits:[/bold]")
        self.config.time_limit = prompt_with_context(
            self.processor, "Time limit (HH:MM:SS)", default=self.config.time_limit,
            module=MOD, description="SLURM wall time",
        )

        # Email notifications
        self.console.print("\n[bold]Email Notifications:[/bold]")
        self.config.email_notifications = confirm_with_context(
            self.processor, "Enable email notifications?", default=False,
            module=MOD, description="Enable SLURM email notifications",
        )

        if self.config.email_notifications:
            self.config.email_address = prompt_with_context(
                self.processor, "Email address",
                module=MOD, description="Email address for notifications",
            )
            self.config.email_types = prompt_with_context(
                self.processor, "Notification types (BEGIN,END,FAIL,ALL)",
                default=self.config.email_types,
                module=MOD, description="Email notification event types",
            )

        # Environment setup
        self.console.print("\n[bold]Environment Setup:[/bold]")

        # Modules
        if confirm_with_context(
            self.processor, "Load modules?", default=True,
            module=MOD, description="Load environment modules",
        ):
            modules_input = prompt_with_context(
                self.processor,
                "Modules to load (comma-separated, e.g., 'amber/22,cuda/11.8')",
                default="amber/22",
                module=MOD, description="Modules to load (comma-separated)",
            )
            self.config.modules = [m.strip() for m in modules_input.split(',') if m.strip()]

        # Conda environment
        if confirm_with_context(
            self.processor, "Activate conda environment?", default=False,
            module=MOD, description="Activate conda environment",
        ):
            self.config.conda_env = prompt_with_context(
                self.processor, "Conda environment name",
                module=MOD, description="Conda environment name",
            )

        # Advanced options
        self.console.print("\n[bold]Advanced Options:[/bold]")

        if confirm_with_context(
            self.processor, "Configure advanced SLURM options?", default=False,
            module=MOD, description="Configure advanced SLURM options",
        ):
            self.config.account = prompt_with_context(
                self.processor, "Account/Allocation (leave blank to skip)", default="",
                module=MOD, description="SLURM account/allocation",
            )

            self.config.qos = prompt_with_context(
                self.processor, "Quality of Service (leave blank to skip)", default="",
                module=MOD, description="SLURM QOS",
            )

            self.config.constraint = prompt_with_context(
                self.processor, "Node constraint (leave blank to skip)", default="",
                module=MOD, description="SLURM node constraint",
            )

            self.config.exclusive = confirm_with_context(
                self.processor, "Request exclusive node access?", default=False,
                module=MOD, description="Request exclusive node access",
            )

        # Display summary
        self._display_config_summary()

        if not confirm_with_context(
            self.processor, "\nUse this configuration?", default=True,
            module=MOD, description="Use SLURM configuration",
        ):
            self.console.print("[yellow]Configuration cancelled[/yellow]")
            return None

        return self.config

    def _display_config_summary(self):
        """Display configuration summary."""
        self.console.print("\n[bold cyan]Configuration Summary:[/bold cyan]\n")

        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("Setting", style="white", width=20)
        table.add_column("Value", style="green", width=40)

        table.add_row("Job Name", self.config.job_name)
        table.add_row("Partition", self.config.partition)
        table.add_row("Nodes", str(self.config.nodes))
        table.add_row("Tasks/Node", str(self.config.ntasks))
        table.add_row("CPUs/Task", str(self.config.cpus_per_task))

        if self.config.gpus > 0:
            gpu_str = f"{self.config.gpus}"
            if self.config.gpu_type:
                gpu_str += f" ({self.config.gpu_type})"
            table.add_row("GPUs", gpu_str)

        table.add_row("Memory", self.config.memory)
        table.add_row("Time Limit", self.config.time_limit)

        if self.config.email_notifications:
            table.add_row("Email", f"{self.config.email_address} ({self.config.email_types})")

        if self.config.modules:
            table.add_row("Modules", ", ".join(self.config.modules))

        if self.config.conda_env:
            table.add_row("Conda Env", self.config.conda_env)

        if self.config.account:
            table.add_row("Account", self.config.account)

        if self.config.qos:
            table.add_row("QoS", self.config.qos)

        self.console.print(table)

    def generate_script(
        self,
        sim_config: Any,  # SimulationConfig
        sim_dir: Path,
        engine: str,
        topology_file: Path,
        coordinate_file: Path,
        dependency_jobid: Optional[str] = None,
        extra_flags: Optional[List[str]] = None,
        topology_arg: Optional[str] = None,
        coord_arg: Optional[str] = None,
        reference_arg: Optional[str] = None,
        output_prefix_override: Optional[str] = None,
        script_path: Optional[Path] = None,
    ) -> str:
        """
        Generate SLURM job script for a single simulation.

        Args:
            sim_config: Simulation configuration
            sim_dir: Simulation directory path
            engine: AMBER engine (sander, pmemd, pmemd.cuda)
            topology_file: Path to topology file
            coordinate_file: Path to coordinate file
            dependency_jobid: Optional SLURM job ID this depends on
            extra_flags: Additional command-line flags (e.g. CpHMD -cpin/-cpout/-cprestrt)

        Returns:
            SLURM script content as string
        """
        script_lines = []

        # Shebang
        script_lines.append("#!/bin/bash")

        # SLURM directives
        script_lines.append(f"#SBATCH --job-name={self.config.job_name}")
        script_lines.append(f"#SBATCH --output={self.config.output_file}")
        script_lines.append(f"#SBATCH --error={self.config.error_file}")
        script_lines.append(f"#SBATCH --partition={self.config.partition}")
        script_lines.append(f"#SBATCH --nodes={self.config.nodes}")
        script_lines.append(f"#SBATCH --ntasks={self.config.ntasks}")
        script_lines.append(f"#SBATCH --cpus-per-task={self.config.cpus_per_task}")

        if self.config.gpus > 0:
            if self.config.gpu_type:
                script_lines.append(f"#SBATCH --gres=gpu:{self.config.gpu_type}:{self.config.gpus}")
            else:
                script_lines.append(f"#SBATCH --gres=gpu:{self.config.gpus}")

        script_lines.append(f"#SBATCH --mem={self.config.memory}")
        script_lines.append(f"#SBATCH --time={self.config.time_limit}")

        # Dependency if provided
        if dependency_jobid:
            script_lines.append(f"#SBATCH --dependency=afterok:{dependency_jobid}")

        # Optional directives
        if self.config.account:
            script_lines.append(f"#SBATCH --account={self.config.account}")

        if self.config.qos:
            script_lines.append(f"#SBATCH --qos={self.config.qos}")

        if self.config.constraint:
            script_lines.append(f"#SBATCH --constraint={self.config.constraint}")

        if self.config.exclusive:
            script_lines.append("#SBATCH --exclusive")

        for directive in self.config.extra_directives:
            line = directive.strip()
            if not line:
                continue
            if not line.startswith("#SBATCH"):
                line = f"#SBATCH {line}"
            script_lines.append(line)

        # Email notifications
        if self.config.email_notifications:
            script_lines.append(f"#SBATCH --mail-type={self.config.email_types}")
            script_lines.append(f"#SBATCH --mail-user={self.config.email_address}")

        script_lines.append("")

        # Header comment
        script_lines.append("# AMBER MD Simulation - Generated by ProPrep")
        script_lines.append(f"# Simulation: {sim_config.name}")
        script_lines.append(f"# Engine: {engine}")
        script_lines.append("")

        # Environment setup
        script_lines.append("# Environment Setup")

        # Load modules (preceded by purge so a stale env can't leak in)
        if self.config.modules:
            script_lines.append("module purge")
            for module in self.config.modules:
                script_lines.append(f"module load {module}")
            script_lines.append("")

        # Activate conda environment
        if self.config.conda_env:
            script_lines.append(f"conda activate {self.config.conda_env}")
            script_lines.append("")

        # Pre-AMBER commands (env vars, scratch-dir setup, staging, ...)
        if self.config.pre_commands:
            script_lines.append("# Pre-AMBER commands")
            for cmd in self.config.pre_commands:
                script_lines.append(cmd)
            script_lines.append("")

        # Resolve the submission-time directory so the generated bundle stays
        # valid after being rsync'd / copied to a cluster. Under sbatch, SLURM
        # copies the script into /var/spool/slurmd/..., so ${BASH_SOURCE[0]}
        # resolves there instead of the user's submission directory — that
        # breaks the "relative to the script" trick. $SLURM_SUBMIT_DIR is
        # populated by SLURM with the directory sbatch was invoked from, and
        # that's the directory we actually want. Fall back to BASH_SOURCE for
        # direct (non-SLURM) shell runs. If we don't know where the script
        # will land (no script_path passed), fall back to the absolute
        # sim_dir — works only if the tree is not moved after generation.
        script_lines.append("# Change into the simulation directory.")
        script_lines.append("# Under sbatch, SLURM_SUBMIT_DIR is the directory where sbatch was invoked")
        script_lines.append("# (BASH_SOURCE would point at the spooled copy under /var, which is wrong).")
        script_lines.append('SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"')
        if script_path is not None:
            try:
                rel = Path(sim_dir).resolve().relative_to(Path(script_path).resolve().parent)
                rel_str = "." if str(rel) == "." else str(rel)
            except ValueError:
                # sim_dir isn't inside the script's directory — fall back
                # to the absolute sim_dir and warn via a comment.
                rel_str = None
        else:
            rel_str = None

        if rel_str is None:
            script_lines.append(f"cd {Path(sim_dir).absolute()}")
        elif rel_str == ".":
            script_lines.append('cd "$SCRIPT_DIR"')
        else:
            script_lines.append(f'cd "$SCRIPT_DIR/{rel_str}"')
        script_lines.append("")

        # Environment variables declared on the resource class.
        if self.config.env_vars:
            script_lines.append("# Environment variables")
            for kv in self.config.env_vars:
                kv = kv.strip()
                if not kv:
                    continue
                if kv.startswith("export "):
                    script_lines.append(kv)
                else:
                    script_lines.append(f"export {kv}")
            script_lines.append("")

        # Print start message
        script_lines.append("# Start message")
        script_lines.append('echo "========================================="')
        script_lines.append(f'echo "Starting AMBER simulation: {sim_config.name}"')
        script_lines.append('echo "Job ID: $SLURM_JOB_ID"')
        script_lines.append('echo "Node: $SLURM_JOB_NODELIST"')
        script_lines.append('echo "Start time: $(date)"')
        script_lines.append('echo "========================================="')
        script_lines.append("")

        # AMBER command
        output_prefix = output_prefix_override or safe_step_slug(sim_config.name).lower()
        mdin_file = sim_dir / "simulation.mdin"
        mdout_file = sim_dir / f"{output_prefix}.mdout"
        restart_file = sim_dir / f"{output_prefix}.rst7"
        trajectory_file = sim_dir / f"{output_prefix}.nc"
        mdinfo_file = sim_dir / f"{output_prefix}.mdinfo"

        # Path strings to hand to pmemd. When the caller supplies explicit
        # overrides (chained workflow with files staged at the parent), use
        # them verbatim; otherwise fall back to the basename of the staged
        # file inside sim_dir, matching legacy behaviour.
        topology_path_str = topology_arg or topology_file.name
        coord_path_str = coord_arg or coordinate_file.name
        reference_path_str = reference_arg or coord_path_str

        script_lines.append("# Run AMBER simulation")
        # MPI binaries need a launcher prefix. The resource class (or profile
        # default) specifies it; fall back to a sane default when absent so
        # legacy manual-wizard runs keep working.
        #
        # The fallback is `srun` (Slurm's native launcher), NOT `mpirun`.
        # Under Slurm with MPICH/Hydra, `mpirun` launches its proxy as a single
        # Slurm task (srun -n1 hydra_pmi_proxy); with TaskPlugin=task/affinity
        # that one task is pinned to one core, and every rank the proxy forks
        # inherits that single-core cpuset — so an N-rank job time-shares one
        # core (~N x slowdown). `srun` makes Slurm bind each rank to its own
        # core. Site-specific PMI bootstrap (e.g. `srun --mpi=pmi2` for
        # MPICH-linked AMBER, `srun --mpi=pmix` for OpenMPI-linked AMBER) is
        # supplied per cluster via the profile's `mpi_launcher`, never here.
        amber_cmd: List[str] = []
        launcher = (self.config.mpi_launcher or "").strip()
        if launcher:
            amber_cmd.extend(launcher.split())
        elif engine.endswith(".MPI") or engine.lower().endswith("mpi"):
            amber_cmd.extend(["srun"])
        amber_cmd.extend([
            engine,
            "-O",
            "-i", str(mdin_file.name),
            "-p", topology_path_str,
            "-c", coord_path_str,
            "-o", str(mdout_file.name),
            "-r", str(restart_file.name),
            "-x", str(trajectory_file.name),
            "-ref", reference_path_str,
        ])

        # Add -inf flag for pmemd variants
        if 'pmemd' in engine:
            amber_cmd.extend(["-inf", str(mdinfo_file.name)])

        # Append extra flags (e.g. CpHMD -cpin/-cpout/-cprestrt)
        if extra_flags:
            amber_cmd.extend(extra_flags)

        # Format command with line continuation for readability
        script_lines.append(" \\\n  ".join(amber_cmd))
        script_lines.append("")

        # Capture AMBER exit status before running any post-AMBER commands so
        # the original exit code survives whatever post-commands return.
        script_lines.append("# Capture AMBER exit status")
        script_lines.append("EXIT_CODE=$?")
        script_lines.append('if [ $EXIT_CODE -eq 0 ]; then')
        script_lines.append('    echo "Simulation completed successfully at $(date)"')
        script_lines.append('else')
        script_lines.append('    echo "Simulation failed with exit code $EXIT_CODE at $(date)"')
        script_lines.append('fi')

        # Post-AMBER commands (e.g. scratch-dir cleanup). Run regardless of
        # success or failure so cleanup is guaranteed; the AMBER exit code is
        # preserved and re-exposed at the final `exit` below.
        if self.config.post_commands:
            script_lines.append("")
            script_lines.append("# Post-AMBER commands (run regardless of AMBER exit code)")
            for cmd in self.config.post_commands:
                script_lines.append(cmd)

        script_lines.append("")
        script_lines.append("exit $EXIT_CODE")

        return "\n".join(script_lines)

    def generate_workflow_scripts(
        self,
        workflow_sims: List[Dict[str, Any]],
        output_dir: Path
    ) -> Dict[str, Path]:
        """
        Generate SLURM scripts for a complete workflow with dependencies.

        Args:
            workflow_sims: List of simulation dicts with configs and paths
            output_dir: Directory to save scripts

        Returns:
            Dictionary mapping simulation names to script paths
        """
        script_paths = {}

        for i, sim_info in enumerate(workflow_sims):
            sim_config = sim_info['config']
            sim_dir = sim_info['sim_dir']
            # step_key drives script filename, --job-name, and the JOBID env
            # var used for dependency chaining. Callers pass 'step1', 'step2',
            # ... so SLURM artifacts match batch-path naming.
            step_key = sim_info.get('step_key') or safe_step_slug(sim_config.name)

            # Determine dependency
            dependency_var = None
            if i > 0:
                prev_info = workflow_sims[i-1]
                prev_key = prev_info.get('step_key') or safe_step_slug(prev_info['config'].name)
                dependency_var = f"${prev_key.upper()}_JOBID"

            # Update job name for this specific simulation
            original_job_name = self.config.job_name
            self.config.job_name = step_key

            # Decide where the script will live so generate_script can emit
            # a SCRIPT_DIR-relative `cd` that survives tree relocation.
            script_filename = f"{step_key}.sh"
            script_path = output_dir / script_filename

            # Generate script
            script_content = self.generate_script(
                sim_config=sim_config,
                sim_dir=sim_dir,
                engine=sim_config.engine,
                topology_file=Path(sim_config.prmtop),
                coordinate_file=Path(sim_config.rst7),
                dependency_jobid=dependency_var,
                extra_flags=sim_info.get('extra_flags'),
                topology_arg=sim_info.get('topology_arg'),
                coord_arg=sim_info.get('coord_arg'),
                reference_arg=sim_info.get('reference_arg'),
                output_prefix_override=sim_info.get('output_prefix'),
                script_path=script_path,
            )

            # Restore original job name
            self.config.job_name = original_job_name

            with open(script_path, 'w') as f:
                f.write(script_content)

            # Make executable
            script_path.chmod(0o755)

            script_paths[sim_config.name] = script_path

        return script_paths

    def generate_submission_script(
        self,
        workflow_sims: List[Dict[str, Any]],
        script_paths: Dict[str, Path],
        output_dir: Path
    ) -> Path:
        """
        Generate master submission script with dependency chaining.

        Args:
            workflow_sims: List of simulation dicts
            script_paths: Dictionary of simulation names to script paths
            output_dir: Output directory

        Returns:
            Path to submission script
        """
        lines = []
        lines.append("#!/bin/bash")
        lines.append("")
        lines.append("# Master SLURM submission script - Generated by ProPrep")
        lines.append("# Submits workflow with dependency chaining")
        lines.append("")

        for i, sim_info in enumerate(workflow_sims):
            sim_config = sim_info['config']
            sim_name = sim_config.name
            step_key = sim_info.get('step_key') or safe_step_slug(sim_name)
            script_path = script_paths[sim_name]
            # Path relative to output_dir (where submit_workflow.sh lives),
            # quoted so names with special chars can't break the sbatch line.
            try:
                rel_script = script_path.relative_to(output_dir)
            except ValueError:
                rel_script = script_path
            script_ref = f'"./{rel_script}"' if not str(rel_script).startswith('/') else f'"{rel_script}"'

            # Variable name for this job ID — step_key based so it's shell-safe
            # and matches script/file naming (STEP1_JOBID, STEP2_JOBID, ...).
            var_name = f"{step_key.upper()}_JOBID"

            if i == 0:
                # First job - no dependency
                lines.append(f"# Submit {sim_name}")
                lines.append(f'{var_name}=$(sbatch --parsable {script_ref})')
                lines.append(f'echo "Submitted {sim_name}: Job ID ${var_name}"')
            else:
                # Dependent job
                prev_info = workflow_sims[i-1]
                prev_key = prev_info.get('step_key') or safe_step_slug(prev_info['config'].name)
                prev_var = f"{prev_key.upper()}_JOBID"
                lines.append("")
                lines.append(f"# Submit {sim_name} (depends on previous job)")
                lines.append(f'{var_name}=$(sbatch --parsable --dependency=afterok:${prev_var} {script_ref})')
                lines.append(f'echo "Submitted {sim_name}: Job ID ${var_name} (depends on ${prev_var})"')

        lines.append("")
        lines.append('echo "All jobs submitted successfully"')

        # Save submission script
        submit_script = output_dir / "submit_workflow.sh"
        with open(submit_script, 'w') as f:
            f.write("\n".join(lines))

        submit_script.chmod(0o755)

        return submit_script

    def save_config(self, output_path: Path):
        """Save SLURM configuration to JSON file."""
        config_dict = asdict(self.config)
        with open(output_path, 'w') as f:
            json.dump(config_dict, f, indent=2)

    def load_config(self, config_path: Path):
        """Load SLURM configuration from JSON file."""
        with open(config_path, 'r') as f:
            config_dict = json.load(f)
        self.config = SlurmConfig(**config_dict)


def build_config_from_profile(
    cluster_profile: Any,
    resource_class_name: str,
    job_name: str,
    time_limit: str = "",
    overrides: Optional[Dict[str, Any]] = None,
) -> SlurmConfig:
    """Resolve a (profile, class, overrides) triple into a concrete SlurmConfig.

    Per-step overrides layer over the class entry, which layers over cluster-
    wide config (modules, account, etc.). Empty / falsy overrides are ignored
    so they never wipe out class defaults.
    """
    cls = cluster_profile.resolve_class(resource_class_name)
    ov = {k: v for k, v in (overrides or {}).items() if v not in (None, "", [])}
    merged = {**cls, **ov}

    cluster_cfg = cluster_profile.cluster or {}
    gpus = int(merged.get("gpus", 0) or 0)
    gres = merged.get("gres", "")
    gpu_type = ""
    if gres:
        # expect 'gpu:<type>:{n}' or 'gpu:<type>:<n>'
        parts = gres.split(":")
        if len(parts) >= 3:
            gpu_type = parts[1]

    effective_time = (
        time_limit
        or merged.get("time_limit")
        or merged.get("default_time")
        or "24:00:00"
    )

    output_pattern = merged.get("output_pattern", "{job_name}.%j.out")
    error_pattern = merged.get("error_pattern", "{job_name}.%j.err")

    return SlurmConfig(
        job_name=job_name,
        output_file=output_pattern.format(job_name=job_name),
        error_file=error_pattern.format(job_name=job_name),
        partition=merged.get("partition", ""),
        nodes=int(merged.get("nodes", 1) or 1),
        ntasks=int(merged.get("ntasks", 1) or 1),
        cpus_per_task=int(merged.get("cpus_per_task", 1) or 1),
        gpus=gpus,
        gpu_type=gpu_type,
        memory=merged.get("memory", ""),
        time_limit=effective_time,
        email_notifications=bool(cluster_cfg.get("email_address")),
        email_address=cluster_cfg.get("email_address", ""),
        email_types=cluster_cfg.get("email_types", "END,FAIL"),
        modules=list(cluster_cfg.get("modules", [])),
        conda_env=cluster_cfg.get("conda_env", ""),
        pre_commands=list(cluster_cfg.get("pre_commands", [])),
        post_commands=list(cluster_cfg.get("post_commands", [])),
        account=cluster_cfg.get("account", ""),
        qos=cluster_cfg.get("qos", ""),
        constraint=merged.get("constraint", ""),
        exclusive=bool(merged.get("exclusive", False)),
        extra_directives=list(cluster_cfg.get("extra_directives", [])),
        binary=merged.get("binary", ""),
        mpi_launcher=merged.get("mpi_launcher", ""),
        env_vars=list(merged.get("env_vars", [])),
    )


def generate_scripts_from_plan(
    cluster_profile: Any,
    run_plan: Any,
    workflow_sims: List[Dict[str, Any]],
    output_dir: Path,
    console: Optional[Console] = None,
) -> Dict[str, Path]:
    """Generate one SLURM script per step using the profile+plan pairing.

    workflow_sims is the list the existing generate_workflow_scripts consumes:
    each item is {config: SimulationConfig, sim_dir: Path, extra_flags?: list}.
    The simulation's `config.step_id` (or `config.name`) is used to look up the
    run plan entry. Scripts are written to output_dir and a master submission
    script with dependency chaining is written alongside them.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    generator = SlurmJobGenerator(console=console)
    script_paths: Dict[str, Path] = {}

    for i, sim_info in enumerate(workflow_sims):
        sim_config = sim_info["config"]
        sim_dir = sim_info["sim_dir"]
        step_id = getattr(sim_config, "step_id", None) or sim_config.name
        if step_id not in run_plan.step_resources:
            raise KeyError(
                f"Run plan '{run_plan.name}' has no assignment for step '{step_id}'"
            )
        sr = run_plan.step_resources[step_id]
        # step_key drives SLURM --job-name, script filename, and JOBID var
        # (step1, step2, ...). Falls back to the sanitized display name when
        # callers haven't supplied one, so older call sites keep working.
        step_key = sim_info.get('step_key') or safe_step_slug(sim_config.name)

        generator.config = build_config_from_profile(
            cluster_profile=cluster_profile,
            resource_class_name=sr.class_name,
            job_name=step_key,
            time_limit=sr.time_limit,
            overrides=sr.overrides,
        )

        engine = generator.config.binary or getattr(sim_config, "engine", "pmemd")
        dependency_var = None
        if i > 0:
            prev_info = workflow_sims[i - 1]
            prev_key = prev_info.get('step_key') or safe_step_slug(prev_info["config"].name)
            dependency_var = f"${prev_key.upper()}_JOBID"

        # Merge run-plan amber_flags with any sim-specific extras
        # (e.g. CpHMD flags from the workflow setup).
        merged_flags = list(sim_info.get("extra_flags") or [])
        if sr.amber_flags:
            merged_flags.extend(str(f) for f in sr.amber_flags)

        script_path = output_dir / f"{step_key}.sh"
        script_content = generator.generate_script(
            sim_config=sim_config,
            sim_dir=sim_dir,
            engine=engine,
            topology_file=Path(sim_config.prmtop),
            coordinate_file=Path(sim_config.rst7),
            dependency_jobid=dependency_var,
            extra_flags=merged_flags or None,
            topology_arg=sim_info.get('topology_arg'),
            coord_arg=sim_info.get('coord_arg'),
            reference_arg=sim_info.get('reference_arg'),
            output_prefix_override=sim_info.get('output_prefix'),
            script_path=script_path,
        )
        with open(script_path, "w") as fh:
            fh.write(script_content)
        script_path.chmod(0o755)
        script_paths[sim_config.name] = script_path

    generator.generate_submission_script(workflow_sims, script_paths, output_dir)
    return script_paths
