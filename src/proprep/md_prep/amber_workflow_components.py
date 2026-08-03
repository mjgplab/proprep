"""
AMBER Workflow Components

Complete components extracted from amber_workflow.py for integration into ProPrep.
Retains all original functionality while adapting to ProPrep architecture.
"""

import os
import subprocess
import time
import json
import re
import threading
# NOTE: Do NOT import readline here. While it enables command history for input(),
# it conflicts with Rich's ANSI escape sequence handling and causes prompts to
# vanish when the user types a character and presses backspace. See issue investigation
# in test_readline_culprit.py for details.
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from collections import defaultdict, deque

from proprep.utils.prompts import prompt_with_context


class AMBERMonitor:
    """Monitor AMBER simulation output for real-time analysis."""
    
    def __init__(self, output_file: str, max_points: int = 50):
        self.output_file = output_file
        self.info_file = output_file.replace('.out', '.info')
        self.max_points = max_points
        self.data = defaultdict(deque)
        self.timing_data = {}
        self.last_position = 0
        self.last_info_position = 0

        # System information
        self.system_info = {}  # Store NATOM, NTYPES, etc.
        self.mdin_params = {}  # Store imin, ntt, ntp, ntb, tempi, temp0
        self.simulation_stage = None  # User-confirmed or inferred stage
        
    @staticmethod
    def create_historical_monitor(output_file: str, max_points: int = 100000):
        """
        Create a monitor for analyzing completed simulations.

        Args:
            output_file: Path to mdout file
            max_points: Maximum data points to store (default 100000 for complete analysis)
                       This allows ~2M steps at 20ps snapshots (40μs simulation)
        """
        monitor = AMBERMonitor(output_file, max_points)
        monitor.last_position = 0  # Start from beginning for complete analysis
        monitor.parse_amber_output()
        return monitor
        
    def parse_amber_output(self):
        """Parse AMBER output file for energy, temperature, pressure data."""
        if not os.path.exists(self.output_file):
            return
            
        try:
            with open(self.output_file, 'r') as f:
                f.seek(self.last_position)
                content = f.read()
                self.last_position = f.tell()
                
            # Parse AMBER output lines
            lines = content.split('\n')
            for i, line in enumerate(lines):
                # Stop parsing at summary sections (these are averages, not actual data points)
                if 'A V E R A G E S' in line or 'FINAL RESULTS' in line or 'TIMINGS' in line:
                    break

                # Handle MD output format
                if 'NSTEP' in line and 'TIME(PS)' in line:
                    # Parse main line: NSTEP = 10000 TIME(PS) = 20.000 TEMP(K) = 6.02 PRESS = 0.0
                    nstep_match = re.search(r'NSTEP\s*=\s*(\d+)', line)
                    time_match = re.search(r'TIME\(PS\)\s*=\s*([\d.]+)', line)
                    temp_match = re.search(r'TEMP\(K\)\s*=\s*([\d.-]+)', line)
                    press_match = re.search(r'PRESS\s*=\s*([\d.-]+)', line)
                    
                    if nstep_match and time_match:
                        step = int(nstep_match.group(1))
                        time_ps = float(time_match.group(1))
                        
                        self.data['step'].append(step)
                        self.data['time'].append(time_ps)
                        
                        if temp_match:
                            self.data['temperature'].append(float(temp_match.group(1)))
                        if press_match:
                            self.data['pressure'].append(float(press_match.group(1)))
                            
                elif 'Etot' in line:
                    # Parse MD energy line: Etot = -266053.2073 EKtot = 824.4644 EPtot = -266877.6717
                    etot_match = re.search(r'Etot\s*=\s*([\d.-]+)', line)
                    ektot_match = re.search(r'EKtot\s*=\s*([\d.-]+)', line)
                    eptot_match = re.search(r'EPtot\s*=\s*([\d.-]+)', line)
                    
                    if etot_match:
                        self.data['total_energy'].append(float(etot_match.group(1)))
                    if ektot_match:
                        self.data['kinetic_energy'].append(float(ektot_match.group(1)))
                    if eptot_match:
                        self.data['potential_energy'].append(float(eptot_match.group(1)))
                        
                # Handle minimization output format
                elif line.strip() and re.match(r'\s*\d+\s+[\d.-]+E[+-]?\d+\s+[\d.E+-]+\s+[\d.E+-]+', line):
                    # Parse minimization table: NSTEP  ENERGY   RMS   GMAX   NAME  NUMBER
                    #                          1      -1.5166E+05  1.0223E+03  2.2397E+05  HD21  5639
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            step = int(parts[0])
                            energy = float(parts[1])
                            rms = float(parts[2])
                            gmax = float(parts[3])

                            self.data['step'].append(step)
                            self.data['time'].append(float(step))  # Use step as pseudo-time for minimization
                            self.data['total_energy'].append(energy)
                            self.data['rms_gradient'].append(rms)
                            self.data['gmax'].append(gmax)  # Track GMAX for minimization

                        except (ValueError, IndexError):
                            continue

                # Parse density and volume (NPT simulations)
                elif 'VOLUME' in line and '=' in line:
                    # VOLUME     =    479978.9347
                    volume_match = re.search(r'VOLUME\s+=\s+([\d.]+)', line)
                    if volume_match:
                        self.data['volume'].append(float(volume_match.group(1)))

                elif 'Density' in line and '=' in line:
                    # Density    =         0.8589
                    density_match = re.search(r'Density\s+=\s+([\d.]+)', line)
                    if density_match:
                        self.data['density'].append(float(density_match.group(1)))

                # Parse system information (NATOM)
                elif 'NATOM' in line and '=' in line and not self.system_info.get('natom'):
                    # NATOM  =   40754 NTYPES =      17 NBONH =   39568 MBONA  =    1178
                    natom_match = re.search(r'NATOM\s+=\s+(\d+)', line)
                    if natom_match:
                        self.system_info['natom'] = int(natom_match.group(1))

                # Parse mdin parameters
                elif 'imin=' in line:
                    imin_match = re.search(r'imin\s*=\s*(\d+)', line)
                    if imin_match:
                        self.mdin_params['imin'] = int(imin_match.group(1))

                elif 'ntt=' in line:
                    ntt_match = re.search(r'ntt\s*=\s*(\d+)', line)
                    if ntt_match:
                        self.mdin_params['ntt'] = int(ntt_match.group(1))

                elif 'ntp=' in line:
                    ntp_match = re.search(r'ntp\s*=\s*(\d+)', line)
                    if ntp_match:
                        self.mdin_params['ntp'] = int(ntp_match.group(1))

                elif 'ntb=' in line:
                    ntb_match = re.search(r'ntb\s*=\s*(\d+)', line)
                    if ntb_match:
                        self.mdin_params['ntb'] = int(ntb_match.group(1))

                elif 'tempi=' in line:
                    tempi_match = re.search(r'tempi\s*=\s*([\d.]+)', line)
                    if tempi_match:
                        self.mdin_params['tempi'] = float(tempi_match.group(1))

                elif 'temp0=' in line:
                    temp0_match = re.search(r'temp0\s*=\s*([\d.]+)', line)
                    if temp0_match:
                        self.mdin_params['temp0'] = float(temp0_match.group(1))
                            
            # Keep only last max_points
            for key in self.data:
                while len(self.data[key]) > self.max_points:
                    self.data[key].popleft()
                    
        except Exception as e:
            print(f"Error parsing output: {e}")

    def parse_info_file(self):
        """Parse AMBER .info file for timing and performance data."""
        if not os.path.exists(self.info_file):
            return
            
        try:
            # Always read the entire .info file since AMBER overwrites timing sections
            with open(self.info_file, 'r') as f:
                content = f.read()
                
            # Look for the most recent timing block
            timing_blocks = content.split('Current Timing Info')
            if len(timing_blocks) > 1:
                # Get the last timing block
                last_block = timing_blocks[-1]
                
                # Parse completion percentage
                completed_match = re.search(r'Completed:\s+(\d+)\s+\(\s*([\d.]+)%\)', last_block)
                if completed_match:
                    completed_steps = int(completed_match.group(1))
                    completion_percent = float(completed_match.group(2))
                    self.timing_data['completed_steps'] = completed_steps
                    self.timing_data['completion_percent'] = completion_percent
                
                # Parse total steps
                total_match = re.search(r'Total steps:\s+(\d+)', last_block)
                if total_match:
                    self.timing_data['total_steps'] = int(total_match.group(1))
                
                # Parse performance data (last N steps)
                last_timing_match = re.search(r'Average timings for last\s+(\d+) steps:.*?Per Step\(ms\)\s*=\s*([\d.]+).*?ns/day\s*=\s*([\d.]+)', last_block, re.DOTALL)
                if last_timing_match:
                    self.timing_data['last_steps_count'] = int(last_timing_match.group(1))
                    self.timing_data['ms_per_step'] = float(last_timing_match.group(2))
                    self.timing_data['ns_per_day'] = float(last_timing_match.group(3))
                
                # Parse overall performance data
                all_timing_match = re.search(r'Average timings for all steps:.*?Per Step\(ms\)\s*=\s*([\d.]+).*?ns/day\s*=\s*([\d.]+)', last_block, re.DOTALL)
                if all_timing_match:
                    self.timing_data['overall_ms_per_step'] = float(all_timing_match.group(1))
                    self.timing_data['overall_ns_per_day'] = float(all_timing_match.group(2))
                
                # Parse estimated time remaining
                eta_match = re.search(r'Estimated time remaining:\s*([\d.]+)\s*(\w+)', last_block)
                if eta_match:
                    eta_value = float(eta_match.group(1))
                    eta_unit = eta_match.group(2).lower()
                    
                    # Convert to minutes for consistency
                    if 'hour' in eta_unit:
                        eta_minutes = eta_value * 60
                    elif 'second' in eta_unit:
                        eta_minutes = eta_value / 60
                    elif 'day' in eta_unit:
                        eta_minutes = eta_value * 24 * 60
                    else:  # assume minutes
                        eta_minutes = eta_value
                        
                    self.timing_data['eta_minutes'] = eta_minutes
                    self.timing_data['eta_original'] = f"{eta_value} {eta_unit}"
                    
        except Exception as e:
            print(f"Error parsing info file: {e}")

    def get_timing_info(self) -> str:
        """Get formatted timing information."""
        self.parse_info_file()
        
        if not self.timing_data:
            return "⏱️  No timing data available yet..."
            
        info = "\n🚀 Performance Monitor\n"
        info += "=" * 50 + "\n"
        
        # Progress
        if 'completion_percent' in self.timing_data:
            percent = self.timing_data['completion_percent']
            completed = self.timing_data.get('completed_steps', 0)
            total = self.timing_data.get('total_steps', 0)
            
            # Progress bar
            bar_width = 30
            filled = int(bar_width * percent / 100)
            bar = "█" * filled + "▒" * (bar_width - filled)
            
            info += f"📊 Progress: {percent:5.1f}% [{bar}]\n"
            info += f"📈 Steps: {completed:,} / {total:,}\n\n"
        
        # Performance metrics
        if 'ns_per_day' in self.timing_data:
            info += f"⚡ Performance (recent {self.timing_data.get('last_steps_count', 0):,} steps):\n"
            info += f"   🕐 {self.timing_data['ms_per_step']:.2f} ms/step\n"
            info += f"   🏃 {self.timing_data['ns_per_day']:.1f} ns/day\n\n"
        
        if 'overall_ns_per_day' in self.timing_data:
            info += f"📊 Overall Performance:\n"
            info += f"   🕐 {self.timing_data['overall_ms_per_step']:.2f} ms/step\n"
            info += f"   🏃 {self.timing_data['overall_ns_per_day']:.1f} ns/day\n\n"
        
        # ETA
        if 'eta_minutes' in self.timing_data:
            eta_min = self.timing_data['eta_minutes']
            if eta_min > 60:
                hours = int(eta_min // 60)
                minutes = int(eta_min % 60)
                eta_str = f"{hours}h {minutes}m"
            else:
                eta_str = f"{eta_min:.1f} min"
                
            info += f"⏰ Estimated Time Remaining: {eta_str}\n"
            info += f"   (Original: {self.timing_data['eta_original']})\n"
        
        return info

    def ascii_plot(self, data_key: str, title: str, units: str = "", height: int = 15, width: int = 60) -> str:
        """Generate ASCII plot of data."""
        if data_key not in self.data or len(self.data[data_key]) == 0:
            return f"No {title} data available yet..."
            
        values = list(self.data[data_key])
        if len(values) < 2:
            return f"Insufficient {title} data (need at least 2 points)..."
            
        # Normalize data to plot height
        min_val = min(values)
        max_val = max(values)
        
        if max_val == min_val:
            return f"{title}: Constant value {max_val:.2f} {units}"
            
        range_val = max_val - min_val
        
        # Create plot grid
        plot = []
        for i in range(height):
            plot.append([' '] * width)
            
        # Plot data points
        for i, val in enumerate(values):
            x = int((i / (len(values) - 1)) * (width - 1))
            y = height - 1 - int(((val - min_val) / range_val) * (height - 1))
            
            if 0 <= x < width and 0 <= y < height:
                plot[y][x] = '*' if i == len(values) - 1 else '●'
                
        # Add axis labels
        result = f"\n📊 {title} {units}\n"
        result += "=" * width + "\n"

        # Max label line - pad to match plot width
        max_label = f"Max: {max_val:8.2f}"
        result += max_label + " " * (width - len(max_label)) + " ┤\n"

        for i, row in enumerate(plot):
            if i == height // 2:
                result += f"{''.join(row)} ┤ {(min_val + max_val)/2:8.2f}\n"
            else:
                result += f"{''.join(row)} ┤\n"

        # Min label line - pad to match plot width
        min_label = f"Min: {min_val:8.2f}"
        result += min_label + " " * (width - len(min_label)) + " ┤\n"
        result += "─" * width + "┘\n"
        
        # Step/time axis for minimization vs MD
        if 'step' in self.data and len(self.data['step']) > 0:
            start_step = self.data['step'][0] if len(self.data['step']) > 0 else 0
            end_step = self.data['step'][-1] if len(self.data['step']) > 0 else 0
            
            # Check if this is minimization (no time data) or MD (has time data)
            if 'time' in self.data and len(self.data['time']) > 0 and self.data['time'][-1] != self.data['step'][-1]:
                # MD simulation - show time
                start_time = self.data['time'][0]
                end_time = self.data['time'][-1]
                result += f"{start_time:8.1f}" + " " * (width - 16) + f"{end_time:8.1f} (ps)\n"
            else:
                # Minimization - show steps
                result += f"{start_step:8.0f}" + " " * (width - 16) + f"{end_step:8.0f} (steps)\n"
            
        return result

    def infer_simulation_stage(self):
        """
        Infer simulation stage from mdin parameters and data trends.
        Returns tuple: (inferred_stage, confidence, evidence_list)
        """
        evidence = []

        # Check for minimization (highest confidence)
        if self.mdin_params.get('imin') == 1:
            return 'minimization', 'high', ['imin=1 (minimization mode)']

        # Determine ensemble type
        ntp = self.mdin_params.get('ntp', 0)
        ntb = self.mdin_params.get('ntb', 1)

        if ntp > 0 and ntb == 2:
            ensemble = 'NPT'
            evidence.append(f'ntp={ntp}, ntb={ntb} (NPT ensemble)')
        else:
            ensemble = 'NVT'
            evidence.append(f'ntp={ntp}, ntb={ntb} (NVT ensemble)')

        # Check for heating (temperature ramping)
        tempi = self.mdin_params.get('tempi', None)
        temp0 = self.mdin_params.get('temp0', None)

        if tempi is not None and temp0 is not None:
            if tempi < temp0:
                evidence.append(f'tempi={tempi}K, temp0={temp0}K (temperature ramping)')
                return f'{ensemble.lower()}_heating', 'high', evidence
            else:
                evidence.append(f'tempi={tempi}K, temp0={temp0}K (at target temperature)')

        # The same assessment applies to both equilibration and production
        # for a given ensemble, so return the unified stage name directly.
        return f'{ensemble.lower()}_equilibration', 'medium', evidence

    def get_natom(self):
        """Return atom count from system info."""
        return self.system_info.get('natom', 0)


class CheckpointManager:
    """Manage workflow checkpoints for resume functionality."""
    
    def __init__(self, checkpoint_file: str = "workflow_checkpoint.json"):
        self.checkpoint_file = checkpoint_file
        
    def save_checkpoint(self, step_index: int, step_name: str, completed_steps: List[str], 
                       prmtop: str, inpcrd: str, current_rst: str, mode: int):
        """Save workflow state to checkpoint file."""
        checkpoint_data = {
            'timestamp': datetime.now().isoformat(),
            'step_index': step_index,
            'step_name': step_name,
            'completed_steps': completed_steps,
            'prmtop': prmtop,
            'inpcrd': inpcrd,
            'current_rst': current_rst,
            'mode': mode
        }
        
        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
            print(f"💾 Checkpoint saved: {step_name}")
        except Exception as e:
            print(f"❌ Failed to save checkpoint: {e}")
            
    def load_checkpoint(self) -> Optional[Dict]:
        """Load workflow state from checkpoint file."""
        if not os.path.exists(self.checkpoint_file):
            return None
            
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint_data = json.load(f)
            return checkpoint_data
        except Exception as e:
            print(f"❌ Failed to load checkpoint: {e}")
            return None
            
    def checkpoint_exists(self) -> bool:
        """Check if checkpoint file exists."""
        return os.path.exists(self.checkpoint_file)
        
    def delete_checkpoint(self):
        """Delete checkpoint file."""
        try:
            if os.path.exists(self.checkpoint_file):
                os.remove(self.checkpoint_file)
                print("🗑️  Checkpoint file deleted")
        except Exception as e:
            print(f"❌ Failed to delete checkpoint: {e}")


class WorkflowExecutor:
    """Execute AMBER workflow steps with monitoring."""
    
    def __init__(self):
        self.current_process = None
        self.workflow_steps = [
            "minimization",
            "heating", 
            "equilibration_1",
            "equilibration_2",
            "production"
        ]
        
    def build_amber_command(self, step: str, use_mpi: bool = True, num_cores: int = 32, 
                           gpu_device: int = 0) -> str:
        """Build AMBER command for the specified step."""
        # Determine input coordinate file
        if step == "minimization":
            coord_file = "system.inpcrd"
        else:
            prev_step_index = self.workflow_steps.index(step) - 1
            prev_step = self.workflow_steps[prev_step_index]
            coord_file = f"{prev_step}.rst7"
        
        # Determine engine and build command
        if step in ["minimization", "equilibration_2"]:
            # CPU with MPI
            if use_mpi:
                cmd_parts = [
                    f"mpirun -np {num_cores} pmemd.MPI",
                    "-O",
                    f"-i {step}.in" if step != "minimization" else "-i min.in",
                    f"-o {step}.out",
                    f"-x {step}.nc",
                    f"-p system.prmtop",
                    f"-c {coord_file}",
                    f"-inf {step}.info",
                    f"-r {step}.rst7"
                ]
            else:
                cmd_parts = [
                    "pmemd",
                    "-O",
                    f"-i {step}.in" if step != "minimization" else "-i min.in",
                    f"-o {step}.out",
                    f"-x {step}.nc",
                    f"-p system.prmtop", 
                    f"-c {coord_file}",
                    f"-inf {step}.info",
                    f"-r {step}.rst7"
                ]
        else:
            # GPU
            cmd_parts = [
                f"export CUDA_VISIBLE_DEVICES={gpu_device} &&",
                "pmemd.cuda",
                "-O",
                f"-i {step}.in" if step != "heating" else "-i heat.in",
                f"-o {step}.out",
                f"-x {step}.nc",
                f"-p system.prmtop",
                f"-c {coord_file}",
                f"-inf {step}.info", 
                f"-r {step}.rst7"
            ]
        
        return " ".join(cmd_parts)
    
    def start_simulation_background(self, command: str, step: str, output_file: str) -> bool:
        """Start simulation in background."""
        try:
            # Handle compound commands (with export)
            if "&&" in command:
                # Split export and main command
                parts = command.split("&&")
                env_cmd = parts[0].strip()
                main_cmd = parts[1].strip()
                
                # Set environment variable
                if "export" in env_cmd:
                    env_parts = env_cmd.replace("export", "").strip().split("=")
                    if len(env_parts) == 2:
                        os.environ[env_parts[0].strip()] = env_parts[1].strip()
                
                cmd_list = main_cmd.split()
            else:
                cmd_list = command.split()
            
            # Start the process
            self.current_process = subprocess.Popen(
                cmd_list,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to start simulation: {e}")
            return False
    
    def stop_simulation(self) -> bool:
        """Stop current simulation."""
        if self.current_process and self.current_process.poll() is None:
            print("⏹️  Stopping simulation...")
            self.current_process.terminate()
            time.sleep(2)
            if self.current_process.poll() is None:
                self.current_process.kill()
            return True
        else:
            print("❌ No simulation running")
            return False
    
    def is_simulation_running(self) -> bool:
        """Check if simulation is currently running."""
        return self.current_process and self.current_process.poll() is None
    
    def get_simulation_status(self) -> Tuple[bool, Optional[int]]:
        """Get simulation status and return code."""
        if not self.current_process:
            return False, None
        
        return_code = self.current_process.poll()
        if return_code is None:
            return True, None  # Still running
        else:
            return False, return_code  # Finished


class AMBERWorkflowCore:
    """Core workflow functionality from the original amber_workflow.py"""
    
    def __init__(self, processor=None):
        self.processor = processor
        self.workflow_steps = [
            "minimization",
            "heating",
            "equilibration_1",
            "equilibration_2",
            "production"
        ]

        self.checkpoint_manager = CheckpointManager()
        self.workflow_executor = WorkflowExecutor()
        self.current_monitor = None
        self.historical_monitors = {}
        
        # Settings
        self.use_mpi = True
        self.num_cores = 32
        self.gpu_device = 0
        
        # State
        self.completed_steps = []
        self.current_step = None

    def display_banner(self):
        """Display welcome banner."""
        print("=" * 70)
        print("    AMBER Molecular Dynamics Workflow Manager v2.0")
        print("    🚀 Now with Real-time Monitoring & Checkpoints!")
        print("=" * 70)
        print("Features:")
        print("  • Interactive workflow with parameter explanations")
        print("  • Real-time ASCII plots of energy, temperature, pressure")
        print("  • Performance monitoring with ETA calculations")
        print("  • Checkpoint/resume functionality")
        print("  • Historical plot analysis")
        print("  • Background simulation monitoring")
        print("  • Step-by-step or continuous execution")
        print("=" * 70)

    def check_for_resume(self) -> bool:
        """Check if user wants to resume from checkpoint."""
        if not self.checkpoint_manager.checkpoint_exists():
            return False
            
        checkpoint = self.checkpoint_manager.load_checkpoint()
        if not checkpoint:
            return False
            
        print(f"\n📄 Found checkpoint from {checkpoint['timestamp']}")
        print(f"   Last completed: {checkpoint['step_name']}")
        print(f"   Completed steps: {', '.join(checkpoint['completed_steps'])}")
        
        resume = prompt_with_context(
            self.processor, "❓ Resume from checkpoint? (y/N)",
            default="N",
            module="AMBER Workflow",
            description="Resume workflow from saved checkpoint",
        ).strip().lower()
        return resume.startswith('y')

    def list_files(self, extensions: List[str]) -> List[str]:
        """List files with given extensions in current directory."""
        files = []
        for ext in extensions:
            files.extend(glob.glob(f"*.{ext}"))
        return sorted(files)

    def select_file(self, prompt: str, extensions: List[str], required: bool = False) -> Optional[str]:
        """Interactive file selection with proper validation."""
        files = self.list_files(extensions)
        if not files:
            print(f"❌ No files found with extensions: {', '.join(extensions)}")
            if required:
                manual = prompt_with_context(
                    self.processor, "Enter filename manually (required)",
                    module="AMBER Workflow",
                    description="Manual filename entry (required input file)",
                ).strip()
                while not manual or not os.path.exists(manual):
                    if not manual:
                        print("❌ This file is required for the simulation!")
                    else:
                        print(f"❌ File '{manual}' not found!")
                    manual = prompt_with_context(
                        self.processor, "Enter valid filename",
                        module="AMBER Workflow",
                        description="Re-enter valid filename (required input file)",
                    ).strip()
                return manual
            return None
            
        print(f"\n{prompt}")
        for i, file in enumerate(files, 1):
            print(f"  {i}. {file}")
        
        if not required:
            print(f"  {len(files) + 1}. Skip")
        print(f"  {len(files) + (2 if not required else 1)}. Manual entry")
        
        max_choice = len(files) + (2 if not required else 1)
        
        while True:
            try:
                file_options_map = {str(i): f for i, f in enumerate(files, 1)}
                if not required:
                    file_options_map[str(len(files) + 1)] = "Skip"
                file_options_map[str(max_choice)] = "Manual entry"
                choice = prompt_with_context(
                    self.processor, f"Select file (1-{max_choice})",
                    module="AMBER Workflow",
                    description="Select input file for workflow step",
                    options_map=file_options_map,
                ).strip()
                
                # Handle skip option (only for non-required files)
                if not required and choice == str(len(files) + 1):
                    return None
                
                # Handle manual entry
                manual_option = len(files) + (2 if not required else 1)
                if choice == str(manual_option):
                    manual = prompt_with_context(
                        self.processor, "Enter filename manually",
                        module="AMBER Workflow",
                        description="Manual filename entry",
                    ).strip()
                    if required:
                        while not manual or not os.path.exists(manual):
                            if not manual:
                                print("❌ This file is required for the simulation!")
                            else:
                                print(f"❌ File '{manual}' not found!")
                            manual = prompt_with_context(
                                self.processor, "Enter valid filename",
                                module="AMBER Workflow",
                                description="Re-enter valid filename",
                            ).strip()
                    return manual if manual else None
                    
                # Handle numbered selection
                idx = int(choice) - 1
                if 0 <= idx < len(files):
                    return files[idx]
                print("Invalid selection. Try again.")
            except (ValueError, KeyboardInterrupt):
                print("Invalid input. Try again.")

    def run_simulation_background(self, command: str, step: str, output_file: str):
        """Run simulation in background with monitoring."""
        return self.workflow_executor.start_simulation_background(command, step, output_file)

    def interactive_monitor(self, step: str):
        """Interactive monitoring while simulation runs."""
        # Define command shortcuts
        shortcuts = {
            'pe': 'plot energy',
            'pt': 'plot temp', 
            'pp': 'plot pressure',
            'pk': 'plot kinetic',
            'ppo': 'plot potential',
            'pr': 'plot rms',
            's': 'status',
            't': 'timing',
            'h': 'help',
            'c': 'continue',
            'st': 'stop'
        }
        
        while True:
            is_running, return_code = self.workflow_executor.get_simulation_status()
            
            if not is_running:
                # Process finished
                if return_code == 0:
                    print(f"\n✅ {step} completed successfully!")
                    self.post_simulation_analysis(step)
                    return True
                else:
                    print(f"\n❌ {step} failed!")
                    return False
            
            try:
                # Get input with retry for empty strings
                user_input = ""
                while not user_input:
                    user_input = prompt_with_context(
                        self.processor, f"[{step}] Monitor>",
                        module="AMBER Workflow - Monitor",
                        description="Interactive monitor command",
                    ).strip()
                    if not user_input:
                        print("⚠️  Empty input detected, please try again.")
                
                # Expand shortcuts
                original_input = user_input
                user_input_lower = user_input.lower()
                if user_input_lower in shortcuts:
                    user_input = shortcuts[user_input_lower]
                    print(f"💡 Shortcut '{original_input}' → '{user_input}'")
                else:
                    user_input = user_input_lower
                
                if user_input == 'status':
                    if self.workflow_executor.is_simulation_running():
                        print("🟢 Simulation running...")
                        self.current_monitor.parse_amber_output()
                        if 'step' in self.current_monitor.data and len(self.current_monitor.data['step']) > 0:
                            last_step = self.current_monitor.data['step'][-1]
                            last_time = self.current_monitor.data['time'][-1] if 'time' in self.current_monitor.data else 0
                            print(f"   Last step: {last_step}, Time: {last_time:.1f} ps")
                    else:
                        print("🔴 Simulation not running")
                        
                elif user_input in ['timing', 'performance', 'eta']:
                    print(self.current_monitor.get_timing_info())
                        
                elif user_input.startswith('plot '):
                    metric = user_input.split(' ', 1)[1]
                    self.current_monitor.parse_amber_output()
                    
                    if metric in ['temp', 'temperature']:
                        print(self.current_monitor.ascii_plot('temperature', 'Temperature', '(K)'))
                    elif metric in ['energy', 'total']:
                        print(self.current_monitor.ascii_plot('total_energy', 'Total Energy', '(kcal/mol)'))
                    elif metric == 'pressure':
                        print(self.current_monitor.ascii_plot('pressure', 'Pressure', '(bar)'))
                    elif metric == 'kinetic':
                        print(self.current_monitor.ascii_plot('kinetic_energy', 'Kinetic Energy', '(kcal/mol)'))
                    elif metric == 'potential':
                        print(self.current_monitor.ascii_plot('potential_energy', 'Potential Energy', '(kcal/mol)'))
                    elif metric in ['rms', 'gradient']:
                        print(self.current_monitor.ascii_plot('rms_gradient', 'RMS Gradient', '(kcal/mol/Å)'))
                    else:
                        print("❌ Unknown metric. Available: temp, energy, pressure, kinetic, potential, rms")
                        
                elif user_input == 'stop':
                    return self.workflow_executor.stop_simulation()
                        
                elif user_input == 'continue':
                    if self.workflow_executor.is_simulation_running():
                        print("⚠️  Simulation still running! Wait for completion or type 'stop'")
                    else:
                        break
                        
                elif user_input in ['help', '?']:
                    print("Available commands:")
                    print("  status (s) - Check simulation status")
                    print("  plot <metric> - Show real-time plot")
                    print("    Shortcuts: pe=energy, pt=temp, pp=pressure, pk=kinetic, ppo=potential, pr=rms")
                    print("  timing/performance/eta (t) - Show performance & ETA")
                    print("  stop (st) - Terminate simulation")
                    print("  continue (c) - Proceed (when simulation done)")
                    print("  help (h) - Show this help")
                    print("💡 Use ↑↓ arrow keys for command history")
                    
                else:
                    print("❓ Unknown command. Type 'help' for available commands.")
                    
            except KeyboardInterrupt:
                print("\n⏹️  Interrupted by user")
                return self.workflow_executor.stop_simulation()
                
        return True

    def post_simulation_analysis(self, step: str):
        """Offer plot analysis after simulation completes."""
        output_file = f"{step}.out"
        if not os.path.exists(output_file):
            return
            
        print(f"\n🎉 {step.replace('_', ' ').title()} completed!")
        view_plots = prompt_with_context(
            self.processor, "📊 View completion plots? (y/N)",
            default="N",
            module="AMBER Workflow",
            description="View completion plots for finished step",
        ).strip().lower()
        
        if view_plots.startswith('y'):
            self.analyze_single_step(step)

    def analyze_single_step(self, step: str):
        """Analyze and show plots for a single completed step."""
        output_file = f"{step}.out"
        if not os.path.exists(output_file):
            print(f"❌ Output file {output_file} not found!")
            return
            
        print(f"\n📊 Analyzing {step.replace('_', ' ').title()}")
        print("-" * 40)
        
        # Create or get cached monitor
        if step not in self.historical_monitors:
            print("📖 Loading simulation data...")
            self.historical_monitors[step] = AMBERMonitor.create_historical_monitor(output_file, max_points=200)
        
        monitor = self.historical_monitors[step]
        
        if not any(len(monitor.data[key]) > 0 for key in monitor.data):
            print(f"❌ No data found in {output_file}")
            return
            
        # Show available metrics
        available_metrics = []
        if len(monitor.data.get('total_energy', [])) > 0:
            available_metrics.append('energy')
        if len(monitor.data.get('temperature', [])) > 0:
            available_metrics.append('temp')
        if len(monitor.data.get('pressure', [])) > 0:
            available_metrics.append('pressure')
        if len(monitor.data.get('rms_gradient', [])) > 0:
            available_metrics.append('rms')
        if len(monitor.data.get('kinetic_energy', [])) > 0:
            available_metrics.append('kinetic')
        if len(monitor.data.get('potential_energy', [])) > 0:
            available_metrics.append('potential')
            
        if not available_metrics:
            print("❌ No plottable data found!")
            return
            
        print(f"📈 Available metrics: {', '.join(available_metrics)}")
        print("💡 Commands: plot <metric>, summary, or 'done' to finish")
        
        while True:
            try:
                cmd = prompt_with_context(
                    self.processor, f"[{step}] Analyze>",
                    module="AMBER Workflow - Analyze",
                    description="Interactive analysis command",
                ).strip().lower()
                
                if cmd == 'done' or cmd == 'exit':
                    break
                elif cmd == 'summary':
                    self.show_step_summary(monitor, step)
                elif cmd.startswith('plot '):
                    metric = cmd.split(' ', 1)[1]
                    if metric in ['energy', 'total']:
                        print(monitor.ascii_plot('total_energy', 'Total Energy', '(kcal/mol)'))
                    elif metric in ['temp', 'temperature']:
                        print(monitor.ascii_plot('temperature', 'Temperature', '(K)'))
                    elif metric == 'pressure':
                        print(monitor.ascii_plot('pressure', 'Pressure', '(bar)'))
                    elif metric in ['rms', 'gradient']:
                        print(monitor.ascii_plot('rms_gradient', 'RMS Gradient', '(kcal/mol/Å)'))
                    elif metric == 'kinetic':
                        print(monitor.ascii_plot('kinetic_energy', 'Kinetic Energy', '(kcal/mol)'))
                    elif metric == 'potential':
                        print(monitor.ascii_plot('potential_energy', 'Potential Energy', '(kcal/mol)'))
                    else:
                        print(f"❌ Unknown metric. Available: {', '.join(available_metrics)}")
                else:
                    print("💡 Commands: 'plot <metric>', 'summary', or 'done'")
                    
            except KeyboardInterrupt:
                break

    def show_step_summary(self, monitor: AMBERMonitor, step: str):
        """Show a summary of key metrics for a step."""
        summary = ""
        
        # Energy summary
        if len(monitor.data.get('total_energy', [])) > 0:
            energies = list(monitor.data['total_energy'])
            summary += f"⚡ Energy: {energies[0]:.1f} → {energies[-1]:.1f} kcal/mol"
            if len(energies) > 1:
                change = energies[-1] - energies[0]
                summary += f" (Δ{change:+.1f})\n"
            else:
                summary += "\n"
        
        # Temperature summary
        if len(monitor.data.get('temperature', [])) > 0:
            temps = list(monitor.data['temperature'])
            summary += f"🌡️  Temperature: {temps[0]:.1f} → {temps[-1]:.1f} K"
            if len(temps) > 1:
                avg_temp = sum(temps) / len(temps)
                summary += f" (avg: {avg_temp:.1f} K)\n"
            else:
                summary += "\n"
        
        # Pressure summary (if available)
        if len(monitor.data.get('pressure', [])) > 0:
            pressures = list(monitor.data['pressure'])
            avg_pressure = sum(pressures) / len(pressures)
            summary += f"📊 Pressure: avg {avg_pressure:.1f} bar\n"
        
        # RMS gradient (for minimization)
        if len(monitor.data.get('rms_gradient', [])) > 0:
            rms_vals = list(monitor.data['rms_gradient'])
            summary += f"📉 RMS Gradient: {rms_vals[0]:.3f} → {rms_vals[-1]:.3f} kcal/mol/Å\n"
        
        # Data points
        if len(monitor.data.get('step', [])) > 0:
            steps = list(monitor.data['step'])
            summary += f"📈 Data points: {len(steps)} (steps {steps[0]}-{steps[-1]})\n"
        
        print(summary if summary else "❌ No data available for analysis")