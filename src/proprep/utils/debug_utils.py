# debug_utils.py
"""Debug utilities for PDB processor modules"""

from typing import Any, Dict


def debug_workspace(module_name: str, action: str, workspace: Dict[str, Any], console: Any):
    """Debug utility to print workspace state"""
    if workspace.get("debug", False):
        console.print(f"[yellow]DEBUG: {module_name} {action} workspace:[/yellow]")
        for key, value in workspace.items():
            if key != "debug":
                value_type = type(value).__name__
                if hasattr(value, "__len__"):
                    value_info = f"{value_type}({len(value)})"
                else:
                    value_info = value_type
                console.print(f"[yellow]  {key}: {value_info}[/yellow]")
