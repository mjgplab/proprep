"""
Debug Workspace Wrapper for MPSA

This module provides a debugging wrapper around the Workspace class
to track all access and modifications to workspace data.
"""

import traceback
from datetime import datetime
from pathlib import Path
import os

class DebugWorkspace:
    """Wrapper around workspace to debug all access"""
    
    def __init__(self, workspace, log_file=None):
        self._workspace = workspace
        self._access_log = []
        self._monitored_keys = ["filtered_structure", "repaired_structure", "structure"]
        self._log_file = log_file
        
        # Write header to log file
        if self._log_file:
            with open(self._log_file, 'w') as f:
                f.write(f"WORKSPACE DEBUG LOG - Started at {datetime.now()}\n")
                f.write("="*80 + "\n\n")
    
    def _log_to_file(self, message):
        """Write message to log file if configured"""
        if self._log_file:
            with open(self._log_file, 'a') as f:
                f.write(message + "\n")
    
    def get(self, key, default=None):
        """Wrapped get method with logging"""
        value = self._workspace.get(key, default)
        
        if key in self._monitored_keys:
            # Get the calling location
            stack = traceback.extract_stack()
            caller = "Unknown"
            caller_detail = ""
            if len(stack) >= 2:
                frame = stack[-2]
                caller = f"{frame.filename.split('/')[-1]}:{frame.lineno} in {frame.name}"
                # Get more context - look up the stack for module methods
                for i in range(len(stack)-3, max(0, len(stack)-6), -1):
                    frame = stack[i]
                    if "handle_menu_option" in frame.name or "transform" in frame.name or "generate" in frame.name:
                        caller_detail = f" via {frame.name}"
                        break
            
            # Log the access
            log_entry = {
                "time": datetime.now().strftime("%H:%M:%S.%f")[:-3],
                "operation": "GET",
                "key": key,
                "value_type": type(value).__name__,
                "is_path": isinstance(value, (str, Path)),
                "caller": caller + caller_detail
            }
            
            self._access_log.append(log_entry)
            
            # Print immediate debug info for monitored keys
            output = []
            output.append(f"\n[DEBUG GET] {log_entry['time']} - {key}")
            output.append(f"  Type: {log_entry['value_type']}")
            output.append(f"  Is Path: {log_entry['is_path']}")
            if isinstance(value, (str, Path)):
                output.append(f"  Path: {value}")
                if os.path.exists(str(value)):
                    output.append(f"  File exists: YES")
                else:
                    output.append(f"  File exists: NO")
            elif hasattr(value, "get_atoms"):
                atom_count = sum(1 for _ in value.get_atoms())
                output.append(f"  BioPython structure: {atom_count} atoms")
            output.append(f"  Caller: {caller + caller_detail}")
            
            # Print to console
            for line in output:
                print(line)
            
            # Log to file
            self._log_to_file("\n".join(output))
        
        return value
    
    def set(self, key, value):
        """Wrapped set method with logging"""
        old_value = self._workspace.get(key)
        self._workspace.set(key, value)
        
        if key in self._monitored_keys:
            # Get the calling location
            stack = traceback.extract_stack()
            caller = "Unknown"
            caller_detail = ""
            if len(stack) >= 2:
                frame = stack[-2]
                caller = f"{frame.filename.split('/')[-1]}:{frame.lineno} in {frame.name}"
                # Get more context
                for i in range(len(stack)-3, max(0, len(stack)-6), -1):
                    frame = stack[i]
                    if "handle_menu_option" in frame.name or "transform" in frame.name or "save" in frame.name:
                        caller_detail = f" via {frame.name}"
                        break
            
            # Log the change
            log_entry = {
                "time": datetime.now().strftime("%H:%M:%S.%f")[:-3],
                "operation": "SET",
                "key": key,
                "old_type": type(old_value).__name__ if old_value is not None else "None",
                "new_type": type(value).__name__,
                "old_is_path": isinstance(old_value, (str, Path)),
                "new_is_path": isinstance(value, (str, Path)),
                "caller": caller + caller_detail
            }
            
            self._access_log.append(log_entry)
            
            # Print immediate debug info for monitored keys
            output = []
            output.append(f"\n[DEBUG SET] {log_entry['time']} - {key}")
            output.append(f"  Old Type: {log_entry['old_type']} (Path: {log_entry['old_is_path']})")
            output.append(f"  New Type: {log_entry['new_type']} (Path: {log_entry['new_is_path']})")
            if isinstance(value, (str, Path)):
                output.append(f"  New Path: {value}")
            elif hasattr(value, "get_atoms") and hasattr(old_value, "get_atoms"):
                old_atoms = sum(1 for _ in old_value.get_atoms()) if old_value else 0
                new_atoms = sum(1 for _ in value.get_atoms())
                output.append(f"  Atom count: {old_atoms} -> {new_atoms}")
            output.append(f"  Caller: {caller + caller_detail}")
            
            # Alert on suspicious changes
            if log_entry['old_is_path'] and not log_entry['new_is_path']:
                output.append(f"  ⚠️  WARNING: Path replaced with non-path object!")
                output.append(f"  ⚠️  This is likely where the issue occurs!")
            
            # Print to console
            for line in output:
                print(line)
            
            # Log to file
            self._log_to_file("\n".join(output))
    
    def has(self, key):
        """Wrapped has method"""
        return self._workspace.has(key)
    
    def clear(self):
        """Wrapped clear method"""
        print("\n[DEBUG CLEAR] Workspace cleared!")
        self._log_to_file("\n[DEBUG CLEAR] Workspace cleared!")
        return self._workspace.clear()
    
    def keys(self):
        """Wrapped keys method"""
        return self._workspace.keys()
    
    def items(self):
        """Wrapped items method"""
        return self._workspace.items()
    
    def print_access_log(self):
        """Print the full access log"""
        output = []
        output.append("\n" + "="*80)
        output.append("WORKSPACE ACCESS LOG SUMMARY")
        output.append("="*80)
        
        for entry in self._access_log:
            if entry['operation'] == 'GET':
                output.append(f"{entry['time']} GET {entry['key']} -> {entry['value_type']} (Path: {entry['is_path']})")
            else:
                output.append(f"{entry['time']} SET {entry['key']}: {entry['old_type']} -> {entry['new_type']}")
                if entry['old_is_path'] and not entry['new_is_path']:
                    output.append(f"           ⚠️  Path replaced with non-path!")
            output.append(f"           Called from: {entry['caller']}")
        
        output.append("="*80)
        
        # Print to console
        for line in output:
            print(line)
        
        # Log to file
        if self._log_file:
            self._log_to_file("\n".join(output))
    
    def __getattr__(self, name):
        """Forward any other attributes to the wrapped workspace"""
        return getattr(self._workspace, name)
    
    # Support dictionary-like access
    def __getitem__(self, key):
        """Support workspace['key'] syntax"""
        return self.get(key)
    
    def __setitem__(self, key, value):
        """Support workspace['key'] = value syntax"""
        self.set(key, value)
    
    def __contains__(self, key):
        """Support 'key in workspace' syntax"""
        return self.has(key)
