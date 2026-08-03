"""
Utility modules for MPSA (Molecular Preparation and Simulation Assistant)

This package provides core utilities used throughout the MPSA package:
- Workspace management for data sharing between modules
- Module registry for plugin-like functionality
- Debugging utilities
- Session recording and replay for reproducibility
"""
from proprep.utils.module_registry import (ModuleRegistry, ProcessingModule,
                                        register_module, registry)
from proprep.utils.mpsa_workspace import MPSAWorkspace
from proprep.utils.workspace import Workspace
from proprep.utils.session_recorder import (
    SessionRecorder,
    SessionReplayer,
    InterceptedPrompt,
    SessionManager,
    integrate_session_manager
)
from .debug_utils import debug_workspace

__all__ = [
    "Workspace",
    "MPSAWorkspace",
    "ProcessingModule",
    "ModuleRegistry",
    "registry",
    "register_module",
    "debug_workspace",
    "SessionRecorder",
    "SessionReplayer",
    "InterceptedPrompt",
    "SessionManager",
    "integrate_session_manager",
]
