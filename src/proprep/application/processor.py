"""
Processor interface for MPSA (Modular PDB Structure Analysis) application.

"""

from abc import ABC, abstractmethod


class Processor(ABC):
    """
    Abstract base class for all processors.
    Used for type checking and interface definition.
    """

    @abstractmethod
    def run_main_menu(self):
        pass

    @abstractmethod
    def cleanup(self):
        pass

    @abstractmethod
    def get_module_instance(self, name: str):
        pass

    @abstractmethod
    def run_workflow_step(self, step_name: str) -> bool:
        pass

    @abstractmethod
    def display_workspace_status(self, detailed=False):
        pass
