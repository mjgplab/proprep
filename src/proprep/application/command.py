"""
Abstract Command Module

"""

from abc import ABC, abstractmethod
from typing import Any, Optional
from datetime import datetime

from .processor import Processor


class Command(ABC):
    """
    Enhanced base command interface with history and metadata support.

    Extends the original Command class from command.py with additional
    functionality for breadcrumb tracking and command metadata.
    """

    def __init__(self, processor: Processor, description: str = ""):
        self.processor = processor
        self.description = description
        self.timestamp = datetime.now()
        self.executed = False
        self.result: Any = None
        self.error: Optional[Exception] = None

    @abstractmethod
    def execute(self) -> Any:
        """Execute the command and return result."""
        pass

    def can_execute(self) -> bool:
        """Check if command can be executed in current state."""
        return True

    def get_breadcrumb(self) -> str:
        """Get breadcrumb text for this command."""
        return self.description or self.__class__.__name__

    def __str__(self) -> str:
        """String representation of the command."""
        status = "✓" if self.executed else "○"
        error_info = f" (Error: {self.error})" if self.error else ""
        return f"{status} {self.description}{error_info}"
