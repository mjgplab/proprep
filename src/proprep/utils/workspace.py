from abc import ABC, abstractmethod
from typing import Any, Dict, List, Tuple


class Workspace(ABC):
    """This is a dictionary like interface for storing and managing the program state"""

    @abstractmethod
    def get(self, key: str, default: Any = None) -> Any:
        pass

    @abstractmethod
    def set(self, key: str, value: Any) -> None:
        pass

    @abstractmethod
    def has(self, key: str) -> bool:
        pass

    @abstractmethod
    def update(self, updates: Dict[str, Any]) -> None:
        pass

    @abstractmethod
    def dump(self) -> None:
        pass

    @abstractmethod
    def clear(self) -> None:
        pass

    @abstractmethod
    def get_keys(self) -> List[str]:
        pass

    @abstractmethod
    def items(self) -> List[Tuple[str, Any]]:
        """Return a list of key-value pairs."""
        pass

    def __iter__(self):
        """Allow iteration over key-value pairs."""
        return iter(self.items())
