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


def unwrap_serialized(value: Any) -> Any:
    """Strip workspace serialization wrappers from a value, recursively.

    Checklist state is written to JSON, and anything that is not a plain scalar
    is stored as ``{"__type__": <name>, "value": <payload>}``. Reading it back
    gives that wrapper, not the payload, so code expecting the original dict
    shape sees ``__type__``/``value`` where it expected its own keys — a
    resumed run handed a ``RedoxSite`` wrapper to ``dict_to_redox_site`` and it
    raised ``KeyError: 'site_id'``.

    Two dict shapes therefore exist for the same object: the flat one written by
    the exporters, and this wrapped one. This makes them interchangeable.

    Enums serialize as a dict carrying ``_value_``; that is unwrapped to the
    underlying value so ``CenterType(...)`` and friends can rebuild from it.
    """
    if isinstance(value, dict):
        if "__type__" in value and "value" in value:
            type_name = value["__type__"]
            inner = value["value"]
            if type_name == "tuple":
                return tuple(unwrap_serialized(v) for v in inner)
            if type_name == "str":
                return str(inner)
            if type_name == "Path":
                return str(inner)
            if isinstance(inner, dict) and "_value_" in inner:
                return inner["_value_"]
            return unwrap_serialized(inner)
        return {k: unwrap_serialized(v) for k, v in value.items()}
    if isinstance(value, list):
        return [unwrap_serialized(v) for v in value]
    return value
