"""
AMBER MDIN Configuration Interface - Exceptions

Custom exceptions for MDIN parsing and keyword lookup.
"""


class MdinSyntaxError(Exception):
    """Raised when an MDIN file has invalid syntax."""

    def __init__(
        self,
        message: str,
        line_number: int | None = None,
        line_text: str | None = None,
    ):
        self.line_number = line_number
        self.line_text = line_text
        if line_number is not None:
            full_message = f"Line {line_number}: {message}"
            if line_text:
                full_message += f"\n  {line_text.strip()}"
        else:
            full_message = message
        super().__init__(full_message)


class AmbiguousKeywordError(Exception):
    """Raised when a keyword name exists in multiple sections and section is not specified."""

    def __init__(self, name: str, sections: list[str]):
        self.name = name
        self.sections = sections
        message = (
            f"Keyword '{name}' exists in multiple sections: {', '.join(sections)}. "
            f"Please specify the section explicitly."
        )
        super().__init__(message)
