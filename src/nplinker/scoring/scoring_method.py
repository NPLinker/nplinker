from enum import Enum
from enum import unique


@unique
class ScoringMethod(Enum):
    """Enum class for scoring methods."""

    METCALF = "metcalf"
    ROSETTA = "rosetta"
    NPLCLASS = "nplclass"

    @classmethod
    def has_value(cls, value: str) -> bool:
        """Check if the enum has a value."""
        return any(value == item.value for item in cls)
