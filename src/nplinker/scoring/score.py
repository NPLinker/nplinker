from __future__ import annotations
from dataclasses import dataclass
from dataclasses import fields
from .scoring_method import ScoringMethod


@dataclass
class Score:
    """A data class to represent score data.

    Attributes:
        name: the name of the scoring method. See `ScoringMethod` for valid values.
        value: the score value.
        parameter: the parameters used for the scoring method.
    """

    name: str
    value: float
    parameter: dict

    def __post_init__(self) -> None:
        """Check if the value of `name` is valid.

        Raises:
            ValueError: if the value of `name` is not valid.
        """
        if ScoringMethod.has_value(self.name) is False:
            raise ValueError(
                f"{self.name} is not a valid value. Valid values are: {[e.value for e in ScoringMethod]}"
            )

    def __getitem__(self, key):
        if key in {field.name for field in fields(self)}:
            return getattr(self, key)
        else:
            raise KeyError(f"{key} not found in {self.__class__.__name__}")

    def __setitem__(self, key, value):
        # validate the value of `name`
        if key == "name" and ScoringMethod.has_value(value) is False:
            raise ValueError(
                f"{value} is not a valid value. Valid values are: {[e.value for e in ScoringMethod]}"
            )

        if key in {field.name for field in fields(self)}:
            setattr(self, key, value)
        else:
            raise KeyError(f"{key} not found in {self.__class__.__name__}")
