from __future__ import annotations
import logging
from abc import ABC
from abc import abstractmethod
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from nplinker.nplinker import NPLinker
    from .link_graph import LinkGraph

logger = logging.getLogger(__name__)


class ScoringBase(ABC):
    """Abstract base class of scoring methods.

    Attributes:
        name: The name of the scoring method.
        npl: The NPLinker object.
    """

    name: str = "ScoringBase"

    def __init__(self, npl: NPLinker):
        """Initialize the scoring method.

        Args:
            npl: The NPLinker object.
        """
        self.npl = npl

    @classmethod
    @abstractmethod
    def setup(cls, npl: NPLinker):
        """Setup class level attributes."""

    @abstractmethod
    def get_links(
        self,
        *objects,
        **parameters,
    ) -> LinkGraph:
        """Get links information for the given objects.

        Args:
            objects: A list of objects to get links for.
            parameters: The parameters used for scoring.

        Returns:
            The LinkGraph object.
        """

    @abstractmethod
    def format_data(self, data) -> str:
        """Format the scoring data to a string."""

    @abstractmethod
    def sort(self, objects, reverse=True) -> list:
        """Sort the given objects based on the scoring data."""
