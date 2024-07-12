from __future__ import annotations
from functools import wraps
from typing import Union
from networkx import Graph
from nplinker.genomics import GCF
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum
from .score import Score
from .scoring_method import ScoringMethod


# Type aliases
Entity = Union[GCF, Spectrum, MolecularFamily]  # using Union to ensure python 3.9 compatibility
LINK_DATA = dict[str, Score]
LINK = tuple[Entity, Entity, LINK_DATA]


def validate_u(func):
    """A decorator to validate the type of the u object."""

    @wraps(func)
    def wrapper(self, u: Entity, *args, **kwargs):
        if not isinstance(u, (GCF, Spectrum, MolecularFamily)):
            raise TypeError(f"{u} is not a GCF, Spectrum, or MolecularFamily object.")

        return func(self, u, *args, **kwargs)

    return wrapper


def validate_uv(func):
    """A decorator to validate the types of the u and v objects."""

    @wraps(func)
    def wrapper(
        self,
        u: Entity,
        v: Entity,
        *args,
        **kwargs,
    ):
        if isinstance(u, GCF):
            if not isinstance(v, (Spectrum, MolecularFamily)):
                raise TypeError(f"{v} is not a Spectrum or MolecularFamily object.")
        elif isinstance(u, (Spectrum, MolecularFamily)):
            if not isinstance(v, GCF):
                raise TypeError(f"{v} is not a GCF object.")
        else:
            raise TypeError(f"{u} is not a GCF, Spectrum, or MolecularFamily object.")

        return func(self, u, v, *args, **kwargs)

    return wrapper


class LinkGraph:
    """Class to represent the links between objects in NPLinker.

    This class wraps the `networkx.Graph` class to provide a more user-friendly interface for
    working with the links.

    The links between objects are stored as edges in a graph, while the objects themselves are
    stored as nodes.

    The scoring data for each link (or link data) is stored as the key/value attributes of the edge.
    """

    def __init__(self) -> None:
        """Initialize a LinkGraph object.

        Examples:
            Create a LinkGraph object:
            >>> lg = LinkGraph()

            Add a link between a GCF and a Spectrum object:
            >>> lg.add_link(gcf, spectrum, metcalf=Score("metcalf", 1.0, {"cutoff": 0.5}))

            Get all links for a given object:
            >>> lg[gcf]
            {spectrum: {"metcalf": Score("metcalf", 1.0, {"cutoff": 0.5})}}

            Get all links:
            >>> lg.links
            [(gcf, spectrum, {"metcalf": Score("metcalf", 1.0, {"cutoff": 0.5})})]

            Check if there is a link between two objects:
            >>> lg.has_link(gcf, spectrum)
            True

            Get the link data between two objects:
            >>> lg.get_link_data(gcf, spectrum)
            {"metcalf": Score("metcalf", 1.0, {"cutoff": 0.5})}
        """
        self._g: Graph = Graph()

    def __str__(self) -> str:
        """Get a short summary of the LinkGraph."""
        return f"{self.__class__.__name__}(#links={len(self.links)}, #objects={len(self)})"

    def __len__(self) -> int:
        """Get the number of objects."""
        return len(self._g)

    @validate_u
    def __getitem__(self, u: Entity) -> dict[Entity, LINK_DATA]:
        """Get all links for a given object.

        Args:
            u: the given object

        Returns:
            A dictionary of links for the given object.

        Raises:
            KeyError: if the input object is not found in the link graph.
        """
        try:
            links = self._g[u]
        except KeyError:
            raise KeyError(f"{u} not found in the link graph.")

        return {**links}  # type: ignore

    @property
    def links(
        self,
    ) -> list[LINK]:
        """Get all links.

        Returns:
            A list of tuples containing the links between objects.
        """
        return list(self._g.edges(data=True))

    @validate_uv
    def add_link(
        self,
        u: Entity,
        v: Entity,
        **data: Score,
    ) -> None:
        """Add a link between two objects.

        The objects `u` and `v` must be different types, i.e. one must be a GCF and the other must be
        a Spectrum or MolecularFamily.

        Args:
            u: the first object, either a GCF, Spectrum, or MolecularFamily
            v: the second object, either a GCF, Spectrum, or MolecularFamily
            data: keyword arguments. At least one scoring method and its data must be provided.
                The key must be the name of the scoring method defined in `ScoringMethod`, and the
                value is a `Score` object, e.g. `metcalf=Score("metcalf", 1.0, {"cutoff": 0.5})`.
        """
        # validate the data
        if not data:
            raise ValueError("At least one scoring method and its data must be provided.")
        for key, value in data.items():
            if not ScoringMethod.has_value(key):
                raise ValueError(
                    f"{key} is not a valid name of scoring method. See `ScoringMethod` for valid names."
                )
            if not isinstance(value, Score):
                raise TypeError(f"{value} is not a Score object.")

        self._g.add_edge(u, v, **data)

    @validate_uv
    def has_link(self, u: Entity, v: Entity) -> bool:
        """Check if there is a link between two objects.

        Args:
            u: the first object, either a GCF, Spectrum, or MolecularFamily
            v: the second object, either a GCF, Spectrum, or MolecularFamily

        Returns:
            True if there is a link between the two objects, False otherwise
        """
        return self._g.has_edge(u, v)

    @validate_uv
    def get_link_data(
        self,
        u: Entity,
        v: Entity,
    ) -> LINK_DATA | None:
        """Get the data for a link between two objects.

        Args:
            u: the first object, either a GCF, Spectrum, or MolecularFamily
            v: the second object, either a GCF, Spectrum, or MolecularFamily

        Returns:
            A dictionary of scoring methods and their data for the link between the two objects, or
            None if there is no link between the two objects.
        """
        return self._g.get_edge_data(u, v)  # type: ignore
