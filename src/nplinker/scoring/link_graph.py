from __future__ import annotations
from collections.abc import Sequence
from functools import wraps
from os import PathLike
from typing import Union
from networkx import Graph
from tabulate import tabulate
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

            Display the empty LinkGraph object:
            >>> lg
            |    | Genomic Object Type   | Genomic Object ID   | Metabolomic Object Type   | Metabolomic Object ID   | Metcalf Score   | Rosetta Score   |
            |----|-----------------------|---------------------|---------------------------|-------------------------|-----------------|-----------------|

            Add a link between a GCF and a Spectrum object:
            >>> lg.add_link(gcf, spectrum, metcalf=Score("metcalf", 1.0, {"cutoff": 0.5}))

            Display all links in LinkGraph object:
            >>> lg
            |    | Genomic Object Type   | Genomic Object ID   | Metabolomic Object Type   | Metabolomic Object ID   | Metcalf Score   | Rosetta Score   |
            |----|-----------------------|---------------------|---------------------------|-------------------------|-----------------|-----------------|
            |  1 | GCF                   | 1                   | Spectrum                  | 1                       | 1.00            | -               |

            Get all links for a given object:
            >>> lg[gcf]
            {spectrum: {"metcalf": Score("metcalf", 1.0, {"cutoff": 0.5})}}

            Get all links in the LinkGraph:
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

    def __repr__(self) -> str:
        """Return a string representation of the LinkGraph."""
        return self._get_table_repr()

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

        Examples:
            >>> lg.links
            [(gcf, spectrum, {"metcalf": Score("metcalf", 1.0, {"cutoff": 0.5})})]
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

        Examples:
            >>> lg.add_link(gcf, spectrum, metcalf=Score("metcalf", 1.0, {"cutoff": 0.5}))
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

        Examples:
            >>> lg.has_link(gcf, spectrum)
            True
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

        Examples:
            >>> lg.get_link_data(gcf, spectrum)
            {"metcalf": Score("metcalf", 1.0, {"cutoff": 0.5})}
        """
        return self._g.get_edge_data(u, v)  # type: ignore

    def filter(self, u_nodes: Sequence[Entity], v_nodes: Sequence[Entity] = [], /) -> LinkGraph:
        """Return a new LinkGraph object with the filtered links between the given objects.

        The new LinkGraph object will only contain the links between `u_nodes` and `v_nodes`.

        If `u_nodes` or `v_nodes` is empty, the new LinkGraph object will contain the links for
        the given objects in `v_nodes` or `u_nodes`, respectively. If both are empty, return an
        empty LinkGraph object.

        Note that not all objects in `u_nodes` and `v_nodes` need to be present in the original
        LinkGraph.

        Args:
            u_nodes: a sequence of objects used as the first object in the links
            v_nodes: a sequence of objects used as the second object in the links

        Returns:
            A new LinkGraph object with the filtered links between the given objects.

        Examples:
            Filter the links for `gcf1` and `gcf2`:
            >>> new_lg = lg.filter([gcf1, gcf2])
            Filter the links for `spectrum1` and `spectrum2`:
            >>> new_lg = lg.filter([spectrum1, spectrum2])
            Filter the links between two lists of objects:
            >>> new_lg = lg.filter([gcf1, gcf2], [spectrum1, spectrum2])
        """
        lg = LinkGraph()

        # exchange u_nodes and v_nodes if u_nodes is empty but v_nodes not
        if len(u_nodes) == 0 and len(v_nodes) != 0:
            u_nodes = v_nodes
            v_nodes = []

        if len(v_nodes) == 0:
            for u in u_nodes:
                self._filter_one_node(u, lg)

        for u in u_nodes:
            for v in v_nodes:
                self._filter_two_nodes(u, v, lg)

        return lg

    @validate_u
    def _filter_one_node(self, u: Entity, lg: LinkGraph) -> None:
        """Filter the links for a given object and add them to the new LinkGraph object."""
        try:
            links = self[u]
        except KeyError:
            pass
        else:
            for node2, value in links.items():
                lg.add_link(u, node2, **value)

    @validate_uv
    def _filter_two_nodes(self, u: Entity, v: Entity, lg: LinkGraph) -> None:
        """Filter the links between two objects and add them to the new LinkGraph object."""
        link_data = self.get_link_data(u, v)
        if link_data is not None:
            lg.add_link(u, v, **link_data)

    def get_table_data(self, display_limit: int | None = None) -> list[dict[str, any]]:
        """Generate the table data for the LinkGraph.

        This method iterates over the links in the LinkGraph and constructs a table
        containing information about genomic and metabolomic objects, as well as their
        associated scores. Each row in the table represents a link between a genomic
        object and a metabolomic object.

        Args:
            display_limit (int | None): The maximum number of rows to include in the
                table. If None, all rows are included.

        Returns:
            list: A list of dictionaries, where each dictionary contains
                the following keys:
                - Index (int)
                - Genomic Object Type (str)
                - Genomic Object ID (str or int)
                - Metabolomic Object Type (str)
                - Metabolomic Object ID (str or int)
                - Metcalf Score (str, formatted to 2 decimal places, or "-")
                - Rosetta Score (str, formatted to 2 decimal places, or "-")
        """
        genomic_object_classes = (GCF,)

        table_data = []

        for index, (u, v, data) in enumerate(self.links, start=1):
            genomic_object = u if isinstance(u, genomic_object_classes) else v
            metabolomic_object = v if isinstance(u, genomic_object_classes) else u
            metcalf_score = data.get("metcalf")
            rosetta_score = data.get("rosetta")

            table_data.append(
                {
                    "index": index,
                    "genomic_object_type": genomic_object.__class__.__name__,
                    "genomic_object_id": genomic_object.id,
                    "metabolomic_object_type": metabolomic_object.__class__.__name__,
                    "metabolomic_object_id": metabolomic_object.id,
                    "metcalf_score": f"{metcalf_score.value:.2f}" if metcalf_score else "-",
                    "rosetta_score": f"{rosetta_score.value:.2f}" if rosetta_score else "-",
                }
            )

            if display_limit is not None and index == display_limit:
                break

        return table_data

    def _get_table_repr(self, display_limit: int | None = 60) -> str:
        """Generate a table representation of the LinkGraph.

        Args:
            display_limit: The maximum number of links to display in the table. Defaults to 60.

        Returns:
            str: A string representation of the table in GitHub-flavored markdown format. If the number of links
            exceeds the display limit, the table is truncated and an additional line indicating the total number
            of links is appended.
        """
        table = tabulate(
            self.get_table_data(display_limit),
            headers="keys",
            tablefmt="github",
            stralign="right",
        )

        if len(self.links) > display_limit:
            truncated_info = f"...\n[ {len(self.links)} links ]"
            table += f"\n{truncated_info}"

        return table

    def export_links(self, file: str | PathLike) -> None:
        """Exports the links in the LinkGraph to a file.

        Args:
            file: the file to write the links to.

        Examples:
            >>> lg.print_links("links.tsv")
        """
        table_data = self.get_table_data()
        headers = table_data[0].keys()
        with open(file, "w") as f:
            f.write("\t".join(headers) + "\n")
            for row in table_data:
                f.write("\t".join(str(row[h]) for h in headers) + "\n")
