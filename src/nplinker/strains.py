from __future__ import annotations
from .logconfig import LogConfig


logger = LogConfig.getLogger(__name__)


class Strain():

    def __init__(self, primary_id: str) -> None:
        """To model the mapping between strain id and its aliases.

        It's recommended to use NCBI taxonomy strain id or name as the primary
        id.

        Args:
            primary_id(str): the representative id of the strain.
        """
        self.id: str = primary_id
        self._aliases: set[str] = set()

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f'Strain({self.id}) [{len(self._aliases)} aliases]'

    def __eq__(self, other) -> bool:
        return (isinstance(other, Strain) and self.id == other.id
                and self._aliases == other._aliases)

    def __hash__(self) -> int:
        return hash(self.id)

    @property
    def aliases(self) -> set[str]:
        """Get the set of known aliases.

        Returns:
            set[str]: A set of aliases associated with the strain.
        """
        return self._aliases

    def add_alias(self, alias: str) -> None:
        """Add an alias to the list of known aliases.

        Args:
            alias(str): The alias to add to the list of known aliases.
        """
        if len(alias) == 0:
            logger.warning(
                'Refusing to add an empty-string alias to strain {%s}', self)
        else:
            self._aliases.add(alias)
