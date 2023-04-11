from __future__ import annotations
from .logconfig import LogConfig
from typing import Iterator


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
        if isinstance(other, Strain):
            return (self.id == other.id
                    and self.aliases == other.aliases)
        return NotImplemented

    def __hash__(self) -> int:
        """Hash function for Strain.

        Note that Strain is a mutable container, so here we hash on only the id
        to avoid the hash value changes when `self._aliases` is updated.
        """
        return hash(self.id)

    def __contains__(self, alias: str) -> bool:
        if not isinstance(alias, str):
            raise TypeError(f'Expected str, got {type(alias)}')
        return alias in self._aliases

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
        if not isinstance(alias, str):
            raise TypeError(f'Expected str, got {type(alias)}')
        if len(alias) == 0:
            logger.warning(
                'Refusing to add an empty-string alias to strain {%s}', self)
        else:
            self._aliases.add(alias)
