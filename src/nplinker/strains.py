from .logconfig import LogConfig

logger = LogConfig.getLogger(__name__)


class Strain():

    def __init__(self, primary_strain_id: str):
        self.id: str = primary_strain_id
        self.aliases: set[str] = set()

    def has_alias(self, alt_id: str) -> bool:
        """Check if strain has an alias.

        Args:
            alt_id(str): Alias to check.

        Returns:
            bool: Whether the alias is registered in the set of aliases or not.
        """
        return alt_id in self.aliases

    def add_alias(self, alt_id: str):
        """Add an alias to the list of known aliases.

        Args:
            alt_id(str): Alternative id to add to the list of known aliases.
        """
        if len(alt_id) == 0:
            logger.warning(
                f'Refusing to add zero-length alias to strain {self}')
            return

        self.aliases.add(name)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f'Strain({self.id}) [{len(self.aliases)} aliases]'
