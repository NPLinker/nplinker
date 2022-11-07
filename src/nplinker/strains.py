from .logconfig import LogConfig

logger = LogConfig.getLogger(__file__)


class Strain():

    def __init__(self, primary_strain_id):
        self.id = primary_strain_id
        self.aliases = set()

    def has_alias(self, name: str) -> bool:
        """Check if alias name exist

        Args:
            name(str): strain alias name

        Returns:
            bool: if given name exist return True
        """
        return name in self.aliases

    def add_alias(self, name: str) -> None:
        """Add one alias name to aliases set

        Args:
            name(str): strain alias name
        """
        if len(name) == 0:
            logger.warning(
                f'Refusing to add zero-length alias to strain {self}')
            return

        self.aliases.add(name)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f'Strain({self.id}) [{len(self.aliases)} aliases]'
