import json
from os import PathLike
from pathlib import Path
from typing import Iterator
from deprecated import deprecated
from .logconfig import LogConfig
from .strains import Strain
from .utils import list_dirs
from .utils import list_files

logger = LogConfig.getLogger(__name__)


class StrainCollection():

    def __init__(self):
        """A collection of Strain objects."""
        self._strains: list[Strain] = []
        # dict of strain name (id and alias) to primary strain object
        self._strain_dict_name: dict[str, Strain] = {}

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if len(self) > 20:
            return f'StrainCollection(n={len(self)})'

        return f'StrainCollection(n={len(self)}) [' + ','.join(
            s.id for s in self._strains) + ']'

    def __len__(self) -> int:
        return len(self._strains)

    def __eq__(self, other) -> bool:
        if isinstance(other, StrainCollection):
            return (self._strains == other._strains
                    and self._strain_dict_name == other._strain_dict_name)
        return NotImplemented

    def __contains__(self, item: str | Strain) -> bool:
        """Check if the strain collection contains the given strain.

        The given strain could be a Strain object, or a strain id or alias.
        """
        if isinstance(item, str):
            return item in self._strain_dict_name
        if isinstance(item, Strain):
            return item.id in self._strain_dict_name
        raise TypeError(f"Expected Strain or str, got {type(item)}")

    def __iter__(self) -> Iterator[Strain]:
        return iter(self._strains)

    def add(self, strain: Strain) -> None:
        """Add strain to the collection.

        If the strain already exists, merge the aliases.

        Args:
            strain(Strain): The strain to add.
        """
        # if the strain exists, merge the aliases
        if strain.id in self._strain_dict_name:
            existing: Strain = self.lookup(strain.id)
            for alias in strain.aliases:
                existing.add_alias(alias)
                self._strain_dict_name[alias] = existing
        else:
            self._strains.append(strain)
            self._strain_dict_name[strain.id] = strain
            for alias in strain.aliases:
                self._strain_dict_name[alias] = strain

    def remove(self, strain: Strain):
        """Remove a strain from the collection.

        Args:
            strain(Strain): The strain to remove.
        """
        if strain.id in self._strain_dict_name:
            self._strains.remove(strain)
            del self._strain_dict_name[strain.id]
            for alias in strain.aliases:
                del self._strain_dict_name[alias]

    def filter(self, strain_set: set[Strain]):
        """
        Remove all strains that are not in strain_set from the strain collection
        """
        # note that we need to copy the list of strains, as we are modifying it
        for strain in self._strains.copy():
            if strain not in strain_set:
                self.remove(strain)

    def lookup(self, name: str) -> Strain:
        """Lookup a strain by name (id or alias).

        Args:

            name(str): Strain name (id or alias) to lookup.

        Returns:
            Strain: Strain identified by the given name.

        Raises:
            KeyError: If the strain name is not found.
        """
        if name in self:
            return self._strain_dict_name[name]
        raise KeyError(f"Strain {name} not found in strain collection.")

    # TODO to move this method to a separate class
    @deprecated(version="1.3.3", reason="This method will be removed")
    def generate_strain_mappings(self, strain_mappings_file: str | PathLike,
                                 antismash_dir: str | PathLike) -> None:
        """Add AntiSMASH BGC file names as strain alias to strain mappings file.

            Note that if AntiSMASH ID (e.g. GCF_000016425.1) is not found in
            existing self.strains, its corresponding BGC file names will not be
            added.

        Args:
            strain_mappings_file(str | PathLike): Path to strain mappings file.
            antismash_dir(str | PathLike): Path to AntiSMASH output directory.
        """
        if Path(strain_mappings_file).exists():
            logger.info('Strain mappings file exist')
            return

        # if not exist, generate strain mapping file with antismash BGC names
        logger.info('Generating strain mappings file')
        subdirs = list_dirs(antismash_dir)
        for d in subdirs:
            antismash_id = Path(d).name

            # use antismash_id (e.g. GCF_000016425.1) as strain name to query
            # TODO: self is empty at the moment, why lookup here?
            try:
                strain = self.lookup(antismash_id)
            except KeyError:
                logger.warning(
                    f'Failed to lookup AntiSMASH strain name: {antismash_id}')
                continue

            # if strain `antismash_id` exist, add all gbk file names as strain alias
            gbk_files = list_files(d, suffix=".gbk", keep_parent=False)
            for f in gbk_files:
                gbk_filename = Path(f).stem
                strain.add_alias(gbk_filename)

        logger.info(f'Saving strains to {strain_mappings_file}')
        self.save_to_file(strain_mappings_file)
