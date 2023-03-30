import csv
import os
from typing import Iterator
from .logconfig import LogConfig
from .strains import Strain
from .utils import list_dirs
from .utils import list_files


logger = LogConfig.getLogger(__name__)


class StrainCollection():

    def __init__(self):
        """A collection of Strain objects."""
        self._strains: list[Strain] = []
        self._strain_dict_id: dict[str, Strain] = {}
        self._strain_dict_index: dict[int, Strain] = {}

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
        result = self._strains == other._strains
        result &= self._strain_dict_id == other._strain_dict_id
        result &= self._strain_dict_index == other._strain_dict_index
        return result

    def __contains__(self, strain: str | Strain) -> bool:
        if isinstance(strain, str):
            value = strain in self._strain_dict_id
        elif isinstance(strain, Strain):
            value = strain.id in self._strain_dict_id
        else:
            raise TypeError(f"Expected Strain or str, got {type(strain)}")
        return value

    def __iter__(self) -> Iterator[Strain]:
        return iter(self._strains)

    def add(self, strain: Strain) -> None:
        """Add strain to the collection.

        If the strain already exists, merge the aliases.

        Args:
            strain(Strain): The strain to add.
        """
        # if the strain exists, merge the aliases
        if strain.id in self._strain_dict_id:
            existing: Strain = self.lookup(strain.id)
            for alias in strain.aliases:
                existing.add_alias(alias)
                self._strain_dict_id[alias] = existing
            return

        self._strain_dict_index[len(self)] = strain
        self._strains.append(strain)
        self._strain_dict_id[strain.id] = strain
        for alias in strain.aliases:
            self._strain_dict_id[alias] = strain

    def remove(self, strain: Strain):
        """Remove a strain from the collection.

        Args:
            strain(Strain): The strain to remove.
        """
        if strain.id in self._strain_dict_id:
            self._strains.remove(strain)
            # remove from dict id
            del self._strain_dict_id[strain.id]
            for alias in strain.aliases:
                del self._strain_dict_id[alias]

    def filter(self, strain_set: set[Strain]):
        """
        Remove all strains that are not in strain_set from the strain collection
        """
        to_remove = [x for x in self._strains if x not in strain_set]
        for strain in to_remove:
            self.remove(strain)

    def lookup_index(self, index: int) -> Strain:
        """Return the strain from lookup by index.

        Args:
            index(int): Position index from which to retrieve the strain

        Returns:
            Strain: Strain identified by the given index.
        """
        return self._strain_dict_index[index]

    def lookup(self, strain_id: str) -> Strain:
        """Check whether the strain id is contained in the lookup table. If so, return the strain, otherwise return `default`.

        Raises:
            Exception if strain_id is not found.

        Args:
            strain_id(str): Strain id to lookup.

        Returns:
            Strain: Strain retrieved during lookup or object passed as default.
        """
        if strain_id not in self._strain_dict_id:
            # logger.error('Strain lookup failed for "{}"'.format(strain_id))
            raise KeyError(f'Strain lookup failed for "{strain_id}"')

        return self._strain_dict_id[strain_id]

    def add_from_file(self, file: str | os.PathLike):
        """Read strains and aliases from file and store in self.

        Args:
            file(str): Path to strain mapping file to load.
        """

        if not os.path.exists(file):
            logger.warning(f'strain mappings file not found: {file}')
            return

        line = 1
        with open(file) as f:
            reader = csv.reader(f)
            for ids in reader:
                if len(ids) == 0:
                    continue
                strain = Strain(ids[0])
                for id in ids[1:]:
                    if len(id) == 0:
                        logger.warning(
                            'Found zero-length strain label in {} on line {}'.
                            format(file, line))
                    else:
                        strain.add_alias(id)
                self.add(strain)

                line += 1

    def save_to_file(self, file: str | os.PathLike):
        """Save this strain collection to file.

        Args:
            file(str): Output file.

        Examples:
            >>>
            """
        with open(file, 'w') as f:
            for strain in self._strains:
                ids = [strain.id] + list(strain.aliases)
                f.write(','.join(ids) + '\n')

    def generate_strain_mappings(self, strain_mappings_file: str,
                                 antismash_dir: str) -> None:
        """Add AntiSMASH BGC file names as strain alias to strain mappings file.

            Note that if AntiSMASH ID (e.g. GCF_000016425.1) is not found in
            existing self.strains, its corresponding BGC file names will not be
            added.

        Args:
            strain_mappings_file(str): Path to strain mappings file
            antismash_dir(str): Path to AntiSMASH directory
        """
        if os.path.exists(strain_mappings_file):
            logger.info('Strain mappings file exist')
            return

        # if not exist, generate strain mapping file with antismash BGC names
        logger.info('Generating strain mappings file')
        subdirs = list_dirs(antismash_dir)
        for d in subdirs:
            antismash_id = os.path.basename(d)

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
                gbk_filename = os.path.splitext(f)[0]
                strain.add_alias(gbk_filename)

        logger.info(f'Saving strains to {strain_mappings_file}')
        self.save_to_file(strain_mappings_file)
