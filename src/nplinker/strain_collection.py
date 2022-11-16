import csv
import os
from .logconfig import LogConfig
from .strains import Strain


logger = LogConfig.getLogger(__file__)


class StrainCollection():

    def __init__(self):
        self._strains: list[Strain] = []
        self._lookup: dict[str, Strain] = {}
        self._lookup_indices: dict[int, Strain] = {}

    def add(self, strain: Strain):
        """Add the strain to the aliases.
        This also adds those strain's aliases to this' strain's aliases.

        Args:
            strain(Strain): Strain to add to self.

        Examples:
            >>> 
            """
        if strain.id in self._lookup:
            # if it already exists, just merge the set of aliases and update
            # lookup entries
            existing = self.lookup(strain.id)
            existing.aliases.update(strain.aliases)
            for alias in strain.aliases:
                self._lookup[alias] = existing
            return

        self._lookup_indices[len(self)] = strain
        self._strains.append(strain)
        # insert a mapping from strain=>strain, plus all its aliases
        self._lookup[strain.id] = strain
        for alias in strain.aliases:
            self._lookup[alias] = strain

    def remove(self, strain: Strain):
        """Remove the specified strain from the aliases.
        TODO: #90 Implement removing the strain also from self._lookup indices.

        Args:
            strain(Strain): Strain to remove.
        """
        if strain.id not in self._lookup:
            return

        self._strains.remove(strain)
        del self._lookup[strain.id]
        for alias in strain.aliases:
            del self._lookup[alias]

    def filter(self, strain_set: set[Strain]):
        """
        Remove all strains that are not in strain_set from the strain collection
        """
        to_remove = [x for x in self._strains if x not in strain_set]
        for strain in to_remove:
            self.remove(strain)

    def __contains__(self, strain_id: str|Strain) -> bool:
        """Check if the strain or strain id are contained in the lookup table.

        Args:
            strain_id(str|Strain): Strain or strain id to look up.

        Returns:
            bool: Whether the strain is contained in the collection.
         """
         
        if isinstance(strain_id, str):
            return strain_id in self._lookup
        # assume it's a Strain object
        return strain_id.id in self._lookup

    def __iter__(self):
        return iter(self._strains)

    def __next__(self):
        return next(self._strains)

    def lookup_index(self, index: int) -> Strain:
        """Return the strain from lookup by index.

        Args:
            index(int): Position index from which to retrieve the strain

        Returns:
            Strain: Strain identified by the given index.
        """
        return self._lookup_indices[index]

    def lookup(self, strain_id: str, default=None) -> Strain|None:
        """Check whether the strain id is contained in the lookup table. If so, return the strain, otherwise return `default`.

        Args:
            strain_id(str): Strain id to lookup.
            default (Strain, optional): Object to return otherwise. Defaults to None.

        Returns:
            Optional[Strain]: Strain retrieved during lookup or object passed as default.
        """
        if strain_id not in self._lookup:
            # logger.error('Strain lookup failed for "{}"'.format(strain_id))
            return default

        return self._lookup[strain_id]

    def add_from_file(self, filename: str):
        """Read strains and aliases from file and store in self.

        Args:
            filename(str): Path to strain mapping file to load.
        """

        if not os.path.exists(filename):
            logger.warning(
                f'strain mappings file not found: {filename}')
            return

        line = 1
        with open(filename) as f:
            reader = csv.reader(f)
            for ids in reader:
                if len(ids) == 0:
                    continue
                strain = Strain(ids[0])
                for id in ids[1:]:
                    if len(id) == 0:
                        logger.warning(
                            'Found zero-length strain label in {} on line {}'.
                            format(filename, line))
                    else:
                        strain.add_alias(id)
                self.add(strain)

                line += 1

    def save_to_file(self, filename: str):
        """Save this strain collection to file.

        Args:
            filename(str): Output filename.

        Examples:
            >>> 
            """
        with open(filename, 'w') as f:
            for strain in self._strains:
                ids = [strain.id] + list(strain.aliases)
                f.write(','.join(ids) + '\n')
    
    def generate_strain_mappings(self, strain_mappings_file: str, antismash_dir: str):
        """Generate stain mappings file from this strain collection, exporting the strains and their aliases.

        Args:
            strain_mappings_file(str): Export path.
            antismash_dir(str): Directory in which to look for gbk files for which to export the aliases.

        Returns:
            StrainCollection: self
        """
        # first time downloading, this file will not exist, should only need done once
        if not os.path.exists(strain_mappings_file):
            logger.info('Generating strain mappings file')
            for root, dirs, files in os.walk(antismash_dir):
                for f in files:
                    if not f.endswith('.gbk'):
                        continue
                    
                    # use the containing folder of the .gbk file as the strain name,
                    # and then take everything but ".gbk" from the filename and use
                    # that as an alias
                    strain_name = os.path.split(root)[1]
                    strain_alias = os.path.splitext(f)[0]
                    if strain_alias.find('.') != -1:
                        strain_alias = strain_alias[:strain_alias.index('.')]
                    if self.lookup(strain_name) is not None:
                        self.lookup(strain_name).add_alias(strain_alias)
                    else:
                        logger.warning(
                            f'Failed to lookup strain name: {strain_name}')
            logger.info(f'Saving strains to {strain_mappings_file}')
            self.save_to_file(strain_mappings_file)
        else:
            logger.info('Strain mappings already generated')
    
        return self

    def __len__(self) -> int:
        return len(self._strains)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self):
        if len(self) > 20:
            return f'StrainCollection(n={len(self)})'

        return f'StrainCollection(n={len(self)}) [' + ','.join(
            s.id for s in self._strains) + ']'
