import csv

from .logconfig import LogConfig
logger = LogConfig.getLogger(__file__)

class Strain(object):

    def __init__(self, primary_strain_id):
        self.id = primary_strain_id
        self.aliases = set()
    
    def has_alias(self, alt_id):
        return alt_id in self.aliases

    def add_alias(self, alt_id):
        self.aliases.add(alt_id)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'Strain({}) [{} aliases]'.format(self.id, len(self.aliases))

class StrainCollection(object):

    def __init__(self):
        self._strains = []
        self._lookup = {}
        self._lookup_indices = {}

    def add(self, strain):
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

    def remove(self, strain):
        if strain.id not in self._lookup:
            print(strain.id)
            return

        self._strains.remove(strain)
        del self._lookup[strain.id]
        for alias in strain.aliases:
            del self._lookup[alias]

    def __contains__(self, strain_id):
        if isinstance(strain_id, str):
            return strain_id in self._lookup
        # assume it's a Strain object
        return strain_id.id in self._lookup

    def __iter__(self):
        return iter(self._strains)

    def __next__(self):
        return next(self._strains)

    def lookup_index(self, index):
        return self._lookup_indices[index]

    def lookup(self, strain_id, default=None):
        if strain_id not in self._lookup:
            # logger.error('Strain lookup failed for "{}"'.format(strain_id))
            return default

        return self._lookup[strain_id]

    def add_from_file(self, filename):
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            for ids in reader:
                if len(ids) == 0:
                    continue
                strain = Strain(ids[0])
                for id in ids[1:]:
                    strain.add_alias(id)
                self.add(strain)

    def __len__(self):
        return len(self._strains)

    def __repr__(self):
        return str(self)

    def __str__(self):
        if len(self) > 20:
            return 'StrainCollection(n={})'.format(len(self))

        return 'StrainCollection(n={}) ['.format(len(self)) + ','.join(s.id for s in self._strains) + ']'
