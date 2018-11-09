import numpy as np

# TODO currently not used
# needed if its only ever going to be a wrapper around a string?
class Strain(object):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return 'Strain({})'.format(self.name)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        # two strains are equal if they have the same name
        return isinstance(other, Strain) and self.name == other.name

    def __hash__(self):
        return hash(self.name)

class StrainManager:
    """Stores a list of all strains present in a dataset (from both genomics and 
    metabolomics sources) and provides some utility methods for dealing with scoring etc"""

    def __init__(self, initial_strains):
        """Initialise the manager object with a collection of strains (should be 
        a list or set of strain IDs (strings))"""
        if isinstance(initial_strains, list):
            # get rid of duplicates
            initial_strains = set(initial_strains)
        self.strains = initial_strains
        self._generate_indexes()

    def add_strains(self, new_strains):
        """Add new strains to the existing set.

        Note that existing lookup tables will need recreated after doing this"""
        self.strains.update(new_strains)
        self._generate_indexes()

    def _generate_indexes(self):
        """Create a dict mapping strain names to indices in a lookup table"""
        self.indexes = {}
        for i, s in enumerate(self.strains):
            self.indexes[s] = i

    @property
    def all_strains_set(self):
        return self.strains

    @property
    def all_strains_np(self):
        return np.array(list(self.strains))

    def generate_lookup_table(self, strains):
        """Given a subset of the full set of known strains, generate a lookup table 
        containing True for the indices corresponding to those strains"""
        lookup = np.full((len(self.strains), ), False)
        indices = [self.indexes[s] for s in strains]
        lookup[indices] = True
        return lookup

    def generate_prob_dict(self, spectra):
        """Generates strain probability dict using code from example notebook"""
        total = 0
        strain_counts = {s: 0 for s in self.strains}
        for spectrum in spectra:
            for strain in strain_counts.keys():
                if spectrum.has_strain(strain):
                    strain_counts[strain] += 1
                    total += 1

        return {s: v/total for s, v in strain_counts.items()}

