class MolecularFamily():

    def __init__(self, family_id):
        self.id = -1
        self.family_id = family_id
        self.spectra = []
        self.family = None
        self.spectra_ids = []

    # def has_strain(self, strain):
    #     for spectrum in self.spectra:
    #         if spectrum.has_strain(strain):
    #             return True

    #     return False

    @property
    def strains(self):
        strains = set()
        for spectrum in self.spectra:
            strains = strains.union(spectrum.strains)
        return strains

    def add_spectrum(self, spectrum):
        self.spectra.append(spectrum)

    def __str__(self):
        return 'MolFam(family_id={}, spectra={})'.format(
            self.family_id, len(self.spectra))

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id
