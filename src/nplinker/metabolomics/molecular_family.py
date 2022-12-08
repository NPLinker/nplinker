from typing_extensions import Self

from nplinker.metabolomics.spectrum import Spectrum
from nplinker.strain_collection import StrainCollection

class MolecularFamily():

    def __init__(self, family_id: int):
        """Class to model molecular families.

        Args:
            family_id(int): Id for the molecular family.
        """
        self.id: int = -1
        self.family_id: int = family_id
        self.spectra: list[Spectrum] = []
        self.family = None
        self.spectra_ids: set[int] = set()

    # def has_strain(self, strain):
    #     for spectrum in self.spectra:
    #         if spectrum.has_strain(strain):
    #             return True

    #     return False

    @property
    def strains(self) -> StrainCollection:
        """Get strains of spectra in the molecular family.

        Returns:
            StrainCollection: StrainCollection of strains from which the spectra in the molecular family are coming.
        """
        strains = set()
        for spectrum in self.spectra:
            strains = strains.union(spectrum.strains)
        return strains

    def add_spectrum(self, spectrum: Spectrum):
        """Add a spectrum to the list of contained spectra.

        Args:
            spectrum(Spectrum): Spectrum to add to the molecular family.
        """
        self.spectra.append(spectrum)

    def __str__(self) -> str:
        return 'MolFam(family_id={}, spectra={})'.format(
            self.family_id, len(self.spectra))

    def __eq__(self, other: Self) -> bool:
        return self.id == other.id

    def __hash__(self) -> int:
        return self.id
