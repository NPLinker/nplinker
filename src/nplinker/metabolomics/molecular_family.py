from typing_extensions import Self

from nplinker.metabolomics.spectrum import Spectrum
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain

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
            set[StrainCollection]: StrainCollection of strains from which the spectra in the molecular family are coming.
        """
        strains: StrainCollection = StrainCollection()
        for spectrum in self.spectra:
            for strain in spectrum.strains:
                strains.add(strain)
        return strains

    def add_spectrum(self, spectrum: Spectrum):
        """Add a spectrum to the spectra list.

        Args:
            spectrum(Spectrum): Spectrum to add to the molecular family.
        """
        self.spectra.append(spectrum)

    def __str__(self) -> str:
        return 'MolFam(family_id={}, spectra={})'.format(
            self.family_id, len(self.spectra))

    def __eq__(self, other: Self) -> bool:
        if isinstance(other, MolecularFamily):
            return (self.id == other.id
                    and self.family_id == other.family_id
                    and set(self.spectra) == set(other.spectra))
        return NotImplemented

    def __hash__(self) -> int:
        """Hash function for MolecularFamily.

        Note that MolecularFamily is a mutable container, so here we hash on
        the id and family_id only to avoid the hash value changing when
        `self.spectra` is updated.
        """
        return hash((self.id, self.family_id))
