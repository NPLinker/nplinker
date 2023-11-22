from __future__ import annotations
from typing import TYPE_CHECKING
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection


if TYPE_CHECKING:
    from .spectrum import Spectrum


class MolecularFamily:
    def __init__(self, family_id: str):
        """Class to model molecular family.

        Args:
            family_id(str): Id for the molecular family.
        """
        self.id: int = -1
        self.family_id: str = family_id
        self.spectra_ids: set[str] = set()
        self._spectra: set[Spectrum] = set()
        self._strains: StrainCollection = StrainCollection()

    def __str__(self) -> str:
        return "MF(family_id={}, #Spectrum_objects={}, #spectrum_ids={}, #strains={})".format(
            self.family_id, len(self._spectra), len(self.spectra_ids), len(self._strains)
        )

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        if isinstance(other, MolecularFamily):
            return self.id == other.id and self.family_id == other.family_id
        return NotImplemented

    def __hash__(self) -> int:
        """Hash function for MolecularFamily.

        Note that MolecularFamily is a mutable container, so here we hash on
        the id and family_id only to avoid the hash value changing when
        `self.spectra` is updated.
        """
        return hash((self.id, self.family_id))

    @property
    def spectra(self) -> set[Spectrum]:
        """Get Spectrum objects in the molecular family."""
        return self._spectra

    @property
    def strains(self) -> StrainCollection:
        """Get strains in the molecular family."""
        return self._strains

    def has_strain(self, strain: Strain) -> bool:
        """Check if the given strain exists.

        Args:
            strain(Strain): `Strain` object.

        Returns:
            bool: True when the given strain exist.
        """
        return strain in self._strains

    # TODO: update the logics, mf should also be added to the spectrum object
    def add_spectrum(self, spectrum: Spectrum):
        """Add a spectrum to the spectra list.

        Args:
            spectrum(Spectrum): Spectrum to add to the molecular family.
        """
        self._spectra.add(spectrum)
