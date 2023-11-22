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
        return "MolecularFamily(family_id={}, #Spectrum_objects={}, #spectrum_ids={}, #strains={})".format(
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

    def add_spectrum(self, spectrum: Spectrum) -> None:
        """Add a Spectrum object to the molecular family.

        Args:
            spectrum(Spectrum): `Spectrum` object to add to the molecular family.
        """
        self._spectra.add(spectrum)
        self.spectra_ids.add(spectrum.spectrum_id)
        self._strains = self._strains + spectrum.strains
        # add the molecular family to the spectrum
        spectrum.family = self

    def detach_spectrum(self, spectrum: Spectrum) -> None:
        """Remove a Spectrum object from the molecular family.

        Args:
            spectrum(Spectrum): `Spectrum` object to remove from the molecular family.
        """
        self._spectra.remove(spectrum)
        self.spectra_ids.remove(spectrum.spectrum_id)
        self._strains = self._update_strains()
        # remove the molecular family from the spectrum
        spectrum.family = None

    def is_singleton(self) -> bool:
        """Check if the molecular family contains only one spectrum.

        Returns:
            bool: True when `MolecularFamily.spectra_ids` contains only one spectrum id.
        """
        return len(self.spectra_ids) == 1

    def _update_strains(self) -> StrainCollection:
        """Update strains in the molecular family.

        The strains is regenerated from the existing spectra in the molecular family. This method
        is mainly used when the spectra are updated (e.g. remove a spectrum from the molecular
        family).

        Returns:
            StrainCollection: Updated StrainCollection object.
        """
        strains = StrainCollection()
        for spectrum in self._spectra:
            strains = strains + spectrum.strains
        return strains
