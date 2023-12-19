from __future__ import annotations
from functools import cached_property
from typing import TYPE_CHECKING
import numpy as np
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection


if TYPE_CHECKING:
    from .molecular_family import MolecularFamily


class Spectrum:
    def __init__(
        self,
        spectrum_id: str,
        mz: list[float],
        intensity: list[float],
        precursor_mz: float,
        rt: float = 0,
        metadata: dict | None = None,
    ) -> None:
        """Class to model MS/MS Spectrum.

        Args:
            spectrum_id (str): the spectrum ID.
            mz (list[float]): the list of m/z values.
            intensity (list[float]): the list of intensity values.
            precursor_mz (float): the precursor m/z.
            rt (float): the retention time in seconds. Defaults to 0.
            metadata (dict, optional): the metadata of the spectrum, i.e. the header infomation
                in the MGF file.

        Attributes:
            spectrum_id (str): the spectrum ID.
            mz (list[float]): the list of m/z values.
            intensity (list[float]): the list of intensity values.
            precursor_mz (float): the m/z value of the precursor.
            rt (float): the retention time in seconds.
            metadata (dict): the metadata of the spectrum, i.e. the header infomation in the MGF
                file.
            gnps_annotations (dict): the GNPS annotations of the spectrum.
            gnps_id (str | None): the GNPS ID of the spectrum.
            strains (StrainCollection): the strains that this spectrum belongs to.
            family (MolecularFamily): the molecular family that this spectrum belongs to.
            peaks (np.ndarray): 2D array of peaks, each row is a peak of (m/z, intensity) values.
        """
        self.spectrum_id = spectrum_id
        self.mz = mz
        self.intensity = intensity
        self.precursor_mz = precursor_mz
        self.rt = rt
        self.metadata = metadata or {}

        self.gnps_annotations: dict = {}
        self.gnps_id: str | None = None
        self.strains: StrainCollection = StrainCollection()
        self.family: MolecularFamily | None = None

    def __str__(self) -> str:
        return f"Spectrum(spectrum_id={self.spectrum_id}, #strains={len(self.strains)})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        if isinstance(other, Spectrum):
            return self.spectrum_id == other.spectrum_id and self.precursor_mz == other.precursor_mz
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.spectrum_id, self.precursor_mz))

    @cached_property
    def peaks(self) -> np.ndarray:
        """Get the peaks, a 2D array with each row containing the values of (m/z, intensity)."""
        return np.array(list(zip(self.mz, self.intensity)))

    def has_strain(self, strain: Strain):
        """Check if the given strain exists in the spectrum.

        Args:
            strain(Strain): `Strain` object.

        Returns:
            bool: True when the given strain exist in the spectrum.
        """
        return strain in self.strains
