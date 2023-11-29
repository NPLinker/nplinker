from __future__ import annotations
from functools import cached_property
from typing import TYPE_CHECKING
import numpy as np
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection


if TYPE_CHECKING:
    from .molecular_family import MolecularFamily

GNPS_KEY = "gnps"


class Spectrum:
    def __init__(
        self,
        spectrum_id: str,
        mz: list[float],
        intensity: list[float],
        precursor_mz: float,
        rt: float | None = None,
        metadata: dict | None = None,
        annotations: dict | None = None,
    ) -> None:
        """Class to model MS/MS Spectrum.

        Args:
            spectrum_id (str): the spectrum ID.
            mz (list[float]): the list of m/z values.
            intensity (list[float]): the list of intensity values.
            precursor_mz (float): the precursor m/z.
            rt (float, optional): the retention time in seconds.
            metadata (dict, optional): the metadata of the spectrum, i.e. the header infomation
                in the MGF file.
            annotations (dict, optional): the annotations of the spectrum, e.g. annotations from
                GNPS.

        Attributes:
            spectrum_id (str): the spectrum ID.
            precursor_mz (float): the precursor m/z.
            rt (float): the retention time in seconds.
            metadata (dict): the metadata of the spectrum, i.e. the header infomation in the MGF
                file.
            annotations (dict): the annotations of the spectrum, e.g. annotations from GNPS.
            gnps_id (str): the GNPS ID of the spectrum.
            strains (StrainCollection): the strains that this spectrum belongs to.
            family (MolecularFamily): the molecular family that this spectrum belongs to.
            peaks (np.ndarray): 2D array of peaks, and each row is a peak of (m/z, intensity).
        """
        self.spectrum_id = spectrum_id
        self.mz = mz
        self.intensity = intensity
        self.precursor_mz = precursor_mz
        self.rt = rt
        self.metadata = metadata or {}
        self.annotations = annotations or {}

        self.gnps_id = None
        self.strains = StrainCollection()
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

    @property
    def gnps_annotations(self):
        if GNPS_KEY not in self.annotations:
            return None

        return self.annotations[GNPS_KEY][0]

    def has_strain(self, strain: Strain):
        return strain in self.strains
