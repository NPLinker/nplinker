from __future__ import annotations
from functools import cached_property
from typing import TYPE_CHECKING
import numpy as np
from nplinker.strain import Strain
from nplinker.strain import StrainCollection


if TYPE_CHECKING:
    from .molecular_family import MolecularFamily


class Spectrum:
    """Class to model MS/MS Spectrum.

    Attributes:
        id: the spectrum ID.
        mz: the list of m/z values.
        intensity: the list of intensity values.
        precursor_mz: the m/z value of the precursor.
        rt: the retention time in seconds.
        metadata: the metadata of the spectrum, i.e. the header information in the MGF
            file.
        gnps_annotations: the GNPS annotations of the spectrum.
        gnps_id: the GNPS ID of the spectrum.
        strains: the strains that this spectrum belongs to.
        family: the molecular family that this spectrum belongs to.
        peaks: 2D array of peaks, each row is a peak of (m/z, intensity) values.
    """

    def __init__(
        self,
        id: str,
        mz: list[float],
        intensity: list[float],
        precursor_mz: float,
        rt: float = 0,
        metadata: dict | None = None,
    ) -> None:
        """Initialize the Spectrum.

        Args:
            id: the spectrum ID.
            mz: the list of m/z values.
            intensity: the list of intensity values.
            precursor_mz: the precursor m/z.
            rt: the retention time in seconds. Defaults to 0.
            metadata: the metadata of the spectrum, i.e. the header information
                in the MGF file.
        """
        self.id = id
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
        return f"Spectrum(id={self.id}, #strains={len(self.strains)})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        if isinstance(other, Spectrum):
            return self.id == other.id and self.precursor_mz == other.precursor_mz
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self.id, self.precursor_mz))

    def __reduce__(self) -> tuple:
        """Reduce function for pickling."""
        return (
            self.__class__,
            (self.id, self.mz, self.intensity, self.precursor_mz, self.rt, self.metadata),
            self.__dict__,
        )

    @cached_property
    def peaks(self) -> np.ndarray:
        """Get the peaks, a 2D array with each row containing the values of (m/z, intensity)."""
        return np.array(list(zip(self.mz, self.intensity)))

    def has_strain(self, strain: Strain) -> bool:
        """Check if the given strain exists in the spectrum.

        Args:
            strain: `Strain` object.

        Returns:
            True when the given strain exist in the spectrum.
        """
        return strain in self.strains
