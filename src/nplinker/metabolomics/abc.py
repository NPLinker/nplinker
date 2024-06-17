from __future__ import annotations
from abc import ABC
from abc import abstractmethod
from .molecular_family import MolecularFamily
from .spectrum import Spectrum


class SpectrumLoaderBase(ABC):
    """Abstract base class for SpectrumLoader."""

    @property
    @abstractmethod
    def spectra(self) -> list[Spectrum]:
        """Get Spectrum objects.

        Returns:
            A sequence of Spectrum objects.
        """


class MolecularFamilyLoaderBase(ABC):
    """Abstract base class for MolecularFamilyLoader."""

    @abstractmethod
    def get_mfs(self, keep_singleton: bool) -> list[MolecularFamily]:
        """Get MolecularFamily objects.

        Args:
            keep_singleton: True to keep singleton molecular families. A
                singleton molecular family is a molecular family that contains
                only one spectrum.

        Returns:
            A sequence of MolecularFamily objects.
        """


class FileMappingLoaderBase(ABC):
    """Abstract base class for FileMappingLoader."""

    @property
    @abstractmethod
    def mappings(self) -> dict[str, list[str]]:
        """Get file mappings.

        Returns:
            A mapping from spectrum ID to the names of files where the spectrum occurs.
        """


class AnnotationLoaderBase(ABC):
    """Abstract base class for AnnotationLoader."""

    @property
    @abstractmethod
    def annotations(self) -> dict[str, dict]:
        """Get annotations.

        Returns:
            A mapping from spectrum ID to its annotations.
        """
