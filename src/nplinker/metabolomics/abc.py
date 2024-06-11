from abc import ABC
from abc import abstractmethod
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from typing import Mapping
    from typing import Sequence
    from .molecular_family import MolecularFamily
    from .spectrum import Spectrum


class SpectrumLoaderBase(ABC):
    """Abstract base class for SpectrumLoader."""

    @property
    @abstractmethod
    def spectra(self) -> Sequence["Spectrum"]:
        """Get Spectrum objects.

        Returns:
            A sequence of Spectrum objects.
        """


class MolecularFamilyLoaderBase(ABC):
    """Abstract base class for MolecularFamilyLoader."""

    @abstractmethod
    def get_mfs(self, keep_singleton: bool) -> Sequence["MolecularFamily"]:
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
    def mappings(self) -> Mapping[str, Sequence[str]]:
        """Get file mappings.

        Returns:
            A mapping from spectrum ID to the names of files where the spectrum occurs.
        """


class AnnotationLoaderBase(ABC):
    """Abstract base class for AnnotationLoader."""

    @property
    @abstractmethod
    def annotations(self) -> Mapping[str, Mapping]:
        """Get annotations.

        Returns:
            A mapping from spectrum ID to its annotations.
        """
