from abc import ABC
from abc import abstractmethod
from collections.abc import Sequence
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from .molecular_family import MolecularFamily
    from .spectrum import Spectrum


class SpectrumLoaderBase(ABC):
    @property
    @abstractmethod
    def spectra(self) -> Sequence["Spectrum"]:
        ...


class MolecularFamilyLoaderBase(ABC):
    @abstractmethod
    def get_mfs(self, keep_singleton: bool) -> Sequence["MolecularFamily"]:
        """Get MolecularFamily objects.

        Args:
            keep_singleton: True to keep singleton molecular families. A
                singleton molecular family is a molecular family that contains
                only one spectrum.

        Returns:
            Sequence[MolecularFamily]: a list of MolecularFamily objects.
        """


class FileMappingLoaderBase(ABC):
    @property
    @abstractmethod
    def mappings(self) -> dict[str, list[str]]:
        ...


class AnnotationLoaderBase(ABC):
    @property
    @abstractmethod
    def annotations(self) -> dict[str, dict]:
        ...
