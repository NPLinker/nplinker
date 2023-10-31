from abc import ABC
from abc import abstractmethod
from collections.abc import Sequence
from nplinker.metabolomics import MolecularFamily
from nplinker.metabolomics import Spectrum


class SpectrumLoaderBase(ABC):

    @property
    @abstractmethod
    def spectra(self) -> Sequence[Spectrum]:
        ...

class MolecularFamilyLoaderBase(ABC):

    @property
    @abstractmethod
    def families(self) -> Sequence[MolecularFamily]:
        ...


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
