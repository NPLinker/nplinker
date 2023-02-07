from abc import ABC, abstractmethod

from nplinker.metabolomics.molecular_family import MolecularFamily
from collections.abc import Sequence
from nplinker.metabolomics.spectrum import Spectrum


class SpectrumLoaderBase(ABC):
    
    @abstractmethod
    def spectra(self) -> Sequence[Spectrum]:
        ...

class MolecularFamilyLoaderBase(ABC):
    
    @abstractmethod
    def families(self) -> Sequence[MolecularFamily]:
        ...


class FileMappingLoaderBase(ABC):

    @abstractmethod
    def mapping(self) -> dict[int, list[str]]:
        ...


class AnnotationLoaderBase(ABC):

    @abstractmethod
    def get_annotations(self) -> dict[int, dict]:
        ...