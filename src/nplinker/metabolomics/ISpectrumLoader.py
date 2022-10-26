from abc import ABC, abstractmethod
from typing import Iterable

from nplinker.metabolomics.spectrum import Spectrum


class ISpectrumLoader(ABC):
    
    @abstractmethod
    def spectra(self) -> Iterable[Spectrum]:
        ...
