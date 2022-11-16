from abc import ABC, abstractmethod
from collections.abc import Sequence

from nplinker.metabolomics.spectrum import Spectrum


class ISpectrumLoader(ABC):
    
    @abstractmethod
    def spectra(self) -> Sequence[Spectrum]:
        ...
