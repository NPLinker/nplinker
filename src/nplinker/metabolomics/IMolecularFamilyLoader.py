from abc import ABC, abstractmethod
from typing import Iterable

from nplinker.metabolomics.molecular_family import MolecularFamily


class IMolecularFamilyLoader(ABC):
    
    @abstractmethod
    def families(self) -> Iterable[MolecularFamily]:
        ...
