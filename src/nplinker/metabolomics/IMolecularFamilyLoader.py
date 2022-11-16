from abc import ABC, abstractmethod
from collections.abc import Sequence

from nplinker.metabolomics.molecular_family import MolecularFamily


class IMolecularFamilyLoader(ABC):
    
    @abstractmethod
    def families(self) -> Sequence[MolecularFamily]:
        ...
