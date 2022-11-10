from abc import ABC, abstractmethod
from typing import Dict


class IFileMappingLoader(ABC):

    @abstractmethod
    def mapping(self) -> Dict[int, list[str]]:
        ...