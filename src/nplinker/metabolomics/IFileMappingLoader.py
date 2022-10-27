from abc import ABC, abstractmethod


class IFileMappingLoader(ABC):

    @abstractmethod
    def mapping(self) -> dict:
        ...