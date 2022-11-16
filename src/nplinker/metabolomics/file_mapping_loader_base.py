from abc import ABC, abstractmethod


class FileMappingLoaderBase(ABC):

    @abstractmethod
    def mapping(self) -> dict[int, list[str]]:
        ...