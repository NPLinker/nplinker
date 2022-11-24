from abc import ABC
from abc import abstractmethod
from .bgc import BGC


class BGCLoaderBase(ABC):

    def __init__(self, data_dir: str):
        """Abstract base class for BGC loader.

        Args:
            data_dir(str): Path to directory that contains BGC metadata files
                (.json) or full data genbank files (.gbk).
        """
        self.data_dir = data_dir

    @abstractmethod
    def get_files(self) -> dict[str, str]:
        """Get path to BGC files

        Returns:
            dict[str, str]: key is BGC name and value is path to BGC file
        """

    @abstractmethod
    def get_bgcs(self) -> dict[str, BGC]:
        """Get BGC objects

        Returns:
            dict[str, BGC]: key is BGC name and value is
                :class:`~nplinker.genomic.BGC` objects
        """
