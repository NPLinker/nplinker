from __future__ import annotations
from abc import ABC
from abc import abstractmethod
from os import PathLike
from .bgc import BGC
from .gcf import GCF


class BGCLoaderBase(ABC):
    """Abstract base class for BGC loader."""

    def __init__(self, data_dir: str | PathLike) -> None:
        """Initialize the BGC loader.

        Args:
            data_dir: Path to directory that contains BGC metadata files
                (.json) or full data genbank files (.gbk).
        """
        self.data_dir = str(data_dir)

    @abstractmethod
    def get_files(self) -> dict[str, str]:
        """Get path to BGC files.

        Returns:
            The key is BGC name and value is path to BGC file
        """

    @abstractmethod
    def get_bgcs(self) -> list[BGC]:
        """Get BGC objects.

        Returns:
            A list of BGC objects
        """


class GCFLoaderBase(ABC):
    """Abstract base class for GCF loader."""

    @abstractmethod
    def get_gcfs(self, keep_mibig_only: bool, keep_singleton: bool) -> list[GCF]:
        """Get GCF objects.

        Args:
            keep_mibig_only: True to keep GCFs that contain only MIBiG
                BGCs.
            keep_singleton: True to keep singleton GCFs. A singleton GCF
                is a GCF that contains only one BGC.

        Returns:
            A list of GCF objects
        """
