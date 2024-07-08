from __future__ import annotations
import logging
from os import PathLike
from pathlib import Path
from nplinker.strain import Strain
from nplinker.utils import list_files
from ..abc import BGCLoaderBase
from ..bgc import BGC
from .mibig_metadata import MibigMetadata


logger = logging.getLogger(__name__)


class MibigLoader(BGCLoaderBase):
    """Parse MIBiG metadata files and return BGC objects.

    MIBiG metadata file (json) contains annotations/metadata information
    for each BGC. See https://mibig.secondarymetabolites.org/download.

    The MiBIG accession is used as BGC id and strain name. The loaded BGC
    objects have Strain object as their strain attribute (i.e. `BGC.strain`).
    """

    def __init__(self, data_dir: str | PathLike):
        """Initialize the MIBiG metadata loader.

        Args:
            data_dir: Path to the directory of MIBiG metadata json files

        Examples:
            >>> loader = MibigLoader("path/to/mibig/data/dir")
            >>> loader.data_dir
            'path/to/mibig/data/dir'
            >>> loader.get_bgcs()
            [BGC('BGC000001', 'NRP'), BGC('BGC000002', 'Polyketide')]
        """
        self.data_dir = str(data_dir)
        self._file_dict = self.parse_data_dir(self.data_dir)
        self._metadata_dict = self._parse_metadata()
        self._bgcs = self._parse_bgcs()

    def get_files(self) -> dict[str, str]:
        """Get the path of all MIBiG metadata json files.

        Returns:
            The key is metadata file name (BGC accession), and the value is path to the metadata
            json file
        """
        return self._file_dict

    @staticmethod
    def parse_data_dir(data_dir: str | PathLike) -> dict[str, str]:
        """Parse metadata directory and return paths to all metadata json files.

        Args:
            data_dir: path to the directory of MIBiG metadata json files

        Returns:
            The key is metadata file name (BGC accession), and the value is path to the metadata
            json file
        """
        file_dict = {}
        json_files = list_files(data_dir, prefix="BGC", suffix=".json")
        for file in json_files:
            fname = Path(file).stem
            file_dict[fname] = file
        return file_dict

    def get_metadata(self) -> dict[str, MibigMetadata]:
        """Get MibigMetadata objects.

        Returns:
            The key is BGC accession (file name) and the value is MibigMetadata object
        """
        return self._metadata_dict

    def _parse_metadata(self) -> dict[str, MibigMetadata]:
        """Parse all metadata files and return MibigMetadata objects.

        Returns:
            The key is BGC accession (file name) and the value is MibigMetadata object
        """
        metadata_dict = {}
        for name, file in self._file_dict.items():
            metadata = MibigMetadata(file)
            metadata_dict[name] = metadata
        return metadata_dict

    def get_bgcs(self) -> list[BGC]:
        """Get BGC objects.

        The BGC objects use MiBIG accession as id and have Strain object as
        their strain attribute (i.e. `BGC.strain`), where the name of the Strain
        object is also MiBIG accession.

        Returns:
            A list of BGC objects
        """
        return self._bgcs

    def _parse_bgcs(self) -> list[BGC]:
        """Parse all metadata files as BGC objects.

        Returns:
            A list of BGC objects
        """
        return [parse_bgc_metadata_json(file) for file in self._file_dict.values()]


def parse_bgc_metadata_json(file: str | PathLike) -> BGC:
    """Parse MIBiG metadata file and return BGC object.

    Note that the MiBIG accession is used as the BGC id and strain name. The BGC
    object has Strain object as its strain attribute.

    Args:
        file: Path to the MIBiG metadata json file

    Returns:
        BGC object
    """
    metadata = MibigMetadata(str(file))
    mibig_bgc = BGC(metadata.mibig_accession, *metadata.biosyn_class)
    mibig_bgc.mibig_bgc_class = metadata.biosyn_class
    mibig_bgc.strain = Strain(metadata.mibig_accession)
    return mibig_bgc
