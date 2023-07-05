import os.path
from nplinker.logconfig import LogConfig
from nplinker.strains import Strain
from nplinker.utils import list_files
from ..abc import BGCLoaderBase
from ..bgc import BGC
from .mibig_metadata import MibigMetadata


logger = LogConfig.getLogger(__name__)


class MibigBGCLoader:

    def __init__(self, data_dir: str):
        """Parse MIBiG metadata files and return BGC objects

        MIBiG metadata file (json) contains annotations/metadata information
        for each BGC. See https://mibig.secondarymetabolites.org/download.

        Args:
            data_dir(str): Path to the directory of MIBiG metadata json files
        """
        self.data_dir = data_dir
        self._file_dict = self.parse_data_dir(self.data_dir)
        self._metadata_dict = self._parse_metadatas()
        self._bgc_dict = self._parse_bgcs()

    def get_bgc_genome_mapping(self) -> dict[str, str]:
        """Get the mapping from BGC to genome.

        Note that for MIBiG BGC, same value is used for BGC id and genome id.
        Users don't have to provide genome id for MIBiG BGCs in the
        `strain_mappings.json` file.

        Returns:
            dict[str, str]: key is BGC id/accession, value is
                genome id that uses the value of BGC accession.
        """
        return {bid: bid for bid in self._file_dict}

    def get_files(self) -> dict[str, str]:
        """Get the path of all MIBiG metadata json files.

        Returns:
            dict[str, str]: key is metadata file name (BGC accession), value is
                path to the metadata json file
        """
        return self._file_dict

    @staticmethod
    def parse_data_dir(data_dir: str) -> dict[str, str]:
        """Parse metadata directory and return pathes to all metadata json
            files.

        Args:
            data_dir(str): path to the directory of MIBiG metadata json files

        Returns:
            dict[str, str]: key is metadata file name (BGC accession), value is
                 path to the metadata json file
        """
        file_dict = {}
        json_files = list_files(data_dir, prefix="BGC", suffix=".json")
        for file in json_files:
            fname = os.path.splitext(os.path.basename(file))[0]
            file_dict[fname] = file
        return file_dict

    def get_metadatas(self) -> dict[str, MibigMetadata]:
        """Get MibigMetadata objects.

        Returns:
            dict[str, MibigMetadata]: key is BGC accession (file name) and
                value is :class:`nplinker.genomics.mibig.MibigMetadata` object
        """
        return self._metadata_dict

    def _parse_metadatas(self) -> dict[str, MibigMetadata]:
        """Parse all metadata files and return MibigMetadata objects.

        Returns:
            dict[str, MibigMetadata]: key is BGC accession (file name) and
                value is :class:`nplinker.genomics.mibig.MibigMetadata` object
        """
        metadata_dict = {}
        for name, file in self._file_dict.items():
            metadata = MibigMetadata(file)
            metadata_dict[name] = metadata
        return metadata_dict

    def get_bgcs(self) -> dict[str, BGC]:
        """Get BGC objects.

        Returns:
            dict[str, BGC]: key is BGC name and value is
                :class:`nplinker.genomics.BGC` object
        """
        return self._bgc_dict

    def _parse_bgcs(self) -> dict[str, BGC]:
        """Parse all metadata files as BGC objects

        Returns:
            dict[str, BGC]: key is BGC accession (file name) and value is
                BGC object
        """
        bgc_dict = {}
        for name, file in self._file_dict.items():
            bgc = parse_bgc_metadata_json(file)
            bgc_dict[name] = bgc
        return bgc_dict


def parse_bgc_metadata_json(file: str) -> BGC:
    """Parse MIBiG metadata file and return BGC object

    Args:
        file(str): Path to the MIBiG metadata json file

    Returns:
        BGC: :class:`nplinker.genomics.BGC` object
    """
    metadata = MibigMetadata(file)
    mibig_bgc = BGC(metadata.mibig_accession, *metadata.biosyn_class)
    mibig_bgc.mibig_bgc_class = metadata.biosyn_class
    mibig_bgc.strain = Strain(metadata.mibig_accession)
    return mibig_bgc


# register as virtual class to prevent metaclass conflicts
BGCLoaderBase.register(MibigBGCLoader)
