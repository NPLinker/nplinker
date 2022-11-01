from nplinker.logconfig import LogConfig
from nplinker.utils import list_files
from nplinker.strains import Strain
from .mibig_bgc import MibigBGC
from .mibig_metadata import MibigMetadata

logger = LogConfig.getLogger(__file__)


class MibigBGCLoader():

    def __init__(self, database: str):
        """Load MIBiG metadata database and return MibigBGC objects.

        MIBiG metadata database is a group of json files, containing the info
        of annotations/metadata for each BGC.

        Args:
            database(str): Path to the directory of MIBiG BGC metadata database
        """
        self.database = database
        self.bgcs = self._load()

    @staticmethod
    def parse_metadata_database(database: str) -> dict:
        """Parse MIBiG BGC metadata database (json files)

        Args:
            database(str): path to the directory of MIBiG BGC metadata database

        Returns:
            dict: key is MIBiG BGC accession, value is MibigMetadata object
        """
        metadatas = {}
        json_files = list_files(database, suffix=".json", prefix=True)
        logger.info("Found %s MIBiG metadata json files", len(json_files))
        for file in json_files:
            metadata = MibigMetadata(file)
            metadatas[metadata.mibig_accession] = metadata
        return metadatas

    def _load(self) -> dict:
        """Load MIBiG metadata as MibigBGC objects

        Returns:
            dict: key is MIBiG BGC accession and value is MibigBGC object
        """
        metadata_dict = MibigBGCLoader.parse_metadata_database(self.database)
        i = 0
        bgcs = {}
        for name, metadata in metadata_dict.items():
            strain = Strain(metadata.mibig_accession)
            mibig_bgc = MibigBGC(i, strain, metadata.mibig_accession,
                                 metadata.biosyn_class)
            bgcs[name] = mibig_bgc
            i += 1
        return bgcs
