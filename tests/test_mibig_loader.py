import shutil
import tempfile
import pytest
from nplinker.genomics.mibig import MibigBGCLoader
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.genomics.mibig.mibig_bgc import MibigBGC
from nplinker.genomics.mibig.mibig_metadata import MibigMetadata


class TestMibigBGCLoader:

    @pytest.fixture
    def database(self):
        download_root = tempfile.mkdtemp()
        extract_path = tempfile.mkdtemp()
        download_and_extract_mibig_metadata(download_root, extract_path)
        yield extract_path
        shutil.rmtree(download_root)

    def test_init(self, database):
        loader = MibigBGCLoader(database)

        assert loader.database == database
        assert isinstance(loader.bgcs, dict)
        assert len(loader.bgcs) == 2502  # MIBiG v3.1 has 2502 BGCs
        assert "BGC0000001" in loader.bgcs
        assert "BGC0000246" not in loader.bgcs
        assert isinstance(loader.bgcs["BGC0000001"], MibigBGC)

    def test_parse_metadata_database(self, database):
        metadatas = MibigBGCLoader.parse_metadata_database(database)
        assert isinstance(metadatas, dict)
        assert len(metadatas) == 2502  # MIBiG v3.1 has 2502 BGCs
        assert "BGC0000001" in metadatas
        assert "BGC0000246" not in metadatas
        assert isinstance(metadatas["BGC0000001"], MibigMetadata)
