import os.path
import pytest
from nplinker.genomics import BGC
from nplinker.genomics import BGCLoaderBase
from nplinker.genomics.mibig import MibigLoader
from nplinker.genomics.mibig import download_and_extract_mibig_metadata
from nplinker.genomics.mibig import parse_bgc_metadata_json
from nplinker.genomics.mibig.mibig_metadata import MibigMetadata
from .. import DATA_DIR


class TestMibigBGCLoader:
    @pytest.fixture
    def data_dir(self, tmp_path):
        download_root = tmp_path / "download"
        extract_path = tmp_path / "metadata"
        download_root.mkdir()
        extract_path.mkdir()
        download_and_extract_mibig_metadata(download_root, extract_path)
        yield str(extract_path)

    @pytest.fixture
    def loader(self, data_dir):
        loader = MibigLoader(data_dir)
        yield loader

    def test_abc(self, loader):
        assert issubclass(MibigLoader, BGCLoaderBase)
        assert isinstance(loader, BGCLoaderBase)

    def test_init(self, loader, data_dir):
        assert loader.data_dir == data_dir

    def test_get_strain_bgc_mapping(self, loader):
        mapping = loader.get_strain_bgc_mapping()
        assert isinstance(mapping, dict)
        assert len(mapping) == 2502
        for bid in mapping:
            assert bid == mapping[bid]

    def test_get_files(self, loader):
        files = loader.get_files()
        assert isinstance(files, dict)
        assert len(files) == 2502  # MIBiG v3.1 has 2502 BGCs
        assert "BGC0000001" in files
        assert "BGC0000246" not in files
        assert isinstance(files["BGC0000001"], str)
        assert os.path.exists(files["BGC0000001"])

    def test_parse_data_dir(self, data_dir):
        files = MibigLoader.parse_data_dir(data_dir)
        assert isinstance(files, dict)
        assert len(files) == 2502  # MIBiG v3.1 has 2502 BGCs
        assert "BGC0000001" in files
        assert "BGC0000246" not in files
        assert isinstance(files["BGC0000001"], str)
        assert os.path.exists(files["BGC0000001"])

    def test_get_metadatas(self, loader):
        metadatas = loader.get_metadatas()
        assert isinstance(metadatas, dict)
        assert len(metadatas) == 2502  # MIBiG v3.1 has 2502 BGCs
        assert "BGC0000001" in metadatas
        assert "BGC0000246" not in metadatas
        assert isinstance(metadatas["BGC0000001"], MibigMetadata)

    def test_get_bgcs(self, loader):
        bgcs = loader.get_bgcs()
        assert isinstance(bgcs, dict)
        assert len(bgcs) == 2502  # MIBiG v3.1 has 2502 BGCs
        assert "BGC0000001" in bgcs
        assert "BGC0000246" not in bgcs
        assert isinstance(bgcs["BGC0000001"], BGC)


def test_parse_bgc_metadata_json():
    json_file = DATA_DIR / "mibig" / "BGC0000001_v3.1.json"
    bgc = parse_bgc_metadata_json(str(json_file))
    assert isinstance(bgc, BGC)
    assert bgc.bgc_id == "BGC0000001"
    assert bgc.mibig_bgc_class == ("Polyketide",)
