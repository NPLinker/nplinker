import pytest
from nplinker.genomics.mibig import MibigMetadata
from . import DATA_DIR


@pytest.mark.parametrize("version", ["v1.4", "v3.1"])
class TestMibigMetadata:

    @pytest.fixture
    def json_file(self, version):
        json_file = DATA_DIR / "mibig" / f"BGC0000001_{version}.json"
        yield json_file

    @pytest.fixture
    def metadata(self, json_file):
        yield MibigMetadata(json_file)

    def test_init(self, metadata, json_file):
        assert metadata.file == json_file
        assert isinstance(metadata.metadata, dict)

    def test_mibig_accession(self, metadata):
        assert metadata.mibig_accession == "BGC0000001"

    def test_biosyn_class(self, metadata):

        assert metadata.biosyn_class == "Polyketide"
