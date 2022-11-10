import pytest

from nplinker.metabolomics.GNPSFileMappingLoader import GNPSFileMappingLoader
from tests import DATA_DIR

@pytest.fixture
def loader() -> GNPSFileMappingLoader:
    filename = DATA_DIR / "nodes.tsv"
    return GNPSFileMappingLoader(filename)

def test_default(loader):
    assert loader is not None


@pytest.mark.parametrize("filename, expected_length, spectrum_id", [
    [ DATA_DIR / "nodes_mwe.csv", 13, 275],
    [ DATA_DIR / "nodes.tsv", 25935, 223]
])
def test_load_mapping_allfiles(filename, expected_length, spectrum_id):
    sut = GNPSFileMappingLoader(filename)
    sut.load_mapping_allfiles()

    actual = sut.mapping()

    assert len(actual) == expected_length
    assert actual[spectrum_id] == ["26c.mzXML", "26c.mzXML", "26c.mzXML"]