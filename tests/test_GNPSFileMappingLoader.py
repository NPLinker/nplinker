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


@pytest.mark.parametrize("filename, expected_length, spectrum_id, samples", [
    [DATA_DIR / "nodes_fbmn_mwe.csv", 9, 301, ["20210623_10_9A_1uL.mzML", "20210623_16_9C_1uL.mzML", "20210623_12_5B_1uL.mzML", "20210623_13_9B_1uL.mzML"]],
    [DATA_DIR / "nodes_fbmn_mwe.csv", 9, 1465, ["20210623_10_9A_1uL.mzML", "20210623_16_9C_1uL.mzML"]],
    [DATA_DIR / "nodes_fbmn.csv", 994, 304, ["20210623_10_9A_1uL.mzML", "20210623_16_9C_1uL.mzML", "20210623_13_9B_1uL.mzML"]]
])
def test_load_mapping_fbmn(filename, expected_length, spectrum_id, samples):
    sut = GNPSFileMappingLoader(filename)
    sut.load_mapping_fbmn()
    actual = sut.mapping()

    assert actual[spectrum_id] == samples
    assert len(actual) == expected_length
