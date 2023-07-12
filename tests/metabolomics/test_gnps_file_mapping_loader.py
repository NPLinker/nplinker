import pytest

from nplinker.metabolomics.gnps.gnps_file_mapping_loader import GNPSFileMappingLoader
from tests import DATA_DIR

@pytest.fixture
def loader() -> GNPSFileMappingLoader:
    filename = DATA_DIR / "nodes.tsv"
    return GNPSFileMappingLoader(filename)

def test_default(loader):
    assert loader is not None

@pytest.mark.parametrize("filename, expected_length, spectrum_id, samples", [
    [DATA_DIR / "nodes_fbmn_mwe.csv", 9, "301", ["20210623_10_9A_1uL.mzML", "20210623_16_9C_1uL.mzML", "20210623_12_5B_1uL.mzML", "20210623_13_9B_1uL.mzML"]],
    [DATA_DIR / "nodes_fbmn_mwe.csv", 9, "1465", ["20210623_10_9A_1uL.mzML", "20210623_16_9C_1uL.mzML"]],
    [DATA_DIR / "nodes_fbmn.csv", 994, "304", ["20210623_10_9A_1uL.mzML", "20210623_16_9C_1uL.mzML", "20210623_13_9B_1uL.mzML"]],
    [DATA_DIR / "nodes_mwe.csv", 13, "275", ["26c.mzXML", "26c.mzXML", "26c.mzXML"]],
    [DATA_DIR / "nodes.tsv", 25935, "223", ["26c.mzXML", "26c.mzXML", "26c.mzXML"]]
])
def test_load_mapping(filename, expected_length, spectrum_id, samples):
    sut = GNPSFileMappingLoader(str(filename))
    actual = sut.mapping()

    assert actual[spectrum_id] == samples
    assert len(actual) == expected_length
