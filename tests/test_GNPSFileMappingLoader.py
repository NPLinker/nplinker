from nplinker.metabolomics.GNPSFileMappingLoader import GNPSFileMappingLoader
from tests import DATA_DIR


def test_default():
    filename = DATA_DIR / "nodes.tsv"
    sut = GNPSFileMappingLoader(filename)

    assert sut is not None