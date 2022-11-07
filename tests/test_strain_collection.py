from nplinker.strain_collection import StrainCollection
from tests import DATA_DIR


def test_default():
    sut = StrainCollection()
    assert sut is not None


def test_add_from_file():
    filename = DATA_DIR / "strain_mappings.csv"
    sut = StrainCollection()
    sut.add_from_file(filename)

    assert len(sut) == 27
    assert len(sut.lookup_index(1).aliases) == 29
