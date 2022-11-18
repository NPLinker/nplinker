import zipfile
from pathlib import Path
import numpy
from typing_extensions import Self
from nplinker.metabolomics.gnps.gnps_loader import GNPSLoader
from . import DATA_DIR


class GNPSLoaderBuilder:

    def __init__(self):
        self._filepath = None

    def with_filepath(self, filepath: Path) -> Self:
        self._filepath = filepath
        return self
    
    def build(self) -> GNPSLoader:
        return GNPSLoader(self._filepath)

def test_default():
    sut = GNPSLoaderBuilder().build()
    assert sut is not None


def test_has_zipfile():
    filepath = DATA_DIR / 'metabolomics_data.zip'
    sut = GNPSLoaderBuilder().with_filepath(filepath).build()
    actual = sut.data()
    
    expected = zipfile.ZipFile(filepath)
    numpy.testing.assert_array_equal(actual.namelist(), expected.namelist())
