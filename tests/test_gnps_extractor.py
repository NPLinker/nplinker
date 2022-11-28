import filecmp
import zipfile
from pathlib import Path
import numpy
from typing_extensions import Self

import pytest
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from . import DATA_DIR


class GNPSExtractorBuilder:

    def __init__(self):
        self._filepath = None
        self._extract_path = None

    def with_filepath(self, filepath: Path) -> Self:
        self._filepath = filepath
        return self

    def with_extract_path(self, extract_path: Path) -> Self:
        self._extract_path = extract_path
        return self
    
    def build(self) -> GNPSExtractor:
        return GNPSExtractor(self._filepath, self._extract_path)

@pytest.fixture(params=[
    "ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip",
    "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip"
])
def extractor(request, tmp_path) -> GNPSExtractor:
    filepath = DATA_DIR / request.param
    sut = GNPSExtractorBuilder().with_filepath(filepath).with_extract_path(tmp_path).build()
    return sut
    
def test_default():
    sut = GNPSExtractorBuilder().build()
    assert sut is not None


def test_has_zipfile():
    filepath = DATA_DIR / 'metabolomics_data.zip'
    sut = GNPSExtractorBuilder().with_filepath(filepath).build()
    actual = sut.data()
    
    expected = zipfile.ZipFile(filepath)
    numpy.testing.assert_array_equal(actual.namelist(), expected.namelist())


def test_has_extract_path(tmp_path):
    filepath = DATA_DIR / 'metabolomics_data.zip'
    sut = GNPSExtractorBuilder().with_filepath(filepath).with_extract_path(tmp_path).build()
    assert sut.target() == tmp_path


def test_creates_spectra(extractor: GNPSExtractor):
    extractor.extract()

    actual = extractor.target() / "spectra.mgf"
    expected = DATA_DIR / "spectra.mgf"

    assert Path.exists(actual)
    assert filecmp.cmp(actual, expected)


def test_creates_molecular_families(extractor: GNPSExtractor):
    extractor.extract()

    actual = extractor.target() / "molecular_families.pairsinfo"
    expected = DATA_DIR / "edges.pairsinfo"

    assert Path.exists(actual)
    assert filecmp.cmp(actual, expected)


def test_creates_file_mappings(extractor: GNPSExtractor):
    extractor.extract()
    actual = extractor.target() / "file_mappings.tsv"
    expected = DATA_DIR / "nodes.tsv"

    assert Path.exists(actual)
    assert filecmp.cmp(actual, expected)