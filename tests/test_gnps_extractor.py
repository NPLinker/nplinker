import filecmp
import os
from tempfile import gettempdir
import zipfile
from pathlib import Path
import numpy
from typing_extensions import Self

import pytest
from nplinker.metabolomics.gnps.gnps_extractor import GNPSExtractor
from nplinker.utils import extract_archive
from . import DATA_DIR


class GNPSExtractorBuilder:

    def __init__(self):
        self._filepath = DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip"
        self._extract_path = gettempdir()

    def with_filepath(self, filepath: Path) -> Self:
        self._filepath = filepath
        return self

    def with_extract_path(self, extract_path: Path) -> Self:
        self._extract_path = extract_path
        return self
    
    def build(self) -> GNPSExtractor:
        return GNPSExtractor(self._filepath, self._extract_path)


def assert_extraction_success(filename: str, outdir: Path, actual: Path):
    expected = outdir / filename
    assert Path.exists(actual)
    assert filecmp.cmp(actual, expected, shallow=False)


def _unpack(archive: Path):
    filepath = DATA_DIR / archive
    outdir = DATA_DIR / filepath.stem
    extract_archive(filepath, outdir)
    return filepath, outdir

    
def test_default():
    sut = GNPSExtractorBuilder().build()
    assert sut is not None


def test_has_zipfile():
    filepath = DATA_DIR / 'ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip'
    sut = GNPSExtractorBuilder().with_filepath(filepath).build()
    actual = sut.data()
    
    expected = zipfile.ZipFile(filepath)
    numpy.testing.assert_array_equal(actual.namelist(), expected.namelist())


def test_has_extract_path(tmp_path):
    filepath = DATA_DIR / 'ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip'
    sut = GNPSExtractorBuilder().with_filepath(filepath).with_extract_path(tmp_path).build()
    assert sut.target() == str(tmp_path)


@pytest.mark.parametrize("archive, filename", [
    ["ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip", "METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf"],
    ["ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip", "spectra/specs_ms.mgf"]
])
def test_creates_spectra(archive: Path, filename: str, tmp_path: Path):
    filepath, outdir = _unpack(archive)

    sut = GNPSExtractorBuilder().with_filepath(filepath).with_extract_path(tmp_path).build()
    sut._extract_spectra()
    actual = Path(sut.target()) / "spectra.mgf"

    assert_extraction_success(filename, outdir, actual)


@pytest.mark.parametrize("archive, filename", [
    ["ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip", "networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo"],
    ["ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip", "networkedges_selfloop/c74fec018736475483e9c8b05e230cce..selfloop"]
])
def test_creates_molecular_families(archive: Path, filename: str, tmp_path: Path):
    filepath, outdir = _unpack(archive)

    sut = GNPSExtractorBuilder().with_filepath(filepath).with_extract_path(tmp_path).build()
    sut._extract_molecular_families()
    actual = Path(sut.target()) / "molecular_families.pairsinfo"
    
    assert_extraction_success(filename, outdir, actual)


@pytest.mark.parametrize("archive, filename", [
    ["ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip", "clusterinfosummarygroup_attributes_withIDs_withcomponentID/d69356c8e5044c2a9fef3dd2a2f991e1.tsv"],
    ["ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-92036537-download_cytoscape_data.zip", "quantification_table_reformatted/1a12f6fbd2ca4e099ec56bdaea56368f.csv"]
])
def test_creates_file_mappings(archive: Path, filename: str, tmp_path: Path):
    filepath, outdir = _unpack(archive)

    sut = GNPSExtractorBuilder().with_filepath(filepath).with_extract_path(tmp_path).build()
    sut._extract_file_mappings()
    actual = Path(sut.target()) / ("file_mappings" + str(Path(filename).suffix))
    
    assert_extraction_success(filename, outdir, actual)