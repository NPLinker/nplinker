
import pathlib
import zipfile

import numpy
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from . import DATA_DIR


def test_has_gnps_task_id():
    sut = GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb")
    assert sut.task_id() == "c22f44b14a3d450eb836d607cb9521bb"


def test_has_url():
    sut = GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb")
    assert sut.url() == 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=c22f44b14a3d450eb836d607cb9521bb&view=download_clustered_spectra'


def test_downloads_file(tmp_path):
    outpath = tmp_path / "data.zip"
    sut = GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb", outpath)
    actual = zipfile.ZipFile(outpath)

    expected = zipfile.ZipFile(DATA_DIR / 'metabolomics_data.zip')
    numpy.testing.assert_array_equal(actual.namelist(), expected.namelist())