import zipfile
from pathlib import Path
import numpy
from typing_extensions import Self
from nplinker.metabolomics.gnps.gnps_downloader import GNPSDownloader
from . import DATA_DIR


class GNPSDownloaderBuilder:
    def __init__(self):
        self._task_id = None
        self._download_root = None
        pass
    
    def with_task_id(self, task_id: str) -> Self:
        self._task_id = task_id
        return self

    def with_download_root(self, download_root: Path) -> Self:
        self._download_root = download_root
        return self
    
    def build(self) -> GNPSDownloader:
        return GNPSDownloader(self._task_id, self._download_root)
    


def test_has_gnps_task_id():
    sut = GNPSDownloaderBuilder().with_task_id("c22f44b14a3d450eb836d607cb9521bb").build()
    assert sut.task_id() == "c22f44b14a3d450eb836d607cb9521bb"


def test_has_url():
    sut = GNPSDownloaderBuilder().with_task_id("c22f44b14a3d450eb836d607cb9521bb").build()
    assert sut.url() == 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=c22f44b14a3d450eb836d607cb9521bb&view=download_clustered_spectra'


def test_downloads_file(tmp_path):
    outpath = tmp_path.joinpath("c22f44b14a3d450eb836d607cb9521bb.zip")
    sut = GNPSDownloader("c22f44b14a3d450eb836d607cb9521bb", tmp_path)
    sut.download()
    actual = zipfile.ZipFile(outpath)

    expected = zipfile.ZipFile(DATA_DIR / 'metabolomics_data.zip')
    numpy.testing.assert_array_equal(actual.namelist(), expected.namelist())