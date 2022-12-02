import filecmp
import os
from pathlib import Path
import zipfile
import numpy
import pytest
from nplinker import utils

from nplinker.pairedomics.downloader import Downloader
from nplinker.pairedomics.downloader import _generate_gnps_download_url
from nplinker.pairedomics.downloader import _execute_download
from . import DATA_DIR


@pytest.fixture
def gnps_url():
    return _generate_gnps_download_url("c22f44b14a3d450eb836d607cb9521bb")



@pytest.mark.parametrize("expected", [
    Path(os.getenv('HOME'), 'nplinker_data', 'pairedomics'),
    pytest.lazy_fixture('tmp_path')
])
def test_default(expected):
    gnps_id = "MSV000079284"

    sut = Downloader(gnps_id, local_cache=str(expected))


    assert sut is not None
    assert sut.gnps_massive_id == gnps_id
    assert sut.local_cache == str(expected)

    assert sut.local_download_cache == str(expected / 'downloads')
    assert sut.project_download_cache == str(expected / 'downloads' / gnps_id)

    assert sut.local_file_cache == str(expected / 'extracted')
    assert sut.project_file_cache == str(expected / 'extracted'/ gnps_id)
    assert sut.strain_mappings_file == str(expected / 'extracted'/ gnps_id / 'strain_mappings.csv')
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'antismash'))
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'bigscape'))

    assert sut.all_project_json_file == str(expected / 'all_projects.json')
    assert sut.project_json_file == str(expected / f"{gnps_id}.json")



def test_download_metabolomics_zipfile(tmp_path):
    sut = Downloader("MSV000079284", local_cache=tmp_path)
    sut._download_metabolomics_zipfile("c22f44b14a3d450eb836d607cb9521bb")
    expected_path = os.path.join(sut.project_download_cache, 'metabolomics_data.zip')

    assert os.path.exists(expected_path)
    assert (Path(sut.project_file_cache) / "networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo").is_file()
    assert (Path(sut.project_file_cache) / "clusterinfosummarygroup_attributes_withIDs_withcomponentID/d69356c8e5044c2a9fef3dd2a2f991e1.tsv").is_file()
    assert (Path(sut.project_file_cache) / "spectra/METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf").is_file()


def test_download_metabolomics_zipfile_v2(tmp_path):
    sut = Downloader("MSV000079284", local_cache=tmp_path)
    sut._download_metabolomics_zipfile_v2("c22f44b14a3d450eb836d607cb9521bb")
    expected_path = os.path.join(sut.project_download_cache, 'c22f44b14a3d450eb836d607cb9521bb.zip')

    assert os.path.exists(expected_path)
    assert (Path(sut.project_file_cache) / "molecular_families.pairsinfo").is_file()
    assert (Path(sut.project_file_cache) / "file_mappings.tsv").is_file()
    assert (Path(sut.project_file_cache) / "spectra.mgf").is_file()


def test_generate_gnps_download_url():
    gnps_task_id = "c22f44b14a3d450eb836d607cb9521bb"
    expected = 'https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=c22f44b14a3d450eb836d607cb9521bb&view=download_clustered_spectra'
    actual = _generate_gnps_download_url(gnps_task_id)
    assert actual == expected


def test_execute_download(gnps_url: str, tmp_path: Path):
    outpath = tmp_path / 'metabolomics_data.zip'
    _execute_download(gnps_url, outpath)
    assert os.path.exists(outpath)


def test_download_gnps_data(tmp_path):
    gnps_task_id = "c22f44b14a3d450eb836d607cb9521bb"
    sut = Downloader("MSV000079284", local_cache=tmp_path / 'actual')
    actual = sut._load_gnps_data(gnps_task_id)
    
    expected = zipfile.ZipFile(DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip")

    actual.extract("networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo", tmp_path / "actual")
    expected.extract("networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo", tmp_path / "expected")

    assert filecmp.cmp(
        tmp_path / "actual/networkedges_selfloop" / "6da5be36f5b14e878860167fa07004d6.pairsinfo",
        tmp_path / "expected/networkedges_selfloop" / "6da5be36f5b14e878860167fa07004d6.pairsinfo",
        shallow=False
    )


def test_extract_metabolomics_data(tmp_path):
    sut = Downloader("MSV000079284", local_cache=tmp_path)
    archive = zipfile.ZipFile(DATA_DIR / "ProteoSAFe-METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra.zip")
    sut._extract_metabolomics_data(archive)

    assert (Path(sut.project_file_cache) / "networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo").is_file()
    assert (Path(sut.project_file_cache) / "clusterinfosummarygroup_attributes_withIDs_withcomponentID/d69356c8e5044c2a9fef3dd2a2f991e1.tsv").is_file()
    assert (Path(sut.project_file_cache) / "spectra/METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf").is_file()

