import filecmp
import os
from pathlib import Path
import zipfile
import numpy
import pytest
from pytest_lazyfixture import lazy_fixture
from nplinker import utils
from nplinker.pairedomics.downloader import PODPDownloader
from nplinker.pairedomics.downloader import STRAIN_MAPPINGS_FILENAME
from .. import DATA_DIR


@pytest.mark.parametrize("expected", [
    Path(os.getenv('HOME'), 'nplinker_data', 'pairedomics'),
    lazy_fixture('tmp_path')
])
def test_default(expected: Path):
    gnps_id = "MSV000079284"

    sut = PODPDownloader(gnps_id, local_cache=str(expected))

    assert sut.gnps_massive_id == gnps_id
    assert sut.local_cache == str(expected)

    assert sut.local_download_cache == str(expected / 'downloads')
    assert sut.project_download_cache == str(expected / 'downloads' / gnps_id)

    assert sut.local_file_cache == str(expected / 'extracted')
    assert sut.project_file_cache == str(expected / 'extracted'/ gnps_id)
    assert sut.strain_mappings_file == str(expected / 'extracted'/ gnps_id / STRAIN_MAPPINGS_FILENAME)
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'antismash'))
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'bigscape'))

    assert sut.all_project_json_file == str(expected / 'all_projects.json')
    assert sut.project_json_file == str(expected / f"{gnps_id}.json")

def test_download_metabolomics_zipfile(tmp_path):
    sut = PODPDownloader("MSV000079284", local_cache=tmp_path)
    sut._download_metabolomics_zipfile("c22f44b14a3d450eb836d607cb9521bb")
    expected_path = os.path.join(sut.project_download_cache, 'metabolomics_data.zip')

    assert os.path.exists(expected_path)
    assert (Path(sut.project_file_cache) / "networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo").is_file()
    assert (Path(sut.project_file_cache) / "clusterinfosummarygroup_attributes_withIDs_withcomponentID/d69356c8e5044c2a9fef3dd2a2f991e1.tsv").is_file()
    assert (Path(sut.project_file_cache) / "spectra/METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf").is_file()


def test_download_metabolomics_zipfile(tmp_path):
    sut = PODPDownloader("MSV000079284", local_cache=tmp_path)
    sut._download_metabolomics_zipfile("c22f44b14a3d450eb836d607cb9521bb")
    expected_path = os.path.join(sut.project_download_cache, 'c22f44b14a3d450eb836d607cb9521bb.zip')

    assert os.path.exists(expected_path)
    assert (Path(sut.project_file_cache) / "molecular_families.pairsinfo").is_file()
    assert (Path(sut.project_file_cache) / "file_mappings.tsv").is_file()
    assert (Path(sut.project_file_cache) / "spectra.mgf").is_file()
