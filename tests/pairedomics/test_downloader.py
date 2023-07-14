import os
from pathlib import Path
import pytest
from pytest_lazyfixture import lazy_fixture
from requests.exceptions import ReadTimeout
from nplinker.globals import STRAIN_MAPPINGS_FILENAME
from nplinker.pairedomics.downloader import PODPDownloader


@pytest.mark.parametrize("expected", [
    Path(os.getenv('HOME'), 'nplinker_data', 'pairedomics'),
    lazy_fixture('tmp_path')
])
def test_default(expected: Path):
    gnps_id = "MSV000079284"

    sut = PODPDownloader(gnps_id, working_dir=str(expected))

    assert sut.gnps_massive_id == gnps_id
    assert sut.working_dir == str(expected)

    assert sut.downloads_dir == str(expected / 'downloads')
    assert sut.project_download_cache == str(expected / 'downloads' / gnps_id)

    assert sut.results_dir == str(expected / 'extracted')
    assert sut.project_file_cache == str(expected / 'extracted'/ gnps_id)
    assert sut.strain_mappings_file == str(expected / 'extracted'/ gnps_id / STRAIN_MAPPINGS_FILENAME)
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'antismash'))
    assert os.path.exists(str(expected / 'extracted'/ gnps_id / 'bigscape'))

    assert sut.all_project_json_file == str(expected / 'all_projects.json')
    assert sut.project_json_file == str(expected / f"{gnps_id}.json")

def test_download_metabolomics_zipfile(tmp_path):
    sut = PODPDownloader("MSV000079284", working_dir=tmp_path)
    try:
        sut._download_metabolomics_zipfile("c22f44b14a3d450eb836d607cb9521bb")
        expected_path = os.path.join(sut.project_download_cache, 'metabolomics_data.zip')

        assert os.path.exists(expected_path)
        assert (Path(sut.project_file_cache) / "networkedges_selfloop/6da5be36f5b14e878860167fa07004d6.pairsinfo").is_file()
        assert (Path(sut.project_file_cache) / "clusterinfosummarygroup_attributes_withIDs_withcomponentID/d69356c8e5044c2a9fef3dd2a2f991e1.tsv").is_file()
        assert (Path(sut.project_file_cache) / "spectra/METABOLOMICS-SNETS-c22f44b1-download_clustered_spectra-main.mgf").is_file()
    except ReadTimeout:
        pytest.skip("GNPS is down")


def test_download_metabolomics_zipfile_scenario2(tmp_path):
    sut = PODPDownloader("MSV000079284", working_dir=tmp_path)
    try:
        sut._download_metabolomics_zipfile("c22f44b14a3d450eb836d607cb9521bb")
        expected_path = os.path.join(sut.project_download_cache, 'c22f44b14a3d450eb836d607cb9521bb.zip')

        assert os.path.exists(expected_path)
        assert (Path(sut.project_file_cache) / "molecular_families.pairsinfo").is_file()
        assert (Path(sut.project_file_cache) / "file_mappings.tsv").is_file()
        assert (Path(sut.project_file_cache) / "spectra.mgf").is_file()
    except ReadTimeout:
        pytest.skip("GNPS is down")
